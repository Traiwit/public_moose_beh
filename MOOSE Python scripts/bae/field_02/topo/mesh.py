# -*- coding: utf-8 -*-
"""Topology classes for unstructured meshes

This is deprecated. Use it for code recycling only...
"""


import numpy as np
from numba import njit, jit
from itertools import izip, chain

from bae.log_01 import msg
from bae.abq_model_02 import Model

from .base import StaticEntities
from ..numpy_geom import (createRTree, rTreeQuery)

###############################################################################
#{Topology Classes


def _findIndexOfSubset(fullSet, subSet):
    '''
    '''
    useUnique = True

    if useUnique:
        nFull = len(fullSet)


        v = np.hstack((fullSet, subSet))

        uni, uniIdx = np.unique(v, return_inverse=True)

        idxFull = uniIdx[:nFull]
        idxSub  = uniIdx[nFull:]

        out = {}
        for ii in np.unique(idxSub):
            out[uni[ii]] = np.where(ii==idxFull)[0]

        out = [out[sub] for sub in subSet]
    else:
        inter = np.intersect1d(fullSet, subSet)
        out = [np.array([]) if val not in inter
               else np.where(fullSet==val)[0]
               for val in subSet]
    msg('done')
    return out


class Mesh(StaticEntities):
    """Topology Class to load, access and save unstructured meshes.

    Structure::

        UnstructuredMesh
            |
            |- nodes: structured numpy-array of length nPts
            |    |
            |    |- nodes['label']: node labels as np.int-array of shape (nPts,)
            |    |- nodes['xyz']: node coordinates as np.float-array of shape (nPts,3)
            |
            |- elNodes: dictionary, keys <- ElementTypes, values <- elNodesArray
                 |
                 |- elNodesArray: structured numpy-array of length nElems
                         |
                         |- elNodes[elType]['label']: element labels as np.int-array of shape (nElems,)
                         |- elNodes[elType]['nodes']: nodesLookup for elements as np.int-array of shape (nElems,nodesPerElem)


    """

    # dataTypes
    indexType = np.uint32
    floatType = np.float32

    # general Indexlookups
    faceIndex ={
            'C3D4' : [[0,1,2], [0,1,3], [1,2,3], [0,2,3]],
            'C3D10': [[0,5,1,5,2,6],
                      [0,4,1,8,3,7],
                      [1,5,2,9,3,8],
                      [0,6,2,9,3,7],],
            'COH3D6' : [[0,1,2], [3,4,5]],
            }
    faceIndex['C3D10M'] = faceIndex['C3D10']

    linearEq = {'C3D10M' : 'C3D4',
                'C3D10' : 'C3D4',}


    def __init__(self, **kwargs):
        """An unstructured mesh can be crated:
            - empty
            - from an existing abq_model_02.Model holding mesh data
            - from a path to an abaqus-input deck holding mesh data

        @kwarg model: abq_model_02.Model with not-empty elNodes, elTypes and
        nodeCoords dictionary

        @kwarg meshPath: path to abaqus-input file
        """
        model = kwargs.pop('model', None)
        if model:
            self.fromAbqModel(model)
            return

        meshPath = kwargs.pop('meshPath', None)
        if meshPath:
            self.fromAbqModel(Model().read(meshPath))


    @classmethod
    def getNodeCoords(cls, elNodes, nodes):
        """Returns the node coordinates for element(labels)

        @param elNodes: an element-label-to-nodes-label-lookup, structured
        numpy-arrayof dtype=[('label',int),('node', int, (nNodesPerElemType,))]
        (e.g. myUnstructuredMesh.elNodes[myElemType])

        @param nodes: node label to node coord lookup, structured numpy-array
        of dtype=[('label',int),('xyz',float, (3,))] (e.g.
        myUnstructuredMesh.nodes), , must contain all node labels requested
        by elNodes

        returns: list of length nNodesPerElemType, each item holds the elemNode
        related coordinates as a (len(elNodes), 3)-shaped numpy array
        """
        elNodes = np.atleast_1d(elNodes)
        pts = []
        for pt in range(elNodes['nodes'].shape[1]):
            eNodes = elNodes['nodes'][:,pt]
            pts.append(nodes['coords'][eNodes])

        return pts


    @classmethod
    def getTetCoordsTransform(cls, elNodes, nodes):
        """Returns the affine Tet-Coordinates-Transformations for elnodes

        @param elNodes: an element-label-to-nodes-label-lookup, structured
        numpy-arrayof dtype=[('label',int),('node', int, (nNodesPerElemType,))]
        (e.g. myUnstructuredMesh.elNodes[myElemType])

        @param nodes: node label to node coord lookup, structured numpy-array
        of dtype=[('label',int),('xyz',float, (3,))]
        (e.g. myUnstructuredMesh.nodes), must contain all node labels requested
        by elNodes

        """

        nodeCoords = cls.getNodeCoords(elNodes, nodes)[:4]

        return tetCoordsTransForm(*nodeCoords)


    @classmethod
    def getEnclosingElementCoords(cls, elNodes, nodes, testPoints):
        """Returns the barycentric coordinates of testPoints in a TetMesh
        described by elNodes and nodes.

        Algorithm:
        ==========
            1) get nodeCoordinates for tetElems (if quadratic only the linear
            counterpart is used)
            2) generate a R-Tree for all elements
            3) for each testPoint find all elements, whos boundingboxes include
            the testPoint
            4) loop these elements and test if testPoint is inside

        @param elNodes: an element-label-to-nodes-label-lookup, structured
        numpy-arrayof dtype=[('label',int),('node', int, (nNodesPerElemType,))]
        (e.g. myUnstructuredMesh.elNodes['C3D6'])

        @param nodes: node label to node coord lookup, structured numpy-array
        of dtype=[('label',int),('xyz',float, (3,))]
        (e.g. myUnstructuredMesh.nodes), must contain all node labels requested
        by elNodes

        @param testPoints: list/numpy-array of coordinates

        @returns: structured numpy array of length = len(testPoints) and
        dtype=[('idx',int), ('coord', float, (4,))], for testpoints not found
        in any given element idx will be np.iinfo(np.uint32).max == 4294967295
        and coord will be [np.nan, np.nan, np.nan]
        """

        nPts = len(testPoints)
        msg('Find surrounding elements for %d points.' % nPts)

        msg('....Collecting nodal coordinates')
        nodeCoords = cls.getNodeCoords(elNodes, nodes)

        nodeCoords = np.asarray(nodeCoords)

        msg('....Finding indices')
        idxs = pointsTetCoordsRTree(nodeCoords, testPoints)

        return idxs


    @classmethod
    def elemFaces(cls, elNodes, elemType, linear=True):
        """Returns the elementFaces

        @param elNodes: an element-label-to-nodes-label-lookup, structured
        numpy-arrayof dtype=[('label',int),('node', int, (nNodesPerElemType,))]
        (e.g. myUnstructuredMesh.elNodes['C3D6'])

        @param elemType: String specifiing the element type

        @param linear: if elemType is quadratic, use the linear counterpart
        """

        if linear:
            try:
                elemType = cls.linearEq[elemType.upper()]
            except KeyError:
                pass

        surfIdxs = np.array(cls.faceIndex[elemType.upper()])

        dType = [('label', cls.indexType),
                 ('surfIdxs', cls.indexType, surfIdxs.shape)]

        surf = np.ones(len(elNodes), dtype=dType)

        surf['label'] = elNodes['label']
        surf['surfIdxs'] = elNodes['nodes'][:,surfIdxs]

        return surf


    @staticmethod
    def indexFromLabel(labels, subLabels):
        inter, ii, jj = np.intersect1d(labels, subLabels,
                                       return_indices=True)

        lookup = dict(zip(inter, ii))
        try:
            return np.array([lookup[subLabel] for subLabel in subLabels])
        except KeyError:
            notFound = np.array([l for l in subLabels
                                 if l not in lookup.keys()])
            raise KeyError(
                "The following labels can't be found: %s" % notFound)


    def fromAbqModel(self, model):
        """Fills UnstructuredMesh object with mesh-data stored in a
        bae.abq_model_02.Model object.

        @param model: abq_model_02.Model with not-empty elNodes, elTypes and
        nodeCoords dictionary

        @returns: self, so you could do:
            >>> uM = UnstructuredMesh().fromAbqModel(myModel)
        """
        try:
            msg('\t rearange nodeCoords from model')
            dTypeNodes = [('label', self.indexType),
                          ('coords', self.floatType, (3,))]
            self.nodes = np.fromiter(model.nodeCoords.iteritems(),
                                     dtype=dTypeNodes,
                                     count=len(model.nodeCoords))
        except AttributeError:
            msg("Can't get nodeCoords from model")

        try:
            msg('\t rearange elementTypes from model')
            dTypeElType = [('label', self.indexType), ('type','S10')]
            elType = np.fromiter(model.elType.iteritems(),
                                 dtype=dTypeElType,
                                 count=len(model.elType))
        except AttributeError:
            msg("Can't get elementTyps from model")

        try:
            msg('\t rearange elementNodes from model')
            elNodes = {}
            for elT in np.unique(elType['type']):
                elNums = elType['label'][elType['type']==elT]

                nodesNum = len(model.elNodes[elNums[0]])
                dTypeElNodes = [('label', self.indexType),
                                ('nodes', self.indexType, (nodesNum,))]
                elNodesIter = ((el, model.elNodes[el]) for el in elNums)
                elNodes[elT] = np.fromiter(elNodesIter,
                                           dtype=dTypeElNodes,
                                           count=len(elNums))


            msg('\t reindex elNodes from labels to np-index')
            for elT, elNode in elNodes.iteritems():
                elNodes[elT]['nodes'] = replaceBy(elNode['nodes'],
                                                  self.nodes['label'],
                                                  np.arange(len(self.nodes)))

            self.elNodes = elNodes

        except AttributeError:
            msg("Can't get elementNodes from model")

        return self

    def toAbqModel(self):
        """Stores the UnstructuredMesh data in a bae.abq_model_02.Model object.
        """
        m = Model()
        m.nodeCoords.update(
            izip(self.nodes['label'], self.nodes['xyz'].tolist()))

        for elType, elNodes in self.elNodes.iteritems():

            labels = elNodes['label']
            nodes = replaceBy(elNodes['nodes'],
                              np.arange(len(self.nodes)),
                              self.nodes['label'])
            dict.update(m.elNodes, izip(labels, nodes.tolist()))
            m.elType.update(dict.fromkeys(labels, elType))

        return m


#} #End Topology Classes


###############################################################################


def tetCoordsTransForm(A, B, C, D):
    '''Returns the node coordinats A, B, C and D of a Tetrahedron to its
    orthogonal system

    Found here U{stackoverflow.com<https://stackoverflow.com/questions/25179693/how-to-check-whether-the-point-is-in-the-tetrahedron-or-not>},
    and vectorized.

    '''
    v1 = B-A
    v2 = C-A
    v3 = D-A

    # mat defines an affine transform from the tetrahedron to the orthogonal
    # system
    swapped = np.moveaxis(np.array((v1,v2,v3,A)), 0, -1)
    toAdd = np.tile([0,0,0,1], (len(A),1,1))
    mat = np.concatenate( (swapped, toAdd), axis=1 )
    # The inverse matrix does the opposite (from orthogonal to tetrahedron)

    M1 = np.linalg.inv(mat)
    return(M1)


@njit
def ptOnSamePlaneSide(v1, v2, v3, v4, p, tol=1E-3):
    """Checks if v1 and p are on the 'positive' side (in direction of normal)
    of plane (v1, v2,v3,v4).

    @param v1: (and v2, v3, v4) points specifiing plane

    @param p: point to check

    @param tol: tolerance

    @returns: True if point is on positive side
    """
    normal = surfaceNormal(v1, v2, v3)

    return (np.dot(normal, v4-v1) * np.dot(normal, p-v1) >= -1*tol)


@njit
def surfaceNormal(v1, v2, v3):
    '''Returns the NormalVector from triangular surface defined by v1, v2, v3.
    The length of normal vector equals the spanned area.
    '''
    rr = v2-v1
    tt = v3-v1

    # np.cross is not implemented in numba yet
    normal = .5 * np.array([[rr[1]*tt[2] - rr[2]*tt[1]],
                            [rr[2]*tt[0] - rr[0]*tt[2]],
                            [rr[0]*tt[1] - rr[1]*tt[0]]]).T
    return normal


@njit
def tetNormals(elNodes):
    permIdxs = [[0,1,2], [0,1,3], [1,2,3], [0,2,3]]

    normals = np.zeros_like(elNodes)
    for ii, permIdx in enumerate(permIdxs):
        normals[ii] = surfaceNormal(elNodes[permIdxs])

    return normals



@jit
def pointToTetCoords(nodes, pts):
    v1, v2, v3, v4 = map(np.atleast_2d, nodes)
    pts = np.atleast_2d(pts)
    # Find the transform matrix from orthogonal to tetrahedron system

    M1 = tetCoordsTransForm(v1, v2, v3, v4)
    # apply the transform to P

    p1s = np.concatenate((pts, np.ones((pts.shape[0], 1))), axis=1)
    coords = M1.dot(p1s.T)

    #last coord from restrictions that sum(coords)=1
    coords[:,-1] = 1 - coords[:,:-1].sum(axis=1)

    return coords


@jit
def _candidateTest(elIdxs, elNodes, points, res, rtol = 1E-3):
    '''Checks if points are in Tet-elems with indexes elIdxs and collects
    elementindex and barycentric coordinates.
    @param elIdxs: index-list of elements to check

    @param elNodes: (all) mesh edge nodes (shape (nEdgeNodes=4,nElems,3d=3))

    @param points: points to find surrounding elements for

    @param res: structured array (dtype: elementIndex, tetCoords) to collect
    result in

    @returns: res

    @note: this subroutine is excluded from pointsTetCoords to:
        1) use numba (straight forward loops)
        2) allow different preselection algorithms (Rtree, kdTree, bruteforce)
    '''
    for ii, (elIdx, point) in enumerate(zip(elIdxs, points)):

        elCoords = pointToTetCoords(elNodes[:4,elIdx], point)

        for jj, elCoord in enumerate(elCoords):
            if np.all((elCoord + rtol) >=0 ) and np.all((elCoord - rtol) <=1):
                res[ii]['idx']  = elIdx[jj]
                res[ii]['coor'] = elCoord.ravel()
                break
    return res


def pointsTetCoordsRTree(elNodes, testPoints, treeBaseName=None,
                         notFoundValue = np.iinfo(np.uint32).max):
    """Calculates the barycentric coordinates of testPoints in tetElements
    which are specified by their element corner coordinates.

    @param elNodes: (all) mesh edge nodes (shape (nEdgeNodes=4,nElems,3d=3))

    @param testPoints: sequence of points to find surrounding elements for

    @param notFoundValue: returned elemnentindex if testPoint is not inside of
    any element, default = np.iinfo(np.uint32).max == 4294967295

    @returns: structured array (dtype: elementIndex, tetCoords)
    """
    msg('Setting up rTree')
    rtree = createRTree(elNodes, treeBaseName=treeBaseName, force=False)

    msg('rTree query')
    candidatesElems = rTreeQuery(rtree, testPoints)

    rtree.close()
    del rtree

    idx = np.ones((testPoints.shape[0],),
                  dtype=[('idx', np.uint32),
                         ('coor', np.float32,(4,))])
    idx[:]['idx'] = notFoundValue
    idx[:]['coor'] = np.nan


    msg('testing candidates')
    idx = _candidateTest(candidatesElems, elNodes, testPoints, idx)

    return idx



def replaceBy(field, oldVals, newVals, strict=True):
    if strict and not len(oldVals) == len(np.unique(oldVals)):
        raise ValueError("If strict=True than the lookupkeys has to be unique")

    field = np.asarray(field)
    oldVals = np.asarray(oldVals)
    newVals = np.asarray(newVals)

    upper = max(field.max(), oldVals.max()) + 1
    lut = np.arange(upper)
    lut[oldVals] = newVals
    """
    @param field:
    @param lookupIn: 
    @param lookupOut: 
    """

    return lut[field.ravel()].reshape(field.shape)


if __name__ == '__main__':
    print 'No syntax errors.'

