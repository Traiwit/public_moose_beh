r"""This module defines classes to create ground support structures like bolts
to be inserted into an Abaqus model.

The classes are subclasses of abq_model_02.Model with all their functionality.

The suggested work flow is to instantiate objects of those classes, then copy,
translate, rotate, renumber and/or merge them to one model which can then be
exported as input file.


Use cases
=========

 - cable bolts: use L{BoltEmbeddedTruss}
 - fully grouted bolts: use L{BoltEmbeddedTrussEqualLen}
 - split sets, dynamic bolts, Yieldlok bolts and so on: use
   L{BoltFrictionSleeve} if you can aford the computational cost.


Glossary
========

 - Collar .. is the end of the bolt that is visible in the tunnel. Also called
   foot point, foot node at times. The opposite end is sometimes called the far
   end of the bolt.

"""

__version__ = "2.13"

_version_history_ = """\
Versions:
=========

1.0 : GP: developed for the yieldlok project 2011. Derived from the Antamina GS
      project in 2010.
gs_mesh_02
2.0 : GP: derived from gs_mesh_01, to be used for the yieldlok project in 2012
2.1 : GP added BoltEmbeddedTruss
2.2 : GP added adjustEmbeddingPos
2.3 : GP removed logfile functionality and introduced bae.log_01
2.4 : GP added BoltEmbeddedTrussEqualLen
2.5 : GP adjustEmbeddingPos accepts (e.g. nodeCoords-) dict as points argument
         and new returnValues argument
2.6 : GP for compatibilty with the BoltFrictionSleeve class added to other bolt
         classes:
         - rotate functions with fake newDirectionId arguments
         - elemCollar, elemBoltEnd ivars (for BoltEmbeddedTrussEqualLen only)
2.7 : GP fixed: adjustEmbeddingPos called wrong partMeshFromElements function,
           Use mesh.getElemCoords to sort out points that don't need to be
           adjusted.
         changed: new correctPointsCallBack argument to adjustEmbeddingPos()
2.8 : GP added addNodeNormalsToMesh
2.9 : GP incompatible change: unified bolt-class attribute name
      bolt_length to boltLength, sleeve_length to sleeveLength,
      bolt_frictionTab to boltFrictionTab,
      sleeve_frictionTab to sleeveFrictionTab
      changed class hierarchy, added BoltBase... class(es)
      added BoltEmbeddedBeamEqualLen class
      changed naming scheme for directions in elsetnames: BOLT_DIR001
2.10 : GP added embedNodesRotConstrained, make BoltEmbeddedTrussEqualLen use it
2.11 : GP make rotary constrained embedding for BoltEmbeddedTrussEqualLen
       optional
2.12 : GP added class BoltDataList, BoltModel.install()
2.13 : GP added class AssemblyModel
"""

from math import log, ceil
from itertools import izip, combinations
from collections import defaultdict
import re
import copy
import csv
import cPickle as pickle

from bae.vecmath_01 import \
    vector, vector_plus, vector_sum,\
    length, norm, vector_scale, vector_modif_scale,\
    vector_rot_x, vector_rot_y, vector_rot_z,\
    dot, cross, trans_vert2dir
from bae.abq_model_02 import Model
from bae.mesh_01 import Mesh
from bae.misc_01 import incrementName, Container, getShortestDistKey,\
    groupIter, findNextRoundNumber

from bae.log_01 import msg


def _parseKwArg(self, kwargs, attrName, default="ParseKwArgsNoDefSpec"):
    try:
        val = kwargs[attrName]
    except KeyError:
        if default=="ParseKwArgsNoDefSpec":
            raise ValueError("%s: mandatory keyword argument %s not"
                             " specified." % (self.__class__, attrName))
        else:
            val = default
    return val


def adjustEmbeddingPos(points, mesh, hostElset=None, safety=1E-3,
                       returnValues="adjustedPtsIds",
                       correctPointsCallBack=None):
    """Adjust the given points such that they all lie within the specified
    host elset.

    @param points: Points to be checked and possibly be adjusted. Might be a
       list of point coordinates or a {point id: point coords}-dictionary like
       abq_model_02.Model.nodeCoords. Note: This list or dict will be updated
       in place, i.e. its values will be exchanged where applicable, i.e. the
       point coordinate triples will be replaced by updated ones.
    @param mesh: abaqus model containing the host elements
    @param hostElset: elset (anything that
       L{getUnionSet<bae.abq_model_02.Model.getUnionSet>} accepts as input)
       defining the host elset. If omitted then the whole mesh is considered.
    @param safety: If a point is being corrected move it into the tet element
       such that the smallest element coordinate will not be smaller than this
       safety value. This is also to account for rounding errors introduced by
       the export of node positions to the abaqus input file.

       Set safety=0.0 if you don't want this feature. With too large values
       like 0.1 the algorithm might fail.

       Note that the value is relative to tet element dimensions because it
       acts on element coordinates in the range of 0..1.
    @param returnValues: A string or tuple of strings specifying what kind of
       diagnostic output the function shall return. Defaults to
       "adjustedPtsIds", in this case the function returns a list of point
       indices/keys that needed adjustment.

       If "maxAdjustDist" then the function returns the maximum adjustment
       distance.

       Might also be a tuple/list containing the aforehead mentioned strings
       in which case the output will be a corresponding tuple as well.
    @param correctPointsCallBack: (optional) This callback function gets a dict
       {point id: new position} of all adjusted points from the current
       iteration as only argument. It's supposed to suggest new starting points
       for some or all of those points. It *must* return a dict {point id: new
       starting point}. Then the adjustment algorithm is restarted with just
       those points. All points not contained in this new dict are considered
       accepted and will be transferred to the points list immediately. All
       points contained in this new dict will be checked again based on their
       updated positions in the new dict. The callback function will be called
       again and again until it returns an empty dict.

    @return: List of point indices/keys that needed adjustment. See
       returnValues argument for other options.

    @Note: Adjusted points are replaced in the points list/dict (in place). If
       you don't want the original list/dict to be modified, make a (shallow)
       copy beforehand.

    @Note: Might easily be modified to also accept a mesh_01.Mesh object as
       mesh argument.

    @Note: This interface is a mess! Needs cleanup in gs_mesh_03!
       Suggestions: ...?
    """

    # strip the mesh to only contain hostElems
    if hostElset is not None:
        hostElset = mesh.getUnionSet("elset", hostElset)
        oldNb = len(mesh.elNodes)
        mesh = mesh.partMeshFromElements(hostElset)
        msg("Stripped the host mesh from %d to %d elements."
            % (oldNb, len(mesh.elNodes)))

    # different initialization of pointCoords depending on type of points
    if isinstance(points, dict):
        pointIds = points.keys()
    else:
        pointIds = range(len(points))
    pointCoords = [points[idx] for idx in pointIds]

    # loop over corrections
    adjustedPtsIds = list()   # list of ids of all adjusted points
    while 1:

        # get element number and element coordinates for each point
        elemCoordsList = mesh.getElemCoords(pointCoords)
        adjustMePtIds = [
            pointIds[i]
            for i, (elem, elCoords) in enumerate(elemCoordsList)
            if elem is None]
        adjustMeCoords = [points[idx] for idx in adjustMePtIds]

        adjustedIdElemCoords = [
            (ptId, elem, elCoords)
            for ptId, (elem, elCoords) in izip(
                adjustMePtIds, mesh.getClosestElemCoords(adjustMeCoords)) ]

        # adjust all points to not lie outside their particular element
        adjustedPts = dict()
        maxAdjustDist = 0.0

        for (ptId, elem, elCoords), oldCoords in izip(
                adjustedIdElemCoords, adjustMeCoords):
            adjusted = False
            for nodeId, xi in enumerate(elCoords):
                if xi >= 0.0:
                    continue
                # reduce other coords in order to have them sum up to 1.0
                factor = (1.0)/(1.0-xi)
                # set current coord to zero = on the tet face
                elCoords[nodeId] = 0.0
                for n2 in range(len(elCoords)):
                    if n2==nodeId:
                        continue
                    elCoords[n2] *= factor

                # set flag: this point has been adjusted
                adjusted = True

            if adjusted:
                # check distance from element face and introduce safety
                if safety>0.0:
                    sumCorrection = 0.0
                    safeIds = list()
                    for nodeId, xi in enumerate(elCoords):
                        if xi<safety:
                            sumCorrection += (safety-xi)
                            elCoords[nodeId] = safety
                        elif xi>3*safety:
                            safeIds.append(nodeId)
                    if len(safeIds)==0:
                        raise ValueError("Safety argument of"
                                         " adjustEmbeddingPos() too large!")
                    factor = (1.0)/(1.0+sumCorrection)
                    for nodeId in safeIds:
                        elCoords[nodeId] *= factor

                msg("Adjusted point <%d>, assumed to be in element %d"
                    % (ptId, elem), debugLevel=10)

                # new point position from adjusted element coords
                newCoords = [
                    sum(weight*x for weight, x in izip(elCoords, allX))
                    for allX in izip(*(mesh.nodeCoords[node]
                                       for node in mesh.elNodes[elem]))
                    ]
                # store in adjustedPts
                adjustedPts[ptId] = newCoords

                # store max adjustment distance
                adjustDist = length(vector(oldCoords, newCoords))
                maxAdjustDist = max(adjustDist, maxAdjustDist)
                msg("... position before adjustment: %s"
                    % oldCoords, debugLevel=10)
                msg("... position after adjustment: %s"
                    % newCoords, debugLevel=10)
                msg("... current adjust distance %g, max so far %g"
                    % (adjustDist, maxAdjustDist), debugLevel=10)


        if len(adjustedPts):
            msg("... %d points need adjustment." % len(adjustedPts),
                debugLevel=10)
        else:
            msg("... no point needs adjustment.", debugLevel=10)

        if correctPointsCallBack and adjustedPts:
            tryAgainPoints = correctPointsCallBack(adjustedPts)
            msg("... callback function accepted %d corrections and delivered"
                " %d new starting points."
                % (len(adjustedPts)-len(tryAgainPoints), len(tryAgainPoints)),
                debugLevel=10)

            # remove points from accepted adjustedPts-dict that got new
            # coordinates to check again from correctPointsCallBack
            for ptId in tryAgainPoints:
                del adjustedPts[ptId]

            # store new coordinates to check again in points and pointCoords
            pointIds = tryAgainPoints.keys()
            pointCoords = [tryAgainPoints[ptId] for ptId in pointIds]

            for ptId in pointIds:
                points[ptId] = tryAgainPoints[ptId]
                msg("... new starting point for point %d: %s"
                    % (ptId, tryAgainPoints[ptId]), debugLevel=10)

        else:
            tryAgainPoints = None

        # store (accepted) adjusted coords in points
        for ptId, newCoords in adjustedPts.iteritems():
            points[ptId] = newCoords
        adjustedPtsIds.extend(adjustedPts)

        if tryAgainPoints:
            continue
        else:
            break

    # return something (or not)
    if returnValues=="adjustedPtsIds":
        return adjustedPtsIds
    elif returnValues=="maxAdjustDist":
        return maxAdjustDist
    else:
        vardict = locals()
        try:
            return tuple(vardict[key] for key in returnValues)
        except TypeError:
            # could not interpret returnValues: no return value!
            pass
    return


###############################################################################

class GSMesh090(Model):
    """Mesh with horizontal and vertical wires, aligned to the tunnel axis
    and perpendicular. Consisting of beam or truss elements.

    The mesh starts at y=0 and ends shortly before or at y=contour.length such
    that it fits an integer number of meshSize units. This is meant in the
    tunnel coordinate system

    For the identifiation of individual nodes and elements there are several
    sets of inidices being used:
     - loop indices (y_lidx, s_lidx) (see L{GSMesh090.getNodesOfLoop()})
     - node indices (y_nidx, s_nidx) (see instance variable L{nodeGrid})
     - elements oriented along the tunnel (y_aidx, s_aidx) (see L{elemsAlong})
     - elements perpendicular to the tunnel (y_cidx, s_cidx), see L{elemsCross}

    The relation between those indices is to be seen in the following two
    sketches of single mesh loops with the index pair (y_lidx, s_lidx). It is
    seen from above,the tunnel stretching from left to right. The first (top)
    sketch being for the neg. side, left hand if you look along the tunnel in
    pos y direction, s_lidx<0. The second (bottom) sketch being for the pos.
    side, right hand if you look along the tunnel in pos y direction,
    s_lidx>0::

         (y_nidx, s_nidx)_0                             (y_nidx, s_nidx)_1
         = (y_lidx, s_lidx)                             = (y_lidx+1, s_lidx)
                    (*)----------(y_aidx, s_aidx)_0---------(*)
                     |           = (y_nidx, s_nidx)_0        |
                     |           = (y_lidx, s_lidx)          |
                     |                                       |
                     |                                       |
                     |                                       |
        (y_cidx, s_cidx)_0        loop:                 (y_cidx, s_cidx)_1
        = (y_lidx, s_lidx)        (y_lidx, s_lidx)      = (y_lidx+1, s_lidx)
                     |            s_lidx<0                   |
                     |                                       |
                     |                                       |
                     |           (y_aidx, s_aidx)_1          |
                     |           = (y_nidx, s_nidx)_3        |
                    (*)----------= (y_lidx, s_lidx+1)-------(*)
         (y_nidx, s_nidx)_3                             (y_nidx, s_nidx)_2
         = (y_lidx, s_lidx+1)                           = (y_lidx+1, s_lidx+1)
     |
     |
     +---> y
     |
     V s
         (y_nidx, s_nidx)_0                             (y_nidx, s_nidx)_1
         = (y_lidx, s_lidx-1)                           = (y_lidx+1, s_lidx-1)
                    (*)----------(y_aidx, s_aidx)_0---------(*)
                     |           = (y_nidx, s_nidx)_0        |
                     |           = (y_lidx, s_lidx-1)        |
                     |                                       |
                     |                                       |
                     |                                       |
        (y_cidx, s_cidx)_0        loop:                 (y_cidx, s_cidx)_1
        = (y_lidx, s_lidx)        (y_lidx, s_lidx)      = (y_lidx+1, s_lidx)
                     |            s_lidx>0                   |
                     |                                       |
                     |                                       |
                     |           (y_aidx, s_aidx)_1          |
                     |           = (y_nidx, s_nidx)_3        |
                    (*)----------= (y_lidx, s_lidx)---------(*)
         (y_nidx, s_nidx)_3                             (y_nidx, s_nidx)_2
         = (y_lidx, s_lidx)                             = (y_lidx+1, s_lidx)

    @ivar meshSize: meshspacing
    @ivar numLoopsAlong: number of mesh loops along the tunnel
    @ivar numLoopsCross: number of mesh loops from the top down along the
        contour

    @ivar nodeGrid: nodeGrid[y_nidx][s_nidx] contains the node number of the
        node at the specified grid index. y_nidx starts from 0 at the front of
        the tunnel (y min). max y_nidx is self.numLoopsAlong.

        s_nidx=0 identifies the row of nodes at the very centre of the tunnel
        back. s_nidx=1,2,3,...,self.numLoopsCross identify the parallel rows
        towards pos s. s_nidx=-1,-2,-3,...,-self.numLoopsCross identify the
        parallel rows towards neg s.

    @ivar elemsAlong: elemsAlong[y_aidx][s_aidx] contains the element numbers
        of the beam elements aligned with the tunnel main axis y. y_aidx starts
        from 0 at the front of the tunnel (y min). max y_aidx is
        self.numLoopsAlong-1.

        s_aidx=0 identifies the elements at the very centre of the tunnel
        back. s_aidx=1,2,3,...,self.numLoopsCross identify the parallel rows
        towards pos s. s_aidx=-1,-2,-3,...,-self.numLoopsCross identify the
        parallel rows towards neg s.

    @ivar elemsCross: elemsCross[y_cidx][s_cidx] contains the element numbers
        of the beam elements perpedicular to the tunnel main axis y. y_cidx
        starts from 0 at the front of the tunnel (y min). max y_cidx is
        self.numLoopsAlong.

        s_cidx identifies the elements at the very centre of the tunnel
        back.

        s_cidx=1,2,3,...,self.numLoopsCross identify the elements on the right
        side when looking along the tunnel in pos y direction, starting just
        right of the centre of the tunnel back proceeding towards pos s.
        s_cidx=-1,-2,-3,... identify the elements on the left side proceeding
        towards neg s.
    """
    def __init__(self, *args, **kwargs):
        """Constructor

        @keyword contour: the contour of the mesh of type Tunnel (mandatory)
        @keyword meshSize: mesh size of the ground support mesh (mandatory)
        @keyword sEnd: Half length of the mesh along the circumference of the
           tunnel cross section. If this is not an integer multiple of
           self.meshSize then the value will be truncated to the next smaller
           multiple of self.meshSize. (mandatory)
        @keyword numNodesBetweenCrossings: for the fem mesh, being ignored at
           the moment.
        """

        if len(args)==1 and isinstance(args[0], self.__class__):
            # copy constructor
            oldmodel = args[0]
            Model.__init__(self, oldmodel)

            self.contour = copy.deepcopy(oldmodel.contour)
            self.length = oldmodel.length
            self.numNodesBetweenCrossings = oldmodel.numNodesBetweenCrossings
            self.elemType = oldmodel.elemType

            self.meshSize = oldmodel.meshSize
            self.numLoopsAlong = oldmodel.numLoopsAlong
            self.numLoopsCross = oldmodel.numLoopsCross

            self.nodeGrid = copy.deepcopy(oldmodel.nodeGrid)
            self.elemsAlong = copy.deepcopy(oldmodel.elemsAlong)
            self.elemsCross = copy.deepcopy(oldmodel.elemsCross)

        elif len(args)==0:
            # create new object from kwargs

            contour = _parseKwArg(self, kwargs, "contour")
            meshSize = _parseKwArg(self, kwargs, "meshSize")
            sEnd = _parseKwArg(self, kwargs, "sEnd")

            numNodesBetweenCrossings = _parseKwArg(
                self, kwargs, "numNodesBetweenCrossings", default=0)
            elemType = _parseKwArg(self, kwargs, "elemType", default='B31')

            Model.__init__(self)

            self.contour = contour
            self.length = float(self.contour.length)
            self.numNodesBetweenCrossings = numNodesBetweenCrossings
            self.elemType = elemType

            self.meshSize = meshSize
            self.numLoopsAlong = int(round(self.length/meshSize))
            self.numLoopsCross = int(sEnd/self.meshSize)

            self.nodeGrid = [
                [None, ]*(1+2*self.numLoopsCross)
                for i in range(1+self.numLoopsAlong)
                ]
            self.elemsAlong = [
                [None, ]*(1+2*self.numLoopsCross)
                for i in range(self.numLoopsAlong)
                ]
            self.elemsCross = [
                # index 0 wont be used,needed anyway
                [None, ]*(1+2*self.numLoopsCross)
                for i in range(1+self.numLoopsAlong)
                ]

            nextNode = 1
            nextElem = 1

            # initialize first row
            nextNode, top_nodes = self._newNodesRow(0, nextNode)
            last_node = top_nodes[0]
            elNodes = dict()
            for new_node in top_nodes[1:]:
                elNodes[nextElem] = [last_node,new_node]
                nextElem += 1
                last_node = new_node
            self.updateElems(elNodes=elNodes, elType=self.elemType)

            # update self.elemsAlong
            for elemList, elem in izip(self.elemsAlong, elNodes):
                elemList[0] = elem

            # pos s side
            nextNode, nextElem = self._createMeshOneSide(
                top_nodes, 1, nextNode, nextElem)

            # neg s side
            nextNode, nextElem = self._createMeshOneSide(
                top_nodes, -1, nextNode, nextElem)

    def _newNodesRow(self, s_idx, nextNode):
        """Add a new row of nodes to model at coordinate s_idx*self.meshSize
        and return the new node numbers. Update self.nodeGrid.

        >>> nextNode, newNodeNums = self._newNodesRow(s, nextNode)
        """

        s = s_idx*self.meshSize

        newMeshRowCoords = [
            self.contour.ys2xyzContour(y=i*self.meshSize, s=s)
            for i in range(self.numLoopsAlong+1)
            ]
        newNodeNums = [i+nextNode for i in range(len(newMeshRowCoords))]
        nextNode += len(newMeshRowCoords)
        self.nodeCoords.update(
            (node, coords)
            for node, coords in izip(newNodeNums, newMeshRowCoords)
            )

        # update self.nodeGrid
        for nodeList, node in izip(self.nodeGrid, newNodeNums):
            nodeList[s_idx] = node

        return nextNode, newNodeNums

    def _createMeshOneSide(self, top_nodes, dir_mesh, nextNode, nextElem):
        """
        @param dir_mesh: Must be -1 or 1. If ==1 -> direction pos s,
          if ==-1 -> direction neg s.
        """

        s_ix_list = [dir_mesh*i for i in range(1, self.numLoopsCross)]

        last_nodes = top_nodes

        elNodesMesh = dict()
        for s_ix in s_ix_list:

            nextNode, newNodeNums = self._newNodesRow(s_ix, nextNode)

            # elements between rows
            for i, nnode in enumerate(newNodeNums):
                elNodesMesh[nextElem] = [last_nodes[i],nnode]
                self.elemsCross[i][s_ix] = nextElem
                nextElem += 1

            # elements in the new row
            last_node = newNodeNums[0]
            for i, new_node in enumerate(newNodeNums[1:]):
                elNodesMesh[nextElem] = [last_node,new_node]
                self.elemsAlong[i][s_ix] = nextElem
                nextElem += 1
                last_node = new_node

            last_nodes = newNodeNums

        # evaluate elNodes, really insert the elements...
        self.updateElems(elNodes=elNodesMesh, elType=self.elemType)

        return nextNode, nextElem


    def getNodesOfLoop(self, y_idx, s_idx):
        """Returns a list of the four corner node numbers of the specified mesh
        loop.

        In the following sketch y and s denote the tunnel coordinate system,
        i.e. we look from above, the tunnel streching from left to right. The
        numbers in the corners of the sketched loop indicate the index of the
        corresponding node in the returned list of node numbers::

         +---> y
         |
         V    (O)------(1)
          s    |        |
               |        |
               |        |
              (3)------(2)

        The relation between the loop indices (y_lidx, s_lidx) and the node
        indices (y_nidx, s_nidx) is as follows::

         +---> y
         |
         V    (y_nidx, s_nidx)_0        (y_nidx, s_nidx)_1
          s   = (y_lidx, s_lidx-1)      = (y_lidx+1, s_lidx-1)
                         (*)-------------------(*)
                          |                     |
                          |                     |
                          |   loop:             |
                          |   (y_lidx, s_lidx)  |
                          |                     |
                          |                     |
                         (*)-------------------(*)
              (y_nidx, s_nidx)_3        (y_nidx, s_nidx)_2
              = (y_lidx, s_lidx)        = (y_lidx+1, s_lidx)

        @param y_idx: Loop index along the tunnel main axis. The first row of
           loops at the front (i.e. min y) end of the tunnel has index 0. No
           negative values allowed.
        @param s_idx: Loop index perpendicular to the tunnel main axis. There
           is a row of nodes at s=0, i.e. the very centre of the back of the
           tunnel. The loops attached to those nodes have index 1 on the pos.
           s side and index -1 on neg. s side. This index may not be zero.
        """
        return [ self.nodeGrid[y_nidx][s_nidx] for (y_nidx, s_nidx) in [
                (y_idx,   s_idx-1),
                (y_idx+1, s_idx-1),
                (y_idx,   s_idx  ),
                (y_idx+1, s_idx  ) ] ]


    def getNearestMeshIndex(self, y, s):
        """Returns the index tuple (y_lidx, s_lidx) of the mesh loop for the
        specified (y,s) coordinates.
        For s_lidx there is no zero index, pos side starts with one, neg side
        with -1. y_lidx is always pos or zero.
        """
        s_idx = int(round((s+0.5)/self.meshSize))
        if s_idx<1:
            s_idx -= 1
        y_idx = int(round((y+0.5)/self.meshSize))
        return (y_idx, s_idx)


###############################################################################

class GSQuadShells090(Model):
    """GS model consisting of shell or membrane elements parallel to the tunnel
    surface. Intended to represent a shotcrete layer.

    The mesh starts at y=0 and ends shortly before or at y=contour.length such
    that it fits an integer number of meshSize units. This is meant in the
    tunnel coordinate system

    Individual elements are being identified by an index pair (y_eidx, s_eidx).
    y_eidx=0 identifies the elements at the front of the tunnel (y min).
    max y_eidx is self.numElemsAlong-1.

    s_eidx=1,2,3,...,self.numElemsCross identify the elements on the right
    side when looking along the tunnel in pos y direction, starting just
    right of the centre of the tunnel back proceeding towards pos s.
    s_eidx=-1,-2,-3,... identify the elements on the left side proceeding
    towards neg s.

    @ivar meshSize: meshspacing
    @ivar numElemsAlong: number of mesh loops along the tunnel
    @ivar numElemsCross: number of mesh loops from the top down along the
        contour

    @ivar nodeGrid: nodeGrid[y_nidx][s_nidx] contains the node number of the
        node at the specified grid index. y_nidx starts from 0 at the front of
        the tunnel (y min). max y_nidx is self.numLoopsAlong.

        s_nidx=0 identifies the row of nodes at the very centre of the tunnel
        back. s_nidx=1,2,3,...,self.numLoopsCross identify the parallel rows
        towards pos s. s_nidx=-1,-2,-3,...,-self.numLoopsCross identify the
        parallel rows towards neg s.

    @ivar elemGrid: elemGrid[y_eidx][s_eidx] contains the element number of the
        element identified by (y_eidx, s_eidx).
    """

    def __init__(self, *args, **kwargs):
        """Constructor

        @keyword contour: the contour of the mesh of type Tunnel
        @keyword meshSize: mesh size of the ground support mesh
        @keyword sEnd: Half length of the mesh along the circumference of the
           tunnel cross section. If this is not an integer multiple of
           self.meshSize then the value will be truncated to the next smaller
           multiple of self.meshSize.

        @todo: add argument numNodesBetweenCrossings: for the fem mesh
        """

        if len(args)==1 and isinstance(args[0], self.__class__):
            # copy constructor
            oldmodel = args[0]
            Model.__init__(self, oldmodel)

            self.contour = copy.deepcopy(oldmodel.contour)
            self.length = copy.deepcopy(oldmodel.length)
            self.elemType = copy.deepcopy(oldmodel.elemType)

            self.meshSize = copy.deepcopy(oldmodel.meshSize)
            self.numElemsAlong = copy.deepcopy(oldmodel.numElemsAlong)
            self.numElemsCross = copy.deepcopy(oldmodel.numElemsCross)

            self.nodeGrid = copy.deepcopy(oldmodel.nodeGrid)
            self.elemGrid = copy.deepcopy(oldmodel.elemGrid)

        else:
            # create new object from kwargs

            contour = _parseKwArg(self, kwargs, "contour")
            meshSize = _parseKwArg(self, kwargs, "meshSize")
            sEnd = _parseKwArg(self, kwargs, "sEnd")
            elemType = _parseKwArg(self, kwargs, "elemType", default='S4R')

            Model.__init__(self)

            self.contour = contour
            self.length = float(self.contour.length)
            self.elemType = elemType

            self.meshSize = meshSize
            self.numElemsAlong = int(round(self.length/meshSize))
            self.numElemsCross = int(sEnd/self.meshSize)

            self.nodeGrid = [
                [None, ]*(1+2*self.numElemsCross)
                for i in range(1+self.numElemsAlong)
                ]
            self.elemGrid = [
                # s_index 0 won"t be used, it's needed anyway
                [None, ]*(1+2*self.numElemsCross)
                for i in range(1+self.numElemsAlong)
                ]

            nextNode = 1
            nextElem = 1

            # initialize first row
            nextNode, top_nodes = self._newNodesRow(0, nextNode)

            # pos s side
            nextNode, nextElem = self._createRowOneSide(
                top_nodes, 1, nextNode, nextElem)

            # neg s side
            nextNode, nextElem = self._createRowOneSide(
                top_nodes, -1, nextNode, nextElem)

    def _newNodesRow(self, s_idx, nextNode):
        """Add a new row of nodes to model at coordinate s_idx*self.meshSize
        and return the new node numbers. Update self.nodeGrid.

        >>> nextNode, newNodeNums = self._newNodesRow(s, nextNode)
        """

        s = s_idx*self.meshSize

        newCoords = [
            self.contour.ys2xyzContour(y=i*self.meshSize, s=s)
            for i in range(self.numElemsAlong+1)
            ]
        newNodeNums = [i+nextNode for i in range(len(newCoords))]
        nextNode += len(newCoords)
        self.nodeCoords.update(
            (node, coords)
            for node, coords in izip(newNodeNums, newCoords)
            )

        # update self.nodeGrid
        for nodeList, node in izip(self.nodeGrid, newNodeNums):
            nodeList[s_idx] = node

        return nextNode, newNodeNums

    def _createRowOneSide(self, top_nodes, dir_mesh, nextNode, nextElem):
        """
        @param dir_mesh: Must be -1 or 1. If ==1 -> direction pos s,
          if ==-1 -> direction neg s.
        """

        s_ix_list = [dir_mesh*i for i in range(1, self.numElemsCross)]

        last_nodes = top_nodes

        elNodes = dict()
        for s_ix in s_ix_list:

            nextNode, newNodeNums = self._newNodesRow(s_ix, nextNode)

            if dir_mesh>0:
                for i in range(len(newNodeNums)-1):
                    elNodes[nextElem] = [last_nodes[i],newNodeNums[i],
                                         newNodeNums[i+1],last_nodes[i+1]]
                    self.elemGrid[i][s_ix] = nextElem
                    nextElem += 1
            else:
                for i in range(len(newNodeNums)-1):
                    elNodes[nextElem] = [last_nodes[i],last_nodes[i+1],
                                         newNodeNums[i+1],newNodeNums[i]]
                    self.elemGrid[i][s_ix] = nextElem
                    nextElem += 1

            last_nodes = newNodeNums

        # evaluate elNodes, really insert the elements...
        self.updateElems(elNodes=elNodes, elType=self.elemType)

        return nextNode, nextElem


###############################################################################
###############################################################################

class BoltData(Container):
    pass

class BoltDataList(list):
    """Stores bolt model data for postprocessing (or other) purposes.

    A list of L{BoltData} objects.
    """

    @classmethod
    def fromPickleFile(cls, fileName):
        return pickle.load(open(fileName, "rb"))

    def writeToCsv(self, fileName):
        boltdataFile = csv.writer(open(fileName, "wb"))

        nbElastElemsPerBolt = max(len(bolt.elastElems) for bolt in self)
        boltdataFile.writerow(['x','y','z', 'type', 'eltype',
                               'footnode','elemCollar','elemBoltEnd']
                              + ['boltElems']*nbElastElemsPerBolt)
        for bolt in self:
            boltdataFile.writerow(
                bolt.footpoint
                + [bolt.groupName, bolt.boltElType,
                   bolt.nodeCollar, bolt.elemCollar, bolt.elemBoltEnd]
                + bolt.elastElems)
        del boltdataFile

    def readFromCsv(self, fileName):
        boltdataFile = csv.reader(open(fileName, "rb"))

        headline = boltdataFile.next()
        # should be: x,y,z, type, eltype, footnode,elemCollar,elemBoltEnd, ...
        # ... boltElems, boltElems, boltElems, ...
        for cnt, line in enumerate(boltdataFile):
            if len(line)<9:
                raise ValueError(
                    "ERROR: illegal data format in %s line %d"
                    % (fileName, cnt+1))
            bolt = BoltData(
                footpoint=map(float, line[:3]),
                groupName=line[3],
                boltElType=line[4],
                nodeCollar=int(line[5]),
                elemCollar=int(line[6]),
                elemBoltEnd=int(line[7]),
                elastElems=[int(x) for x in line[8:]])
            self.append(bolt)
        del boltdataFile

    def writeToPickleFile(self, fileName):
        pickle.dump(self, open(fileName, "wb"))

    def getGroups(self):
        """Return a set of all groupName values in self.
        """
        return set(bolt.groupName for bolt in self)

    def filterGroup(self, groupName):
        """Return a new L{BoltDataList} containing only bolts of the
        specified group (Rhino layer).
        """
        res = BoltDataList(
            bolt for bolt in self
            if bolt.groupName == groupName)
        return res


###############################################################################
###############################################################################

class BoltModel(Model):
    """Base class for all bolts.

    @ivar boltLength: total length in meters = "length" constructor-arg
    @ivar nodeCollar: node (number) at the collar
    @ivar nodeBoltEnd: node (number) at the far end
    @ivar nodesElast: a list of node numbers
    @ivar elemsElast: a list of element numbers
    @ivar nodesEmbedded: list of numbers of embedded nodes
    @ivar elemCollar: element number of the element at the collar
    @ivar elemBoltEnd: element number of the element at the far end
    @ivar namePrefix: (exist only after call to L{initMeshEquallen})
    """

    def copyConstructor(self, oldmodel):
        Model.__init__(self, oldmodel)
        self.boltLength = oldmodel.boltLength
        self.nodeCollar = oldmodel.nodeCollar
        self.nodeBoltEnd = oldmodel.nodeBoltEnd
        self.nodesElast = list(oldmodel.nodesElast)
        self.elemsElast = list(oldmodel.elemsElast)
        self.nodesEmbedded = list(oldmodel.nodesEmbedded)
        if hasattr(oldmodel, "elemCollar"):
            self.elemCollar = oldmodel.elemCollar
        if hasattr(oldmodel, "elemBoltEnd"):
            self.elemBoltEnd = oldmodel.elemBoltEnd

    def initMeshEquallen(self, length, footpoint=[0,0,0], direction=[0,0,1],
                         elementLength=0.2, elType=None,
                         namePrefix="", installerNsetName=None,
                         **kwargs):
        """Initialize bolt mesh of equal size elements.
        Used for fully grouted bolts. All nodes are embedded in rock.
        All elements are of the specified length. The actual bolt length will
        be adjusted.

        This function is called by __init__ of derived "Equallen"- subclasses.

        There will be an ELSET "<namePrefix>BOLT" containing all bolt elements.

        There will be an ELSET "<namePrefix>EMBEDDED_ELEMS". Those elements
        need to be embedded in the bulk rock elements:
         >>> *EMBEDDED ELEMENT,HOST ELSET=...,ABSOLUTE EXTERIOR TOLERANCE=...
         >>> <namePrefix>EMBEDDED_ELEMS"

        If elType="B31" then there will be an NSET "<namePrefix>INSTALLER_ROT"
        for all nodes that you might want to fix in their rotation. Note that
        this might have to be stated for each step separately.
         >>> *BOUNDARY, TYPE=VELOCITY
         >>> <namePrefix>INSTALLER_ROT, 4, 6, 0.0

        @param footpoint: Position of the bolts foot node, usually close to
        the tunnel wall. Default: [0,0,0]

        @param direction: vector (of arbitrary length) pointing away from the
        foot point to the other end of the bolt. Default: [0,0,1]

        @param length: total length of the bolt (approximately)

        @param elementLength: Length of the individual elements of the elastic
        bolt.

        @param namePrefix: name prefix for elset, section, material names

        @param installerNsetName: Create an nset with this name to apply
        temperature to at installation time.

        @param elType: "T3D2" for truss elements, "B31" for beam
        elements

        @Note: The footpoint and direction can come from a
        bae.generatemesh.tunnel_01.StraightTunnel object like that:
         >>> from bae.generatemesh.tunnel_01 import XSectCircBack,StraightTunnel
         >>> tunnelXSect = XSectCircBack(width=5.0, height=5.0,backRadius = 3.0)
         >>> tunnelGeom = StraightTunnel(
         ...     crossSection=tunnelXSect, length=12.0, heading=0.0)
         >>> for y in [1.0, 2.5, 4.0, 5.5]:
         >>>     for s in [-5.0, -2.5, 0.0, 2.5, 5.0]:
         >>>         pos = tunnelGeom.ys2xyzContour(y, s)
         >>>         drn = tunnelGeom.ys2normal(y, s)
         >>>         currentBolt = BoltEmbeddedTrussEqualLen(
         ...             pos, drn, **boltsparameter)
        """

        if installerNsetName is None:
            installerNsetName="%sINSTALLER" % namePrefix

        Model.__init__(self)

        self.boltLength = float(length)
        direction = norm(direction)

        nbElastBoltElems = int(round(self.boltLength/elementLength))

        # create nodes and elems
        firstNode = 1
        firstElem = 1

        # define nodes along the elastic steel bolt
        nodePositions = [i*elementLength for i in range(nbElastBoltElems+1)]
        coords = [vector_plus(footpoint,
                              vector_scale(direction, pos))
                  for pos in nodePositions]
        self.nodesElast = self.insertNodes(coords, firstNode=firstNode)
        firstNode = self.nodesElast[-1] + 1  # for next iteration

        # nodes embedded in rock
        self.nodesEmbedded = list(self.nodesElast)

        # installer nset (for compatibility with other bolt types)
        self.forceNset(installerNsetName).update(self.nodesElast)

        # elastic bolt elements (T3D2 truss or B31 beam)
        elNodes = [[self.nodesElast[i], self.nodesElast[i+1]]
                   for i in range(nbElastBoltElems)]
        self.elemsElast = self.insertElems(
            elNodes, elType, firstElem=firstElem)
        self.forceElset("%sBOLT" % namePrefix).update(self.elemsElast)
        embeddedElset = self.forceElset("%sEMBEDDED_ELEMS" % namePrefix)
        embeddedElset.update(self.elemsElast)

        if elType == "B31":
            nsetRotBC = self.forceNset("%sINSTALLER_ROT" % namePrefix)
            nsetRotBC.update(self.nodesElast)

        self.nodeCollar = self.nodesElast[0]
        self.nodeBoltEnd = self.nodesElast[-1]
        self.elemCollar = self.elemsElast[0]
        self.elemBoltEnd = self.elemsElast[-1]

        self.namePrefix = namePrefix  # just store
        return

    def renumber(self, nodeStart=1, elemStart=1):
        """renumber elements and nodes

        Examples:
          >>> # renumber nodes and elements to start at 500001
          >>> bolt.renumber(nodeStart=500001, elemStart=500001)

        @param nodeStart: New start for node numbers. If None, nodes are not
        renumbered.
        @param elemStart: New start for element numbers. If None, elements are
        not renumbered.

        @Note: For further details look at L{bae.abq_model_02.Model.renumber}.

        @Returns: Tuple of two dictionaries (nodesOldToNew, elemsOldToNew)
        maping old node/element numbers to the new ones.
        """

        (nodesOldToNew, elemsOldToNew) = Model.renumber(
            self, nodeStart, elemStart)

        self.nodeCollar = nodesOldToNew[self.nodeCollar]
        self.nodeBoltEnd = nodesOldToNew[self.nodeBoltEnd]
        self.nodesElast = [nodesOldToNew[node] for node in self.nodesElast]
        self.elemsElast = [elemsOldToNew[elem] for elem in self.elemsElast]
        self.nodesEmbedded = [nodesOldToNew[node]
                              for node in self.nodesEmbedded]
        if hasattr(self, "elemCollar"):
            self.elemCollar = elemsOldToNew[self.elemCollar]
        if hasattr(self, "elemBoltEnd"):
            self.elemBoltEnd = elemsOldToNew[self.elemBoltEnd]
        if hasattr(self, "nodesBoreHoles"):
            self.nodesBoreHoles = [nodesOldToNew[node]
                                   for node in self.nodesBoreHoles]

        return (nodesOldToNew, elemsOldToNew)

    def install(self, footpoint, direction, groupName,
                model, dirVectors, boltDataList, embeddedElsets,
                firstNode=None, firstElem=None):
        """Install bolt in model updating relevant data containers.

        @param footpoint: position of the first node
        @param direction: vector (list of floats). Needn't be unit length.
        @param groupName: Identifier string for this group of bolts.
            (Layer name in the Rhino file)
        @param model: L{bae.abq_model_02.Model}-object to insert the new bolt
        @param dirVectors: dictionary {direction name: direction unit vector}
            of all bolts installed so far. Direction names will be D001, D002
            and so forth. Directions within about 1 deg of each other are
            considered identical and get the same name.
        @param boltDataList: L{BoltDataList}-object to be updated with this
            new bolt
        @param embeddedElsets: Set of elset names, containing all elsets to
            be embedded. Will be updated.
        @param firstNode: renumber nodes before inserting into model
        @param firstElem: renumber elements before inserting into model

        @Note: This will most likely change self: rotate, translate,
            rename elsets and orientations, renumber elements and nodes.
        """

        # determine, sort and store directions
        direction = norm(direction)

        if len(dirVectors)==0:
            # tweak next condition to make this first direction be stored
            oldDir = [0,0,0]
        else:
            # find old direction in dirVector that is closest to new direction
            oldDirName = getShortestDistKey(direction, dirVectors)
            oldDir = dirVectors[oldDirName]

        if length(vector(oldDir, direction)) > 0.02:  # that's about tol = 1 deg
            dirName = "%03d" % (len(dirVectors)+1)
            dirVectors[dirName] = direction
        else:
            dirName = oldDirName

        # rotate and move bolt into position
        self.rotate([0,0,0], trans_vert2dir(direction), newDirectionId=dirName)
        self.translate(footpoint)

        # insert into common model
        self.renumber(nodeStart=firstNode, elemStart=firstElem)
        (nodesNew, elemsNew) = model.insertModel(
            self, orientationTolerance=0.05)

        # find some elset names
        if self.hasDirection:
            elastBoltElset = "%sBOLT_DIR%s" % (self.namePrefix, dirName)
            elastBoltElsetWoDir = "%sBOLT" % self.namePrefix
        else:
            elastBoltElset = "%sBOLT" % self.namePrefix
            elastBoltElsetWoDir = elastBoltElset

        # update all bolts elsets
        elastBoltsElems = set(elemsNew[elem]
                              for elem in self.elset[elastBoltElset])
        model.forceElset(elastBoltElsetWoDir).update(elastBoltsElems)
        model.forceElset(elastBoltElset).update(elastBoltsElems)

        # add to embedded elements elset list
        if self.hasEmbeddedRubber:
            embeddedElsets.add("%sEMBEDDED_ELEMS" % self.namePrefix)
        else:
            embeddedElsets.add(elastBoltElset)

        # store some bolt data (basic types only: float, int, string)
        boltData = Container(
            footpoint = footpoint,
            directionName = dirName,
            groupName = groupName,
            boltElType = self.boltElType,
            elastElems = [elemsNew[elem]
                          for elem in sorted(self.elset[elastBoltElset])],
            elemCollar = elemsNew[self.elemCollar],
            elemBoltEnd = elemsNew[self.elemBoltEnd],
            nodeCollar = nodesNew[self.nodeCollar],
            nodeBoltEnd = nodesNew[self.nodeBoltEnd],
            nodesElast = [nodesNew[node] for node in self.nodesElast],
            # matName = matName,
        )
        # try:
        #     boltData.crossSectArea = boltsparameter["crossSectArea"]
        # except (AttributeError, KeyError):
        #     pass
        try:
            boltData.nodesBoreHoles = [
                nodesNew[node] for node in self.nodesBoreHoles]
        except AttributeError:
            pass
        try:
            boltData.nodesEmbedded = [
                nodesNew[node] for node in self.nodesEmbedded]
        except AttributeError:
            pass
        boltDataList.append(boltData)


###############################################################################

class BoltBaseTrussElems(BoltModel):
    """additional base class for bolts with truss elements"""

    hasDirection = False
    boltElType = "truss"

    def initSteelBoltSection(self, elsetName=None, material="BOLT_STEEL",
                             crossSectArea=1.0):
        """Initialize a beam section definition for the steel bolt.

        @param elsetName: Elset containing the steel bolt elements. If not
           specified it's being guessed from what happens in the constructor.
           This will work if you do this straight away before modifying this
           model any further.
        @param material: A string: XXX in *BEAM SECTION, MATERIAL=XXX
        @param crossSectArea: cross sectional area of the beam: The number on
           the data line of the *SOLID SECTION command.
        """
        if elsetName is None:
            elsetNames = [els for els in self.elset
                          if els.endswith("BOLT")]
            if len(elsetNames)<1:
                raise ValueError("Could not guess the steel bolts elset name,"
                                 " no elset present in the model.")
            if len(elsetNames)>1:
                msg("WARNING: There are multiple elsets that might be the"
                    " steel bolts elset:%s.\n <%s> chosen."
                    % (",".join(elsetNames), elsetName))
            elsetName = elsetNames[0]

        try:
            crossSectArea = float(crossSectArea)
        except (ValueError, TypeError):
            raise ValueError("Invalid crossSectArea %s for elset %s."
                             % (crossSectArea, elsetName))

        self.properties.updateSection(
            type="SOLIDSECTION", elset=elsetName,
            material=material, data=[[crossSectArea, ], ] )

    def rotate(self, refPoint, rotMat, newDirectionId=None):
        """Rotate the mesh according to the specified transformation matrix
        rotMat relative to the given point refPoint.

        @param refPoint: point around which to rotate the mesh
        @param rotMat: transformation matrix. No check is performed but the
           matrix should perform a pure rotation.
        @param newDirectionId: This parameter is being ignored. It's there for
           compatibility to BoltFrictionSleeve.rotate()

        @Note: See L{rotateX} for further arguments and notes on limitations
        and features.
        """
        # rotate nodes
        Mesh.rotate(self, refPoint, rotMat)

    def rotateX(self, refPoint, alpha, newDirectionId=None):
        """Rotate the model around a line parallel to the x axis.

        @param refPoint: point around which to rotate the mesh
        @param alpha: angle in degree, pos == clockwise when looking in
        x-direction.
        @param newDirectionId: This parameter is being ignored. It's there for
        compatibility to BoltFrictionSleeve.rotate()
        """
        Mesh.rotateX(self, refPoint, alpha)

    def rotateY(self, refPoint, alpha, newDirectionId=None):
        """Rotate the mesh around a line parallel to the y axis.

        @param alpha: angle in degree, pos == clockwise when looking in
        y-direction.

        @Note: See L{rotateX} for further arguments and notes on limitations
        and features.
        """
        Mesh.rotateY(self, refPoint, alpha)

    def rotateZ(self, refPoint, alpha, newDirectionId=None):
        """Rotate the mesh around a line parallel to the z axis.

        @param alpha: angle in degree, pos == anticlockwise when seen from
        above (i.e. looking in neg. z-direction).

        @Note: See L{rotateX} for further arguments and notes on limitations
        and features.
        """
        Mesh.rotateZ(self, refPoint, alpha)


###############################################################################

class BoltBaseBeamElems(BoltModel):
    """additional base class for bolts with beam elements"""

    hasDirection = True
    boltElType = "beam"

    def initSteelBoltSection(self, elsetName=None, material="BOLT_STEEL",
                             section="CIRC", sectionPara=[0.0086],
                             sectionOrientation=[0,1,0]):
        """Initialize a beam section definition for the steel bolt.

        @param elsetName: Elset containing the steel bolt elements. If not
           specified it's being guessed from what happens in the constructor.
           This will work if you do this straight away before modifying this
           model any further.
        @param material: A string: XXX in *BEAM SECTION, MATERIAL=XXX
        @param section: A string: XXX in *BEAM SECTION, SECTION=XXX
        @param sectionPara: A list of floats defining the first data line of
           the *BEAM SECTION command.
            - In case of section="CIRC" it's the radius of the circular
              section: sectionPara=[myradius]
            - For pipe it's the outer radius and the pipe thickness:
              sectionPara=[myradius, mythickness]
        @param sectionOrientation: a list of three floats for the second data
           line of the *BEAM SECTION command. Defines the first beam section
           axis. Must be roughly perpendicular to the beam orientation as given
           by the direction argument to the constructor.
        """
        if elsetName is None:
            elsetNames = [els for els in self.elset
                          if els.endswith("BOLT")]
            if len(elsetNames)<1:
                raise ValueError("Could not guess the steel bolts elset name,"
                                 " no elset present in the model.")
            if len(elsetNames)>1:
                msg("WARNING: There are multiple elsets that might be the"
                    " steel bolts elset:%s.\n <%s> chosen."
                    % (",".join(elsetNames), elsetName))
            elsetName = elsetNames[0]

        self.properties.updateSection(
            type="BEAMSECTION", elset=elsetName, section=section,
            # very important: make a copy of the data lists, otherwise
            material=material, data=[list(sectionPara), list(sectionOrientation)])

    def _rotateRenameElset(self, elsetName, newDirectionId):
        """Rename elset associated with a direction.
        In the process of rotating the model elsets associated with a direction
        or orientation shall be renamed in order to distinguish them from
        elsets with other orientations.

        E.g. create a bolt and rotate it, then create a similar bolt and rotate
        it differently. Then merge those two models. In this case we need to
        distinguish the elsets from the two models.

        Direction identifiers are identified by their "DIR"-prefix and separated
        from the rest by underscore "_". As in "BOLT_DIRNORTH", "ORI_DIR05_BLT".

        This function is only called internally from _rotateRename and this in
        turn from the various rotate functions.

        @param elsetName: old name of the elset
        @param newDirectionId: new direction part of the elset name excluding
           "DIR". I.e. "06" for "BOLT_DIR06".
        """
        if newDirectionId is None:
            return

        # determine new elset name
        if newDirectionId == "autoRename":
            newElsetName = incrementName(
                elsetName,
                numberpattern=re.compile(r"_DIR(\d+)"),
                defaultnumber="_DIR001")
        else:
            (newElsetName, nbSub) = re.subn(
                r"(.*)_(DIR.*)((?:_.*)?)", "\\1_DIR%s\\3" % newDirectionId,
                elsetName)
            if nbSub==0:
                newElsetName="%s_DIR%s" % (elsetName, newDirectionId)

        if newElsetName!=elsetName:
            # make sure this elset name is not already defined
            while ((newElsetName in self.elset)
                   or (newElsetName in self.properties)):
                newElsetName = incrementName(newElsetName)

            # rename elset
            self.properties[newElsetName] = self.properties[elsetName]
            del self.properties[elsetName]
            self.elset[newElsetName] = self.elset[elsetName]
            del self.elset[elsetName]

    def _rotateRename(self, refPoint, newDirectionId):
        """Rename orientations (if present) and elsets used in section
        definitions that are affected by rotation.

        This function is only called internally from the various rotate
        functions.

        @param refPoint: currently not used.
        @param newDirectionId: New direction part of the name of the
        orientation and elset excluding the "DIR" prefix.
        """

        if newDirectionId is None:
            return

        newOriName = None
        if len(self.orientation) and hasattr(self, "orientationName"):

            # determine new name for orientation
            if len(self.orientation) > 1:
                msg("WARNING: There are more than just one orientation"
                    " definitions in the model. Only <%s> will be renamed when"
                    " rotating the model." % self.orientationName)
                if newDirectionId == "autoRename":
                    newOriName = incrementName(
                        self.orientationName,
                        numberpattern=re.compile(r"_DIR(\d+)"),
                        defaultnumber="_DIR01")
            else:
                (newOriName, nbSubs) = re.subn(
                    r"ORIENTATION_.*", "ORIENTATION_DIR%s" % newDirectionId,
                    self.orientationName)
                if nbSubs==0:
                    newOriName = "%s_DIR%s" % (
                        self.orientationName, newDirectionId)

            # rename orientation-dict key
            if newOriName!=self.orientationName:
                self.orientation.renameItem(self.orientationName, newOriName)

        # rename in connector section
        # loop over a copy of the key, value pairs because the dict will be
        # modified!
        for elset, prop in self.properties.items():
            if prop["TYPE"]=="CONNECTORSECTION":
                if len(self.orientation):
                    try:
                        # thisOri is the second data line, i.e. a list of items
                        # the first (and only) item -- thisOri[0] -- is expected
                        # to be the name of the orientation
                        thisOri = prop["DATA"][1]
                        thisOri[0]  # just a check it's there
                    except (KeyError, IndexError):
                        # no orientation defined in this section
                        continue

                    if hasattr(self, "orientationName"):
                        if thisOri[0] == self.orientationName:
                            # rename the orientation
                            thisOri[0] = newOriName
                        else:
                            msg("WARNING: Expected orientation <%s> but found"
                                " <%s> in connector section definition for"
                                " elset %s. Will not be renamed. The result is"
                                " unlikely to be correct."
                                % (self.orientationName, thisOri[0], elset))
                            continue

                # rename the corresponding elset as well
                self._rotateRenameElset(elset, newDirectionId)

            elif prop["TYPE"]=="BEAMSECTION":
                if prop['SECTION'] in ("BOX", "CIRC", "HEX", "I", "L", "PIPE",
                                       "RECT", "TRAPEZOID", "ARBITRARY"):

                    # rename the corresponding elset
                    self._rotateRenameElset(elset, newDirectionId)

        # store the new orientation name
        if newOriName:
            self.orientationName = newOriName

    def rotate(self, refPoint, rotMat, newDirectionId="autoRename"):
        """Rotate the mesh according to the specified transformation matrix
        rotMat relative to the given point refPoint.

        @param refPoint: point around which to rotate the mesh
        @param rotMat: transformation matrix. No check is performed but the
           matrix should perform a pure rotation.

        @Note: See L{rotateX} for further arguments and notes on limitations
        and features.
        """
        Model._rotateXXX(
            self, refPoint, rotMat, Mesh.rotate, Model._vector_rot_mat)
        self._rotateRename(refPoint, newDirectionId)

    def rotateX(self, refPoint, alpha, newDirectionId="autoRename"):
        """Rotate the model around a line parallel to the x axis.

        @param refPoint: point around which to rotate the mesh
        @param alpha: angle in degree, pos == clockwise when looking in
        x-direction.
        @param newDirectionId: Elset and orientation names will be changed.
        This parameter states the new name of the direction part of the
        orientation and elset excluding the "DIR" prefix.

        The result is as if the constructor would have been supplied with
        this string as directionID argument.

        Can be "autoRename" (the default), in this case a default will be used.

        Can be None, in this case elsets and orientations will not be renamed.

        Be careful with those two options: If you copy one template several
        times and then rotate each of the copies only once to its final
        position and then assemble this then you will have the same name for
        different orientations in the final model. So: this will not work!

        @Note: This works reliably only if no new items have been added
        to the model since initialization or if the following conditions are
        met:
         - no new orientations
        """
        Model._rotateXXX(self, refPoint, alpha, Mesh.rotateX, vector_rot_x)
        self._rotateRename(refPoint, newDirectionId)

    def rotateY(self, refPoint, alpha, newDirectionId="autoRename"):
        """Rotate the mesh around a line parallel to the y axis.

        @param alpha: angle in degree, pos == clockwise when looking in
        y-direction.

        @Note: See L{rotateX} for further arguments and notes on limitations
        and features.
        """
        Model._rotateXXX(self, refPoint, alpha, Mesh.rotateY, vector_rot_y)
        self._rotateRename(refPoint, newDirectionId)

    def rotateZ(self, refPoint, alpha, newDirectionId="autoRename"):
        """Rotate the mesh around a line parallel to the z axis.

        @param alpha: angle in degree, pos == anticlockwise when seen from
        above (i.e. looking in neg. z-direction).

        @Note: See L{rotateX} for further arguments and notes on limitations
        and features.
        """
        Model._rotateXXX(self, refPoint, alpha, Mesh.rotateZ, vector_rot_z)
        self._rotateRename(refPoint, newDirectionId)


###############################################################################

class BoltFrictionSleeve(BoltBaseBeamElems):
    r"""
    The bolt can be switched on by switching the temperature of a certain node
    set between certain values. (See the __init__() function for further
    descriptions of the arguments to the constructor mentioned below.)

    The connection between the bolt and the rock is assumed to consist of
    two regions with different pull out resistence / friction. The far end
    region is referred to as sleeve, the rest as bolt (assuming that the
    latter represents more or less the friction of the naked bolt in the
    bore hole). sleeveLength is the length of the former region and
    length-sleeveLength the length of the latter.

    Each of the two regions generate a certain integral friction force
    specified by sleeve_frictionForce and bolt_frictionForce respectively.
    Those forces are distributed onto the individual nodes connected to the
    rock such that the sum of sleeve_frictionForce and bolt_frictionForce
    is equal to the total force all connecting nodes transmit to the rock.
    The sleeve_frictionForce is distributed onto the nodes at the far end
    region of length sleeveLength and bolt_frictionForce is distributed
    onto the nodes at the near end region. In between there is usually one
    node that gets contributions from both.

    The friction forces can depend on the pull out distance of each of the
    nodes. This can be used to model degradation of the connection.

    Don't forget to install the bolts by setting the temperature to the
    "On"-value specified in tempControlOffOnValues. The name of the nset
    for the *TEMPERATURE command is specified as a parameter.

    Note that some nodes need a boundary condition to constrain their
    rotation. In case of a step with newly specified boundary conditions
    (*BOUNDARY, OP=NEW) you have to also respecify those boundary conditions.
    For this purpose you can create an input file to include into the *STEP
    definition by means of self.writeBoundaryConditions(). See also
    self.writeConstraints().

    There will be an NSET "<namePrefix>INSTALLER_ROT" for all nodes that need
    to be fixed in their rotation. Note that this might have to be stated for
    each step separately.
     >>> *BOUNDARY, TYPE=VELOCITY
     >>> <namePrefix>INSTALLER_ROT, 4, 6, 0.0

    There will be an ELSET "<namePrefix>BOLT_DIRXXX" containing all bolt
    elements (the actual steel bolt). SECTION properties and material needs
    to be specified separately and will not be predefined. XXX is replaced by
    the directionId argument of the constructor.

    There will be an ELSET "<namePrefix>EMBEDDED_ELEMS". Those T3D2 elements
    need to be embedded in the bulk rock elements:
     >>> *EMBEDDED ELEMENT, HOST ELSET= ... , ABSOLUTE EXTERIOR TOLERANCE= ...
     >>> <namePrefix>EMBEDDED_ELEMS"

    There will be an ELSET "<namePrefix>BOREHOLE" identifying all bore hole
    connector elements.

    There will be an ELSET "<namePrefix>BOREHOLE_COLLAR" and an ELSET
    "<namePrefix>_FAREND" identifying the first and last bore hole connector
    elements.

    There will be ELSETs "<namePrefix>BOREHOLE_X_DIRZZZ". ZZZ is replaced by
    the directionId argument of the constructor, X is replaced by A, B, C, ...
    identifying sections of different friction force along the bolt starting
    from the collar. Those elsets are required for internal purposes. Don't
    bother.

    There will be an ELSET "<namePrefix>INSTALLER" containing all installer
    connector elements.

    There will be an NSET as specified by the constructor argument
    installerNsetName. It contains all nodes on the "installer" connector
    elements. The bolts must be switched on or off by specifying corresponding
    temperature values to those nodes.
     >>> *AMPLITUDE, NAME=INSTALL_TIMING
     >>> ...
     >>> *STEP
     >>> ...
     >>> *TEMPERATURE, AMPLITUDE=INSTALL_TIMING
     >>> INSTALLER, 1.0

    @ivar boltLength: total length in meters = "length" constructor-arg
    @ivar boltFrictionTab: list of (relative displacement, total friction
        force)-tuples specifying the (possibly degrading) friction force
        of the non sleeve part of the bolt (the inner section of the bolt)
    @ivar sleeveLength: Length of the plastic sleeve at the top/far end
        of the Yieldlok bolt. In this region the friction force specified by
        sleeveFrictionTab applies. May be zero.
    @ivar sleeveFrictionTab: A list of (plastic displacement, friction
        force) tuples. This friction force will be distributed over the nodes
        within the region specified by sleeveLength on the far end of the
        bolt. Does not need to be sensible if sleeveLength is zero.
    @ivar tempControlOffOnValues: An (off value, on value)-tuple
        specifying the temperatur values that indicate the state of the bolts.
        (off means: not yet installed, on means: installed)
    @ivar orientationName: The key of the supposedly only orientation in
        self.orientation.
        Default: "%sORIENTATION_DIR%s"%(namePrefix, directionId)
        Will usually be changed if the model is rotated

    @ivar elemCollar: element number of the connector element at the collar
    @ivar elemBoltEnd: element number of the connector element at the far end
    @ivar nodeCollar: node (number) fixed to the rock at the collar
    @ivar nodeBoltEnd: bore hole node (i.e. fixed to the rock only after
        installation) at the far end
    @ivar nodesElast: nodes along the elastic steel bolt
    @ivar elemsElast: element numbers of the elastic steel bolt
    @ivar nodesBoreHoles: Bore hole nodes. They will get fixed when the
        installer elements lock. There is no such node on the collar (first
        station), here the bore hole is always fixed to the rock.
    @ivar nodesEmbedded: nodes embedded in rock

    @Note: The following boundary conditions are not implemented in this model,
    they must be specified in the input file independently from this model
    definition. Possibly has to be defined for each step separately:
     >>> *BOUNDARY, TYPE=VELOCITY
     >>> <namePrefix>INSTALLER_ROT, 4,6, 0.0

    @Note: *EMBEDDED ELEMENTS card not implemented in this model, must be
    specified in the input file independently from this model definition.
    Only lasts for two consecutive steps. Might need an import analysis. In
    this case has to be specified again.
     >>> *EMBEDDED ELEMENT, HOST ELSET=..., ABSOLUTE EXTERIOR TOLERANCE=...
     >>> <namePrefix>EMBEDDED_ELEMS

    @Todo: Remove self.orentationName. In __init__ make it a local variable.
    Later (for rotations) just treat all existing orientations accordingly...
    """

    hasEmbeddedRubber = True

    def __init__(self, *args, **kwargs):
        """
        @kwarg footpoint: Position of the bolts foot node, usually close to
        the tunnel wall. Default: [0,0,0]

        @kwarg direction: vector (of arbitrary length) pointing away from the
        foot point to the other end of the bolt. Default: [0,0,1]

        @kwarg directionId: name for the direction to identify orientations
        and elsets associated with this direction. This string will be
        prepended by "DIR" to generate names like "BOLT_DIRNORTH",
        "ORI_DIR05_BOLT" or so. Defaults to "00".
        It's suggested to skip this argument if you later use self.install().

        @kwarg length: total length of the bolt assembly

        @kwarg installerNsetName: nodes on the "installer" connector elements
        of this bolt will be added to this nset. If an nset of this name
        does not exist yet, it will be created. If not specified it defaults to
        "<namePrefix>INSTALLER". This nset will be initialized
        with temperature 0.0 (off) automatically. The bolts can be switched on
        by specifying temperature 1 to those nodes. The temperature values for
        on and off can be specified, see self.__init__()::

         *AMPLITUDE, NAME=INSTALL_TIMING
         ...
         *STEP
         ...
         *TEMPERATURE, AMPLITUDE=INSTALL_TIMING
         INSTALLER, 1.0

        @kwarg bolt_frictionForce: Either just a float to specify the constant
        friction force or a list of (plastic displacement, friction force)
        tuples. This friction force will be distributed over the innermost
        nodes within the region (length-sleeveLength) of the bolt.
        If both bolt_frictionForce and sleeve_frictionForce are tables the
        plastic displacement values must be identical in both tables!

        @kwarg sleeveLength: Length of the plastic sleeve at the top/far end
        of the Yieldlok bolt. In this region the friction force specified by
        sleeveFrictionTab applies. May be zero.

        @kwarg sleeve_frictionForce: Either just a float to specify the
        constant friction force or a list of (plastic displacement, friction
        force) tuples. This friction force will be distributed over the nodes
        within the region specified by sleeveLength on the far end of the
        bolt. Ignored if sleeveLength is zero.
        If both bolt_frictionForce and sleeve_frictionForce are tables the
        plastic displacement values must be identical in both tables!

        @kwarg elementLengthMax: Either the approximate length of the
        individual elements of the elastic bolt section or None. If None there
        will be exactly one elastic element per bolt.

        @kwarg tempControlOffOnValues: An (off value, on value)-tuple
        specifying the temperatur values that indicate the state of the bolts.
        (off means: not yet installed, on means: installed)

        @kwarg namePrefix: name prefix for elset, section, material names

        @Note: The footpoint and direction can come from a
        bae.generatemesh.tunnel_01.StraightTunnel object like that:
         >>> from bae.generatemesh.tunnel_01 import XSectCircBack,StraightTunnel
         >>> tunnelXSect = XSectCircBack(width=5.0, height=5.0,backRadius = 3.0)
         >>> tunnelGeom = StraightTunnel(
         ...     crossSection=tunnelXSect, length=12.0, heading=0.0)
         >>> assembly = BoltsAssembly(**boltsparameter)
         >>> for y in [1.0, 2.5, 4.0, 5.5]:
         >>>     for s in [-5.0, -2.5, 0.0, 2.5, 5.0]:
         >>>         pos = tunnelGeom.ys2xyzContour(y, s)
         >>>         drn = tunnelGeom.ys2normal(y, s)
         >>>         assembly.addBolt(pos, drn)
        """

        if len(args)==1 and isinstance(args[0], self.__class__):
            # copy constructor
            oldmodel = args[0]
            self.copyConstructor(oldmodel)

            self.boltFrictionTab = copy.deepcopy(oldmodel.boltFrictionTab)
            self.sleeveLength = oldmodel.sleeveLength
            self.sleeveFrictionTab = copy.deepcopy(oldmodel.sleeveFrictionTab)
            self.tempControlOffOnValues = copy.deepcopy(
                oldmodel.tempControlOffOnValues)
            self.orientationName = oldmodel.orientationName
            self.nodesBoreHoles = list(oldmodel.nodesBoreHoles)
        else:
            length = _parseKwArg(self, kwargs, "length")
            footpoint = _parseKwArg(self, kwargs, "footpoint", default=[0,0,0])
            direction = _parseKwArg(self, kwargs, "direction", default=[0,0,1])
            directionId = _parseKwArg(self, kwargs, "directionId", default="00")
            installerNsetName = _parseKwArg(
                self, kwargs, "installerNsetName", default=None)
            bolt_frictionForce = _parseKwArg(
                self, kwargs, "bolt_frictionForce", default=50000.0)
            sleeveLength = _parseKwArg(
                self, kwargs, "sleeveLength", default=0.0)
            sleeve_frictionForce = _parseKwArg(
                self, kwargs, "sleeve_frictionForce", default=None)
            elementLengthMax = _parseKwArg(
                self, kwargs, "elementLengthMax", default=None)
            tempControlOffOnValues = _parseKwArg(
                self, kwargs, "tempControlOffOnValues", default=(0.0, 1.0))
            namePrefix = _parseKwArg(self, kwargs, "namePrefix", default="")

            Model.__init__(self)

            if installerNsetName is None:
                installerNsetName="%sINSTALLER" % namePrefix

            self.boltLength = float(length)

            # bolt
            try:
                float(bolt_frictionForce)
            except TypeError:
                self.boltFrictionTab = bolt_frictionForce
            else:
                self.boltFrictionTab = [[0.0, bolt_frictionForce]]

            # sleeve
            self.sleeveLength = sleeveLength

            try:
                float(sleeve_frictionForce)
            except TypeError:
                self.sleeveFrictionTab = sleeve_frictionForce
            else:
                self.sleeveFrictionTab = [[0.0, sleeve_frictionForce]]

            # other
            if elementLengthMax is None:
                elementLengthMax = self.boltLength

            self.tempControlOffOnValues = tempControlOffOnValues

            direction = norm(direction)

            # create nodes and elems
            firstNode = 1
            firstElem = 1

            nbElastBoltElems = int(ceil(self.boltLength
                                        / elementLengthMax))
            elementLength = self.boltLength / nbElastBoltElems
            nbSigDigs = int(log(nbElastBoltElems,10)+0.5) + 4

            # friction force fractions for the individual nodes
            lastBoreholeFrictFracs = None
            boreholeFrictFracList = list()
            boreholeFrictFracFirstEl = list()
            bareBoltLength = self.boltLength-self.sleeveLength
            for i in range(nbElastBoltElems + 1):
                boreholeFrictFracs = [None, None]

                bareBoltSectStart = max(0.0, (i-0.5)*elementLength)
                bareBoltSectEnd = min(bareBoltLength, (i+0.5)*elementLength)
                boreholeFrictFracs[0] = max(
                    0.0, (bareBoltSectEnd-bareBoltSectStart)/bareBoltLength)

                if self.sleeveLength>0:
                    sleeveSectStart = max(bareBoltLength, (i-0.5)*elementLength)
                    sleeveSectEnd = min(self.boltLength, (i+0.5)*elementLength)
                    boreholeFrictFracs[1] = max(
                        0.0, ((sleeveSectEnd-sleeveSectStart)
                              /self.sleeveLength))
                else:
                    boreholeFrictFracs[1] = 0.0

                # round to significant digits in order to treat (almost) same
                # fractions as same (because different fractions induce
                # different connector sections)
                boreholeFrictFracs = [round(x, nbSigDigs)
                                      for x in boreholeFrictFracs]
                # The first element of each bolt must get a separate section
                # ("A") because there is a connector stop. The ...or i==1...
                # makes sure element one and two are in different sections!
                if boreholeFrictFracs!=lastBoreholeFrictFracs or i==1:
                    boreholeFrictFracFirstEl.append(i)
                    boreholeFrictFracList.append(boreholeFrictFracs)
                    lastBoreholeFrictFracs = boreholeFrictFracs

            # in boreholeFrictFracList store:
            # [("A", [frac-bolt, frac-sleeve]), ("B", [frac-bolt, frac-sleeve]),
            # ... ]
            # frac-bolt is the fraction of bolt friction beeing applied to the
            # particular node, frac-sleeve the corresponding fraction of the
            # sleeve friction
            boreholeFrictFracList = [
                (chr(ord("A") + i), boreholeFrictFracs)
                for i, boreholeFrictFracs in enumerate(boreholeFrictFracList)]

            # define nodes along the elastic steel bolt
            coords = [vector_plus(footpoint,
                                  vector_scale(direction, i*elementLength))
                      for i in range(nbElastBoltElems+1)]
            self.nodesElast = self.insertNodes(coords, firstNode=firstNode)
            firstNode = self.nodesElast[-1] + 1  # for next items

            # bore hole nodes: will get fixed when the installer elements lock
            # No such node on the collar (first station), here the bore hole is
            # always fixed to the rock
            self.nodesBoreHoles = self.insertNodes(
                coords[1:], firstNode=firstNode)
            firstNode = self.nodesBoreHoles[-1] + 1  # for next items
            self.forceNset("%sINSTALLER_ROT"%namePrefix).update(
                self.nodesBoreHoles)

            # nodes embedded in rock
            self.nodesEmbedded = self.insertNodes(coords, firstNode=firstNode)
            firstNode = self.nodesEmbedded[-1] + 1  # for next items
            self.forceNset("%sINSTALLER_ROT"%namePrefix).add(
                self.nodesEmbedded[0])

            # elastic bolt elements (beam)
            elNodes = [[self.nodesElast[i], self.nodesElast[i+1]]
                       for i in range(nbElastBoltElems)]
            self.elemsElast = self.insertElems(
                elNodes, "B31", firstElem=firstElem)
            firstElem = self.elemsElast[-1] + 1  # for next items
            self.forceElset("%sBOLT_DIR%s"%(namePrefix, directionId)
                            ).update(self.elemsElast)

            # embedded truss elements (T3D2 truss)
            elNodes = [[self.nodesEmbedded[i], self.nodesEmbedded[i+1]]
                       for i in range(nbElastBoltElems)]
            embeddedTrussElems = self.insertElems(
                elNodes, "T3D2", firstElem=firstElem)
            firstElem = embeddedTrussElems[-1] + 1  # for next items
            self.forceElset("%sEMBEDDED_ELEMS"%namePrefix).update(
                embeddedTrussElems)

            # bore hole at collar is connected directly to embedded rubber
            elNodes = [ [ self.nodesEmbedded[0], self.nodesElast[0] ] ]

            # bore hole connector elements (all but the first)
            elNodes.extend(
                [n1, n2]
                for n1, n2 in izip(self.nodesBoreHoles, self.nodesElast[1:]))

            boreHoleElems = self.insertElems(
                elNodes, "CONN3D2", firstElem=firstElem)
            firstElem = boreHoleElems[-1] + 1  # for next items

            for i, (secName, dummy) in enumerate(boreholeFrictFracList):
                if i==0:
                    # first section
                    i0=0
                    i1=boreholeFrictFracFirstEl[i+1]
                elif i==len(boreholeFrictFracList)-1:
                    # last section
                    i0=boreholeFrictFracFirstEl[i]
                    i1=len(boreHoleElems)
                else:
                    # intermediate section
                    i0=boreholeFrictFracFirstEl[i]
                    i1=boreholeFrictFracFirstEl[i+1]

                # elset for the connectors of this section
                self.forceElset( "%sBOREHOLE_%s_DIR%s" % ( namePrefix, secName, \
                    directionId)).update(boreHoleElems[i0:i1])

                # special bore hole elset: COLLAR
                if secName=="A":
                    self.forceElset("%sBOREHOLE_COLLAR"%namePrefix).update(
                        boreHoleElems[i0:i1])

            # special bore hole elset: FAREND (data from last loop iteration)
            self.forceElset("%sBOREHOLE_FAREND"%namePrefix).update(
                boreHoleElems[i0:i1])

            # common bore hole elsets
            self.forceElset("%sBOREHOLE"%namePrefix).update(boreHoleElems)

            # installer elements (none at the collar)
            elNodes = [[n1, n2] for n1, n2 in izip(self.nodesBoreHoles,
                                                   self.nodesEmbedded[1:])]
            installerElems = self.insertElems(
                elNodes, "CONN3D2", firstElem=firstElem)
            firstElem = installerElems[-1] + 1  # for next iteration
            self.forceElset("%sINSTALLER"%namePrefix).update(installerElems)
            self.forceNset(installerNsetName).update(self.nodesBoreHoles)
            self.forceNset(installerNsetName).update(self.nodesEmbedded)
            self.forceNset(installerNsetName).add(self.nodesElast[0])

            self.elemCollar = boreHoleElems[0]
            self.elemBoltEnd = boreHoleElems[-1]
            self.nodeCollar = self.nodesEmbedded[0]
            self.nodeBoltEnd = self.nodesBoreHoles[-1]

            # connector orientations (for bore hole connectors)
            if direction[1]>0.7071:
                # direction closer to y-axes than 45 deg
                dirSecond = [1.0, 0.0, 0.0]
            else:
                # otherwise: default 2-axis is y-axis
                dirSecond = [0.0, 1.0, 0.0]

            self.orientationName = "%sORIENTATION_DIR%s"%(namePrefix, directionId)

            self.orientation.updateItem(self.orientationName,
                                        direction, dirSecond)

            # connector lock datalines (for installers)
            offValue = (0.55*self.tempControlOffOnValues[0]
                        + 0.45*self.tempControlOffOnValues[1])
            onValue = (0.45*self.tempControlOffOnValues[0]
                       + 0.55*self.tempControlOffOnValues[1])

            dataLineOff = [ None,None,None,None, -99999.9, 99999.9, offValue ]
            dataLineOn = [ None,None,None,None, 0.0, 0.0, onValue ]
            if offValue<onValue:
                connLockDataLines = [dataLineOff, dataLineOn]
            else:
                connLockDataLines = [dataLineOn, dataLineOff]

            # connector lock datalines for blocking the borehole at the collar
            # before installation
            dataLineUnBlk = [ None,None,None,None, -99999.9, 99999.9, onValue ]
            dataLineBlock = [ None,None,None,None, 0.0, 0.0, offValue ]
            if offValue<onValue:
                connBlockCollarDataLines = [dataLineBlock, dataLineUnBlk]
            else:
                connBlockCollarDataLines = [dataLineUnBlk, dataLineBlock]

            # material and section properties for embedded rubber
            self.material.updateItem("%sEMBEDDED_RUBBER"%namePrefix,{
                "Alpha": 0.5,
                "Density": 1000.0,
                "EMod": 1E6,
                "Nue": 0.3,
                })
            self.properties.updateSection(
                type="SOLIDSECTION", elset="%sEMBEDDED_ELEMS"%namePrefix,
                material="%sEMBEDDED_RUBBER"%namePrefix, data=[[1E-4],] )

            # connector section and behaviour for boreholes
            for secCnt,(secName,frictFrac) in enumerate(boreholeFrictFracList):

                #-- init behavior, elasticity
                connBehavior = {"elasticity": [(1, "rigid"),]}

                #-- prepare friction tab
                # frictTab = [[force1,displacement1],[force2,displacement2],...]
                if frictFrac[0]==0.0:
                    # sleeve only section
                    frictTab = [
                        [f*frictFrac[1], x]
                        for x, f in self.sleeveFrictionTab]
                elif frictFrac[1]==0.0:
                    # bolt only section
                    frictTab = [
                        [f*frictFrac[0], x]
                        for x, f in self.boltFrictionTab]
                else:
                    # mixed section
                    if len(self.sleeveFrictionTab)==1:
                        protoFrictTab = [
                            [x, f, self.sleeveFrictionTab[0][1]]
                            for x, f in self.boltFrictionTab]
                    elif len(self.boltFrictionTab)==1:
                        protoFrictTab = [
                            [x, self.boltFrictionTab[0][1], f]
                            for x, f in self.sleeveFrictionTab]
                    else:
                        protoFrictTab = [
                            [x1, f1, f2]
                            for (x1, f1), (x2, f2) in izip(
                                self.sleeveFrictionTab, self.boltFrictionTab)]
                    frictTab = [
                        [f1*frictFrac[0]+f2*frictFrac[1], x]
                        for x, f1, f2 in protoFrictTab]
                # condense friction tab if constant friction
                if len(frictTab)==1:
                    assert frictTab[0][1]==0.0, (
                        "Single friction force for non-zero displacement!\n%s"
                        % frictTab)
                    frictTab = [[frictTab[0][0]], ]
                # store
                connBehavior["plasticity"] = [(1, frictTab),]

                #-- first section gets a connector stop
                if secCnt==0:
                    connBehavior["stop"] = [( 1, [[None, 0.0], ] ),]
                    connBehavior["lock"] = [
                        ( 1, "ALL", connBlockCollarDataLines ),]

                #-- store behavior data in self
                self.connectorBehavior.updateItem(
                    "%sBEHAV_BOREHOLE_%s"
                    % (namePrefix, secName), connBehavior )

                #-- store connector sections for boreholes
                self.properties.updateSection(
                    type="CONNECTORSECTION",
                    elset="%sBOREHOLE_%s_DIR%s" % (namePrefix,secName,directionId),
                    behavior="%sBEHAV_BOREHOLE_%s" % (namePrefix, secName),
                    data=[["SLOT"], [self.orientationName]] )

            # connector section and behaviour for installer
            self.connectorBehavior.updateItem(
                "%sbehav_installer" % namePrefix, {
                    "lock": [ (comp, "ALL", connLockDataLines)
                              for comp in [1,2,3] ],
                    } )
            self.properties.updateSection(
                type="CONNECTORSECTION", elset="%sINSTALLER" % namePrefix,
                behavior="%sbehav_installer" % namePrefix,
                data=[["CARTESIAN"], ] )

            # initial conditions for installer nsets
            newInitCond = self.initCond.addNewList("TEMPERATURE")
            newInitCond.append(
                installerNsetName, self.tempControlOffOnValues[0])

        return


###############################################################################

class BoltEmbeddedTruss(BoltBaseTrussElems):
    r"""
    Simple bolt model made of truss elements. All nodes are embedded in rock.
    The first element of each bolt has a specified length (elementLengthMax)
    and at the far end a section of a specified length (sleeveLength) has as
    many nodes as fit into it given the specified elements maximum length of
    elementLengthMax. In between there is just one element.

    This setup is especially suitable for cable bolts that are grouted at the
    far end over a length of sleeveLength.

    There will be an ELSET "<namePrefix>BOLT" containing all bolt elements.
    SECTION properties should be defined using L{initSteelBoltSection}.
    Material data needs to be specified separately.

    There will be an ELSET "<namePrefix>EMBEDDED_ELEMS". Those T3D2 elements
    need to be embedded in the bulk rock elements:
     >>> *EMBEDDED ELEMENT, HOST ELSET= ... , ABSOLUTE EXTERIOR TOLERANCE= ...
     >>> <namePrefix>EMBEDDED_ELEMS"

    @ivar boltLength: total length in meters = "length" constructor-arg
    @ivar sleeveLength: Length of the embedded region at the far end of the
        bolt.
    @ivar nodeBoltEnd: bore hole node (i.e. fixed to the rock only after
        installation) at the far end
    @ivar nodesElast: nodes along the elastic steel bolt
    @ivar elemsElast: element numbers of the elastic steel bolt
    @ivar nodesEmbedded: nodes embedded in rock

    @Note: *EMBEDDED ELEMENTS card not implemented in this model. It must be
    specified in the input file independently from this model definition.
    Only lasts for two consecutive steps. Might need an import analysis. In
    this case has to be specified again.
     >>> *EMBEDDED ELEMENT, HOST ELSET=..., ABSOLUTE EXTERIOR TOLERANCE=...
     >>> <namePrefix>EMBEDDED_ELEMS
    """

    hasEmbeddedRubber = False

    def __init__(self, *args, **kwargs):
        """
        @kwarg footpoint: Position of the bolts foot node, usually close to
        the tunnel wall. Default: [0,0,0]

        @kwarg direction: vector (of arbitrary length) pointing away from the
        foot point to the other end of the bolt. Default: [0,0,1]

        @kwarg length: total length of the bolt assembly

        @kwarg sleeveLength: Length of the embedded region at the far end of
        the bolt.

        @kwarg elementLengthMax: Either the approximate length of the
        individual elements of the elastic bolt section or None. If None there
        will be exactly three elements per bolt.

        @kwarg namePrefix: name prefix for elset, section, material names

        @Note: The footpoint and direction can come from a
        bae.generatemesh.tunnel_01.StraightTunnel object like that:
         >>> from bae.generatemesh.tunnel_01 import XSectCircBack,StraightTunnel
         >>> tunnelXSect = XSectCircBack(width=5.0, height=5.0,backRadius = 3.0)
         >>> tunnelGeom = StraightTunnel(
         ...     crossSection=tunnelXSect, length=12.0, heading=0.0)
         >>> for y in [1.0, 2.5, 4.0, 5.5]:
         >>>     for s in [-5.0, -2.5, 0.0, 2.5, 5.0]:
         >>>         pos = tunnelGeom.ys2xyzContour(y, s)
         >>>         drn = tunnelGeom.ys2normal(y, s)
         >>>         currentBolt = BoltEmbeddedTruss(pos, drn, **boltsparameter)
        """

        if len(args)==1 and isinstance(args[0], self.__class__):
            oldmodel = args[0]
            self.copyConstructor(oldmodel)
            self.sleeveLength = oldmodel.sleeveLength
            self.nbOfSleeveElems = oldmodel.nbOfSleeveElems
        else:
            length = _parseKwArg(self, kwargs, "length")
            footpoint = _parseKwArg(self, kwargs, "footpoint", default=[0,0,0])
            direction = _parseKwArg(self, kwargs, "direction", default=[0,0,1])
            sleeveLength = _parseKwArg(
                self, kwargs, "sleeveLength", default=0.0)
            elementLengthMax = _parseKwArg(
                self, kwargs, "elementLengthMax", default=0.2)
            namePrefix = _parseKwArg(self, kwargs, "namePrefix", default="")

            Model.__init__(self)

            self.boltLength = float(length)
            self.sleeveLength = float(sleeveLength)
            direction = norm(direction)

            self.nbOfSleeveElems = int(ceil(
                self.sleeveLength/elementLengthMax))
            nbElastBoltElems = self.nbOfSleeveElems+2

            # create nodes and elems
            firstNode = 1
            firstElem = 1

            # define nodes along the elastic steel bolt
            nodePositions = (
                [0.0, elementLengthMax] +
                [ (self.boltLength
                   + (float(i)/self.nbOfSleeveElems-1.0)*self.sleeveLength )
                  for i in range(self.nbOfSleeveElems+1) ] )
            coords = [vector_plus(footpoint,
                                  vector_scale(direction, pos))
                      for pos in nodePositions]
            self.nodesElast = self.insertNodes(coords, firstNode=firstNode)
            firstNode = self.nodesElast[-1] + 1  # for next iteration

            # nodes embedded in rock
            self.nodesEmbedded = sorted(set(self.nodesElast[0:2]).union(
                self.nodesElast[-self.nbOfSleeveElems-1:]))
            # In our case this is a complicated expression for
            # self.nodesEmbedded = list(self.nodesElast)
            # It's done this way for historical reasons. And maybe some day
            # we want more than one node in between the collar element and
            # the elements of the sleeve at the far end which would not be
            # embedded. Then this formular still works.

            # elastic bolt elements (T3D2 truss)
            elNodes = [[self.nodesElast[i], self.nodesElast[i+1]]
                       for i in range(nbElastBoltElems)]
            self.elemsElast = self.insertElems(
                elNodes, "T3D2", firstElem=firstElem)
            firstElem = self.elemsElast[-1] + 1  # for next iteration
            self.forceElset("%sBOLT"%namePrefix).update(self.elemsElast)
            embeddedElset = self.forceElset("%sEMBEDDED_ELEMS"%namePrefix)
            embeddedElset.update(self.elemsElast)

            self.nodeCollar = self.nodesElast[0]
            self.nodeBoltEnd = self.nodesElast[-1]
        return


#######################################################################

class BoltEmbeddedTrussEqualLen(BoltBaseTrussElems):
    r"""
    Simple bolt model for fully grouted bolts made of truss elements. All nodes
    are embedded in rock. All elements are of the specified length. The actual
    bolt length will be adjusted.

    There will be an ELSET "<namePrefix>BOLT" containing all bolt elements.
    SECTION properties should be defined using L{initSteelBoltSection}.
    Material data needs to be specified separately.

    There will be an ELSET "<namePrefix>EMBEDDED_ELEMS". Those T3D2 elements
    need to be embedded in the bulk rock elements:
     >>> *EMBEDDED ELEMENT, HOST ELSET= ... , ABSOLUTE EXTERIOR TOLERANCE= ...
     >>> <namePrefix>EMBEDDED_ELEMS"

    @ivar boltLength: total length in meters = "length" constructor-arg
    @ivar nodeCollar: node (number) at the collar
    @ivar nodeBoltEnd: node at the far end
    @ivar nodesElast: nodes along the elastic steel bolt
    @ivar elemsElast: element numbers of the elastic steel bolt
    @ivar nodesEmbedded: nodes embedded in rock
    @ivar elemCollar: element number of the connector element at the collar
    @ivar elemBoltEnd: element number of the connector element at the far end

    @Note: *EMBEDDED ELEMENTS card not implemented in this model, must be
    specified in the input file independently from this model definition.
    Only lasts for two consecutive steps. Might need an import analysis. In
    this case has to be specified again.
     >>> *EMBEDDED ELEMENT, HOST ELSET=..., ABSOLUTE EXTERIOR TOLERANCE=...
     >>> <namePrefix>EMBEDDED_ELEMS
    """

    hasEmbeddedRubber = False

    def __init__(self, *args, **kwargs):
        """
        @kwarg footpoint: Position of the bolts foot node, usually close to
        the tunnel wall. Default: [0,0,0]

        @kwarg direction: vector (of arbitrary length) pointing away from the
        foot point to the other end of the bolt. Default: [0,0,1]

        @kwarg length: total length of the bolt (approximately)

        @kwarg elementLength: Length of the individual elements of the elastic
        bolt.

        @kwarg namePrefix: name prefix for elset, section, material names

        @kwarg installerNsetName: just for compatibility with other bolt types:
        Create an nset with this name to apply temperature to at installation
        time.

        @Note: The footpoint and direction can come from a
        bae.generatemesh.tunnel_01.StraightTunnel object like that:
         >>> from bae.generatemesh.tunnel_01 import XSectCircBack,StraightTunnel
         >>> tunnelXSect = XSectCircBack(width=5.0, height=5.0,backRadius = 3.0)
         >>> tunnelGeom = StraightTunnel(
         ...     crossSection=tunnelXSect, length=12.0, heading=0.0)
         >>> for y in [1.0, 2.5, 4.0, 5.5]:
         >>>     for s in [-5.0, -2.5, 0.0, 2.5, 5.0]:
         >>>         pos = tunnelGeom.ys2xyzContour(y, s)
         >>>         drn = tunnelGeom.ys2normal(y, s)
         >>>         currentBolt = BoltEmbeddedTrussEqualLen(
         ...             pos, drn, **boltsparameter)
        """

        if len(args)==1 and isinstance(args[0], self.__class__):
            self.copyConstructor(oldmodel=args[0])
        else:
            kwargs["elType"] = "T3D2"
            self.initMeshEquallen(**kwargs)


#######################################################################

class BoltEmbeddedBeamEqualLen(BoltBaseBeamElems):
    r"""
    Simple bolt model for fully grouted bolts made of beam elements. All nodes
    are embedded in rock. All elements are of the specified length. The actual
    bolt length will be adjusted.

    There will be an ELSET "<namePrefix>BOLT" containing all bolt elements.
    SECTION properties should be defined using L{initSteelBoltSection}.
    Material data needs to be specified separately.

    There will be an ELSET "<namePrefix>EMBEDDED_ELEMS". Those B31 elements
    need to be embedded in the bulk rock elements:
     >>> *EMBEDDED ELEMENT, HOST ELSET= ... , ABSOLUTE EXTERIOR TOLERANCE= ...
     >>> <namePrefix>EMBEDDED_ELEMS"

    There will be an NSET "<namePrefix>INSTALLER_ROT" for all nodes that you
    might want to fix in their rotation. Note that this might have to be stated
    for each step separately.
     >>> *BOUNDARY, TYPE=VELOCITY
     >>> <namePrefix>INSTALLER_ROT, 4, 6, 0.0

    @ivar boltLength: total length in meters = "length" constructor-arg
    @ivar nodeCollar: node (number) at the collar
    @ivar nodeBoltEnd: node at the far end
    @ivar nodesElast: nodes along the elastic steel bolt
    @ivar elemsElast: element numbers of the elastic steel bolt
    @ivar nodesEmbedded: nodes embedded in rock
    @ivar elemCollar: element number of the connector element at the collar
    @ivar elemBoltEnd: element number of the connector element at the far end

    @Note: *EMBEDDED ELEMENTS card not implemented in this model, must be
    specified in the input file independently from this model definition.
    Only lasts for two consecutive steps. Might need an import analysis. In
    this case has to be specified again.
     >>> *EMBEDDED ELEMENT, HOST ELSET=..., ABSOLUTE EXTERIOR TOLERANCE=...
     >>> <namePrefix>EMBEDDED_ELEMS
    """

    hasEmbeddedRubber = False

    def __init__(self, *args, **kwargs):
        """
        @kwarg footpoint: Position of the bolts foot node, usually close to
        the tunnel wall. Default: [0,0,0]

        @kwarg direction: vector (of arbitrary length) pointing away from the
        foot point to the other end of the bolt. Default: [0,0,1]

        @kwarg length: total length of the bolt (approximately)

        @kwarg elementLength: Length of the individual elements of the elastic
        bolt.

        @kwarg namePrefix: name prefix for elset, section, material names

        @kwarg installerNsetName: just for compatibility with other bolt types:
        Create an nset with this name to apply temperature to at installation
        time.

        @kwarg embedRotConstrained: If True then embed the nodes using
        L{embedNodesRotConstrained}. Default: False

        @Note: The footpoint and direction can come from a
        bae.generatemesh.tunnel_01.StraightTunnel object like that:
         >>> from bae.generatemesh.tunnel_01 import XSectCircBack,StraightTunnel
         >>> tunnelXSect = XSectCircBack(width=5.0, height=5.0,backRadius = 3.0)
         >>> tunnelGeom = StraightTunnel(
         ...     crossSection=tunnelXSect, length=12.0, heading=0.0)
         >>> for y in [1.0, 2.5, 4.0, 5.5]:
         >>>     for s in [-5.0, -2.5, 0.0, 2.5, 5.0]:
         >>>         pos = tunnelGeom.ys2xyzContour(y, s)
         >>>         drn = tunnelGeom.ys2normal(y, s)
         >>>         currentBolt = BoltEmbeddedTrussEqualLen(
         ...             pos, drn, **boltsparameter)
        """

        if len(args)==1 and isinstance(args[0], self.__class__):
            self.copyConstructor(oldmodel=args[0])
        else:

            try:
                embedRotConstrained = kwargs["embedRotConstrained"]
            except KeyError:
                embedRotConstrained = False
            else:
                del kwargs["embedRotConstrained"]
                self.hasEmbeddedRubber = True

            kwargs["elType"] = "B31"
            self.initMeshEquallen(**kwargs)

            if embedRotConstrained:
                # embed nodes rotconstrained
                embeddedRubber = embedNodesRotConstrained(
                    self, self.nodesEmbedded)
                embeddedElset = self.forceElset(
                    "%sEMBEDDED_ELEMS" % kwargs["namePrefix"])
                embeddedElset.update(embeddedRubber)


#######################################################################

def addNodeNormalsToMesh(model):
    """Add node normals to nodes attached to a GS-mesh such that they
    point vertical to the mesh surface (tunnel surface). This will be the
    beam section 2-axis.
    """

    msg("Calculating node normals: Retrieving attached line elements.")
    nodeDirections = defaultdict(list)
    for elem in sorted(model.elNodes):
        elNodes = model.elNodes[elem]
        if len(elNodes) != 2:
            continue
        dir_ = norm(vector(model.nodeCoords[elNodes[0]],
                           model.nodeCoords[elNodes[1]]))
        nodeDirections[elNodes[0]].append(dir_)
        nodeDirections[elNodes[1]].append(vector_scale(dir_, -1))

    msg("Calculating node normals: Condensing tangent directions and"
        " storing normals.")
    for node, directions in nodeDirections.iteritems():
        if len(directions)<2:
            msg("ERROR: Less than two attached lines on node %d."
                " No normal here." % node)
            continue
        if len(directions)>2:
            # more than two directions need to be condensed to only two

            # calculate cos of angles between any two directions
            angleIds = sorted(
                ((dot(dir1, dir2), [i,j])
                 for (i, dir1), (j, dir2)
                 in combinations(enumerate(directions), 2)),
                key=lambda x:abs(x[0]), reverse=True)

            # condense to two directions...
            # while condensing ids change. idTranslate gives the new (current)
            # index for any original one. angleIds-values contain original ids.
            idTranslate = dict( (i,i) for i in range(len(angleIds)) )
            # groups of vectors referring to the same direction (to be averaged)
            directionLists = [ [x] for x in directions ]
            nbJoins = len(directions)-2
            for angle, idTuple in angleIds[:nbJoins]:
                idTuple.sort()
                if angle < 0:
                    # invert directions to be merged
                    for dir_ in directionLists[idTranslate[idTuple[1]]]:
                        vector_modif_scale(dir_, -1)
                directionLists[idTranslate[idTuple[0]]].extend(
                    directionLists[idTranslate[idTuple[1]]])
                del directionLists[idTranslate[idTuple[1]]]
                # update idTranslate
                newTranslate = dict(
                    (i, int(j>=idTranslate[idTuple[1]])*(j-1) or j)
                    for i, j in idTranslate.iteritems())
                idTranslate = newTranslate

            directions[0] = norm(vector_sum(*directionLists[idTranslate[0]]))
            directions[1] = norm(vector_sum(*directionLists[idTranslate[1]]))

        # storing normal
        normal = cross(directions[0], directions[1])
        model.nodeCoords[node][4:6] = normal

    return


def embedNodesRotConstrained(
        model, nodes, offset1=[0.01,0,0], offset2=[0,0.01,0]):
    """Embed the listed nodes including their rotational degrees of freedom
    using connector elements and additional nodes attached to (very light and
    soft) truss elements.

    The given nodes have to be already defined in the model.

    For each listed node there will be two additional nodes offset from the
    listed nodes by vectors given by parameters offset1 and offset2 and
    connected by connectors. In order to be able to embed those new nodes
    they'll be attached to a truss element. Those truss elements are assigned
    to the material "EMBEDDED_RUBBER" that will automatically be added to the
    model if it does not exist yet.

    @param model: L{bae.abq_model_02.Model} object
    @param nodes: iterable of node numbers. Those nodes must already be defined
      in model.
    @param offset1: offset vector from each listed node to its first additional
      node
    @param offset2: offset vector from each listed node to its second additional
      node

    @Note: The given nodes attached to the structural elements have to be
    embedded as well. This is not accomplished by this function.

    @Note: This embedding technique leeads to large implicitly solved
    constrains. Abaqus refuses to start if there are too many. And presumably
    it's very costly anyway...
    """

    msg("embedNodesRotConstrained: inserting new nodes", debugLevel=1)
    newnodes = model.insertNodes(
        vector_plus(model.nodeCoords[node0],x)
        for node0 in nodes
        for x in (offset1, offset2))

    msg("embedNodesRotConstrained: adding connector elements", debugLevel=1)
    iterNN = iter(newnodes)
    newConnectors = model.insertElems(
        ([node0, iterNN.next()] for node0 in nodes for i in range(2)),
        elType="CONN3D2")
    model.forceElset("CONNECTOR_EMBED_RC1").update(newConnectors[0::2])
    model.forceElset("CONNECTOR_EMBED_RC2").update(newConnectors[1::2])

    msg("embedNodesRotConstrained: adding extra embedded elements",
        debugLevel=1)
    iterNN = iter(newnodes)
    rubberElems = model.insertElems(
        ([node1, node2] for node1, node2 in izip(iterNN, iterNN)),
        elType="T3D2")
    elsetRubber = model.forceElset("EMBEDDED_RUBBER_RC")
    elsetRubber.update(rubberElems)

    msg("embedNodesRotConstrained: adding properties", debugLevel=1)
    model.properties.updateSection(
        type="CONNECTORSECTION",
        elset="CONNECTOR_EMBED_RC1",
        data=[["SLOT",],["ORI_EMBED_RC1",]],)
    model.properties.updateSection(
        type="CONNECTORSECTION",
        elset="CONNECTOR_EMBED_RC2",
        data=[["SLIDE-PLANE",],["ORI_EMBED_RC2",]],)
    model.orientation.updateItem(
        "ORI_EMBED_RC1", norm(offset1), norm(offset2))
    model.orientation.updateItem(
        "ORI_EMBED_RC2", norm(cross(offset1, offset2)), norm(offset2))
    model.properties.updateSection(
        type="SOLIDSECTION",
        elset="EMBEDDED_RUBBER_RC",
        material="EMBEDDED_RUBBER",
        data=[[1e-4],])
    if "EMBEDDED_RUBBER" not in model.material:
        model.material.updateItem("EMBEDDED_RUBBER", {
            "DampingAlpha": 0.5,
            "Density": 1.0,
            "EMod": 10,
            "Nue": 0.3,
        })
    return elsetRubber


def checkSectionFirstAxis(model):
    """To be implemented... Check default approximate n1 axis (0,0,-1) against
    node normals.

     - t ... element tangent vector
     - n_1 ... approximate first axis of the beam cross section
     - n_2 ... node normal and second axis of the beam cross section

    If t x n_1 is more than 90 deg from the normal n_2 then n_2 will be flipped
    by Abaqus. To prevent this checkSectionFirstAxis suggests a different n_1
    for affected beam elements.
    """
    pass


############################################################################

class AssemblyModel(Model):
    """Subclass of L{bae.abq_model_02.Model} to help assembling
    ground support models consisting of various parts.
    """

    def __init__(self, firstNode=5000001, firstElem=5000001):
        Model.__init__(self)
        self.nextNode = firstNode
        self.nextElem = firstElem

    def addPart(self, newPart):
        """
        @param newPart: a file name to an Abaqus input file or a
            L{bae.abq_model_02.Model}-object that is to be inserted into self.

        @returns: Tuple of two sets (newNodes, newElems) with newly added node
           and element numbers.
        """

        if isinstance(newPart, basestring):
            m = Model().read(newPart)
        else:
            m = newPart.copy()

        (nodesOldToNew, elemsOldToNew) = m.renumber(
            nodeStart=self.nextNode, elemStart=self.nextElem)
        newNodes = set(nodesOldToNew.itervalues())
        newElems = set(elemsOldToNew.itervalues())
        self.insertModel(m)

        # update firstNode and firstElem for next part
        self.nextNode = findNextRoundNumber(max(self.nodeCoords), 100000) + 1
        self.nextElem = findNextRoundNumber(max(self.elNodes), 100000) + 1

        return (newNodes, newElems)

    def writeEmbeddedElementCard(
            self, resultFile,
            hostElsetName="GS_HOST", absExtTol=0.02, embeddedElsets=None):
        """As long as the *EMBEDDED ELEMENT card is not implemented in the
        Abaqus model use this method to write this card to an output file
        @param resultFile: open file, destination Abaqus input file
        @param hostElsetName: parameter of the HOST ELSET option
        @param absExtTol:parameter of the ABSOLUTE EXTERIOR TOLERANCE option
        @param embeddedElsets: iterable of elset names: which elsets shall be
           embedded in the host elset. If not given then take all elsets
           defined in self with names that end in "EMBEDDED_ELEMS".
        """

        if embeddedElsets is None:
            embeddedElsets = self.regexpNames("elset", ".*EMBEDDED_ELEMS")

        resultFile.write("**\n** embedded elements\n**\n")
        resultFile.write(
            "*EMBEDDED ELEMENT, HOST ELSET=%s, ABSOLUTE EXTERIOR TOLERANCE=%g\n"
            % (hostElsetName, absExtTol))
        for line in groupIter(embeddedElsets):
            resultFile.write("%s\n" % ", ".join(line))
