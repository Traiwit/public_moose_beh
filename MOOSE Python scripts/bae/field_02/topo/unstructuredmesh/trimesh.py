"""module for the class L{TriMesh}.
"""

__version__ = "0.1"

_version_history_ = """\
Versions:
=========

0.1 GP new, based on fieldData_00.topo.mesh.py and surface_03.py v3.11
"""

from copy import deepcopy
import numpy as np
import struct

from .meshbase import HomoMesh
from .elementType import elTypeLinTri, ElementType

from bae.log_01 import msg, MsgTicker
from bae.abq_model_02 import Model as AbqModel, checkModelModuleVersion
from bae.mesh_01 import Mesh as PurePyMesh


class TriMesh(HomoMesh):
    """A surface made of triangles. Can be used to just represent an arbitrary
    surface approximated by linear triangles or for topo classes L{MeshNodes}
    L{MeshElems} and L{MeshGaussPts}.

    @ivar nodeCoords: ... node coords, ndarray Nx3
    @ivar elNodes: ... connectivity, ndarray M x 3
    @ivar elType: L{ElementType}-object that describes properties of the
      elements
    """

    def __init__(self, *args, **kwargs):
        """Initialize a triangle surface.

        The first (positional) argument might be another triangle surface
        (L{TriMesh})-object. In this case a *deep* copy is performed.

        @kwarg nodeCoords: optional, initialize self.nodeCoords with this
        @kwarg elNodes: optional, initialize self.elNodes with this
        @kwarg elType: optional, initialize self.elType with this
        """

        if args and isinstance(args[0], TriMesh):
            # Copy-constructor
            other = args[0]
            self.elNodes = deepcopy(other.elNodes)
            self.nodeCoords = deepcopy(other.nodeCoords)
            self.edgeToTri = deepcopy(other.edgeToTri)
            self.elType = other.elType
        elif args:
            #... need to differentiate:
            # ... and isinstance(args[0], Whatever-Surface):
            raise NotImplementedError(
                "Not yet there: TriMesh from other Surface classes.")
        else:
            try:
                self.nodeCoords = np.asarray(kwargs["nodeCoords"])
            except KeyError:
                self.nodeCoords = np.array([])  # empty ndarray
            try:
                self.elNodes = np.asarray(kwargs["elNodes"])
            except KeyError:
                self.elNodes = np.array([])  # empty ndarray

            try:
                self.elType = kwargs["elType"]
            except KeyError:
                self.elType = elTypeLinTri

            # the following are None as long as they are not calculated yet
            # don't query these variables, use the corresponding get...()
            # method, ie. for self.edgeToTri use self.getEdgeToTri()
            # if any operation renders those data invalid, set the correspondig
            # variable (eg. self.edgeToTri) to None and it will be recalculated
            # automatically
            self.edgeToTri = None


    ######################################################
    #{ inititialize new surface
    @classmethod
    def fromAbqModelTris(cls, abqModel, elset=None):
        """Create a surface from some tri elements of a L{bae.mesh_01.Mesh}- or
        an L{abaqus model<bae.abq_model_02.Model>}-object.

        All kinds of triangle elements in the abqModel (or in the given elsets)
        are considered, non triangle elements are ignored.

        Elements in elset that are not found in abqModel.elNodes are ignored, a
        warning is posted.

        @param abqModel: a L{bae.mesh_01.Mesh}- or an
        L{abaqus model<bae.abq_model_02.Model>}-object or a filename (a string)
        or an open file of an abaqus input file.

        @param elset: The surface shall consist of the elements in this elset.
        If abqModel is an AbqModel then elset might be anything that
        abqModel.getUnionSet() accepts as input: an elset name, a set of element
        numbers, a list of elset names and element numbers. Otherwise it must
        be an iterable of element numbers.

        @note: This creates additional attributes for the TriMesh object:
        elemLabel and nodeLabel.
        """

        # process the abqModel argument: determine the source mesh
        if isinstance(abqModel, (basestring, file)):
            # basestring = any kind of string
            m = AbqModel()
            m.read(abqModel)
            abqModel = m
        elif isinstance(abqModel, AbqModel):
            checkModelModuleVersion(abqModel, "%s.%s.fromAbqModelTris"
                                    % (__name__, cls.__name__))
        elif not isinstance(abqModel, PurePyMesh):
            ValueError(
                "%s.%s.fromAbqModelTris expects some sort of Mesh as first"
                " argument, istead got a %s."
                % (__name__, cls.__name__, type(abqModel)))

        # process elset argument: take which elements from the mesh?
        alltris = (abqModel.shapeEl.get('TRI', set())
                   | abqModel.shapeEl.get('TRI_L', set())
                   | abqModel.shapeEl.get('TRI_Q', set()))
        if elset is None:
            elset = alltris.intersection(abqModel.elNodes)
        else:
            try:
                # see if we have a getUnionSet-method (if it's an AbqModel)
                elset = abqModel.getUnionSet("elset", elset)
            except AttributeError:
                # if not, elset must be a set of element numbers, make a copy
                elset = set(elset)
                pass
            elset.intersection_update(alltris)
            trisNotFound = elset.difference(abqModel.elNodes)
            if trisNotFound:
                msg("WARNING from %s.%s.fromAbqModelTris():"
                    " Could not find element connectivity for %d of the %d tri"
                    " elements in the specified elset."
                    % (__name__, cls.__name__, len(trisNotFound), len(elset)))
                elset.difference_update(trisNotFound)

        # find connected nodes and check that all nodes are defined
        connectedNodes = abqModel.getConnectedNodes(elset)
        if connectedNodes.difference(abqModel.nodeCoords):
            trisMissingNodes = [
                e for e in elset
                if any(n not in abqModel.nodeCoords
                       for nodes in abqModel.elNodes[e]
                       for n in nodes)]
            msg("WARNING from %s.%s.fromAbqModelTris():"
                " For %d of %d tri elements some nodes where not defined in"
                " the model. Ignoring those triangles."
                % (__name__, cls.__name__, len(trisMissingNodes), len(elset)))
            elset.difference_update(trisMissingNodes)

        # determine elType
        types = set(abqModel.elType[e] for e in elset)
        if len(types)==1:
            elType = ElementType(abqName=types.pop())
        else:
            elType = elTypeLinTri

        # initialize data members
        elset = sorted(elset)
        elemLabel = np.array(elset, dtype=int)
        connectedNodes = sorted(connectedNodes)
        nodeLabel = np.array(connectedNodes, dtype=int)

        nodeCoords = np.array(
            [abqModel.nodeCoords[n][:3] for n in connectedNodes], dtype=float)
        nodeLabelToIdx = dict((n, i) for i, n in enumerate(connectedNodes))
        elNodes = np.array([
            [nodeLabelToIdx[n] for n in abqModel.elNodes[e][:3]]
            for e in elset], dtype=int)

        if len(elNodes)==0:
            msg('WARNING: %s.%s.fromAbqModelTris initialized without elements.'
                % (__name__, cls.__name__))

        res = cls(nodeCoords=nodeCoords, elNodes=elNodes)
        res.elemLabel = elemLabel
        res.nodeLabel = nodeLabel
        res.elType = elType

        return res


    @classmethod
    def fromAbqModelSurface(cls, abqModel, surfName):
        """Create a TriangleSurface from the surface named surfName in a
        given ABAQUS model.

        The surface must consist of faces of tet elements (C3D4 or C3D10M).
        Only the corner nodes of the tets are considered.

        The elsets surfName+"_S1", surfName+"_S2", ... "+_S4" define the
        surface.

        The node coords of this surface are deep copied.

        @note: This creates additional attributes for the TriMesh object:
        elemLabel (size N), nodeLabel (size M) and faceType (size N).
        faceType is "S1", "S2, "S3" or "S4" for tet element faces. N is the
        number of faces, M is the number of nodes/vertices.

        @note: There is room for optimization if it's too slow. But mind the
        various consistency checks.
        """

        checkModelModuleVersion(abqModel, "%s.%s.fromAbqModelSurface"
                                % (__name__, cls.__name__))

        triElementFaces = {'S1':[0,1,2],'S2':[0,1,3],'S3':[1,2,3],'S4':[0,2,3]}
        validElements = abqModel.shapeEl.get("TET", set()).union(
            abqModel.shapeEl.get("TETHT", set()))

        # a dict {faceType: set of element labels}
        abqModelSurf = abqModel.surface[surfName]

        # collect faces with element labels and face-id
        elNodes = list()
        elemLabel = list()
        faceTypeList = list()
        ignoredElements = set()
        nbTotalFaces = 0  # nb of potential faces (including those ignored)
        for faceType in set(abqModelSurf).intersection(triElementFaces):
            nodeIds = triElementFaces[faceType]
            elems = abqModel.elset[abqModelSurf[faceType]]
            nbTotalFaces += len(elems)
            invalidElems = elems.difference(validElements)
            ignoredElements.update(invalidElems)  # store for diagnostic output
            elems = sorted(elems.difference(invalidElems))
            elemLabel.extend(elems)
            faceTypeList.extend( [faceType,]*len(elems) )
            elNodes.extend([abqModel.elNodes[el][i] for i in nodeIds]
                           for el in elems)

        # diagnostic output
        if ignoredElements:
            ignoredTypes = ", ".join(sorted(set(
                abqModel.elType[el] for el in ignoredElements)))
            msg("WARNING from TriMesh.fromAbqModelSurface(): %d/%d faces of the"
                " surface %s could not be considered because they are"
                " non-triangular or the element type is not implemented."
                " Ignored element types: %s."
                % (len(elNodes), nbTotalFaces, surfName, ignoredTypes))

        # find connected nodes and check that all nodes are defined
        connectedNodes = set(n for nodes in elNodes for n in nodes)
        if connectedNodes.difference(abqModel.nodeCoords):
            idsMissingNodes = set(
                i for i, nodes in enumerate(elNodes)
                if any(n not in abqModel.nodeCoords for n in nodes))
            msg("WARNING from TriMesh.fromAbqModelSurface():"
                " For %d of %d faces some nodes where not defined in"
                " the model. Ignoring those faces."
                % (len(idsMissingNodes), len(elNodes)))
            N = len(elNodes)
            elNodes = [elNodes[i] for i in range(N)
                       if i not in idsMissingNodes]
            elemLabel = [elemLabel[i] for i in range(N)
                         if i not in idsMissingNodes]
            faceTypeList = [faceTypeList[i] for i in range(N)
                            if i not in idsMissingNodes]

        # create numpy-arrays
        connectedNodes = sorted(connectedNodes)
        nodeCoords = np.array(
            [abqModel.nodeCoords[n][:3] for n in connectedNodes], dtype=float)
        nodeLabelToIdx = dict((n, i) for i, n in enumerate(connectedNodes))
        elNodes = np.array([
            [nodeLabelToIdx[n] for n in nodes] for nodes in elNodes], dtype=int)
        elemLabel = np.array(elemLabel, dtype=int)
        nodeLabel = np.array(connectedNodes, dtype=int)
        faceTypeList = np.array(faceTypeList)

        if len(elNodes)==0:
            msg("WARNING: TriMesh.fromAbqModelSurface initialized"
                " without elements.")

        res = cls(nodeCoords=nodeCoords, elNodes=elNodes)
        res.elemLabel = elemLabel
        res.nodeLabel = nodeLabel
        res.faceType = faceTypeList
        res.elType = elTypeLinTri
        return res


    @classmethod
    def fromSTL(cls, filename):
        """Create a TriangleSurface from some tri elements read from a .stl
        (stereo lithography) file.

        Derived from raw2inp.py and something from the internet
        @param filename: File name of a .stl file to read surface data from.
        """

        coordsToNodes = dict()
        nodes = []
        elements = []

        # read start of file to determine if its a binay stl file or an ascii
        # stl file
        fp=open(filename,'rb')
        h=fp.read(80)
        type=h[0:5]
        fp.close()

        if type=='solid':
            msg("reading text file %s" % filename)

            # read text stl match keywords to grab the points to build the model
            fp=open(filename,'r')

            for line in fp:
                words=line.split()
                if len(words)>0:
                    if words[0]=='solid':
                        msg("Reading an object called %s from %s"
                            % (words[1], filename))
                    if words[0]=='facet':
                        triangle=[]
                        # ignoring the normal information
                        # normal=(eval(words[2]),eval(words[3]),eval(words[4]))

                    if words[0]=='vertex':
                        coords = tuple(float(x) for x in words[1:4])
                        try:
                            node = coordsToNodes[coords]
                        except KeyError:
                            node = len(nodes)
                            nodes.append(coords)
                            coordsToNodes[coords] = node
                        triangle.append(node)

                    if words[0]=='endloop':
                        # make sure we got the correct number of values before
                        # storing
                        if len(triangle)==3:
                            elements.append(triangle)
                        else:
                            msg("WARNING: found facet with %d vertices in %s"
                                % (len(triangle), filename))
            fp.close()

        else:
            msg("reading binary stl file %s" % filename)

            # Load binary stl file:
            # Check wikipedia for the binary layout of the file.
            # We use the struct library to read in and convert binary data
            # into a format we can use.
            fp=open(filename,'rb')
            h=fp.read(80)

            l=struct.unpack('I',fp.read(4))[0]
            def parse_pt(bytes):
                if len(bytes)==12:
                    return tuple(struct.unpack('f',bytes[slice(i, i+4)])[0]
                                 for i in (0, 4, 8))
                else:
                    raise EOFError("Read less than 12 bytes for a coordinate"
                                   " tripel in %s." % filename)

            count=0
            while True:
                try:
                    p=fp.read(12)
                    n=parse_pt(p)

                    triangle=[]
                    for i in range(3):
                        p=fp.read(12)
                        coords = parse_pt(p)
                        try:
                            node = coordsToNodes[coords]
                        except KeyError:
                            node = len(nodes)
                            nodes.append(coords)
                            coordsToNodes[coords] = node
                        triangle.append(node)

                    elements.append(triangle)
                    count+=1
                    fp.read(2)

                except EOFError:
                    break
            fp.close()

        # copy data into TriMesh object
        if len(elements)==0:
            msg('WARNING: TriangleSurface.fromSTL initialized'
                ' without elements.')
        return cls(nodeCoords=np.array(nodes, dtype=float),
                   elNodes=np.array(elements, dtype=int),
                   elType=elTypeLinTri)


    @classmethod
    def fromCohesives(cls, abqModel, elset=None, elType="COH3D6"):
        """Create a L{TriMesh} from some cohesive elements COH3D6 of an
        Abaqus model

        Cohesive elements are considered, other elements are ignored.

        Elements in elset that are not found in abqModel.elNodes are ignored, a
        warning is posted.

        @param abqModel: an L{abq_model_02.Model} object or a filename (a string)
        or an open file of an abaqus input file.

        @param elset: The surface shall consist of the elements in this elset.
        Might be anything that abqModel.getUnionSet() accepts as input:
        An elset name, a set of element numbers, a list of elset names
        and element numbers.
        Elset might also be None, then all COH3D6 elements of abqModel are used.

        @param elType: Can be used to override the element type to be
        considered instead of the COH3D6 - cohesive elements.
        """

        if isinstance(abqModel, (basestring, file)):
            # basestring = any kind of string
            m = AbqModel()
            m.read(abqModel)
            abqModel = m
        else:
            checkModelModuleVersion(abqModel, "%s.%s.fromCohesives"
                                    % (__name__, cls.__name__))

        if elset is None:
            elset = abqModel.typeEl[elType]
        else:
            elset = abqModel.getUnionSet("elset", elset)
            elset.intersection_update(abqModel.typeEl[elType])

        # find connected nodes and check that all nodes are defined
        connectedNodes = set(
            n for elem in elset for n in abqModel.elNodes[elem][:3])
        if connectedNodes.difference(abqModel.nodeCoords):
            trisMissingNodes = [
                e for e in elset
                if any(n not in abqModel.nodeCoords
                       for nodes in abqModel.elNodes[e]
                       for n in nodes)]
            msg("WARNING from %s.%s.fromCohesives():"
                " For %d of %d elements some nodes where not defined in"
                " the model. Ignoring those elements."
                % (__name__, cls.__name__, len(trisMissingNodes), len(elset)))
            elset.difference_update(trisMissingNodes)

        # initialize data members
        elset = sorted(elset)
        elemLabel = np.array(elset, dtype=int)
        connectedNodes = sorted(connectedNodes)
        nodeLabel = np.array(connectedNodes, dtype=int)

        nodeCoords = np.array(
            [abqModel.nodeCoords[n][:3] for n in connectedNodes], dtype=float)
        nodeLabelToIdx = dict((n, i) for i, n in enumerate(connectedNodes))
        elNodes = np.array([
            [nodeLabelToIdx[n] for n in abqModel.elNodes[e][:3]]
            for e in elset], dtype=int)

        if len(elNodes)==0:
            msg('WARNING: %s.%s.fromCohesives initialized without elements.'
                % (__name__, cls.__name__))

        res = cls(nodeCoords=nodeCoords, elNodes=elNodes)
        res.elemLabel = elemLabel
        res.nodeLabel = nodeLabel
        res.elType = elTypeLinTri
        return res


    @classmethod
    def convexHullFromPoints(cls, points):
        """Create a TriMesh that is the convex hull of the given points.
        """
        try:
            import scipy.spatial
        except ImportError:
            raise NotImplementedError(
                "%s.%s.convexHullFromPoints requires scipy.spatial"
                " to be installed." % (__name__, cls.__name__))

        qhull = scipy.spatial.ConvexHull(points)
        nodeCoords = np.array(
            [points[idx] for idx in qhull.vertices], dtype=float)
        elNodes = np.array(qhull.simplices, dtype=int)
        return cls(nodeCoords=nodeCoords, elNodes=elNodes)

    #{ end of inititialize new surface

    ######################################################
    #{ file I/O

    ###### CONSTRUCTION SITE AHEAD ####################
    #### functions still missing

    #{ end of file I/O


    ######################################################
    #{ methods for data extraction / getting info
    def __str__(self):
        """Give a short description of the surface.

        >>> surf = TriMesh.fromSTL("example.stl")
        >>> print surf
        TriMesh with 214 nodes, 201 triangles, ...
        ... min: [0.0,0.4,1.2], max: [5.3,8.4,3.2]

        """
        bbox = self.getBoundingBox()
        return ("TriMesh with %d nodes, %d triangles, min: [%s], max: [%s]"
                % (len(self.nodeCoords), len(self.elNodes),
                   ",".join(self.coordFormat % x for x in bbox[0]),
                   ",".join(self.coordFormat % x for x in bbox[1])))

    def getCentroid(self):
        """Returns the area-centroid of the surface.

        See also: L{getArea} to get the area of the surface.
        And L{getTriAreaNormal} for the area of individual faces.
        """
        # array of all triangle-vertex-coords: nodeCoords[i,j,k]=...
        # i-th triangle, j-th vertex (0..2) of the triangle, k-th component
        nodeCoords = self.nodeCoords[self.elNodes]

        # centroid for each triangle
        centroids = nodeCoords.sum(axis=1) / 3.0

        # area of each triangle
        a1 = nodeCoords[:,1,:] - nodeCoords[:,0,:]
        a2 = nodeCoords[:,2,:] - nodeCoords[:,0,:]
        areas = 0.5 * np.linalg.norm(np.cross(a1, a2), axis=1)

        # sum up centroids weighted by area and divide by accumulated area
        centroid = (centroids * areas[:,np.newaxis]).sum(axis=0) / areas.sum()

        return centroid

    def getElCentroids(self):
        """
        Calculate element centroids for all elements
        """
        # array of all triangle-vertex-coords: nodeCoords[i,j,k]=...
        # i-th triangle, j-th vertex (0..2) of the triangle, k-th component
        nodeCoords = self.nodeCoords[self.elNodes]

        # centroid for each triangle
        return np.average(nodeCoords, axis=1)
    
    ###### CONSTRUCTION SITE AHEAD ####################
    #### functions still missing

    #} end of methods for data extraction / getting info


    ######################################################
    #{ methods for projecting points and searching points

    ###### CONSTRUCTION SITE AHEAD   ####################
    # This brute force version is too slow.
    #
    # check out fieldData_00.unstructured-topo-mesh... rtree
    #    ... this needs to call the find_point function for each single point
    #    - this implementation uses the rtree module
    #      https://pypi.org/project/Rtree/ which is a ctype-wrapper force
    #      libspatialindex (https://libspatialindex.org/en/latest/)
    #      Currently I'm very reluctant to make field_02 dependent on
    #      libspatialindex. But maybe it makes sense?
    #    - the rtree interface in fieldData_00.numpy_geom seems not well
    #      designed to me (?)
    #
    # or should we reimplement the bae.surface_03 initializePointSearch /
    #    ... pointIsBelow / projectPoints methods?
    #
    # or use scipy KDTree?
    #
    # there is a pure python implementation of an rtree:
    #    ... https://pypi.org/project/rtreelib/
    #    this seems a bit experimental to me...

    def projectPoints(self, points, dir=None, maxChunkSize=10000000.0):
        """projects arbitrary points vertically on the surface

        Returns a list of as many items as there are items in points. Each of
        those items is a list of projected points. It may be empty if the point
        is not below or above the surface or it may contain one or many points
        on the surface.

        If two of the points are closer than self.pointSearchTol only the first
        is recognized.

        @param points: points to project, an Nx3 array of point coordinates.
        @param dir: projection direction. Defaults to the z-axis.
            Can be 0, 1, 2 for x-,y-,z-axis. Or a vector. Not yet recognized.
        @param maxChunkSize: don't touch, for testing purposes only:
            max allowed matrix size in bytes *as float*.
        @return: list of (N_i x 3) arrays. One item/array per given point.
        """

        # Todo: rotate system to account for dir argument
        if dir is not None and dir!=2:
            raise NotImplementedError(
                "ERROR: dir argument for %s.%s.projectPoints not implemented,"
                " yet." % (__name__, cls.__name__))

        # array of all triangle-vertex-coords: nodeCoords[i,j,k]=...
        # i-th triangle, j-th vertex (0..2) of the triangle, k-th component
        nodeCoords = self.nodeCoords[self.elNodes]

        # base vectors of each triangle as matrix
        # base vec 1: a1 = nodeCoords[:,0] - nodeCoords[:,2]
        # base vec 2: a2 = nodeCoords[:,1] - nodeCoords[:,2]
        # we need the 2D part of it for ainv and the full matrix later
        amat = nodeCoords[:,0:2] - nodeCoords[:,2:3]  # Nx2x3

        # filter tris with determinant almost zero
        valid = np.abs(np.linalg.det(amat[:,:,:2]))>1e-6
        nodeCoords = nodeCoords[valid]
        amat = amat[valid]
        msg("Filtered too small triangles, %d remaining of %d."
            % (len(nodeCoords), len(valid)), debugLevel=10)

        # invert matrix of x-y-components of base vectors
        # amat_ji^2D xi_j = x_i  <=>  ainv_ji x_j = xi_i
        ainv = np.linalg.inv(amat[:,:,:2])[np.newaxis,:,:,:]  # 1xNx2x2

        # base point (third node) of each triangle, only x,y components
        basePts = nodeCoords[:,2,:2]  # Nx2
        ### print "### basePts :\n%s" % basePts   ##### DEBUG

        # result list, will be a list of (N_i x 3) arrays corresponding to the
        # items in the points array, initialize with empty array
        result = [np.array([]),] * len(points)

        # divide points into chunks such that the size of the intermediate
        # matrix stays below 10MB
        N = len(amat)
        bytesN = ainv.nbytes/4  # nb of triangles * sizof(float)
        M = len(points)  # nb of points
        bytesRes = bytesN * M * 2
        nbChunks = max(1.0, np.ceil(bytesRes / maxChunkSize))  # float constants!
        MM = int(np.ceil(M / nbChunks))  # MM = chunk size
        nbChunks = int(nbChunks)
        msg("Processing %d points in %d chunks of size %d."
            % (M, nbChunks, MM), debugLevel=10)

        ### print "amat.shape %s, ainv.shape %s" % (amat.shape, ainv.shape)   ##### DEBUG
        ticker = MsgTicker("...processing chunk %%d of %d" % nbChunks)
        tempMat = np.zeros((MM, N, 2))
        elemCoords = np.zeros((MM, N, 2))
        for iChunk in range(nbChunks):
            ticker.tick()
            testPts = points[iChunk*MM:(iChunk+1)*MM,:2]  # MM x 2
            ### print "### testPts :\n%s" % testPts   ##### DEBUG
            np.subtract(  # store the point position relative to base nodes
                testPts[:,np.newaxis,:], basePts[np.newaxis,:,:],
                out=tempMat)  # MM x N x 2
            ### print "### tempMat MMxNx2 point pos relative to base nodes :\n%s" % tempMat   ##### DEBUG
            ### print "### tempMat.shape, ainv.shape, elemCoords.shape:", tempMat.shape, ainv.shape, elemCoords.shape    ##### DEBUG
            # matmul(a, b)[i,j,k] = sum(a[i,j,:] * b[i,j,:,k])
            elemCoords = np.einsum(  # MM x N x 2
                "ijk,ijkm->ijm",
                tempMat,  # MM x N x 2
                ainv,  # MM x N x 2 x 2
                out=elemCoords)
            ### print "### elemCoords:\n%s" % elemCoords    ##### DEBUG

            # store xi+eta in tempMat[:,:,0]
            # this must be <=1 in order for the point to sit in the triangle
            # the third element coordinate would be 1-xi-eta.
            np.add(elemCoords[:,:,0], elemCoords[:,:,1], out=tempMat[:,:,0])
            ### print "### third elemCoord:\n%s" % tempMat[:,:,0]    ##### DEBUG

            # create IndexMat
            # np.nonzero(M_ij) returns two vectors (size = nb of non-zero items
            # of M_ij): the first vector holds the i-indexes, the second the
            # j-indexes
            ptIds, triIds = np.nonzero(np.logical_and(np.logical_and(
                elemCoords[:,:,0]>=0, elemCoords[:,:,1]>=0), tempMat[:,:,0]<=1))
            ### print "### ptIds:%s" % ptIds    ##### DEBUG
            ### print "### triIds:%s" % triIds    ##### DEBUG

            # now collect items for each point
            for ptId in np.unique(ptIds):
                tIds = triIds[ptIds==ptId]  # tri ids only for this point
                ptIdGlobal = ptId + iChunk*MM  # id in points
                # 
                ### print "### ptId:%s, tIds: %s" % (ptId, tIds)    ##### DEBUG
                ### print "### ptIdGlobal:%s" % ptIdGlobal    ##### DEBUG
                ### print "### amat[tIds,:]=%s" %  amat[tIds,:]   ##### DEBUG
                ### print "### elemCoords[ptId,tIds,:]=%s" % elemCoords[ptId,tIds,:]    ##### DEBUG
                result[ptIdGlobal] = (
                    np.einsum("ijk,ij->ik", amat[tIds,:], elemCoords[ptId,tIds,:])
                    + nodeCoords[tIds,2,:])

        del ticker

        # Todo: rotate results to account for dir argument
        return result
    
    #} end of methods for projecting points and searching points

    ###### CONSTRUCTION SITE AHEAD   ####################
    ###### CONSTRUCTION SITE AHEAD   ####################


    ######################################################
    #{ file I/O
    def exportAsSTL(self, fileName):
        """Export the triangle surface to the new file with the specified file
        name in .stl (stereolithography) format. This is suitable for import as
        mesh into Rhino or as geometry into icem.

        The binary format is exported. The normals are not exported.

        It might still be advisable to call self.unifyNormals() in advance.
        This will ensure that the vertices (nodes) are ordered in a consistent
        way. Not necessarily with normals pointing outwards as STL requires,
        though.

        @Note: Existing files with the same name are silently overwritten.
        """

        msg("Writing binary stl file %s" % fileName)

        # Write binary stl file:
        # Check wikipedia for the binary layout of the file.
        # Uses the struct library to convert into binary form.
        fp=open(fileName,'wb')

        header = struct.pack('80s', "binary stl")
        fp.write(header)

        nbTris = struct.pack('I', len(self.elNodes))
        fp.write(nbTris)

        emptynormal = struct.pack("<3f", 0.0, 0.0, 0.0)
        emptyattrib = struct.pack("2s", "")
        ptpacker = struct.Struct("<3f")

        for nodes in self.elNodes.itervalues():
            assert(len(nodes)==3)
            fp.write(emptynormal)
            for node in nodes:
                coords = self.nodeCoords[node]
                fp.write(ptpacker.pack(*coords))
            fp.write(emptyattrib)

        fp.close()
        msg("Finished writing binary stl file %s" % fileName)
        return

    def exportAsIv(self, fileName, colour=None, transparency=0):
        """Export the triangle surface to the new file with the specified file
        name in OpenInventor (.iv) format. This is suitable for import as
        shapes into voxler.

        @param fileName: ... guess
        @param colour: RGB-tuple, three floats between 0 and 1
        @param transparency: float between 0 and 1
        """

        msg("Writing OpenInventor ascii file %s" % fileName)

        # find all connected nodes
        allnodes = set()
        for nodes in self.elNodes.itervalues():
            allnodes.update(nodes)
        allnodes = sorted(allnodes)

        # point list
        pointList = list()
        nodeToPtIdx = dict()
        for cnt, node in enumerate(allnodes):
            coords = self.nodeCoords[node]
            pointList.append(coords)
            nodeToPtIdx[node] = cnt

        # cell list
        triangles = [
            [nodeToPtIdx[node] for node in nodes]
            for nodes in self.elNodes.itervalues()
            ]

        output = open(fileName, "w")
        output.write("#Inventor V2.0 ascii\n")
        output.write("Separator {\n")

        # write material
        if colour is not None or transparency>0:
            if colour is None:
                colour = (1,1,1)
            else:
                colour = tuple(colour)
            output.write("  Material {\n")
            output.write("    ambientColor 0 0 0\n")
            output.write("    diffuseColor %g %g %g\n" % colour)
            output.write("    specularColor 0 0 0\n")
            output.write("    emissiveColor 0 0 0\n")
            output.write("    shininess 1\n")
            output.write("    transparency %g\n" % transparency)
            output.write("  }\n")
            output.write("  ShapeHints {\n")
            output.write("    vertexOrdering COUNTERCLOCKWISE\n")
            output.write("    shapeType UNKNOWN_SHAPE_TYPE\n")
            output.write("    faceType CONVEX\n")
            output.write("    creaseAngle 0\n")
            output.write("  }\n")

        # write point / node data
        output.write("  Coordinate3 {\n    point [\n")
        for coords in pointList:
            output.write("%g %g %g,\n" % tuple(coords))
        output.write("    ]\n  }\n")

        # write face / triangle data
        output.write("  IndexedFaceSet {\n    coordIndex [\n")
        for pts in triangles:
            output.write("%d,%d,%d,-1,\n" % tuple(pts))
        output.write("    ]\n  }\n")

        # fini
        output.write("}\n")
        output.close()
        msg("Finished writing OpenInventor file %s"% fileName)
        return

    def exportAsVtk(self, fileName, objectName="Object1", format='binary'):
        """Export the triangle surface to the specified file in lagacy .vtk
        format. This is suitable for paraview and others. Not for voxler.

        @param fileName: ... guess!
        @param objectName: Name parameter in the vtk file. This parameter
        accepts a string.
        @param format: "binary" or "ascii"
        @Note: Existing files with that name are being overwritten silently.
        """
        # point list
        pointList = list()
        nodeToPtIdx = dict()
        for cnt, (node, coords) in enumerate(self.nodeCoords.iteritems()):
            pointList.append(coords)
            nodeToPtIdx[node] = cnt

        # cell list
        triangles = [
            [nodeToPtIdx[node] for node in nodes]
            for nodes in self.elNodes.itervalues()]

        # export to vtk
        vtkGrid = UnstructuredGrid( pointList, triangle=triangles )
        v = VtkData(vtkGrid, objectName)
        v.tofile(fileName, format)

        # fini
        msg("Finished writing vtk file %s"% fileName)

    #{ end of file I/O


    ######################################################
    #{ methods for data extraction / getting info

    def getTriNormal(self):
        """Returns a dict of normal vectors of length one:
        {element number: unit normal vector}.

        You may use L{unifyNormals} beforehand.
        """
        triNormal = dict(
            (el, getTriNormalSc(self.nodeCoords, nodes))
            for el, nodes in self.elNodes.iteritems())
        return triNormal

    def getTriAreaNormal(self):
        """Returns a dict: {element number: normal vector}.
        The length of each normal vector corresponds to the area of the
        triangle.

        Get the surface area of the whole TriangleSurface object:
         >>> from bae.vecmath_01 import length
         >>> area = sum(length(v)
         >>>            for e,v in surf.getTriAreaNormal().iteritems())

        See also: L{getArea}
        """
        triNormal = dict(
            (el, [0.5*x for x in self.getTriNormalUsc(nodes)])
            for el, nodes in self.elNodes.iteritems())
        return triNormal

    def getNodeNormal(self):
        """Calculate nodal B{unit} normals as arithmetic mean of the normals of
        all attached triangle elements. A warning will be issued if
        self.L{unifyNormals}() has not been called beforehand.

        @Note: For a node on an edge in self.flipEdges (see L{unifyNormals})
        the normal will be oriented according to the triangle element with the
        lowest number. Elements on the other side of flipEdges will contribute
        to the node normal with their inverse element normal.

        @Note: If the tri elements connected to a particular node do not form
        exactly one group of edge-connected tri elements -i.e. if one or more
        tri elements are separate from the rest of the tri elements on this
        particular node- then the nodal normal on this particular node might
        be unreasonable. This is especially true if this node lies on one of
        the flipEdges.

        @Returns: dict {node number -> averaged normal of the surface}
        """
        triNormal = self.getTriNormal()
        nodeElems = self.getNodeElems()
        nodeNormalDict = dict()
        try:
            flipEdges = self.flipEdges
        except AttributeError:
            msg("WARNING: Calculating node normals before tri element normals"
                " have been unified is very likely to give unreasonable"
                " results. If element normals are guaranteed to be consistent"
                " for other reasons then consider to add an empty set as"
                " flipEdges attribute to the TriangleSurface object as marker"
                " to suppress this warning.")
            flipEdges = set()

        # nodes on flipEdges: only needed to check if the flipEdges-algorithm
        # is required on a particular node
        flipNodes = set(node for edge in flipEdges for node in edge)

        for node, elements in nodeElems.iteritems():
            if node in flipNodes:
                # walk around node considering flipEdges
                # outer nodes on each element
                elOuterNodes = dict(
                    (el, set(self.elNodes[el]).difference((node,)))
                    for el in elements)
                # outer node to connected elements
                nodeTwoElems = defaultdict(set)
                for (el, nodes) in elOuterNodes.iteritems():
                    for node2 in nodes:
                        nodeTwoElems[node2].add(el)

                # initialize start at element with smallest number
                startElNum = min(elements)
                el = startElNum
                node2, endNode = elOuterNodes[el]
                rightSide = True
                nodeNormal = list(triNormal[el])  # make a copy, will add more
                while node2!=endNode:

                    # switch side if crossing flipEdges
                    if frozenset((node2, node)) in flipEdges:
                        rightSide = not(rightSide)

                    # proceed to next element
                    try:
                        el = nodeTwoElems[node2].difference((el,)).pop()
                        node2 = elOuterNodes[el].difference((node2,)).pop()
                    except KeyError:
                        # reached end (border) of surface
                        if startElNum is None:
                            # already went in both directions stop it
                            break

                        # continue from start in opposite direction
                        el = startElNum
                        startElNum = None  # don't turn around a second time!
                        node2 = endNode

                        # switch side if crossing flipEdges
                        if frozenset((node2, node)) in flipEdges:
                            rightSide = not(rightSide)

                        try:
                            el = nodeTwoElems[node2].difference((el,)).pop()
                            node2 = elOuterNodes[el].difference((node2,)).pop()
                        except KeyError:
                            break

                    # sum up tri normals
                    if rightSide:
                        vector_modif_add(nodeNormal, triNormal[el])
                    else:
                        vector_modif_add_scaled(nodeNormal, triNormal[el], -1)
            else:
                nodeNormal = [0,0,0]
                for el in elements:
                    nodeNormal = vector_plus(nodeNormal,triNormal[el])
            scale = length(nodeNormal)
            if scale>0:
                nodeNormalDict[node] = vector_scale(nodeNormal, 1.0/scale)
            else:
                nodeNormalDict[node] = nodeNormal

        # return the dict
        return nodeNormalDict

    def getArea(self):
        """Returns the accumulated area of the surface.

        See also: L{getTriAreaNormal} for the area of individual faces.
        """
        return 0.5 * sum(length(self.getTriNormalUsc(nodes))
                         for nodes in self.elNodes.itervalues())

    def getVolume(self):
        """Returns the volume of the surface.

        The surface obviously must be closed and simply connected -i.e. it must
        be a manifold. Additionally all normals must be unified. Otherwise the
        returned result will be meaningless.

        This subroutine does not check anything! Use L{findBorderEdges}
        and L{unifyNormals} to check.

        The returned volume will be positive if the normals point outwards and
        negative if the normals point inwards. Use abs() if in doubt.
        """

        # *** z component of normal ***
        # ab = vector(coords[0], coords[1])
        # ac = vector(coords[0], coords[2])
        # normal[2] = cross(ab, ac)[2] = -ab[1]*ac[0]+ab[0]*ac[1]
        #  = -(coords[1][1]-coords[0][1])*(coords[2][0]-coords[0][0])
        #   + (coords[1][0]-coords[0][0])*(coords[2][1]-coords[0][1])
        # *** z component of centroid ***
        # centroid[2] = (coords[0][2] + coords[1][2] + coords[2][2])/3
        # sum up the product of (0.5*normal[2]*centroid[2]) for all faces

        return sum(
            (-(coords[1][1]-coords[0][1])*(coords[2][0]-coords[0][0])
             +(coords[1][0]-coords[0][0])*(coords[2][1]-coords[0][1]))
            * (coords[0][2] + coords[1][2] + coords[2][2]) / 6.0
            for coords in (
                    [self.nodeCoords[n] for n in nodes]
                    for nodes in self.elNodes.itervalues()) )

    #} end of methods for data extraction / getting info
