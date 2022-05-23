"""surface_03.py

surface classes
"""

__version__ = "3.14"

_version_history_ = """\
Versions:
=========

3.01 GP new from surface_02.py v2.11
3.02 GP moved getConnectedNodes() to mesh_01.Mesh
3.03 GP added getTriAreaNormal()
3.04 GP added getCentroid(), getArea(), getVolume()
3.05 GP added fromPointsTriFaces(), convexHullFromPoints(), offset()
3.06 GP added ElemFaceSurface.invert() and ElemFaceSurface.addSurfFromElsetHole
3.07 GP added ElemFaceSurface-property faceIdShapeCornerNodeIds and two methods
        differenceUpdateNset and intersectionUpdateNset
3.08 GP added ElemFaceSurface.addFreeSurfFromSplitMesh(), getTetFaceNormal()
3.09 GP added ElemFaceSurface.difference() and ElemFaceSurface.intersection()
3.10 AF added ElemFaceSurface.getHexFaceNormal()
        improved storeInAbqModel() to only consider non-empty elsets
3.11 GP added ElemFaceSurface.updateFromSequence()
        fixed ElemFaceSurface.difference_update() to not have empty sets in
        faceEl
3.12 GP added pointIsBelow() axis-argument, changed internal variables:
        removed pointSearchInitialized and separate point-search-data variables
        boundingBox, cellNumber, cellSize and cells
3.13 GP changed ElemFaceSurface.invert() now recognizes arbitrary elements,
        not only tets
3.14 GP added returnVolRef option to ElemFaceSurface.getTriangleSurface()
        and ElemFaceSurface.getTetFaceCentroid()
"""

## ---------------------------------------------------------------------------
_bugs_ = """\
known bugs/problems:
====================
- exportAsVtk() issues a warning from the underlying pyvtk library that can be
  ignored. It complains about the mesh having no data attribute, which
  apparently is fine for paraview. Suppressing this warning is not easily
  possible. Would have to modify the pyvtk library.
"""

_todo_ = """\
ideas / todo:
=============

- make simplyConnectedParts really return only simply connected parts
  (currently it ignores connections by T-edges but also returns surface parts
  that include T-edges if all triangles are also connected by normal edges on
  different paths)
  . incorporate findConnectedParts into simplyConnectedParts.
    findConnectedParts is only being called by simplyConnectedParts
  . replace the algorithm in findConnectedParts by the one in unifyNormals
    so we get unified normals and flipEdges when collecting the simply
    connected parts
  . should this feature be optional?
  ...



- more test routines
  . for getNodeNormal

- from getNodeNormal extract walk-around-the-split-node to separate function
  that creates a structure {node : (start-elem, list of this-side elems, list
  of other-side elems)}
  This structure can later be used to assign tet nodes to either side of the
  surface even on splitEdges.
   . We'd need an additional dict for tets on flip nodes that are
     "late-assigned" (don't share a face with the surface) stating from which
     "early assigned" tet their assignment stems.
   . I.e.: for every early-assigned elem in the node: (start elem, ...)
     structure store this

- rename splitForSurfToContact

- replace simplyConnectedParts() with new routine: newsimplyConnectedParts ?

- simplyConnectedParts(): at edges with more than two tris instead of
  splitting it up in any case accept a connection with an angle of less then
  5 deg if all other intersections are > 30 deg.

- fromVRML: has been started some years ago but never finished. See
  surface_02.SurfaceFromVRML - class.
"""

from collections import defaultdict, deque  # deque = double ended queue
from itertools import izip
from math import cos, pi, sqrt
import struct
from copy import deepcopy

from bae.future_01 import classproperty
from bae.mesh_01 import Mesh
from bae.abq_model_02 import Model as AbqModel, checkModelModuleVersion

from bae.vecmath_01 import vector, vector_plus, vector_scale, vector_sum, \
    vector_modif_add, vector_modif_add_scaled,\
    norm, dot, cross, length, dist
from bae.misc_01 import BoundingBox, Container
from pyvtk import UnstructuredGrid, VtkData
from bae.log_01 import msg, MsgTicker
# IMPORTANT: look at the end of the file, there is yet another import


## ---------------------------------------------------------------------------

def getTriNormalSc(nodeCoords, nodes):
    """Computes the scaled ("Sc") normal of a triangle face. I.e. this normal
    vector is of length 1.

    @param nodeCoords: a dict {node number: coords} like
       L{bae.mesh_01.Mesh.nodeCoords}

    @param nodes: A list or tuple of node numbers identifying the triangle face.
       Nodes come in counter-clockwise succession if you look from "above"
       (positive side), against the direction of the normal.
    """
    vectorAB=vector(nodeCoords[nodes[0]], nodeCoords[nodes[1]])
    vectorAC=vector(nodeCoords[nodes[0]], nodeCoords[nodes[2]])
    return norm(cross(vectorAB,vectorAC))

class SurfaceConnectionError(Exception):
    """raised when there is an error in the connection of the elements
    forming the surface"""
    pass

## ---------------------------------------------------------------------------


class Surface(object):
    "Base class for surfaces."
    pass

## ---------------------------------------------------------------------------


class TriangleSurface(Mesh, Surface):
    """
    A surface made of triangles.

    "smooth mode": When finding points above or below the surface with
    L{pointIsBelow} or L{projectPoints} "smooth mode" additionally considers
    points to be below the surface that are in fact outside but within the
    L{smoothFindRadius} to the closest centre point of any of the triangles
    of the surface (in the x-y-plane, or the plane perpendicular to the
    chosen axis).

    @ivar smoothFindRadius: Can be set to a radius>0.0 in order to activate
      "smooth mode"
    @ivar elNodes: is a dictionary: {element number:  list of node numbers}.
      It contains only the elements that actually form the surface.
    @ivar nodeCoords: is a reference to the dictionary: {node number:
      coordinates of the node} as given to the constructor.
      If you move the original nodes, the surface moves also.
    """
    def __init__(self, *args, **kwargs):
        """Initialize a triangle surface.

        The first (positional) argument might be:
         - Another triangle surface object. In this case a *deep* copy is
           performed.
         - A L{bae.abq_model_02.Model}-object or a filename or an open file
           object of an Abaqus input file. A keyword argument elset can be
           provided that will be passed through
           L{model.getUnionSet()<bae.abq_model_02.Model.getUnionSet>}
           to specify a subset to be considered.

        @kwarg nodeCoords: optional, initialize self.nodeCoords with this
        @kwarg elNodes: optional, initialize self.elNodes with this
        @kwarg elset: If the first (positional) argument identifies an
           L{bae.abq_model_02.Model}-object then this keyword argument might
           be given additionally to specify elements to be considered:
           The surface shall consist of the elements in this elset. Might be
           anything that model.getUnionSet() accepts as input:
           an elset name, a set of element numbers, a list of elset names
           and element numbers

        @note: The triangles are assumed to be simply connected
        (no T-junctions, all connected, no separate parts).
        """

        if args and isinstance(args[0], TriangleSurface):
            # Copy-constructor
            other = args[0]
            self.elNodes = deepcopy(other.elNodes)
            self.nodeCoords = deepcopy(other.nodeCoords)
            self.edgeToTri = deepcopy(other.edgeToTri)
        elif args and isinstance(args[0], Surface):
            #... need to differentiate
            raise NotImplementedError(
                "Not yet there: TriangleSurface from other Surface classes.")
        else:
            try:
                self.nodeCoords = kwargs["nodeCoords"]
            except KeyError:
                self.nodeCoords = dict()
            try:
                self.elNodes = kwargs["elNodes"]
            except KeyError:
                self.elNodes = dict()

            # the following are None as long as they are not calculated yet
            # don't query these variables, use the corresponding get...()
            # method, ie. for self.edgeToTri use self.getEdgeToTri()
            # if any operation renders those data invalid, set the correspondig
            # variable (eg. self.edgeToTri) to None and it will be recalculated
            # automatically
            self.edgeToTri = None

        # for compatibility with Mesh
        self.shapeEl = FakeShapeEl(self.elNodes, shape="TRI_L")
        self.elShape = FakeElShape(self.elNodes, shape="TRI_L")

        # other flag
        self.pointSearchData = [False, False, False]
        self.smoothFindRadius = False

    ######################################################
    #{ inititialize new surface
    @classmethod
    def fromMeshTris(cls, model, elset=None):
        """Create a surface from some tri elements of a L{bae.mesh_01.Mesh}- or
        an L{abaqus model<bae.abq_model_02.Model>}-object.

        All kinds of triangle elements in the model (or in the given elsets)
        are considered, non triangle elements are ignored.

        The self.nodeCoords attribute is taken from the model, no copy is made.
        You may use self.L{cleanNodeCoords}() to make a (shallow) copy and
        strip nodes not needed for the surface.

        The node connectivity from model.elNodes however is deep copied: The
        node lists referenced by the dict self.elNodes are others then those
        referenced by model.elNodes.

        Elements in elset that are not found in model.elNodes are ignored, a
        warning is posted.

        @param model: a L{bae.mesh_01.Mesh}- or an
        L{abaqus model<bae.abq_model_02.Model>}-object or a filename (a string)
        or an open file of an abaqus input file.

        @param elset: The surface shall consist of the elements in this elset.
        If model is an AbqModel then elset might be anything that
        model.getUnionSet() accepts as input: an elset name, a set of element
        numbers, a list of elset names and element numbers. Otherwise it must
        be an iterable of element numbers.
        """

        if isinstance(model, (basestring, file)):
            # basestring = any kind of string
            m = AbqModel()
            m.read(model)
            model = m
        elif isinstance(model, AbqModel):
            checkModelModuleVersion(model, "%s.%s.fromMeshTris"
                                    % (__name__, cls.__name__))
        elif not isinstance(model, Mesh):
            ValueError(
                "%s.%s.fromMeshTris expects some sort of Mesh as first"
                " argument, istead got a %s."
                % (__name__, cls.__name__, type(model)))

        alltris = (model.shapeEl.get('TRI', set())
                   | model.shapeEl.get('TRI_L', set())
                   | model.shapeEl.get('TRI_Q', set()))
        if elset is None:
            elset = alltris.intersection(model.elNodes)
        else:
            try:
                # see if we have a getUnionSet-method (if it's an AbqModel)
                elset = model.getUnionSet("elset", elset)
            except AttributeError:
                # if not, elset must be a set of element numbers, make a copy
                elset = set(elset)
                pass
            elset.intersection_update(alltris)

        nb_trisnotfound = 0
        elNodes = dict()
        for elem in elset:
            try:
                nodes = model.elNodes[elem][:3]
                elNodes[elem] = nodes
            except KeyError:
                nb_trisnotfound += 1
                pass
        if nb_trisnotfound>0:
            msg("WARNING from TriangleSurface.fromMeshTris():"
                " Could not find element connectivity for %d"
                " of the %d tri elements in the specified"
                " elset." % (nb_trisnotfound, len(elset)))
        if len(elNodes)==0:
            msg('WARNING: TriangleSurface.fromMeshTris initialized without'
                ' elements.')

        return cls(nodeCoords=model.nodeCoords, elNodes=elNodes)

    @classmethod
    def fromAbqModelSurface(cls, model, surfName):
        """Create a TriangleSurface from the surface named surfName in a
        given ABAQUS model.

        The surface must consist of faces of tet elements (C3D4 or C3D10M).
        Only the corner nodes of the tets are considered.

        The elsets model+"_S1", model+"_S2", ... "+_S4" define the
        surface.

        The node coords of this surface are a reference (no copy) to those of
        the abaqus model. If the abaqus model moves, so does this surface.
        You may use self.L{cleanNodeCoords}() to make a (shallow) copy and
        strip nodes not needed for the surface.
        """

        checkModelModuleVersion(model, "%s.%s.fromAbqModelSurface"
                                % (__name__, cls.__name__))

        triElementFaces={'S1':[0,1,2],'S2':[0,1,3],'S3':[1,2,3],'S4':[0,2,3]}

        elnum = 1
        elNodes = dict()
        for faceType in set(model.surface[surfName]).intersection(
                triElementFaces):
            nodeIds = triElementFaces[faceType]
            for el in model.elset[
                    model.surface[surfName][faceType]]:
                if model.elType[el] not in ('C3D4', 'C3D10M'):
                    raise NotImplementedError(
                        "%s.%s.__init__() can not create a surface from"
                        " elements of type %s."
                        % (__name__, cls.__name__, model.elType[el]))
                nodes = [model.elNodes[el][i] for i in nodeIds]
                elNodes[elnum] = nodes
                elnum += 1

        if len(elNodes)==0:
            msg('WARNING: TriangleSurface.fromAbqModelSurface initialized'
                ' without elements.')

        return cls(nodeCoords=model.nodeCoords, elNodes=elNodes)

    @classmethod
    def fromPointsTriFaces(cls, points, faces):
        """Create a TriangleSurface from a list of points and a list of faces.
        A face is a list of indexes in the points list.
        """        
        nodeCoords = dict((cnt+1, coords)
                          for cnt, coords in enumerate(points))
        elNodes = dict((cnt+1, [node+1 for node in nodes])
                       for cnt, nodes in enumerate(faces))
        if len(elNodes)==0:
            msg('WARNING: TriangleSurface object initialized'
                ' without elements.')
        return cls(nodeCoords=nodeCoords, elNodes=elNodes)

    @classmethod
    def fromDXF(cls, filename, faceID='3DFACE'):
        """Create a TriangleSurface from some tri elements read from a DXF file

        Derived from dxv2inp.py

        @Note: This does not read DXF files created by rhino (default settings).
        Don't know why.

        @param filename: File name of a DXF file to read surface data from.

        @param faceID: Face identifier in the DXF file.
        """

        inputFile = open(filename)
        lineNb = 0

        nodes = []
        elements = []
        while 1:
            line = inputFile.readline().strip()
            lineNb += 1
            if not line: break

            if line.split()[0]==faceID:

                line = inputFile.readline()  # 2 lines junk
                line = inputFile.readline()
                lineNb += 2

                points = [0,0,0]
                for i in range(4):

                    coords = [0,0,0]
                    for j in range(3):

                        line = inputFile.readline()  # point/xyz counter
                        line = inputFile.readline()  # coord value
                        lineNb += 2
                        coords[j] = float(line)

                    try:
                        nodesIndex = nodes.index(coords)  # node ID
                    except ValueError:
                        nodesIndex = len(nodes)
                        nodes.append(coords)

                    try:
                        if nodesIndex not in points:
                            points[i]=nodesIndex  # element connectivity
                    except IndexError:
                        msg("WARNING: Found more than three different points"
                            " for triangle %d on line %d of the DXF input"
                            " file. Ignoring this point."
                            % (len(elements)+1, lineNb))

                elements.append(points)

        inputFile.close()

        nodeCoords = dict((cnt+1, coords)
                          for cnt, coords in enumerate(nodes))
        elNodes = dict((cnt+1, [node+1 for node in nodes])
                       for cnt, nodes in enumerate(elements))
        if len(elNodes)==0:
            msg('WARNING: TriangleSurface.fromDXF initialized'
                ' without elements.')
        return cls(nodeCoords=nodeCoords, elNodes=elNodes)

    @classmethod
    def fromRaw(cls, filename):
        """Create a TriangleSurface from some tri elements read from a .raw
        file like Rhino creates them.

        Derived from raw2inp.py

        @param filename: File name of a .raw file to read surface data from.
        """

        inputFile = open(filename)
        lineNb = 0

        nodes = []
        elements = []

        # ignore first line
        line = inputFile.readline().strip()
        msg("Reading an object called %s from %s"
            % (line, filename))

        while 1:
            line = inputFile.readline().strip()
            lineNb += 1
            if not line: break

            try:
                raw=map(float, line.split())
            except ValueError:
                # if something is not convertible to a float then
                # it must be another object, ignore this line
                msg("Reading an object called %s from %s"
                    % (line, filename))
                continue

            points = [0,0,0]
            for i in range(3):

                coords = raw[i*3:(i+1)*3]

                try:
                    nodesIndex = nodes.index(coords)  # node ID
                except ValueError:
                    nodesIndex = len(nodes)
                    nodes.append(coords)

                points[i]=nodesIndex  # element connectivity

            elements.append(points)

        inputFile.close()

        nodeCoords = dict((cnt+1, coords)
                          for cnt, coords in enumerate(nodes))
        elNodes = dict((cnt+1, [node+1 for node in nodes])
                       for cnt, nodes in enumerate(elements))
        if len(elNodes)==0:
            msg('WARNING: TriangleSurface.fromRaw initialized'
                ' without elements.')
        return cls(nodeCoords=nodeCoords, elNodes=elNodes)

    @classmethod
    def fromSTL(cls, filename):
        """Create a TriangleSurface from some tri elements read from a .stl
        (stereo lithography or what) file.

        Derived from raw2inp.py and something from the internet
        @param filename: File name of a .raw file to read surface data from.
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

        # copy data into triangle object
        nodeCoords = dict((cnt+1, list(coords))
                          for cnt, coords in enumerate(nodes))
        elNodes = dict((cnt+1, [node+1 for node in nodes])
                       for cnt, nodes in enumerate(elements))
        if len(elNodes)==0:
            msg('WARNING: TriangleSurface.fromSTL initialized'
                ' without elements.')
        return cls(nodeCoords=nodeCoords, elNodes=elNodes)

    @classmethod
    def fromCohesives(cls, model, elset=None, elType="COH3D6"):
        """Create a TriangleSurface from some cohesive elements COH3D6 of an
        abaqus model

        Cohesive elements are considered, other elements are ignored.

        The self.nodeCoords attribute is taken from the model, no copy is made.
        You may use self.L{cleanNodeCoords}() to make a (shallow) copy and
        strip nodes not needed for the surface.

        The node connectivity from model.elNodes however is deep copied: The
        node lists referenced by the dict self.elNodes are others then those
        referenced by model.elNodes.

        Elements in elset that are not found in model.elNodes are ignored, a
        warning is posted.

        @param model: an L{abq_model_02.Model} object or a filename (a string)
        or an open file of an abaqus input file.

        @param elset: The surface shall consist of the elements in this elset
        might be anything that model.getUnionSet() accepts as input:
        an elset name, a set of element numbers, a list of elset names
        and element numbers
        elset might also be None, then all COH3D6 elements of model are used.

        @param elType: Can be used to override the element type to be
        considered instead of the COH3D6 - cohesive elements.
        """

        if isinstance(model, (basestring, file)):
            # basestring = any kind of string
            m = AbqModel()
            m.read(model)
            model = m
        else:
            checkModelModuleVersion(model, "%s.%s.fromCohesives"
                                    % (__name__, cls.__name__))

        if elset is None:
            elset = model.typeEl[elType]
        else:
            elset = model.getUnionSet("elset", elset)
            elset.intersection_update(model.typeEl[elType])

        nb_trisnotfound = 0
        elNodes = dict()
        for elem in elset:
            try:
                nodes = model.elNodes[elem][:3]
                elNodes[elem] = nodes
            except KeyError:
                nb_trisnotfound += 1
                pass
        if nb_trisnotfound>0:
            msg("WARNING from TriangleSurface.fromCohesives():"
                " Could not find element connectivity for %d"
                " of the %d tri elements in the specified"
                " elset." % (nb_trisnotfound, len(elset)))

        if len(elNodes)==0:
            msg('WARNING: TriangleSurface.fromCohesives initialized without'
                ' elements.')

        return cls(nodeCoords=model.nodeCoords, elNodes=elNodes)


    @classmethod
    def convexHullFromPoints(cls, points):
        """Create a TriangleSurface that is the convex hull of the given
        points.
        """
        try:
            import scipy.spatial
        except ImportError:
            raise NotImplementedError(
                "TriangleSurface.convexHullFromPoints requires scipy.spatial"
                " to be installed.")

        qhull = scipy.spatial.ConvexHull(points)
        nodeCoords = dict((idx+1, points[idx])
                          for idx in qhull.vertices)
        elNodes = dict((cnt+1, [node+1 for node in nodes])
                       for cnt, nodes in enumerate(qhull.simplices))
        return cls(nodeCoords=nodeCoords, elNodes=elNodes)

    #{ end of inititialize new surface

    ######################################################
    #{ split surface
    def newsimplyConnectedParts(self):
        """NOT TESTED YET, might eventually replace L{simplyConnectedParts}()
        or will be renamed! Or maybe we don't need this algorithm at all
        because L{simplyConnectedParts}() has all features we need.

        CONSTRUCTION SITE: Docs missing, testing missing, needs different name

        LAST IDEA: Morph algorithm into an alternate unifyNormals-function

        @Note: You might want to L{filterDegeneratedTris}() before calling this
        method.
        """
        edgeToTri = self.getEdgeToTri()
        triNormal = self.getTriNormal()

        sortedEdgeList = [
            (abs(dot(triNormal[tris[0]], triNormal[tris[1]])), edge)
            for edge, tris in edgeToTri.iteritems()
            if len(tris)==2]
        sortedEdgeList.sort(reverse=True)
        sortedEdgeList = [y for x,y in sortedEdgeList]

        # surfaces: {surface id : list of tri element numbers}
        nextSurfId = 1
        surfaces = dict()

        # dict {element number: surface index}
        triSurf = dict()

        # dict {surface id: list of flipped edges}
        # over a flipped edge tri orientation is not consistent
        flipEdges = defaultdict(list)

        for edge in sortedEdgeList:

            # the following creates:
            # twoTris: element numbers of the two tris on the edge
            # connected: surface id, or None if not assigned to a surface yet
            # The tri not assigned to any surface comes first in both lists
            # (if only one tri is assigned).
            twoTris, connected = zip(*sorted(
                ((tri, triSurf[tri])
                 for tri in edgeToTri[edge]),
                key=lambda (x,y): y))

            # nodeIndexes: ...
            # [ [tri_1 node id 1, tri_1 node id 2], [tri_2 n-id 1, t2 n-id 2] ]
            # tri_1 is the possibly unassigned first tri, tri_2 is the other
            # nodeIndexes[i][1]-nodeIndexes[i][0] can be 1 or -2 meaning the
            # tri is "upright" on the edge. Or it can be -1 or 2 meaning the
            # tri is "upside-down" on the edge.
            # If both tris are "upright" or both are "upside-down" then the
            # orientation of both tris is consistent. Otherwise one tri must be
            # flipped.
            nodeIndexes = [
                [self.elNodes[tri].index(edgeNode) for edgeNode in edge]
                for tri in twoTris]
            nodeOrder = [(idx[1]-idx[0]) % 3 for idx in nodeIndexes]
            sameOrientation = (nodeOrder[0] == nodeOrder[1])

            # both tris already assigned to a surface
            if connected[0] is not None:
                if connected[0]==connected[1]:
                    # both tris belong to the same surface already
                    # if different orientation then this is a twisted /
                    # one-sided surface
                    if not(sameOrientation):
                        flipEdges[connected[0]].append(edge)
                else:  # if connected[0]!=connected[1]:
                    # tris on different surfaces so far
                    if not(sameOrientation):
                        # flip one entire surface, choose the smallest
                        surfId = min((len(surfaces[i]), i)
                                     for i in connected)[1]
                        for el in surfaces[surfId]:
                            self.flipTri(el, (0,1))

                    # merge the two surfaces
                    surfIds = sorted((len(surfaces[i]), i) for i in connected)
                    surfIds = [j for i,j in surfIds]  # smaller surface first
                    smallSurf = surfaces[surfIds[0]]
                    surfaces[surfIds[1]].extend(smallSurf)
                    triSurf.update( (el, surfIds[1]) for el in smallSurf)
                    del surfaces[surfIds[0]]
                    if surfIds[0] in flipEdges:
                        flipEdges[surfIds[1]].extend(flipEdges[surfIds[0]])
                        del flipEdges[surfIds[0]]

                # finished with both tris already assigned to a surf before
                continue

            # if no tri assigned to a surface yet then initialize new surface
            elif connected[1] is None:
                surfaces[nextSurfId] = list(twoTris)
                triSurf[twoTris[0]] = nextSurfId
                triSurf[twoTris[1]] = nextSurfId
                nextSurfId += 1

            # if one tri assigned to a surface then add the other
            else:
                surfaces[connected[1]].append[twoTris[0]]
                triSurf[twoTris[0]] = connected[1]

            # possibly flip (one) newly assigned tri
            if not(sameOrientation):
                # flip first tri
                self.flipTri(twoTris[0], nodeIndexes[0])

        # prepare result list with flipEdges-attribute
        surfaceList = list()
        for i in surfaces:
            surf = TriangleSurface(
                nodeCoords=self.nodeCoords,
                elNodes=dict(
                    (elem, self.elNodes[elem])
                    for elem in surfaces[i]))
            surf.flipEdges = flipEdges[i]
            surfaceList.append(surf)
        return surfaceList

    def simplyConnectedParts(self, splitangle=None):
        """Creates a list of parts of the surface that have all triangles
        connected. Edges with more than two connected elements or with a
        sharp angle are ignored, i.e. those edges don't constitute an tri to
        tri connection within this surface.

        @param splitangle: Stop at edges with an angle greater than
        splitangle (in degrees)

        @return: list of sets of element numbers. Each set represents a single
        part.

        @Note: You might want to L{filterDegeneratedTris}() before calling this
        method.
        """

        # collect T-junctions as splitting edges
        multiTriEdges = self.findNonManyfoldEdges()
        splitEdges = set(multiTriEdges)
        msg("Found %d non-manifold edges. They are not considered as"
            " connection." % len(multiTriEdges), debugLevel=1)

        # collect sharp edges as splitting edges
        if splitangle is not None:
            sharpEdges = self.findSharpEdges(splitangle)
            splitEdges.update(sharpEdges)
            msg("Found %d edges with angle greater than %s deg. They are not"
                " considered as connection."
                % (len(sharpEdges),splitangle), debugLevel=1)

        # collect single tris to exclude from the list of potentially connected
        # tri elements
        singleTris = self.findSingleTris()
        allElems = set(self.elNodes)
        allElems.difference_update(singleTris)
        msg("Found %d unconnected surface parts with only one triangle (single"
            " tris)." % len(singleTris), debugLevel=1)

        # find connected parts
        surflist = self.findConnectedParts(
            elems=allElems, splitEdges=splitEdges)
        msg("Found %d unconnected surface parts with more than one triangle."
            % len(surflist), debugLevel=1)

        # append single tris to the list
        surflist.extend( set((elem,)) for elem in singleTris )
        msg("TriangleSurface.simplyConnectedParts found %d individual surface"
            " parts." % len(surflist), debugLevel=1)

        return surflist

    @staticmethod
    def _joinPartsAngleWeightingFunction(angleCos):
        """Not meant to be called from "outside"
        angle weighting function used by getSimpleSurfacesFromMeshTris():
        at alpha = 0 => weighting function f = 1.0
        at alpha = 1 => weighting function f = 0.0
        => f = sum_i a_i cos^i alpha
        a_i = aList
        """
        # alpha max = 60 deg
        # cosMaxAngle = cos(60.0 *pi/180.0)
        # aList=[-3.52478665,24.41584112,-61.15783254,64.43440953,-23.16763147]

        # alpha max = 15 deg
        cosMaxAngle = cos(15.0 *pi/180.0)
        aList = [1673924.77719293, -6884242.83380361, 10617036.49685724,
                 -7277187.95816989, 1870470.51792333]

        if angleCos<cosMaxAngle:
            return 0.0
        else:
            return sum([a_i * angleCos**i
                        for i, a_i in enumerate(aList)])

    def splitForSurfToContact(self, splitangle=None, joinPartsAngle=15.0):
        """Split surface into separate parts with unified normals.

        Splitting will be done where there is no edge-to-edge connection,
        where the angle between adjacent tris is larger than splitangle and
        at non-manifold edges.

        Calls self.unifyNormals() and will therefore create
        self.flipEdges-attribute if self is not a double-sided surface (like
        the Moebius-strip)

        Tries to join parts over adges with an tri-to-tri-angle of less than
        joinPartsAngle.

        @param splitangle: split at edges with an angle larger than this.

        @param joinPartsAngle: Only join parts over an edge with a smaller
        angle than this. (In deg)

        @Returns: a list of TriangleSurface objects.
        """

        surfParts = self.simplyConnectedParts(splitangle)
        surflist = [
            TriangleSurface(
                nodeCoords=self.nodeCoords,
                elNodes=dict((el, self.elNodes[el])
                             for el in elems))
            for elems in surfParts]
        if len(surflist)==1:
            msg("Surface has not been split into pieces (nonmanifold and"
                " sharp edges were not considerd as connections).")
        else:
            msg("Split the surface into %d surface parts (nonmanifold and"
                " sharp edges were not considerd as connections)."
                % len(surflist))

        edgeToTri = self.getEdgeToTri()
        multiTriEdges = self.findNonManyfoldEdges()
        mteTris = set(elem for edge in multiTriEdges
                      for elem in edgeToTri[edge])
        mteTriPartIdx = dict(
            (elem, i)
            for i, elems in enumerate(surfParts)
            for elem in elems
            if elem in mteTris)
        triNormal = dict(
            (elem, getTriNormalSc(self.nodeCoords, self.elNodes[elem]))
            for elem in mteTris)

        # would like to but cannot delete mteTris: the dicts get generator
        # expressions that might evaluate the expression containing the
        # variables at a later time, after mteTris has already been deleted...
        # del mteTris, surfParts

        ticker = MsgTicker("Unifying normals for surface part %%d/%d"
                           % len(surflist))
        for surf in surflist:
            ticker.tick()
            surf.unifyNormals()
        del ticker
        msg("Finished unifying normals for surface parts.")

        ## now try to join parts together again
        msg("Try to reconnect parts at edges with more than two triangles."
            " There are %d such edges. Minimum join angle: %g."
            % (len(multiTriEdges), joinPartsAngle), debugLevel=1)

        joinPartsAngleCos = cos(joinPartsAngle*pi/180)

        ## possibleConnections is a dict { connection key : value list }
        # connection key is a set of two part indexes in connectedTris
        # value list components:
        #  . sum of edge connection values = values of the angle weighting
        #    function _joinPartsAngleWeightingFunction()
        #  . list of contradicting connections
        #  . min number of tris in parts
        #  . twist normals factor: -1 when the normals have to be twisted and
        #    1 if not
        possibleConnections = dict()
        for edgeCnt, thisEdge in enumerate(multiTriEdges):

            # ignore all tris on this edge that belong to a part with more than
            # one tris on this edge
            # {connected tri : part index in connectedTris}
            connTriToPart = dict([
                (tri, mteTriPartIdx[tri])
                for tri in edgeToTri[thisEdge]])
            connPartToNbTri = defaultdict(int)  # {part index: }
            for pn in connTriToPart.values():
                connPartToNbTri[pn] += 1

            # list of tris connected to this edge
            connectedTris = [tri for tri in edgeToTri[thisEdge]
                             if connPartToNbTri[connTriToPart[tri]]==1]
            if len(connectedTris)<2:
                continue

            # determine the orientation of the tris with respect to the
            # connecting edge
            orderedEdge = tuple(thisEdge)
            orderedEdgeTwisted = tuple(orderedEdge[::-1])
            triTwistNormalsFactor = list()
            for tri in connectedTris:
                triOrderedEdges = set([
                    ee[1] for ee in self.triEdgesOrderedNodesIter(
                        self.elNodes[tri])])
                if orderedEdge in triOrderedEdges:
                    triTwistNormalsFactor.append(1.0)
                elif orderedEdgeTwisted in triOrderedEdges:
                    triTwistNormalsFactor.append(-1.0)
                else:
                    raise Exception(
                        "Multiple connected edge %s not on tri %d with nodes"
                        " %s anymore."
                        % (orderedEdge, tri, self.elNodes[tri]))

            # all possible part connections
            for i1,t1 in enumerate(connectedTris):
                part1 = mteTriPartIdx[t1]
                normal1 = vector_scale(triNormal[t1],
                                       triTwistNormalsFactor[i1])
                for i2 in range(i1+1, len(connectedTris)):
                    t2 = connectedTris[i2]
                    part2 = mteTriPartIdx[t2]
                    twistNormalsFactor2 = -triTwistNormalsFactor[i2]
                    normal2 = vector_scale(triNormal[t2],
                                           twistNormalsFactor2)
                    connTwistNormalsFactor = (triTwistNormalsFactor[i1]
                                              * twistNormalsFactor2)

                    # ignore all connections with an angle greater than
                    # _joinPartsAngleDeg
                    t1t2AngleCos = dot(normal1, normal2)
                    if t1t2AngleCos<joinPartsAngleCos:
                        continue

                    # get or initialize the value list for this connection
                    thisConn = frozenset((part1, part2))
                    if thisConn in possibleConnections:
                        thisConnVal = possibleConnections[thisConn]
                    else:
                        thisConnVal = [0.0, set(), 0, connTwistNormalsFactor]
                        possibleConnections[thisConn] = thisConnVal

                    # Check 4th component of connection value list:
                    # twist normals factor: -1 when the normals have to be
                    # twisted and 1 if not. 0 if found to be contradicting.
                    # For such a connection we don't need further information
                    if thisConnVal[3]*connTwistNormalsFactor < 0.5:
                        thisConnVal[3] = 0
                        continue

                    # Update 3rd component of connection value list:
                    # min number of tris in parts
                    thisConnVal[2] = min(thisConnVal[2],
                                         len(surflist[i1].elNodes),
                                         len(surflist[i2].elNodes))

                    # Update 1st component of connection value list:
                    # sum of edge connection values = values of the angle
                    # weighting function
                    thisConnVal[0] += (
                        self._joinPartsAngleWeightingFunction(t1t2AngleCos))

                    # Update 2nd component of connection value list:
                    # list of contradicting connections
                    otherConnTris = set(connectedTris)
                    otherConnTris.remove(t1)
                    otherConnTris.remove(t2)
                    for t3 in otherConnTris:
                        thisConnVal[1].add(
                            frozenset((connTriToPart[t1], connTriToPart[t3])))
                        thisConnVal[1].add(
                            frozenset((connTriToPart[t2], connTriToPart[t3])))
        ## finished possibleConnections dict
        msg("Found %d possible part connections in the first place."
            % len(possibleConnections))

        # remove impossible connections
        # 4th component: twist normals factor == 0 if found to be contradicting
        for key in possibleConnections.keys():
            if possibleConnections[key][3] < 0.5:
                del possibleConnections[key]

        # calculate rating, 1st component (sum of edge connection values)
        # rating = rating of self*nbPossibleConnects minus sum of contradicting
        #          ratings
        nbPossibleConnects = len(possibleConnections)
        ratingConnect = list()
        for key, connVal in possibleConnections.iteritems():
            connVal[1].intersection_update(possibleConnections)
            rating = connVal[0] * (nbPossibleConnects-1)
            for contraKey in connVal[1]:
                rating -= possibleConnections[contraKey][0]
            ratingConnect.append((rating,key))

        # delete contradicting connections from top to bottom of rating list
        ratingConnect.sort()
        for rating, key in ratingConnect:
            try:
                for contraKey in possibleConnections[key][1]:
                    del possibleConnections[contraKey]
            except KeyError:
                # this connection has already been deleted
                pass

        # sort possibleConnections, result in list joinConnections:
        # transform frozenset keys to tuples (high number, low number)
        joinConnections = [list(key) for key in possibleConnections]
        for partList in joinConnections:
            partList.sort(reverse=True)
        joinConnections.sort(reverse=True)

        # Remove duplicates, each surface is at most joined into one other
        # (the second, target surface having a lower number).
        # Because joinConnections is sorted, if a surface could be joined with
        # several others, the lowest surface number is chosen.
        joinConnections = dict(joinConnections)
        joinConnections = list(joinConnections.iteritems())
        joinConnections.sort(reverse=True)

        # join the surfaces
        deleteParts = set()
        for partList in joinConnections:
            deleteParts.add(partList[0])

            msg("Merging surface part %d into part %d."
                % (partList[0]+1, partList[1]+1))

            # merging partList[0] into partList[1]
            surf = surflist[partList[1]]
            tomerge = surflist[partList[0]]
            # need only update elNodes and possibly flipEdges
            # nodeCoords is identical for all surfaces!
            surf.elNodes.update(tomerge.elNodes)
            if hasattr(tomerge, "flipEdges"):
                try:
                    surf.flipEdges.update(tomerge.flipEdges)
                except AttributeError:
                    # surf.flipEdges does not exist yet
                    surf.flipEdges = tomerge.flipEdges
            del surf, tomerge

        # delete the merged in parts
        surflist = [surf
                    for i, surf in enumerate(surflist)
                    if i not in deleteParts]

        msg("Joined parts of the surface %d times, leaving a total of %d"
            " parts." % (len(joinConnections), len(surflist)))

        return surflist

    #{ end of split surface

    ######################################################
    #{ file I/O
    def exportAsRaw(self, outputFile, objectName="Object1", numberFormat="%g"):
        """Export the triangle surface to the specified file in .raw format.
        This is suitable for import as mesh into Rhino.

        @param outputFile: If outputFile has got a write method (i.e. if
        outputFile is a file object or something similar) the output is being
        written to that file and the file is left open. Otherwise outputFile
        is treated as a file name and a corresponding file is being opened
        for write and closed after this procedure. If a file with that name
        existed before it is silently overwritten.
        @param objectName: Each mesh object in the raw file may be given a
        name that will become the mesh objects name in Rhino. This parameter
        accepts a string for this name.
        @param numberFormat: You may specify a format specifier that'll be used
        for each number (float) in the resulting raw file output. This can be
        used to increase the precision. E.g. "%1.3f" will print three decimal
        places.
        """
        try:
            writemethod = outputFile.write
        except AttributeError:
            output = open(outputFile, "w")
            writemethod = output.write

        writemethod("%s\n" % objectName)
        elems = self.elNodes.keys()
        elems.sort()
        for el in elems:
            nodes = self.elNodes[el]
            coords = [self.nodeCoords[node][ix]
                      for node in nodes
                      for ix in xrange(3)]
            line = " ".join(numberFormat%x for x in coords)
            writemethod("%s\n" % line)

        try:
            output.close()
        except NameError:
            pass
        return

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
    def __str__(self):
        """Give a short description of the surface.

        >>> surf = TriangleSurface.fromSTL("example.stl")
        >>> print surf
        TriangleSurface with 214 nodes, 201 triangles, ...
        ... min: [0.0,0.4,1.2], max: [5.3,8.4,3.2]

        """
        return ("TriangleSurface with %d nodes, %d triangles, %s"
                % (len(self.getConnectedNodes()),
                   len(self.elNodes), self.getBoundingBox()))

    def getBoundingBox(self):
        """
        Returns the bounding box of the surface

        Example:
          >>> # bounding box of that surface
          >>> bb = surf.getBoundingBox()
          >>> print "lower left front corner:", bb[0]
          >>> print "upper right back corner:", bb[1]

        @Returns: a BoundingBox object which is basically list of two
        coordinate tuples (rather lists).
        The first coord tuple states the min value, the last the max.
        I.e. self.getBoundingBox()[i][j]:
        i \in {0,1} for min/max value, j \in {0,1,2} coordinate index
        """

        # nodes set (contains unique nodes)
        nodes = self.getConnectedNodes()

        # update bounding box
        boundingBox = BoundingBox()
        boundingBox.update(nodeCoords=self.nodeCoords, elNodes=nodes)
        return boundingBox

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

    def getCentroid(self):
        """Returns the area-centroid of the surface.

        See also: L{getArea} to get the area of the surface.
        And L{getTriAreaNormal} for the area of individual faces.
        """
        # (centroid, area)-tuple for each triangle
        centroidAreaList = [
            # triangle centroid:
            (vector_scale(vector_sum(
                *(self.nodeCoords[n] for n in nodes)), 1.0/3),
             # triangle area:
             0.5*length(self.getTriNormalUsc(nodes)))
            # ... for each triangle
            for nodes in self.elNodes.itervalues()]

        # sum up centroids weighted by area and divide by accumulated area
        centroid = vector_scale(
            vector_sum(*(vector_scale(c,a) for c,a in centroidAreaList)),
            1.0 / sum(a for c,a in centroidAreaList))

        return centroid

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

    ######################################################
    #{ methods for connectivity and consistency info
    def getEdgeToTri(self):
        """Returns a dict {edge: set of tri-element-ids}. And edge is a
        frozenset of the node numbers of the two attached nodes.

        @Note: Currently this method stores a copy of its result in case it is
        called again. But this is an implementation detail that might change
        in future versions. Never rely on the stored copy directly. Instead
        just call this function repeatedly if required.
        """
        if self.edgeToTri is None:
            self.edgeToTri = defaultdict(set)
            for element, nodesOfThisTri in self.elNodes.iteritems():
                for thisEdge in self.triEdgesIter(nodesOfThisTri):
                    self.edgeToTri[thisEdge].add(element)
        return self.edgeToTri

    def getNodesToElem(self):
        """Return two dicts nodesToElem and nodesToDoubles.

        nodesToElem is a dict {frozenset(three node numbers) : element nb}

        If there are double tris they are placed in nodesToDoubles
        {frozenset(node numbers) : list of element nbs}
        """
        nodeToElem = dict()
        nodesToDoubles = defaultdict(list)
        for element, nodesOfThisTri in self.elNodes.iteritems():
            face = frozenset(nodesOfThisTri)
            if face in self.nodeToElem:
                nodesToDoubles[face].append(element)
            else:
                nodeToElem[face] = element
        return nodeToElem, nodesToDoubles

    def getNodeElems(self):
        """returns a dict {node id: set of attached element numbers}

        If you need a list of all connected nodes use L{getConnectedNodes}
        """
        nodeTris = defaultdict(set)
        for element, nodes in self.elNodes.iteritems():
            for node in nodes:
                nodeTris[node].add(element)
        nodeTris.default_factory = None  # finalize nodeTris
        return nodeTris

    def findNonManyfoldEdges(self):
        """
        Check if the surface triangles are simply connected. I.e. if all edges
        have exactly one or two triangles connected to it.

        If an edges has more than two triangles attached to it then this is a
        non-manifold edge or a 'surface junction'.

        Usage showing how to also get all nodes at those edges and connected
        tri elements:
         >>> nonManyfoldEdges = surf.findNonManyfoldEdges()
         >>> edgeToTri = surf.getEdgeToTri()
         >>> wrongNodes = set()
         >>> wrongTris = set()
         >>> for edge in nonManyfoldEdges:
         >>>     wrongTris.update(edgeToTri[edge])
         >>>     wrongNodes.update(edge)

        @returns: non-manifold edges. An edge is a frozenset of two node
        numbers.
        """
        nonManyfoldEdges = [
            edge
            for edge, tris in self.getEdgeToTri().iteritems()
            if len(tris) > 2]
        return nonManyfoldEdges

    def findBorderEdges(self):
        """Returns a set of edges ( =set(node1,node2) ) with only one tri
        connected to any other of this surface
        """
        edgeToTri = self.getEdgeToTri()
        borderEdges = set([edge
                           for edge,tris in edgeToTri.iteritems()
                           if len(tris)==1])
        return borderEdges

    def findBorderNodes(self):
        """Return the set of all nodes anywhere on the border.
        """
        borderNodes = set()
        for edge in self.findBorderEdges():
            borderNodes.update(edge)
        return borderNodes

    def findSingleTris(self):
        """Returns a list of tri elements without neighbours
        """

        edgeToTri = self.getEdgeToTri()
        singleTris = set(
            elem for elem, nodes in self.elNodes.iteritems()
            if all( len(edgeToTri[edge])==1
                    for edge in self.triEdgesIter(nodes)) )
        return singleTris

    def findConnectedParts(self, elems=None, splitEdges=None):
        """
        @param elems: Optional iterable of element numbers. If given then only
        consider these elements in the search for connected parts.

        @param splitEdges: an iterable of edges to be ignored when traversing
        connected elements.

        @Returns: A list of element-lists each list containing only
        elements that are connected over an edge.
        """
        if splitEdges:
            splitEdges = set(splitEdges)
        else:
            splitEdges = set()

        if elems is None:
            unassignedTris = set(self.elNodes)
        else:
            # make it a set and make a copy, don't modify the argument
            unassignedTris = set(elems)

        surflist = list()
        edgeToTri = self.getEdgeToTri()

        while len(unassignedTris)>0:

            # get any one element as a starting point
            startTriId = iter(unassignedTris).next()

            # collect tris of this surface-part in connectedTris
            connectedTris = set()
            connectedTris.add(startTriId)

            # trisToFindNeighbours is a list of (TriId, ElNodes, Normal)-tuples
            # of those tris for which the neighbours have to be examined yet.
            # The normals of the tris in this list have already been
            # calculated. the ElNodes-item contains the node numbers in
            # counter-clockwise succession if you look from "above" (pos side)
            trisToFindNeighbours = deque()  # a double-ended queue
            trisToFindNeighbours.append(startTriId)

            # main loop traversing all attached tris
            while len(trisToFindNeighbours) > 0:
                # get the first (oldest) tri in the list
                startTriId = trisToFindNeighbours.popleft()

                for thisEdge, edgeNodes in self.triEdgesOrderedNodesIter(
                        self.elNodes[startTriId]):

                    if thisEdge in splitEdges:
                        continue

                    trisOnThisEdge = edgeToTri[thisEdge].difference(
                        connectedTris)
                    for newTri in trisOnThisEdge:
                        connectedTris.add(newTri)
                        trisToFindNeighbours.append(newTri)

            # add to list
            surflist.append(connectedTris)

            # remove new surface from unassignedTris
            unassignedTris.difference_update(connectedTris)
            continue  # while len(unassignedTris)>0

        return surflist

    def findSharpEdges(self, minAngle, edgeToTri=None, triNormal=None):
        """
        @param minAngle: Find edges with an angle greater than minAngle
        (in degrees)

        @param edgeToTri: use this {edge: set of attached tris}-dict, otherwise
        get one from self.getEdgeToTri()

        @param triNormal: use this {element number: unit normal vector}-dict,
        otherwise initialize an empty dict. Note that this dictionary will be
        updated with normals as required by the algorithm.

        @returns: list of sharp edges. An edge is a frozenset of two node
        numbers.

        @Note: Edges with more or less than two attached tris are silently
        ignored.
        """
        if edgeToTri is None:
            edgeToTri = self.getEdgeToTri()
        if triNormal is None:
            triNormal = dict()
        anglecosThreshold = cos(float(minAngle)*pi/180.)

        def getNormal(elem, triNormal):
            try:
                normal = triNormal[elem]
            except KeyError:
                normal = getTriNormalSc(self.nodeCoords, self.elNodes[elem])
                triNormal[elem] = normal
            return normal

        sharpEdges = list()
        for edge, attachedTris in edgeToTri.iteritems():
            if len(attachedTris)!=2:
                continue

            anglecos = dot(*[getNormal(elem, triNormal)
                             for elem in attachedTris])
            if not(self.equallyOrientedTris(attachedTris, edge)):
                anglecos *= -1

            if anglecos < anglecosThreshold:
                sharpEdges.append(edge)

        return sharpEdges

    def equallyOrientedTris(self, twoTris, edge):
        """Check whether the two given tri elements attached to the given
        edge are oriented consistently.

        @Note: the two tris must be connected over the given edge.
        """
        # nodeIndexes: ...
        # [ [tri_1 node id 1, tri_1 node id 2], [tri_2 n-id 1, t2 n-id 2] ]
        # tri_1 is the possibly unassigned first tri, tri_2 is the other
        # nodeIndexes[i][1]-nodeIndexes[i][0] can be 1 or -2 meaning the
        # tri is "upright" on the edge. Or it can be -1 or 2 meaning the
        # tri is "upside-down" on the edge.
        # If both tris are "upright" or both are "upside-down" then the
        # orientation of both tris is consistent. Otherwise one tri must be
        # flipped.
        nodeIndexes = [
            [self.elNodes[tri].index(edgeNode) for edgeNode in edge]
            for tri in twoTris]
        nodeOrder = [(idx[1]-idx[0]) % 3 for idx in nodeIndexes]
        sameOrientation = (nodeOrder[0] != nodeOrder[1])
        return sameOrientation

    #} end of methods for connectivity and consistency info

    ######################################################
    #{ edit / modify

    def insertSurface(self, otherSurface, mergeNodes={}):
        """Add the contents of another surface to self.

        This has been taken from abq_model_02.Model.insertModel().

        Node and element number of the otherSurface will be conserved if there
        is no node resp. element in self falling in between the those new
        numbers. I.e. if there is node 1, 2, 10, 11 in self and nodes 3, 4, 5
        in otherSurface then they will retain their numbers upon import to self
        and will be added to self as nodes 3, 4, 5. If there is again node 1,
        2, 10, 11 in self and now nodes with numbers 7, 12, 17 in otherSurface
        then they will get new numbers and become the new nodes 12, 13, 14 in
        self. (The same applies to element numbers accordingly.)

        Nodes that are not defined in otherSurface.nodeCoords but used in the
        element definitions otherSurface.elNodes will keep their old numbers.
        This is a feature to enable import of elements connected to already
        defined nodes with this function. You might want to check beforehand
        if all nodes are well defined if you don't need this feature.

        Node coordinates are deep copied, i.e. if you later change individual
        node coordinates in otherSurface then the node coordinates of self are
        not affected. Same with the element connectivity elNodes.

        @param otherSurface: Surface-instance to be added to self.

        @param mergeNodes: A dictionary {nodeOther : nodeThis}. nodeOther
          is the node number in otherSurface of the node to be replaced by
          the nodeThis-node in self.

          I.e. say mergeNodes={1000:10, ...}. Node 10 is a node in self. All
          occurences of node 1000 in otherSurface are replaced by 10 during the
          process of adding otherSurface to self. The node 1000 from
          otherSurface won't be imported to self.

        @Return: A tuple (nodesOldToNew, elemsOldToNew):
          - nodesOldToNew: {node number in otherSurface, new node num}
          - elemsOldToNew: {elNum in otherSurface, new elNum}
        """

        otherNodesToAdd = sorted( node for node in otherSurface.nodeCoords
                                  if node not in mergeNodes )

        # determine new node numbers
        if len(otherNodesToAdd)==0:
            newNodes = []
        elif set(self.nodeCoords).intersection(
                range(otherNodesToAdd[0], otherNodesToAdd[-1]+1)):
            # conflicting node number ranges
            firstNode = max(self.nodeCoords)+1
            newNodes = range(firstNode, firstNode+len(otherNodesToAdd))
        else:
            newNodes = otherNodesToAdd

        # insert nodes
        self.nodeCoords.update(
            (newN, otherSurface.nodeCoords[oldN][:])
            for newN, oldN in izip(newNodes, otherNodesToAdd)
            )

        # initialize nodesOldToNew: {node number in otherSurface, new node num}
        nodesOldToNew = dict(mergeNodes)
        nodesOldToNew.update(izip(otherNodesToAdd, newNodes))
        del otherNodesToAdd, newNodes

        # determine new element numbers
        oldElems = sorted(otherSurface.elNodes)
        if len(oldElems) and set(self.elNodes).intersection(
                range(oldElems[0], oldElems[-1]+1)):
            # conflicting element number ranges
            firstElem = max(self.elNodes)+1
            newElems = range(firstElem, firstElem+len(oldElems))
        else:
            newElems = oldElems

        # initialize elemsOldToNew: {elNum in otherSurface, new elNum}
        elemsOldToNew = dict(izip(oldElems, newElems))
        del oldElems, newElems

        # insert elements
        elNodesDict = dict(
            (elemsOldToNew[elem],
             [nodesOldToNew.get(node, node) for node in nodes])
            for elem, nodes in otherSurface.elNodes.iteritems())
        # update Surface object attributes
        dict.update(self.elNodes, elNodesDict)

        # The following attributes are reset to None.
        # They are supposed to be None as long as they are not calculated yet
        # don't query these variables, use the corresponding get...() method,
        # ie. for self.edgeToTri use self.getEdgeToTri()
        # This operation renders this data invalid / incomplete therefore
        # the corresponding variable (eg. self.edgeToTri) is reset to None.
        # It will then be recalculated automatically
        self.edgeToTri = None
        self.pointSearchData = [False, False, False]

        # return some of the transfer parameters
        return (nodesOldToNew, elemsOldToNew)

    def removeElems(self, elems, removeNodes=False):
        """remove tris from the surface

        update self.elNodes, self.edgeToTri

        @param elems: A list or a set of element numbers (also accepts a dict,
        in that case takes its keys.)

        @param removeNodes: if True remove unused nodes. This is not very
        efficient if many elements are beeing deleted, in this case rather
        clean up after all elements have been removed using
        getNodesFromElset():
          >>> surf.removeElems(elems=myBigElset)
          >>> allNodes = surf.getNodesFromElset()
          >>> removeNodes = set(surf.nodeCoords).difference(allNodes)
          >>> for node in removeNodes:
          ...    del surf.nodeCoords[node]

        @note: Silently ignores missing elements.
        @note: It is save to do the following (it is not only save but the
        recommended way to accomplish those tasks):
          >>> surf.removeElems(elems=surf.elNodes)  # remove all elements
        """

        # special treatment if you want to remove all, wouldn't work otherwise
        if elems is self.elNodes:
            self.elNodes.clear()
            self.edgeToTri = None
            if removeNodes:
                self.nodeCoords.clear()
            return

        # update self.edgeToTri
        if self.edgeToTri is not None:
            removedEdges = list()
            for elNum in elems:
                try:
                    nodes = self.elNodes[elNum]
                except KeyError:
                    continue
                for thisEdge in self.triEdgesIter(nodes):
                    self.edgeToTri[thisEdge].remove(elNum)
                    removedEdges.append(thisEdge)
                    if len(self.edgeToTri[thisEdge]) == 0:
                        del self.edgeToTri[thisEdge]

        # node list of element nodes
        if removeNodes:
            removeNodesSet = set()
            for elNum in elems:
                try:
                    removeNodesSet.update(self.elNodes[elNum])
                except KeyError:
                    # elNum not in self.elNodes
                    pass

        # remove elements from self.elNodes
        for elNum in elems:
            try:
                del self.elNodes[elNum]
            except KeyError:
                pass

        # remove nodes
        if removeNodes:
            for nodes in self.elNodes.itervalues():
                removeNodesSet.difference_update(nodes)
                if len(removeNodesSet)==0:
                    break
            for node in removeNodesSet:
                del self.nodeCoords[node]

        # point search has to be initialized again
        self.pointSearchData = [None, None, None]

    def addSingleTri(self, nodes, elNum=None):
        """Inserts a new tri into the surface attached to already existing
        nodes.

        Swap node ordering in order to have its normal consistent with its
        neighbours.

        update self.elNodes, self.edgeToTri

        @param nodes: node number list of the new element
        """
        if elNum is None:
            elNum = max(self.elNodes)+1

        edgeToTri = self.getEdgeToTri()

        # update edgeToTri and
        # find a neighbour to determine normal
        # ... only if edgeToTri already defined
        neighbour = None
        for edge in self.triEdgesIter(nodes):
            try:
                newneighbours = edgeToTri[edge]

            except KeyError:
                # this edges is not yet in edgeToTri, its a new edge
                # a new edge is an edge on the border
                pass

            else:
                if len(newneighbours)>1:
                    raise SurfaceConnectionError(
                        "There are already more than one elements"
                        " connected to the edge %s when trying to insert"
                        " new tri element %d on nodes %s."
                        % (str(edge), elNum, str(nodes)))

                # if we don't have a neighbour yet, but need one...
                if neighbour is None:
                    neighbour = iter(newneighbours).next()  # get any

            # anyway update edgeToTri
            edgeToTri[edge].add(elNum)

        if (not neighbour):
            raise SurfaceConnectionError(
                "New tri element %d (nodes: %s) has no neighbours."
                % (elNum, str(nodes)))
        # sort the nodes in an order such that its normal fits to
        # the neighbour
        newNodes = [0,0,0]
        NodesUnsorted = set(nodes)
        for idx, node in enumerate(self.elNodes[neighbour]):
            if node in NodesUnsorted:
                newNodes[-idx] = node
                NodesUnsorted.remove(node)
            else:
                thirdIdx = -idx
        newNodes[thirdIdx] = NodesUnsorted.pop()

        # update elNodes
        self.elNodes[elNum] = newNodes

        # point search has to be initialized again
        self.pointSearchData = [None,None,None]

    def updateElems(self, elNodes):
        """Inserts new tris into the surface replacing those that already exist
        with the same element number.

        updates self.elNodes

        removes self.edgeToTri

        @param elNodes: dictionary {element number: node number
        list} of elements to be updated or inserted
        """

        dict.update(self.elNodes, elNodes)
        self.edgeToTri = None
        # point search has to be initialized again
        self.pointSearchData = [None,None,None]

    def replaceNode(self, triId, oldNode, newNode):
        """Change surface topology: Reconnect triangle to other node.

        Reconnect the corner of the triangle triId which was originally
        connected to oldNode to newNode. That is used when the surface is split
        by surfToContact.

        Update self.edgeToTri accordingly. Don't remove the oldNode and don't
        replace the node in other triangles of the surface.
        """
        try:
            idx = self.elNodes[triId].index(oldNode)
        except ValueError:
            # ... .index(oldNode) failed
            raise ValueError("Node %d not in tri element %d in this surface"
                             " when trying to replace it."
                             % (oldNode, triId))
        except IndexError:
            # self.elNodes[triId] failed
            raise IndexError("No tri element %d in this surface when trying"
                             " to replace node %d on it."
                             % (triId, oldNode))
        # update elNodes
        self.elNodes[triId][idx] = newNode

        # update edgeToTri
        if self.edgeToTri is not None:
            for edge,triList in self.edgeToTri.items():
                if (oldNode in edge) and (triId in triList):
                    triList.remove(triId)
                    newEdge = set(edge)
                    newEdge.remove(oldNode)
                    newEdge.add(newNode)
                    newEdge = frozenset(newEdge)
                    triList = self.edgeToTri[newEdge]
                    triList.add(triId)

        # point search has to be initialized again
        self.pointSearchData = [None,None,None]
    #{ end of edit / modify

    ######################################################
    #{ transform geometry
    def translate(self, move_vector):
        """Move the surface in space by move_vector

        Example:

        >>> # move this surface up by 10
        >>> surf.translate(move_vector=[0,0,10])

        Translating is done by simply moving the point coordinates for this
        surface. Other surfaces or the original mesh from which the surface
        originates are not affected.

        This creates a new node-dictionary, which is refered to from now on.
        It contains only those nodes connected to any element in self.elNodes.

        The point search initialization is modified and doesn't have to be done
        again.
        """

        all_nodes = self.getConnectedNodes()

        new_coords = dict()
        for node_id in all_nodes:
            new_coords[node_id] = vector_plus(
                self.nodeCoords[node_id],move_vector)

        self.nodeCoords = new_coords

        if len(self.elNodes)==0:
            return

        # adapt point search properties
        for axis, psData in enumerate(self.pointSearchData):
            if not psData:
                continue

            # axis index for the first and second in-plane axes
            # (defaults to 0, 1 for the default axis==2,
            # (or else [1,2] or [2,0] for axis==0 or 1)
            aix = [1, 2, 0, 1, 2][axis:axis+2]

            # tet Centroids
            for centroid in psData.triCentroids.itervalues():
                centroid[0] += move_vector[aix[0]]
                centroid[1] += move_vector[aix[1]]

            # BoundingBox
            psData.boundingBox[0][0] += move_vector[aix[0]]
            psData.boundingBox[0][1] += move_vector[aix[1]]
            psData.boundingBox[1][0] += move_vector[aix[0]]
            psData.boundingBox[1][1] += move_vector[aix[1]]

    def scale(self, scale_factor, scale_origin):
        """Scale the surface by a certain factor relative to scale_origin

        scale_factor may be a single number in which case it is applied to all
        three dimensions or a list of three numbers with one factor for each
        dimension.

        Example:

        >>> # scale this surface by 2 in each direction
        >>> surf.scale(scale_factor=2, scale_origin=[0,0,0])

        Scaling is done by simply moving the point coordinates for this
        surface. Other surfaces or the original mesh from which the surface
        originates are not affected.

        This creates a new node-dictionary, which is refered to from now on.
        It contains only those nodes connected to any element in self.elNodes.

        After scaling point search has to be initialized again.
        """

        all_nodes = self.getConnectedNodes()

        if isinstance(scale_factor, (float, int)):
            scale_factor = [scale_factor,scale_factor,scale_factor]

        new_coords = dict()
        for node_id in all_nodes:
            orig2old = vector(scale_origin, self.nodeCoords[node_id])
            orig2new = [orig2old[0]*scale_factor[0],
                        orig2old[1]*scale_factor[1],
                        orig2old[2]*scale_factor[2]]
            new_coords[node_id] = vector_plus(scale_origin, orig2new)

        self.nodeCoords = new_coords

        # after scaling point search has to be initialized again
        self.pointSearchData = [None,None,None]

    def offset(self, distance):
        """Move each node by the specified distance perpendicular to the
        surface into the direction of the node normal.

        For a closed surface (to check: assert not(self.L{findBorderEdges}()))
        you may want to call self.L{unifyNormals}() and then
        self.L{getVolume}() to check if the normals point outwards.
         >>> dist = 10.0  # move it 10m out
         >>> surf.unifyNormals()
         >>> if surf.getVolume()>0:
         >>>     dd = dist
         >>> else:
         >>>     dd = -dist
         >>> self.offset(dd)
        """
        nodeNormal = self.getNodeNormal()  # delivers unit normals
        for node in self.getConnectedNodes():
            vector_modif_add_scaled(
                self.nodeCoords[node], nodeNormal[node], distance)
        return

    #{ end of transform geometry

    ######################################################
    #{ edit / modify
    def cleanNodeCoords(self, checkOnly=False):
        """Replace self.nodeCoords with a copy of itself that only contains
        nodes that are attached to any triangle (as listed in self.elNodes)
        Note that it's a shallow copy: the individual lists of node
        coordinates are identical to those of the old nodeCoords-dict.

        @param checkOnly: Skip the clean operation and only report whether
        self.nodeCoords is "clean".

        @Returns: True if nodeCoords is clean (already), i.e. all nodes are
        connected to triangles and no change had to be performed on self.
        Returns False if nodeCoords contained nodes not connected to any
        triangle. Note that the result only reflects the situation *before*
        the possible clean operation.

        @Note: No check is performed whether all nodes in self.elNodes are
        actually defined in self.nodeCoords.
        """
        connectedNodes = set()
        for nodes in self.elNodes.itervalues():
            connectedNodes.update(nodes)
        clean = not(set(self.nodeCoords).difference(connectedNodes))
        if checkOnly:
            return clean
        self.nodeCoords = dict(
            (node, coords)
            for node, coords in self.nodeCoords.iteritems()
            if node in connectedNodes)
        return clean

    def flipTri(self, elNum, swapNodeIds):
        """Flip the node ordering of a single tri element.

        If you have already calculated a triNormal-dictionary like the one you
        get from self.L{getTriNormal}() it makes a lot of sense to update that:
         >>> vector_modif_scale(triNormal[elNum], -1)

        @param elNum: tri element to be flipped over
        @param swapNodeIds: tuple of the indexes of the nodes that shall be
        exchanged. E.g. (0,1) or (2,0)

        @Note: Contrary to the poor conception of the simple minded author of
        these lines this function does not affect the edgeToTri-structure that
        you get from self.L{getEdgeToTri}().
        """
        elNodes = self.elNodes[elNum]  # just an abbreviation
        i1, i2 = swapNodeIds
        swap = elNodes[i1]
        elNodes[i1] = elNodes[i2]
        elNodes[i2] = swap

    def unifyNormals(self, startElNum=None, skipMoebiusCheck=False):
        """Make the orientation of the tri elements consistent over connecting
        edges.

        Creates the set self.flipEdges if skipMoebiusCheck is False: If the
        surface does not have two sides (like the Moebius strip) then
        self.flipEdges contains edges at which the surface normals flip (after
        the unifying process).

        @param skipMoebiusCheck: If True no check is performed, whether the
        surface has two sides and the normals can be made consistent.

        @Returns: a list of all triangles that have been flipped.
        """

        if not(len(self.elNodes)):
            return list()

        edgeToTri = self.getEdgeToTri()

        flippedTris = list()
        trisToCheck = set(self.elNodes)

        if startElNum is None:
            startElNum = trisToCheck.pop()
        else:
            trisToCheck.remove(startElNum)

        # collect moebius-stop-edges in flipEdges
        if not skipMoebiusCheck:
            self.flipEdges = set()

        while len(trisToCheck)>0:

            # trisToFindNeighbours is a list of (elNum, skipEdge)-tuples for
            # each tri for which the neighbours have to be examined yet.
            # skipEdge is the edge over which the algorithm approached this tri
            # so this edge has already been processed before and can be
            # skipped. The tris in this list already have the correct
            # orientation:
            # node numbers are in counter-clockwise succession if you look from
            # "above" (positive side)
            trisToFindNeighbours = deque()  # a double-ended queue
            trisToFindNeighbours.append((startElNum, None))

            while len(trisToFindNeighbours) > 0:
                # get the first (oldest) tri in the list
                startElNum, skipEdge = trisToFindNeighbours.popleft()

                for thisEdge, edgeNodes in self.triEdgesOrderedNodesIter(
                        self.elNodes[startElNum]):

                    if thisEdge == skipEdge:
                        continue

                    trisOnThisEdge = edgeToTri[thisEdge]
                    if len(trisOnThisEdge)!=2:
                        # more than two tris on this edge: ignore this edge
                        # less than two: edge on border, no other tri connected
                        continue

                    # newTri := other tri on this edge, not startElNum
                    for newTri in trisOnThisEdge:
                        if newTri != startElNum: break

                    # if no check and normal already calculated, skip the rest
                    triCheckedEarlier = newTri not in trisToCheck
                    if skipMoebiusCheck and triCheckedEarlier:
                        continue

                    # append this neighbour to trisToFindNeighbours-list
                    if not triCheckedEarlier:
                        trisToFindNeighbours.append((newTri, thisEdge))
                        trisToCheck.remove(newTri)

                    # check orientation:
                    # indexes of edgeNodes in newTri.elNodes (i1, i2) must be
                    # ... (0,1) or (1,2) or (2,0)
                    newTriElNodes = self.elNodes[newTri]  # just an abbreviation
                    ids = [newTriElNodes.index(edgeNodes[i]) for i in 1,0]
                    if (ids[1]-ids[0]) % 3 == 1:
                        # newTri has right orientation, nothing to be done
                        continue

                    # orientation of newTri different
                    if triCheckedEarlier:
                        # here that implies: skipMoebiusCheck is False
                        # because of earlier "if skipMoebiusCheck and ..."
                        self.flipEdges.add(thisEdge)
                    else:
                        # reverse orientation
                        self.flipTri(newTri, ids)
                        flippedTris.append(newTri)

                    continue  # while len(trisToCheck)>0

            # prepare next iteration
            if len(trisToCheck)>0:
                startElNum = trisToCheck.pop()

        # end of the procedure
        return flippedTris

    def filterDegeneratedTris(self):
        """Remove tris with two ore more corners on the same node.

        This is done by counting the elements in set(self.elNodes[i]). So tri
        elements that don't have exactly three nodes in the first place are
        being removed as well.

        This does not create new holes if everything else is ok. It does not
        account for duplicate nodes: Nodes with different node numbers but
        same coordinates.

        @Returns: a list of removed tri element numbers.
        """
        remove = [el for el, nodes in self.elNodes.iteritems()
                  if len(set(nodes))!=3]
        for el in remove:
            del self.elNodes[el]
        return remove

    #} end of edit / modify

    ######################################################
    #{ methods for searching points in the x-y-projection

    # Constant for point search: cell Size for array relative to 'avgSize'
    cellSizeFactor = 0.5

    # Constant for point search:
    # in method pointIsBelow() pointSearchTol is compared to values 0..1
    pointSearchTol = 1E-4

    def initializePointSearch(self, smooth_findradius=0.0, axis=2):
        """
        Usually you don't need to call this function.

        You may call this function once before L{pointIsBelow} or
        L{projectPoints}, otherwise it is done automatically.

        smooth mode: Additionally consider points to be below the surface that
        are in fact outside but within the smooth_findradius to the
        closest centre point of any of the triangles of the surface (in the
        x-y-plane). Note that the argument smooth_findradius is deprecated.
        Instead set the instance attribute smoothFindRadius directly.

        @param smooth_findradius: DEPRECATED argument. Set the instance
           variable L{smoothFindRadius} instead.

        @param axis: normal axis index. axis=2 means vertical projection i.e.
           pointIsBelow actually checks if the vertical projection of a point
           is on the surface. axis=0/1 means pointIsBelow checks if the
           projection in x-/y-direction is on the surface.
        """

        # list of centroid-coord-tuples
        triCentroids = dict()

        # list of [[xmin,ymin],[xmax,ymax]] for each tri element
        elementBox = dict()

        boundingBox = BoundingBox(dim=2)
        avgSize = 0

        # axis index for the first and second in-plane axes
        # (defaults to 0, 1 for the default axis==2,
        # (or else [1,2] or [2,0] for axis==0 or 1)
        aix = [1, 2, 0, 1, 2][axis:axis+2]

        # size of the projection on the x-z plane
        # pos facing one direction, neg facing the other
        elSignedSize = dict()
        for element,nodes in self.elNodes.iteritems():

            node_pts = [self.nodeCoords[node]
                        for node in nodes]
            v12 = vector(node_pts[0], node_pts[1])
            v13 = vector(node_pts[0], node_pts[2])
            detA = v12[aix[0]]*v13[aix[1]] - v12[aix[1]]*v13[aix[0]]
            elSignedSize[element] = 0.5*detA

            # calculate 2D centroid
            x=y=0
            for coords in node_pts:
                x=x+coords[0]/3.0
                y=y+coords[1]/3.0
            triCentroids[element]=[x,y]

            # 2D bounding box
            xy = zip(*node_pts)
            elementBox[element]=[[min(xy[aix[0]]),min(xy[aix[1]])],
                                 [max(xy[aix[0]]),max(xy[aix[1]])]]

            # update bounding box of all elements
            try:
                for i in xrange(2):  # x,y
                    if elementBox[element][0][i]<boundingBox[0][i]:  # min
                        boundingBox[0][i]=elementBox[element][0][i]
                    if elementBox[element][1][i]>boundingBox[1][i]:  # max
                        boundingBox[1][i]=elementBox[element][1][i]
            except IndexError:
                boundingBox[0]=elementBox[element][0][:]  # min, make a copy
                boundingBox[1]=elementBox[element][1][:]  # max, make a copy

            avgSize+=self.dist2D(
                elementBox[element][1],elementBox[element][0])

        avgSize /= len(self.elNodes)

        msg('bounding box of the surface: %s' % boundingBox)
        msg('avg. element size: %g' % avgSize)


        ## -----------------------------------------------------------------
        ## map elements to cells

        cellSize = avgSize*self.cellSizeFactor

        # x,y,z number of cells
        cellDimensions = [boundingBox[1][0]-boundingBox[0][0],
                          boundingBox[1][1]-boundingBox[0][1]]
        msg('cell array dimensions: %s' % " ".join(map(str, cellDimensions)))

        cellNumber = [int(x/cellSize)+1 for x in cellDimensions]
        msg('cell numbers: %s cells: %d'
            % (" ".join(map(str, cellNumber)), cellNumber[0]*cellNumber[1]))

        # create empty array
        for i in range(cellNumber[0]):
            if i==0: cells=[]
            for j in range(cellNumber[1]):
                if j==0: LJ=[]
                LJ.append([])
            cells.append(LJ)

        # tris with a too small projected size will be ignored
        zerosize = (self.pointSearchTol*avgSize*0.1)**2

        # map elements
        for element in self.elNodes:

            # ignore too small tris
            if abs(elSignedSize[element]) < zerosize:
                continue

            # get index ranges for the cells that the box of the current
            # element spans
            cellOrigin = boundingBox[0]
            ebox = elementBox[element]  # abbreviation
            cellIndex = [
                # lower bound of cell indexes
                [int((ebox[0][0] - cellOrigin[0])/cellSize),
                 int((ebox[0][1] - cellOrigin[1])/cellSize)],
                # high bound of cell indexes
                [int((ebox[1][0] - cellOrigin[0])/cellSize),
                 int((ebox[1][1] - cellOrigin[1])/cellSize)]]

            # append element ID to cells that overlap the element's bounding box
            for i in range(cellIndex[0][0],cellIndex[1][0]+1):
                for j in range(cellIndex[0][1],cellIndex[1][1]+1):
                    cells[i][j].append(element)

        # permanently store objects as object attributes
        pointSearchData = Container()
        pointSearchData.triCentroids = triCentroids
        pointSearchData.boundingBox = boundingBox
        pointSearchData.cellNumber = cellNumber
        pointSearchData.cellSize = cellSize
        pointSearchData.cells = cells

        # for smooth mode store radius
        self.smoothFindRadius = float(smooth_findradius)

        # finished initializing
        self.pointSearchData[axis] = pointSearchData
        return

    def pointIsBelow(self, point, returnElems=False, axis=2):
        """
        Test if point is exactly below or above the surface in the sense of
        a vertical projection on the x-y plane (see axis-argument). The
        function does not decide whether a point is actually below or above
        the surface. The z coordinates are ignored completely, if axis is
        specified then those of the corresponding axis.

        @param point: a tuple of point coordinates

        @param returnElems: If True the function returns a dictionary
        {triNb:(a0, a1)}.
        The vertical projection of the point on the tri element triNb is
        pt_0 + a0*vec_01 + a1*vec_02, with pt_0 being
        self.nodeCoords[self.elNodes[triNb][0]], vec_01 being
        the vector from node self.elNodes[triNb][0] to self.elNodes[triNb][1]
        and vec_02 from self.elNodes[triNb][0] to self.elNodes[triNb][2]:

        >>> from bae.vecmath_01 import *
        >>> point = [...some coords...]
        >>> surface = TriangleSurface(....)
        >>> triA0A1 = surface.pointIsBelow(point, returnElems=True)
        ... for triNb, (a0,a1) in triA0A1.iteritems():
        ...     coords = [surface.nodeCoords[node]
        ...               for node in surface.elNodes[triNb]]:
        ...                  # newPt = pt_0 + a0*vec_01 + a1*vec_02
        ...                   newPt = list(coords[0])
        ...                   vector_modif_add(newPt, vector_scale(vec_01,a0))
        ...                   vector_modif_add(newPt, vector_scale(vec_02,a1))

        Elements that are only found by means of the "smooth mode" (see
        self.L{smoothFindRadius}) will have (a0, a1) = (0,0).

        @param axis: normal axis index. axis=2 means vertical projection i.e.
           pointIsBelow actually checks if the vertical projection of a point
           is on the surface. axis=0/1 means pointIsBelow checks if the
           projection in x-/y-direction is on the surface.

        @Note: On the first invocation of pointIsBelow() the search algorithm
        is initialized by automatically calling self.initializePointSearch().
        If you now change the coordinates or the connectivity of the surface
        subsequent calls to pointIsBelow() will not work properly. To force a
        new initialization based on the updated data either call
        self.initializePointSearch() manually or set
        self.pointSearchInitialized = False.

        @Bug: tris with a too small projected size will be ignored
        """
        if returnElems:
            notFound = dict()
        else:
            notFound = False

        if len(self.elNodes)==0:
            return notFound

        if not self.pointSearchData[axis]:
            self.initializePointSearch(axis=axis)
        psData = self.pointSearchData[axis]

        # axis index for the first and second in-plane axes
        # (defaults to 0, 1 for the default axis==2,
        # (or else [1,2] or [2,0] for axis==0 or 1)
        aix = [1, 2, 0, 1, 2][axis:axis+2]

        # cell check (includes bounding box test from cell index)
        cellOrigin = psData.boundingBox[0]
        idx0 = int((point[aix[0]]-cellOrigin[0])/psData.cellSize)
        if idx0<0 or idx0>=psData.cellNumber[0]:
            return notFound
        idx1 = int((point[aix[1]]-cellOrigin[1])/psData.cellSize)
        if idx1<0 or idx1>=psData.cellNumber[1]:
            return notFound

        # if cell check passed ... go on
        elems_in_cell = psData.cells[idx0][idx1]
        # sorted by distance
        centroid_to_point_distance = sorted(
            [self.dist2D(point, psData.triCentroids[element]), element]
            for element in elems_in_cell)

        # full check if point inside any tri of elems_in_cell
        inside = dict()

        if (self.smoothFindRadius and centroid_to_point_distance):
            for L, element in centroid_to_point_distance:
                if L >= self.smoothFindRadius:
                    break
                inside[element] = (0,0)
                if not returnElems:
                    break

        for [L, element] in centroid_to_point_distance:
            node_pts = [self.nodeCoords[node]
                        for node in self.elNodes[element]]
# Explanation of the following algorithm:

#       2
#      ^
#     /  .
#  b /     .
#   /        .
# 0 ----------> 1
#        a

# How to find out, if a point is inside the triangle 0 - 1 - 2:
# or, more precise, if the projection of that point is inside
# the triangle:

# Vectors starting on 0 end ending on the dotted line are:
#  v. = delta a + (1-delta) b              with 0 < delta < 1
# Vectors shorter than v. but pointing in the same direction are:
#  v = gamma v.                          with 0 < gamma < 1
#    = gamma ( delta a + (1 - delta) b)
#    = gamma delta a + gamma ( 1 - delta) b
#
# Define:
# alpha := gamma delta                     with 0 < alpha < 1
#  beta := gamma (1 - delta)               with 0 < beta  < 1
#
# So that we have:
#  v = alpha a + beta b

# Looking at which values alpha and beta can take:
# alpha + beta = gamma (delta + ( 1 - delta) ) = gamma < 1
# <=> alpha + beta < 1   =>  alpha < 1 AND beta < 1

# Algorithm:
# Given a point p, which could be inside the triangle 0 - 1 - 2, by
# its position vector x_p. Find alpha, beta by:

# 1) Compute v = x_p - x_0, a vector pointing from 0 to the point.
# 2) Compute alpha and beta from the first 2 equations of:
#    v = alpha a + beta b
# 3) Check if alpha + beta < 1 and if alpha > 0 and if beta > 0.

            # vectors
            v12=vector(node_pts[0], node_pts[1])
            v13=vector(node_pts[0], node_pts[2])
            detA = float(v12[aix[0]]*v13[aix[1]] - v12[aix[1]]*v13[aix[0]])
            adjA = [[v13[aix[1]], -v13[aix[0]]], [-v12[aix[1]], v12[aix[0]]]]
            invA = [[x/detA for x in row] for row in adjA]

            v1p=vector(node_pts[0], point)
            a0 = invA[0][0]*v1p[aix[0]] + invA[0][1]*v1p[aix[1]]
            a1 = invA[1][0]*v1p[aix[0]] + invA[1][1]*v1p[aix[1]]

            # if a0>=0 and a1>=0 and (a0+a1)<=1:
            #     inside = 1
            #     break
            if ((a0 > -self.pointSearchTol)
                and (a1 > -self.pointSearchTol)
                and ((a0+a1) < (1+self.pointSearchTol))):
                inside[element] = (a0, a1)
                if not returnElems:
                    break

        if returnElems:
            return inside
        else:
            return bool(len(inside))

    def projectPoints(self, points, axis=2):
        """projects arbitrary points vertically on the surface

        Returns a list of as many items as there are items in points. Each of
        those items is a list of projected points. It may be empty if the point
        is not below or above the surface or it may contain one ore many points
        on the surface.

        If two of the points are closer than self.pointSearchTol only the first
        is recognized.

        @param axis: normal axis index. axis=2 means vertical projection i.e.
           pointIsBelow actually checks if the vertical projection of a point
           is on the surface. axis=0/1 means pointIsBelow checks if the
           projection in x-/y-direction is on the surface.

        @Note: Uses self.pointIsBelow(). See its documentation in case the
        coordinates of the surface change between calls.
        """
        result = list()
        for point in points:
            triA0A1 = self.pointIsBelow(point, returnElems=True, axis=axis)
            projectedPts = list()
            for triNb, (a0,a1) in triA0A1.iteritems():
                coords = [self.nodeCoords[node]
                          for node in self.elNodes[triNb]]
                vec_01 = vector(coords[0], coords[1])
                vec_02 = vector(coords[0], coords[2])

                # newPt = pt_0 + a0*vec_01 + a1*vec_02
                newPt = list(coords[0])
                vector_modif_add(newPt, vector_scale(vec_01,a0))
                vector_modif_add(newPt, vector_scale(vec_02,a1))

                tolDist = self.pointSearchTol*(
                    0.5*(length(vec_01)+length(vec_02)))
                for oldPt in projectedPts:
                    if dist(newPt, oldPt)<tolDist:
                        newPt = None
                        break
                if newPt:
                    projectedPts.append(newPt)
            result.append(projectedPts)
        return result

    #} end of methods for searching points in the x-y-projection

    ######################################################
    #{ some service functions

    @staticmethod
    def dist2D(a,b):
        "2D distance only considering x and y coordinates"
        return sqrt((b[0]-a[0])**2 + (b[1]-a[1])**2)

    @staticmethod
    def triEdgesOrderedNodesIter(nodesOfThisTri):
        """Returns an iterator of the edges of the given triangle in two forms

        nodesOfThisTri must be a list of node ids

        The return value is a tuple of two elements. Each of those states the
        edge by its node numbers. The first element of the return tuple is a
        frozenset, that means the nodes are not ordered. The second element is
        a tuple of the node ids in the correct order.

        similar to TriangleSurface.triEdgesIter()
        """
        for i1, i2 in ((0,1), (1,2), (2,0)):
            orderedEdge = (nodesOfThisTri[i1], nodesOfThisTri[i2])
            yield (frozenset(orderedEdge), orderedEdge)

    @staticmethod
    def triEdgesIter(nodesOfThisTri):
        """Returns an iterator of the edges of the given triangle

        nodesOfThisTri will be converted to a tuple of node ids

        >>> myTri = (3, 43, 12)
        >>> for edges in triEdgesIter(myTri):
        >>> # edges will be:
        >>> # frozenset((3, 43)), frozenset((43, 12)), frozenset((12, 3))

        @note: See triEdgesOrderedNodesIter if you also need the nodes in
        the correct order.

        """
        for i1, i2 in ((0,1), (1,2), (2,0)):
            yield frozenset((nodesOfThisTri[i1], nodesOfThisTri[i2]))

    def getTriNormalUsc(self, nodes):
        """Not meant to be called from "outside"

        Computes the unscaled ("Usc") normal of a triangle. I.e. the length of
        the normal vector is twice its area.

        Nodes come in counter-clockwise succession if you look from "above"
        (positive side), against the direction of the normal.
        """
        vectorAB=vector(self.nodeCoords[nodes[0]],self.nodeCoords[nodes[1]])
        vectorAC=vector(self.nodeCoords[nodes[0]],self.nodeCoords[nodes[2]])
        return cross(vectorAB,vectorAC)

    def getTriNormalSc(self, nodes):
        """DEPRECATED! Use
        L{bae.surface_03.getTriNormalSc}(self.nodeCoords, nodes) instead.

        Not meant to be called from "outside"

        Computes the scaled ("Sc") normal of a triangle. I.e. this normal
        vector is of length 1.

        Nodes come in counter-clockwise succession if you look from "above"
        (positive side), against the direction of the normal.
        """
        vectorAB=vector(self.nodeCoords[nodes[0]],self.nodeCoords[nodes[1]])
        vectorAC=vector(self.nodeCoords[nodes[0]],self.nodeCoords[nodes[2]])
        return norm(cross(vectorAB,vectorAC))

    def calcAngleToNeighbours(self, elNodes, default=None):
        """calculates the cosines of the angles between the normals of the tri
        given by elNodes and of all of its neighbours.

        There does not have to be a tri on the surface connecting elNodes,
        this is meant to check the angle before actually inserting a new tri.

        The sign of the angle is arbitrary, since node ordering is not taken
        into account.

        The result is always a list of three members. Where there is no
        neighbour, the value of the argument "default" is returned in that
        list.
        """
        edgeToTri = self.getEdgeToTri()
        elNodes = tuple(elNodes)
        n1 = getTriNormalSc(self.nodeCoords, elNodes)

        angles = [default, default, default]
        for edgeIdx, thisEdge in enumerate(self.triEdgesIter(elNodes)):
            nodes = set()
            for thisTri in edgeToTri[thisEdge]:
                nodes.update(self.elNodes[thisTri])
            nodes.difference_update(elNodes)
            if len(nodes)!=1:
                continue
            nodes.update(thisEdge)
            n2 = getTriNormalSc(self.nodeCoords, tuple(nodes))
            angles[edgeIdx] = dot(n1,n2)

        return angles

    #} end of some service functions

######################################################
#{ Classes for making TriangleSurface a subclass of Mesh


class FakeElShape(object):
    """Objects of this class look like a dict but they always
    return the same value.

    It is used in TriangleSurface class.
    """

    def __init__(self, elNodes, shape="TRI_L"):
        self.elNodes = elNodes
        self.shape = shape

    def __getitem__(self, el):
        if el not in self.elNodes:
            raise KeyError("Element %d not there." % el)
        return self.shape

    def __iter__(self):
        return self.elNodes.__iter__()

    def iterkeys(self):
        return self.__iter__()

    def iteritems(self):
        return ((el, self.shape) for el in self.elNodes)


class FakeShapeEl(object):
    """Objects of this class look like a dict but they always
    return the same value.

    It is used in TriangleSurface class.
    """

    def __init__(self, elNodes, shape="TRI_L"):
        self.elNodes = elNodes
        self.shape = shape
        simpleShape = shape.rsplit("_",1)[0]
        self.shapes = set((shape, simpleShape))

    def __getitem__(self, shape):
        """get set of element numbers of the specified shape.
        Both specifically linear or quad versions (like "TRI_L") as well as
        the general version (like "TRI") might be queried.

        @raises KeyError: if this shape is not present in the model.
        """
        if shape not in self.shapes:
            raise KeyError("Shape %s does not exist." % shape)
        return set(self.elNodes)

    def get(self, shape, default=None):
        """get set of element numbers of the specified shape.
        Both specifically linear or quad versions (like "TRI_L") as well as
        the general version (like "TRI") might be queried.

        @returns: a set of elements or the value of the default argument if this
        shape is not present in the model.
        """
        if shape in self.shapes:
            return set(self.elNodes)
        else:
            return default

    def __iter__(self):
        yield self.shape

    def iteritems(self):
        yield (self.shape, set(self.elNodes))
#} end of classes for making TriangleSurface a subclass of Mesh


## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------
class ElemFaceSurface(Surface):
    """
    A surface made of faces of (typically volume) elements of a mesh. This is
    very much related to Abaqus meshes. It uses Abaqus' element face numbering.

    In other words: This is an Abaqus element based surface, a surface on
    -typically continuum- elements in an Abaqus mesh.


    How to create an instance
    =========================
     - exterior surface of some elements:
       surf = ElemFaceSurface().L{addFreeSurfFromElset}()
     - exterior surface of the remainder of an elset: L{addSurfFromElsetHole}
     - top surface of a mesh: L{bae.mesh_01.Mesh.getTopSurface}()
     - load it from an Abaqus input file:
       surf = ElemFaceSurface().L{updateFromModel}(model, "mySurf")


    @ivar faceEl: {face key: set of element numbers} face key like "S1", "S2"
    @ivar mesh: a L{mesh_01.Mesh} or L{abq_model_02.Model} instance defining
        the elements referenced in faceEl
    """
    # Todo:
    #   - union, difference, intersection functions returning a new surface
    #     object

    @property
    def faceIdShapeCornerNodeIds(cls):
        """A dictionary storing node ids for all faces of all element shapes.
        Stores only the corner nodes for each face.

        In particular it's {faceId: elShapeFaceNodes-dict} with:
         - faceId a string like "S1", "S2", ...
         - elShapeFaceNodes a dict {elShape: node index tuple}
         - elShape a string like "TET_Q", "TRI_L", ....
         - node index tuple: e.g. for the S1 face of a tet it's (0,1,2)

        This is a property (function). It creates its value only on demand
        (when it's accessed for the first time) and then overwrites the
        property function with the actual value thus preventing the function
        from being called a second time.
        """
        faceIdShapeCornerNodeIds = defaultdict(dict)
        for shape in Mesh.elShapeFaceNodes:
            # get the faces of the corner nodes only
            faces = Mesh.elShapeFaceNodes[shape.split("_")[0]]
            for i, face in enumerate(faces):
                faceIdShapeCornerNodeIds["S%d" % (i+1)][shape] = face
        # turn into ordinary dict
        faceIdShapeCornerNodeIds.default_factory = None

        # store and return result
        # (Overwrite the property-function with the value just determined, so
        # this expensive calculation is only done once and only if needed at
        # all.)
        ElemFaceSurface.faceIdShapeCornerNodeIds = faceIdShapeCornerNodeIds
        return faceIdShapeCornerNodeIds

    def __init__(self, faceEl=None, mesh=None):
        """Can also serve as a copy constructor if only the first argument is
        supplied and of type ElemFaceSurface. In this case a deep copy of the
        faceEl item is performed, for self.mesh we simply store the same
        reference.

        Whereas if you supply faceEl as a dictionary then it is just stored as
        is!

        @param faceEl: Optional, initialize self.faceEl with this.
           If this argument is of type ElemFaceSurface and if it's the only
           argument then this function serves as a copy constructor and
           performes a deep copy of the faceEl attribute.
        @param mesh: optional, initialize self.mesh with this
        """
        if isinstance(faceEl, ElemFaceSurface) and mesh is None:
            # copy constructor
            # if the first argument is another ElemFaceSurface-instance
            self.mesh = faceEl.mesh
            self.faceEl = dict(
                (typ, set(elset))
                for typ, elset in faceEl.faceEl.iteritems())
        elif faceEl is None:
            self.faceEl = dict()
            self.mesh = mesh
        else:
            self.faceEl = faceEl
            self.mesh = mesh

    def __str__(self):
        """Give a short description of the surface.

        >>> surf = ElemFaceSurface().updateFromModel(model, "mySurf")
        >>> print surf
        ElemFaceSurface with 214 nodes, 201 triangles, ...
        ... min: [0.0,0.4,1.2], max: [5.3,8.4,3.2]

        """
        return ("ElemFaceSurface with %d faces (triangles)"
                % self.countFaces())

    def countFaces(self):
        "Returns the total number of faces (triangles) of this surface."
        return sum(map(len, self.faceEl.itervalues()))

    def _checkMeshArg(self, mesh, funcName):
        """Service function: Check if mesh-argument == self.mesh if both are
        given. Either could be None initially. If self.mesh is None then update
        with mesh-argument.

        Usage:
         >>> model = self._checkMeshArg(model, "addFreeSurfFromElset")
        """
        if self.mesh is None:
            self.mesh = mesh
        elif mesh is not None and self.mesh is not mesh:
            raise ValueError(
                "ElemFaceSurface.%s(): Given mesh or model argument"
                " must be identical to self.mesh and it's not." % funcName)
        if not isinstance(self.mesh, Mesh):
            raise ValueError(
                "ElemFaceSurface.%s(): Given mesh or model argument or"
                " self.mesh must be a bae.mesh_01.Mesh instance. It's not."
                % funcName)
        return self.mesh

    def _checkElsetArg(self, elset, funcName):
        """Service function: Check and convert elset-argument for
        addFreeSurfFromElset, addSurfFromElsetHole

        Usage:
         >>> elset = self._checkElsetArg(elset, "addFreeSurfFromElset")
        """
        checkElems = False
        if elset is None:
            elset = set(self.mesh.elNodes)
        elif isinstance(self.mesh, AbqModel):
            if not isinstance(elset, set):
                # if we have an abq_model_02.Model as self.mesh
                # then we can use Model.getUnionSet()
                elset = self.mesh.getUnionSet("elset", elset)
                checkElems = True
        else:
            # if self.mesh is a mesh_01.Mesh object then the elset
            # argument must be a set of element numbers
            assert(isinstance(self.mesh, Mesh))
            if not isinstance(elset, set):
                # convert lists or so to a set
                elset = set(elset)
                checkElems = True

        if checkElems:
            if not elset:
                raise ValueError(
                    "Elset argument to ElemFaceSurface.%s() is empty. It must"
                    " be an iterable of elset numbers defined in the mesh."
                    % funcName)
            foundElems = elset.intersection(self.mesh.elNodes)
            if len(foundElems) < len(elset):
                # at least some elements in elset not defined in self.mesh!
                if not foundElems:
                    raise ValueError(
                        "Elset argument to ElemFaceSurface.%s() must be an"
                        " iterable of elset numbers defined in the mesh. Did"
                        " not find any valid element numbers. Instead got"
                        " e.g.:\n%s"
                        % (funcName, sorted(elset.difference(foundElems))[:20]))
                else:
                    msg("WARNING: ElemFaceSurface.%s() is ignoring %d elements"
                        " from the given elset of %d elements because they are"
                        " not defined in the mesh."
                        % (funcName, len(elset)-len(foundElems), len(elset)))
                    elset = foundElems
        return elset

    def addFreeSurfFromElset(self, model, elset=None, ndConstraints=dict()):
        """Adds to self the perimeter surface / exterior surface
        of the given elset.

        See also: L{addFreeSurfFromSplitMesh} if you want cohesive elements
        considered as internal part of the mesh and not identifying (part of)
        the exterior surface.

        @param model: A L{bae.abq_model_02.Model} or L{bae.mesh_01.Mesh}
           instance.
        @param elset: may be a set of element numbers or an element set name
           or a list of element set names, anything that Model.getUnionSet
           accepts as input.
           If not given or None then take the free surface of the whole model.
        @param ndConstraints: may be a dict of {slaveNdLabel:masterNdLabel, ...}
           if for example two nodes are connected via a MPC.
        @returns: self, so that surf=ElemFaceSurface().addFreeSurfFromElset(...)
           is possible.
        """

        model = self._checkMeshArg(model, "addFreeSurfFromElset")
        elset = self._checkElsetArg(elset, "addFreeSurfFromElset")

        faceToElem = defaultdict(list)
        # for performance reasons take the following if-then-else construct
        # outside and do the loops inside
        if not ndConstraints:
            for thisElem in elset:
                for faceNum, face in enumerate(model.elemFacesIter(thisElem)):
                    faceToElem[face].append([thisElem,faceNum])
        else:
            for thisElem in elset:
                for faceNum, face in enumerate(model.elemFacesIter(thisElem)):
                    faceNew = frozenset(ndConstraints.get(nd,nd) for nd in face)
                    faceToElem[faceNew].append([thisElem,faceNum])

        freeSurfElems = defaultdict(set)
        for face, elems in faceToElem.iteritems():
            if len(elems)==1:
                freeSurfElems['S%d' %(elems[0][1]+1)].add(elems[0][0])
            else:
                pass

        self.faceEl.update(freeSurfElems)
        return self

    def addFreeSurfFromSplitMesh(
            self, model, elset=None,
            volElemShape="TET", cohsvElShape="WEDGE"):
        """Adds to self the perimeter surface / exterior surface
        of the given elset in the given model. Cohesive elements are considered
        to connect adjacent tets in a way as if the mesh had not been split:

        Contrary to this method (L{addFreeSurfFromSplitMesh}) the method
        L{addFreeSurfFromElset} finds quadratic tet element faces that are
        connected to cohesive elements as part of the free surface.
        This method (L{addFreeSurfFromSplitMesh}) avoids this by sort of
        temporarily reconnecting the tet faces on either side of each cohesive
        element.

        @param model: A L{bae.abq_model_02.Model} or L{bae.mesh_01.Mesh}
           instance.
        @param elset: may be a set of element numbers or an element set name
           or a list of element set names, anything that Model.getUnionSet
           accepts as input.
           If not given or None then take the free surface of the whole model.

        @param volElemShape: Filter: only faces on elements of this shape are
           considered to form the resulting free surface.
        @param cohsvElShape: Element shape of elements to be considered as
           cohesive elements.

        @returns: self, so that surf=ElemFaceSurface().addFreeSurfFromElset(...)
           is possible.
        """

        model = self._checkMeshArg(model, "addFreeSurfFromSplitMesh")

        msg("ElemFaceSurface.addFreeSurfFromSplitMesh: Determining cohsNodes.",
            debugLevel=1)
        elNodes = model.elNodes  # abbreviation
        # cohsNodes = (elNodes[cohs] for cohs in model.shapeEl[cohsvElShape])
        try:
            cohsNodes = (elNodes[cohs] for cohs in model.shapeEl[cohsvElShape])
            msg("ElemFaceSurface.addFreeSurfFromSplitMesh: Determined cohsNodes"
            " for %d cohesive elements."
            % len(model.shapeEl[cohsvElShape]), debugLevel=1)
        except:
            msg("The mesh contain no cohesive elements!")
            cohsNodes = []
        
        msg("ElemFaceSurface.addFreeSurfFromSplitMesh: Preparing oldToNew.")
        # create node pairs, sorted ascending
        # Some notes on the algorithm:
        # - cohsNodePairs will be iterated from small node numbers to large
        # - assignments oldToNew will assign from n2/large/old to n1/small/new
        #   node numbers
        # - Therefore in most cases the assignment destination (n1/new) in
        #   oldToNew will not have to be changed later in the process of
        #   assembling it because we won't find n2->n1 with smaller n2 anymore.
        #   Because of the ascending order of cohsNodePairs and the node tuples
        #   in it.
        # - There is one exception: If a new n2->n1 appears with n2/large/old
        #   already reassigned earlier, i.e. we have an earlier assignment
        #   n2->n1' already stored in oldToNew. In this case we store n1->n1'
        #   (or the inverse n1'->n1 such that we maintain the rule to assign
        #   from large to small). Now we also have to check all earlier
        #   oldToNew contents for assignments to a destination that is just
        #   about to be reassigned. I.e. find assignments n2"->n1 and replace
        #   their destinations n1 by n1'.
        cohsNodePairs = sorted(set(
            tuple(sorted(twoNodes))
            for twoNodes in (
                (nodes[i+len(nodes)/2], nodes[i])
                for nodes in cohsNodes
                for i in range(len(nodes)/2))
            if twoNodes[0]!=twoNodes[1]))
        msg("ElemFaceSurface.addFreeSurfFromSplitMesh: ... found %d node pairs"
            "..." % len(cohsNodePairs), debugLevel=1)

        # higher node number shall be replaced by lower
        oldToNew = dict()
        newToOld = defaultdict(list)
        for n1, n2 in cohsNodePairs:
            n1 = oldToNew.get(n1, n1)
            try:
                # check if the assignment destination n1/new has already been
                # reassigned earlier.
                n2 = oldToNew[n2]
            except KeyError:
                pass
            else:
                # higher node already assigned earlier to oldToNew
                n1, n2 = sorted((n1, n2))
                # n2 might be already assigned as lower number (destination in
                # oldToNew). Check all earlier assignment destinations:
                try:
                    nn2List = newToOld[n2]
                except KeyError:
                    pass
                else:
                    for nn2 in nn2List:
                        oldToNew[nn2] = n1
                    del newToOld[n2]
                    newToOld[n1].extend(nn2List)
            # in any case, store result: assign higher node to lower node
            oldToNew[n2] = n1
            newToOld[n1].append(n2)
        msg("ElemFaceSurface.addFreeSurfFromSplitMesh: ... found %d pairs of"
            " duplicate nodes." % len(oldToNew), debugLevel=1)

        # determine volume elements
        if elset is None:
            elset = model.shapeEl[volElemShape]
        else:
            if isinstance(self.mesh, AbqModel):
                if not isinstance(elset, set):
                    # if we have an abq_model_02.Model as self.mesh
                    # then we can use Model.getUnionSet()
                    elset = self.mesh.getUnionSet("elset", elset)
            # intersect with volElemShape
            elset = model.shapeEl[volElemShape].intersection(elset)
        msg("Considering %d volume elements." % len(elset), debugLevel=10)

        # prepare node ids for faces
        # tuple of face-number and such that normal points out
        faceNodeIds = model.elShapeFaceNodes[
            volElemShape.replace("_Q", "_L")]
        faceNodeIds = tuple(
            (cnt, ids[::-1]) for cnt, ids in enumerate(faceNodeIds))

        # collect elements on faces
        # faces contain node numbers, modified by oldToNew where applicable
        msg("ElemFaceSurface.addFreeSurfFromSplitMesh: Examining face"
            " connectivity.", debugLevel=1)
        faceToElem = defaultdict(list)
        for thisElem in elset:
            nodes = [oldToNew.get(node, node) for node in elNodes[thisElem]]
            for faceNum, nodeIds in faceNodeIds:
                faceToElem[frozenset(nodes[i] for i in nodeIds)].append(
                    [thisElem,faceNum])

        # create exteriour surface
        msg("ElemFaceSurface.addFreeSurfFromSplitMesh: Identifying exterior"
            " faces.", debugLevel=1)
        freeSurfElems = defaultdict(set)
        for face, elems in faceToElem.iteritems():
            if len(elems)==1:
                freeSurfElems['S%d' %(elems[0][1]+1)].add(elems[0][0])

        # store result in ElemFaceSurface object
        msg("ElemFaceSurface.addFreeSurfFromSplitMesh: Storing results.",
            debugLevel=1)
        self.faceEl.update(freeSurfElems)
        return self

    def addSurfFromElsetHole(self, model, elset, ndConstraints=dict(),
                             smallElsetRatio=0.5):
        """Adds to self the surface of the hole in the model when the given
        elset is extracted from the model. This is the opposite surface of the
        exterior / free surface of the elset.

        @param model: A L{bae.abq_model_02.Model} or L{bae.mesh_01.Mesh}
           instance.
        @param elset: The hole-elset. May be a set of element numbers or an
           element set name or a list of element set names, anything that
           Model.getUnionSet accepts as input.

        @param ndConstraints: (optional) May be a dict of
           {slaveNdLabel: masterNdLabel, ...} if for example two nodes are
           connected via a MPC. Each of the masterNdLabels must not appear
           as any of the slaveNdLabels.

        @param smallElsetRatio: Don't specify this value unless you have a
           special reason to do so. The default value should be a reasonable
           choice and may change when we learn more or when the implemented
           algorithm changes.

           Explanation: There are two alternative algorithms implemented: The
           first determines all faces on the hole-elset and the second
           determines all faces on the remaining elements. As this is an
           expensive part of the procedure we select the algorithm based on the
           ratio of numbers of elements in the given (hole-) elset vs the
           number of remaining elements.

           If set to 0 then the first algorithm for smaller elsets is always
           chosen. If set to 1 (or larger) then the second algorithm is always
           chosen. This can be used for testing purposes.

        @returns: self, so that surf=ElemFaceSurface().addSurfFromElsetHole(...)
           is possible.
        """

        model = self._checkMeshArg(model, "addSurfFromElsetHole")
        elset = self._checkElsetArg(elset, "addSurfFromElsetHole")

        # node ids for specific faces
        faceNodes = dict(
            ("S%d" % (i+1), nodeIds)
            for i, nodeIds in enumerate(self.mesh.elShapeFaceNodes["TET_L"]))
        nbFaceNodes = 3
        nbCornerNodes = self.mesh.elShapeToNbCornerNodes["TET"]

        if len(elset) < smallElsetRatio*len(model.elNodes):
            # If elset is considerably smaller than the rest of the model then
            # it's probably (not tested yet) more efficient to determine the
            # exterior surface of elset, then find all nodes on this surface
            # and finally determine elements connected to these nodes and faces.

            # faces (keys in faceToElem) are frozensets of node numbers
            faceToElem = defaultdict(list)
            # for performance reasons take the following if-then-else construct
            # outside and do the loops inside
            if not ndConstraints:
                for thisElem in elset:
                    for face in model.elemFacesIter(thisElem):
                        faceToElem[face].append(thisElem)
            else:
                for thisElem in elset:
                    for face in model.elemFacesIter(thisElem):
                        faceNew = frozenset(ndConstraints.get(nd,nd)
                                            for nd in face)
                        faceToElem[faceNew].append(thisElem)

            # collect all faces with only one element attached
            faces = set(face for face, elems in faceToElem.iteritems()
                        if len(elems)==1)

            # find search set: all nodes on self
            allNodes = set().union(*faces)

            # find search set: all elements attached to nodes on self
            # but not in elset
            connectedElems = [
                (el, nodes)
                for el, nodes in model.elNodes.iteritems()
                if el not in elset
                and len(allNodes.intersection(nodes[:nbCornerNodes]))
                >= nbFaceNodes]

            # test each element in connectedElems if and how it's connected
            faceEl = [
                (faceId2, elset2)
                for faceId2, elset2 in (
                    (faceId,
                     set(elem for elem, nodes in connectedElems
                         if frozenset(nodes[i] for i in nodeIds) in faces))
                    for faceId, nodeIds in faceNodes.iteritems())
                if elset2]

        else:
            # If elset is about half of the model or even bigger then it's
            # probably (not tested yet) more efficient to determine all nodes
            # on elset, then all faces of the remainder of model and then
            # identify all faces that consist entirely of nodes attached to
            # elset as the faces we are looking for --the faces of the
            # elset-hole in model.
            elsetNodes = model.getConnectedNodes(elset)
            elemsRemaining = set(model.elNodes).difference(elset)

            # faces (keys in faceToElem) are frozensets of node numbers
            faceToElem = defaultdict(list)
            # for performance reasons take the following if-then-else construct
            # outside and do the loops inside
            if not ndConstraints:
                for el in elemsRemaining:
                    for faceNum, face in enumerate(model.elemFacesIter(el)):
                        if len(face.intersection(elsetNodes))==nbFaceNodes:
                            faceToElem[face].append([el,faceNum])
            else:
                for el in elemsRemaining:
                    for faceNum, face in enumerate(model.elemFacesIter(el)):
                        faceNew = frozenset(ndConstraints.get(nd,nd)
                                            for nd in face)
                        if len(faceNew.intersection(elsetNodes))==nbFaceNodes:
                            faceToElem[faceNew].append([el,faceNum])

            faceEl = defaultdict(set)
            for face, elems in faceToElem.iteritems():
                if len(elems)==1:
                    faceEl['S%d' %(elems[0][1]+1)].add(elems[0][0])

        # store result
        self.faceEl.update(faceEl)
        return self

    def updateFromModel(self, model, surfaceName):
        """Update self from the definition in an L{bae.abq_model_02.Model}
        instance, i.e. from an Abaqus input file.

        Example:
        ========
         >>> model = Model().read("MyModelWithSurfaces.inp")
         >>> surf = ElemFaceSurface().updateFromModel(model, "mySurf")

        @param model: an L{bae.abq_model_02.Model} instance
        @param surfaceName: Surface name of type string

        @returns: self, so that surf=ElemFaceSurface().updateFromModel(...)
           is possible.

        @note: This stores references of the element sets. If you modify this
          data you may also modify the elset the surface referes to.
        """

        model = self._checkMeshArg(model, "updateFromModel")

        try:
            surfElemSets = model.surface[surfaceName]
        except KeyError:
            raise KeyError("Surface %s is not defined in model." % surfaceName)

        for faceId, elsetName in surfElemSets.iteritems():
            self.faceEl[faceId] = model.elset[elsetName]
        return self

    def updateFromSequence(
            self, model=None, elsetsExcav=None, elsetsFill=None,
            ignoreElsets=set(), ignoreElems=set()):
        """Return an iterator that step-by-step updates self in place by adding
        sequence elsets according to elsetsFill and removing others according
        to elsetsExcav.

        B{PLEASE} use keyword arguments only because the order of arguments may
        change!


        Example:
        ========
         >>> from bae.abq_model_02 import Model
         >>> from bae.makeHistory_03 import sequenceElsetsIter
         >>> model = Model().read("mymodel.inp")
         >>> topSurf = model.getTopSurface()
         >>> topSurfIter = topSurf.updateFromSequence(
         >>>     model=model,
         >>>     elsetsExcav=sequenceElsetsIter(
         >>>         odb.postData, frameList=frameList,
         >>>         seqElsets=sequenceSetDictNameExcav,
         >>>         incremental=True),
         >>>     elsetsFill=sequenceElsetsIter(
         >>>         odb.postData, frameList=frameList,
         >>>         seqElsets=sequenceSetDictNameFill,
         >>>         incremental=True))
         >>> for (step, frame) in frameList:
         >>>     topSurfIter.next()
         >>>     surf = topSurf.getTriangleSurface()
         >>>     surf.exportAsSTL("export_%2d.stl" % frame)


        Algorithm:
        ==========

        We start with the current configuration of self. For each step the
        parameter elsetsExcav provides a set of elements to be removed. Its
        exterior surface is being divided into one part intersecting the
        current configuration of self and the rest being disjoint from self.
        The first part is then being subtracted from self and the second part
        is being added to self.

        Then elsetsFill is being processed as the set of elements to be added:
        Firstly the part of the current surface (self) that will be covered by
        the new fill gets subtracted from self. Secondly the remaining part of
        the new fill is being added to self.


        @param model: L{abq_model_02.Model} object supplying the elset
           definitions. Note that the mesh itself is taken from self.
           If the items in elsetsExcav and elsetsFill are set of elements
           then this is not required and defaults to None.

        @param elsetsExcav: An iterable supplying a list for each time step
           that defines what's going to be extracted for this step.

           If the model argument has been supplied (i.e. if model is not None)
           then items can be anything that
           L{model.getUnionSet<bae.abq_model_02.Model.getUnionSet>} accepts,
           i.e. lists of elset names or element numbers.
           If the model argument has not been supplied or if model is None
           then those items must be sets of element numbers.
           The elsetsExcav argument can also be None if nothing is to be
           extracted.

        @param elsetsFill: Like the elsetsExcav-argument but for elements to be
           added: An iterable supplying a list of elset names or element
           numbers to be added to the model configuration for each time step.
           Or None if nothing is to be added.

        @param ignoreElsets: set of elset names to be ignored, i.e. disregarded
           where they appear (by their names) in elsetsExcav and elsetsFill.

        @param ignoreElems: set of element numbers to be ignored, i.e. removed
           from any set of elements before it's being added or subtracted from
           the current configuration.

        @returns: self (the updated ElemFaceSurface) at each step.

        @Note: self.mesh must contain the complete mesh, at least all elements
           that ever become part of the surface during the course of this
           procedure.

        @Note: Provided sequence elsets must connect seamlessly and be disjoint
           for the algorithm to work. See above for an explanation of the
           algorithm.
        """
        # convert sequence provider to iterators
        if elsetsExcav:
            elsetsExcav = iter(elsetsExcav)
        if elsetsFill:
            elsetsFill = iter(elsetsFill)

        while 1:
            if elsetsExcav:
                elsetNamesExcav = elsetsExcav.next()
                msg("New excavation: got %d items to extract."
                    % len(elsetNamesExcav), debugLevel=5)
                if model and isinstance(elsetNamesExcav, (tuple, list)):
                    elsetNamesExcav = set(
                        elsetNamesExcav).difference(ignoreElsets)
                    newExcavElems = model.getUnionSet("elset", elsetNamesExcav)
                elif not isinstance(elsetNamesExcav, set):
                    raise ValueError(
                        "bae.surface_03.ElemFaceSurface.updateFromSequence: If"
                        " the model argument has not been supplied then the"
                        " items in elsetsExcav must be sets of element numbers."
                        " Instead we got: %s." % type(elsetNamesExcav))
                else:
                    newExcavElems = elsetNamesExcav
                msg("New excavation: after filtering and resolving items for"
                    " actual elements got %d elements to extract."
                    % len(newExcavElems), debugLevel=5)
                if len(newExcavElems)==0:
                    msg("No elements to excavate.\nWe have %d elsets: %s"
                        "\nWe have a getUnionSet-attribute: %s."
                        "\nOf the elsets to excavate %d are actually *not* in"
                        " the model."
                        % (len(elsetNamesExcav), sorted(elsetNamesExcav)[:10],
                           hasattr(model, "getUnionSet"),
                           (not model and -99) or len(
                               elsetNamesExcav.difference(model.elset))),
                        debugLevel=5)
                # filter only elements in self.mesh
                newExcavElems.intersection_update(self.mesh.elNodes)
                newExcavElems.difference_update(ignoreElems)
                msg("New excavation for this time step of %d elements..."
                    % len(newExcavElems), debugLevel=1)

                # modify self
                if newExcavElems:
                    # ... take away free surface of new excavation
                    newExcavFreeSurf = ElemFaceSurface().addFreeSurfFromElset(
                        self.mesh, newExcavElems)
                    toRemoveFromTopsurf = self.intersection(newExcavFreeSurf)
                    self.difference_update(toRemoveFromTopsurf)

                    # ... add hole (opposite side) of new excavation
                    newExcavFreeSurf.difference_update(toRemoveFromTopsurf)
                    newExcavHoleSurf = newExcavFreeSurf.invert()
                    self.update(newExcavHoleSurf)
                    msg("Updated top surface. Now has %d faces."
                        % self.countFaces(), debugLevel=5)

            if elsetsFill:
                elsetNamesFill = elsetsFill.next()
                msg("New refill: got %d items to add."
                    % len(elsetNamesFill), debugLevel=5)
                if model and isinstance(elsetNamesFill, (tuple, list)):
                    elsetNamesFill = set(
                        elsetNamesFill).difference(ignoreElsets)
                    newFillElems = model.getUnionSet("elset", elsetNamesFill)
                elif not isinstance(elsetNamesFill, set):
                    raise ValueError(
                        "bae.surface_03.ElemFaceSurface.updateFromSequence: If"
                        " the model argument has not been supplied then the"
                        " items in elsetsFill must be sets of element numbers."
                        " Instead we got: %s." % type(elsetNamesFill))
                else:
                    newFillElems = elsetNamesFill
                msg("New refill: after filtering and resolving items for"
                    " actual elements got %d elements to add."
                    % len(newFillElems), debugLevel=5)
                # filter only elements in self.mesh
                newFillElems.intersection_update(self.mesh.elNodes)
                newFillElems.difference_update(ignoreElems)
                msg("New refill for this time step of %d elements..."
                    % len(newFillElems), debugLevel=1)

                # modify self
                if newFillElems:
                    # ... create exterior surface and hole (opposite side) of
                    # fill
                    newFillFreeSurf = ElemFaceSurface().addFreeSurfFromElset(
                        self.mesh, newFillElems)
                    newFillHoleSurf = newFillFreeSurf.invert()

                    # ... take away hole (opposite side) of new fill
                    toRemoveFromTopsurf = self.intersection(newFillHoleSurf)
                    self.difference_update(toRemoveFromTopsurf)

                    # ... add free surface of new fill
                    newFillFreeSurf.difference_update(
                        toRemoveFromTopsurf.invert())
                    self.update(newFillFreeSurf)
                    msg("Updated top surface. Now has %d faces."
                        % self.countFaces(), debugLevel=5)

            if not elsetsExcav and not elsetsFill:
                msg("WARNING: There is nothing to update for"
                    " ElemFaceSurface.updateFromSequence().")

            # ready for this time step: updated self
            yield self

    def storeInAbqModel(self, model, name):
        """Add this surface to abq_model with a specified name, i.e. add
        corresponding elsets and model.surface - dictionary

        @param model: abq_model_02.Model instance
        @param name: name of the new surface
        """
        if name in model.surface:
            raise ValueError("Surface with the specified name %s already exists"
                             " in the model." % name)
        surfdict = dict()
        for typ in sorted(self.faceEl):
            elsetName = "%s_%s" % (name, typ)
            if elsetName in model.elset:
                cnt = 0
                template = "%s_%%d" % elsetName
                while elsetName in model.elset:
                    cnt += 1
                    elsetName = template % cnt
            if len(self.faceEl[typ])>0:
                model.elset[elsetName] = self.faceEl[typ]
                surfdict[typ] = elsetName
        model.surface[name] = surfdict


    def getTriangleSurface(self, quadTetTo4Tris=False, returnVolRef=False):
        """Return a L{TriangleSurface} created from self.

        The surface must consist of faces of tet elements (C3D4 or C3D10M).
        Only the corner nodes of the tets are considered if the parameter
        quadTetTo4Tris is not set to True.

        The node coords of this surface are a reference (no copy) to those of
        the abaqus model self.mesh. If self.mesh moves, so does this surface.

        Normals of the resulting triangle surface point outward.

        @param quadTetTo4Tris: If True then create four tris for each quadratic
        tet face. If False then only the corner nodes of the tets are
        considered and one (linear) tri is created for each tet face.

        @param returnVolRef: If True then additionally return a reference to
        the original volume element and face. The result will be a tuple
        (tri-surf, vol-ref) with tri-surf being the TriangleSurface object and
        vol-ref being a dict {tri elem label: (vol elem label, faceId,
        Gauss-pt-interpol-data)}.
        The Gauss-pt-interpol-data is a list of lists of weights for the
        interpolation of Gauss-pt values from the original volume elements to
        the triangle elements representing the surface.

        Example for a surface on linear tet elements with only one resulting
        triangle per tet face:
        vol-ref = {1: (3245931, "S2", [[1.0,],]), 2: (3558741, "S4", [[1.0,],]),
        ...}.
        The Gauss-pt-interpol data is a list of values for all Gauss points on
        the target element (triangle). Each item is a list of weights to
        multiply with the respective Gauss point of the source element.

        Example for a surface on quadratic tet elements with four resulting
        triangles per tet face:
        vol-ref = {1: (3245931, "S2", [[0, 1, 0, 0],]), 2: (3245931, "S2",
        [[0, 0, 0, 1],]), 3: (3245931, "S2", [[1, 0, 0, 0],]), 4: (3245931,
        "S2", [[1/3, 1/3, 0, 1/3]],]), 5: (3558741, "S4", [[1, 0, 0, 0],]),
        6: (3558741, "S4", [[0, 0, 0, 1],]), ...}.

        @Note: This implementation has to be cleaned up. Currently you should
        only use the quadTetTo4Tris parameter if the mesh consists of quadratic
        tets only.

        @Note: Only faces on tets and brick elements considered/implemented so
        far.

        @Note: Elements (triangles) in the resulting surface are always ordered
        in the same way, two calls to this method are guaranteed to produce the
        same output. (Resulting tri order results from alphabetic order of face
        identifiers followed by volume element numbers.)
        """

        volMesh = self.mesh  # abbreviation

        # connectivity data
        elShapeFaceType2TriNodeIds = {
            "TET_L": dict(
                # yields: {'S1':[[0,1,2],], 'S2':[[1,0,3],], 'S3':[[2,1,3],],
                #          'S4':[[0,2,3],]}
                # Note: nodeIds-list inverted in order to have outward pointing
                # normals.
                ("S%d" % (i+1), [nodeIds[::-1],])
                for i, nodeIds in enumerate(
                        volMesh.elShapeFaceNodes["TET_L"])),
            "TET_Q": dict(
                # tetNodeIdsOnFace are the six indexes in volMesh.elNode[el] of
                # the six nodes on one tet-face.
                # faceNodeIds4Tri is a list of the three indexes in
                # tetNodeIdsOnFace for each of the four tris on one tet face.
                ("S%d" % (i+1), [
                    [tetNodeIdsOnFace[j] for j in faceNodeIds4Tri]
                    for faceNodeIds4Tri in [
                            # index-triples for each triangle
                            # in list of tet's face nodes
                            [0,5,3], [5,2,4], [4,1,3], [5,4,3]]])
                for i, tetNodeIdsOnFace in enumerate(
                        volMesh.elShapeFaceNodes["TET_Q"])),
            "HEX_L": {
                "S1": [[1,0,2], [0,3,2]],  # two tri elements for face S1
                "S2": [[4,5,7], [5,6,7]],
                "S3": [[0,1,4], [1,5,4]],
                "S4": [[1,2,5], [2,6,5]],
                "S5": [[2,3,6], [3,7,6]],
                "S6": [[3,0,7], [0,4,7]],
                },
            }
        elShapeFaceType2gpInterpol = {
            "TET_L": dict(
                ("S%d" % (i+1), [[[1.0,],],]) for i in range(4)),
            "TET_Q": {
                "S1": [
                    [[1.0, 0.0, 0.0, 0.0],],
                    [[0.0, 0.0, 1.0, 0.0],],
                    [[0.0, 1.0, 0.0, 0.0],],
                    [[1.0/3, 1.0/3, 1.0/3, 0.0]],],
                "S2": [
                    [[0.0, 1.0, 0.0, 0.0],],
                    [[0.0, 0.0, 0.0, 1.0],],
                    [[1.0, 0.0, 0.0, 0.0],],
                    [[1.0/3, 1.0/3, 0.0, 1.0/3]],],
                "S3": [
                    [[0.0, 1.0, 0.0, 0.0],],
                    [[0.0, 0.0, 1.0, 0.0],],
                    [[0.0, 0.0, 0.0, 1.0],],
                    [[0.0, 1.0/3, 1.0/3, 1.0/3]],],
                "S4": [
                    [[1.0, 0.0, 0.0, 0.0],],
                    [[0.0, 0.0, 0.0, 1.0],],
                    [[0.0, 0.0, 1.0, 0.0],],
                    [[1.0/3, 0.0, 1.0/3, 1.0/3]],],},
            "HEX_L": dict(
                ("S%d" % (i+1), [[[1.0,],[1.0,]],]) for i in range(6)),
            }

        if not quadTetTo4Tris:
            elShapeFaceType2TriNodeIds["TET_Q"] = (
                elShapeFaceType2TriNodeIds["TET_L"])
            elShapeFaceType2gpInterpol["TET_Q"] = dict(
                (faceId, [weights[-1],])
                for faceId, weights in
                elShapeFaceType2gpInterpol["TET_Q"].iteritems())

        # collect tri nodes in this list
        elNodesList = list()
        if returnVolRef:
            volRefList = list()

        # main loop over faceTypes
        for faceType in sorted(self.faceEl):
            for el in sorted(self.faceEl[faceType]):

                # ... for each face ...
                elShape = volMesh.elShape[el]
                try:
                    # list of node ids for all tris for current face
                    # [ tri 1 - [node id 1, node id 2, node id 3],
                    #   tri 2 - [node id 1, node id 2, node id 3],
                    #   tri 3 - [...], tri 4 ...]
                    nodeIds4FaceTris = (
                        elShapeFaceType2TriNodeIds[elShape][faceType])
                    gpWeights = (
                        elShapeFaceType2gpInterpol[elShape][faceType])
                except KeyError:
                    raise NotImplementedError(
                        "%s.%s.getTriangleSurface() can not create a"
                        " surface from elements of shape %s."
                        % (__name__, self.__class__.__name__,
                           volMesh.elShape[el]))

                # update elNodesList
                volNodes = volMesh.elNodes[el]
                elNodesList.extend(
                    [volNodes[i] for i in triNodeIds]
                    for triNodeIds in nodeIds4FaceTris)

                # optionally store vol elem label, face id and Gauss-pt ids and
                # weights for (possibly multiple) triangles created from the
                # current vol elem face
                if returnVolRef:
                    volRefList.extend(
                        (el, faceType, i) for i in gpWeights)

        # store tri nodes in newsurf
        elNodes = dict(
            (i+1, nodes)
            for i, nodes in enumerate(elNodesList))

        # copy attached nodes
        nodestocopy = set(
            i for nodes in elNodesList for i in nodes)
        nodeCoords = dict(
            (node, volMesh.nodeCoords[node])
            for node in nodestocopy
            if node in volMesh.nodeCoords)

        newsurf = TriangleSurface(nodeCoords=nodeCoords, elNodes=elNodes)
        if len(newsurf.elNodes)==0:
            msg('WARNING: SurfaceFromMeshTris initialized with'
                ' empty element set.')

        if returnVolRef:
            volRef = dict((i+1,x) for i, x in enumerate(volRefList))
            return (newsurf, volRef)
        else:
            return newsurf


    def getTetFaceNormal(self):
        """Returns a dict of normal vectors of length one:
        {faceId: {element number: unit normal vector}}.

        Normals point outward.
        """
        # abbreviations
        nodeCoords = self.mesh.nodeCoords
        elNodes = self.mesh.elNodes
        allTets = self.mesh.shapeEl["TET"]

        triNormal = dict()
        for faceId, elems in self.faceEl.iteritems():
            elems = elems.intersection(allTets)
            # reverse nodeIds to get outward pointing normal
            nodeIds = self.faceIdShapeCornerNodeIds[faceId]["TET_L"][::-1]
            triNormal[faceId] = dict(
                (el,
                 getTriNormalSc(nodeCoords, [elNodes[el][i] for i in nodeIds]))
                for el in elems)

        return triNormal


    def getTetFaceCentroid(self):
        """Returns a dict of face centroids:
        {faceId: {element number: [x, y, z]}}.
        """
        # abbreviations
        nodeCoords = self.mesh.nodeCoords
        elNodes = self.mesh.elNodes
        allTets = self.mesh.shapeEl["TET"]

        centroids = dict()
        for faceId, elems in self.faceEl.iteritems():
            elems = elems.intersection(allTets)
            nodeIds = self.faceIdShapeCornerNodeIds[faceId]["TET_L"]
            weight = 1.0 / len(nodeIds)
            centroids[faceId] = dict(
                (el,
                 vector_scale(
                     vector_sum(*(nodeCoords[elNodes[el][i]] for i in nodeIds)),
                     weight))
                for el in elems)

        return centroids


    def getHexFaceNormal(self):
        """Returns a dict of normal vectors of length one:
        {faceId: {element number: unit normal vector}}.

        Normals point outward.

        Note: the face of a Hex is a Quad. The normal of the quad (with
        nodes 0,1,2,3) is determined as average of four tris (0,1,2),(0,2,3)
        and (0,1,3),(1,2,3).
        """
        # abbreviations
        nodeCoords = self.mesh.nodeCoords
        elNodes = self.mesh.elNodes
        allHexs = self.mesh.shapeEl["HEX"]

        triNormal = dict()
        for faceId, elems in self.faceEl.iteritems():
            elems = elems.intersection(allHexs)
            # reverse nodeIds to get outward pointing normal
            nodeIds = self.faceIdShapeCornerNodeIds[faceId]["HEX_L"][::-1]

            triNormal[faceId] = dict()
            for el in elems:
                n1 = getTriNormalSc(nodeCoords, [elNodes[el][nodeIds[0]],
                                                 elNodes[el][nodeIds[1]],
                                                 elNodes[el][nodeIds[2]]])
                n2 = getTriNormalSc(nodeCoords, [elNodes[el][nodeIds[0]],
                                                 elNodes[el][nodeIds[1]],
                                                 elNodes[el][nodeIds[2]]])
                n3 = getTriNormalSc(nodeCoords, [elNodes[el][nodeIds[0]],
                                                 elNodes[el][nodeIds[1]],
                                                 elNodes[el][nodeIds[2]]])
                n4 = getTriNormalSc(nodeCoords, [elNodes[el][nodeIds[0]],
                                                 elNodes[el][nodeIds[1]],
                                                 elNodes[el][nodeIds[2]]])
                triNormal[faceId][el] = [
                    .25*(nn1+nn2+nn3+nn4)
                    for nn1,nn2,nn3,nn4 in zip(n1,n2,n3,n4)]

        return triNormal

    def update(self, other):
        """
        Constructs a boolean add operator for Abaqus continuum element surfaces.
        @param other: ElemFaceSurface object to be added.
        """
        assert isinstance(other, ElemFaceSurface)

        self._checkMeshArg(other.mesh, "update")
        for face, elems in other.faceEl.iteritems():
            try:
                oldElems = self.faceEl[face]
            except KeyError:
                oldElems = set()
                self.faceEl[face] = oldElems
            oldElems.update(elems)

    def difference_update(self, other):
        """
        Constructs a boolean subtract operator for Abaqus continuum element
        surfaces.

        @param other: ElemFaceSurface object to be subtracted.

        @Note: The name difference_update based on the corresponding
        Python builtin set operation. differenceUpdate is a camelCase alias.
        """
        assert isinstance(other, ElemFaceSurface)

        self._checkMeshArg(other.mesh, "difference_update")
        for face, elems in other.faceEl.iteritems():
            try:
                oldElems = self.faceEl[face]
            except KeyError:
                continue
            else:
                oldElems.difference_update(elems)
                if not oldElems:
                    del self.faceEl[face]

    # camelCase alias for L{difference_update}
    differenceUpdate = difference_update

    def intersection_update(self, other):
        """Constructs a boolean intersect operator for Abaqus continuum element
        surfaces. I.e. remove all surface parts from self that are not found in
        other.

        @param other: ElemFaceSurface object to be intersected with self.

        @Note: The name intersection_update based on the corresponding
        Python builtin set operation. intersectionUpdate is a camelCase alias.
        """
        assert isinstance(other, ElemFaceSurface)

        self._checkMeshArg(other.mesh, "intersection_update")
        for face, elems in other.faceEl.iteritems():
            try:
                oldElems = self.faceEl[face]
            except KeyError:
                del self.faceEl[face]
            else:
                oldElems.intersection_update(elems)
                if not oldElems:
                    del self.faceEl[face]

    # camelCase alias for L{intersection_update}
    intersectionUpdate = intersection_update

    def intersection(self, other):
        """Constructs a boolean intersect operator for Abaqus continuum element
        surfaces.

        @param other: ElemFaceSurface object to be intersected with self.

        @returns: The intersection of self with other, a new
        L{ElemFaceSurface}-object.
        """
        assert isinstance(other, ElemFaceSurface)

        self._checkMeshArg(other.mesh, "intersection")

        newFaceEl = dict()
        for face, elems in other.faceEl.iteritems():
            try:
                intersectElems = self.faceEl[face].intersection(elems)
            except KeyError:
                continue
            else:
                if intersectElems:
                    newFaceEl[face] = intersectElems
        return ElemFaceSurface(faceEl=newFaceEl, mesh=self.mesh)

    def differenceUpdateNset(self, nodes):
        """Removes all faces that have B{all} corner nodes in the given node
        set.

        @param nodes: set of nodes
        """
        assert isinstance(nodes, set)

        # abbreviations
        elShape = self.mesh.elShape
        elNodes = self.mesh.elNodes

        for face, elems in self.faceEl.iteritems():
            shapeToNodeIds = self.faceIdShapeCornerNodeIds[face]

            # the following command is identical (but faster) than these three:
            # elNodesIdsIter = (
            #     (el, elNodes[el], shapeToNodeIds[elShape[el]])
            #     for el in elems)
            # elFaceIter = (
            #     (el, set(nodes[i] for i in nodeIds))
            #     for el, nodes, nodeIds in elNodesIdsIter)
            # foundElems = [
            #     el for el, faceNodes in elFaceIter
            #     if faceNodes.issubset(nodes)]
            foundElems = [
                el
                for el in elems
                if set(elNodes[el][i] for i in shapeToNodeIds[elShape[el]])
                .issubset(nodes)]

            elems.difference_update(foundElems)

    def difference(self, other):
        """Constructs a boolean difference operator for Abaqus continuum element
        surfaces.

        @param other: ElemFaceSurface object to be subtracted from self.

        @returns: The difference of self with other, a new
        L{ElemFaceSurface}-object.
        """
        assert isinstance(other, ElemFaceSurface)

        self._checkMeshArg(other.mesh, "difference")

        newFaceEl = dict()
        for face, elems in self.faceEl.iteritems():
            try:
                diffElems = elems.difference(other.faceEl[face])
            except KeyError:
                newFaceEl[face] = set(elems)
            else:
                if diffElems:
                    newFaceEl[face] = diffElems
        return ElemFaceSurface(faceEl=newFaceEl, mesh=self.mesh)

    def intersectionUpdateNset(self, nodes):
        """Removes all faces that don't have B{all} corner nodes in the given
        node set.

        @param nodes: set of nodes
        """
        assert isinstance(nodes, set)

        # abbreviations
        elShape = self.mesh.elShape
        elNodes = self.mesh.elNodes

        for face, elems in self.faceEl.iteritems():
            shapeToNodeIds = self.faceIdShapeCornerNodeIds[face]
            # see description in method differenceUpdateNset
            # for the following expression...
            foundElems = [
                el
                for el in elems
                if set(elNodes[el][i] for i in shapeToNodeIds[elShape[el]])
                .issubset(nodes)]
            elems.intersection_update(foundElems)

    def invert(self, wholeMesh=None):
        """Get the opposite element based surface: All opposing faces to those
        on self. Those faces belong to other (opposing) volume (i.e. tet-)
        elements. Obviously this will only return something for internal
        surfaces and not for the external surface of self.mesh.

        @param wholeMesh: (optional) a L{bae.abq_model_02.Model} or
            L{bae.mesh_01.Mesh} object that contains all elements of self and
            of the new inverted surface. If not given then self.mesh is used.

            This option might be useful if self.mesh does not contain all
            elements of the opposing side.

        @returns: a new L{ElemFaceSurface} object.

        @Note: "Degenerate" cases: If two opposing faces belong to self then
            the result will not contain either of those faces.
            If an element belongs to self with one face it can belong to the
            result with another face. This can happen if the face normals are
            not unified, i.e. if the "outside" direction swaps across an edge
            (for two adjacent faces).
        """

        # use self.mesh if wholeMesh not given
        if wholeMesh:
            mesh = wholeMesh
        else:
            mesh = self.mesh

        # node ids for specific faces
        faceNodes = dict(
            # { (elShape, faceId) : [nodeIds] }
            ((elShape, "S%d" % (i+1)), nodeIds)
            for elShape in mesh.shapeEl
            for i, nodeIds in enumerate(mesh.elShapeFaceNodes[elShape]))

        # collect all faces as frozensets of node numbers
        faces = set(
            frozenset(nodes[i] for i in faceNodes[(elShape, faceType)])
            for faceType, elset in self.faceEl.iteritems()
            for nodes, elShape in (
                    (mesh.elNodes[elem], mesh.elShape[elem])
                    for elem in elset))

        # find search set: all nodes on self
        allNodes = set().union(*faces)

        # minimum number of nodes per face for a particular element shape
        # to later compare to the number of nodes actually connected to the
        # interface --i.e. self
        shapeNbFaceNodes = dict(
            (elShape, min(len(nodeIds)
                          for nodeIds in mesh.elShapeFaceNodes[elShape]))
            for elShape in mesh.shapeEl)

        # Find search set: all elements attached to nodes on self.
        # Serves as a pre-select filter, not all those elements must really
        # share a face with the interface / with self...
        connectedElems = [
            (el, nodes)
            for el, nodes in mesh.elNodes.iteritems()
            if (len(allNodes.intersection(nodes))
                >= shapeNbFaceNodes[mesh.elShape[el]])]
        # alternative, may be faster or slower... (can possibly be simplified?) 
        # connectedElems = [
        #     (el, nodes)
        #     for (el, nodes, nbFaceNodes) in (
        #             (el2, mesh.elNodes[el2], nbFaceNodes1)
        #             for elset, nbFaceNodes1 in (
        #                     (elset2, shapeNbFaceNodes[elShape])
        #                     for elShape, elset2 in mesh.shapeEl.iteritems())
        #             for el2 in elset)
        #     if (len(allNodes.intersection(nodes)) >= nbFaceNodes)]

        # Test each element in connectedElems if and how it's connected.
        # After this step faceEl contains all faces connected to the surface
        # self on either sides. I.e. faceEl contains both: self and the
        # opposing side that we're looking for.
        faceEl = [(faceId,
                   set(elem
                       for elem, nodes in connectedElems
                       if frozenset(nodes[i] for i in nodeIds) in faces))
                  for (elShape, faceId), nodeIds in faceNodes.iteritems()]

        # remove faces of self
        faceEl = dict(
            (faceId2, elset2)
            for faceId2, elset2 in (
                    (faceId, elset.difference(self.faceEl.get(faceId, set())))
                    for faceId, elset in faceEl)
            if elset2)

        # prepare result
        result = ElemFaceSurface(faceEl=faceEl, mesh=mesh)
        return result


## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------
