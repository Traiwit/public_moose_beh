"""volume_02 - module

Some functionality concerning a distinct volume in space

Used primarily to find elements by their location in space.

Usage:
 - create a volume (instantiate any of the Volume... classes)
 - found_elements = volume.intersectionElset(model, searchset)
"""

__version__ = "2.23"
_version_history_ = """
Versions:

1.01 new   - pt_is_inside() and initialize_point_search() derived
             from findAtMesh 1.17
             scale, copy, getboundingbox methods
1.02 added - translate method, more doc, added smooth mode for point finding
1.03 added - intersection_elset,
             don't need to specify elset for create_from_tet_elem
     fixed: all from future
     added: VolumeBase class, VolumeExtrudeOutline,
              Volume.get_exterior_surface() method
     changed: introduced tolerance in point_is_inside methods

2.01 changed: use abq_model_02, base class is Volume, old Volume class ->
         VolumeFromTetMesh, change most names to CamelCase
     added: VolumeSphere
2.02 changed: rename VolumeFromTetMesh to VolumeFromMesh and add element
         conversion from arbitrary volume elements to tets
2.03 added: VolumeBoxRotZ, getCentroid() for the box volumes
2.04 GP added: VolumeBetweenSurfaceAndZPlane, added intersectionIdPointTuples,
         intersectionIdSomePointsTuples(), gauss point method for
         intersectionElset()
2.05 GP added: means for element interpolation from nodal values:
         VolumeFromMesh.pointIsInside() returnLinShapeFunc argument added
2.06 GP changed: INCOMPATIBLE INTERFACE CHANGE! renamed first argument of
         VolumeFromMesh.initializePointSearch() and added second argument.
         Changed the order of arguments to initializePointSearchInBox()
2.07 GP added: VolumeInSurfZPrj
2.08 GP added: VolumeFromMesh.getCentroidVolume()
2.09 GP INCOMPATIBLE INTERFACE CHANGE! added: elementConversion option to
        VolumeFromMesh.__init__() which is now off by default. Before
        elements (other than tets) were always converted. To reduce memory
        consumption and procesing time.
2.10 GP changed to bae.log_01
2.11 GP changed pointIsInside uses 6D KDTree
2.12 GP added: accept Mesh instance as input to VolumeFromMesh
2.13 GP reinstall old pointIsInside
2.14 GP fixed: in VolumeFromMesh.__init__(): Model.shapeEl.__getitem__ recently
        changed to throws KeyError, therefore use ...get() instead
2.15 GP added surface based volumes should now also accept TriangleSurface
        objects from bae.surface_03 instead of solely surface_02.
2.16 GP added: VolumeInSurfCheckRays
2.17 GP added: VolumeVoxel
2.18 GP added: intersectionGrid as general method for all Volumes and
        streamlined the special version of VolumeInSurfCheckRays.
2.19 GP added: VolumeVoxel.getConnectedVolumePart(), fromFieldOnPtGrid()
        and getCellCentroidsGrid()
2.20 GP changed: getExteriorSurface() yields a surface_03-class by default
2.21 GP incompatible change to VolumeVoxel.fromVtkFieldOnCellCentr()
        insideFilterForFirstField argument: all non-zero values indicate
        "inside". Before all values>0.5 indicated "inside".
2.22 GP added projectionAxis to VolumeInSurfCheckRays,
        improved VolumeInSurfCheckRays.intersectionGrid()
2.23 GP added VolumeBox.intersectionGrid

todo:
* VolumeInSurfCheckRays: make pointIsInside aware of self.projectionAxis.
  handle projectionAxis-tuple: check first axis as usual, then check dubious
  points on the remaining axes.
* rotate method
* pointIsInside ... faster! (which one?)
"""

from itertools import izip
from math import pi, sin, cos, sqrt, ceil, floor
import bisect

from bae.vecmath_01 import \
    vector, vector_minus, vector_plus, vector_scale,\
    vector_modif_add,\
    dist
from bae.future_01 import defaultdict
from bae.misc_01 import BoundingBox, KDTree
from bae.abq_model_02 import Model, checkModelModuleVersion
from bae.mesh_01 import Mesh, MeshStructuredPoints
from bae.field_01 import Field
from bae.surface_03 import TriangleSurface
from bae.vtk_02 import FieldsOnStructPointGrid as Vtk

from bae.log_01 import msg, MsgTicker

class Volume(object):
    """Base class for all different types of volumes.

    all derived classes must provide the following member functions:
     - self.pointIsInside(point)
     - self.getBoundingBox()
     - self.scale(scale_factor, scale_origin)
     - self.translate(move_vector)
    and might provide
     - self.getExteriorSurface()
    """

    def __init__(self, logfile=None):
        "@param logfile: deprecated for compatibility, ignored"
        pass

    def intersectionIdPointTuples(self, idPointTuples):
        """
        @param idPointTuples: A sequence (list or iterator) of tuples
          (point id, coordinates tuple).
          Or a dictionary {point id: coordinates tuple}.

        @returns: A set of point ids for which the corresponding points have
          been found inside the volume.
        """
        # check if all checkpoints are in the volume
        ids_inside = set()
        if isinstance(idPointTuples, dict):
            idPointTuples = idPointTuples.iteritems()
        ticker = MsgTicker("... intersecting point list: tested %d points so"
                           " far and already found %d in the volume.")
        for cnt, (id, point) in enumerate(idPointTuples):
            if self.pointIsInside(point):
                ids_inside.add(id)
            ticker.msg(cnt+1, len(ids_inside))
        del ticker
        return ids_inside


    def intersectionIdSomePointsTuples(self, idPointsTuples):
        """
        Find which groups of points lie completely (all individual points of
        the group) in the volume self.

        @param idPointsTuples: A sequence (list or iterator) of tuples
          (group id, point list) with point list being a list of points
          (coordinates tuples). Or it may be a dictionary
          {point id: coordinates tuple}.

        @returns: A set of group ids for which the corresponding points have
          all been found inside the volume.

        @note: The same group id may occur multiple times. If any of the
          corresponding point groups are completely within the volume, this
          group id is in the resulting set.

        @note: If all groups consist of only one item, use
          self.intersectionIdPointTuples().
        """
        # check if all checkpoints are in the volume
        ids_inside = set()
        if isinstance(idPointsTuples, dict):
            idPointsTuples = idPointsTuples.iteritems()
        ticker = MsgTicker("... intersecting point list: tested %d points so"
                           " far and already found %d in the volume.")
        for cnt, (id, points) in enumerate(idPointsTuples):
            if all([self.pointIsInside(pt) for pt in points]):
                ids_inside.add(id)
            ticker.msg(cnt+1, len(ids_inside))
        del ticker
        return ids_inside


    def intersectionElset(self, model, elset=None, method='centroid'):
        """find all elements from the elset in model that are within the
        volume

        @param model: an abq_model_02.Model object or a filename (a string) or
              an open file of an abaqus input file.

        @param elset: states the "searchset". Ignore the rest of model.
              might be anything that model.getUnionSet() accepts as input:
              an elset name, a set of element numbers, a list of elset names
              and element numbers
              elset might also be None, then all elements of model are tested.
              Element numbers in elset that are not defined in the model are
              silently ignored.

        @param method: optional method parameter
              May be 'centroid'/'c' or 'gauss point'/'gp'

        @returns: a set of element numbers.

        @note: bugs/limitations:
         - uses model.L{getElCentroids<bae.abq_model_02.model.Model.getElCentroids>}()
           to calculate the centroid of each element. See its limitations.
         - gauss point method not tested yet (?)

        @note: If the same model and elset is applied to different volumes or
        the volume changes between different invocations of
        self.intersectionElset(), consider the following aproach for the
        centroid method:
         >>> elemCentroid = model.getElCentroids(elset)
         >>> output = Model()
         >>> for cnt, vol in enumerate(volumesToCheck):
         >>>     output.elset["ELSET_%02d" % (cnt+1)] = (
         ...         vol.intersectionIdPointTuples(elemCentroid))
         >>> output.write("myresults.inp")

        ...and for the gauss point method:
         >>> elemGPts = model.getQuadGaussPointsDict(elset)
         >>> output = Model()
         >>> for cnt, vol in enumerate(volumesToCheck):
         >>>     output.elset["ELSET_%02d" % (cnt+1)] = (
         ...         vol.intersectionIdSomePointsTuples(elemGPts))
         >>> output.write("myresults.inp")

        @note: multivol...intersectionElset() might be faster than calling this
        function multiple times for different volumes with the same model and
        elset to intersect.
        """

        method = method.replace(' ', '').lower()

        if isinstance(model, (basestring, file)):
            m = Model()
            m.read(model)
            model = m

        if method in ('centroid', 'c'):
            el_pts = model.getElCentroids(elems=elset)
            els_inside = self.intersectionIdPointTuples(el_pts)
        elif method in ('gausspoint', 'gp'):
            el_pts = model.getQuadGaussPointsDict(elset=elset)
            els_inside = self.intersectionIdSomePointsTuples(el_pts)
        else:
            raise Exception('Volume.intersectionElset: method <%s> not'
                            ' implemented so far.' % method)

        return els_inside


    def intersectionElsetOld(self, model, elset=None, method='centroid'):
        """find all elements from the elset in model that are within the
        volume

        @param model: an abq_model_02.Model object or a filename (a string) or
              an open file of an abaqus input file.

        @param elset: states the "searchset". Ignore the rest of model.
              might be anything that model.getUnionSet() accepts as input:
              an elset name, a set of element numbers, a list of elset names
              and element numbers
              elset might also be None, then all elements of model are tested.
              Element numbers in elset that are not defined in the model are
              silently ignored.

        @param method: optional method parameter
              May be 'centroid'/'c' or 'gauss point'/'gp'

        @returns: a set of element numbers.

        @note: bugs/limitations:
         - uses model.getCentroid() to calculate the centroid of each element
           see its limitations
         - gauss point method not implemented yet.

        @note: If the same model and elset is applied to different volumes or
        the volume changes between different invocations of
        self.intersectionElset(), consider the following aproach:
         >>> elemCentroid = model.getElCentroidDict(elset)
         >>> output = Model()
         >>> for cnt, vol in enumerate(volumesToCheck):
         >>>     output.elset["ELSET_%02d" % (cnt+1)] = (
         ...         vol.intersectionIdPointTuples(elemCentroid))
         >>> output.write("myresults.inp")

        @note: multivol...intersectionElset() might be faster than calling this
        function multiple times for different volumes with the same model and
        elset to intersect.
        """

        method = method.replace(' ', '')

        if isinstance(model, (basestring, file)):
            m = Model()
            m.read(model)
            model = m
        else:
            checkModelModuleVersion(model, "%s.%s.intersectionElset"
                                    % (__name__, self.__class__.__name__))

        if elset is None:
            # caution, this is fast but dangerous if you modify the code
            # in this case here, elset is a dictionary, not a set
            elset = model.elNodes
        else:
            elset = model.getUnionSet("elset", elset, onlyValid=True)

        # create list with tuples (el numer, [points to check])
        el_pts = list()
        if method in ('centroid', 'c'):
            ticker = MsgTicker("... calculated %s/%d element"
                               " centroids so far." % ("%d", len(elset)))
            el_pts = [None]*len(elset)
            for cnt, el in enumerate(elset):
                el_pts[cnt] = (el, [model.getCentroid(el), ])
                ticker.msg(cnt+1)
            del ticker
        else:
            raise Exception('Volume.intersectionElset: method <%s> not'
                            ' implemented so far.' % method)

        # check if all checkpoints are in the volume
        els_inside = set()
        ticker = MsgTicker("... intersecting model: tested %s/%d"
                           " and found %s elements in the volume."
                           % ("%d", len(elset), "%d"))
        for cnt, (el, points) in enumerate(el_pts):
            if all([self.pointIsInside(pt) for pt in points]):
                els_inside.add(el)
            ticker.msg(cnt+1, len(els_inside))
        del ticker

        return els_inside

    def intersectionGrid(self, grid, fieldName="insideVolume"):
        """Find points of the regular grid inside self.

        Example: Get a field indicating points in or outside the volume:
         >>> from bae.mesh_01 import MeshStructuredPoints
         >>> from bae.volume_02 import Volume....
         >>> from bae.vtk_02 import FieldsOnStructPointGrid as Vtk
         >>>
         >>> vol = Volume....(...)
         >>> grid = MeshStructuredPoints(....)
         >>> insideField = vol.intersectionGrid(grid)
         >>>
         >>> vtk = Vtk(mesh=grid)
         >>> vtk.updateFields(insideField)
         >>> vtk.toVtk("inVolume.vtk")

        @param grid: L{bae.mesh_01.MeshStructuredPoints}-object or a
            L{MeshStructuredPointsRot} object.

        @param fieldName: fieldName attribute of the field to be returned

        @returns: L{bae.field_01.Field}-object (position="structPt",
            dataType="scalar") containing 1 for "in", 0 for "out"

        @Note: This is a fallback function. It might be very slow and there
        may be a simple way to implement a much faster algorithm for a
        particular Volume-subclass. (If it's not there already...)
        """
        result = Field.classFromPosType(fieldName, "structPt", "scalar")(
            int(self.pointIsInside(point)) for point in grid.getPointsIter())
        return result


###################################################################
###################################################################
###################################################################
###################################################################

class VolumeFromMesh(Volume):
    """Volume class represented by the volume elements of a given Abaqus model.

    Example:

    >>> # initialization
    >>> from bae.volume_02 import VolumeFromMesh
    >>> VolumeFromMesh.setPointSearchMethod(method="cellGrid")
    >>> v = VolumeFromMesh(model)
    >>>
    >>> # test if nodes are in that volume
    >>> v.initializePointSearch() # only needed if specifying smooth_findradius
    >>> for point in [[0,0,0], [1.0, 0.4, 1.6], [-0.1, 2.0, 0.3]]:
    >>>    print point, "is inside:", v.pointIsInside(point)

    VolumeFromTetMesh is deprecated. It's a synonym for VolumeFromMesh and
    provided for compatibility reasons.

    @Note: that VolumeFromMesh.pointIsInside may yield unexpected behaviour if
    there is an element number of zero in the model, so please avoid that.
    """

    cellSizeFactor=0.5      # cell Size for array relative to 'avgSize'

    def __init__(self, model, elset=None, logfile=None,
                 elementConversion=False):
        """a volume created from some elements of an abaqus model

        @param model: an L{abq_model_02.Model} or L{mesh_01.Mesh} object or a
              filename (a string) or an open file of an abaqus input file.

        @param elset: the volume shall consist of the elements in this elset
              might be anything that model.getUnionSet() accepts as input:
              an elset name, a set of element numbers, a list of elset names
              and element numbers
              elset might also be None, then all elements of model are used.

        @param logfile: deprecated argument, only for compatibility reasons.
              No effect.

        @param elementConversion: If model contains volume elements other than
              tets those elements can be converted to tets. This takes longer
              and needs more memory and is therefore disabled by default.

        @note: No deep copy of the nodeCoords member of model is performed. If
        model.nodeCoords changes, it will affect this object. On the other
        hand from model.elNodes, the elements in elset are being deep copied
        and possibly converted to tets.

        @note: Elements in elset that are not found in model.elNodes are
        silently ignored.
        """

        Volume.__init__(self)

        # internal variable, use getBoundigBox method! (initialize it)
        self.boundingBox = None

        # variables used for finding points in the volume already set?
        self.point_search_initialized_6DTree = False
        self.point_search_initialized_cellGrid = False

        if isinstance(model, (basestring, file)):
            # basestring = any kind of string, or file
            m = Model()
            m.read(model)
            model = m
        elif not isinstance(model, Mesh):
            raise ValueError(
                "VolumeFromMesh.__init__: The supplied model is an <%s> object"
                " from the module %s. It's got to be an instance of Mesh."
                % (model.__class__.__name__, model.__class__.__module__))

        self.nodeCoords = model.nodeCoords
        if not len(self.nodeCoords):
            raise RuntimeError(
                "VolumeFromMesh.__init__: The supplied model does not have"
                " any nodes, nodeCoords is empty.")

        if elset is None:
            elset = set(model.elNodes)
        elif hasattr(model, "getUnionSet"):
            elset = model.getUnionSet("elset", elset)

        if elementConversion:
            if not(hasattr(model, "convertToTets")):
                msg("WARNING: elementConversion argument will be ignored if"
                    " the model supplied is not an abq_model_02.Model"
                    " instance")
            else:
                submodel = model.copySubmodelFromElsets(elset, "ELEMENT")
                submodel.convertToTets()

                # remove anything that is not a tet
                # note: it only updates elNodes because that's all we need but
                # submodel is corrupt after that
                for typ, elems in submodel.typeEl.iteritems():
                    if submodel.elTypeToShape[typ] != "TET":
                        for el in elems:
                            del submodel.elNodes[el]

                self.elNodes = submodel.elNodes

        else:
            elemsToCopy = model.shapeEl.get("TET_L", set()).union(
                model.shapeEl.get("TET_Q",set()))
            elemsToCopy.intersection_update(elset)
            if len(elemsToCopy)>100000:
                msg("Copying elements to the new VolumeFromMesh object.")
            self.elNodes = dict(
                (elem, model.elNodes[elem][:4])
                for elem in elemsToCopy )
            if len(elemsToCopy)>100000:
                msg("...Copied %d elements." % len(self.elNodes))

        if len(self.elNodes)==0:
            msg('WARNING: VolumeFromMesh initialized with empty element set.')

        return

    @classmethod
    def setPointSearchMethod(cls, method):
        """Choose which of the point search methods should be used.

        @param method: "6DTree" or "cellGrid"
        """
        if method == "6DTree":
            cls.initializePointSearch = \
                cls.initializePointSearch_6DTree
            cls.initializePointSearchInBox = \
                cls.initializePointSearchInBox_6DTree
            cls.pointIsInside = \
                cls.pointIsInside_6DTree
        elif method == "cellGrid":
            cls.initializePointSearch = \
                cls.initializePointSearch_cellGrid
            cls.initializePointSearchInBox = \
                cls.initializePointSearchInBox_cellGrid
            cls.pointIsInside = \
                cls.pointIsInside_cellGrid
        else:
            raise ValueError(
                "Got %s as argument for setPointSearchMethod where only"
                " '6DTree' or 'cellGrid' would be applicable."
                % method)

    ## ------------------------------------------------------------------------
    ## some info

    def getBoundingBox(self):
        """Bounding Box of the whole volume.

        Returns a list of two coordinate tuples (rather lists)
        The first coord tuple states the min value, the last the max.
        I.e. self.getBoundingBox()[i][j]:
        i \in {0,1} for min/max value, j \in {0,1,2} coordinate index

        Example:

        >>> # bounding box of that volume
        >>> bb = v.getBoundingBox()
        >>> print "lower left front corner:", bb[0]
        >>> print "upper right back corner:", bb[1]
        """
        if not(self.boundingBox):
            self.boundingBox = BoundingBox()
            self.boundingBox.update(
                nodeCoords=self.nodeCoords, elNodes=self.elNodes)

        return self.boundingBox

    def getCentroidVolume(self):
        """Calculate centroid and volume.

        Usage:
         >>> vol = VolumeFromMesh(mymodel, "myelset")
         >>> centroid, size  = vol.getCentroidVolume()[0]
         >>> print "Elset myelset at %s, size %g" % (centroid, size)

        To get just the centroid:
         >>> centroid  = vol.getCentroidVolume()[0]

        @Returns: A tuple: the coordinate tuple (rather list) of the centre of
        this volume and the volume itself
        """

        if len(self.elNodes)==0:
            msg("WARNING: %s.%s.getCentroidVolume(): No tets in this"
                " volume!" % (__name__, self.__class__.__name__))
            return [0.0, 0.0, 0.0], 0.0

        allVol = 0.0
        allCentre = [0.0, 0.0, 0.0]
        avg_n=4
        negVolWarn = False
        for el, nodes in self.elNodes.iteritems():

            points = [self.nodeCoords[node] for node in nodes]

            # centroid
            x=y=z=0.0
            for coords in points:
                x=x+coords[0]/avg_n
                y=y+coords[1]/avg_n
                z=z+coords[2]/avg_n

            # volume
            vec = [vector_minus(pt, points[0])
                   for pt in points[1:]]
            vol = ((
                  vec[0][0]*vec[1][1]*vec[2][2]
                + vec[0][1]*vec[1][2]*vec[2][0]
                + vec[0][2]*vec[1][0]*vec[2][1]
                - vec[0][2]*vec[1][1]*vec[2][0]
                - vec[0][1]*vec[1][0]*vec[2][2]
                - vec[0][0]*vec[1][2]*vec[2][1])
                   / 6.0)

            if vol<0:
                negVolWarn = True
                vol = -vol
            allVol += vol
            vector_modif_add(allCentre,
                             vector_scale([x,y,z], vol))
            continue

        if negVolWarn:
            msg("WARNING: %s.%s.getCentroidVolume(): Tets with"
                " negative size in this volume!"
                % (__name__, self.__class__.__name__))

        if allVol==0.0:
            msg("WARNING: %s.%s.getCentroidVolume(): The given volume"
                " has zero size. Cannot determine centroid properly."
                % (__name__, self.__class__.__name__))
            return [0.0, 0.0, 0.0], 0.0

        allCentre = vector_scale(allCentre, 1.0/allVol)
        return allCentre, allVol


    ## ------------------------------------------------------------------------
    ## get the exterior surface of this volume as TriangleSurface object

    def getExteriorSurface(self, nodeCoordsByRef=False, surfaceClass=None):
        """Return a L{TriangleSurface<bae.surface_03.TriangleSurface>} object
        with all outer triangles of this volume.

        Outer triangles are those with exactly one tet connected to.

        IMPORTANT NOTE: Node coordinates are passed by reference, if the volume
        is moved or scaled, the surface changes accordingly. The nodeCoords
        dictionary however by default is created new containing only connected
        nodes. See the parameter nodeCoordsByRef.

        @param nodeCoordsByRef: if this is True then take a reference of
           self.nodeCoords as the nodeCoords-attribute of the resulting surface
           object. This should speed up the task.

           IMPORTANT NOTE: Only specify this as keyword argument to ensure
           compatibility with future versions.

        @param surfaceClass: Optional class of the result. If not specified
           the result will be of class L{bae.surface_03.TriangleSurface}.
           This is intended for old scripts that require a
           L{bae.surface_02.TriangleSurface}-object.

           IMPORTANT NOTE: Only specify this as keyword argument to ensure
           compatibility with future versions.
        """

        faceOrder = dict()
        faceToNbTets = defaultdict(int)
        ticker = MsgTicker("computing exterior surface: faces of"
                           " tet elements %s/%d." % ("%d", len(self.elNodes)))
        for cnt, nodesOfThisTet in enumerate(self.elNodes.itervalues()):
            ticker.msg(cnt+1)
            for face in self.orientedFacesIter(nodesOfThisTet):
                id_ = frozenset(face)
                faceOrder[id_] = face
                faceToNbTets[id_] += 1
        del ticker

        ticker = MsgTicker()
        extFaces = [id_ for id_, nb in faceToNbTets.iteritems() if nb==1]
        ticker.msg("computing exterior surface: found %d exterior faces."
                   % len(extFaces))
        del ticker, faceToNbTets

        # now add the actual faces to the surf object
        ticker = MsgTicker()
        elNodes = dict(
            (cnt+1, list(faceOrder[id_]))
            for cnt, id_ in enumerate(extFaces) )
        ticker.msg("computing exterior surface: generated %d triangles."
                   % len(extFaces))
        del ticker, extFaces

        if nodeCoordsByRef:
            nodeCoords = self.nodeCoords

        else:

            # find all nodes
            ticker = MsgTicker()
            allNodes = set(n for face in elNodes.itervalues() for n in face)
            ticker.msg("computing exterior surface: identified %d connected"
                       " nodes." % len(allNodes))
            del ticker

            # update nodeCoords
            ticker = MsgTicker()
            nodeCoords = dict(
                (node, self.nodeCoords[node]) for node in allNodes)
            ticker.msg("computing exterior surface: copied %d nodes."
                       % len(allNodes))
            del ticker, allNodes

        # create the actual TriangleSurface object
        if isinstance(surfaceClass, type):
            moduleTail = surfaceClass.__module__.rsplit(".")[-1]
            if moduleTail=="surface_02":
                surf = surfaceClass(nodeCoords=nodeCoords, triNodes=elNodes)
            elif moduleTail=="surface_03":
                surf = surfaceClass(nodeCoords=nodeCoords, elNodes=elNodes)

        # default result type
        else:
            surf = TriangleSurface(nodeCoords=nodeCoords, elNodes=elNodes)

        return surf


    ## ------------------------------------------------------------------------
    ## methods for testing if a point lies inside the volume or outside
    ## Those are genereric placeholders that do nothing. Look at
    ## setPointSearchMethod() and the specific variants for details

    def initializePointSearch_6DTree(self, smoothRadius=0.0, smoothFactor=0.0):
        """Call this function once before pointIsInside.
        (Otherwise it is done automatically py pointIsInside)
        Builds a kd-tree for O(n log n) nearest neighbour search.

        Use setPointSearchMethod() to choose the actual implementation

        @Note: The volume-object self must not be empty, must contain at least
        one tet. In the current implementation...
        """
        pass

    def initializePointSearchInBox_6DTree(
            self, box=None, smoothRadius=0.0, smoothFactor=0.0):
        """This is a variant of self.initializePointSearch()

        Use setPointSearchMethod() to choose the actual implementation
        """
        pass

    def pointIsInside(
            self, point, returnVolElem=False, returnLinShapeFunc=False):
        """test if point is inside the volume

        Use setPointSearchMethod() to choose the actual implementation

        point is tuple of coordinates

        If returnVolElem==True the function returns the element number of the
        tet element the point is in or False if the point is not inside the
        volume at all. Otherwise it just returns True or False.

        If returnLinShapeFunc==True the function returns a tuple
        (element number, [SF1, SF2, SF3, SF4]) SF1...SF4 being the values of
        the (linear tetrahedron) shape function associated with the four nodes
        of the tet at the given point. Each of the SF1...SF4 values lies
        between 0 and 1. They can be used to linearily interpolate values at
        the corner nodes. Smoothmode must be off (that's the case unless you
        specify a nonzero smooth_findradius argument to the function
        initializePointSearch().

        @Note: If the volume is created from elements other then tets which
        have been converted to tets, this method returns the id of the
        converted tet.

        @Note: If there is an element number of zero, the typical term
        "if vol.pointIsInside(pt, True)" may yield unexpected behaviour.
        Instead you could test like that:
           >>> if vol.pointIsInside(pt, True)!=False:

        Warning: Don't test the type of the return value with the built in
        isinstance() function. The type bool is a subclass of int!
        """
        pass


    ## ------------------------------------------------------------------------
    ## methods for testing if a point lies inside the volume or outside
    ## this is the classic method that stores a list of intersecting elements
    ## for each cell of a structured homogenous grid

    def initializePointSearch_cellGrid(
            self, smoothRadius=0.0, smoothFactor=0.0):
        """call this function once before pointIsInside
        (it is done automatically)

        If either smoothRadius or smoothFactor are given, run "smoothmode":
        Additionally consider points to be inside the volume that are in fact
        outside the volume but within the "smooth find radius" to the closest
        centre point of any of the elements of the volume.

        The "smooth find radius" is set to smoothRadius if given, this is the
        priority if smoothFactor is given as well. Otherwise it is set to
        (smoothFactor * avarage element size) of self.
        """

        # dict of centroid-coord-tuples, should be initialized in any case
        self.tetCentroids=dict()

        if len(self.elNodes)==0:
            return


        msg("Initializing point search, element boxes...")
        ticker = MsgTicker("... processed %s/%d tet elements"
                           %("%d", len(self.elNodes)))

        # list of [[xmin,ymin,zmin],[xmax,ymax,zmax]] for each element
        elementBox = dict()

        nodeCoords = self.nodeCoords    # local copy for efficiency
        volBoundingBox = BoundingBox()  # later copied to self.boundingBox
        avgSize=0
        for cnt, (element, nodes) in enumerate(self.elNodes.iteritems()):

            # calculate centroid
            x=y=z=0
            for node in nodes:
                coords=nodeCoords[node]
                x=x+coords[0]/4.0
                y=y+coords[1]/4.0
                z=z+coords[2]/4.0
            self.tetCentroids[element]=[x,y,z]

            # calc element bounding box
            xyz = zip(*(self.nodeCoords[node] for node in nodes))
            newElemBox = [[min(xyz[0]),min(xyz[1]),min(xyz[2])],
                          [max(xyz[0]),max(xyz[1]),max(xyz[2])]]

            # update bounding box of all elements
            try:
                for i in xrange(3):  # x,y,z
                    if newElemBox[0][i]<volBoundingBox[0][i]:  # min
                        volBoundingBox[0][i]=newElemBox[0][i]
                    if newElemBox[1][i]>volBoundingBox[1][i]:  # max
                        volBoundingBox[1][i]=newElemBox[1][i]
            except IndexError:
                volBoundingBox[0]=newElemBox[0][:3]  # min, make a copy
                volBoundingBox[1]=newElemBox[1][:3]  # max, make a copy

            avgSize+=dist(newElemBox[1],newElemBox[0])
            elementBox[element] = newElemBox
            ticker.msg(cnt+1)

        del ticker
        avgSize /= len(self.elNodes)

        # in method pointIsInside() tolerance is compared to
        # dot(cross(...)...) of vectors of avgSize magnitude
        self.tolerance = (1E-4 * avgSize)**3

        msg('bounding box of the volume: %s' % volBoundingBox)
        msg('avg. element size: %g' % avgSize)


        ## -----------------------------------------------------------------
        ## map elements to cells

        msg("Initializing point search, mapping elements to cells...")

        cellOrigin=volBoundingBox[0]
        cellSize=avgSize*self.cellSizeFactor

        # x,y,z number of cells
        cellDimensions = [volBoundingBox[1][i]-volBoundingBox[0][i]
                          for i in range(3)]
        msg('cell array dimensions: %g %g %g'
            % tuple(cellDimensions))

        cellNumber = [int(cdim/cellSize)+1 for cdim in cellDimensions]
        msg('cell numbers: %g %g %g' % tuple(cellNumber))

        # create empty array
        self.cells=[]
        ticker = MsgTicker(
            "... initialized %s of %d cells"
            % ("%3.1f%%", cellNumber[0]*cellNumber[1]*cellNumber[2]))
        tickerscale = 100.0/(cellNumber[0]*cellNumber[1])
        for i in xrange(cellNumber[0]):
            LJ = []
            self.cells.append(LJ)
            for j in xrange(cellNumber[1]):
                LK = [[] for k in xrange(cellNumber[2])]
                LJ.append(LK)
                ticker.msg((i*cellNumber[1]+j)*tickerscale)
        del ticker

        # map elements
        ticker = MsgTicker("... processing tet element %s of %d"
                           %("%d", len(self.elNodes)))
        cellIndex=[[0,0,0],[0,0,0]]  # low,high
        for cnt, element in enumerate(self.elNodes):
            ticker.msg(cnt+1)
            # get index ranges
            for l in xrange(2):
                for i in xrange(3):
                    cellIndex[l][i]=int(
                        (elementBox[element][l][i]-cellOrigin[i])/cellSize)
            # append element ID to cells that overlap the element's bounding box
            for i in xrange(cellIndex[0][0],cellIndex[1][0]+1):
                for j in xrange(cellIndex[0][1],cellIndex[1][1]+1):
                    for k in xrange(cellIndex[0][2],cellIndex[1][2]+1):
                        self.cells[i][j][k].append(element)

        ticker.msgTemplate = ("... finished processing all %d tet elements."
                              % len(self.elNodes))
        ticker.msg()
        del ticker

        # store local values as object attributes
        # they have been local in the first place for performance reasons
        self.cellSize = cellSize
        self.boundingBox = volBoundingBox
        self.cellNumber = cellNumber

        # for smooth mode store radius
        self.smoothmode = (smoothRadius>0.0 or smoothFactor!=0.0)
        self.smooth_findradius = float(smoothRadius)
        if smoothRadius==0.0 and smoothFactor!=0.0:
            self.smooth_findradius = float(smoothFactor)*avgSize

        # finished initializing
        self.point_search_initialized_cellGrid = True
        return

    def initializePointSearchInBox_cellGrid(
            self, box=None, smoothRadius=0.0, smoothFactor=0.0):
        """This is a variant of self.initializePointSearch()

        Call this function once before pointIsInside, if you only want to find
        points that are guaranteed to lie with in a certain boundingBox that
        does not include the whole volume (self).

        This variant is especially useful when you want to find the tet
        elements that contain some points within a certain region of the
        volume.

        If either smoothRadius or smoothFactor are given, run "smoothmode":
        Additionally consider points to be inside the volume that are in fact
        outside the volume but within the "smooth find radius" to the closest
        centre point of any of the elements of the volume.

        The "smooth find radius" is set to smoothRadius if given, this is the
        priority if smoothFactor is given as well. Otherwise it is set to
        (smoothFactor * avarage element size) of self.

        @param smoothRadius: Additionally consider points to be inside the
        volume that are in fact outside the volume but within the smoothRadius
        to the closest centre point of any of the elements of the volume.

        @param smoothFactor: Additionally consider points to be inside the
        volume that are in fact outside the volume but within (smoothFactor *
        average element size of the mesh self consists of) to the closest
        centre point of any of the elements of the volume.

        @param box: All points to be tested with self.pointIsInside() certainly
        lie within in this box. box = [[xmin, ymin, zmin], [xmax, ymax, zmax]]

        @note: If a point is tested with the self.pointIsInside() method that
        does not lie within the specified box but in the volume (self), the
        result is unreliable, most probably self.pointIsInside() will state
        that this point is not in the volume.
        """

        # dict of centroid-coord-tuples, should be initialized in any case
        self.tetCentroids=dict()

        if len(self.elNodes)==0:
            return

        msg("Initializing point search, element boxes...")
        ticker = MsgTicker("... processed %s/%d tet elements"
                           %("%d", len(self.elNodes)))

        # list of [[xmin,ymin,zmin],[xmax,ymax,zmax]] for each element
        elementBox=dict()

        nodeCoords = self.nodeCoords    # local copy for efficiency
        volBoundingBox = BoundingBox()  # later copied to self.boundingBox
        relevantElements = list()
        avgSize=0
        for cnt, (element, nodes) in enumerate(self.elNodes.iteritems()):

            # calculate centroid
            x=y=z=0
            for node in nodes:
                coords=nodeCoords[node]
                x=x+coords[0]/4.0
                y=y+coords[1]/4.0
                z=z+coords[2]/4.0
            self.tetCentroids[element]=[x,y,z]

            # calc element bounding box
            xyz = zip(*node_pts)
            newElemBox = [[min(xyz[0]),min(xyz[1]),min(xyz[2])],
                          [max(xyz[0]),max(xyz[1]),max(xyz[2])]]

            # check if this elements bounding box intersects the box argument
            disjoint = False
            for e1,e2,b1,b2 in izip(newElemBox[0],newElemBox[1],box[0],box[1]):
                if min(e2,b2)<=max(e1,b1):
                    disjoint = True
                    break
            # ignore this element if it does not intersect box
            if disjoint:
                ticker.msg(cnt+1)
                continue

            # update bounding box
            try:
                for i in xrange(3):  # x,y,z
                    if newElemBox[0][i]<volBoundingBox[0][i]:  # min
                        volBoundingBox[0][i]=newElemBox[0][i]
                    if newElemBox[1][i]>volBoundingBox[1][i]:  # max
                        volBoundingBox[1][i]=newElemBox[1][i]
            except IndexError:
                volBoundingBox[0]=newElemBox[0][:3]
                volBoundingBox[1]=newElemBox[1][:3]

            relevantElements.append(element)
            avgSize += dist(newElemBox[1],newElemBox[0])
            elementBox[element] = newElemBox
            ticker.msg(cnt+1)

        del ticker
        if not len(relevantElements):
            avgSize = dist(*box)
            volBoundingBox.update(box)
        else:
            avgSize /= len(relevantElements)

        # in method pointIsInside() tolerance is compared to
        # dot(cross(...)...) of vectors of avgSize magnitude
        self.tolerance = (1E-4 * avgSize)**3

        # boundingBox = intersection(boundingBox, box)
        volBoundingBox[0] = map(max, volBoundingBox[0], box[0])  # min
        volBoundingBox[1] = map(min, volBoundingBox[1], box[1])  # max

        msg('bounding box of the relevant volume: %s' % volBoundingBox)
        msg('%d of the %d elements in the volume are within or'
            ' close to the limiting box.'
            % (len(relevantElements), len(self.elNodes)))
        msg('avg. element size: %g' % avgSize)


        ## -----------------------------------------------------------------
        ## map elements to cells

        msg("Initializing point search, mapping elements to cells...")

        cellOrigin=volBoundingBox[0]
        cellSize = avgSize*self.cellSizeFactor

        # x,y,z number of cells
        cellDimensions = [volBoundingBox[1][i]-volBoundingBox[0][i]
                          for i in range(3)]
        msg('cell array dimensions: %g %g %g' % tuple(cellDimensions))

        cellNumber = [int(cdim/cellSize)+1 for cdim in cellDimensions]
        msg('cell numbers: %g %g %g' % tuple(cellNumber))

        # create empty array
        self.cells=[]
        ticker = MsgTicker("... initialized %s of %d cells"
            % ("%3.1f%%", cellNumber[0]*cellNumber[1]*cellNumber[2]))
        tickerscale = 100.0/(cellNumber[0]*cellNumber[1])
        for i in xrange(cellNumber[0]):
            LJ = []
            self.cells.append(LJ)
            for j in xrange(cellNumber[1]):
                LK = [[] for k in xrange(cellNumber[2])]
                LJ.append(LK)
                ticker.msg((i*cellNumber[1]+j)*tickerscale)
        del ticker

        # map elements
        ticker = MsgTicker("... processing tet element %s of %d"
                           %("%d", len(relevantElements)))
        cellIndex=[[0,0,0],[0,0,0]]  # low,high
        for cnt, element in enumerate(relevantElements):
            ticker.msg(cnt+1)
            # get index ranges
            for i in xrange(3):
                cellIndex[0][i]=max(0, int(
                    (elementBox[element][0][i]-cellOrigin[i])/cellSize))
                cellIndex[1][i]=min(cellNumber[i]-1, int(
                    (elementBox[element][1][i]-cellOrigin[i])/cellSize))
            # append element ID to cells that overlap the element's bounding box
            for i in xrange(cellIndex[0][0],cellIndex[1][0]+1):
                for j in xrange(cellIndex[0][1],cellIndex[1][1]+1):
                    for k in xrange(cellIndex[0][2],cellIndex[1][2]+1):
                        self.cells[i][j][k].append(element)

        ticker.msgTemplate = ("... finished processing all %d tet elements."
                              % len(relevantElements))
        ticker.msg()
        del ticker

        # store local values as object attributes
        # they have been local in the first place for performance reasons
        self.cellSize = cellSize
        self.boundingBox = volBoundingBox
        self.cellNumber = cellNumber


        # for smooth mode store radius
        self.smoothmode = (smoothRadius>0.0 or smoothFactor!=0.0)
        self.smooth_findradius = float(smoothRadius)
        if smoothRadius==0.0 and smoothFactor!=0.0:
            self.smooth_findradius = float(smoothFactor)*avgSize

        # finished initializing
        self.point_search_initialized_cellGrid = True
        return

    def pointIsInside_cellGrid(
            self, point, returnVolElem=False, returnLinShapeFunc=False):
        """test if point is inside the volume
        point is tuple of coordinates

        If returnVolElem==True the function returns the element number of the
        tet element the point is in or False if the point is not inside the
        volume at all. Otherwise it just returns True or False.

        If returnLinShapeFunc==True the function returns a tuple
        (element number, [SF1, SF2, SF3, SF4]) SF1...SF4 being the values of
        the (linear tetrahedron) shape function associated with the four nodes
        of the tet at the given point. Each of the SF1...SF4 values lies
        between 0 and 1. They can be used to linearily interpolate values at
        the corner nodes. Smoothmode must be off (that's the case unless you
        specify a nonzero smooth_findradius argument to the function
        initializePointSearch().

        @Note: If the volume is created from elements other then tets which
        have been converted to tets, this method returns the id of the
        converted tet.

        @Note: If there is an element number of zero, the typical term
        "if vol.pointIsInside(pt, True)" may yield unexpected behaviour.
        Instead you could test like that:
           >>> if vol.pointIsInside(pt, True)!=False:

        Warning: Don't test the type of the return value with the built in
        isinstance() function. The type bool is a subclass of int!
        """

        if len(self.elNodes)==0:
            return False

        if not self.point_search_initialized_cellGrid:
            self.initializePointSearch_cellGrid()

        # cell check (includes bounding box test from cell index)
        index=[0,0,0]
        for i in xrange(3):
            index[i]=int((point[i]-self.boundingBox[0][i])/self.cellSize)
            if index[i]<0 or index[i]>=self.cellNumber[i]:
                return False

        # if cell check passed ... go on
        elems_in_cell=self.cells[index[0]][index[1]][index[2]]
        centroid_to_point_distance=[]
        # sort by distance
        for element in elems_in_cell:
            l=dist(point,self.tetCentroids[element])
            centroid_to_point_distance.append([l,element])
        centroid_to_point_distance.sort()

        # full check if point inside any tet of elems_in_cell
        inside=False

        if self.smoothmode:
            if (centroid_to_point_distance and
                (centroid_to_point_distance[0][0]<self.smooth_findradius)):
                inside = [centroid_to_point_distance[0], None]

        for [l,element] in centroid_to_point_distance:

            # coordinates of the four corner nodes
            x = [self.nodeCoords[node] for node in self.elNodes[element]]
            # coordinates relative to the first node
            # dp = [pi - x0i for pi, x0i in zip(point, x[0])]
            # dx = [[xxi-x0i for xxi, x0i in zip(xx, x[0])]
            #       for xx in x[1:]]
            dp = [point[0]-x[0][0], point[1]-x[0][1], point[2]-x[0][2]]
            dx = [[x[1][0]-x[0][0], x[1][1]-x[0][1], x[1][2]-x[0][2]],
                  [x[2][0]-x[0][0], x[2][1]-x[0][1], x[2][2]-x[0][2]],
                  [x[3][0]-x[0][0], x[3][1]-x[0][1], x[3][2]-x[0][2]]]
            # determinant of the "Jakobian"
            # cylcePos = [[0,1,2], [1,2,0], [2,0,1]]
            # detJ = sum(
            #     dx[0][ip[0]]*dx[1][ip[1]]*dx[2][ip[2]]
            #     - dx[2][ip[0]]*dx[1][ip[1]]*dx[0][ip[2]]
            #     for ip in cylcePos)
            detJ = (dx[0][0]*dx[1][1]*dx[2][2]
                    + dx[0][1]*dx[1][2]*dx[2][0]
                    + dx[0][2]*dx[1][0]*dx[2][1]
                    - dx[2][0]*dx[1][1]*dx[0][2]
                    - dx[2][1]*dx[1][2]*dx[0][0]
                    - dx[2][2]*dx[1][0]*dx[0][1])
            if detJ <= 0.0:
                # negative or zero volume, ignore this element
                continue

            # inverse * determinant
            # invJ = [ [a[ip2[1]]*b[ip2[2]] - a[ip2[2]]*b[ip2[1]] for ip2 in cylcePos]
            #          for a, b in ((dx[ip1[1]], dx[ip1[2]]) for ip1 in cylcePos)]
            invJ = [[dx[1][1]*dx[2][2] - dx[1][2]*dx[2][1],
                     dx[1][2]*dx[2][0] - dx[1][0]*dx[2][2],
                     dx[1][0]*dx[2][1] - dx[1][1]*dx[2][0]],
                    [dx[2][1]*dx[0][2] - dx[2][2]*dx[0][1],
                     dx[2][2]*dx[0][0] - dx[2][0]*dx[0][2],
                     dx[2][0]*dx[0][1] - dx[2][1]*dx[0][0]],
                    [dx[0][1]*dx[1][2] - dx[0][2]*dx[1][1],
                     dx[0][2]*dx[1][0] - dx[0][0]*dx[1][2],
                     dx[0][0]*dx[1][1] - dx[0][1]*dx[1][0]]]
            # parametric coordinates: xi_p := invJ * dp
            # xip = [float(sum(invJ[i][j]*dp[j] for j in range(3)))/detJ
            #        for i in range(3)]
            xip = [float(invJ[0][0]*dp[0]+invJ[0][1]*dp[1]+invJ[0][2]*dp[2])/detJ,
                   float(invJ[1][0]*dp[0]+invJ[1][1]*dp[1]+invJ[1][2]*dp[2])/detJ,
                   float(invJ[2][0]*dp[0]+invJ[2][1]*dp[1]+invJ[2][2]*dp[2])/detJ]

            shapeFuncs = [1-sum(xip), xip[0], xip[1], xip[2]]
            if all(xi>-self.tolerance for xi in shapeFuncs):
                inside = (element, shapeFuncs)
                break

        if not(returnLinShapeFunc) and inside is not False:
            if returnVolElem:
                inside = inside[0]
            else:
                inside = True

        return inside


    ## ------------------------------------------------------------------------
    ## methods for testing if a point lies inside the volume or outside
    ## using 6D KDTree

    def initializePointSearch_6DTree(self, smoothRadius=0.0, smoothFactor=0.0):
        """Call this function once before pointIsInside.
        (Otherwise it is done automatically py pointIsInside)
        Builds a kd-tree for O(n log n) nearest neighbour search.

        @param smoothRadius: no function, for compatibility only
        @param smoothFactor: no function, for compatibility only

        @Note: The volume-object self must not be empty, must contain at least
        one tet. In the current implementation...
        """

        assert len(self.elNodes)>0, "VolumeFromMesh.initializePointSearch" \
            "_KDTree needs at least one element!"

        labels = list()             # element labels
        points = list()             # 6D points == bounding boxes of tets
        centroids = list()          # normalized centroid of the tet in its box
        avgSize=0
        bbox = BoundingBox()             # just for diagnostic output
        nodeCoords = self.nodeCoords     # abbreviation for speed
        msg("Calculating bounding boxes for all %d tets." % len(self.elNodes),
            debugLevel=10)
        for elem, nodes in self.elNodes.iteritems():
            labels.append(elem)
            coordsT = zip(*(nodeCoords[n] for n in nodes))

            pt = [min(coordsT[0]), min(coordsT[1]), min(coordsT[2]),
                  -max(coordsT[0]), -max(coordsT[1]), -max(coordsT[2])]
            points.append( pt )
            # note: pt[3:6] := (-1) * tet bounding box upper limit
            avgSize += sqrt(
                ((pt[0]+pt[3])**2)+((pt[1]+pt[4])**2)+((pt[2]+pt[5])**2))

            bbox.update([pt[:3], [-x for x in pt[3:]]])
            centroids.append([
                    0.25*(coordsT[0][0]+coordsT[0][1]+coordsT[0][2]+coordsT[0][3]),
                    0.25*(coordsT[1][0]+coordsT[1][1]+coordsT[1][2]+coordsT[1][3]),
                    0.25*(coordsT[2][0]+coordsT[2][1]+coordsT[2][2]+coordsT[2][3]) ])
        msg("...finished. Bounding box: %s" % bbox, debugLevel=10)

        if len(self.elNodes):
            avgSize /= len(self.elNodes)
        # in method pointIsInside() tolerance is compared to
        # dot(cross(...)...) of vectors of avgSize magnitude
        self.tolerance = (1E-4 * avgSize)**3

        msg("Preparing KDTree...", debugLevel=10)
        self.tree = KDTree(points, leafsize=100, dim=6)
        self.labels = labels
        self.centroids = centroids
        self.point_search_initialized_6DTree = True
        msg("...finished.", debugLevel=10)


    def initializePointSearchInBox_6DTree(
            self, box=None, smoothRadius=0.0, smoothFactor=0.0):
        """
        Deprecated, for compatibility only, identical to
        initializePointSearch

        @param box: no function, for compatibility only
        @param smoothRadius: no function, for compatibility only
        @param smoothFactor: no function, for compatibility only
        """
        self.initializePointSearch()

    def pointIsInside_6DTree(
            self, point, returnVolElem=False, returnLinShapeFunc=False):
        """test if point is inside the volume
        point is tuple of coordinates

        If returnVolElem==True the function returns the element number of the
        tet element the point is in or False if the point is not inside the
        volume at all. Otherwise it just returns True or False.

        If returnLinShapeFunc==True the function returns a tuple
        (element number, [SF1, SF2, SF3, SF4]) SF1...SF4 being the values of
        the (linear tetrahedron) shape function associated with the four nodes
        of the tet at the given point. Each of the SF1...SF4 values lies
        between 0 and 1. They can be used to linearily interpolate values at
        the corner nodes. Smoothmode must be off (that's the case unless you
        specify a nonzero smooth_findradius argument to the function
        initializePointSearch().

        @Note: If the volume is created from elements other then tets which
        have been converted to tets, this method returns the id of the
        converted tet.

        @Note: If there is an element number of zero, the typical term
        "if vol.pointIsInside(pt, True)" may yield unexpected behaviour.
        Instead you could test like that:
           >>> if vol.pointIsInside(pt, True)!=False:

        Warning: Don't test the type of the return value with the built in
        isinstance() function. The type bool is a subclass of int!
        """

        if len(self.elNodes)==0:
            return False

        if not self.point_search_initialized_6DTree:
            self.initializePointSearch_6DTree()

        pL = self.tree.pointList  # abbreviation
        cL = self.centroids       # abbreviation
        pp = [x+self.tolerance for x in
              point[0], point[1], point[2], -point[0], -point[1], -point[2]]
        enclosingBoxIds = self.tree.searchBelow(pp)

        normDistIdxList = sorted(
            (sum(((point[i]-cL[idx][i])/(-pL[idx][i+3]-pL[idx][i]))**2
                 for i in (0,1,2) ),    idx )
            for idx in enclosingBoxIds )

        # full check if point inside any tet of elems_in_cell
        inside=False

        for [l,idx] in normDistIdxList:

            element = self.labels[idx]

            # coordinates of the four corner nodes
            x = [self.nodeCoords[node] for node in self.elNodes[element]]
            # coordinates relative to the first node
            # dp = [pi - x0i for pi, x0i in zip(point, x[0])]
            # dx = [[xxi-x0i for xxi, x0i in zip(xx, x[0])]
            #       for xx in x[1:]]
            dp = [point[0]-x[0][0], point[1]-x[0][1], point[2]-x[0][2]]
            dx = [[x[1][0]-x[0][0], x[1][1]-x[0][1], x[1][2]-x[0][2]],
                  [x[2][0]-x[0][0], x[2][1]-x[0][1], x[2][2]-x[0][2]],
                  [x[3][0]-x[0][0], x[3][1]-x[0][1], x[3][2]-x[0][2]]]
            # determinant of the "Jakobian"
            # cylcePos = [[0,1,2], [1,2,0], [2,0,1]]
            # detJ = sum(
            #     dx[0][ip[0]]*dx[1][ip[1]]*dx[2][ip[2]]
            #     - dx[2][ip[0]]*dx[1][ip[1]]*dx[0][ip[2]]
            #     for ip in cylcePos)
            detJ = (dx[0][0]*dx[1][1]*dx[2][2]
                    + dx[0][1]*dx[1][2]*dx[2][0]
                    + dx[0][2]*dx[1][0]*dx[2][1]
                    - dx[2][0]*dx[1][1]*dx[0][2]
                    - dx[2][1]*dx[1][2]*dx[0][0]
                    - dx[2][2]*dx[1][0]*dx[0][1])
            if detJ <= 0.0:
                # negative or zero volume, ignore this element
                continue

            # inverse * determinant
            # invJ = [ [a[ip2[1]]*b[ip2[2]] - a[ip2[2]]*b[ip2[1]] for ip2 in cylcePos]
            #          for a, b in ((dx[ip1[1]], dx[ip1[2]]) for ip1 in cylcePos)]
            invJ = [[dx[1][1]*dx[2][2] - dx[1][2]*dx[2][1],
                     dx[1][2]*dx[2][0] - dx[1][0]*dx[2][2],
                     dx[1][0]*dx[2][1] - dx[1][1]*dx[2][0]],
                    [dx[2][1]*dx[0][2] - dx[2][2]*dx[0][1],
                     dx[2][2]*dx[0][0] - dx[2][0]*dx[0][2],
                     dx[2][0]*dx[0][1] - dx[2][1]*dx[0][0]],
                    [dx[0][1]*dx[1][2] - dx[0][2]*dx[1][1],
                     dx[0][2]*dx[1][0] - dx[0][0]*dx[1][2],
                     dx[0][0]*dx[1][1] - dx[0][1]*dx[1][0]]]
            # parametric coordinates: xi_p := invJ * dp
            # xip = [float(sum(invJ[i][j]*dp[j] for j in range(3)))/detJ
            #        for i in range(3)]
            xip = [float(invJ[0][0]*dp[0]+invJ[0][1]*dp[1]+invJ[0][2]*dp[2])/detJ,
                   float(invJ[1][0]*dp[0]+invJ[1][1]*dp[1]+invJ[1][2]*dp[2])/detJ,
                   float(invJ[2][0]*dp[0]+invJ[2][1]*dp[1]+invJ[2][2]*dp[2])/detJ]

            shapeFuncs = [1-sum(xip), xip[0], xip[1], xip[2]]
            if all(xi>-self.tolerance for xi in shapeFuncs):
                inside = (element, shapeFuncs)
                break

        if not(returnLinShapeFunc) and inside is not False:
            if returnVolElem:
                inside = inside[0]
            else:
                inside = True

        return inside


    ## ------------------------------------------------------------------------
    ## ------------------------------------------------------------------------
    ## transforming the volume

    def scale(self, scale_factor, scale_origin):
        """Scale the volume by a certain factor relative to scale_origin

        scale_factor may be a single number in which case it is applied to all
        three dimensions or a list of three numbers with one factor for each
        dimension.

        Example:

        >>> # scale this volume by 2 in each direction (volume := volume*8)
        >>> v.scale(scale_factor=2, scale_origin=[0,0,0])

        Scaling is done by simply moving the point coordinates for this
        volume. Other volumes or the original mesh from which the volume
        originates are not affected.

        This creates a new node-dictionary, which is refered to from now on.
        It contains only those nodes connected to any element in self.elNodes.

        After scaling point search has to be initialized again,
        self.point_search_initialized is set to False.

        Possible improvement: To avoid initializing the point search again
        after scaling it would be necessary to also scale self.boundingBox and
        self.cellSize. In case of scale_factor being a list (different for the
        three dimensions) self.cellSize would have to become a list, too.
        """

        all_nodes = set()
        for elemNodes in self.elNodes.itervalues():
            all_nodes.update(elemNodes)

        if type(scale_factor) in (float, int):
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
        self.point_search_initialized_6DTree = False
        self.point_search_initialized_cellGrid = False


    def translate(self, move_vector):
        """Move the volume in space by move_vector

        Example:

        >>> # move this volume up by 10
        >>> vol.translate(move_vector=[0,0,10])

        Translating is done by simply moving the point coordinates for this
        volume. Other volumes or the original mesh from which the volume
        originates are not affected.

        This creates a new node-dictionary, which is refered to from now on.
        It contains only those nodes connected to any element in self.elNodes.

        The point search initialization is modified and doesn't have to be done
        again.
        """

        all_nodes = set()
        for elemNodes in self.elNodes.itervalues():
            all_nodes.update(elemNodes)

        new_coords = dict()
        for node_id in all_nodes:
            new_coords[node_id] = vector_plus(
                self.nodeCoords[node_id],move_vector)

        self.nodeCoords = new_coords

        # after translating point search has to be initialized again
        # this could be changed: just translate all coordinates and section
        # points in the KDTree as well...
        self.point_search_initialized_6DTree = False
        self.point_search_initialized_cellGrid = False


    ## ------------------------------------------------------------------------
    ## some service functions: iterator function

    @staticmethod
    def tetFacesIter(nodesOfThisTet):
        """Returns an iterator of the faces of the given tet element
        nodesOfThisTet gives the node ids, it is converted to a tuple

        use as in

        >>> myTet = mdb.elNodes[5731]  # list of node ids, e.g. [3, 43, 12, 36]
        >>> for faces in tetFacesIter(myTet):
        ...     print faces
        ...
        frozenset((3, 43, 12))
        frozenset((43, 3, 36))
        ...

        @Note: Deprecated, use self.orientedFacesIter()
        """

        # An extra comment on the index ordering used in the for loop below: If
        # the points order would stay that of the used tuple, each faces normal
        # points inwardly.
        nodesOfThisTet = tuple(nodesOfThisTet)
        for i1,i2,i3 in ((0,1,2), (1,0,3), (2,1,3), (0,2,3)):
            yield frozenset((nodesOfThisTet[i1], nodesOfThisTet[i2],
                             nodesOfThisTet[i3]))

    @staticmethod
    def orientedFacesIter(nodesOfThisTet):
        """Returns an iterator of the oriented faces of the given tet element.

        use as in

        >>> myTet = mdb.elNodes[5731]  # list of node ids, e.g. [3, 43, 12, 36]
        >>> for faces in tetFacesIter(myTet):
        ...     print faces
        ...
        (3, 12, 43)
        (3, 43, 36)
        ...

        @param nodesOfThisTet: list of ids of the four corner nodes

        @returns: An iterator that yields three-node-number-tuples. The nodes
        in each of those triples are ordered anticlockwise when looking from
        outside into the tet. The normal vector associated with this sense of
        rotation points outwards.
        """

        # An extra comment on the index ordering used in the for loop below: If
        # the points order would stay that of the used tuple, each faces normal
        # points inwardly.
        nodesOfThisTet = tuple(nodesOfThisTet)
        for i1,i2,i3 in ((0,2,1), (0,1,3), (1,2,3), (0,3,2)):
            yield (nodesOfThisTet[i1], nodesOfThisTet[i2],
                   nodesOfThisTet[i3])

# default search method: old method
VolumeFromMesh.setPointSearchMethod(method="cellGrid")

# for compatibility reasons provide the name VolumeFromTetMesh
VolumeFromTetMesh = VolumeFromMesh


###################################################################
###################################################################
###################################################################
###################################################################


class VolumeExtrudeOutline(Volume):
    """Volume class representing a volume between two horizontal planes with an
    outline in the x-y-plane defined by projecting a surface on the x-y-plane.

    >>> from bae.volume_02 import VolumeExtrudeOutline
    >>> vol = VolumeExtrudeOutline(surface=mysurf, zminmax=[0, 10])

    Volumes of this class are generated the following way: The given surface is
    projected on the x-y-plane and extruded in z direction from zminmx[0] to
    zminmax[1].

    This volume object keeps a reference of the triangle nodes and coordinates
    of the surface. So don't modify the surface.
    """

    # Todo: Should store or copy the surface object given in the constructor and
    # use its pointIsBelow(pt, returnElems=False) function for pointIsInside().
    # Remove VolumeExtrudeOutline.initializePointSearch(). Adapt self.translate().

    cellSizeFactor=0.5      # cell Size for array relative to 'avgSize'

    def __init__(self, surface, zminmax=[0,0], logfile=None):
        """
        @param surface: A L{bae.surface_02.TriangleSurface} object.
              L{bae.surface_03.TriangleSurface} should be possible too.

        @param zminmax: A tuple: the zmin and zmax coordinates. May be
              specified or modified later.

        @param logfile: deprecated argument, only for compatibility reasons.
              No effect.
        """
        Volume.__init__(self)

        if surface.__module__.endswith("surface_02"):
            self.triNodes = surface.triNodes
        elif surface.__module__.endswith("surface_03"):
            self.triNodes = surface.elNodes

        self.nodeCoords = surface.nodeCoords
        self.zminmax = list(zminmax)

        # internal variable, use getBoundigBox method! (initialize it)
        self.boundingBox = None

        # variables used for finding points in the volume already set?
        self.point_search_initialized = False
        return


    def setZminmax(self, *zminmax):
        """change z-range of the volume
        accepts two values: vol.setZminmax(zmin, zmax)
        or a list or tuple: vol.setZminmax([zmin, zmax])
        """
        if len(zminmax)==1:
            self.zminmax = map(float, zminmax[0])
        elif len(zminmax)==2:
            self.zminmax = map(float, zminmax)
        else:
            raise ValueError("setZminmax() expects two values or one tuple!")


    ## ------------------------------------------------------------------------
    ## some info

    def getBoundingBox(self):
        """Bounding Box of the whole volume.

        Example:

        >>> # bounding box of that volume
        >>> bb = v.getBoundingBox()
        >>> print "lower left front corner:", bb[0]
        >>> print "upper right back corner:", bb[1]

        @Returns: a list of two coordinate tuples (rather lists)
        The first coord tuple states the min value, the last the max.
        I.e. self.getBoundingBox()[i][j]:
        i \in {0,1} for min/max value, j \in {0,1,2} coordinate index
        """
        if not self.boundingBox:
            # initialize bounding box
            self.boundingBox=BoundingBox(dim=2)

            self.boundingBox.update(
                nodeCoords=self.nodeCoords, elNodes=self.elNodes)

        # copy 2D bounding box
        boundingBox = BoundingBox(list(xy) for xy in self.boundingBox)
        boundingBox[0].append(self.zminmax[0])  # min z
        boundingBox[1].append(self.zminmax[1])  # max z
        return boundingBox


    ## ------------------------------------------------------------------------
    ## methods for testing if a point lies inside the volume or outside

    def initializePointSearch(self, smooth_findradius=0.0):
        """call this function once before pointIsInside
        (it is done automatically)

        smoothmode: Additionally consider points to be inside the volume that
        are in fact outside the volume but within the smooth_findradius to the
        closest centre point of any of the elements of the volume.
        """

        # list of centroid-coord-tuples
        self.triCentroids=dict()

        # list of [[xmin,ymin],[xmax,ymax]] for each tri element
        elementBox=dict()

        self.boundingBox = None
        avgSize=0

        # size of the projection on the x-z plane
        # pos facing one direction, neg facing the other
        elSignedSize = dict()
        for element,nodes in self.triNodes.iteritems():

            node_pts = [self.nodeCoords[node]
                        for node in nodes]
            v12=vector(node_pts[0], node_pts[1])
            v13=vector(node_pts[0], node_pts[2])
            detA = v12[0]*v13[1] - v12[1]*v13[0]
            elSignedSize[element] = 0.5*detA

            xy=[[],[]]
            x=y=0
            for coords in node_pts:
                for i in xrange(2):  # for min max
                    xy[i].append(coords[i])
                x=x+coords[0]/3.0  # for centroid
                y=y+coords[1]/3.0
            self.triCentroids[element]=[x,y]

            elementBox[element]=[[min(xy[0]),min(xy[1])],
                                 [max(xy[0]),max(xy[1])]]

            # update bounding box
            if self.boundingBox:
                for i in xrange(2):  # x,y
                    if elementBox[element][0][i]<self.boundingBox[0][i]:  # min
                        self.boundingBox[0][i]=elementBox[element][0][i]
                    if elementBox[element][1][i]>self.boundingBox[1][i]:  # max
                        self.boundingBox[1][i]=elementBox[element][1][i]
            else:
                self.boundingBox=[[0,0],[0,0]]
                for i in xrange(2):  # x,y
                    self.boundingBox[0][i]=elementBox[element][0][i]
                    self.boundingBox[1][i]=elementBox[element][1][i]

            avgSize+=self.dist2D(
                elementBox[element][1],elementBox[element][0])

        avgSize /= len(self.triNodes)

        # in method pointIsInside() tolerance is compared to values 0..1
        self.tolerance = 1E-4

        msg('bounding box of the surface: [(%.1f,%.1f),(%.1f,%.1f)]'
            % (self.boundingBox[0][0],self.boundingBox[0][1],
               self.boundingBox[1][0],self.boundingBox[1][1]))
        msg('avg. element size: %g' % avgSize)


        ## -----------------------------------------------------------------
        ## map elements to cells

        self.cellOrigin=self.boundingBox[0]
        cellDimensions=[0,0]
        self.cellSize=avgSize*self.cellSizeFactor
        self.cellNumber=[0,0]

        # x,y,z number of cells
        for i in xrange(2):
            cellDimensions[i]=self.boundingBox[1][i]-self.boundingBox[0][i]
        msg('cell array dimensions: %s' % (" ".join(map(str, cellDimensions))))

        for i in xrange(2):
            self.cellNumber[i]=int(cellDimensions[i]/self.cellSize)+1
            msg(" "+str(self.cellNumber[i]))
        msg('cell numbers: %s cells: %d'
            % (" ".join(map(str, self.cellNumber)),
               self.cellNumber[0]*self.cellNumber[1]))

        # create empty array
        for i in xrange(self.cellNumber[0]):
            if i==0: self.cells=[]
            for j in xrange(self.cellNumber[1]):
                if j==0: LJ=[]
                LJ.append([])
            self.cells.append(LJ)

        # tris with a too small projected size will be ignored
        zerosize = (self.tolerance*avgSize*0.1)**2

        # map elements
        for element in self.triNodes:

            # ignore too small tris
            if abs(elSignedSize[element]) < zerosize:
                continue

            # get index ranges
            cellIndex=[[0,0],[0,0]]  # low,high
            for l in xrange(2):  # low,high - index
                for i in xrange(2):  # dimension (x,y) - index
                    cellIndex[l][i]=int((elementBox[element][l][i]
                                         -self.cellOrigin[i])/self.cellSize)
            # append element ID to cells that overlap the element's bounding box
            for i in xrange(cellIndex[0][0],cellIndex[1][0]+1):
                for j in xrange(cellIndex[0][1],cellIndex[1][1]+1):
                    self.cells[i][j].append(element)

        # for smooth mode store radius
        self.smoothmode = smooth_findradius>0
        self.smooth_findradius = float(smooth_findradius)

        # finished initializing
        self.point_search_initialized = True
        return

    def pointIsInside(self, point):
        """test if point is inside the volume
        point is tuple of coordinates
        """

        if len(self.triNodes)==0:
            return False

        if not self.point_search_initialized:
            self.initializePointSearch()

        if (point[2]<self.zminmax[0]
            or point[2]>self.zminmax[1]):
            # point outside z-region
            return False

        # cell check (includes bounding box test from cell index)
        index=[0,0]
        for i in xrange(2):
            index[i]=int((point[i]-self.cellOrigin[i])/self.cellSize)
            if index[i]<0 or index[i]>=self.cellNumber[i]:
                return False

        # if cell check passed ... go on
        elems_in_cell=self.cells[index[0]][index[1]]
        centroid_to_point_distance=[]
        # sort by distance
        for element in elems_in_cell:
            l=self.dist2D(point,self.triCentroids[element])
            centroid_to_point_distance.append([l,element])
        centroid_to_point_distance.sort()

        # full check if point inside any tri of elems_in_cell
        inside=0

        if self.smoothmode:
            if (centroid_to_point_distance and
                (centroid_to_point_distance[0][0]<self.smooth_findradius)):
                inside = 1

        for [l,element] in centroid_to_point_distance:
            node_pts = [self.nodeCoords[node]
                        for node in self.triNodes[element]]
            # vectors
            v12=vector(node_pts[0], node_pts[1])
            v13=vector(node_pts[0], node_pts[2])
            detA = v12[0]*v13[1] - v12[1]*v13[0]
            adjA = [[v13[1], -v13[0]], [-v12[1], v12[0]]]
            invA = [[x/detA for x in row] for row in adjA]

            v1p=vector(node_pts[0], point)
            a0 = invA[0][0]*v1p[0] + invA[0][1]*v1p[1]
            a1 = invA[1][0]*v1p[0] + invA[1][1]*v1p[1]

#             if a0>=0 and a1>=0 and (a0+a1)<=1:
#                 inside = 1
#                 break

            if ((a0 > -self.tolerance) and
                (a1 > -self.tolerance) and
                ((a0+a1) < (1+self.tolerance))):
                inside = 1
                break

        return inside


    ## ------------------------------------------------------------------------
    ## ------------------------------------------------------------------------
    ## transforming the volume

    def translate(self, move_vector):
        """Move the volume in space by move_vector

        Example:

        >>> # move this volume up by 10
        >>> vol.translate(move_vector=[0,0,10])

        Translating is done by simply moving the point coordinates for this
        volume. Other volumes or the original mesh from which the volume
        originates are not affected.

        This creates a new node-dictionary, which is refered to from now on.
        It contains only those nodes connected to any element in self.elNodes.

        The point search initialization is modified and doesn't have to be done
        again.
        """

        all_nodes = set()
        for nodes in self.triNodes.itervalues():
            all_nodes.update(nodes)

        new_coords = dict()
        for node_id in all_nodes:
            new_coords[node_id] = vector_plus(
                self.nodeCoords[node_id],move_vector)

        self.nodeCoords = new_coords
        for i in xrange(2):  # min, max
            self.zminmax[i] += move_vector[2]


        # adapt point search properties
        if self.point_search_initialized and len(self.triNodes)>0:

            # tet Centroids
            for centroid in self.triCentroids.itervalues():
                vector_modif_add(centroid, move_vector)

            # BoundingBox
            for i in xrange(2):  # x,y
                self.boundingBox[0][i] += move_vector[i]
                self.boundingBox[1][i] += move_vector[i]

            # search cells: nothing to modify here
            # self.cellOrigin is just an alias for self.boundingBox[0]


    ## ------------------------------------------------------------------------
    ## some service functions: 2D distance

    @staticmethod
    def dist2D(a,b):
        return sqrt((b[0]-a[0])**2 + (b[1]-a[1])**2)


###################################################################

class VolumeExtrudeSurfVert(Volume):
    """
    Volume class representing a volume vertically extruded from a given
    triangle surface. The lower boundary of the volume is defined by moving the
    given surface vertically by zRelMinMax[0]. The upper boundary of the volume
    is defined by moving the given surface vertically by zRelMinMax[1].

    >>> from bae.volume_02 import VolumeBetweenSurfaceAndZPlane
    >>> vol = VolumeBetweenSurfaceAndZPlane(surface=mysurf, zplane=0.0,
    ...                                     volumePosition="below")

    This volume object keeps a reference of the triangle nodes and coordinates
    of the surface. So don't modify the surface or the volume will be modified
    too. Note that if you translate or scale the volume, the surface will also
    be changed accordingly!
    """

    def __init__(self, surface, zRelMinMax=[0,0], logfile=None):
        """
        @param surface: A L{bae.surface_02.TriangleSurface} object.
          L{bae.surface_03.TriangleSurface} should be possible too.

        @param zRelMinMax: A tuple specifying the offset of the lower and upper
          boundary of the volume to the given surface. May be specified or
          modified later.

        @param logfile: deprecated argument, only for compatibility reasons.
              No effect.
        """
        Volume.__init__(self)

        self.surface = surface
        self.zRelMinMax = list(zRelMinMax)

        return


    ## ------------------------------------------------------------------------
    ## some info

    def getBoundingBox(self):
        """Bounding Box of the whole volume.

        Example:

        >>> # bounding box of that volume
        >>> bb = vol.getBoundingBox()
        >>> print "lower left front corner:", bb[0]
        >>> print "upper right back corner:", bb[1]

        @Returns: a list of two coordinate tuples (rather lists)
        The first coord tuple states the min value, the last the max.
        I.e. self.getBoundingBox()[i][j]:
        i \in {0,1} for min/max value, j \in {0,1,2} coordinate index
        """
        # copy the bounding box of the surface: to be cautious, don't modify it
        boundingBox = BoundingBox(
            [list(xyz) for xyz in self.surface.getBoundingBox()])
        boundingBox[0][2] += self.zRelMinMax[0]
        boundingBox[1][2] += self.zRelMinMax[1]
        return boundingBox


    ## ------------------------------------------------------------------------
    ## methods for testing if a point lies inside the volume or outside

    def initializePointSearch(self, smooth_findradius=0.0):
        """
        Usually you don't need to call this function.

        You may call this function once before pointIsInside, otherwise it is
        done automatically.

        smoothmode: Additionally consider points to be inside the volume that
        are in fact outside the volume but within the smooth_findradius to the
        closest centre point of any of the elements of the volume.
        """
        self.surface.initializePointSearch(smooth_findradius)
        return

    def pointIsInside(self, point):
        """test if point is inside the volume
        point is tuple of coordinates
        """

        try:
            if len(self.surface.triNodes)==0:
                # no surface (bae.surface_02.TriangleSurface object)
                return False
        except AttributeError:
            if len(self.surface.elNodes)==0:
                # no surface (bae.surface_03.TriangleSurface object)
                return False


        # list of projected points
        # if parts of the surface lie on top of each other there may be more
        # than one projection for each point
        projectedPts = self.surface.projectPoints([point])[0]
        if len(projectedPts)==0:
            return False

        # check if any of the projections is close enough to the actual point
        for pt in projectedPts:
            zRel = point[2] - pt[2]
            if zRel>=self.zRelMinMax[0] and zRel<=self.zRelMinMax[1]:
                return True
        return False


    ## ------------------------------------------------------------------------
    ## ------------------------------------------------------------------------
    ## transforming the volume

    def translate(self, move_vector):
        """Move the volume in space by move_vector

        Example:

        >>> # move this volume up by 10
        >>> vol.translate(move_vector=[0,0,10])

        Translating is done by simply moving the point coordinates for this
        volume. The surface that was used to create this volume is moved as
        well.

        The point search initialization is modified and doesn't have to be done
        again.
        """

        self.surface.translate(move_vector)


    def scale(self, scale_factor, scale_origin):
        """Scale the volume by a certain factor relative to scale_origin

        scale_factor may be a single number in which case it is applied to all
        three dimensions or a list of three numbers with one factor for each
        dimension.

        Example:

        >>> # scale this volume by 2 in each direction
        >>> vol.scale(scale_factor=2, scale_origin=[0,0,0])

        Scaling is done by simply moving the point coordinates for this
        volume. The surface that was used to create this volume is scaled as
        well.

        After scaling point search has to be initialized again.
        """

        if isinstance(scale_factor, (float, int)):
            scale_factor = [scale_factor,scale_factor,scale_factor]

        self.surface.scale(scale_factor, scale_origin)
        self.zRelMinMax[0] *= scale_factor[2]
        self.zRelMinMax[1] *= scale_factor[2]


###################################################################

class VolumeBetweenSurfaceAndZPlane(Volume):
    """
    Volume class representing a volume between an arbitrary surface made of
    triangles and a horizontal plane at a given z coordinate.

    >>> from bae.volume_02 import VolumeBetweenSurfaceAndZPlane
    >>> vol = VolumeBetweenSurfaceAndZPlane(surface=mysurf, zplane=0.0,
    ...                                   volumePosition="below")

    This volume object keeps a reference of the triangle nodes and coordinates
    of the surface. So don't modify the surface or the volume will be modified
    too. Note that if you translate or scale the volume, the surface will also
    be changed accordingly!

    In case of surface parts lying above each other: The bounding surface
    of the volume is taken as the set of those points on the given triangle
    surface that are furthest away from the given zPlane into the direction
    indicated by the parameter volumePosition. (Try to say that in a simpler
    way!)

    @ivar zplane: The z coordinate of the bounding horizontal plane.
    @ivar surfOutsideNormalZ: 1 if the volume is below the surface, -1 if the
      volume is above the surface
    """

    def __init__(self, surface, zplane=0.0, volumePosition=None,
                 logfile=None):
        """
        @param surface: A surface_02.TriangleSurface object.
          L{bae.surface_03.TriangleSurface} should be possible too.

        @param zplane: The z coordinate of the bounding plane. You may change
          this later by just changing self.zplane.

        @param volumePosition: Either "below" or "above" the surface. E.g. if
          the surface represents a pit shape the excavated volume would be
          "above". If not specified, it is determined automatically: The volume
          is considered to be below the surface if most of its nodes lie above
          the given zplane. The volume is considered to be above the surface if
          most of its nodes lie below the given zplane.

        @param logfile: deprecated argument, only for compatibility reasons.
              No effect.
        """
        Volume.__init__(self)

        # determine self.surfOutsideNormalZ
        if volumePosition and volumePosition.lower() == "below":
            self.surfOutsideNormalZ = 1
        elif volumePosition and volumePosition.lower() == "above":
            self.surfOutsideNormalZ = -1
        else:
            avgSurfZ = 0.0

            if surface.__module__.endswith("surface_02"):
                allNodes = set(surface.getNodeTris())
            elif surface.__module__.endswith("surface_03"):
                allNodes = set(surface.getNodeElems())

            nodeCoords = surface.nodeCoords
            for node in allNodes:
                avgSurfZ += nodeCoords[node][2]
            avgSurfZ /= len(allNodes)
            if avgSurfZ>zplane:
                self.surfOutsideNormalZ = 1
                posToSurf, posToZPlane = "below", "above"
            else:
                self.surfOutsideNormalZ = -1
                posToSurf, posToZPlane = "above", "below"
            msg("The volume of class VolumeBetweenSurfaceAndZPlane is"
                " considered to be %s the surface at about z=%g and %s the"
                " zplane at z=%g."
                % (posToSurf, avgSurfZ, posToZPlane, zplane))
            del avgSurfZ, allNodes, nodeCoords, node, posToSurf, posToZPlane

        self.surface = surface
        self.zplane = zplane
        return


    ## ------------------------------------------------------------------------
    ## some info

    def getBoundingBox(self):
        """Bounding Box of the whole volume.

        Example:

        >>> # bounding box of that volume
        >>> bb = vol.getBoundingBox()
        >>> print "lower left front corner:", bb[0]
        >>> print "upper right back corner:", bb[1]

        @Returns: a list of two coordinate tuples (rather lists)
        The first coord tuple states the min value, the last the max.
        I.e. self.getBoundingBox()[i][j]:
        i \in {0,1} for min/max value, j \in {0,1,2} coordinate index
        """
        # copy the bounding box of the surface: to be cautious, don't modify it
        boundingBox = BoundingBox()
        boundingBox[:] = self.surface.getBoundingBox()
        if self.surfOutsideNormalZ>0:
            boundingBox[0][2] = self.zplane
        else:
            boundingBox[1][2] = self.zplane
        return boundingBox


    ## ------------------------------------------------------------------------
    ## methods for testing if a point lies inside the volume or outside

    def initializePointSearch(self, smooth_findradius=0.0):
        """
        Usually you don't need to call this function.

        You may call this function once before pointIsInside, otherwise it is
        done automatically.

        smoothmode: Additionally consider points to be inside the volume that
        are in fact outside the volume but within the smooth_findradius to the
        closest centre point of any of the elements of the volume.
        """
        self.surface.initializePointSearch(smooth_findradius)
        return

    def pointIsInside(self, point):
        """test if point is inside the volume
        point is tuple of coordinates
        """

        try:
            if len(self.surface.triNodes)==0:
                # no surface (bae.surface_02.TriangleSurface object)
                return False
        except AttributeError:
            if len(self.surface.elNodes)==0:
                # no surface (bae.surface_03.TriangleSurface object)
                return False

        # pointDistAboveZPlane: the distance from the zplane to the given point
        # in the direction of self.surfOutsideNormalZ
        # (from the zplane into the volume == from the surface out of the vol.)
        pointDistAboveZPlane = (point[2]-self.zplane) * self.surfOutsideNormalZ
        if pointDistAboveZPlane<0:
            # point on the wrong side of self.zplane
            return False

        # list of projected all points
        projectedPts = self.surface.projectPoints([point])[0]
        if len(projectedPts)==0:
            return False

        # Of those projected points find the point farest above/below the
        # zplane. Store the greatest projPtDistAboveZPlane: the distance from
        # the zplane to the projected point in the direction of
        # self.surfOutsideNormalZ
        # (from the zplane into the volume == from the surface out of the vol.)
        projPtDistAboveZPlane = max(
            (pt[2]-self.zplane) * self.surfOutsideNormalZ
            for pt in projectedPts )
        return (pointDistAboveZPlane <= projPtDistAboveZPlane)

    def intersectionGrid(self, grid, fieldName="insideVolume"):
        """Find points of the regular grid inside self.

        Example: Get a field indicating points in or outside the volume:
         >>> from bae.mesh_01 import MeshStructuredPoints
         >>> from bae.volume_02 import VolumeInSurfCheckRays
         >>> from bae.surface_03 import TriangleSurface
         >>> from bae.vtk_02 import FieldsOnStructPointGrid as Vtk
         >>>
         >>> surf = TriangleSurface.fromSTL("myVolume.stl")
         >>> vol = VolumeInSurfCheckRays(surf)
         >>> grid = MeshStructuredPoints(....)
         >>> insideField = vol.intersectionGrid(grid)
         >>>
         >>> vtk = Vtk(mesh=grid)
         >>> vtk.updateFields(insideField)
         >>> vtk.toVtk("inVolume.vtk")

        @param grid: L{bae.mesh_01.MeshStructuredPoints}-object.
            A L{MeshStructuredPointsRot} object works as well if the rotation
            is about the z-axis only. Otherwise use the general method
            L{Volume.intersectionGrid}.
        @param fieldName: fieldName attribute of the field to be returned

        @returns: L{bae.field_01.Field}-object (position="structPt",
            dataType="scalar") containing 1 for "in", 0 for "out"
        """
        # counting pos z in the direction from the zplane to the surface
        # get *this* z-coordinate of the zplane
        zPlaneZPlane = self.zplane * self.surfOutsideNormalZ

        result = Field.classFromPosType(fieldName, "structPt", "scalar")(
            [0]*len(grid))

        points = [
            grid.getPoint((i,j,0))
            for j in range(grid.gridPtNb[1])
            for i in range(grid.gridPtNb[0])]

        intersections = self.surface.projectPoints(points)

        ticker = MsgTicker("...processing line %%d/%d"
                           % grid.gridPtNb[0]*grid.gridPtNb[1])
        for j in range(grid.gridPtNb[1]):
            i0 = j*grid.gridPtNb[0]
            for i in range(grid.gridPtNb[0]):
                ticker.tick()
                idxBase = i*grid.strides[0]+j*grid.strides[1]
                projectedPts = intersections[i+i0]

                if len(projectedPts)==0:
                    # all pts out if no intersection
                    for k in range(grid.gridPtNb[2]):
                        result[idxBase+k*grid.strides[2]] = 0
                    continue

                # counting pos z in the direction from the zplane to the surface
                # get outmost surface intersection
                projPtZPlane = max(pt[2]*self.surfOutsideNormalZ
                                   for pt in projectedPts)
                if projPtZPlane < zPlaneZPlane:
                    # surface "below" zplane at this x-y-location, all pts out
                    for k in range(grid.gridPtNb[2]):
                        result[idxBase+k*grid.strides[2]] = 0
                    continue

                if self.surfOutsideNormalZ>0:
                    z0 = self.zplane
                    z1 = projPtZPlane
                else:
                    z0 = -projPtZPlane
                    z1 = self.zplane

                for k in grid.getIdsInIntervalOverZ(z0, z1):
                    result[idxBase+k*grid.strides[2]] = 1

        del ticker
        return result


    ## ------------------------------------------------------------------------
    ## ------------------------------------------------------------------------
    ## transforming the volume

    def translate(self, move_vector):
        """Move the volume in space by move_vector

        Example:

        >>> # move this volume up by 10
        >>> vol.translate(move_vector=[0,0,10])

        Translating is done by simply moving the point coordinates for this
        volume. The surface that was used to create this volume is moved as
        well.

        The point search initialization is modified and doesn't have to be done
        again.
        """

        self.surface.translate(move_vector)
        self.zplane += move_vector[2]


    def scale(self, scale_factor, scale_origin):
        """Scale the volume by a certain factor relative to scale_origin

        scale_factor may be a single number in which case it is applied to all
        three dimensions or a list of three numbers with one factor for each
        dimension.

        Example:

        >>> # scale this volume by 2 in each direction
        >>> vol.scale(scale_factor=2, scale_origin=[0,0,0])

        Scaling is done by simply moving the point coordinates for this
        volume. The surface that was used to create this volume is scaled as
        well.

        After scaling point search has to be initialized again.
        """

        if isinstance(scale_factor, (float, int)):
            scale_factor = [scale_factor,scale_factor,scale_factor]

        self.surface.scale(scale_factor, scale_origin)

        orig2new = (self.zplane - scale_origin[2]) * scale_factor[2]
        self.zplane = orig2new + scale_origin[2]


###################################################################

class VolumeInSurfZPrj(Volume):
    """
    Volume class representing a volume between the lowest and highest triangle
    of a given surface.

    Points are tested by projecting the point vertically to any appropriate
    triangle of the surface and then testing whether the point is between the
    highest and lowest of those projected points.

    This is suitable for volumes defined by their exteriour surface if those
    surfaces are convex in the sense that any vertical line has no more than
    two intersections with the boundary of the volume. (There may be additional
    internal parts of the defining surface that don't define a volume boundary.
    They are simply ignored.)

    >>> from bae.volume_02 import VolumeInSurfZPrj
    >>> vol = VolumeInSurfZPrj(surface=mysurf)

    This volume object keeps a reference of the triangle nodes and coordinates
    of the surface. So don't modify the surface or the volume will be modified
    too. Note that if you translate or scale the volume, the surface will also
    be changed accordingly!
    """

    def __init__(self, surface, logfile=None):
        """
        @param surface: A surface_02.TriangleSurface object.
              L{bae.surface_03.TriangleSurface} should be possible too.

        @param logfile: deprecated argument, only for compatibility reasons.
              No effect.
        """
        Volume.__init__(self)

        self.surface = surface

        return


    ## ------------------------------------------------------------------------
    ## some info

    def getBoundingBox(self):
        """Bounding Box of the whole volume.

        Example:

        >>> # bounding box of that volume
        >>> bb = vol.getBoundingBox()
        >>> print "lower left front corner:", bb[0]
        >>> print "upper right back corner:", bb[1]

        @Returns: a list of two coordinate tuples (rather lists)
        The first coord tuple states the min value, the last the max.
        I.e. self.getBoundingBox()[i][j]:
        i \in {0,1} for min/max value, j \in {0,1,2} coordinate index
        """
        # copy the bounding box of the surface: to be cautious, don't modify it
        boundingBox = BoundingBox(
            [list(xyz) for xyz in self.surface.getBoundingBox()])
        return boundingBox


    ## ------------------------------------------------------------------------
    ## methods for testing if a point lies inside the volume or outside

    def initializePointSearch(self, smooth_findradius=0.0):
        """
        Usually you don't need to call this function.

        You may call this function once before pointIsInside, otherwise it is
        done automatically.

        smoothmode: Additionally consider points to be inside the volume that
        are in fact outside the volume but within the smooth_findradius to the
        closest centre point of any of the elements of the volume.
        """
        self.surface.initializePointSearch(smooth_findradius)
        return

    def pointIsInside(self, point):
        """test if point is inside the volume
        point is tuple of coordinates
        """

        try:
            if len(self.surface.triNodes)==0:
                # no surface (bae.surface_02.TriangleSurface object)
                return False
        except AttributeError:
            if len(self.surface.elNodes)==0:
                # no surface (bae.surface_03.TriangleSurface object)
                return False

        # list of projected points
        # if parts of the surface lie on top of each other there may be more
        # than one projection for each point
        projectedPts = self.surface.projectPoints([point])[0]
        if len(projectedPts)==0:
            return False

        # check if any of the projections is close enough to the actual point
        projectedZs = [pt[2] for pt in projectedPts]
        return (min(projectedZs) <= point[2] <= max(projectedZs))


    ## ------------------------------------------------------------------------
    ## ------------------------------------------------------------------------
    ## transforming the volume

    def translate(self, move_vector):
        """Move the volume in space by move_vector

        Example:

        >>> # move this volume up by 10
        >>> vol.translate(move_vector=[0,0,10])

        Translating is done by simply moving the point coordinates for this
        volume. The surface that was used to create this volume is moved as
        well. (There is actually no difference between moving the initial
        surface and moving this volume.)

        The point search initialization is modified and doesn't have to be done
        again.
        """

        self.surface.translate(move_vector)


    def scale(self, scale_factor, scale_origin):
        """Scale the volume by a certain factor relative to scale_origin

        scale_factor may be a single number in which case it is applied to all
        three dimensions or a list of three numbers with one factor for each
        dimension.

        Example:

        >>> # scale this volume by 2 in each direction
        >>> vol.scale(scale_factor=2, scale_origin=[0,0,0])

        Scaling is done by simply moving the point coordinates for this
        volume. The surface that was used to create this volume is scaled as
        well. (There is actually no difference between scaling the initial
        surface and scaling this volume.)

        After scaling point search has to be initialized again.
        """

        self.surface.scale(scale_factor, scale_origin)


###################################################################

class VolumeInSurfCheckRays(Volume):
    """
    Volume class representing a volume enclosed by a given surface.

    Points are tested by casting rays from the point through the surface and
    then counting the number of intersections.

    The current implementation casts vertical rays.

    This is especially suitable for water tight closed surfaces. Holes in the
    connectivity (free edges) are no concern as long as no ray can fall through.

    Every point will be decided upon whether it's in or out. This may be wrong.
    The class maintains a list of points which are uncertain whether they are
    in or out. It's good practice to check the contents of that list.

    >>> from bae.volume_02 import VolumeInSurfCheckRays
    >>> vol = VolumeInSurfCheckRays(mysurf)
    >>> elsetInVol = vol.intersectionElset(mesh)
    >>> if vol.dubiousPoints:
    >>>     msg("WARNING: The assignments of %d centroids are questionable and"
    >>>         " insecure." % len(vol.dubiousPoints))

    This volume object keeps a reference of the triangle nodes and coordinates
    of the surface. So don't modify the surface or the volume will be modified
    too. Note that if you translate or scale the volume, the surface will also
    be changed accordingly!

    @ivar dubiousPoints: list of coordinates identifying those points which have
       been found hard to decide whether they are in or out. Every call to
       pointIsInside or intersectionGrid can add one point to this list.
       Clearing the list is done like that:
        >>> vol.dubiousPoints = []

       It's advisable to check this list after points have been evaluated.
    """

    def __init__(self, surface, projectionAxis=2):
        """
        @param surface: A L{bae.surface_02.TriangleSurface} or
              L{bae.surface_03.TriangleSurface} object.

        @param projectionAxis: normal axis index. axis=2 means vertical
           projection i.e. pointIsBelow actually checks if the vertical
           projection of a point is on the surface. axis=0/1 means pointIsBelow
           checks if the projection in x-/y-direction is on the surface.
        """
        Volume.__init__(self)
        self.surface = surface
        self.dubiousPoints = []
        self.projectionAxis = projectionAxis
        return


    ## ------------------------------------------------------------------------
    ## some info

    def getBoundingBox(self):
        """Bounding Box of the whole volume.

        Example:

        >>> # bounding box of that volume
        >>> bb = vol.getBoundingBox()
        >>> print "lower left front corner:", bb[0]
        >>> print "upper right back corner:", bb[1]

        @Returns: a list of two coordinate tuples (rather lists)
        The first coord tuple states the min value, the last the max.
        I.e. self.getBoundingBox()[i][j]:
        i \in {0,1} for min/max value, j \in {0,1,2} coordinate index
        """
        # copy the bounding box of the surface: to be cautious, don't modify it
        boundingBox = BoundingBox(
            [list(xyz) for xyz in self.surface.getBoundingBox()])
        return boundingBox


    ## ------------------------------------------------------------------------
    ## methods for testing if a point lies inside the volume or outside

    def pointIsInside(self, point):
        """Test if point is inside the volume by shooting a vertical ray
        (ignoring the projectionAxis argument to the constructor).

        point is a tuple of coordinates.

        Check self.dubiousPoints whether the decision seemed firm. This list
        will be updated with point if not. (self.dubiousPoints will not be
        cleared by this function and can be used to collect points from
        subsequent calls to this function.)
        """

        # list of projected points
        # if parts of the surface lie on top of each other there may be more
        # than one projection for each point
        projectedPts = self.surface.projectPoints([point])[0]
        projectedZs = sorted(pt[2] for pt in projectedPts)
        if len(projectedZs) == 0:
            return False

        # if number of intersections is one then point is out and dubious
        if len(projectedZs) == 1:
            self.dubiousPoints.append(point)
            return False

        # if two intersections are found then the result is clear/definite
        if len(projectedZs) == 2:
            return (projectedZs[0] <= point[2] <= projectedZs[-1])

        # if number of intersections is odd then this is dubious
        # The point is considered inside if it's between the first and
        # last intersection of the rays.
        if len(projectedZs) % 2 != 0:
            self.dubiousPoints.append(point)
            return (projectedZs[0] <= point[2] <= projectedZs[-1])

        # assign rest to in / out
        else:
            cntBelow = bisect.bisect_left(projectedZs, point[2])
            return (cntBelow % 2 != 0)

    def intersectionGrid(self, grid, fieldName="insideVolume"):
        """Find points of the regular grid inside self.

        This method casts rays parallel to the axis specified by the
        projectionAxis argument to the constructor. If this leaves "dubious"
        points, i.e. points that yield an odd number of intersections then
        those points are checked with a ray in a perpendicular direction.

        Example: Get a field indicating points in or outside the volume:
         >>> from bae.mesh_01 import MeshStructuredPoints
         >>> from bae.volume_02 import VolumeInSurfCheckRays
         >>> from bae.surface_03 import TriangleSurface
         >>> from bae.vtk_02 import FieldsOnStructPointGrid as Vtk
         >>>
         >>> surf = TriangleSurface.fromSTL("myVolume.stl")
         >>> vol = VolumeInSurfCheckRays(surf)
         >>> grid = MeshStructuredPoints(....)
         >>> insideField = vol.intersectionGrid(grid)
         >>>
         >>> vtk = Vtk(mesh=grid)
         >>> vtk.updateFields(insideField)
         >>> vtk.toVtk("inVolume.vtk")

        Additionally to the resulting field containing in (1) and out (0)
        values the instance variable L{dubiousPoints} will be set to contain
        grid indexes of grid points that have been found hard to decide
        whether they are in or out.

        @Note: The aim of this procedure is foremost to quickly return a result
        that possibly is complete. No effort is made to check if casting rays
        in different directions yields consistent results or if a contradiction
        occurs.

        @param grid: L{bae.mesh_01.MeshStructuredPoints}-object.
            A L{MeshStructuredPointsRot} object works as well if the rotation
            is about the z-axis only. Otherwise use the general method
            L{Volume.intersectionGrid}.
        @param fieldName: fieldName attribute of the field to be returned

        @returns: L{bae.field_01.Field}-object (position="structPt",
            dataType="scalar") containing 1 for "in", 0 for "out"
        """
        # result initialized with 0 (out) values
        result = Field.classFromPosType(fieldName, "structPt", "scalar")(
            [0]*len(grid))
        self.dubiousPoints = []

        axis = self.projectionAxis  # abbreviation

        # axis index for the first and second axes perpendicular to the
        # projection direction.
        # (defaults to 0, 1 for the default axis==2,
        # (or else [1,2] or [2,0] for axis==0 or 1)
        aix = [1, 2, 0, 1, 2][axis:axis+2]

        # determine points for which we'll search projections
        NN = grid.gridPtNb  # abbreviation
        if axis==0:
            gridIdIter = ((0,j,k) for k in range(NN[2]) for j in range(NN[1]))
        elif axis==1:
            gridIdIter = ((i,0,k) for k in range(NN[2]) for i in range(NN[0]))
        else:  # axis==2
            gridIdIter = ((i,j,0) for j in range(NN[1]) for i in range(NN[0]))
        points = [grid.getPoint(ijk) for ijk in gridIdIter]

        intersections = self.surface.projectPoints(points, axis=axis)

        ijkBase = [0,0,0]  # grid index of dubious points to be filled
        ticker = MsgTicker("...processing %s-ray %%d/%d"
                           % ("xyz"[axis], len(points)))
        for j in range(NN[aix[1]]):
            i0 = j*NN[aix[0]]
            for i in range(NN[aix[0]]):
                ticker.tick()
                interZ = sorted(pt[axis] for pt in intersections[i+i0])
                # # check only working for axis=2
                # basePt = grid.getPoint((i,j,0))
                # print ("for i,j=(%d,%d) x,y=(%5.1f,%5.1f) intersections %s"
                #        % (i,j, basePt[0], basePt[1],
                #           ", ".join("%5.1f" % x for x in interZ)))
                idxBase = i*grid.strides[aix[0]]+j*grid.strides[aix[1]]

                # if number of intersections is odd then points are dubious...
                if len(interZ) % 2 != 0:

                    # initialize grid index for dubious points for this ray
                    ijkBase[aix[0]] = i
                    ijkBase[aix[1]] = j

                    if len(interZ)==1:
                        # if there is only one intersection then the complete
                        # column will be assumed to be "out" and dubious
                        for k in range(NN[axis]):
                            ijkBase[axis] = k
                            self.dubiousPoints.append(list(ijkBase))

                    else:
                        # if there are three or more intersections (and the
                        # number is odd) then assume "in" and dubious between
                        # the first and last intersection
                        kIds = grid.getIdsInInterval(
                            interZ[0], interZ[-1], axis=axis)
                        for k in kIds:
                            result[idxBase + k*grid.strides[axis]] = 1
                            ijkBase[axis] = k
                            self.dubiousPoints.append(list(ijkBase))  # append a copy

                # rest is out of question, make assignments
                else:
                    if len(interZ)>1:
                        # everything in any of the intervals is "inside"
                        for z0, z1 in zip(*((iter(interZ),)*2)):
                            # print "... intersections z=%5.1f, %5.1f" % (z0,z1)
                            for k in grid.getIdsInInterval(z0, z1, axis=axis):
                                result[idxBase+k*grid.strides[axis]] = 1

        del ticker
        msg("intersectionGrid finished pass one along axis %d. %d of %d points"
            " dubious. Including...:\n%s"
            % (axis, len(self.dubiousPoints), len(result),
               self.dubiousPoints[:15]),
            debugLevel=5)

        # no dubious points: end here and return result
        if not self.dubiousPoints:
            return result

        # otherwise check and treat dubious points
        axis = aix[0]  # choose new projection axis

        # axis index for the first and second axes perpendicular to the
        # projection direction.
        # (defaults to 0, 1 for the default axis==2,
        # (or else [1,2] or [2,0] for axis==0 or 1)
        aix = [1, 2, 0, 1, 2][axis:axis+2]

        # group dubious points by new-axis-rays
        gridIdxDict = defaultdict(list)
        for ijk in self.dubiousPoints:
            gridIdxDict[(ijk[aix[0]], ijk[aix[1]])].append(ijk[axis])

        # determine points for which we'll search projections
        rayStartGridIds = sorted(gridIdxDict)
        NN = grid.gridPtNb  # abbreviation
        if axis==0:
            gridIdIter = ((0,j,k) for j,k in rayStartGridIds)
        elif axis==1:
            gridIdIter = ((i,0,k) for k,i in rayStartGridIds)
        else:  # axis==2
            gridIdIter = ((i,j,0) for i,j in rayStartGridIds)
        rayStartPoints = [grid.getPoint(ijk) for ijk in gridIdIter]

        self.dubiousPoints = []  # clear dubious points list for the second loop
        axisStrides = grid.strides[axis]  # abbreviation
        msg("intersectionGrid starting second pass along axis %d. checking %d"
            " points. Including...:\n%s"
            % (axis, len(rayStartGridIds), rayStartGridIds[:10]), debugLevel=5)

        # project points in the second direction
        intersections = self.surface.projectPoints(rayStartPoints, axis=axis)

        ticker = MsgTicker("...processing %s-ray %%d/%d"
                           % ("xyz"[axis], len(rayStartPoints)))
        for (i,j), intersPts in izip(rayStartGridIds, intersections):
            ticker.tick()
            idxBase = i*grid.strides[aix[0]]+j*grid.strides[aix[1]]
            interZ = sorted(pt[axis] for pt in intersPts)
            dubPtsK = gridIdxDict[(i,j)]
            dubPtsZ = [grid.origin[axis] + k*grid.spacing[axis]
                       for k in dubPtsK]
            if (len(interZ) % 2) == 0:
                # if number of intersections is even then assign dubious points

                if len(interZ)==0:
                    # no intersection in this ray: all dubious pt are out
                    for k in dubPtsK:
                        result[idxBase+k*axisStrides] = 0
                else:
                    for k, z in izip(dubPtsK, dubPtsZ):
                        # index of the dubious point in sorted intersection list
                        interIdx = bisect.bisect_left(interZ, z)
                        if interIdx>=len(interZ):  # point outside intersections
                            result[idxBase+k*axisStrides] = 0
                        elif (interZ[interIdx]==z  # point exactly on surface =>in
                            or interIdx%2==1):  # point iside of an interval
                            result[idxBase+k*axisStrides] = 1
                        else:
                            result[idxBase+k*axisStrides] = 0
            
            else:
                # if number of intersections is odd then points are dubious...

                # very simple algorithm, programmed rather inefficiently
                # only consider first and last intersection
                # possible improvements:
                # - sort k-indexes of dubious points, find interZ[0] in this
                #   list and interZ[-1] and assign all dubious points in between
                #   to "in" and rest to "out". Still same result, maybe faster.
                # - better results: check neighbours of dubious points if they
                #   are dubious as well and if there is an intersection in
                #   between (or multiple intersections).
                
                for k, z in izip(dubPtsK, dubPtsZ):
                    if interZ[0] <= z <= interZ[-1]:
                        result[idxBase+k*axisStrides] = 1
                        self.dubiousPoints.append(ijkBase)
                    else:
                        result[idxBase+k*axisStrides] = 0

        del ticker
        msg("intersectionGrid finished second pass along axis %d. %d points"
            " still dubious." % (axis, len(self.dubiousPoints)),
            debugLevel=5)
        return result

    ## ------------------------------------------------------------------------
    ## ------------------------------------------------------------------------
    ## transforming the volume

    def translate(self, move_vector):
        """Move the volume in space by move_vector

        Example:

        >>> # move this volume up by 10
        >>> vol.translate(move_vector=[0,0,10])

        Translating is done by simply moving the point coordinates for this
        volume. The surface that was used to create this volume is moved as
        well. (There is actually no difference between moving the initial
        surface and moving this volume.)

        The point search initialization is modified and doesn't have to be done
        again.
        """

        self.surface.translate(move_vector)


    def scale(self, scale_factor, scale_origin):
        """Scale the volume by a certain factor relative to scale_origin

        scale_factor may be a single number in which case it is applied to all
        three dimensions or a list of three numbers with one factor for each
        dimension.

        Example:

        >>> # scale this volume by 2 in each direction
        >>> vol.scale(scale_factor=2, scale_origin=[0,0,0])

        Scaling is done by simply moving the point coordinates for this
        volume. The surface that was used to create this volume is scaled as
        well. (There is actually no difference between scaling the initial
        surface and scaling this volume.)

        After scaling point search has to be initialized again.
        """

        self.surface.scale(scale_factor, scale_origin)


###################################################################
###################################################################
###################################################################
###################################################################


class VolumeBox(Volume):
    """Volume class representing a box aligned to the coordinate axes.

    >>> from bae.volume_02 import VolumeBox
    >>> vol = VolumeBox(box=[xyz_min, xyz_max])
    """

    def __init__(self, box, logfile=None):
        """
        @param box: [[xmin, ymin, zmin], [xmax, ymax, zmax]]
        @param logfile: deprecated argument, only for compatibility reasons.
              No effect.
        """
        Volume.__init__(self)

        # initialize internal variable, use getBoundigBox method to query it
        assert len(box[0]) and len(box[1])
        self.boundingBox = box

        return

    def getBoundingBox(self):
        """Bounding Box of the whole volume.
        Returns a list of two coordinate tuples (rather lists)
        The first coord tuple states the min value, the last the max.
        I.e. self.getBoundingBox()[i][j]:
        i \in {0,1} for min/max value, j \in {0,1,2} coordinate index

        Example:

        >>> # bounding box of that volume
        >>> bb = v.getBoundingBox()
        >>> print "Bounding box: %s" % bb

        @Note: The result is a deep copy of the volume's box itself. It's
        save to modify the contents of the bounding box result. You can change
        the volume itself by modifying its self.boundingBox attribute.
        """
        return BoundingBox(box=(list(xyz) for xyz in self.boundingBox))

    def getCentroid(self):
        """Centroid of the volume.
        Returns the coordinate tuple (rather list) of the centre of this box
        """
        cnr2ctr = vector_scale(vector(*self.boundingBox),0.5)
        return vector_plus(self.xyzOrigin, cnr2ctr)

    def getExteriorSurface(self):
        """return a TriangleSurface object with all outer triangles of this
        volume.
        """

        raise Exception("VolumeBox.getExteriorSurface() not implemented yet.")

    def pointIsInside(self, point):
        """test if point is inside the volume
        point is tuple of coordinates
        """
        inside = True
        for i in xrange(3):
            if (point[i]<self.boundingBox[0][i]
                or point[i]>self.boundingBox[1][i]):
                inside = False
                break
        return inside

    def intersectionGrid(self, grid, fieldName="insideVolume"):
        """Find points of the regular grid inside self.

        Example: Get a field indicating points in or outside the volume:
         >>> from bae.mesh_01 import MeshStructuredPoints
         >>> from bae.volume_02 import VolumeBox
         >>> from bae.vtk_02 import FieldsOnStructPointGrid as Vtk
         >>>
         >>> vol = VolumeBox(...)
         >>> grid = MeshStructuredPoints(....)
         >>> insideField = vol.intersectionGrid(grid)
         >>>
         >>> vtk = Vtk(mesh=grid)
         >>> vtk.updateFields(insideField)
         >>> vtk.toVtk("inVolume.vtk")

        @param grid: L{bae.mesh_01.MeshStructuredPoints}-object or a
            L{MeshStructuredPointsRot} object.

        @param fieldName: fieldName attribute of the field to be returned

        @returns: L{bae.field_01.Field}-object (position="structPt",
            dataType="scalar") containing 1 for "in", 0 for "out"
        """
        if hasattr(grid, "rotmat"):
            # is a MeshStructuredPointsRot...
            # Not Implemented
            return super(self).intersectionGrid(grid, fieldName=fieldName)
        else:
            # transform lower and upper corner of self into grid-coordinates
            x0 = [float(x)/r for x, r in izip(
                vector(grid.origin, self.boundingBox[0]), grid.spacing)]
            x1 = [float(x)/r for x, r in izip(
                vector(grid.origin, self.boundingBox[1]), grid.spacing)]

            # initialize resulting field with 0
            result = Field.classFromPosType(fieldName, "structPt", "scalar")(
                [0]*len(grid))

            # set points inside bounds to 1
            for k in range(int(ceil(x0[2])), 1+int(floor(x1[2]))):
                for j in range(int(ceil(x0[1])), 1+int(floor(x1[1]))):
                    i0 = int(ceil(x0[0]))
                    i1 = 1+int(floor(x1[0]))
                    offs = j*grid.strides[1] + k*grid.strides[2]
                    result[offs+i0:offs+i1] = [1]*(i1-i0)
            
        return result
    
    def scale(self, scale_factor, scale_origin=None):
        """Scale the volume relative to scale_origin by scale_factor.

        Example:
          >>> # scale this volume by 2 in each direction (volume := volume*8)
          >>> v.scale(scale_factor=2, scale_origin=[0,0,0])

        @param scale_factor: May be a single number in which case it is applied
        to all three dimensions or a list of three numbers with one factor for
        each dimension.
        @param scale_origin: The scaling in done relative to this point.
        Defaults to the centre of the volume.
        """
        if scale_origin is None:
            scale_origin = vector_scale(
                vector_plus(*self.boundingBox), 0.5)

        if type(scale_factor) in (float, int):
            scale_factor = [scale_factor,scale_factor,scale_factor]

        for i in xrange(2):
            orig2old = vector(scale_origin, self.boundingBox[i])
            orig2new = [orig2old[0]*scale_factor[0],
                        orig2old[1]*scale_factor[1],
                        orig2old[2]*scale_factor[2]]
            self.boundingBox[i] = vector_plus(scale_origin, orig2new)

    def translate(self, move_vector):
        """Translate the volume by move_vector.

        Example:

        >>> # move this volume up by 10
        >>> vol.translate(move_vector=[0,0,10])
        """
        for i in xrange(2):
            self.boundingBox[i] = vector_plus(self.boundingBox[i], move_vector)


class VolumeBoxRotZ(Volume):
    """A box volume with sides originally aligned to the coordinate axes and
    then rotated about the vertical axis (parallel to z axis) through the first
    point of the box argument.

    (The rotation is about the lower left edge of the box when seen from above).

    >>> from bae.volume_02 import VolumeBoxRotZ
    >>> vol = VolumeBoxRotZ(box=[xyz_min, xyz_max], rotZ=60.0)
    """

    def __init__(self, box, rotZ=0.0, logfile=None):
        """
        @param box: coordinates before rotation
        [[xmin, ymin, zmin], [xmax, ymax, zmax]]

        @param rotZ: angle in degree to rotate about the [xmin, ymin]-z-axis

        @param logfile: deprecated argument, only for compatibility reasons.
              No effect.
        """
        Volume.__init__(self)

        self.xyzOrigin = list(box[0])
        self.span = vector(box[0], box[1])

        alpha = rotZ*pi/180
        self.ca=cos(alpha)
        self.sa=sin(alpha)
        # rotated vec: [ca*v[0]+sa*v[1], -sa*v[0]+ca*v[1], v[2]]

        # due to the transformation routines involved a tolerance seems
        # necessary
        self.tolerance = 1E-4

        return

    def getCentroid(self):
        """Centroid of the volume.
        Returns the coordinate tuple (rather list) of the centre of this box
        """
        cnr2ctr = vector_scale(self.span,0.5)
        cnr2ctr = [self.ca*cnr2ctr[0]-self.sa*cnr2ctr[1],
                   +self.sa*cnr2ctr[0]+self.ca*cnr2ctr[1],
                   cnr2ctr[2]]
        return vector_plus(self.xyzOrigin, cnr2ctr)

    def getBoundingBox(self):
        """Bounding Box of the whole volume.
        Returns a list of two coordinate tuples (rather lists)
        The first coord tuple states the min value, the last the max.
        I.e. self.getBoundingBox()[i][j]:
        i \in {0,1} for min/max value, j \in {0,1,2} coordinate index

        Example:

        >>> # bounding box of that volume
        >>> bb = v.getBoundingBox()
        >>> print "lower left front corner:", bb[0]
        >>> print "upper right back corner:", bb[1]

        @Note: Not implemenmted yet. Raises an exception.
        """
        raise Exception("VolumeBoxRotZ.getBoundingBox() not implemented yet.")

    def getExteriorSurface(self):
        """return a TriangleSurface object with all outer triangles of this
        volume.
        """
        raise Exception("VolumeBox.getExteriorSurface() not implemented yet.")

    def pointIsInside(self, point):
        """test if point is inside the volume box
        point is tuple of coordinates
        """
        v = [point[0]-self.xyzOrigin[0],
             point[1]-self.xyzOrigin[1],
             point[2]-self.xyzOrigin[2]]
        v = [self.ca*v[0]+self.sa*v[1],
             -self.sa*v[0]+self.ca*v[1],
             v[2]]
        inside = all(
            x>=-self.tolerance and x<=y+self.tolerance
            for x, y in izip(v, self.span))
        return inside

    def scale(self, scale_factor, scale_origin):
        """Scale the volume relative to scale_origin by scale_factor.

        Example:
          >>> # scale this volume by 2 in each direction (volume := volume*8)
          >>> v.scale(scale_factor=2, scale_origin=[0,0,0])
        """
        self.span = vector_scale(self.span, scale_factor)
        # new xyzOrigin
        old2orig = vector(self.xyzOrigin, scale_origin)
        old2new = vector_scale(old2orig, 1.0-scale_factor)
        vector_modif_add(self.xyzOrigin, old2new)

    def translate(self, move_vector):
        """Translate the volume by move_vector.

        Example:

        >>> # move this volume up by 10
        >>> vol.translate(move_vector=[0,0,10])
        """
        self.xyzOrigin = vector_plus(self.xyzOrigin, move_vector)


###################################################################
###################################################################
###################################################################
###################################################################


class VolumeSphere(Volume):
    """Volume class representing a sphere.

    >>> from bae.volume_02 import VolumeSphere
    >>> vol = VolumeSphere(centre=[101.0,734.0,1137.0], radius=5.0)
    """

    def __init__(self, centre, radius, logfile=None):
        """
        @param centre: [x,y,z]
        @param radius:
        @param logfile: deprecated argument, only for compatibility reasons.
              No effect.
        """
        Volume.__init__(self)
        self.centre = centre
        self.radius = radius
        return

    def getBoundingBox(self):
        """Bounding Box of the whole volume.
        Returns a list of two coordinate tuples (rather lists)
        The first coord tuple states the min value, the last the max.
        I.e. self.getBoundingBox()[i][j]:
        i \in {0,1} for min/max value, j \in {0,1,2} coordinate index

        Example:

        >>> # bounding box of that volume
        >>> bb = v.getBoundingBox()
        >>> print "lower left front corner:", bb[0]
        >>> print "upper right back corner:", bb[1]
        """
        rvec = [self.radius, self.radius, self.radius]
        return BoundingBox(box=[vector_plus(self.centre, rvec),
                                vector_minus(self.centre, rvec)])

    def getExteriorSurface(self):
        """return a TriangleSurface object with all outer triangles of this
        volume.
        """
        raise Exception(
            "VolumeSphere.getExteriorSurface() not implemented yet.")

    def pointIsInside(self, point):
        """test if point is inside the volume
        point is tuple of coordinates
        """
        d = dist(self.centre, point)
        return d <= self.radius

    def scale(self, scale_factor, scale_origin):
        """Scale the volume relative to scale_origin by scale_factor.

        Example:

        >>> # scale this volume by 2 in each direction (volume := volume*8)
        >>> v.scale(scale_factor=2, scale_origin=[0,0,0])
        """
        self.radius *= scale_factor
        orig2old = vector(scale_origin, self.centre)
        orig2new = vector_scale(orig2old, scale_factor)
        self.centre = vector_plus(scale_origin, orig2new)

    def translate(self, move_vector):
        """Translate the volume by move_vector.

        Example:

        >>> # move this volume up by 10
        >>> vol.translate(move_vector=[0,0,10])
        """
        self.centre = vector_plus(self.centre, move_vector)


###################################################################
###################################################################
###################################################################
###################################################################

class VolumeVoxel(Volume):
    """
    Volume represented by a regular grid of cells of rectangular cuboid shape.

    @ivar origin: lower corner of first cell (not its centroid)
    @ivar cellSize: a float-triple giving the size of a single cell (voxel)
        in x,y,z direction respectively.
    @ivar cellNb: an integer-triple, number of cells in each direction
    @ivar strides: A triple of index offsets to the next cell in x,y,z
      direction. I.e. In the list self.inside strides[0] is the offset to the
      next neighbouring cell in x direction, strides[1] is next in y direction
      and strides[2] in z direction.
    @ivar inside: A list of bools. True means the cell belongs to the volume.

    @note: This class is somehow similar and related to
    L{bae.mesh_01.MeshStructuredPoints}. But there are also significant
    differences you should be aware of. For example the instance attribute
    origin has a slightly different meaning in both classes.
    """
    def __init__(self, *args, **kwargs):
        """Initialize a dummy empty VolumeVoxel object. If given another
        VolumeVoxel argument then serves as copy constructor.

        @param args: If given a single positional argument of type VolumeVoxel
           then initialize a duplocate of this given object.
        @kwarg origin: If keyword arguments origin, cellSize and cellNb are
           given then initialize a VolumeVoxel-object with corresponding
           attributes. See the corresponding instance attributes for a
           description.
        @kwarg cellSize: see keyword argument origin
        @kwarg cellNb: see keyword argument origin
        @kwarg inside: If keyword arguments origin, cellSize and cellNb are
           given then inside may optionally be given to initialise the
           corresponding instance attribute as well.
        """
        if not args and not kwargs:
            # initialize empty dummy values
            self.origin = [0.0, 0.0, 0.0]
            self.cellSize = [1.0, 1.0, 1.0]
            self.cellNb = [0, 0, 0]
            self.strides = [1, 1, 1]
            self.inside = []
        elif len(args)==1 and isinstance(args[0], VolumeVoxel):
            # copy constructor
            other = args[0]
            self.origin = list(other.origin)
            self.cellSize = list(other.cellSize)
            self.cellNb = list(other.cellNb)
            self.strides = list(other.strides)
            self.inside = list(other.inside)
        elif all(x in kwargs for x in ("origin", "cellSize", "cellNb")):
            self.origin = list(kwargs["origin"])
            self.cellSize = list(kwargs["cellSize"])
            self.cellNb = list(kwargs["cellNb"])
            self.strides = [1, self.cellNb[0], self.cellNb[0]*self.cellNb[1]]
            N = self.cellNb[0]*self.cellNb[1]*self.cellNb[2]
            try:
                self.inside = list(kwargs["inside"])
                if len(self.inside) != N:
                    raise ValueError(
                        "The inside argument has %d items but there are %d"
                        " cells to supply values for."
                        % (len(self.inside), N))
            except KeyError:
                self.inside = [False]*N
        else:
            raise NotImplementedError(
                "Don't know what kind of VolumeVoxel to create from these"
                " parameters: args: %s; kwargs: %s"
                % (args, kwargs))

    @classmethod
    def fromFieldOnPtGrid(cls, grid, field):
        """Initialize VolumeVoxel object based on a field defined on a regular
        points grid. Those points are being considered the cell/voxel centroids.

        Usage:
         >>> from bae.volume_02 import VolumeVoxel
         >>> vol = VolumeVoxel.fromVtkFieldOnCellCentr("XFER_10.vtk")
         >>> elset = vol.intersectionElset(mesh)

        Or:
         >>> from bae.vtk_02 import FieldsOnStructPointGrid as Vtk
         >>> from bae.volume_02 import VolumeVoxel
         >>> vtk = Vtk().fromVtk(vtkFileName)
         >>>
         >>> # apply modifyData() function to data in vtk
         >>> firstFieldName = vtk.data.iterkeys().next()
         >>> firstField = vtk.data[firstFieldName]
         >>> vtk.data[firstFieldName] = type(firstField)(
         >>>    modifyData(x) for x in firstField)
         >>>
         >>> # get elements inside that volume
         >>> vol = VolumeVoxel.fromVtkFieldOnCellCentr(vtk)
         >>> elset = vol.intersectionElset(mesh)
        """
        self = cls()

        # origin is corner of first cell; grid.origin is centroid of first cell
        self.origin = vector(vector_scale(grid.spacing, 0.5), grid.origin)
        self.cellSize = grid.spacing
        self.cellNb = grid.gridPtNb
        self.strides = grid.strides
        self.inside = [bool(x) for x in field]

        nbCellsInside = sum(int(x) for x in self.inside)
        vol = nbCellsInside*self.cellSize[0]*self.cellSize[1]*self.cellSize[2]
        msg("Initialized volume. Overall volume size: %g" % vol)

        return self

    @classmethod
    def fromVtkFieldOnCellCentr(
            cls, vtk, insideFilterForFirstField=None):
        """Initialize VolumeVoxel object based on fields defined on a regular
        points grid. Those points are being considered the cell/voxel centroids.

        Usage:
         >>> from bae.volume_02 import VolumeVoxel
         >>> vol = VolumeVoxel.fromVtkFieldOnCellCentr("XFER_10.vtk")
         >>> elset = vol.intersectionElset(mesh)

        Or:
         >>> from bae.vtk_02 import FieldsOnStructPointGrid as Vtk
         >>> from bae.volume_02 import VolumeVoxel
         >>> vtk = Vtk().fromVtk(vtkFileName)
         >>>
         >>> # apply modifyData() function to data in vtk
         >>> firstFieldName = vtk.data.iterkeys().next()
         >>> firstField = vtk.data[firstFieldName]
         >>> vtk.data[firstFieldName] = type(firstField)(
         >>>    modifyData(x) for x in firstField)
         >>>
         >>> # get elements inside that volume
         >>> vol = VolumeVoxel.fromVtkFieldOnCellCentr(vtk)
         >>> elset = vol.intersectionElset(mesh)

        @param vtk: Can be the name of a vtk file (a string). In this case the
        file will be read. Or it can be a L{bae.vtk_02.FieldsOnStructPointGrid}
        object.
        @param insideFilterForFirstField: a function that will be called for
        each single grid point value. It gets the value of the first field in
        the given vtk-object and is expected to return True for "inside" and
        False for "outside".

        If not given (or None) then defaults to "bool", i.e. all non zero
        values indicate "inside".

        @Note: The plan is to add another argument insideFilter that get's the
        vtk.data attribute as a whole and is expected to return a list of bools
        to be used to initialize self.inside. As a more general procedure...
        """

        # default arguments
        if insideFilterForFirstField is None:
            def insideFilterForFirstField(x):
                return bool(x)

        self = cls()

        # if filename specified then load that file
        if isinstance(vtk, basestring):
            vtkFileName = vtk
            vtk = Vtk().fromVtk(vtkFileName)
        else:
            vtkFileName = ""  # for diagnostic output

        # abbreviation
        grid = vtk.mesh

        # origin is corner of first cell; grid.origin is centroid of first cell
        self.origin = vector(vector_scale(grid.spacing, 0.5), grid.origin)
        self.cellSize = grid.spacing
        self.cellNb = grid.gridPtNb
        self.strides = grid.strides

        if insideFilterForFirstField:
            firstFldName = vtk.data.iterkeys().next()
            msgText = (
                "from field %s in vtk file %s" % (firstFldName, vtkFileName))
            fld = vtk.data[firstFldName]
            self.inside = [insideFilterForFirstField(x) for x in fld]
        else:
            self.inside = [False,]*len(grid)

        nbCellsInside = sum(int(x) for x in self.inside)
        vol = nbCellsInside*self.cellSize[0]*self.cellSize[1]*self.cellSize[2]
        msg("Initialized volume %s. Overall volume size: %g" % (msgText, vol))

        return self

    def getBoundingBox(self):
        """Get the bounding box of the volume.

        @Note: The current implementation is not correct and has been switched
        off!
        """
        raise NotImplementedError(
            "VolumeVoxel.getBoundingBox() not implemented yet.")
        # only consider voxels with self.inside == true! Not easy!
        return BoundingBox([
            list(self.origin),  # make a copy
            vector_plus(self.origin,
                        [n*d for n, d in izip(self.cellNb, self.cellSize)])
            ])

    def getCellCentroidsGrid(self):
        """Get a L{bae.mesh_01.MeshStructuredPoints}-object of the cell
        centroids of the voxel-grid, no matter whether they are in the volume
        or not.

        This is merely for internal purposes.
        """
        return MeshStructuredPoints(
            origin=[(x+0.5*d) for x, d in izip(self.origin, self.cellSize)],
            spacing=self.cellSize, gridPtNb=self.cellNb)

    def getCellIdx(self, point):
        """If the given point is in any of the cells of self then return its
        indexes (i,j,k integer triple). This function is merely for internal
        purposes.

        Usage:
         >>> idx = vol.getCellIdx(point)
         >>> inside = idx and vol.inside[
         >>>     sum(i*n for i, n in izip(idx, vol.strides))]
        """
        point0 = [(x-d) for x, d in izip(point, self.origin)]
        # treated first because following formula is wrong on negative values
        # ... int(-0.5)==0
        if any(x<0 for x in point0):
            return False

        idx = [int(x/d) for x, d in izip(point0, self.cellSize)]
        if any(i>=n for i, n in izip(idx, self.cellNb)):
            return False
        else:
            return idx

    def pointIsInside(self, point):
        """Determine whether the given point is in the volume (self).
        """
        ptIdx = self.getCellIdx(point)
        return ptIdx and self.inside[
            sum(i*n for i, n in izip(ptIdx, self.strides))]

    def getConnectedVolumePart(self, point, connectionType="face"):
        """Find the volume part of self that is connect to the given point

        @param point: list of x,y,z values
        @param connectionType: the kind of connection to be recognized: "face",
           "edge" or "corner". (Obviously "corner" implies "edge" and "face"
           and "edge" implies "face".)
        @returns: a L{VolumeVoxel}-object
        """
        if connectionType=="face":
            neighbourIdxOffs = [
                (-1,0,0), (1,0,0),
                (0,-1,0), (0,1,0),
                (0,0,-1), (0,0,1)]
        elif connectionType=="edge":
            neighbourIdxOffs = [
                (-1,0,0), (1,0,0),
                (0,-1,0), (0,1,0),
                (0,0,-1), (0,0,1),

                (-1,-1,0), (1,-1,0),
                (-1,1,0),  (1,1,0),
                (-1,0,-1), (1,0,-1),
                (-1,0,1),  (1,0,1),
                (0,-1,-1), (0,1,-1),
                (0,-1,1),  (0,1,1),
                (-1,-1,0), (-1,1,0),
                (1,-1,0),  (1,1,0),
                (-1,0,-1), (-1,0,1),
                (1,0,-1),  (1,0,1),
                (0,-1,-1), (-1,0,1),
                (0,1,-1),  (1,0,1),
            ]
        if connectionType=="corner":
            raise NotImplementedError(
                "VolumeVoxel.getConnectedVolumePart(connectionType='corner')"
                " not implemented yet. But is very easy to achieve. See source"
                " code.")

        # identify and check first point
        ptIdx = self.getCellIdx(point)
        if ptIdx:
            ptII = sum(i*n for i, n in izip(ptIdx, self.strides))
        if not(ptIdx and self.inside[ptII]):
            # point is not in the volume!
            return VolumeVoxel()

        inside = [False]*len(self.inside)
        inside[ptII] = True
        border = [ptIdx,]
        while border:
            ptIdx0 = border.pop()
            for idxOffs in neighbourIdxOffs:
                ptIdx1 = [i+j for i, j in izip(ptIdx0, idxOffs)]
                if not(all(0 <= i < n for i, n in izip(ptIdx1, self.cellNb))):
                    continue
                ptII1 = sum(i*n for i, n in izip(ptIdx1, self.strides))
                if not self.inside[ptII1] or inside[ptII1]:
                    continue
                # this neighbour cell is "in" and new (and of course connected)
                inside[ptII1] = True
                border.append(ptIdx1)

        # end, all neighbours exhausted...
        result = VolumeVoxel(self)
        result.inside = inside
        return result

    def getVolume(self):
        """Calculate the geometrical volume. (number of "inside"-cells times
        cell volume)
        """
        return (self.cellSize[0]*self.cellSize[1]*self.cellSize[2])*sum(
            int(i) for i in self.inside)

    def getCentroid(self):
        """Centroid of the volume.
        Returns the coordinate tuple (rather list) of the centre of gravity of
        this volume.
        """
        raise NotImplementedError("VolumeVoxel.getCentroid() not implemented yet."
                             " But it's very simple. See comments in code.")
        # 1. sum up individually:
        #   [idx for flag, idx in izip(
        #        self.inside, ((i,j,k) for k in range(self.cellNb[2])
        #                              for j in range(self.cellNb[1])
        #                              for i in range(self.cellNb[0]))
        #    if flag]
        # 2. scale by cellSize/len(self.inside)
        # 3. translate by origin+0.5*cellsize

    def getExteriorSurface(self):
        """Returns a TriangleSurface object with all outer triangles of this
        volume.
        """

        raise NotImplementedError(
            "VolumeVoxel.getExteriorSurface() not implemented yet.")
