# -*- coding: utf-8 -*-
"""Fundamental base classes for other unstructured mesh topo classes.
"""

import numpy as np

from bae.log_01 import msg
from bae.abq_model_02 import Model as AbqModel, checkModelModuleVersion
from bae.mesh_01 import Mesh as PurePyMesh

from ..base import TopoBase
from .elementType import ElementType

class HomoMesh(object):
    """Base class for homogeneous *unstructured* meshes that have only one type
    of elements. This serves as base class for subclasses for certain types
    like L{TetMesh}, L{TriMesh}. But will also be used directly.

    @ivar elType: element type information of type L{ElementType}
    @ivar nodeCoords: ... node coords, ndarray Nx3
    @ivar elNodes: ... connectivity, ndarray M x elType.nodesPerElem
    """

    # output format for coordinates
    coordFormat = "%.1f"

    def __init__(self, *args, **kwargs):
        """Initialize a new HomoMesh object.

        The first (positional) argument might be another L{HomoMesh}-object.
        In this case a *deep* copy is performed.

        @kwarg nodeCoords: optional, initialize self.nodeCoords with this
        @kwarg elNodes: optional, initialize self.elNodes with this
        @kwarg elType: optional, initialize self.elType with this
        """

        if args and isinstance(args[0], HomoMesh):
            # Copy-constructor
            other = args[0]
            self.elNodes = deepcopy(other.elNodes)
            self.nodeCoords = deepcopy(other.nodeCoords)
            self.elType = other.elType
        elif args:
            #... need to differentiate:
            # ... and isinstance(args[0], Whatever-Surface):
            raise ValueError(
                "HomoMesh copy constructor called with object of other type"
                " than HomoMesh. args: %s" % args)
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
                self.elType = ElementType()

                # provide minimal information in case we have some elements
                if self.elNodes:
                    self.elType.nodesPerElem = self.elNodes.shape[1]

    ######################################################
    #{ inititialize new mesh object
    @classmethod
    def fromAbqModel(cls, abqModel, elset=None):
        """Create a HomoMesh instance from some elements of a
        L{bae.mesh_01.Mesh}- or an L{abaqus model<bae.abq_model_02.Model>}-
        object.

        Elements in elset that are not found in abqModel.elNodes are ignored, a
        warning is posted.

        @param abqModel: a L{bae.mesh_01.Mesh}- or an
        L{abaqus model<bae.abq_model_02.Model>}-object or a filename (a string)
        or an open file of an abaqus input file.

        @param elset: only consider elements of this elset.
        If abqModel is an AbqModel then elset might be anything that
        abqModel.getUnionSet() accepts as input: an elset name, a set of
        element numbers, a list of elset names and element numbers. Otherwise
        it must be an iterable of element numbers.

        @note: This creates additional attributes for the HomoMesh object:
        elemLabel and nodeLabel.
        """

        # process the abqModel argument: determine the source mesh
        if isinstance(abqModel, (basestring, file)):
            # basestring = any kind of string
            m = AbqModel()
            m.read(abqModel)
            abqModel = m
        elif isinstance(abqModel, AbqModel):
            checkModelModuleVersion(abqModel, "%s.%s.fromAbqModel"
                                    % (__name__, cls.__name__))
        elif not isinstance(abqModel, PurePyMesh):
            ValueError(
                "%s.%s.fromAbqModel expects some sort of Mesh as first"
                " argument, istead got a %s."
                % (__name__, cls.__name__, type(abqModel)))

        # process elset argument: take which elements from the mesh?
        if elset is None:
            elset = abqModel.elNodes
        else:
            try:
                # see if we have a getUnionSet-method (if it's an AbqModel)
                elset = abqModel.getUnionSet("elset", elset)
            except AttributeError:
                # if not, elset must be a set of element numbers, make a copy
                elset = set(elset)
                pass
            elemsNotFound = elset.difference(abqModel.elNodes)
            if elemsNotFound:
                msg("WARNING from %s.%s.fromAbqModel():"
                    " Could not find element connectivity for %d of the %d"
                    " elements in the specified elset."
                    % (__name__, cls.__name__, len(elemsNotFound), len(elset)))
                elset.difference_update(elemsNotFound)

        # find connected nodes and check that all nodes are defined
        connectedNodes = abqModel.getConnectedNodes(elset)
        if connectedNodes.difference(abqModel.nodeCoords):
            elemsMissingNodes = [
                e for e in elset
                if any(n not in abqModel.nodeCoords
                       for nodes in abqModel.elNodes[e]
                       for n in nodes)]
            msg("WARNING from %s.%s.fromAbqModel():"
                " For %d of %d elements some nodes where not defined in"
                " the model. Ignoring those elements."
                % (__name__, cls.__name__, len(elemsMissingNodes), len(elset)))
            elset.difference_update(elemsMissingNodes)

        # determine elType
        try:
            elType = abqModel.elType
        except AttributeError:
            # abqModel has no elType, i.e. it is an mesh_01.Mesh
            shapes = set(abqModel.elShape[e] for e in elset)
            if len(shapes)==1:
                elShape = shapes.pop()
                elType = AbqModel.defaultMapShapeToType[elShape]
                elType = ElementType(abqName=elType)
            else:
                raise ValueError(
                    "HomoMesh.fromAbqModel() can only initialize elements of one"
                    " element shape. Instead found: %s" % elShape)
            
        else:
            # abqModel has elType, i.e. it is an abq_model_02.Model
            types = set(elType[e] for e in elset)
            if len(types)==1:
                elType = ElementType(abqName=types.pop())
            else:
                raise ValueError(
                    "HomoMesh.fromAbqModel() can only initialize elements of one"
                    " type. Instead found: %s" % elType)

        # initialize data members
        elset = sorted(elset)
        elemLabel = np.array(elset, dtype=int)
        connectedNodes = sorted(connectedNodes)
        nodeLabel = np.array(connectedNodes, dtype=int)

        nodeCoords = np.array(
            [abqModel.nodeCoords[n][:elType.nodesPerElem]
             for n in connectedNodes], dtype=float)
        nodeLabelToIdx = dict((n, i) for i, n in enumerate(connectedNodes))
        elNodes = np.array([
            [nodeLabelToIdx[n]
             for n in abqModel.elNodes[e][:elType.nodesPerElem]]
            for e in elset], dtype=int)

        if len(elNodes)==0:
            msg('WARNING: %s.%s.fromAbqModel initialized without elements.'
                % (__name__, cls.__name__))

        res = cls(nodeCoords=nodeCoords, elNodes=elNodes, elType=elType)
        res.elemLabel = elemLabel
        res.nodeLabel = nodeLabel
        res.elType = elType

        return res

    #{ end of inititialize new mesh object

    ######################################################
    #{ methods for data extraction / getting info
    def getBoundingBox(self):
        """
        Returns the bounding box of the surface

        Example:
          >>> # bounding box of that surface
          >>> bb = surf.getBoundingBox()
          >>> print "lower left front corner:", bb[0]
          >>> print "upper right back corner:", bb[1]

        @Returns: An array of shape (2,3) of the min and max coordinates.
        I.e. self.getBoundingBox()[i,j]:
        i \in {0,1} for min/max value, j \in {0,1,2} coordinate index
        """
        a = self.nodeCoords
        boundingBox = np.array((a.min(axis=0), a.max(axis=0)))
        return boundingBox

    def getElCentroids(self):
        """
        Calculate element centroids for all elements
        """

        msg("WARNING: The general HomoMesh.getElCentroids() method calulates"
            " element centroids by just averaging the node coordinates per"
            " element. This is does not yield the correct centroid for all"
            " types of elements.")

        # array of all triangle-vertex-coords: nodeCoords[i,j,k]=...
        # i-th triangle, j-th vertex (0..2) of the triangle, k-th component
        nodeCoords = self.nodeCoords[self.elNodes]

        # centroid for each triangle
        return np.average(nodeCoords, axis=1)

    def getConnectedNodes(self):
        """Get all nodes that are connected to any of the elements.

        The inverse --all nodes not connected to any element-- might be of more
        interest and can be derived from this:
         >>> notConnected = np.setdiff1d(
         >>>     np.arange(len(mesh.nodeCoords)), mesh.getConnectedNodes(),
         >>>     assume_unique=True)
        """
        return np.unique(self.elNodes)

    def getElemFaces(self):
        """Returns the element faces as NxMxP array of corner node indexes, i.e.
        the values in the resulting array are suitable as indexes in
        nodeCoords. N is the number of elements, len(elNodes). M is the number
        of faces per element, e.g. 4 for tet elements. P is the number of nodes
        identifying the face, i.e. three for tet elements --linear or quadratic.
        """
        return self.elNodes.take(self.elType.faceNodes, axis=1)

    #} end of methods for data extraction / getting info

    ######################################################
    #{ transform geometry
    def scale(self, scale_factor, scale_origin):
        """Scale the mesh by a certain factor relative to scale_origin

        scale_factor must be a single number that is applied to all
        three dimensions.

        Example:
         >>> # scale this mesh by 2 in each direction (volume := volume*8)
         >>> mesh.scale(scale_factor=2, scale_origin=[0,0,0])

        Scaling is done by simply moving the coordinates of all nodes.
        """
        self.nodeCoords = (
            scale_factor*(self.nodeCoords-scale_origin) + scale_origin)
    #} end of transform geometry group


class HeteroMesh(object):
    """Base class for unstructured meshes with mixed element types.

    @ivar nodeCoords: ... all nodes of the mesh
    @ivar elTypeNodes: ... list of (elType, elNodes)-tuples. One tuple for each
        element type.
    @ivar elNodes: ... possibly a HeteroElNodesProxy object
        *if we really need that*
    """

    ######################################################
    #{ methods for data extraction / getting info
    def getBoundingBox(self):
        """
        Returns the bounding box of the surface

        Example:
          >>> # bounding box of that surface
          >>> bb = surf.getBoundingBox()
          >>> print "lower left front corner:", bb[0]
          >>> print "upper right back corner:", bb[1]

        @Returns: An array of shape (2,3) of the min and max coordinates.
        I.e. self.getBoundingBox()[i,j]:
        i \in {0,1} for min/max value, j \in {0,1,2} coordinate index
        """
        a = self.nodeCoords
        boundingBox = np.array((a.min(axis=0), a.max(axis=0)))
        return boundingBox
    #} end of methods for data extraction / getting info
    
    ######################################################
    #{ transform geometry
    def scale(self, scale_factor, scale_origin):
        """Scale the mesh by a certain factor relative to scale_origin

        scale_factor must be a single number that is applied to all
        three dimensions.

        Example:
         >>> # scale this mesh by 2 in each direction (volume := volume*8)
         >>> mesh.scale(scale_factor=2, scale_origin=[0,0,0])

        Scaling is done by simply moving the coordinates of all nodes.
        """
        self.nodeCoords = (
            scale_factor*(self.nodeCoords-scale_origin) + scale_origin)
    #} end of transform geometry group


class HeteroElNodesProxy(object):
    """Service class for L{HeteroMesh}-objects which have an elNodes attribute
    of this type.

    As the connectivity data in HeteroMesh-objects is spread over different
    arrays --one for each element type-- the __getitem__ method of the elNodes
    attribute (i.e. of this class) retrieves the node indexes of the given
    element(s) from the right item in the HeteroMesh.elTypeNodes list.

    ... stores indices if suitable to accelerate the access...
    __getitem__ method somehow does the magic...
    """
    pass


class MeshTopoBase(TopoBase):
    """Base class for L{MeshNodes}, L{MeshElems}, L{MeshGaussPts}. Delivers
    functions common to all those classes.

    These three classes take an unstructured mesh object as single argument
    on construction:
     >>> topo = MeshNodes(TriMesh())

    @cvar position: can be "node", "element", "elemIP" depending on the
    position of the datapoints in the mesh.
    """
    def __init__(self, mesh):
        self.mesh = mesh

    ######################################################
    #{ methods for data extraction / getting info
    def getBoundingBox(self):
        """
        Returns the bounding box of the surfacetopography, i.e. of the mesh

        Example:
          >>> bb = topo.getBoundingBox()
          >>> print "lower left front corner:", bb[0]
          >>> print "upper right back corner:", bb[1]

        @Returns: An array of shape (2,3) of the min and max coordinates.
        I.e. self.getBoundingBox()[i,j]:
        i \in {0,1} for min/max value, j \in {0,1,2} coordinate index
        """
        return self.mesh.getBoundingBox()
    #} end of methods for data extraction / getting info

class MeshNodes(MeshTopoBase):
    """Objects of this type serve as FieldsCollection.topo for nodal fields.

    @ivar mesh: is a HomoMesh or a HeteroMesh object. The same object can be
        used (referenced) by other MeshElems and MeshGaussPts objects.
    @cvar position: "node". For L{FieldsCollection}.topo objects the
      position attribute identifies the position of values ("node", "element",
      "elemIP", "point" or "structPt").
    """
    position = "node"

    def getPoints(self):
        return self.mesh.nodeCoords
        
    def __len__(self):
        return len(self.mesh.nodeCoords)

    def getSubTopo(self, index):
        """Service method for FieldsCollection.__getitem__() selecting from the
        space dimension.

        Returns a new L{MeshNodes} object with a new mesh object of only the
        specified nodes and corresponding elements.

        @param index: A serial integer index --anything suitable to have numpy
           choose from one dimension.
        """
        implementation_notes = """
  If the getSubTopo() method is expected to accept arbitrary indexes then
  MeshNodes needs some kind of filter array: which of the nodes are actually
  part of the topology. Because then it might be that we have elements (which
  can only be part of the mesh with all their nodes) with only some of their
  nodes being selected as part of the sub-topo. Same problem for MeshGaussPts
  but there we might be able to connect this to the problem of grouping
  (assigning) Gauss pts to elements.

  Or we throw an error if getSubTopo gets indexes that would tear apart
  elements.
        """
        raise NotImplementedError(
            "getSubTopo not yet implemented for class MeshNodes")


class MeshElems(MeshTopoBase):
    """Objects of this type serve as FieldsCollection.topo for whole element
    fields.

    @ivar mesh: is a HomoMesh or a HeteroMesh object. The same object can be
        used (referenced) by other MeshElems and MeshGaussPts objects.
    @cvar position: "element". For L{FieldsCollection}.topo objects the
      position attribute identifies the position of values ("node", "element",
      "elemIP", "point" or "structPt").
    """
    position = "element"

    def getPoints(self):
        """Returns the cell/element centroid coordinates
        """
        return self.mesh.getElCentroids()

    def __len__(self):
        return len(self.mesh.elNodes)

    def getSubTopo(self, index):
        """Service method for FieldsCollection.__getitem__() selecting from the
        space dimension.

        Returns a new L{MeshElems} object with a new mesh object of only the
        specified elements.

        @param index: A serial integer index --anything suitable to have numpy
           choose from one dimension.
        """
        if not isinstance(self.mesh, HomoMesh):
            raise NotImplementedError(
                "getSubTopo not yet implemented for meshes of type %s"
                % type(self.mesh))

        # all nodes connected to the given element set (index)
        filteredElNodes = self.mesh.elNodes[index,:]
        nodes = np.unique(filteredElNodes)

        # initialize new mesh object
        # ... later to be fed to the new MeshElems object
        subMesh = type(self.mesh)()
        subMesh.elType = self.mesh.elType
        subMesh.nodeCoords = self.mesh.nodeCoords[nodes,:]

        # create lookup array to look up a new node index for each old node id
        # ...an array containing the new node index at the position
        # corresponding to the old nodeCoords array
        # initialize with -1 values for nodes not needed anymore.
        oldNodesToNew = np.full(len(self.mesh.elNodes), -1)
        oldNodesToNew[nodes] = np.arange(len(nodes))

        subMesh.elNodes = oldNodesToNew[filteredElNodes]

        return MeshElems(subMesh)


class MeshGaussPts(MeshTopoBase):
    """Objects of this type serve as FieldsCollection.topo for integration
    point fields.

    MeshGaussPts must somehow deliver a means to identify / group the
    values/indexes for each element. So that we can bring together whole element
    and GaussPt values. For example supply an array "elGPIds" Mx2, dtype int,
    first index is element index in mesh "i_m", then elGPIds[i_m, 0] would be
    the first index in MeshGaussPts for that element and elGPIds[i_m, 1] would
    be the number of Gauss pts for this element.  (?)
    We need this for averaging per element or something like that...
    We also need the element index for a particular GaussPt index (the global
    index in pts = MeshGaussPts.getPoints()). Could be an array (vector) to be
    created by a particular method MeshGaussPts.getElemIds() => vector with
    nb-of-Gauss-pts index items.

    We might also need a filter which Gauss-pts are part of the topo.
    This would be necessary if the getSubTopo() method is expected to accept
    arbitrary indexes. Then it might be that we have elements with only some
    of their Gauss pts being selected as part of the sub-topo. Or we throw an
    error if getSubTopo gets indexes that would tear apart elements.

    @ivar mesh: is a HomoMesh or a HeteroMesh object. The same object can be
        used (referenced) by other MeshElems and MeshGaussPts objects.
    @cvar position: "elemIP". For L{FieldsCollection}.topo objects the
      position attribute identifies the position of values ("node", "element",
      "elemIP", "point" or "structPt").
    """
    position = "elemIP"

    def getPoints(self):
        raise NotImplementedError("Not yet there...")

    def __len__(self):
        raise NotImplementedError("Not yet there...")

    def getSubTopo(self, index):
        """Service method for FieldsCollection.__getitem__() selecting from the
        space dimension.

        Returns a new L{MeshGaussPts} object with a new mesh object of only the
        specified elements.

        @param index: A serial integer index --anything suitable to have numpy
           choose from one dimension.
        """
        implementation_notes = """
  If the getSubTopo() method is expected to accept arbitrary indexes then
  MeshNodes needs some kind of filter array: which of the nodes are actually
  part of the topology. Because then it might be that we have elements (which
  can only be part of the mesh with all their nodes) with only some of their
  nodes being selected as part of the sub-topo. Same problem for MeshGaussPts
  but there we might be able to connect this to the problem of grouping
  (assigning) Gauss pts to elements.

  Or we throw an error if getSubTopo gets indexes that would tear apart
  elements.
        """
        raise NotImplementedError(
            "getSubTopo not yet implemented for class MeshNodes")
