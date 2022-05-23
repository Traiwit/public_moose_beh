"""package to access Rhino 3dm files

Intended to be imported like a module.

Usage
=====
 >>> from bae.rhino_02 import RhinoModel
 >>> m = RhinoModel()
 >>> help(m) for more information

 >>> from bae.surface_03 import TriangleSurface
 >>> from bae.rhino_02 import RhinoModel
 >>>
 >>> surf = TriangleSurface.fromSTL("torus.stl")
 >>> m = RhinoModel()
 >>> m.layers.append(["torus-layer", (255,10,10), None])
 >>> m.importTriangleSurface(surf, objName="torus", layerIdx=0)
 >>> m.write("torus_out.3dm")
"""

__version__ = "2.3"

_version_history_ = """
rhino_01.py was unversioned and is completely different: It was meant to
create a script to be pasted to the Rhino command line. We're not using
this technique anymore so the rhino-name was captured for this new module.

2.0 GP new
2.1 GP change from fearo.createShrinkElsets to createShrinkElset,
       added importAbqModelElsetAsShrinkTets
2.2 GP added ColorsIter-class and objects
2.3 GP added RhinoModel.addPointCloud
"""

from itertools import izip, cycle as itertools_cycle

# local import from submodules of this package
import rhino_binary as binary
from fearo import createShrinkElset

from bae.log_01 import msg

#-----------------------------------------------------------------------------

class RhinoModel(object):
    """
    Usage:
    ======
     >>> from bae.rhino_02 import RhinoModel
     >>> from bae.abq_model_02 import Model as AbqModel
     >>> from bae.surface_03 import TriangleSurface
     >>>
     >>> model = AbqModel().read("FromIcem.inp")
     >>> rhModel = RhinoModel()
     >>> allS3 = model.typeEl["S3"]
     >>> for cnt, eln in enumerate(sorted(model.elset)):
     >>>     els = model.elset[eln].intersection(allS3)
     >>>     if not els:
     >>>         continue
     >>>     rhModel.layers.append([eln, (0,0,0), None])
     >>>     surf = TriangleSurface.fromMeshTris(m, els)
     >>>     rhModel.importTriangleSurface(surf, objName=eln, layerIdx=cnt)
     >>> rhModel.write("FromIcem_onlyS3.3dm")

    @ivar layers: list of layers, each layer being a [name, colour,
       parent_layer_id]-list.

       colour is an RGB-tuple of integers between 0 and 255.

       parent_layer_id if it is given and an integer then it's the index of
       the parent layer in this list. Otherwise the layer has no parent layer.
       (I.e. specifying None --or other non-integers-- means: no parent.)

       Note that parent layers must always be defined earlier in the list than
       their children.

    @ivar objects: list of geometry objects, for exampole L{rhino_binary.Mesh}
    @ivar model: the L{rhino_binary.RhinoModel} object (internal use)
    """

    def __init__(self):
        self.model = binary.RhinoModel()
        self.layers = self.model.layers
        self.objects = self.model.objects

    def write(self, filename, version=0):
        self.model.write(filename, version)

    # def read(self, filename):
    #     pass

    def importAbqModelAsShrinkTets(
            self, mesh, elsetNames=None, shrinkfactor=0.05):
        """Add shrinked tets mesh object to the RhinoModel self.

        @param mesh: a L{bae.abq_model_02.Model} object
        @param elsetNames: List of elset names to be imported into the self.
          Only those elsets will be imported. If not given then all elsets in
          the FE model will be imported. Each elset corresponds to one mesh
          object in the RhinoModel.
        @param shrinkfactor: ratio by which the tets in Rhino will be smaller
          than in the Abaqus mesh.
          I.e. lengthInRhino = (1.0-tetShrinkfactor) * lengthInAbaqus

        @Note: You might want to use
          L{model.regexpNames("elset", ...)<bae.abq_model_02.Model.regexpNames>}
          to get a list of elset names from a regular expression.
        """

        if elsetNames is None:
            elsetNames = sorted(mesh.elset)

        cnt = 0
        for eln in elsetNames:
            self.importAbqModelElsetAsShrinkTets(
                mesh, eln, objectName=eln, shrinkfactor=shrinkfactor)
            cnt += 1
        msg("Added %d shrinked tets objects." % cnt, debugLevel=1)
        return

    def importAbqModelElsetAsShrinkTets(
            self, mesh, elset, objectName="", layerIdx=0, shrinkfactor=0.05):
        """Add shrinked tets mesh object to the RhinoModel self.

        @param mesh: A L{bae.mesh_01.Mesh} object or a L{bae.abq_model_02.Model}
          object containing nodes, elements and (possibly) elsets.
        @param elset: set of element numbers or the name of an elset (a string).
          If it's a string (an elset name) then the mesh argument must be a
          L{bae.abq_model_02.Model} object.
        @param objectName: Optional Rhino object name.
        @param layerIdx: Optional Layer index. If specified (and not 0) then
           make sure there is a corresponding item in L{self.layers}.
        @param shrinkfactor: ratio by which the tets in Rhino will be smaller
          than in the Abaqus mesh.
          I.e. lengthInRhino = (1.0-tetShrinkfactor) * lengthInAbaqus

        @returns: The new L{bae.rhino_02.binary.Mesh}-object (self.objects[-1]).
           Or None if the elset contains no elements.
        """
        (vertices, faces, texCoords) = createShrinkElset(
            mesh, elset, shrinkfactor)
        if vertices:
            rhMesh = binary.Mesh(
                vertices, faces, texCoords, objectName, layerIdx)
            self.model.objects.append(rhMesh)
            msg("Added elset as shrinked tets object, objectName %s, layerIdx"
                " %d." % (objectName, layerIdx), debugLevel=10)
        else:
            rhMesh = None
            msg("Elset contains no elements for objectName %s, layerIdx %d."
                % (objectName, layerIdx), debugLevel=10)
        return rhMesh

    def importTriangleSurface(
            self, surface, objName="", layerIdx=0, welded=False):
        """Add a TriangleSurface object as mesh to the RhinoModel self.
        Node and element nubers are stored in texture coordinates.

         >>> surf = TriangleSurface...
         >>> model = RhinoModel()
         >>> rhMesh = model.importTriangleSurface(surf, ...)
         >>> rhMesh.objectLayerIndex = 2  # later change the layer
         >>> rhMesh.objectName = "CAVE"
         >>> rhMesh.setUserString("key", "value")

        @param surface: a L{bae.surface_03.TriangleSurface} object
        @param objName: Optional Rhino object name.
        @param layerIdx: Optional Layer index. If specified (and not 0) then
           make sure there is a corresponding item in L{self.layers}.
        @param welded: if True then create a mesh object with each vertex /
           node stored only once: a welded mesh. Otherwise each node /
           vertex is stored once for each facet / triangle / cohesive
           element. Node and element numbers will only be stored if this
           argument is false.

        @returns: The new L{bae.rhino_02.binary.Mesh}-object (self.objects[-1]).
           Or None if the surface contains no elements.
        """

        if welded:

            # find all connected nodes
            allnodes = set()
            for nodes in surface.elNodes.itervalues():
                allnodes.update(nodes)
            allnodes = sorted(allnodes)

            # vertices
            vertices = list()
            nodeToPtIdx = dict()
            for cnt, node in enumerate(allnodes):
                vertices.append(surface.nodeCoords[node])
                nodeToPtIdx[node] = cnt

            # faces
            faces = [ [nodeToPtIdx[node] for node in nodes]
                      for nodes in surface.elNodes.itervalues() ]

            texCoords = None

        else:

            nodeIdx = 0
            vertices = list()
            faces = list()
            texCoords = list()

            for elem, nodes in surface.elNodes.iteritems():
                nbnodes = len(nodes)
                facenodes = range(nodeIdx, nodeIdx+nbnodes)

                for i, node in izip(facenodes, nodes):
                    # vertices
                    vertices.append(surface.nodeCoords[node])
                    # element and node number in texture coordinates
                    texCoords.append([elem, node])

                # add face to faces list
                faces.append(facenodes)

                nodeIdx += nbnodes

        if vertices:
            rhMesh = binary.Mesh(vertices, faces, texCoords, objName, layerIdx)
            self.model.objects.append(rhMesh)
        else:
            rhMesh = None
        return rhMesh

    def addPointCloud(self, points, normals=None, objectName="", layerIdx=0):
        """Add a pointcloud object to the RhinoModel self.

        @param normals: optional list of normals (triples of floats).

        @returns: the PointCloud object. You can change layer and object name
          later.
        """
        rhPtCloud = binary.PointCloud(points, normals)
        del points, normals  # save memory

        rhPtCloud.objectLayerIndex = layerIdx
        rhPtCloud.objectName = objectName

        self.model.objects.append(rhPtCloud)
        return rhPtCloud
        
    
    def importAbqModelAsPointCloud(self, mesh):
        """Add tet mesh as pointcloud object to the RhinoModel self.

        Usage:
         >>> model = RhinoModel()
         >>> ptCloud = model.importAbqModelAsPointCloud(mesh)
         >>> ptCloud.objectLayerIndex = 2  # later change the layer
         >>> ptCloud.objectName = "G02"
         >>> ptCloud.setUserString("key", "value")

        @param mesh: a L{bae.abq_model_02.Model} object

        @returns: the PointCloud object. You can change layer and object name
          later.

        @Note: Use L{importAbqModelElsetsAsPointClouds} to import elsets as
          pointcloud objects.
        """

        elCentroids = mesh.getElCentroids()
        if not elCentroids:
            return None

        elements = sorted(elCentroids)
        points = [elCentroids[e] for e in elements]
        normals = [[0.0, 0.0, float(e)] for e in elements]
        rhPtCloud = binary.PointCloud(points, normals)
        del points, normals  # save memory
        self.model.objects.append(rhPtCloud)
        msg("Added %d element centroids as point cloud." % len(elements),
            debugLevel=1)
        return rhPtCloud

    def importAbqModelElsetsAsPointClouds(
            self, mesh, elsetNames=None):
        """Add elset centroids as pointcloud objects to the RhinoModel self.

        @param mesh: a L{bae.abq_model_02.Model} object

        @param elsetNames: List of elset names to be imported into self.
          Only those elsets will be imported. Each elset corresponds to one
          mesh object in the RhinoModel.

          If not given then all elsets will be imported.

        @Note: Use L{importAbqModelAsPointCloud} to import a complete mesh as
          one pointcloud object.

        @Note: You might want to use
          L{model.regexpNames("elset", ...)<bae.abq_model_02.Model.regexpNames>}
          to get a list of elset names from a regular expression.
        """

        if elsetNames is None:
            elsetNames = sorted(mesh.elset)

        elCentroids = mesh.getElCentroids(elsetNames)
        
        cnt = 0
        for eln in elsetNames:
            elements = sorted(mesh.elset[eln])
            points = [elCentroids[e] for e in elements]
            normals = [[0.0, 0.0, float(e)] for e in elements]
            rhPtCloud = binary.PointCloud(points, normals, objectName=eln)
            self.model.objects.append(rhPtCloud)
            cnt += 1
        msg("Added %d shrinked tets objects." % cnt, debugLevel=1)

#-----------------------------------------------------------------------------

class ColorsIter(itertools_cycle):
    r"""Loop infinitely over a set of rgb colors.

    Usage of the predefined ColorsIter instances as layer colours:

    >>> from bae.rhino_02 import RhinoModel, \
    >>>     colorsIterRhinoDefault as layerColors
    >>>
    >>> for name, col in izip(layerNames, layerColors):
    >>>     m.layers.append([name, col, None])


    How to define another colour sequence:

    >>> colorsRGB = ColorsIter([(255,0,0),(0,255,0),(0,0,255)])
    """
    pass

# ColorsIter: brown, red, orange, green, dark green, cyan, blue, purple,
# magenta, pink, white
colorsIterBrownToWhite = ColorsIter([
    (127,  0,0),   # brown
    (255,  0,0),   # red
    (255,127,0),   # orange
    (  0,255,0),   # green
    (  0,127,0),   # dark green
    (  0,255,255), # cyan
    (  0,  0,255), # blue
    (127,  0,255), # purple
    (255,  0,255), # magenta
    (255,191,191), # pink
    (255,255,255), # white
])

# ColorsIter: Colours from the Rhino Colour dialog except yellow
colorsIterRhinoDefault = ColorsIter([
    (255, 0, 0), (191, 63, 63), (255, 127, 0), (255, 191, 0),
    (127, 255, 0), (0, 255, 0), (0, 127, 0), (63, 191, 127),
    (127, 255, 191), (0, 255, 255), (63, 191, 191), (191, 191, 255),
    (0, 0, 255), (0, 0, 191), (191, 63, 255), (255, 0, 255),
    (255, 127, 255), (255, 191, 191),
])

# ColorsIter: Colours from the default initial five layers
colorsIterRhinoInitialLayers = ColorsIter([
    (200,0,0), (125,38,205), (0,0,255), (0,127,0), (255,255,255),
])
