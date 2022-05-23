"""bae.rhino_02.fearo module

This module supplies service functions for fearos (Finite Element Augmented
Rhino Objects) or more generally speaking Rhino files that contain mesh objects
that repressent FE mesh structures.
"""

###############################################################################
###############################################################################
###############################################################################

__version__ = "1.01"

_version_history_ = """
Versions:
=========

1.0 : GP new
1.01: GP incompatible change from createShrinkElsets to createShrinkElset
"""

from bae.vecmath_01 import vector_plus, vector_scale

# shrinkedTetFaces defines how the Rhino shrinked tet object is created from
# the four corner nodes of an Abaqus tet element.
# shrinkedTetFaces contains the node indexes for the four faces of the shrinked
# tet object
# This is the inverse to cornerNodeFaceVertIds below.
shrinkedTetFaces = [ [0,2,1], [0,1,3], [0,3,2], [1,2,3] ]

def createShrinkElset(mesh, elset, shrinkfactor=0.05):
    """Determine properties needed to generate Rhino shrinked tet mesh object
    for a particular set of elements.

    The items returned can directly be passed on to the
    L{bae.rhino_02.rhino_binary.Mesh} constructor:
     >>> from bae.abq_model_02 import Model
     >>> from bae.rhino_02.fearo import createShrinkElset
     >>> from bae.rhino_02.rhino_binary import Mesh
     >>> from bae.rhino_02 import RhinoModel
     >>> model = Model().read("myAbqModel.inp")
     >>> (vertices, faces, texCoords) = createShrinkElset(model, "STOPE")
     >>> rhMesh = Mesh(vertices, faces, texCoords,
     >>>               objectName="someTets", layerIndex=0)
     >>> rhModel = RhinoModel()
     >>> rhModel.model.objects.append(rhMesh)
     >>> rhModel.write("shrinked.3dm")

    @param mesh: A L{bae.mesh_01.Mesh} object or a L{bae.abq_model_02.Model}
      object containing nodes, elements and (possibly) elsets.
    @param elset: set of element numbers or the name of an elset (a string).
      If it's a string (an elset name) then the mesh argument must be a
      L{bae.abq_model_02.Model} object.
    @param shrinkfactor: ratio by which the tets in Rhino will be smaller than
      in the Abaqus mesh.
      I.e. lengthInRhino = (1.0-tetShrinkfactor) * lengthInAbaqus

    @returns: a (vertices, faces, texCoords)-tuple.
    """

    if isinstance(elset, basestring):
        elset = mesh.elset[elset]

    vertices = list()
    faces = list()
    texCoords = list()

    elCentroid = mesh.getElCentroids(elset)
    baseVertId = 0
    for el in elset:

        centroid = elCentroid[el]
        shrCentr = vector_scale(centroid, shrinkfactor)
        elNodes = mesh.elNodes[el][:4]  # just the corner nodes

        # vertex coords
        vertices.extend(
            vector_plus(
                vector_scale(mesh.nodeCoords[node], 1.0-shrinkfactor),
                shrCentr)
            for node in elNodes)

        # element and node number in texture coordinates
        texCoords.extend(
            (el, node)
            for node in elNodes)

        # add four faces
        faces.extend(
            [i+baseVertId for i in face]
            for face in shrinkedTetFaces)

        # four vertices added, prepare next base id
        baseVertId += 4

    # return properties of this Rhino mesh object
    return (vertices, faces, texCoords)
