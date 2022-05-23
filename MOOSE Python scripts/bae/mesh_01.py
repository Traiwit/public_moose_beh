"""Module to store mesh information (e.g. nodes and elements for unstructured
meshes and also mesh data for structured grids.).
"""

__version__ = "1.46"

_version_history_ = """
Versions:
=========

1.0 : GP initial setup. derived from abq_model_02 v2.4 and VolumeFromMesh
1.1 : GP added: renumber(), rotateZ(), findNodes()
1.2 : GP added: rotateX(), rotateY(), rotate()
1.3 : GP added: MeshStructuredPoints (pointsgrid)
         changed: use bae.log_01 now
1.4 : GP added: Mesh.getClosestElemCoords(), Mesh.checkElVolPos()
1.5 : GP added: Mesh.warp()
1.6 : GP added: Mesh.getNodalAveragingIPCoeffs() for nodal interpolation
1.7 : GP added: new class MeshUnstructuredPoints
1.8 : GP added: MeshStructuredPoints.getGridIdxIter() and len()
1.9 : GP fixed: getClosestElemCoords() made use of bug in shapeEl.__getitem__
                that's been corrected in the meantime
1.10: GP added: new class InterpolMeshToPoints,
         added MeshStructuredPoints.__cmp__()
1.11: GP changed: Mesh.partMeshFromElements() with fewer diagnostic output and
               more generator extressions. Hoping to speed things up.
             - new .interpol (pickle) format, weighting factors not stored
               .interpol-format stored in the file for compatibilty check
             - using msg.tick() more often
             - in some places: nbPoints now really an optional argument
1.12: GP added new interpolation format v2.0
1.13: GP added argument shapesConsidered to Mesh.initializePointSearch()
         defaulting to tet-only shapes. (I.e. default is to only look for
         points in tet elements.)
1.14: GP added Mesh.read(), Mesh.write() and co.
1.15: GP changed logfile argument to Mesh() is now being ignored
1.16: GP added InterpolMeshToPoints.writeRemapParam
1.17: GP incompatible change: renamed getMeshCallback argument to
         getMeshCallBack for consistency in InterpolMeshToPoints.__init__
1.18: GP added MeshStructuredPointsRot
1.19: GP fixed InterpolMeshToPoints._gridToFile to write floats with
         sufficient precision (using "%r" instead of "%g")
1.20: GP added getMesh() function
1.21: GP added MeshStructuredPoints.__str__() and
         MeshStructuredPointsRot.__str__()
1.22 GP fixed MeshStructuredPointsRot.getPointsIter()
1.23 GP added getConnectedNodes()
1.24 GP added __eq__ and __ne__ to MeshUnstructuredPoints.
    To prevent interpolation pickle files from becoming invalid due to rounding
    errors.
1.25 GP added make writeRemapParam optionally write to strings,
        modified getMesh to also return MeshUnstructuredPoints
        incompatibel API change: getMesh now takes kwargs not a dict anymore
1.26 GP added getElemCoordsForGrid, use this instead of old getElemCoords
1.27 GP changed InterpolMeshToPoints.writeRemapParam(): works with different
        element types now, before: only 10-node 4-IP quad tets.
1.28 GP added notDefinedKey to writeRemapParam
1.29 GP for MeshStructuredPoints added getGridIdxFromSerialIndex,
        getClosestGridPointIdx, getPointAlignedGrid, getPoint accepts neg index
        and fixed bug in getPtIdsWeights
1.30 GP fixed getElemVolField() for thick wedges with specified thickness
1.31 GP added MeshStructuredPoints.getIdsInIntervalOverZ for volume_02 classes
1.32 GP added MeshStructuredPointsRot.getPoint accepts single negative index.
1.33 AF added splitElset_DisjointParts
1.34 GP fixed getIdsInBox() for integer grid dimensions and box values
1.35 GP fixed InterpolMeshToPoints._tetElemCoordsToFile(): earlier if the first
        grid point was outside the mesh then all element coords got zero length
1.36 TR let getMesh call MeshStructuredPointsRot for any rotation related kwargs
1.37 GP added InterpolMeshToCentroids, InterpolMeshToFaceCentroids (and
    InterpolMeshBaseClass), added interpolation data for wedge elements,
    made elShapeFunct for TET_Q use all four element coords and restructured
    this adding constant dicts elShapeElemShapeFunc. Added dicts
    elShapeFaceCentrCoords and elShapeFaceCentrWeights.
1.38 GP changed Mesh.shapeEl into an object of the new type ShapeElDict that
    also finds the general shapes like "TET"
        added Mesh.getExteriorSurface
1.39 GP added __str__
1.40 GP changed Gauss point locations from text book values to measured Abaqus
    element C3D10M values.
1.41 GP added Mesh.getTopSurface()
1.42 GP changed InterpolMeshToFaceCentroids arguments, prefered now:
    ElemFaceSurface
1.43 GP added: InterpolMeshToPoints now checks if the directory for the
    interpolation data pickle file exists and if not creates it.
1.44 GP added: getGridDataString() for MeshStructuredPoints and
    MeshStructuredPointsRot.
1.45 GP added new elshape TETHT_L for DC3D4 linear heat tranfer tet element with
    four Gauss points
1.46 GP generalized getIdsInIntervalOverZ to getIdsInInterval
"""

todo = """
derivatives of nodal fields:
dU_i / dx_j
   = sum_a U_ai dN_a / dx_j
   = sum_a U_ai  dN_a / d xi_k  J^-1 _kj
J_ij := dx_j / d xi_i
   = sum_a X_aj  dN_a / d xi_i
Notes:
 - a is the node index, N_a is the shape function that equals one at node a
   and zero at all others. U_ai is nodal displacement at node a component i.
   xi_j is the element coordinate j.
 - J_ij for quadratic elements is not constant. Inverse needs to be calculated
   separately for each point. dN_a / d xi_j is not constant either.
 - What is the derivative of a Gauss-Pt field? Might be problematic.

for v2.0:
- remove deprecated
- rename getElemCoordsForGrid to getElemCoords
- change InterpolMeshToPoints-interface (?)
  . constructor does nothing (if called without arguments)
  . readFromFile, writeToFile, calculate - methods
  . autoConstruct method: does what __init__ does today
- rename classes:
  . base class = Topography or better Topo (some call it Topology, so we don't
    have to argue about that...)
  . Mesh -> UnstructMesh(Topo)
- Make Mesh.getExteriorSurface return an ElemFaceSurface because this would be
  more general and easy to convert to TriangleSurface (or have an option).
- Move rejoinSplitMesh from bae.utils.splitting_01 here, maybe as method of
  Mesh? Currently we can only import it when the function is being called
  to prevent cyclic import. It seems to me that it belongs here. Maybe the
  whole splitting business should be an integral part of the Mesh class?
"""


from math import sqrt, pi, sin, cos, ceil, floor
from itertools import izip, count
from collections import defaultdict
import re
from array import array
import cPickle as pickle
from struct import Struct
from cStringIO import StringIO

from bae.vecmath_01 import vector, vector_plus, vector_minus, vector_sum,\
    vector_modif_add, vector_scale, cross, dot, dist, length,\
    mat_multvec, mat_determinant, mat_transpose, mat_orthoNormBaseGS,\
    mat_toRodrigues, trans_rodrigues
from bae.misc_01 import BoundingBox, KDTree

from bae.log_01 import msg, MsgTicker, log

# ... from bae.misc_01 import ... , selectattrversion
# try:
#     import numpy as np
# except ImportError:
#     np = None


#---------------------------------------------------------
# Error handling
#---------------------------------------------------------
class ModelInconsistent(Exception): pass
class UnknownElShape(Exception): pass
class ZeroOrNegVolumeElements(Exception): pass
class BadPickleFile(Exception): pass


def getMesh(*args, **kwargs):
    """Return a suitable object of the appropriate subclass of L{MeshBaseType}
    depending on the type and contents of the arguments.

    To get a L{MeshStructuredPoints} or L{MeshStructuredPointsRot} you must
    supply keyword arguments. In this case any valid combination of arguments
    for the subclass that you wish to instantiate must be provided:
     >>> gridData = dict(firstPoint=[...], lastPoint=[...], ...)
     >>> grid = getMesh(**gridData)

    If a positional argument is supplied it is interpreted as list of points
    and a corresponding L{MeshUnstructuredPoints} object is returned:
     >>> points = [[0,0,0], [0.5,1,0.5], ...]
     >>> grid = getMesh(points)
    """
    if args:
        # it's assumed to be a list of points
        return MeshUnstructuredPoints(args[0])
    elif any(kw in kwargs for kw in [
            "box", "firstPoint", "origin", "originBackRotated"]):
        # it's one of the MeshStructuredPoints classes
        if any(kw in kwargs for kw in [
            "rotString","originBackRotated","deg",
            "baseVectors","edgeVectors","rotVector"]):
            return MeshStructuredPointsRot(**kwargs)
        else:
            return MeshStructuredPoints(**kwargs)
    else:
        raise NotImplementedError(
            "Don't know what kind of mesh to create from these parameters:"
            " args: %s; kwargs: %s"
            % (args, kwargs))


#---------------------------------------------------------
# mesh base class, service classes
#---------------------------------------------------------
class MeshBaseType(object):
    """Base class for different mesh types
    """
    pass


class ShapeElDict(dict):
    """Variant of dict that will also find elements for tghe general shape
    like "TET" if you have a corresponding specific / separate shape like
    "TET_L" or "TET_Q" in the mesh.
    """

    def __getitem__(self, shape):
        """Get set of element numbers of the specified shape.
        Both specifically linear or quad versions (like "TET_L") as well as
        the general version (like "TET") might be queried.

        @raises KeyError: if this shape is not present in the mesh.
        """
        if "_" in shape:
            # shape is a specific linear or quad version
            return dict.__getitem__(self, shape)
        else:
            # shape is a general version
            res = dict.get(self, "%s_L" % shape, set()).union(
                dict.get(self, "%s_Q" % shape, set()))
            if not res:
                raise KeyError(
                    "No elements neither of shape %s_L nor %s_Q exist."
                    % (shape, shape))
            return res

    def get(self, shape, default=None):
        """Get set of element numbers of the specified shape.
        Both specifically linear or quad versions (like "TET_L") as well as
        the general version (like "TET") might be queried.
        """
        if "_" in shape:
            # shape is a specific linear or quad version
            return dict.get(self, shape, default)
        else:
            # shape is a general version
            res = dict.get(self, "%s_L" % shape, set()).union(
                dict.get(self, "%s_Q" % shape, set()))
            if res:
                return res
            else:
                return default


#---------------------------------------------------------
# mesh class for unstructured meshes like abaqus meshes
#---------------------------------------------------------
class Mesh(MeshBaseType):
    """Class to hold the data for an unstructured mesh of possibly varying
    element types.

    @ivar nodeCoords: dict {node number: [x,y,z]},
    @ivar elNodes: dict {element number: [node numbers]}, B{do not modify
        directly},
        use L{updateElems()<bae.mesh_01.Mesh.updateElems>},
        L{removeElems()<bae.mesh_01.Mesh.removeElems>}.
    @ivar elShape: dict {element number: element type string}, B{do not modify
        directly},
        use L{updateElems()<bae.mesh_01.Mesh.updateElems>},
        L{removeElems()<bae.mesh_01.Mesh.removeElems>}.
    @ivar shapeEl: {element type string: set of element numbers} of type
        L{ShapeElDict}, will retrieve elements for general shapes like "TET"
        as well,
        B{do not modify directly},
        use L{updateElems()<bae.mesh_01.Mesh.updateElems>},
        L{removeElems()<bae.mesh_01.Mesh.removeElems>}.

    @ivar centroidsKDTree: a L{KDTree<bae.misc_01.KDTree>}-object initialized
        with all tet element centroids. This attribute is not there initially.
        It is being initialized (and used) by getClosestElemCoords(). Can
        savely be removed to save memory. It's being automatically recreated
        when needed by getClosestElemCoords().

    @Note: If you import the abq_model_02 module you will get a
    Mesh.asAbqModel() method to convert the mesh to an Abaqus model, e.g for
    output into an Abaqus input file.

    @Note: See L{getElemShapeFuncs} for a description of the element
    coordinates convention.
    """

    #: Node indexes for each face of a particular element shape.
    #: The index ordering of each face is such that each faces normal
    #: points inwardly.
    #: The general shapes like "TRI", "TET", ... have only the corner nodes.
    elShapeFaceNodes = {
        'TRI_L': ((0,1), (1,2), (0,2)),
        'QUAD_L': ((0,1), (1,2), (2,3), (3,0)),
        'TET_L': ((0,1,2), (1,0,3), (2,1,3), (0,2,3)),
        'TETHT_L': ((0,1,2), (1,0,3), (2,1,3), (0,2,3)),
        'WEDGE_L': ((0,1,2), (3,5,4), (0,3,4,1), (1,4,5,2), (0,2,5,3)),
        'HEX_L': ((0,1,2,3), (4,7,6,5), (0,4,5,1), (1,5,6,2), (2,6,7,3),
                  (0,3,7,4)),

        'TRI_Q': ((0,1,3), (1,2,4), (0,2,5)),
        'QUAD_Q': ((0,1,4), (1,2,5), (2,3,6), (3,0,7)),
        'TET_Q': ((0,1,2,4,5,6), (1,0,3,4,7,8), (2,1,3,5,8,9), (0,2,3,6,9,7)),
        'WEDGE_Q': ((0,1,2,6,7,8), (3,5,4,9,10,11), (0,3,4,1,12,9,13,8),
                    (1,4,5,2,13,10,14,7), (0,2,5,3,8,14,11,12)),
        'HEX_Q': ((0,1,2,3,8,9,10,11), (4,7,6,5,12,13,14,15),
                  (0,4,5,1,16,11,17,8), (1,5,6,2,17,13,18,9),
                  (2,6,7,3,18,14,19,10), (0,3,7,4,11,19,15,16)),

        'TRI': ((0,1), (1,2), (0,2)),
        'QUAD': ((0,1), (1,2), (2,3), (3,0)),
        'TET': ((0,1,2), (1,0,3), (2,1,3), (0,2,3)),
        'TETHT': ((0,1,2), (1,0,3), (2,1,3), (0,2,3)),
        'WEDGE': ((0,1,2), (3,5,4), (0,3,4,1), (1,4,5,2), (0,2,5,3)),
        'HEX': ((0,1,2,3), (4,7,6,5), (0,4,5,1), (1,5,6,2), (2,6,7,3),
                (0,3,7,4)),
        }

    #: mid nodes for quadratic elements on edge, a dict {shape : edges list}
    #: shape is something like 'TRI', 'TET', ...
    #: the first item of the edges list is a list of node indices identifying
    #: the corner nodes between which the first quadratic mid node lies and so
    #: on
    quadMidNodeEdges = {
        'TRI': ((0,1),(1,2),(0,2)),
        'TET': ((0,1),(1,2),(0,2),(0,3),(1,3),(2,3)),
        'TETHT': ((0,1),(1,2),(0,2),(0,3),(1,3),(2,3)),
        'WEDGE': ((0,1),(1,2),(0,2),(3,4),(4,5),(3,5),
                  (0,3),(1,4),(2,5),),
        }

    #: number of corner nodes for each of the element shapes
    elShapeToNbCornerNodes = dict(
        (t+ext, nb)
        for t, nb in [
                ('LINE', 2),
                ('TRI', 3),
                ('QUAD', 4),
                ('TET', 4),
                ('TETHT', 4),
                ('WEDGE', 6),
                ('HEX', 8),
                ]
        for ext in ("_L", "_Q", "") )

    #: Element shape functions: given element coordinates return the shape
    #: functions (i.e. nodal weights)
    #: Element coordinates are:
    #:  - for tet elements: four barycentric coordinates, the sum of all four
    #:    always equals one. Inside the tet element each of the coordinates is
    #:    between 0 and 1. They are identical to the shape functions of the
    #:    linear tet element.
    #:  - for wedge elements: three barycentric coordinates on the top and
    #:    bottom triangle and the forth coordinate going from -1 to 1 from
    #:    bottom (node 1..3) to top (node 4..6).
    #:  - for bricks/HEX elements: three coordinates going from -1 to 1.
    elShapeElemShapeFunc = {
        "TET_L": lambda c:c,
        "TETHT_L": lambda c:c,
        "TET_Q": lambda c:[
            2.*c[0]**2-c[0],
            2.*c[1]**2-c[1],
            2.*c[2]**2-c[2],
            2.*c[3]**2-c[3],
            4.*c[0]*c[1],
            4.*c[1]*c[2],
            4.*c[0]*c[2],
            4.*c[0]*c[3],
            4.*c[1]*c[3],
            4.*c[2]*c[3]],
        # FYI only: alternative formulation using only the last three element
        # coords[1:4] == [xi, eta, zeta]:
# node  1.0, xi, eta, zeta, xi**2, eta**2, zeta**2, xi*eta, xi*zeta, eta*zeta
# 1   [[ 1.  -3.  -3.  -3.    2.     2.       2.      4.      4.       4.]
# 2    [ 0.  -1.   0.   0.    2.     0.       0.      0.      0.       0.]
# 3    [ 0.   0.  -1.   0.    0.     2.       0.      0.      0.       0.]
# 4    [ 0.   0.   0.  -1.    0.     0.       2.      0.      0.       0.]
# 5    [ 0.   4.   0.   0.   -4.     0.       0.     -4.     -4.       0.]
# 6    [ 0.   0.   0.   0.    0.     0.       0.      4.      0.       0.]
# 7    [ 0.   0.   4.   0.    0.    -4.       0.     -4.      0.      -4.]
# 8    [ 0.   0.   0.   4.    0.     0.      -4.      0.     -4.      -4.]
# 9    [ 0.   0.   0.   0.    0.     0.       0.      0.      4.       0.]
# 10   [ 0.   0.   0.   0.    0.     0.       0.      0.      0.       4.]]
        "WEDGE_L": lambda c:[
            0.5*c[0]*(1-c[3]),
            0.5*c[1]*(1-c[3]),
            0.5*c[2]*(1-c[3]),
            0.5*c[0]*(1+c[3]),
            0.5*c[1]*(1+c[3]),
            0.5*c[2]*(1+c[3])],
        # "WEDGE_Q":
        # Not implemented because it's not clear which form we need
        # see Abaqus Theory Guide, 3 Elements, 3.2 Continuum elements,
        # 3.2.6 Triangular, tetrahedral, and wedge elements.
        # Note that:
        # coords[0]==(1-g-h), coords[1]==g, coords[2]==h, coords[3]==r
        "HEX_L": lambda c:[
            0.125*(1-c[0])*(1-c[1])*(1-c[2]),
            0.125*(1+c[0])*(1-c[1])*(1-c[2]),
            0.125*(1+c[0])*(1+c[1])*(1-c[2]),
            0.125*(1-c[0])*(1+c[1])*(1-c[2]),
            0.125*(1-c[0])*(1-c[1])*(1+c[2]),
            0.125*(1+c[0])*(1-c[1])*(1+c[2]),
            0.125*(1+c[0])*(1+c[1])*(1+c[2]),
            0.125*(1-c[0])*(1+c[1])*(1+c[2]),
            ],
        }

    #: position of centroid in element coordinates
    shapeCentrCoords = dict(
        (t+ext, c)
        for t, c in [
            ('TET', [0.25, 0.25, 0.25, 0.25]),
            ('TETHT', [0.25, 0.25, 0.25, 0.25]),
            ('WEDGE', [1./3, 1./3, 1./3, 0.0]),
            ('HEX', [0.0, 0.0, 0.0]),
            ]
        for ext in ("_L", "_Q"))

    #: Coefficients for the matrix elShapeGaussPoints['TET_Q']
    #: This is according to some text book.
    _gpTetA_tb = (5+3*sqrt(5))/20.
    _gpTetB_tb = (5-sqrt(5))/20.
    # _gpTetA = _gpTetA_tb
    # _gpTetB = _gpTetB_tb

    #: Coefficients for the matrix elShapeGaussPoints['TET_Q']
    #: These are measured values for Abaqus C3D10M.
    _gpTetA = 0.46875
    _gpTetB = 0.53125/3

    # #: Coefficients for the matrix elShapeGaussPoints['TET_Q']
    # #: These are measured values for Abaqus C3D10  (senza 'M').
    # #: These are rather close (about 1E-8 relative deviation) to the
    # #: "text book" values, see above.
    # _gpTetA = 0.585410177707672
    # _gpTetB = 0.138196602463722  # corrected: 0.138196607430776

    #: position of integration points in element coordinates
    elShapeGaussPoints = {
        'TET_L': ((0.25,0.25,0.25,0.25),),
        'TETHT_L': ((_gpTetA_tb,_gpTetB_tb,_gpTetB_tb,_gpTetB_tb),
                    (_gpTetB_tb,_gpTetA_tb,_gpTetB_tb,_gpTetB_tb),
                    (_gpTetB_tb,_gpTetB_tb,_gpTetA_tb,_gpTetB_tb),
                    (_gpTetB_tb,_gpTetB_tb,_gpTetB_tb,_gpTetA_tb)),
        'HEX_L': ((0,0,0),),
        'WEDGE_L': ((1,0,0,0), (0,1,0,0), (0,0,1,0)),

        'TET_Q': ((_gpTetA,_gpTetB,_gpTetB,_gpTetB),
                  (_gpTetB,_gpTetA,_gpTetB,_gpTetB),
                  (_gpTetB,_gpTetB,_gpTetA,_gpTetB),
                  (_gpTetB,_gpTetB,_gpTetB,_gpTetA)),
        }

    #: Coefficients for the matrix _gpQTetInterpolation
    #: Here are analytically calculated values for the "text book"-locations
    #: of the Gauss points.
    #: I.e. Those values are *NOT* consistent with the actual Gauss point
    #: locations in the Abaqus C3D10M elements!
    #: We might use them anyway to mitigate the effect of overshooting
    #: values close to the elements borders from the linear extrapolation.
    #: And we use them for the linear heat transfer tet elements DC3D4
    _gpTetC_tb = (1+3*sqrt(5))/4.
    _gpTetD_tb = (1-sqrt(5))/4.
    # _gpTetC = _gpTetC_tb
    # _gpTetD = _gpTetD_tb

    #: Coefficients for the matrix _gpQTetInterpolation
    #: Here are values that Abaqus/CAE apparently uses for C3D10M elements.
    #: (Actually measured _gpTetD=-0.44230772395219103, later corrected.)
    #: Those values are *NOT* consistent with the actual Gauss point
    #: locations in the Abaqus C3D10M elements!
    #: These coefficients would be consistent with Gauss point positions
    #: you would find if you used different coefficients for
    #: elShapeGaussPoints["TET_Q"], namely:
    #: _gpTetA=0.5208333656458779,_gpTetB=0.15972221145137402
    #: Using these inconsistent coeficients mitigates the effect of overshooting
    #: values close to the elements borders from the linear extrapolation. But
    #: less so than the values derived from the "text book"-locations.
    _gpTetC = 2.3269228291298654
    _gpTetD = (1-_gpTetC)/3

    # #: Coefficients for the matrix _gpQTetInterpolation
    # #: Alternatively you can use the "correct" consistent values by using the
    # #: formula below.
    # #: Invert 4x4 matrix with values:
    # #: gpTetA on the main diagonal and gpTetB elsewhere.
    # #: Result is gpTetC on the diagonal and getTetD elsewhere:
    # _gpTetD = 1.0/(3.0*_gpTetB - _gpTetA*_gpTetA/_gpTetB - 2.0*_gpTetA)
    # _gpTetC = -_gpTetD*(_gpTetA/_gpTetB + 2.0)

    #: Values and a conversion matrix for linear interpolation of integration
    #: point values in a quadratic tet element:
    #:  - for the point with element coordinates (alias linear shape function
    #:    values) xi_j (j=1..4) the linear interpolation of the values at the
    #:    integration points f_k (k=1..4) equals:
    #:    sum_j sum_k xi_j f_k _gpQTetInterpolation[k][j]
    #:  - mesh.getElemIPInterpolFuncs(elemCoords) yields the weighting factors
    #:    for integration point k: sum_j xi_j _gpQTetInterpolation[k][j]
    _gpQTetInterpolation = ((_gpTetC,_gpTetD,_gpTetD,_gpTetD),
                            (_gpTetD,_gpTetC,_gpTetD,_gpTetD),
                            (_gpTetD,_gpTetD,_gpTetC,_gpTetD),
                            (_gpTetD,_gpTetD,_gpTetD,_gpTetC))
    _gpQTetInterpolation_tb = ((_gpTetC_tb,_gpTetD_tb,_gpTetD_tb,_gpTetD_tb),
                               (_gpTetD_tb,_gpTetC_tb,_gpTetD_tb,_gpTetD_tb),
                               (_gpTetD_tb,_gpTetD_tb,_gpTetC_tb,_gpTetD_tb),
                               (_gpTetD_tb,_gpTetD_tb,_gpTetD_tb,_gpTetC_tb))

    #: For extrapolating integration point values to the nodes of an element
    #: for qudratic tet elements this is a linear interpolation of the values
    #: at four integration points
    #:  >>> elShapeIPInterpolCoeff = { element shape:
    #:  >>>     [ [ node 1 weight for val at IP 1, IP 2, ... ],
    #:  >>>       [ node 2 ... ], ... ] }
    elShapeIPInterpolCoeff = {
        "TRI_L": [[1.0]]*3,
        "TET_L": [[1.0]]*4,
        "TETHT_L": [  # for DC3D4 lin tets for heat transfer: text book values
            [sum(convVal*linSFVal
                 for convVal, linSFVal in izip(convRow, coords))
             for convRow in _gpQTetInterpolation_tb]
            for coords in [
                [1.0,0.0,0.0,0.0],
                [0.0,1.0,0.0,0.0],
                [0.0,0.0,1.0,0.0],
                [0.0,0.0,0.0,1.0],
                ]],
        "TET_Q": [
            [sum(convVal*linSFVal
                 for convVal, linSFVal in izip(convRow, coords))
             for convRow in _gpQTetInterpolation]
            for coords in [
                [1.0,0.0,0.0,0.0],
                [0.0,1.0,0.0,0.0],
                [0.0,0.0,1.0,0.0],
                [0.0,0.0,0.0,1.0],
                [0.5,0.5,0.0,0.0],
                [0.0,0.5,0.5,0.0],
                [0.5,0.0,0.5,0.0],
                [0.5,0.0,0.0,0.5],
                [0.0,0.5,0.0,0.5],
                [0.0,0.0,0.5,0.5],
                ]],
        "HEX_L": [[1.0]]*8,
        "WEDGE_L": [
            [ 1.0, 0.0, 0.0 ],
            [ 0.0, 1.0, 0.0 ],
            [ 0.0, 0.0, 1.0 ],
            [ 1.0, 0.0, 0.0 ],
            [ 0.0, 1.0, 0.0 ],
            [ 0.0, 0.0, 1.0 ]
            ],
        }
    del coords, convRow

    _f3 = 1./3
    #: location of face centroids given in element coordinates
    elShapeFaceCentrCoords = dict(
        (t+ext, c)
        for t, c in [
                ('TET', {
                    "S1": (_f3,_f3,_f3,0.0),
                    "S2": (_f3,_f3,0.0,_f3),
                    "S3": (0.0,_f3,_f3,_f3),
                    "S4": (_f3,0.0,_f3,_f3),
                    }),
                ('TETHT', {
                    "S1": (_f3,_f3,_f3,0.0),
                    "S2": (_f3,_f3,0.0,_f3),
                    "S3": (0.0,_f3,_f3,_f3),
                    "S4": (_f3,0.0,_f3,_f3),
                    }),

                ("HEX", {
                    "S1": ( 0., 0.,-1.),
                    "S2": ( 0., 0., 1.),
                    "S3": ( 0.,-1., 0.),
                    "S4": ( 1., 0., 0.),
                    "S5": ( 0., 1., 0.),
                    "S6": (-1., 0., 0.),
                    }),
                ("WEDGE", {
                    "S1": (_f3,_f3,_f3,-1.),
                    "S2": (_f3,_f3,_f3, 1.),
                    "S3": ( 0., 0., 1., 0.),
                    "S4": ( 1., 0., 0., 0.),
                    "S5": ( 0., 1., 0., 0.),
                    }),
                ]
        for ext in ("_L", "_Q"))

    #: shape functions of face centroids
    elShapeFaceCentrWeights = dict(
        (t+ext, c)
        for t, c in [
                ('TET', {
                    "S1": (_f3,_f3,_f3,0.0),
                    "S2": (_f3,_f3,0.0,_f3),
                    "S3": (0.0,_f3,_f3,_f3),
                    "S4": (_f3,0.0,_f3,_f3),
                    }),
                ('TETHT', {
                    "S1": (_f3,_f3,_f3,0.0),
                    "S2": (_f3,_f3,0.0,_f3),
                    "S3": (0.0,_f3,_f3,_f3),
                    "S4": (_f3,0.0,_f3,_f3),
                    }),
                ("HEX", {
                    "S1": (0.25, 0.25, 0.25, 0.25, 0.00, 0.00, 0.00, 0.00),
                    "S2": (0.00, 0.00, 0.00, 0.00, 0.25, 0.25, 0.25, 0.25),
                    "S3": (0.25, 0.25, 0.00, 0.00, 0.25, 0.25, 0.00, 0.00),
                    "S4": (0.00, 0.25, 0.25, 0.00, 0.00, 0.25, 0.25, 0.00),
                    "S5": (0.00, 0.00, 0.25, 0.25, 0.00, 0.00, 0.25, 0.25),
                    "S6": (0.25, 0.00, 0.00, 0.25, 0.25, 0.00, 0.00, 0.25),
                    }),
                ("WEDGE", {
                    "S1": (_f3,_f3,_f3, 0., 0., 0.),
                    "S2": ( 0., 0., 0.,_f3,_f3,_f3),
                    "S3": (0.25, 0.25, 0.00, 0.25, 0.25, 0.00),
                    "S4": (0.00, 0.25, 0.25, 0.00, 0.25, 0.25),
                    "S5": (0.25, 0.00, 0.25, 0.25, 0.00, 0.25),
                    }),
                ]
        for ext in ("_L", "_Q"))

    def __init__(self, logfile="default_logfile"):
        """
        @param logfile: DEPRECATED. Ignored. Use bae.log_01
        """

        # data members, they are always there
        self.nodeCoords = dict()
        self.elNodes = dict()
        self.elShape = dict()
        self.shapeEl = ShapeElDict()

        self.point_search_initialized = False


    #{ diagnostic output
    def __str__(self):
        """Some diagnostic output. What's in the mesh?
        """
        res = []

        nb = len(self.nodeCoords)
        if nb:
            res.append("%d nodes" % nb)

        nbElems = len(self.elNodes)
        if nbElems:

            if len(self.shapeEl)==1:
                res.append("%d elements of shape %s"
                           % (nbElems, self.shapeEl.iterkeys().next()))
            else:
                res.append("%d elements in total" % nbElems)

                for shape in sorted(self.shapeEl):
                    elems = self.shapeEl[shape]
                    if not elems:
                        continue
                    res.append("%d elements of shape %s" % (len(elems), shape))

        return ("Unstructured mesh (type %s): %s" % (
            self.__class__.__name__, (", ".join(res))))
    #} end of diagnostic output


    ## ------------------------------------------------------------------------
    #{ file I/O

    def writeNodeCoords(self, pickleFile):
        nodeNbs = sorted(self.nodeCoords)
        msg("Writing mesh: %d node coordinates" % len(nodeNbs), debugLevel=5)
        pickleFile.write("POINTS %d\n" % len(nodeNbs))
        a = array("f", (x for node in nodeNbs
                        for x in self.nodeCoords[node]))
        a.tofile(pickleFile)
        pickleFile.write("NODE_NUMBERS %d\n" % len(nodeNbs))
        a = array("I", nodeNbs)
        a.tofile(pickleFile)

    def readNodeCoords(self, pickleFile):
        """
        @Note: Nodes from the pickle file are added to existing nodes
        overwriting those with identical node numbers.
        """
        # read nodeCoords ... 1. coords
        line = pickleFile.readline().rstrip().upper().split(None)
        if line[0] != "POINTS":
            raise BadPickleFile(
                "Expected POINTS option. Instead found %s." % line)
        nbNodes = int(line[1])
        arrCoords = array("f")
        arrCoords.fromfile(pickleFile, nbNodes*3)
        msg("Reading mesh: got %d node coordinates" % nbNodes,
            debugLevel=5)

        # read nodeCoords ... 2. node numbers
        line = pickleFile.readline().rstrip().upper().split(None)
        if line[0] != "NODE_NUMBERS":
            raise BadPickleFile(
                "Expected NODE_NUMBERS option. Instead found %s." % line)
        if int(line[1]) != nbNodes:
            raise BadPickleFile(
                "Incompatible number of nodes in option POINTS (%d)"
                " and option NODE_NUMBERS (%d)." % (nbNodes, int(line[1])))
        arrNodeNb = array("I")
        arrNodeNb.fromfile(pickleFile, nbNodes)
        msg("Reading mesh: got %d node numbers" % nbNodes, debugLevel=5)

        # update nodeCoords dictionary
        #  - for the izip(*[iter]*N)-phrase see ...
        #    ... https://docs.python.org/2/library/ ...
        #    ... itertools.html#itertools.izip
        iterList = [arrNodeNb,] + [iter(arrCoords)]*3
        self.nodeCoords.update( (int(x[0]), list(x[1:4]))
                                for x in izip(*iterList) )

    def writeElNodesOneShape(self, pickleFile, shape):
        """
        @param shape: element shape to be written to file. Example: "TET_Q"
        """
        elemNbs = sorted(self.shapeEl[shape])
        if not elemNbs:
            return
        nbNodes = len(self.elNodes[elemNbs[0]])
        msg("Writing mesh: shape %s (%d nodes) with %d elements"
            % (shape, nbNodes, len(elemNbs)), debugLevel=5)
        pickleFile.write("ELSHAPE %s\n" % shape)
        pickleFile.write("ELNODES %d %d\n" % (nbNodes, len(elemNbs)))
        a = array("I", (x for elem in elemNbs
                        for x in self.elNodes[elem]))
        a.tofile(pickleFile)
        pickleFile.write("ELEMENT_NUMBERS %d\n" % len(elemNbs))
        a = array("I", elemNbs)
        a.tofile(pickleFile)

    def readElNodesOneShape(self, pickleFile):
        """
        @Note: Elements from the pickle file are added to existing elements
        overwriting those with identical numbers.
        """
        line = pickleFile.readline().rstrip().upper().split(None)
        if line[0] != "ELSHAPE":
            raise BadPickleFile(
                "Expected ELSHAPE option. Instead found %s." % line)
        shape = line[1]
        msg("Reading mesh: reading elements of type %s" % shape,
            debugLevel=5)

        line = pickleFile.readline().rstrip().upper().split(None)
        if line[0] != "ELNODES":
            raise BadPickleFile(
                "Expected ELNODES option. Instead found %s." % line)
        nbNodes = int(line[1])
        nbElems = int(line[2])
        arrNodes = array("I")
        arrNodes.fromfile(pickleFile, nbNodes*nbElems)
        msg("Reading mesh: read connectivity of %d elements with"
            " %d nodes each" % (nbElems, nbNodes), debugLevel=5)

        line = pickleFile.readline().rstrip().upper().split(None)
        if line[0] != "ELEMENT_NUMBERS":
            raise BadPickleFile(
                "Expected ELEMENT_NUMBERS option. Instead found %s." % line)
        if nbElems!=int(line[1]):
            raise BadPickleFile(
                "Incompatible number of elements in option ELNODES (%d)"
                " and option ELEMENT_NUMBERS (%d)."
                % (nbElems, int(line[1])))
        arrElemNb = array("I")
        arrElemNb.fromfile(pickleFile, nbElems)
        msg("Reading mesh: read %d element numbers" % nbElems,
            debugLevel=5)

        # create elNodes dictionary
        #  - for the izip(*[iter]*N)-phrase see ...
        #    ... https://docs.python.org/2/library/ ...
        #    ... itertools.html#itertools.izip
        iterList = [arrElemNb,] + [iter(arrNodes)]*nbNodes
        elNodes = dict( (int(x[0]), [int(xx) for xx in x[1:nbNodes+1]])
                        for x in izip(*iterList) )
        self.updateElems(elNodes, shape)

    def write(self, pickleFile):
        """Write the mesh to the given open file object."""
        # write self.nodeCoords
        if not hasattr(pickleFile, "write"):
            raise ValueError(
                "Mesh.write() needs a file-like object as argument. Instead"
                " got %s: %s" % (type(pickleFile), pickleFile))
        self.writeNodeCoords(pickleFile)

        # write self.elnodes
        elShapes = sorted(
            shape for shape, elems in self.shapeEl.iteritems() if elems)
        pickleFile.write("UNSTRUCTURED_MESH %d shapes\n" % len(elShapes))
        for shape in elShapes:
            self.writeElNodesOneShape(pickleFile, shape)

    def read(self, pickleFile):
        """Read the mesh from the given open file object.
        @Note: Existing mesh data (nodes, elements) will be overwritten in
        case of identical labels.
        """
        # read nodeCoords
        self.readNodeCoords(pickleFile)

        # element shapes
        line = pickleFile.readline().rstrip().upper().split(None)
        if line[0] != "UNSTRUCTURED_MESH":
            raise BadPickleFile(
                "Expected UNSTRUCTURED_MESH option. Instead found %s." % line)
        nbElShapes = int(line[1])
        msg("Reading mesh: reading UNSTRUCTURED_MESH with %d shapes"
            % nbElShapes, debugLevel=5)
        for _ in range(nbElShapes):
            self.readElNodesOneShape(pickleFile)

    #} end of file I/O
    ## ------------------------------------------------------------------------

    ## ------------------------------------------------------------------------
    #{ consistency checks

    def checkElShape(self):
        "Check if set(elShape) == set(elNodes) and check shapeEl"
        if set(self.elShape) != set(self.elNodes):
            raise ModelInconsistent(
                "Different elements in Mesh.elNodes (#=%d) and Mesh.elShape"
                " (#=%d)." % (len(self.elNodes), len(self.elShape)))

        in_shapeEl = set()
        for typ, els in self.shapeEl.iteritems():
            if len(in_shapeEl.intersection(els))!=0:
                raise ModelInconsistent(
                    "%d element(s) found more than once in Mesh.shapeEl"
                    % len(in_shapeEl.intersection(els)))
            in_shapeEl.update(els)
            if not all([typ == self.elShape[el] for el in els]):
                wrongElTypeEls = [el for el in els
                                  if typ != self.elShape[el]]
                nbWrong = len(wrongElTypeEls)
                wrongElTypeEls.sort()
                wrongElTypeEls = wrongElTypeEls[:4]
                raise ModelInconsistent(
                    "Inconsistency between shapeEl and elShape: %d elements in"
                    " shapeEl['%s'] have other types: Some examples:\n%s"
                    % (nbWrong, typ,
                       "; ".join(["element %d: %s" % (el,self.elShape[el])
                                  for el in wrongElTypeEls])))

        if (in_shapeEl != set(self.elShape)):
            raise ModelInconsistent(
                "all(shapeEl.values()) != set(elShape)")
        return

    def checkElVolPos(self, elems=None, tol=1E-300, raiseExc=True):
        """Check that all elements have a positive non negative volume

        @param elems: list or iterable of element numbers to check. Defaults
           to all elements of the model. (See note below.)

        @param tol: the volume must not be smaller than tol. Specify 0.0
           if you want to accept elements with zero volume (e.g. Cohesives).

        @param raiseExc: If true raise ZeroOrNegVolumeElements exception if
           bad elements are being found.

        @raises ZeroOrNegVolumeElements: ... see raiseExc argument.

        @Note: Only tet elements considered so far.
        """

        badElems = list()

        # elemsAll = self.shapeEl["TET"]
        elemsAll = set.union(*(self.shapeEl.get(x, set())
                               for x in ("TET_L", "TET_Q", "TETHT_L")))
        if elems is not None:
            elemsAll.intersection_update(elems)

        ticker = MsgTicker(
            "Checking tet element %%d/%d for bad volume." % len(elemsAll))
        for element in elemsAll:

            ticker.tick()
            x = [self.nodeCoords[node]
                 for node in self.elNodes[element][:4]]
            # coordinates relative to the first node
            # dx = [[xxi-x0i for xxi, x0i in zip(xx, x[0])]
            #       for xx in x[1:]]
            dx = [[x[1][0]-x[0][0], x[1][1]-x[0][1], x[1][2]-x[0][2]],
                  [x[2][0]-x[0][0], x[2][1]-x[0][1], x[2][2]-x[0][2]],
                  [x[3][0]-x[0][0], x[3][1]-x[0][1], x[3][2]-x[0][2]]]
            # determinant of the "Jacobian"
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
            if detJ < tol:
                # negative or zero volume
                badElems.append(element)

        del ticker
        if raiseExc and len(badElems):
            raise ZeroOrNegVolumeElements(
                "Mesh.checkElVolPos() found %d bad elements" % len(badElems))

        return badElems
    #} end of consistency checks
    ## ------------------------------------------------------------------------

    ## ------------------------------------------------------------------------
    #{ some service functions: iterator function
    def elemFacesIter(self, elNum, sepLQ=False):
        """Returns an iterator of the faces of the given element

        >>> myElem = 5731
        >>> for faces in elemFacesIter(myElem):
        >>>    ...

        faces will be:
          - frozenset((3, 43, 12)), then
          - frozenset((43, 3, 36)) etc.

        @param elNum: element label
        @param sepLQ: "use seperate Linear/Quadratic element shape?"
           This option only makes a difference for quadratic elements.
           If True then return frozensets of all nodes on the respective face.
           If False then return frozensets of only the linear, i.e. corner-
           nodes on each face.

        @returns: An iterator that yields all faces of the element. Each face
           is represented by a frozenset of the node labels on this face. For
           a linear tet element this would be all combinations of three out of
           the four corner nodes of this tet.
        """
        elShape = self.elShape[elNum]
        if not sepLQ:
            elShape = elShape.replace("_Q", "_L")
        try:
            thisNodeIdx = self.elShapeFaceNodes[elShape]
        except KeyError:
            raise KeyError("elShape %s not equivalent to one of the known:"
                           " %s." % (elShape, self.elShapeFaceNodes.keys()))

        nodesOfThisElem = self.elNodes[elNum]
        for i in thisNodeIdx:
            yield frozenset([nodesOfThisElem[j] for j in i])
    #} end of some service functions: iterator function
    ## ------------------------------------------------------------------------

    ## ------------------------------------------------------------------------
    #{ edit / modify

    def updateElems(self, elNodes, elShape):
        """Inserts new elements into the mesh or replaces existing.

        @param elNodes: dictionary {element number: node number
        list} of elements to be updated or inserted
        @param elShape: Is the string of the element shape.
        """

        # possibly remove elements from self.shapeEl
        overwriteElems = set(elNodes).intersection(self.elNodes)
        for elem in overwriteElems:
            self.shapeEl[self.elShape[elem]].remove(elem)

        # update Mesh object attributes
        dict.update(self.elNodes, elNodes)
        self.elShape.update(dict.fromkeys(elNodes, elShape))
        try:
            elems = self.shapeEl[elShape]
        except KeyError:
            elems = set()
            self.shapeEl[elShape] = elems
        elems.update(elNodes)

    def removeElems(self, elems, removeNodes=False):
        """removes elements from the mesh
        updates elNodes, elShape, shapeEl dicts

        @param elems: A list or a set of element numbers (also accepts a dict,
        in that case takes its keys.)

        @param removeNodes: if True remove unused nodes. This is not very
        efficient if many elements are beeing deleted, in this case rather
        clean up after all elements have been removed using
        L{getConnectedNodes}:
          >>> mesh.removeElems(elems=myBigElset)
          >>> allNodes = mesh.getConnectedNodes()
          >>> removeNodes = set(mesh.nodeCoords).difference(allNodes)
          >>> for node in removeNodes:
          ...    del mesh.nodeCoords[node]

        @note: Silently ignores missing elements.
        @note: It is save to do the following (it is not only save but the
        recommended way to accomplish those tasks):
          >>> mesh.removeElems(elems=mesh.elNodes)  # remove all elements
          >>> mesh.removeElems(elems=mesh.shapeEl['TRI_L']) # remove one type
        """

        # special treatment if you want to remove all, wouldn't work otherwise
        if elems is self.elNodes:
            self.elNodes.clear()
            self.shapeEl.clear()
            self.elShape.clear()
            return

        # node list of element nodes
        if removeNodes:
            removeNodesSet = set()
            for elNum in elems:
                try:
                    removeNodesSet.update(self.elNodes[elNum])
                except KeyError:
                    # elNum not in self.elNodes
                    pass

        # remove elements from self.elNodes and self.elShape
        for elNum in elems:
            try:
                del self.elShape[elNum]
            except KeyError:
                pass

            try:
                del self.elNodes[elNum]
            except KeyError:
                pass

        # remove elements from self.shapeEl
        # use a copy (self.shapeEl.items()) in the loop
        for thisshape, els in self.shapeEl.items():
            els.difference_update(elems)
            if len(els)==0:
                del self.shapeEl[thisshape]

        # remove nodes
        if removeNodes:
            for nodes in self.elNodes.itervalues():
                removeNodesSet.difference_update(nodes)
                if len(removeNodesSet)==0:
                    break
            for node in removeNodesSet:
                del self.nodeCoords[node]

    def renumber(self, nodeStart=1, elemStart=1):
        """renumber elements and nodes

        Examples:
          >>> # renumber node and element numbers consecutively starting at 1
          >>> mesh.renumber()

          >>> # renumber only the nodes to start at 500001
          >>> mesh.renumber(nodeStart=500001, elemStart=None)

        @param nodeStart: New start for node numbers. If None, nodes are not
        renumbered.
        @param elemStart: New start for element numbers. If None, elements are
        not renumbered.

        @Note: The relative order of nodes and elements is not touched.

        @Note: If the element numbers are touched (elemStart!=None), the
        point search must eventually be initialized again. This could be changed
        (i.e. such that this method also updates the point search initialization
        but this is not done.

        @Returns: Tuple of two dictionaries (nodesOldToNew, elemsOldToNew)
        maping old node/element numbers to the new ones.

        @Raises KeyError: If there are nodes defined in the model but at least
        one node used in an element's definition is *not*. Or if the elShape
        or shapeEl dictionaries are not consistent, i.e. they contain an
        element that is not defined in elNodes.
        """

        if nodeStart is None or len(self.nodeCoords)==0:
            nodesOldToNew = dict()
        else:
            nodesOld = sorted(self.nodeCoords)
            nodesOldToNew = dict(izip(nodesOld, count(nodeStart)))

            newNodeCoords = dict(
                (nodesOldToNew[oldNode], coords)
                for oldNode, coords in self.nodeCoords.iteritems())
            self.nodeCoords = newNodeCoords

            try:
                newElNodes = dict(
                    (el, [nodesOldToNew[oldNode] for oldNode in nodes])
                    for el, nodes in self.elNodes.iteritems())
            except KeyError:
                for el, nodes in self.elNodes.iteritems():
                    notFound = sorted(
                        set(nodes).difference(nodesOldToNew))
                    if notFound:
                        break
                text = (
                    "Old node number(s) %s found in the connectivity list of"
                    " old element number %d is not defined in the model.\n"
                    "This would leave the model in a corrupted state. Bailing"
                    " out." % (notFound, notFound and el))
                # ... includes additional check for notFound not empty ...
                msg(text)
                raise KeyError(text)
            else:
                self.elNodes = newElNodes

        if elemStart is None or len(self.elNodes)==0:
            elemsOldToNew = dict()
        else:
            elemsOld = self.elNodes.keys()
            elemsOld.sort()
            elemsOldToNew = dict(izip(elemsOld, count(elemStart)))

            newElNodes = dict(
                (elemsOldToNew[oldEl], nodes)
                for oldEl, nodes in self.elNodes.iteritems())
            self.elNodes = newElNodes

            if isinstance(self.elShape, dict):
                try:
                    newElShape = dict(
                        (elemsOldToNew[oldEl], shape)
                        for oldEl, shape in self.elShape.iteritems())
                except KeyError:
                    notFound = sorted(
                        set(self.elShape).difference(elemsOldToNew))
                    text = (
                        "%d old element number(s) %d found in self.elShape"
                        " is/are not defined in the model. E.g.:\n%s\n"
                        "This would leave the model in a corrupted state."
                        " Bailing out." % (len(notFound), notFound[:20]))
                    msg(text)
                    raise KeyError(text)
                else:
                    self.elShape = newElShape

            # this is treated differently because it's a defaultdict
            if isinstance(self.shapeEl, dict):
                try:
                    newShapeEl = [
                        (shape, set(elemsOldToNew[oldEl] for oldEl in elems))
                        for shape, elems in self.shapeEl.iteritems()]
                except KeyError:
                    for shape, elems in self.shapeEl.iteritems():
                        notFound = sorted(
                            elems.difference(elemsOldToNew))
                        if notFound:
                            break
                    text = (
                        "%d old element number(s) found in self.shapeEl['%s']"
                        " is/are not defined in the model. Some of them:\n%s\n"
                        "This would leave the model in a corrupted state."
                        " Bailing out."
                        % (len(notFound), shape, notFound[:20]))
                    msg(text)
                    raise KeyError(text)
                else:
                    self.shapeEl.clear()
                    self.shapeEl.update(newShapeEl)

            # element numbers changed => invalidates point search initialization
            self.point_search_initialized = False

        return (nodesOldToNew, elemsOldToNew)

    #} end of edit / modify
    ## ------------------------------------------------------------------------
    ## ------------------------------------------------------------------------
    #{ query mesh properties, findAtMesh stuff

    def getBoundingBox(self):
        """Bounding Box of the whole mesh. Ignore nodes not connected to any
        element.

        Example:

        >>> # bounding box of that volume
        >>> bb = v.getBoundingBox()
        >>> print "lower left front corner:", bb[0]
        >>> print "upper right back corner:", bb[1]

        @Returns: A list of two coordinate tuples (rather lists)
        The first coord tuple states the min value, the last the max.
        I.e. self.getBoundingBox()[i][j]:
        i \in {0,1} for min/max value, j \in {0,1,2} coordinate index
        """
        try:
            dummy = self.boundingBox[0]
        except (AttributeError, IndexError):
            self.boundingBox = BoundingBox()
            self.boundingBox.update(
                nodeCoords=self.nodeCoords, elNodes=self.elNodes)

        return self.boundingBox

    def getElemVolField(self, elems=None, fieldName="elemVol",
                        wedgeThickness=None):
        """Returns a Field containing the element volumes, which is basically
        a dict {element number: element volume}

        @param elems: A sequence of element numbers or None for all elements.
           Elements not in mesh are silently ignored. Can also be a dict, then
           its keys are taken.
        @param fieldName: name of the field for output purposes
        @param wedgeThickness: If not None: specified thickness for cohesive
           elements.

        @Note: The volume of Tet elements is approximated by the volume of the
           tetrahedron with the four corner nodes and flat faces.
        """

        if wedgeThickness==0:
            wedgeThickness = None

        # we need bae.field_01.Field but can't import globally because we'd get
        # a cyclic import: bae.mesh_01 is needed in bae.field_01 as well
        from bae.field_01 import Field
        volFld = Field.classFromPosType(fieldName, "element", "scalar")()

        if elems is None:
            nbElems = len(self.elNodes)
        else:
            nbElems = len(elems)
        ticker = MsgTicker(
            "... calculated %%d/%d element volumes so far." % nbElems)

        cnt = 0
        unknownShapes = list()
        nbUnknownShape = 0
        for shape, elemsAll in self.shapeEl.iteritems():
            if elems is not None:
                # note: don't use intersection_update, copy elemsAll!
                # (otherwise you would modify self.shapeEl)
                elemsAll = elemsAll.intersection(elems)
            if len(elemsAll)==0:
                continue
            if shape[:-2] == "TET" or shape[:-2] == "TETHT":
                for elementNumber in elemsAll:
                    cnt += 1
                    try:
                        x = [self.nodeCoords[node]
                             for node in self.elNodes[elementNumber][:4]]
                    except KeyError:
                        msg('WARNING: Mesh.getElemVolField():'
                            ' For elementNumber=%d node %d is not in the mesh.'
                            % (elementNumber, node))
                        continue
                    # volume
                    vec = [vector_minus(pt, x[0])
                           for pt in x[1:]]
                    vol = ((
                          vec[0][0]*vec[1][1]*vec[2][2]
                        + vec[0][1]*vec[1][2]*vec[2][0]
                        + vec[0][2]*vec[1][0]*vec[2][1]
                        - vec[0][2]*vec[1][1]*vec[2][0]
                        - vec[0][1]*vec[1][0]*vec[2][2]
                        - vec[0][0]*vec[1][2]*vec[2][1])
                           / 6.0)
                    volFld[elementNumber] = vol
                    ticker.msg(cnt)
            elif shape[:-2] == "WEDGE":
                for elementNumber in elemsAll:
                    cnt += 1
                    try:
                        x = [self.nodeCoords[node]
                             for node in self.elNodes[elementNumber][:6]]
                    except KeyError:
                        msg('WARNING: Mesh.getElemVolField():'
                            ' For elementNumber=%d node %d is not in the mesh.'
                            % (elementNumber, node))
                        continue
                    # volume
                    area4 = vector_plus(
                        cross(vector(x[0],x[1]), vector(x[0],x[2])),
                        cross(vector(x[3],x[4]), vector(x[3],x[5]))
                        )
                    if wedgeThickness is None:
                        height3 = vector_minus(
                            vector_sum(x[3],x[4],x[5]),
                            vector_sum(x[0],x[1],x[2])
                        )
                        vol = dot( area4, height3 ) / 12.0
                    else:  # if wedgeThickness is not None:
                        vol = 0.25 * length(area4) * wedgeThickness
                    volFld[elementNumber] = vol
                    ticker.msg(cnt)
            else:
                unknownShapes.append(shape)
                nbUnknownShape += len(elemsAll)
                cnt += len(elemsAll)
        del ticker

        if nbUnknownShape>0:
            msg("WARNING: Can't compute element volume for %d elements of the"
                " following unrecognized shapes: %s. "
                % (nbUnknownShape, unknownShapes))

        return volFld

    def getElCentroids(self, elems=None):
        """
        Calculate element centroids for all specified elements and return a
        dict {element: centroid}

        @param elems: A sequence of element numbers or None for all elements.
        Elements not in mesh are silently ignored. Can also be a dict, then
        its keys are taken.

        @note: IMPORTANT! This is function returns the mean of the coordinates
        of the corner nodes of the elements. This is intended and exactly what
        you need for certain use cases. It is not the centre of gravity of the
        element in all cases, however!

        @note: The message ticker is wrong when not all elements in elems are
        actually present in the mesh, but besides that it should also work in
        that case.
        """
        if elems is None:
            nbElems = len(self.elNodes)
        else:
            nbElems = len(elems)
        ticker = MsgTicker(
            "... calculated %%d/%d element centroids so far." % nbElems)
        el_pts = dict()
        cnt = 0
        for shape, elemsAll in self.shapeEl.iteritems():
            if elems is not None:
                # note: don't use intersection_update, copy elemsAll!
                # (otherwise you would modify self.shapeEl)
                elemsAll = elemsAll.intersection(elems)
            if len(elemsAll)==0:
                continue
            try:
                avg_n = self.elShapeToNbCornerNodes[shape]
            except KeyError:
                raise UnknownElShape(
                    'Element shape %s not implemented in %s. Cannot calculate'
                    ' centroid.' % (shape, self.__module__))
            for el in elemsAll:
                cnt += 1
                try:
                    nodes = self.elNodes[el]
                except KeyError:
                    continue
                try:
                    x=y=z=0.0
                    for node in nodes[:avg_n]:
                        coords=self.nodeCoords[node]
                        x=x+coords[0]/avg_n
                        y=y+coords[1]/avg_n
                        z=z+coords[2]/avg_n
                except KeyError:
                    raise KeyError(
                        'Mesh.getElCentroids(elementNumber=%d): node %d'
                        ' is not in the mesh.' % (el, node))
                el_pts[el] = [x,y,z]
                ticker.msg(cnt)
        del ticker
        return el_pts

    def getElGaussPoints(self, elems=None):
        """
        Calculate the coordinates of the Gauss points for the specified elements
        and return a dict {element: [list of gauss points]}

        @param elems: A sequence of element numbers or None for all elements.
        Elements not in the mesh are silently ignored. Can also be a dict, then
        its keys are taken.

        @note: Only tested for tet elements so far. Only considers the linear
        shape functions. This is correct for plane element faces in the
        undeformed reference.

        @note: See also L{bae.abq_model_02.Model.getQuadGaussPointsDict}.
        """
        if elems is None:
            nbElems = len(self.elNodes)
        else:
            nbElems = len(elems)
        ticker = MsgTicker(
            "... calculated %%d/%d element Gauss points so far." % nbElems)
        el_pts = dict()
        cnt = 0
        for shape, elemsAll in self.shapeEl.iteritems():
            if elems is not None:
                # note: don't use intersection_update, copy elemsAll!
                # (otherwise you would modify self.shapeEl)
                elemsAll = elemsAll.intersection(elems)
            if len(elemsAll)==0:
                continue
            try:
                gpElCoords = self.elShapeGaussPoints[shape]
            except KeyError:
                raise UnknownElShape(
                    "Element shape %s not fully implemented in %s. Gauss"
                    " point coordinates are missing in elShapeGaussPoints."
                    " Cannot calculate Gauss points."
                    % (shape, self.__module__))
            try:
                elemShapeFunc = self.elShapeElemShapeFunc[shape]
            except KeyError:
                raise UnknownElShape(
                    "Element shape %s not fully implemented in %s. Element"
                    " shape functions are missing in elShapeElemShapeFunc."
                    " Cannot calculate Gauss points."
                    % (shape, self.__module__))
            if shape=="TET_Q":
                # special: for TET_Q only consider corner nodes
                # ... for efficiency
                gpShapeFuncts = gpElCoords
            else:
                # the general case
                gpShapeFuncts = [elemShapeFunc(c) for c in gpElCoords]
            for el in elemsAll:
                cnt += 1
                try:
                    nodes = self.elNodes[el]
                except KeyError:
                    continue
                # coords for all relevant nodes
                nodeCoords = [self.nodeCoords[node]
                              for node in nodes[:len(gpShapeFuncts[0])]]
                # sum up x*sf, y*sf, z*sf for all nodes and shape functions sf
                el_pts[el] = [map(sum, izip(*(
                                [x*sf for x in pt]
                                for pt, sf in izip(nodeCoords, shapeFuncts))))
                              for shapeFuncts in gpShapeFuncts]
                ticker.msg(cnt)

        del ticker
        return el_pts

    _cellSizeFactor=0.5      # cell Size for array relative to 'avgSize'
    def initializePointSearch(self, box=None, elems=None,
                              shapesConsidered=["TET_L", "TET_Q", "TETHT_L"]):
        """DEPRECATED: Use L{getElemCoordsForGrid} without separate
        initialization instead.

        Call this function once before point lookup with
        self.getElemCoords().

        If it's not explicitly called this is done automatically by
        self.getElemCoords() with the default arguments.

        @param box: All points to be supplied to self.getElemCoords() certainly
        lie within or on the border of this box.
        box = [[xmin, ymin, zmin], [xmax, ymax, zmax]]

        @param elems: A sequence of elements that are to be considered

        @param shapesConsidered: An iterable of element shapes to be
        considered.

        @note: If the self.getElemCoords() method is called for points that
        do not lie within the specified box but in the mesh (self), the
        result is unreliable, most probably self.pointIsInside() will state
        that this point is not in the mesh. Same if a point is in an element
        not contained in the elems argument.
        """

        assert (box is None) or all(x1<=x2 for x1, x2 in izip(*box))

        shapesConsidered = set(shapesConsidered).intersection(self.shapeEl)

        # dict of centroid-coord-tuples, should be initialized in any case
        # elCentroids is a local reference for speed
        elCentroids = dict()
        self.elCentroids = elCentroids

        if elems is None:
            nbElems = sum(len(self.shapeEl.get(shape, set()))
                          for shape in shapesConsidered)
        else:
            elems = set(elems).intersection(self.elNodes)
            nbElems = len(elems)

        msg("Initializing point search, element boxes...")

        # list of [[xmin,ymin,zmin],[xmax,ymax,zmax]] for each element
        elementBox=dict()

        nodeCoords = self.nodeCoords    # local copy for efficiency
        volBoundingBox = BoundingBox()  # later copied to self.boundingBox
        relevantElements = list()
        avgSize=0

        cnt = 0
        ticker = MsgTicker(
            "... processed %%d/%d tet elements" % nbElems, printLast=True)
        ticker.msg(cnt)  # at least one message as diagnostic output
        for shape in shapesConsidered:
            elemsAll = self.shapeEl.get(shape, set())
            try:
                avg_n = self.elShapeToNbCornerNodes[shape]
            except KeyError:
                raise UnknownElShape(
                    'Element shape %s not implemented in %s. Cannot calculate'
                    ' centroid.' % (shape, self.__module__))
            if elems is not None:
                # note: don't use intersection_update, copy elemsAll!
                # (otherwise you would modify self.shapeEl)
                elemsAll = elemsAll.intersection(elems)
            for el in elemsAll:
                cnt += 1
                try:
                    nodes = self.elNodes[el]
                except KeyError:
                    continue
                try:
                    # xyz = [[x1, x2, x3, x4], [y1, y2, y3, y4], [z1, ...]]
                    xyz = zip(*(nodeCoords[node]
                                for node in nodes[:avg_n]))
                except KeyError:
                    # not all nodes in nodeCoords
                    continue

                newElemBox = [[min(xyz[0]),min(xyz[1]),min(xyz[2])],
                              [max(xyz[0]),max(xyz[1]),max(xyz[2])]]

                # check if this elements bounding box intersects the box
                # argument
                disjoint = False
                if box is not None:
                    for e1,e2,b1,b2 in izip(newElemBox[0],newElemBox[1],
                                            box[0],box[1]):
                        if min(e2,b2)<max(e1,b1):
                            disjoint = True
                            break
                # ignore this element if it does not intersect box
                if not disjoint:

                    elCentroids[el] = [sum(xlist)/avg_n for xlist in xyz]

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

                    relevantElements.append(el)
                    avgSize += dist(newElemBox[1],newElemBox[0])
                    elementBox[el] = newElemBox

                # end if
                ticker.msg(cnt)
        del ticker

        if not len(relevantElements):
            # this is not sensible but it prevents raising errors
            if box is not None:
                avgSize = dist(*box)
                volBoundingBox.update(box)
            else:
                avgSize = 1.0
                volBoundingBox.update([[1,1,1]])
            if avgSize<1E-200:
                avgSize = 1.0
        else:
            avgSize /= len(relevantElements)

        # boundingBox = intersection(boundingBox, box)
        if box is not None:
            volBoundingBox[0] = map(max, volBoundingBox[0], box[0])  # min
            volBoundingBox[1] = map(min, volBoundingBox[1], box[1])  # max

        msg('... bounding box of the relevant volume: %s' % volBoundingBox)
        msg('... %d of the %d elements are within or close to the limiting'
            ' box.' % (len(relevantElements), nbElems))
        msg('... avg. element size: %g' % avgSize)

        ## -----------------------------------------------------------------
        ## map elements to cells

        msg("Initializing point search, initializing cells...")

        cellOrigin=volBoundingBox[0]
        cellSize = avgSize*self._cellSizeFactor

        # x,y,z number of cells
        cellDimensions = [volBoundingBox[1][i]-volBoundingBox[0][i]
                          for i in range(3)]
        msg('cell array dimensions: %g %g %g' % tuple(cellDimensions))

        cellNumber = [int(cdim/cellSize)+1 for cdim in cellDimensions]
        msg('cell numbers: %g %g %g' % tuple(cellNumber))

        # create empty array
        self.cells=[]
        ticker = MsgTicker("... initialized %%3.1f%%%% of %d cells"
                           % (cellNumber[0]*cellNumber[1]*cellNumber[2]))
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
        msg("Initializing point search, mapping elements to cells...")
        ticker = MsgTicker(
            "... processing tet element %%d of %d" % len(relevantElements))
        cellIndex=[[0,0,0],[0,0,0]]  # low,high
        for element in relevantElements:
            ticker.tick()
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

        # finished initializing
        self.point_search_initialized = True
        return

    def getElemCoords(self, points, nbPoints=None, tolerance=1E-4):
        """DEPRECATED: Use L{getElemCoordsForGrid} instead.

        Return a list of tuples (element number, element coordinates) for
        the given list of points.

        If a point is not within any of the elements of the mesh or more
        precisely in any of the elements considered by the last call to
        self.initializePointSearch() then (None, []) is returned.

        Elements implemented so far:

        TET_L, TET_Q, TETHT_L: The element coordinates are identical with the
        (linear tetrahedron) shape function associated with the four nodes
        SF1...SF4.  Each of the SF1...SF4 values lies between 0 and 1. SF1==1
        at the first node, SF2==1 at the second node and so on. All of those
        must be between 0 and 1 for a point to be inside the tet.

        If the point is not in any of the elements of an already implemented
        type it is considered to be outside the mesh.

        @param points: A list of point coordinates.

        @param nbPoints: Optional number of points to process. For diagnostic
        output if the points argument is an iterator.

        @param tolerance: assume points to be inside a particular element if
        none of the four element coordinates (usually valued between 0 and 1)
        have a negative value smaller than (-tolerance). Spurious tests
        seem to indicate that a value of 0.0 still finds all points on element
        borders.

        @Returns: list of tuples (element number, element coordinates) for
        the given list of points (for tets we have four element coords)

        """

        if not self.point_search_initialized:
            self.initializePointSearch()

        try:
            nbPoints = len(points)
        except TypeError:
            pass

        if len(self.elNodes)==0 or not self.point_search_initialized:
            if not nbPoints:
                nbPoints = sum(1 for _ in points)
            return [(None, []),]*nbPoints

        if nbPoints is None:
            ticker = MsgTicker("Calculating element coordinates nb %d")
        else:
            ticker = MsgTicker(
                "Calculating element coordinates %%d/%d" % nbPoints)
        elCoords = list()
        for point in points:
            ticker.tick()

            # default: not found
            thisElCoords = (None, [])

            # cell check (includes bounding box test from cell index)
            index=[0,0,0]
            for i in xrange(3):
                index[i]=int((point[i]-self.boundingBox[0][i])/self.cellSize)
            if any(index[i]<0 or index[i]>=self.cellNumber[i]
                   for i in xrange(3)):
                elCoords.append(thisElCoords)
                continue

            # if cell check passed ... go on
            elems_in_cell=self.cells[index[0]][index[1]][index[2]]
            centroid_to_point_distance=[]
            # sort by distance
            for element in elems_in_cell:
                l=dist(point,self.elCentroids[element])
                centroid_to_point_distance.append([l,element])
            centroid_to_point_distance.sort()

            # full check if point inside any tet of elems_in_cell
            for [l,element] in centroid_to_point_distance:

                elShape = self.elShape[element]

                # Tet element
                if elShape=="TET_L" or elShape=="TET_Q" or elShape=="TETHT_L":
                    # coordinates of the four corner nodes
                    x = [self.nodeCoords[node]
                         for node in self.elNodes[element][:4]]
                    # coordinates relative to the first node
                    # dp = [pi - x0i for pi, x0i in zip(point, x[0])]
                    # dx = [[xxi-x0i for xxi, x0i in zip(xx, x[0])]
                    #       for xx in x[1:]]
                    dp = [point[0]-x[0][0], point[1]-x[0][1], point[2]-x[0][2]]
                    dx = [[x[1][0]-x[0][0], x[1][1]-x[0][1], x[1][2]-x[0][2]],
                          [x[2][0]-x[0][0], x[2][1]-x[0][1], x[2][2]-x[0][2]],
                          [x[3][0]-x[0][0], x[3][1]-x[0][1], x[3][2]-x[0][2]]]
                    # determinant of the "Jacobian"
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
                        thisElCoords = (None, [])
                        continue

                    # inverse * determinant
                    # invJ = [ [a[ip2[1]]*b[ip2[2]] - a[ip2[2]]*b[ip2[1]]
                    #           for ip2 in cylcePos]
                    #          for a, b in ((dx[ip1[1]], dx[ip1[2]])
                    #                       for ip1 in cylcePos)]
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
                    if all(xi>=-tolerance for xi in shapeFuncs):
                        thisElCoords = (element, shapeFuncs)
                        break
                # other shapes not implemented yet, ignore those elements
                else:
                    thisElCoords = (None, [])
                    continue

            elCoords.append(thisElCoords)

        del ticker
        return elCoords

    def getClosestElemCoords(self, points, nbPoints=None):
        """Return a list of tuples (element number, element coordinates) for
        the given list of points.

        For any point in the points-list this function searches for the closest
        (tet-) element in terms of the distance to its centroid. Then it
        calculates the element coordinates of the point with respect to this
        closest element.

        This is similar to getElemCoords() but significantly different!

        The point does not need to lie inside any element, so the element
        coordinates need not make sense in the usual sense of showing the
        position within a particular element. At the same time this particular
        point might well be within another (larger) element whose centroid
        happens to be further away from the point. Suppose for example the
        situation sketched below. "1" and "2" indicate the centroids of the
        two tets (sketched as triangles) and "X" is the point to find the
        closest element for. This function would find element 2 and the
        element coordiantes would indicate that the point is outside of the
        element. Whereas it's obvious that the point X lies within the element
        nr 1::

            +----------------------------------+----+
                ''---__                        |   /
                       ''---__         1      X|2 /
                              ''---__          | /
                                     ''---__   |/
                                            ''-+

        Elements implemented so far:

        TET_L, TET_Q, TETHT_L: The element coordinates are identical with the
        (linear tetrahedron) shape function associated with the four nodes
        SF1...SF4.  Each of the SF1...SF4 values lies between 0 and 1. SF1==1
        at the first node, SF2==1 at the second node and so on. All of those
        must be between 0 and 1 for a point to be inside the tet.

        @param points: A list of point coordinates.

        @param nbPoints: Optional number of points to process. For diagnostic
        output if the points argument is an iterator.

        @Returns: list of tuples (element number, element coordinates) for
        the given list of points

        """

        try:
            nbPoints = len(points)
        except TypeError:
            pass

        if len(self.elNodes)==0:
            # empty mesh!
            if not nbPoints:
                nbPoints = sum(1 for _ in points)
            return [(None, []),]*nbPoints

        if not hasattr(self,"centroidsKDTree"):
            try:
                # alltets = self.shapeEl["TET"]
                alltets = set.union(*(self.shapeEl.get(x, set())
                                      for x in ("TET_L", "TET_Q", "TETHT_L")))

            except KeyError:
                raise KeyError('No Element of shape TET_L or TET_Q or TETHT_L')

            # checking element volumes
            badElems = self.checkElVolPos(elems=alltets, raiseExc=False)
            if len(badElems):
                msg("WARNING: Ignoring %d tet elements because of zero or"
                    " small volumes." % len(badElems))
                alltets = alltets.difference(badElems)

            msg("Calculating centroids for %d tet elements." % len(alltets))
            tetCentroids = self.getElCentroids(alltets)

            msg("Initializing kd-tree.")
            self.centroidsKDTree = KDTree(tetCentroids)

        if nbPoints:
            ticker = MsgTicker(
                "Calculating nearest element coordinates %%d/%d" % nbPoints)
        else:
            ticker = MsgTicker("Calculating nearest element coordinates nb %d")

        elCoords = list()
        for point in points:
            ticker.tick()

            # find closest element
            sqd, element = self.centroidsKDTree.knnSearch(point, K=1)[0]

            # coordinates of the four corner nodes
            x = [self.nodeCoords[node]
                 for node in self.elNodes[element][:4]]
            # coordinates relative to the first node
            # dp = [pi - x0i for pi, x0i in zip(point, x[0])]
            # dx = [[xxi-x0i for xxi, x0i in zip(xx, x[0])]
            #       for xx in x[1:]]
            dp = [point[0]-x[0][0], point[1]-x[0][1], point[2]-x[0][2]]
            dx = [[x[1][0]-x[0][0], x[1][1]-x[0][1], x[1][2]-x[0][2]],
                  [x[2][0]-x[0][0], x[2][1]-x[0][1], x[2][2]-x[0][2]],
                  [x[3][0]-x[0][0], x[3][1]-x[0][1], x[3][2]-x[0][2]]]
            # determinant of the "Jacobian"
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
            # negative or zero volume, this should not happen:
            # see earlier call to self.checkElVolPos()
            assert detJ>0.0

            # inverse * determinant
            # invJ = [ [a[ip2[1]]*b[ip2[2]] - a[ip2[2]]*b[ip2[1]]
            #           for ip2 in cylcePos]
            #          for a, b in ((dx[ip1[1]], dx[ip1[2]])
            #                       for ip1 in cylcePos)]
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
            thisElCoords = (element, shapeFuncs)
            elCoords.append(thisElCoords)

        del ticker
        return elCoords

    def getElemCoordsForGrid(self, grid,
                             elems=None,
                             shapesConsidered=["TET_L", "TET_Q", "TETHT_L"],
                             tolerance=1E-6):
        """Return a list of tuples (element number, element coordinates) for
        the given points grid.

        If a point is not within any of the elements of the mesh then
        (None, []) is returned.

        Elements implemented so far:

        TET_L, TET_Q, TETHT_L: The element coordinates are identical with the
        (linear tetrahedron) shape function associated with the four nodes
        SF1...SF4.  Each of the SF1...SF4 values lies between 0 and 1. SF1==1
        at the first node, SF2==1 at the second node and so on. All of those
        must be between 0 and 1 for a point to be inside the tet.

        If the point is not in any of the elements of an already implemented
        type it is considered to be outside the mesh.

        @param grid: L{MeshStructuredPoints}, L{MeshStructuredPointsRot}
        or L{MeshUnstructuredPoints} - object

        @param elems: A sequence of elements that are to be considered

        @param shapesConsidered: An iterable of element shapes to be
        considered.

        @param tolerance: assume points to be inside a particular element if
        none of the four element coordinates (usually valued between 0 and 1)
        have a negative value smaller than (-tolerance). Spurious tests
        seem to indicate that a value of 0.0 still finds all points on element
        borders. But for safety in favour of spending time on further
        investigating the matter the default is a small number.

        @Returns: list of tuples (element number, element coordinates) for
        the given list of points (for tets we have four element coords)

        """
        shapesConsidered = set(shapesConsidered).intersection(self.shapeEl)

        if elems is None:
            nbElems = sum(len(self.shapeEl.get(shape, set()))
                          for shape in shapesConsidered)
        else:
            elems = set(elems).intersection(self.elNodes)
            nbElems = len(elems)
        msg("Initializing point search."
            " %d elements in total (nb of element shapes: %d) will be"
            " considered for the point search of %d grid points."
            % (nbElems, len(shapesConsidered), len(grid)))

        if isinstance(grid, MeshStructuredPointsRot):
            # special for rotated box: use rotated mesh coords
            msg("Rotating coordinates of %d nodes of the mesh into the"
                " rotated grids coordinate system..." % len(self.nodeCoords))
            orig = grid.origin
            rotMat = mat_transpose(grid.rotmat)
            nodeCoords = dict(
                (n, vector_plus(orig, mat_multvec(
                    rotMat, vector(orig, coords))))
                for n, coords in self.nodeCoords.iteritems())
            # del orig, rotMat  # can't delete because of nested scope blabla
            msg("... finished.")

            # take some methods and the box from MeshStructuredPoints
            getIdsInBox = MeshStructuredPoints.getIdsInBox
            getPoint = MeshStructuredPoints.getPoint
            box = MeshStructuredPoints.getBoundingBox(grid)
            msg("For the point search will only consider elements in or"
                " intersecting with the bounding box of the rotated grid,"
                " ***taking the rotated coordinates of the grid points*** %s."
                % box, debugLevel=1)
        else:
            # mesh node coords, some grid methods: just the plain ordinary
            # items in this case. (For what? Compare to previous section)
            nodeCoords = self.nodeCoords
            getIdsInBox = type(grid).getIdsInBox
            getPoint = type(grid).getPoint
            box = grid.getBoundingBox()
            msg("For the point search will only consider elements in or"
                " intersecting with the bounding box of the grid points: %s."
                % box, debugLevel=1)

        # initialize resulting list of (elemNr, elemCoords)-tuples with
        # None for each grid point
        elCoords = [None]*len(grid)

        msg("Calculating element coordinates.")
        cnt = 0
        ticker = MsgTicker(
            "... processed %%d/%d tet elements" % nbElems, printLast=True)
        ticker.msg(cnt)  # at least one message as diagnostic output
        for shape in shapesConsidered:
            elemsAll = self.shapeEl.get(shape, set())
            try:
                avg_n = self.elShapeToNbCornerNodes[shape]
            except KeyError:
                raise UnknownElShape(
                    'Element shape %s not implemented in %s. Cannot calculate'
                    ' centroid.' % (shape, self.__module__))
            if elems is not None:
                # note: don't use intersection_update, copy elemsAll!
                # (otherwise you would modify self.shapeEl)
                elemsAll = elemsAll.intersection(elems)
            for elemNr in elemsAll:
                cnt += 1
                try:
                    nodes = self.elNodes[elemNr]
                except KeyError:
                    continue
                try:
                    # xyz = [[x1, x2, x3, x4], [y1, y2, y3, y4], [z1, ...]]
                    xyz = zip(*(nodeCoords[node]
                                for node in nodes[:avg_n]))
                except KeyError:
                    # not all nodes in nodeCoords
                    continue

                newElemBox = [[min(xyz[0]),min(xyz[1]),min(xyz[2])],
                              [max(xyz[0]),max(xyz[1]),max(xyz[2])]]

                # check if this elements bounding box intersects the box
                # argument
                disjoint = False
                for e1,e2,b1,b2 in izip(newElemBox[0],newElemBox[1],
                                        box[0],box[1]):
                    if min(e2,b2)<max(e1,b1):
                        disjoint = True
                        # msg("Elem #%d not considered, outside bb" % elemNr,
                        #     debugLevel=10)   ### DEBUG
                        break
                # ignore this element if it does not intersect box
                if disjoint:
                    continue

                # grid point ids that are inside the element's bounding box
                ptIds = getIdsInBox(grid, newElemBox)
                # msg("Elem #%d: %d pts in bb" % (elemNr, len(ptIds)),
                #     debugLevel=10)   ### DEBUG

                # ignore grid point if it's already got a
                # (elemNr, elemCoords)-tuple
                ptIds = [pt for pt in ptIds if not elCoords[pt]]
                if not ptIds:
                    continue

                # msg("Elem #%d: %d not yet assigned pts for exact check"
                #     % (elemNr, len(ptIds)), debugLevel=10)   ### DEBUG
                if shape=="TET_L" or shape=="TET_Q" or shape=="TETHT_L":
                    # preparation for calculation of elem coords
                    # for Tet element:

                    # coordinates relative to the first node
                    # (This is the Jacobian matrix dx_j / d xi_i assuming only
                    # linear shape functions.)
                    # dx = [[xxi-x0i for xxi, x0i in zip(xx, x[0])]
                    #       for xx in x[1:]]
                    dx = [[xyz[0][1]-xyz[0][0], xyz[1][1]-xyz[1][0],
                           xyz[2][1]-xyz[2][0]],
                          [xyz[0][2]-xyz[0][0], xyz[1][2]-xyz[1][0],
                           xyz[2][2]-xyz[2][0]],
                          [xyz[0][3]-xyz[0][0], xyz[1][3]-xyz[1][0],
                           xyz[2][3]-xyz[2][0]]]
                    # determinant of the Jacobian
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

                    # inverse * determinant of the Jacobian
                    # invJ = [ [a[ip2[1]]*b[ip2[2]] - a[ip2[2]]*b[ip2[1]]
                    #           for ip2 in cylcePos]
                    #          for a, b in ((dx[ip1[1]], dx[ip1[2]])
                    #                       for ip1 in cylcePos)]
                    invJ = [[dx[1][1]*dx[2][2] - dx[1][2]*dx[2][1],
                             dx[1][2]*dx[2][0] - dx[1][0]*dx[2][2],
                             dx[1][0]*dx[2][1] - dx[1][1]*dx[2][0]],
                            [dx[2][1]*dx[0][2] - dx[2][2]*dx[0][1],
                             dx[2][2]*dx[0][0] - dx[2][0]*dx[0][2],
                             dx[2][0]*dx[0][1] - dx[2][1]*dx[0][0]],
                            [dx[0][1]*dx[1][2] - dx[0][2]*dx[1][1],
                             dx[0][2]*dx[1][0] - dx[0][0]*dx[1][2],
                             dx[0][0]*dx[1][1] - dx[0][1]*dx[1][0]]]

                    for ptId in ptIds:
                        point = getPoint(grid, ptId)

                        # coordinates relative to the first node
                        # dp = [pi - x0i for pi, x0i in zip(point, x[0])]
                        dp = [point[0]-xyz[0][0],
                              point[1]-xyz[1][0],
                              point[2]-xyz[2][0]]

                        # parametric coordinates: xi_p := invJ * dp
                        # xip = [float(sum(invJ[i][j]*dp[j]
                        #                  for j in range(3)))/detJ
                        #        for i in range(3)]
                        xip = [float(invJ[0][0]*dp[0]+invJ[0][1]*dp[1]
                                     +invJ[0][2]*dp[2])/detJ,
                               float(invJ[1][0]*dp[0]+invJ[1][1]*dp[1]
                                     +invJ[1][2]*dp[2])/detJ,
                               float(invJ[2][0]*dp[0]+invJ[2][1]*dp[1]
                                     +invJ[2][2]*dp[2])/detJ]
                        shapeFuncs = [1-sum(xip), xip[0], xip[1], xip[2]]

                        # if inside: store elem number, elem coords
                        if all(xi>=-tolerance for xi in shapeFuncs):
                            elCoords[ptId] = (elemNr, shapeFuncs)

                # other shapes not implemented yet, ignore those elements
                else:  # if shape not in ("TET_L", "TET_Q", "TETHT_L"):
                    msg("WARNING: So far ony tet elements implemented."
                        " Ignoring elements of shape %s." % shape)
                    continue

                # end if
                ticker.msg(cnt)
        del ticker

        # convert values for not-found-points to (None, [])
        elCoords = [x or (None, []) for x in elCoords]
        return elCoords

    def splitElset_DisjointParts(self, elems, connectedBy="sharedFace"):
        """Returns a list of disjoint element-sets, each set containing
        connected elements.

        @param elems: A set of element numbers

        @param connectedBy: One of [ "sharedFace", "sharedEdge", "sharedNode" ]
          describing how elements are considered to be connected.
        """
        # assert(connectedBy in ["sharedFace", "sharedEdge", "sharedNode"])
        if connectedBy not in ["sharedFace", "sharedEdge", "sharedNode"]:
            raise ValueError(
                "Wrong optional connectedBy argument '%s' in method"
                " splitElset_DisjointParts. Must be one of sharedFace,"
                " sharedEdge or sharedNode." % connectedBy)

        #-- initialize return variable
        disjointElsets = []
        #spltElsets = defaultdict(set)

        #-- OPTIONAL to place here: loop over all elsets (here, just one)

        #-- build sharedPiece dict { sharedKey: set(elNumbers), ...}
        sharedDict = defaultdict(set)
        for el in elems:

            try:
                elShape = self.elShape[el]
            except KeyError:
                continue
            if not elShape.startswith("TET"):
                raise NotImplementedError(
                    "splitElset_DisjointParts: elShape '%s' not implemented."
                    " Only elShape 'TET' implemented so far." % elShape)
            else:
                elShape = 'TET'

            if connectedBy=="sharedFace":
                for faceElNds in self.elShapeFaceNodes[elShape]:
                    faceKey = frozenset([self.elNodes[el][idx]
                                         for idx in faceElNds])
                    sharedDict[faceKey].add(el)
            elif connectedBy=="sharedEdge":
                for edgeElNds in self.quadMidNodeEdges[elShape]:
                    edgeKey = frozenset([self.elNodes[el][idx]
                                         for idx in edgeElNds])
                    sharedDict[edgeKey].add(el)
            elif connectedBy=="sharedNode":
                for nodeKey in self.elNodes[el]:
                    sharedDict[nodeKey].add(el)

        #-- now split up
        testElems = set(elems)
        disjointPart = 0
        while len(testElems)>0:
            el = list(testElems)[0]
            lstElset = set()
            curElset = set([el,])
            incElset = set(curElset-lstElset)
            disjointPart += 1
            while len(incElset)>0:
                lstElset = set(curElset)
                for el in incElset:
                    if connectedBy=="sharedFace":
                        for faceElNds in self.elShapeFaceNodes[elShape]:
                            faceKey = frozenset([self.elNodes[el][idx]
                                                 for idx in faceElNds])
                            try:
                                connectedElems = set(sharedDict[faceKey])
                            except:
                                pass
                            else:
                                curElset.update(connectedElems)
                    elif connectedBy=="sharedEdge":
                        for edgeElNds in self.quadMidNodeEdges[elShape]:
                            edgeKey = frozenset([self.elNodes[el][idx]
                                                 for idx in edgeElNds])
                            try:
                                connectedElems = set(sharedDict[edgeKey])
                            except:
                                pass
                            else:
                                curElset.update(connectedElems)
                    elif connectedBy=="sharedNode":
                        for nodeKey in self.elNodes[el]:
                            try:
                                connectedElems = set(sharedDict[nodeKey])
                            except:
                                pass
                            else:
                                curElset.update(connectedElems)
                incElset = set(curElset-lstElset)
            #newName = eln+"P%d"%disjointPart
            #spltElsets[newName] = set(curElset)
            disjointElsets.append(set(curElset))

            testElems.difference_update(curElset)

        del sharedDict

        #-- OPTIONAL to place here: end of loop over all elsets (here, just one)

        return disjointElsets

    #} end of query mesh properties, findAtMesh stuff

    ## ------------------------------------------------------------------------
    ## ------------------------------------------------------------------------
    #{ other data extraction methods
    def getConnectedNodes(self, elems=None):
        """Returns a set of nodes that are connected to the elements given as
        argument.

        @param elems: A list or a set of element numbers (also accepts a dict,
          in that case takes its keys.) If elems==None, return all nodes
          connected to any element in the model. Elements not defined are
          silently ignored.

        @Note: it does not check if the nodes are actually defined in the model.
        """
        if elems is None:
            elemIter = self.elNodes.itervalues()
        else:
            elemIter = (self.elNodes[el]
                        for el in set(elems).intersection(self.elNodes))

        connectedNodes = set()
        for nodes in elemIter:
            connectedNodes.update(nodes)
        return connectedNodes

    def getExteriorSurface(self, elems=None, pruneNodeCoords=False):
        """Return a L{Mesh} object with all outer faces of this mesh as shell
        elements. Outer faces are those with exactly one element connected to.

        IMPORTANT NOTE: The nodeCoords attribute of the resulting surface is
        taken as a reference of self.nodeCoords. That means that changes of
        the position (like scaling, translation, ...) of either self or the
        resulting mesh object will affect both instances. But: see argument
        pruneNodeCoords below.

        In order to get an L{ElemFaceSurface<bae.surface_03.ElemFaceSurface>}
        use L{bae.surface_03.ElemFaceSurface.addFreeSurfFromElset}()
        or L{bae.surface_03.ElemFaceSurface.addFreeSurfFromSplitMesh}():
         >>> from bae.surface_03 import ElemFaceSurface
         >>> extSurf = ElemFaceSurface().addFreeSurfFromElset(mesh)
         >>> extSurf2 = ElemFaceSurface().addFreeSurfFromSplitMesh(meshQuad)

        @param elems: A list or a set of element numbers (also accepts a dict,
          in that case takes its keys.) If elems==None, return all nodes
          connected to any element in the model. Elements not defined are
          silently ignored.

        @param pruneNodeCoords: If True then the result gets a shallow copy
          of self.nodeCoords with only the nodes required for the mesh.
          Otherwise self.nodeCoords is used for the results as well.
        """
        if elems is None:
            nbElems = len(self.elNodes)
        else:
            nbElems = len(elems)

        # volume elements for which this method works.
        surfaceShapes = {
            "TET_L": "TRI_L",
            "TETHT_L": "TRI_L",
            "TET_Q": "TRI_Q"}

        shapeFaces = list()
        faceToNbTets = defaultdict(int)
        ticker = MsgTicker("Computing exterior surface: faces of"
                           " tet elements %%d/%d." % nbElems)
        cnt = 0
        for shape, elemsInShape in self.shapeEl.iteritems():
            if elems:
                elemsInShape = elemsInShape.intersection(elems)
            if not elemsInShape:
                continue
            try:
                faceShape = surfaceShapes[shape]
            except KeyError:
                msg("WARNING: %d elements of shape %s are not considered for"
                    " the exterior surface. Feature not implemented."
                    % (len(elemsInShape), shape))
                continue

            # process elements of current shape
            # reverse order of nodes for outward pointing normal
            orientedFaceIds = self.elShapeFaceNodes[shape][::-1]
            faceOrder = dict()
            shapeFaces.append((faceShape, faceOrder))
            for el in elemsInShape:
                ticker.msg(cnt+1)
                cnt += 1
                nodesOfThisTet = self.elNodes[el]
                for faceNodeIds in orientedFaceIds:
                    face = [nodesOfThisTet[i] for i in faceNodeIds]
                    id_ = frozenset(face)
                    faceOrder[id_] = face  # ordered face
                    faceToNbTets[id_] += 1
        del ticker

        ticker = MsgTicker()
        extFaces = set((id_ for id_, nb in faceToNbTets.iteritems() if nb==1))
        ticker.msg("Computing exterior surface: found %d exterior faces."
                   % len(extFaces))
        del ticker, faceToNbTets

        # create the new mesh object
        mesh = Mesh()

        # now add the actual faces to the surf object
        ticker = MsgTicker()
        cnt = 0
        for elShape, faceOrder in shapeFaces:
            faces = extFaces.intersection(faceOrder)
            elNodes = dict(
                (cnt+i+1, list(faceOrder[id_]))
                for i, id_ in enumerate(faces))
            mesh.updateElems(elNodes, elShape)
            cnt += len(faces)
            ticker.msg("Computing exterior surface: generated %d face elements."
                       % cnt)
        del ticker, shapeFaces, elShape, elNodes, faces

        # nodeCoords
        if pruneNodeCoords:
            allNodes = mesh.getConnectedNodes()
            mesh.nodeCoords = dict(
                (node, self.nodeCoords[node]) for node in allNodes)
        else:
            mesh.nodeCoords = self.nodeCoords

        return mesh

    def getTopSurface(self, surfZTol=0.01):
        """Get the top surface of self as
        L{ElemFaceSurface<bae.surface_03.ElemFaceSurface>}.

        This method assumes a quadratic tet mesh with cohesives. For a linear
        mesh or a tet-only mesh it will work as well but most of its work would
        be unnecessary.

        The resulting ElemFaceSurface keeps a linearized variant of self as
        mesh-attribute. I.e. all tet elements are there, not only those listed
        in the faceEl attribute. This makes it possible to get a time-sequence
        of top surfaces by means of the
        L{updateFromSequence.<bae.surface_03.ElemFaceSurface.updateFromSequence>}
        method of the returned ElemFaceSurface object:
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

        @param surfZTol: minimum z-component of the unit normal vector a.k.a.
           z-direction cosinus for a face to be considered "looking up" and
           therefore belonging to the top surface. The default 0.01 corresponds
           to all faces having at least an angle of 0.5 deg to a vertical line.
        """

        # can't import this upfront because it would cause a cyclic import
        from bae.utils.splitting_01 import rejoinSplitMesh
        from bae.surface_03 import ElemFaceSurface

        # create linear mesh without splitting
        ticker = MsgTicker("Identifying top surface: %s")
        ticker.msg("Removing cohesive elements from mesh, stitching the gaps.")
        linMesh = Mesh()
        linMesh.nodeCoords = self.nodeCoords
        linMesh.updateElems(
            dict( (el, self.elNodes[el][:4])
                  # for el in self.shapeEl["TET"]
                  for el in set.union(*(self.shapeEl.get(x, set())
                                        for x in ("TET_L", "TET_Q", "TETHT_L")))

            ), "TET_L")
        linMesh.updateElems(dict(
            (el, self.elNodes[el]) for el in self.shapeEl["WEDGE"]), "WEDGE_L")
        rejoinSplitMesh(linMesh)
        ticker.msg("Finished stitching cohesive gaps.")
        msg("Finished stitching cohesive gaps. Resulting mesh: %s" % linMesh,
            debugLevel=1)
        linMesh.checkElShape()

        # determine top surface of self
        ticker.msg("Creating exterior surface of whole mesh.")
        topSurf = ElemFaceSurface().addFreeSurfFromElset(linMesh)
        ticker.msg("Determining face normals for %d element faces."
                   % topSurf.countFaces())
        faceIdElemNormal = topSurf.getTetFaceNormal()
        ticker.msg("Searching for top faces by checking the %d face normals."
                   % sum(len(x) for x in faceIdElemNormal.itervalues()))
        for faceId, elemNormal in faceIdElemNormal.iteritems():
            upFaceElems = [el for el, normal in elemNormal.iteritems()
                           if normal[2]>surfZTol]
            topSurfElems = topSurf.faceEl[faceId]
            topSurfElems.intersection_update(upFaceElems)

            # remove empty set from the dict if applicable
            if not topSurfElems:
                del topSurf.faceEl[faceId]

        # finalize
        topSurf.mesh = linMesh
        ticker.msg("Finished.")
        msg("Created top surface of %d faces." % topSurf.countFaces(),
            debugLevel=1)
        del ticker
        return topSurf

    #} end of other data extraction methods

    ## ------------------------------------------------------------------------
    ## ------------------------------------------------------------------------
    #{ interpolation

    def getElemShapeFuncs(self, elemCoords, nbPoints=None):
        """The function returns a tuple (element number, [SF1, SF2, ..., SFN])
        for each item in elemCoords; N is the number of nodes on the particular
        element type (elShape). SF1...SFN are the values of the shape function
        associated with the (e.g.) ten nodes of the tet at the given point.

        If any of the points is not in any element of the mesh which would be
        indicated by an element number None in elemCoords then (None, []) is
        returned for that particular point.

        This can be used for interpolating nodal values. The following is for
        illustration purposes only! Use the Field.interpolateToPoints()
        methods in the bae.field_01 module for this.

        Suppose we have a mesh and a scalar field defined through values at the
        nodes:
         >>> mesh = ...                                       # a Mesh object
         >>> nodeToValue = {1:5.4, 2:2.1, 3:6.4, 4:3.8, ...}  # nodal values

        Now we want to interpolate the values to some arbitrary points:
         >>> gridPoints = [[1.0, 0.3, 0.6], [4.0, 1.6, 2.8]]  # some points
         >>> elemCoords = mesh.getElemCoords(gridPoints)
         >>> elemSFuncs = mesh.getElemShapeFuncs(elemCoords)
         >>> pointValues = [
         ...     sum(weight*nodeToValue[node]
         ...         for node, weight in zip(mesh.elNodes[elem], weights))
         ...     for elem, weights in elemSFuncs]

        Write point coordinates and corresponding values to a csv file:
         >>> import csv
         >>> output = csv.writer(open("myoutput.csv", "wb"))
         >>> output.writerow(['x', 'y', 'z', 'value'])
         >>> for pt, val in zip(gridPoints, pointValues):
         >>>     output.writerow(pt+[val])  # val is a scalar

        Element coordinates are:
         - for tet elements: four barycentric coordinates, the sum of all four
           always equals one. Inside the tet element each of the coordinates is
           betwenn 0 and 1. They are identical to the shape functions of the
           linear tet element.
         - for wedge elements: three barycentric coordinates on the top and
           bottom triangle and the forth going from -1 to 1 from bottom (node
           1..3) to top (node 4..6).
         - for bricks/HEX elements: three coordinates going from -1 to 1.
           (Not yet implemented.)

        @param elemCoords: list of tuples (element number, element coordinates)
        identifying the points for which the element shape function values
        should be calculated. The list might as well be an iterator, in this
        case the additional argument nbPoints must be specified. (Iterators do
        not have a length, which is needed for the progress indicator.)

        @param nbPoints: Optional number of points to process. For diagnostic
        output if the elemCoords argument is an iterator.

        @Note: The function returns an iterator. You can only get each value
        once. If you need it several times, convert to (i.e. store in) a list:
         >>> elemSFuncs = list(mesh.getElemShapeFuncs(elemCoords))

        @Note: The element number of each item in elemCoords will be checked
        whether it's None (indicating a point outside the mesh). Aside from
        that it will only be passed through to the result. As a result the
        element number in fact does not need to be a number. This can be
        (mis-) used if each item in elemCoords corresponds to more than one
        point and more than one element.
        See InterpolMeshToFaceCentroids.shapeFaceCentrWeights for an example.
        """

        # finding elements and shape function values for each point
        try:
            nbPoints = len(elemCoords)
        except TypeError:
            pass
        if nbPoints is None:
            ticker = MsgTicker("Calculating shape function nb %d")
        else:
            ticker = MsgTicker("Calculating shape functions %%d/%d" % nbPoints)
        shapeNotImplemented = defaultdict(int)
        for element, coords in elemCoords:
            ticker.tick()
            if element is None:
                yield (None, [])
                continue

            # determine element shape
            shape = self.elShape[element]
            try:
                func = self.elShapeElemShapeFunc[shape]
            except KeyError:
                shapeNotImplemented[shape] += 1
                yield (None, [])
                continue
            else:
                yield (element, func(coords))
        del ticker
        for shape, cnt in shapeNotImplemented.iteritems():
            msg("WARNING: Element shape %s is not implemented in"
                " Mesh.elShapeElemShapeFunc. This affects %d elements."
                % (shape, cnt))
        return

    def getElemIPInterpolFuncs(self, elemCoords, nbPoints=None):
        """The function returns a tuple (element number, [N1, N2, ..., Nn])
        for each item in points; n is the number of integration points (IPs) on
        the particular element type (shape). N1...Nn are the values of the
        interpolating functions associated with the integration points of the
        element at the given point.

        If any of the points is not in any element of the mesh which would be
        indicated by a element number None in elemCoords, (None, []) is
        returned for that particular point.

        This can be used for interpolating integration point values. The
        following is for illustration purposes only! Use the
        Field.interpolateToPoints() methods in the bae.field_01 module for
        this.

        Suppose we have a mesh and a scalar field defined through values at
        four integrations points per element:
         >>> mesh = ...                                         # a Mesh object
         >>> elemToValues = {1:[5.4,2.1,6.4,3.8], 2:[...], ...} # IP values

        Now we want to interpolate the values to some arbitrary points:
         >>> gridPoints = [[1.0, 0.3, 0.6], [4.0, 1.6, 2.8]]  # some points
         >>> elemCoords = mesh.getElemCoords(gridPoints)
         >>> elemIPFuncs = mesh.getElemIPInterpolFuncs(elemCoords)
         >>> pointValues = [
         ...     sum(weight*ipVal
         ...         for ipVal, weight in zip(elemToValues[elem], weights))
         ...     for elem, weights in elemIPFuncs]

        Write point coordinates and corresponding values to a csv file:
         >>> import csv
         >>> output = csv.writer(open("myoutput.csv", "wb"))
         >>> output.writerow(['x', 'y', 'z', 'value'])
         >>> for pt, val in zip(gridPoints, pointValues):
         >>>     output.writerow(pt+[val])  # val is a scalar

        See L{getElemShapeFuncs} for a description of the element coordinates
        convention.

        @param elemCoords: list of tuples (element number, element coordinates)
        identifying the points for which the element shape function values
        should be calculated.

        @param nbPoints: Optional number of points to process. For diagnostic
        output if the elemCoords argument is an iterator.

        @Note: The function returns an iterator. You can only get each value
        once. If you need it several times, convert to (i.e. store in) a list:
         >>> elemSFuncs = list(mesh.getElemIPInterpolFuncs(elemCoords))
        """

        # some constants (access to local variables is fastest)
        tetQConvMat = self._gpQTetInterpolation

        # finding elements and shape function values for each point
        try:
            nbPoints = len(elemCoords)
        except TypeError:
            pass
        if nbPoints is None:
            ticker = MsgTicker("Calculating shape function (IP) nb %d")
        else:
            ticker = MsgTicker(
                "Calculating shape functions (IP) %%d/%d" % nbPoints)
        shapeNotImplemented = defaultdict(int)
        for element, coords in elemCoords:
            ticker.tick()
            if element is None:
                yield (None, [])
                continue

            # determine element shape
            #... linSFVal = linear shape function value (for specific location
            # in element)
            shape = self.elShape[element]
            if shape=="TET_L":
                yield (element, [1])
                continue
            if shape=="TETHT_L":
                yield (element, [
                    sum(convVal*linSFVal
                        for convVal, linSFVal in izip(convRow, coords))
                    for convRow in self._gpQTetInterpolation_tb])
                continue
            elif shape=="TET_Q":
                yield (element, [
                    sum(convVal*linSFVal
                        for convVal, linSFVal in izip(convRow, coords))
                    for convRow in tetQConvMat])
                continue
            elif shape=="WEDGE_L":
                # valid for cohesive elements with 3 integration points
                yield (element, coords[:3])
            else:
                shapeNotImplemented[shape] += 1
                yield (None, [])
                continue
        del ticker
        for shape, cnt in shapeNotImplemented.iteritems():
            msg("WARNING: Element shape %s is not implemented in"
                " Mesh.getElemIPInterpolFuncs(). This affects %d"
                " elements.\n" % (shape, cnt))
        return

    def getNodalAveragingIPCoeffs(self, elems=None):
        """Provides the data structure needed for nodal interpolation of fields
        defined at element integration points. This is used with
        L{field_01.FieldPosMeshElemIP.getNodeField()}
        to create a field object with position="node" that contains the field
        continuously extrapolated to the nodes. The values at the integration
        points are being extrapolated to each node for each element. Then the
        arithmetic mean of the values coming from all attached elements to a
        particular node is taken as the nodal value.

        This method returns weight factors for all intergration points of
        all elements attached to a particular node. Those weight factors will
        be multiplied to the field values at each integration point and then
        summed up to yield the averaged nodal value.

        @param elems: List or other iterable of element numbers to be
           considered. Note that nodes that are not connected to any of the
           elements in elems won't get a value at all. If None then take all
           elements in self. Defaults to None.

        @Return: A dictionary { node label : [elemIpWeights list]} where
           elemIpWeights contains a [element label, weights list] tuple for
           each element contributing to the respective nodal value. The weights
           list contains the weight for each integration point value (IP)::
             { node : [ [elem, [weight for IP1, weight for IP2, ....]],
                        [elem, [weights ...]], ],   ... }

        @Note: To get the effect of the "don't average on region boundaries"
        option in Abaqus cae supply region by region to the elems argument at
        successive calls to this method.
        """
        # alternative approach:
        #
        # Find the nodal values yielding the smallest global error at the
        # integration points or between the (usually) elementwise linear
        # interpolation from integration points and the FEM function space
        # interpolation that is provided by this method.
        #
        #
        # ideas for improvements:
        #
        # problem: large data structure. This is necessary for compatibility
        # with a yet to write alternative method described above.
        #
        # field class provides iterator yielding all ip values and mesh class
        # provides two sets of functions. One set for averaging nodal values,
        # one set for finding the minimum error solution. Each set contains
        # one function for building up the data structure and one for doing the
        # actual interpolation (IP values -> nodal values). Or put each set
        # into a distinct subclass of a nodal-interpolation-class. This could
        # then also hold the data structure

        if elems is None:
            nbElems = len(self.elNodes)
        else:
            nbElems = len(elems)
        msg("Accumulating connected elements for each node.")
        ticker = MsgTicker("Processing element %%d/%d." % nbElems)

        # nodeElIpWt: { nodelabel : [ [elemlabel, [weight IP1, weight IP2]],...
        nodeElIpWt = defaultdict(list)
        for thisShape, currentElems in self.shapeEl.iteritems():
            if elems is not None:
                currentElems = currentElems.intersection(elems)
            if not currentElems:
                continue
            ipInterpolCoeff = self.elShapeIPInterpolCoeff[thisShape]
            for elem in currentElems:
                ticker.tick()
                for node, ipWeights in izip(self.elNodes[elem],ipInterpolCoeff):
                    nodeElIpWt[node].append([elem, ipWeights])
        del ticker

        nodeElIpWt.default_factory = None    # switch off defaultdict function
        msg("Finalizing weight factors for interpolation.")
        ticker = MsgTicker("Processing node %%d/%d." % len(nodeElIpWt))
        for elemWeightsList in nodeElIpWt.itervalues():
            ticker.tick()
            scale = 1./len(elemWeightsList)
            for elemWeights in elemWeightsList:
                elemWeights[1] = [ w*scale for w in elemWeights[1] ]
        del ticker

        return nodeElIpWt

    #} end of interpolation

    ## ------------------------------------------------------------------------
    ## ------------------------------------------------------------------------
    #{ query mesh properties, other than findAtMesh stuff

    def findNodes(self, pointList, tolerance=0.01):
        """
        returns the node numbers to a list of points assumed to be coordinates
        of the nodes of the mesh

        @param pointList: iterable of points, each point being a list of
        coordinates
        @param tolerance: point finding tolerance. Points closer to a node than
        the specified tolerance are guaranteed to be found and points further
        away than twice the tolerance are not associated with that node. Any
        distance in between leads to an arbitrary result.

        @returns: a list of node numbers, None for each point that could not
        be identified as node

        @note: If the distance between two nodes is close to or smaller than
        tolerance one of them is returned for a queried point nearby. The
        outcome is arbitrary and does not depend on which of the node is closer
        to the queried point. There is no advantage in speed for choosing a
        larger or smaller tolerance.

        @note: For larger meshes most of the work of this method is its
        initialization, so try to avoid to call it multiple times. Rather first
        collect all the points then call it once, if speed is a concern.
        """
        factor = 1.0/tolerance
        # this is a dict
        coordToNode = dict(
            ( tuple(int(round(x*factor)) for x in coord), node)
            for node, coord in self.nodeCoords.iteritems())
        result = list()
        class LocalBreak(Exception): pass
        for pt in pointList:
            pt = tuple(int(round(x*factor)) for x in pt)
            try:
                # check the immidiate block
                try:
                    id = coordToNode[pt]
                except KeyError:
                    pass
                else:
                    raise LocalBreak
                # check neighbouring blocks, over edges and corner points
                for ioffs in [
                    # over faces
                    (1,0,0), (-1,0,0),
                    (0,1,0), (0,-1,0),
                    (0,0,1), (0,0,-1),
                    # over edges
                    (1,1,0), (1,-1,0), (-1,-1,0), (-1,1,0),
                    (1,0,1), (1,0,-1), (-1,0,-1), (-1,0,1),
                    (0,1,1), (0,1,-1), (0,-1,-1), (0,-1,1),
                    # over corner points
                    (1,1,1), (-1,1,1), (-1,-1,1), (1,-1,1),
                    (1,1,-1), (-1,1,-1), (-1,-1,-1), (1,-1,-1),
                    ]:
                    ptnext = (pt[0]+ioffs[0], pt[1]+ioffs[1], pt[2]+ioffs[2])
                    try:
                        id = coordToNode[ptnext]
                    except KeyError:
                        pass
                    else:
                        raise LocalBreak
                id = None
            except LocalBreak:
                pass
            result.append(id)
        return result
    #} end of query mesh properties, other than findAtMesh stuff

    ## ------------------------------------------------------------------------
    ## ------------------------------------------------------------------------
    #{ transform geometry

    def scale(self, scale_factor, scale_origin):
        """Scale the mesh by a certain factor relative to scale_origin

        scale_factor may be a single number in which case it is applied to all
        three dimensions or a list of three numbers with one factor for each
        dimension.

        Example:
         >>> # scale this mesh by 2 in each direction (volume := volume*8)
         >>> mesh.scale(scale_factor=2, scale_origin=[0,0,0])

        Scaling is done by simply moving the coordinates for all nodes.

        After scaling point search has to be initialized again,
        self.point_search_initialized is set to False.

        Possible improvement: To avoid initializing the point search again
        after scaling it would be necessary to also scale self.boundingBox and
        self.cellSize. In case of scale_factor being a list (different for the
        three dimensions) self.cellSize would have to become a list, too.
        """

        if type(scale_factor) in (float, int):
            scale_factor = [scale_factor,scale_factor,scale_factor]

        for coords in self.nodeCoords.itervalues():
            for i in xrange(3):
                coords[i] = ((coords[i]-scale_origin[i])*scale_factor[i]
                             +scale_origin[i])

        # after scaling point search has to be initialized again
        self.point_search_initialized = False

    def translate(self, move_vector):
        """Move the mesh in space by move_vector

        Example:

        >>> # move this mesh up by 10
        >>> mesh.translate(move_vector=[0,0,10])

        Translating is done by simply moving the coordinates of all nodes.

        The point search initialization is modified and doesn't have to be done
        again.
        """

        for coords in self.nodeCoords.itervalues():
            coords[0] += move_vector[0]
            coords[1] += move_vector[1]
            coords[2] += move_vector[2]

        # adapt point search properties
        if self.point_search_initialized and len(self.elNodes)>0:

            # tet Centroids
            for centroid in self.tetCentroids.itervalues():
                vector_modif_add(centroid, move_vector)

            # BoundingBox
            vector_modif_add(self.boundingBox[0], move_vector)
            vector_modif_add(self.boundingBox[1], move_vector)

    def rotate(self, refPoint, rotMat):
        """Rotate the mesh according to the specified transformation matrix
        rotMat relative to the given point refPoint.

        @param refPoint: point around which to rotate the mesh
        @param rotMat: transformation matrix. No check is performed but the
           matrix should perform a pure rotation.
        """
        for coords in self.nodeCoords.itervalues():
            coords[:] = vector_plus(
                refPoint,
                mat_multvec(rotMat, vector(refPoint, coords))
                )

    def rotateX(self, refPoint, alpha):
        """Rotate the mesh around a line parallel to the x axis.

        @param refPoint: point around which to rotate the mesh
        @param alpha: angle in degree, pos == clockwise when looking in
        x-direction.
        """
        alpha = alpha*pi/180
        ca=cos(alpha)
        sa=sin(alpha)
        for coords in self.nodeCoords.itervalues():
            y = coords[1] - refPoint[1]
            z = coords[2] - refPoint[2]
            coords[1] = ca*y-sa*z + refPoint[1]
            coords[2] = sa*y+ca*z + refPoint[2]

    def rotateY(self, refPoint, alpha):
        """Rotate the mesh around a line parallel to the y axis.

        @param refPoint: point around which to rotate the mesh
        @param alpha: angle in degree, pos == clockwise when looking in
        y-direction.
        """
        alpha = alpha*pi/180
        ca=cos(alpha)
        sa=sin(alpha)
        for coords in self.nodeCoords.itervalues():
            x = coords[0] - refPoint[0]
            z = coords[2] - refPoint[2]
            coords[2] = ca*z-sa*x + refPoint[2]
            coords[0] = sa*z+ca*x + refPoint[0]

    def rotateZ(self, refPoint, alpha):
        """Rotate the mesh around a line parallel to the z axis.

        @param refPoint: point around which to rotate the mesh
        @param alpha: angle in degree, pos == anticlockwise when seen from
        above (i.e. looking in neg. z-direction).
        """
        alpha = alpha*pi/180
        ca=cos(alpha)
        sa=sin(alpha)
        for coords in self.nodeCoords.itervalues():
            x = coords[0] - refPoint[0]
            y = coords[1] - refPoint[1]
            coords[0] = ca*x-sa*y + refPoint[0]
            coords[1] = sa*x+ca*y + refPoint[1]

    def warp(self, displacement, scaleFactor=1.0):
        """move each node in self by a certain displacement:

        x_i,new := x_i,old + scaleFactor * displacement_i

        @param displacement: {node number: displacement vector}
        @param scaleFactor: constant scale factor for the displacements
        @note: Nodes not listed in displacement are not touched, nodes listed
          in displacement that don't exist in self are being ignored.
        """
        for node, disp in displacement.iteritems():
            try:
                coords = self.nodeCoords[node]
            except KeyError:
                continue
            coords[0] += disp[0]*scaleFactor
            coords[1] += disp[1]*scaleFactor
            coords[2] += disp[2]*scaleFactor
    #} end of transform geometry

    ## ------------------------------------------------------------------------
    ## ------------------------------------------------------------------------
    #{ other methods extracting parts of the mesh

    def partMeshFromElements(self, elems):
        """Return another mesh with just all the given elements. And only
        attached nodes.

        nodeCoords and elNodes are shallow copies the rest is deep copied

        If elems is empty an empty mesh is returned. No warning!

        If there are no element or no node definitions to any given element
        this element is silently ignored.
        """
        mo = Mesh()

        ticker = MsgTicker(
            "Creating a mesh part of %d out of %d elements: %%s"
            % (len(elems), len(self.elNodes)))

        nodestocopy = set(
            n
            for elem in set(elems).intersection(self.elNodes)
            for n in self.elNodes[elem])
        ticker.msg("Collected %d nodes to copy." % len(nodestocopy))

        for thisShape, thisElems in self.shapeEl.iteritems():
            elNodes = dict(
                (elem, self.elNodes[elem])
                for elem in thisElems.intersection(elems))
            mo.updateElems(elNodes, thisShape)
        ticker.msg("Copied %d elements to new mesh." % len(mo.elNodes))

        mo.nodeCoords.update(
            (node, self.nodeCoords[node])
            for node in nodestocopy.intersection(self.nodeCoords))
        ticker.msg("Copied %d nodes to new mesh." % len(mo.nodeCoords))

        del ticker
        return mo
    #} end of other methods extracting parts of the mesh


#---------------------------------------------------------
# mesh classes for a "structured points" grid
#---------------------------------------------------------
class MeshStructuredPoints(MeshBaseType):
    """Class to hold the data for a "structured points" grid aligned to the
    coordinate axes.

    @ivar gridPtNb: A triple of integers specifying the number of grid
      points in x, y, z direction
    @ivar origin: first point (min x, y, z)
    @ivar spacing: A triple of floats: Point spacing in x, y, z direction
    @ivar strides: A triple of index offsets to the next point in x,y,z
      direction. I.e. In the order of grid points supplied by getPointsIter()
      strides[0] is the offset to the next neighbouring point in x direction,
      strides[1] is next in y direction and strides[2] in z direction.

      Suppose you've got index as a triple of (zero based) x,y,z-indices.
      I.e. index=[i,j,k] for the point at the i-th position in x direction
      j-th position in y direction and k-th in z. Then the following formula
      would give you the grid point index in the order of getPointsIter():
       >>> sum( i*s for i, s in zip(index, self.strides) )

      Note that strides is for information only. Don't set it to arbitrary
      values and expect the class to work correctly!
    """

    def __init__(self, **kwargs):
        """Constructor.

        You may specify either of the following argument combinations:
         - firstPoint, lastPoint and spacing
         - gridPtNb, origin and spacing
         - box and spacing

        @kwarg gridPtNb: A triple of integers specifying the number of grid
          points in x, y, z direction.
          You need to supply spacing as well. And origin might be useful.
        @kwarg origin: First point (min x, y, z). Optional, defaults to [0,0,0]
          You need to supply gridPtNb and spacing as well.
        @kwarg box: A bounding box [[xmin,ymin,zmin], [xmax, ymax, zmax]].
          You need to supply spacing as well.
        @kwarg spacing: A triple of floats: Point spacing in x, y, z direction
          Additionally you need to supply either box or gridPtNb and
          (optionally) origin. May also be just a single number that will then
          be used for x, y and z directions.
        @kwarg firstPoint: first point (min x, y, z) (synonym to origin)
        @kwarg lastPoint: last point (max x, y, z)

        @Note: Might add a single positional argument later. This would then be
          used as initializer. E.g. for copy constructor.
        """

        # replace alias firstPoint, lastPoint -> box
        if "firstPoint" in kwargs and "lastPoint" in kwargs:
            kwargs["box"] = [kwargs["firstPoint"], kwargs["lastPoint"]]
            # note: it's save to delete from the kwargs dict because even if
            # __init__() is called with a dictionary argument, i.e. as
            # __init__(self, **myparams) kwargs will be a copy of (not a
            # reference to) myparams. Automatically...
            del kwargs["firstPoint"], kwargs["lastPoint"]

        # make spacing a list
        try:
            spacing = kwargs["spacing"]
        except KeyError:
            pass
        else:
            if isinstance(spacing, (int, float)):
                kwargs["spacing"] = [spacing,spacing,spacing]

        # now process the different variants of arguments
        if ("box" in kwargs and "spacing" in kwargs):
            gridBox = kwargs["box"]
            if any(x1>x2 for x1, x2 in izip(*gridBox)):
                raise ValueError(
                    "Grid box values don't make sense, lower bound larger than"
                    " upper bound: %s" % gridBox)

            self.spacing = kwargs["spacing"]
            self.origin = gridBox[0]
            gridBoxSize = [x1-x0 for x0, x1 in zip(*gridBox)]
            self.gridPtNb = [
                int(round(x/dx))+1 for x, dx in zip(gridBoxSize, self.spacing)]
        elif all((x in kwargs) for x in ("gridPtNb", "spacing")):
            # origin defaults to [0,0,0]
            self.origin = kwargs.get("origin", [0.0, 0.0, 0.0])
            # all other arguments are mandatory
            # (as long as the copy constructor has not been implemented)
            self.gridPtNb = kwargs["gridPtNb"]
            self.spacing = kwargs["spacing"]
        else:
            raise ValueError(
                "MeshStructuredPoints.__init__ needs at least box and spacing"
                " or gridPtNb and spacing arguments. Got %s." % kwargs.keys())

        # make sure that origin, spacing and gridPtNb are lists
        if not isinstance(self.spacing, list):
            self.spacing = list(self.spacing)
        if not isinstance(self.origin, list):
            self.origin = list(self.origin)
        if not isinstance(self.gridPtNb, list):
            self.gridPtNb = list(self.gridPtNb)

        self.strides = [1, self.gridPtNb[0], self.gridPtNb[0]*self.gridPtNb[1]]
        assert all(x>0 for x in self.gridPtNb)
        assert all(x>0 for x in self.strides)
        assert all(x>0 for x in self.spacing)

    def __cmp__(self, other):
        """Compare to other object of same class based on gridPtNb, spacing and
        origin.
        """
        if type(self)!=type(other):
            msg("MeshStructuredPoints.__cmp__: different type of self (%s)"
                " and other (%s)." % (type(self), type(other)))
            return 1
        if self.gridPtNb != other.gridPtNb:
            msg("MeshStructuredPoints differ in gridPtNb: self %s, other %s"
                % (self.gridPtNb, other.gridPtNb), debugLevel=10)
            return (
                cmp(self.gridPtNb[0]*self.gridPtNb[1]*self.gridPtNb[2],
                    other.gridPtNb[0]*other.gridPtNb[1]*other.gridPtNb[2])
                or 1 )  # make sure it's not equal by chance

        if (length(vector(self.spacing, other.spacing))
            / length(self.spacing) > 1E-6):
            msg("MeshStructuredPoints differ in spacing: self %s, other %s"
                % (self.spacing, other.spacing), debugLevel=10)
            return (
                cmp(self.spacing[0]*self.spacing[1]*self.spacing[2],
                    other.spacing[0]*other.spacing[1]*other.spacing[2])
                or 1 )  # make sure it's not equal by chance

        endPt = [(n-1)*x for n, x in zip(self.gridPtNb, self.spacing)]
        if (length(vector(self.origin, other.origin))
            / length(vector(self.origin,endPt)) > 1E-6):
            msg("MeshStructuredPoints differ in origin: self %s, other %s"
                % (self.origin, other.origin), debugLevel=10)
            return (
                cmp(self.origin[0]*self.origin[1]*self.origin[2],
                    other.origin[0]*other.origin[1]*other.origin[2])
                or 1 )  # make sure it's not equal by chance

        # if all checks succeed then they are equal
        return 0

    def __len__(self):
        """len(mygrid) returns the number of grid points"""
        return self.gridPtNb[0]*self.gridPtNb[1]*self.gridPtNb[2]

    def __str__(self):
        """return a descriptive string of self.
        """
        return (
            "MeshStructuredPoints %d points: %d, %d, %d points in"
            " x,y,z-direction, origin %s, spacing %s"
            % (len(self), self.gridPtNb[0], self.gridPtNb[1], self.gridPtNb[2],
               self.origin, self.spacing))

    def getGridDataString(self):
        """Return a Python command string that can be used as arguments for the
        constructor or for L{getMesh} to rebuild this object.

         >>> kwargStr = box.getGridDataString()
         >>> print kwargStr
         dict(
             firstPoint=[10.00,3.00,1.00],
             lastPoint=[110.00,103.00,11.00],
             spacing=[5,5,5])
         >>> kwargs = eval(kwargStr)
         >>> mesh = getMesh(**kwargs)
        """
        lastPoint = [x0 + (n-1)*dx
                     for x0, n, dx in izip(
                             self.origin, self.gridPtNb, self.spacing)]
        gridDataString = (
            'dict(\n'
            '    firstPoint=[%.2f,%.2f,%.2f],\n'
            '    lastPoint=[%.2f,%.2f,%.2f],\n'
            '    spacing=[%s,%s,%s])'
            % (self.origin[0], self.origin[1], self.origin[2],
               lastPoint[0], lastPoint[1], lastPoint[2],
               self.spacing[0], self.spacing[1], self.spacing[2]))
        return gridDataString

    def getPoint(self, idx):
        """
        @param idx: Might be the integer index in the order of
           self.getPointsIter(). I.e.
            >>> self.getPoint(i)==list(self.getPointsIter())[i]
            True

           Or it might be an i,j,k tuple giving the grid index. I.e.
            >>> idx = list(self.getGridIdxIter())[i]
            >>> self.getPoint(idx)==list(self.getPointsIter())[i]
            True

           Grid indexes are zero based and the first changes fastest.
           In the first of the above cases --when idx is a single number--
           this index might be negative with the common semantics.
           I.e. self.getPoint(-1) returns the point in the corner opposite
           to the origin.
        """
        try:
            i, j, k = idx
        except TypeError:
            # accept negative indexes
            while idx<0:
                idx += len(self)
            k = int(idx / self.strides[2])
            idx -= k*self.strides[2]
            j = int(idx / self.strides[1])
            i = idx - j*self.strides[1]

        return [self.origin[0]+i*self.spacing[0],
                self.origin[1]+j*self.spacing[1],
                self.origin[2]+k*self.spacing[2]]

    def getPointsIter(self):
        """
        Returns an iterator yielding all grid points one after the other. Each
        grid point is a float valued coordinate triple.
        """

        z = self.origin[2]
        for k in xrange(self.gridPtNb[2]):
            y = self.origin[1]
            for j in xrange(self.gridPtNb[1]):
                x = self.origin[0]
                for i in xrange(self.gridPtNb[0]):
                    yield [x, y, z]
                    x += self.spacing[0]
                y += self.spacing[1]
            z += self.spacing[2]

        # astonishingly the following takes twice as long as the
        # yield-loop-version above!
        # return (
        #     [self.origin[0]+i*self.spacing[0],
        #      self.origin[1]+j*self.spacing[1],
        #      self.origin[2]+k*self.spacing[2]]
        #     for k in xrange(self.gridPtNb[2])
        #     for j in xrange(self.gridPtNb[1])
        #     for i in xrange(self.gridPtNb[0]) )

    def getGridIdxIter(self):
        """
        Returns an iterator yielding all grid indexes one after the other. Each
        grid index is a tuple of three integers. The first will be (0,0,0).
        """
        for k in xrange(self.gridPtNb[2]):
            for j in xrange(self.gridPtNb[1]):
                for i in xrange(self.gridPtNb[0]):
                    yield (i, j, k)

    def getGridIdxFromSerialIndex(self, idx):
        """Return a grid index / (i,j,k) tuple for the given serial index.
        I.e:
         >>> grid = MeshStructuredPoints(...)
         >>> for serIdx, gridIdx in enumerate(grid.getGridIdxIter()):
         >>>     assert gridIdx == grid.getGridIdxFromSerialIndex(serIdx)

        Note: The reverse operation is performed by the following expression:
         >>> idx = sum( i*s for i, s in zip((i,j,k), self.strides) )
         >>> # ...or presumably faster:
         >>> idx = i*self.strides[0]+j*self.strides[1]+k*self.strides[2]

        @param idx: a single integer, the zero based serial index of the
           grid point
        @returns: the corresponding i,j,k triple
        """
        k = int(idx / self.strides[2])
        idx -= k*self.strides[2]
        j = int(idx / self.strides[1])
        i = idx - j*self.strides[1]
        return (i,j,k)

    def getIdsInBox(self, box):
        """Returns an iterable (currently implemented: a list) of serial point
        ids for all points inside the given box. Points on the border included.

        The returned serial point ids are single integers that can be converted
        into i,j,k-tuples by means of L{getGridIdxFromSerialIndex}.

        @param box: [[min x, min y, min z], [max x, max y, max z]]. If self is
            a L{MeshStructuredPointsRot}-object then the given box coords must
            be backrotated into the rotated coordinate system. This is exactly
            what we have in the use case of L{Mesh.getElemCoordsForGrid()}.
            So don't overload this function in the MeshStructuredPointsRot
            class.

        @Note: WARNING: uses the rotated coordinate system for
            L{MeshStructuredPointsRot}-objects. And this is the intention!
        """
        ximin = [
            max(0, int(ceil(float(box[0][j] - self.origin[j])/self.spacing[j])))
            for j in 0,1,2]
        ximax = [
            min(self.gridPtNb[j] - 1,
                int(floor(float(box[1][j] - self.origin[j])/self.spacing[j])))
            for j in 0,1,2]

        # Note on speed of the algorithm:
        # Earlier we had a version that concatenated lists created by the range
        # function like that:
        # sum((range(ximin[0]+i1+i2, ximax[0]+i1+i2+1)
        #      for i1 in range(...) for i2 in range(...)), [])
        # This was considerable (>50x) slower!
        return [i0
                for i2 in xrange(ximin[2]*self.strides[2],
                                 ximax[2]*self.strides[2]+1,
                                 self.strides[2])
                for i1 in xrange(ximin[1]*self.strides[1]+i2,
                                 ximax[1]*self.strides[1]+i2+1,
                                 self.strides[1])
                for i0 in xrange(ximin[0]*self.strides[0]+i1,
                                 ximax[0]*self.strides[0]+i1+1,
                                 self.strides[0])]

    def getIdsInInterval(self, x0, x1, axis):
        """Returns an iterable (currently implemented: a list) of components
        of a grid index (one component of the i,j,k-tuple) for grid points
        inside the given interval along the specified axis. Points on the border
        included.

        Service function for
        L{bae.volume_02.VolumeInSurfCheckRays.intersectionGrid}
        and L{bae.volume_02.VolumeBetweenSurfaceAndZPlane.intersectionGrid}
        ...and maybe others...

        @param axis: 0, 1 or 2 for x-, y- or z-axis
        """
        ximin = max(0, int(ceil((x0 - self.origin[axis])/self.spacing[axis])))
        ximax = min(self.gridPtNb[axis] - 1,
                    int(floor((x1 - self.origin[axis])/self.spacing[axis])))
        return range(ximin, ximax+1)

    def getIdsInIntervalOverZ(self, z0, z1):
        """Returns an iterable (currently implemented: a list) of k-indexes
        (the last of the i,j,k-tuple) for grid points inside the given
        z-interval. Points on the border included.

        Deprecated alias for getIdsInInterval(..., axis=2)
        """
        return self.getIdsInInterval(z0, z1, axis=2)

    def getPtIdsWeights(self, points, nbPoints=None):
        """Return point indices and weights for interpolation of values
        defined at the points of the grid.

        This is a generator function. It returns point indices and weights
        for each of the supplied points. Each point that lies inside the grid
        has eight grid points surrounding it. Those grid points have a certain
        index in the order returned by self.getPointsIter(). In order to
        interpolate from values given at the grid points this function delivers
        the indices of the (up to) eight grid points and the weights to be
        multiplied to the values at those grid points. Summing up those weight-
        grid point products you would get the trilinear interpolation of the
        values at the grid points for the specified point.

        @param points: A list or other iterable of point coordinates. If this
        object does not have a length --i.e. len(points)-- then you have to
        supply nbPoints as well for diagnostic output.

        @param nbPoints: If you specify the number of points to process (just
        for diagnostic output, may as well be zero), the points argument may as
        well be an iterator.

        @returns: (ptIds, weights) tuples for each point. ptIds is a list of
        1, 2, 4  or 8 indices in a list constructed from the output of
        self.getPointsIter(). weights is a list of equal length supplying the
        weights. For points found to be outside the grid (None, None) is
        returned.

        @Note: Points on the border might get less than eight point indices.
        Or they might be considered outside the grid. This needs to be
        adjusted when needed.
        """

        if nbPoints is None:
            nbPoints = len(points)

        ticker = MsgTicker(
            "calculating point ids and weights %%d/%d"%nbPoints)
        for point in points:
            ticker.tick()

            # ijk-index of the grid cell (float valued)
            index = [(point[i]-self.origin[i])/self.spacing[i]
                     for i in xrange(3)]

            # check not outside
            if any(index[i]<0 or index[i]>=self.gridPtNb[i]-1
                   for i in xrange(3)):
                yield (None, None)
                continue

            ptIdx0 = sum( int(i)*s for i, s in zip(index, self.strides) )
            ptIds = [ ptIdx0,
                      ptIdx0+self.strides[0],
                      ptIdx0                +self.strides[1],
                      ptIdx0+self.strides[0]+self.strides[1],
                      ptIdx0                                +self.strides[2],
                      ptIdx0+self.strides[0]                +self.strides[2],
                      ptIdx0                +self.strides[1]+self.strides[2],
                      ptIdx0+self.strides[0]+self.strides[1]+self.strides[2] ]

            xis = [i-int(i) for i in index]
            weights = [
                (1-xis[0])*(1-xis[1])*(1-xis[2]),
                xis[0]    *(1-xis[1])*(1-xis[2]),
                (1-xis[0])*xis[1]    *(1-xis[2]),
                xis[0]    *xis[1]    *(1-xis[2]),
                (1-xis[0])*(1-xis[1])*xis[2]    ,
                xis[0]    *(1-xis[1])*xis[2]    ,
                (1-xis[0])*xis[1]    *xis[2]    ,
                xis[0]    *xis[1]    *xis[2]     ]
            yield (ptIds, weights)

        del ticker
        return

    def getPtIdsWeightsFromIdxFloat(self, idxFloat):
        """Return point indices and weights for interpolation of values
        defined at the grid-points. The point to interpolate to is identified
        by its i-j-k-index tuple.


        @returns: (ptIds, weights) tuples for each point. ptIds is a list of
        1, 2, 4  or 8 indices in a list constructed from the output of
        self.getPointsIter(). weights is a list of equal length supplying the
        weights. For points found to be outside the grid (None, None) is
        returned.
        """
        # check not outside
        if any(idxFloat[i]<0 or idxFloat[i]>=self.gridPtNb[i]-1
               for i in xrange(3)):
            return (None, None)

        ptIdx0 = sum( int(i)*s for i, s in zip(idxFloat, self.strides) )
        ptIds = [ ptIdx0,
                  ptIdx0+self.strides[0],
                  ptIdx0                +self.strides[1],
                  ptIdx0+self.strides[0]+self.strides[1],
                  ptIdx0                                +self.strides[2],
                  ptIdx0+self.strides[0]                +self.strides[2],
                  ptIdx0                +self.strides[1]+self.strides[2],
                  ptIdx0+self.strides[0]+self.strides[1]+self.strides[2] ]

        xis = [i-int(i) for i in index]
        weights = [
            (1-xis[0])*(1-xis[1])*(1-xis[2]),
            xis[0]    *(1-xis[1])*(1-xis[2]),
            (1-xis[0])*xis[1]    *(1-xis[2]),
            xis[0]    *xis[1]    *(1-xis[2]),
            (1-xis[0])*(1-xis[1])*xis[2]    ,
            xis[0]    *(1-xis[1])*xis[2]    ,
            (1-xis[0])*xis[1]    *xis[2]    ,
            xis[0]    *xis[1]    *xis[2]     ]
        return (ptIds, weights)

    def getBoundingBox(self):
        """Returns a L{BoundingBox<bae.misc_01.BoundingBox>} with self.origin
        being one corner and the last point of the point grid being the
        opposite.
        """
        return BoundingBox([
            self.origin,
            [x0 + (n-1)*dx for x0, n, dx in izip(
                self.origin, self.gridPtNb, self.spacing)]
            ])

    def getClosestGridPointIdx(self, point):
        """Find the grid point closest to the given x,y,z coordinates.
        Return a [i,j,k]-list.

        See also L{self.getPtIdsWeights} if you also want the actual position
        relative to neighpouring grid points (i.e. weight faktors for
        interpolation.

        Note: This function returns the closest grid point even if the given
        point is far away. You may want to check if the point is inside the
        grid box like this: (Or add 0.5*spacing[i]...?)
         >>> grid = MeshStructuredPoints(...)
         >>> point = [x, y, z]
         >>> if grid.getBoundingBox().pointIsInside(point):
         >>>     ijk = grid.getClosestGridPointIdx(point)
        """
        # ijk-index of the grid cell, rounded to the nearest int
        index = [int(round((point[i]-self.origin[i])/self.spacing[i]))
                 for i in xrange(3)]
        index = [max(i, 0) for i in index]
        index = [min(i, m-1) for i, m in izip(index, self.gridPtNb)]
        return index

    def findNeighbourPoints(self, pointList, radius, maxNbPoints=None):
        """
        @param pointList: list of xyz-lists
        @param radius: search radius
        @param maxNbPoints: only evaluate the (approximately) closest
           N points (N=maxNbPoints).

        @return: list of (ptIndex in self, [ptIds in pointsList])-tuples
        """
        R = float(radius)

        #--- deltaIdxXXXX: ijk-offsets for cells inside radius
        # Suppose we have a sphere of the given radius around each grid point.
        # Now consider a grid cell identified by it's lower-left-point/origin.
        # Find all grid points whose sphere intersects or fully contains the
        # grid cell.
        # deltaIdx is a prerequisite for the two variables deltaIdxPartial
        # and deltaIdxCompl.
        # Stores the maximum indexes, i.e the circumference like that:
        # [[i_max for j=0,1,..,j_max]_k=0 , [i_max for j=0,1,..,j_max]_k=1,
        #   ..., [i_max for j=0,1,..,j_max]_k=k_max]
        deltaIdx = []
        for k in range(int(R/self.spacing[2])+1):
            z2 = (k*self.spacing[2])**2
            deltaIdx.append([
                int(sqrt(R**2-z2-(j*self.spacing[1])**2)/self.spacing[0])
                for j in range(int(sqrt(R**2-z2)/self.spacing[1])+1)])
        msg("MeshStructuredPoints.findNeighbourPoints, radius=%g,"
            " deltaIdx:\n%s" % (radius, deltaIdx), debugLevel=10)

        #--- deltaIdxPartial: ijk-offsets for cells *partially* inside radius
        # With the spheres of the given radius around each grid point:
        # Find all grid points whose sphere intersects the grid cell and store
        # their index offsets relative to the cells origin in deltaIdxPartial.
        # Only store the minimum and maximum indexes, i.e the circumference
        # so that the following nested loops yields minimum and maximum
        # i-indexes with corresponding j and k indexes:
        #  >>> for k, jList in deltaIdxPartial:
        #  >>>     for j, iMin, iMax in jList:
        #  >>>         print "%d<=i<=%d; j=%d, k=%d" % (-iMin, iMax, j, k)
        deltaIdxPartial = []
        for k, jList in enumerate(
                deltaIdx[::-1]+deltaIdx, start=1-len(deltaIdx)):
            deltaIdxPartial.append((
                k, [(j, -iMax, iMax+1)
                    for j, iMax in enumerate(
                            jList[::-1]+jList, start=1-len(jList))]))
        msg("MeshStructuredPoints.findNeighbourPoints, radius=%g,"
            " deltaIdxPartial:\n%s" % (radius, deltaIdxPartial), debugLevel=10)

        #--- deltaIdxCompl: ijk-offsets for cells *completely* inside radius
        # With the spheres of the given radius around each grid point:
        # Find all grid points whose sphere fully contains the grid cell and
        # store their index offsets relative to the cells origin in
        # deltaIdxCompl.
        # Only store the minimum and maximum indexes, i.e the circumference
        # so that the following nested loops yields minimum and maximum
        # i-indexes with corresponding j and k indexes:
        #  >>> for k, jList in deltaIdxCompl:
        #  >>>     for j, iMin, iMax in jList:
        #  >>>         print "%d<=i<=%d; j=%d, k=%d" % (-iMin, iMax, j, k)
        deltaIdxCompl = []
        for k, jList in enumerate(
                deltaIdx[:0:-1]+deltaIdx[1:], start=2-len(deltaIdx)):
            deltaIdxCompl.append((
                k, [(j, 1-iMax, iMax)
                    for j, iMax in enumerate(
                            jList[:0:-1]+jList[1:], start=2-len(jList))]))
        msg("MeshStructuredPoints.findNeighbourPoints, radius=%g,"
            " deltaIdxCompl:\n%s" % (radius, deltaIdxCompl), debugLevel=10)

        ##### CONSTRUCTION SITE AHEAD ####

        # only for checking intermediate results:
        from bae.vecmath_01 import length, vector
        print "#### deltaIdxPartial:", deltaIdxPartial
        print "partially in:"
        for k, jList in deltaIdxPartial:
            for j, iMin, iMax in jList:
                print "%d<=i<=%d; j=%d, k=%d" % (-iMin, iMax, j, k)
                a1 = [iMin*self.spacing[0], j*self.spacing[1], k*self.spacing[2]]
                a2 = [iMax*self.spacing[0], j*self.spacing[1], k*self.spacing[2]]
                for pt in pointList:
                    print "... d1=%g, d2=%g" % (
                        length(vector(pt, a1)), length(vector(pt, a2)))

        print "#### deltaIdxCompl:", deltaIdxCompl
        print "completely in:"
        for k, jList in deltaIdxCompl:
            for j, iMin, iMax in jList:
                print "%d<=i<=%d; j=%d, k=%d" % (-iMin, iMax, j, k)
                a1 = [iMin*self.spacing[0], j*self.spacing[1], k*self.spacing[2]]
                a2 = [iMax*self.spacing[0], j*self.spacing[1], k*self.spacing[2]]
                for pt in pointList:
                    d1 = length(vector(pt, a1))
                    d2 = length(vector(pt, a2))
                    if d1<=radius and d2<=radius:
                        ok = "ok"
                    else:
                        ok = "NOT OK!!!"
                    print "... d1=%g, d2=%g, %s" % (d1, d2, ok)

        raise NotImplementedError()
        return
        """
        ##### CONSTRUCTION SITE AHEAD ####

        # some explanantion
        gridIdToPtIdsComplete = defaultdict(list)
        #.... same with partial


        # loop over pointList
        for ptIdx, point in enumerate(pointList):

            # ijk-index of the grid cell, rounded to the nearest int
            index = [int(round((point[i]-self.origin[i])/self.spacing[i]))
                     for i in xrange(3) ]

            # loop over deltaIdxCompl
            for dijk in deltaIdxCompl:
                currentIdx = .... index + dijk
                # filter points outside grid
                if any((i<0 or i>imax) for i, imax in izip(currentidx, self.gridPtNb):
                    continue
                gridIdToPtIdsComplete[currentIdx].append(ptIdx)
            # same for partial


        # loop over grid points
        # loop over grid points in gridIdToPtIdsComplete and gridIdToPtIdsPartial


        # calc distances
        # get maxNbPoints closest within radius

        """

        return []

    def getPointAlignedGrid(self, newGrid):
        """Adjust the given new grid with the same spacing as self.
        The resulting grid has points that coincide exactly.

        The translation is as small as possible to achieve this goal.

        @param newGrid: a L{MeshStructuredPoints} object, the grid to be
           slightly translated.

        @returns: approximately the same as given by newGrid with a small shift
           such that (at least some) grid points coincide. (Or would coincide
           if the grids where continued further over their actual bounds.)
        """

        # check: spacing is the same
        if newGrid.spacing != self.spacing:
            raise ValueError(
                "ERROR in MeshStructuredPoints.getPointAlignedGrid:"
                " spacings differ! Old %s, new %s"
                % (self.spacing, newGrid.spacing))

        # get the closest grid point in new grid for the origin of self
        oldToNew = newGrid.getClosestGridPointIdx(self.origin)
        msg("Offset for the old grid in the new one: %s" % oldToNew,
            debugLevel=1)

        # check origin of newGrid
        newGridOrigAligned = [
            x0-ii*dx
            for x0, ii, dx in zip(self.origin, oldToNew, self.spacing)]

        return MeshStructuredPoints(
            origin=newGridOrigAligned,
            gridPtNb=newGrid.gridPtNb,
            spacing=newGrid.spacing)


class MeshStructuredPointsRot(MeshStructuredPoints):
    """Class to hold the data for a "structured points" grid that is rotated
    against the coordinate system.

    A right handed cartesian coordinate system is defined by three orthogonal
    base unit vectors u,v,w and the origin. This is the rotated "physical"
    coordinate system of the grid. ("Physical" means it has the global
    "physical" length unit and is still orthogonal and right handed.)

    The base unit vectors u,v,w form the columns of the transformation matrix
    L{rotmat}. To transform the rotated coordinates [x*, y*, z*] of a point
    into the corresponding global coordinates [x,y,z] (left-)multiply this
    rotation matrix and add the origin:
    [x,y,z] = rotmat * [x*, y*, z*] + origin

    Multiplying each base unit vectors u,v,w with the corresponding grid
    spacing gives non-unit grid base vectors a,b,c. They form the columns of
    the grid transformation matrix L{gridmat}.

    To get the global coordinates [x,y,z] of a point identified by it's
    grid-index-triple [i,j,k] (left-)multiply the grid transformation matrix
    and add the origin: [x,y,z] = gridmat * [i,j,k] + origin

    Note that the last paragraphs are intended for descriptive purposes.
    To actually get coordinates of grid points use the L{getPoint}- or
    L{getPointsIter}-methods.

    @ivar origin: Global x,y,z coordinates of point with index (0,0,0).
    @ivar originBackRotated: The grid can be constructed along global x,y,z
      coordinates from this alternative origin and then rotated around
      the centre of the grid.

    @ivar rotmat: Transformation matrix rotating from the u,v,w-coordinate
      system to the global x,y,z system. rotmat[i][0] is the unit vector u,
      rotmat[i][1] is v and rotmat[i][2] is w.
    @ivar gridmat: Transformation matrix rotating and scaling from grid
      indexes to the global x,y,z system. gridmat[i][0] is the grid vector a,
      gridmat[i][1] is b and gridmat[i][2] is c. To get the global coordinates
      x,y,z of the grid point identified by its grid indexes i,j,k do:
      [x,y,z] = origin + gridmat * [i,j,k]
      See the formula in L{getPoint}.

    @ivar gridPtNb: A triple of integers specifying the number of grid
      points in a,b,c (u,v,w) direction
    @ivar spacing: A triple of floats: Point spacing in a,b,c direction
    @ivar strides: A triple of index offsets to the next point in a,b,c
      direction. I.e. In the order of grid points supplied by getPointsIter()
      strides[0] is the offset to the next neighbouring point in direction a,
      strides[1] is next in direction b and strides[2] in direction c.

      Suppose you've got index as a triple of (zero based) a,b,c-indices.
      I.e. index=[i,j,k] for the point at the i-th position in direction a
      j-th position in direction b and k-th in c. Then the following formula
      would give you the grid point index in the order of getPointsIter():
       >>> sum( i*s for i, s in zip(index, self.strides) )

      Note that strides is for information only. Don't set it to arbitrary
      values and expect the class to work correctly!
    """

    # self.__len__ and self.getGridIdxIter inherited from base class
    # MeshStructuredPoints

    def __init__(self, **kwargs):
        """Constructor.

        You may specify either of the following argument combinations:
         - gridPtNb, origin, baseVectors and spacing
         - origin, edgeVectors and spacing
         - gridPtNb, origin, spacing, rotVector and deg
         - gridPtNb, origin, spacing, rotString
         - firstPoint, lastPoint, spacing, rotString
        In any case originBackRotated can be provided instead of origin.

        @kwarg gridPtNb: A triple of integers specifying the number of grid
          points in a, b, c direction.
        @kwarg origin: Global x,y,z coordinates of point with index (0,0,0).
        @kwarg originBackRotated: The grid can be constructed along global
          x,y,z coordinates from this alternative origin and then rotated
          around the centre of the grid.
        @kwarg firstPoint: synonym for originBackRotated. Make a box from
          firstPoint to lastPoint in the gloabl x,y,z-coordinate system then
          rotate it around the centre of this box and you'll get the real final
          rotated grid box.
        @kwarg lastPoint: back-rotated final grid point; see description of the
          firstPoint argument.
        @kwarg baseVectors: A list of the three base vectors u,v,w. They will
          be scaled to unit length and made orthogonal and then form the
          columns of the instance attribute rotmat.
        @kwarg spacing: A triple of floats: Point spacing in u, v, w direction.
          May also be just a single number that will then be used for all
          directions.
        @kwarg edgeVectors: A list of the three vectors A,B,C spanning the
          grid box from the origin. They should be (at least approximately)
          orthogonal. A,B,C will be scaled to unit length and made orthogonal
          to get u,v,w which then form the columns of the instance attribute
          rotmat. Their original lengths are the edge lengths of the grid box.
          Divided by the corresponding spacing that gives the number of grid
          points along the axes.
        @kwarg rotVector: Rotation vector gives the rotaion axis and
          orientation of the rotation angle deg. Will be scaled to unit
          length.
        @kwarg deg: Rotation angle.
        @kwarg rotString: A string describing the rotation of the box. Format:
          Rxxx+yyy-zzzRaa_a, R-xx+yy+zzR-aa. The R's stay verbatim, xxx, yyy,
          zzz give the vector of the rotaion axis; they are integer numbers
          with arbitrary lengths; xxx can contain a minus sign; + and - in
          front of yyy and zzz according to sign of the component. aaa after
          the second R give the angle; optional decimal mpoint replaced by
          underscore '_'; optional minus sign.

        @Note: Might add a single positional argument later. This would then be
          used as initializer. E.g. for copy constructor.
        """

        # get spacing, make it a list
        try:
            spacing = kwargs["spacing"]
        except KeyError:
            raise ValueError(
                "We need a spacing parameter for"
                " MeshStructuredPointsRot.__init__().")
        else:
            if isinstance(spacing, (int, float)):
                self.spacing = [spacing,spacing,spacing]
            else:
                self.spacing = list(spacing)  # make a copy
            del spacing

        # replace alias firstPoint -> originBackRotated
        # note: it's save to delete from the kwargs dict because even if
        # __init__() is called with a dictionary argument, i.e. as
        # __init__(self, **myparams) kwargs will be a copy of (not a
        # reference to) myparams. Automatically...
        if "firstPoint" in kwargs:
            if "originBackRotated" in kwargs:
                raise ValueError(
                    "MeshStructuredPointsRot.__init__(): wrong arguments."
                    " firstPoint (=%s) is an alias for originBackRotated (=%s),"
                    " you can't supply both."
                    % (kwargs["firstPoint"], kwargs["originBackRotated"]))
            kwargs["originBackRotated"] = kwargs["firstPoint"]
            del kwargs["firstPoint"]

        # process rotString argument if present
        if "rotString" in kwargs:
            if "rotVector" in kwargs or "deg" in kwargs:
                raise ValueError(
                    "MeshStructuredPointsRot.__init__(): wrong arguments."
                    " You can't have a rotString argument and at the same time"
                    " rotVector or deg.")
            rotVector, deg = self.parseRotString(kwargs["rotString"])
            kwargs["rotVector"] = rotVector
            kwargs["deg"] = deg
            del kwargs["rotString"]

        # now process the different variants of arguments
        if all((x in kwargs) for x in (
                "gridPtNb", "baseVectors", "spacing")):
            # normalize and orthogonalize
            self.rotmat = mat_orthoNormBaseGS(kwargs["baseVectors"])
            self.rotmat = mat_transpose(self.rotmat)
            self.gridmat = [ [x*s for x, s in izip(a, self.spacing)]
                             for a in self.rotmat ]
            self.gridPtNb = list(kwargs["gridPtNb"])  # make a copy
        elif all((x in kwargs) for x in (
                "edgeVectors", "spacing")):
            gridBoxSize = [length(x) for x in kwargs["edgeVectors"]]
            self.rotmat = mat_orthoNormBaseGS(kwargs["edgeVectors"])
            self.rotmat = mat_transpose(self.rotmat)
            self.gridmat = [ [x*s for x, s in izip(a, self.spacing)]
                             for a in self.rotmat ]
            self.gridPtNb = [
                int(round(x/dx))+1 for x, dx in zip(gridBoxSize, self.spacing)]
        elif all((x in kwargs) for x in (
                "gridPtNb", "spacing", "rotVector", "deg")):
            self.rotmat = trans_rodrigues(kwargs["rotVector"], kwargs["deg"])
            self.gridmat = [ [x*s for x, s in izip(a, self.spacing)]
                             for a in self.rotmat ]
            self.gridPtNb = list(kwargs["gridPtNb"])  # make a copy
        elif all((x in kwargs) for x in (
                 "originBackRotated", "lastPoint", "spacing",
                 "rotVector", "deg")):
            # note: originBackRotated is a synonym for firstPoint
            self.rotmat = trans_rodrigues(kwargs["rotVector"], kwargs["deg"])
            self.gridmat = [ [x*s for x, s in izip(a, self.spacing)]
                             for a in self.rotmat ]
            gridBoxSize = [
                (x1-x0)
                for x0, x1 in zip(
                    kwargs["originBackRotated"], kwargs["lastPoint"])]
            self.gridPtNb = [
                int(round(x/dx))+1 for x, dx in zip(gridBoxSize, self.spacing)]
        else:
            raise ValueError(
                "MeshStructuredPointsRot.__init__ needs a proper set of"
                " arguments. Got %s." % kwargs.keys())

        # get origin and originBackRotated (calculate the other)
        gridBoxSize = [(n-1)*s for n, s in izip(self.gridPtNb, self.spacing)]
        relCentre = vector_scale(gridBoxSize, 0.5)  # ... of the unrotated grid
        offsetOrig = vector_minus(mat_multvec(self.rotmat,relCentre),relCentre)
        if "origin" in kwargs:
            self.origin = list(kwargs["origin"])
            # calculate originBackRotated
            # move the origin such that you can first create the grid along
            # global coordinate axes and then rotate around the centre of the
            # grid
            self.originBackRotated = vector_plus(self.origin, offsetOrig)
        elif "originBackRotated" in kwargs:
            self.originBackRotated = list(kwargs["originBackRotated"])
            # calculate origin
            # move originBackRotated such that you can first create the grid
            # along global coordinate axes and then rotate around the origin
            self.origin = vector_minus(self.originBackRotated, offsetOrig)
        else:
            raise ValueError(
                "We either need a parameter origin or originBackRotated for"
                " MeshStructuredPointsRot.__init__().")
        del gridBoxSize, relCentre, offsetOrig

        # calculate strides, check values to be positive
        self.strides = [1, self.gridPtNb[0], self.gridPtNb[0]*self.gridPtNb[1]]
        assert all(x>0 for x in self.gridPtNb)
        assert all(x>0 for x in self.strides)
        assert all(x>0 for x in self.spacing)

    def __cmp__(self, other):
        """Compare to other object of same class based on gridPtNb, origin and
        gridmat.
        """
        if type(self)!=type(other):
            msg("MeshStructuredPointsRot.__cmp__: different type of self (%s)"
                " and other (%s)." % (type(self), type(other)))
            return 1
        if self.gridPtNb != other.gridPtNb:
            msg("MeshStructuredPointsRot differ in gridPtNb: self %s, other %s"
                % (self.gridPtNb, other.gridPtNb), debugLevel=10)
            return (
                cmp(self.gridPtNb[0]*self.gridPtNb[1]*self.gridPtNb[2],
                    other.gridPtNb[0]*other.gridPtNb[1]*other.gridPtNb[2])
                or 1)  # make sure it's not equal by chance

        endPt = [(n-1)*x for n, x in zip(self.gridPtNb, self.spacing)]
        if (length(vector(self.origin, other.origin))
            / length(vector(self.origin,endPt)) > 1E-6):
            msg("MeshStructuredPointsRot differ in origin: self %s, other %s"
                % (self.origin, other.origin), debugLevel=10)
            return (
                cmp(self.origin[0]*self.origin[1]*self.origin[2],
                    other.origin[0]*other.origin[1]*other.origin[2])
                or 1)  # make sure it's not equal by chance

        if any((length(vector(selfGridVec, otherGridVec))
                / length(self.spacing) > 1E-6)
               for selfGridVec, otherGridVec in izip(
                       self.gridmat, other.gridmat)):
            msg("MeshStructuredPointsRot differ in gridmat: self %s, other %s"
                % (self.gridmat, other.gridmat), debugLevel=10)
            return (
                cmp(mat_determinant(self.gridmat),
                    mat_determinant(other.gridmat))
                or 1)  # make sure it's not equal by chance

        # if all checks succeed then they are equal
        return 0

    def __str__(self):
        """return a descriptive string of self.
        """
        gridVec = mat_transpose(self.gridmat)

        return (
            "MeshStructuredPointsRot %d points: %d points in direction %s,"
            " %d in %s, %d in %s, origin %s, spacing %s"
            % (len(self),
               self.gridPtNb[0], gridVec[0],
               self.gridPtNb[1], gridVec[1],
               self.gridPtNb[2], gridVec[2],
               self.origin, self.spacing))

    def getGridDataString(self):
        """Return a Python command string that can be used as arguments for the
        constructor or for L{getMesh} to rebuild this object.

         >>> kwargStr = box.getGridDataString()
         >>> print kwargStr
         dict(
             firstPoint=[10.00,3.00,1.00],
             lastPoint=[110.00,103.00,11.00],
             spacing=[5,5,5])
         >>> kwargs = eval(kwargStr)
         >>> mesh = getMesh(**kwargs)
        """
        lastPoint = [
            x0 + (n-1)*dx
            for x0, n, dx in izip(
                    self.originBackRotated, self.gridPtNb, self.spacing)]
        gridDataString = (
            'dict(\n'
            '    firstPoint=[%.2f,%.2f,%.2f],\n'
            '    lastPoint=[%.2f,%.2f,%.2f],\n'
            '    spacing=[%s,%s,%s],\n'
            '    rotString="%s")'
            % (self.originBackRotated[0], self.originBackRotated[1],
               self.originBackRotated[2],
               lastPoint[0], lastPoint[1], lastPoint[2],
               self.spacing[0], self.spacing[1], self.spacing[2],
               self.getRotString()))
        return gridDataString

    def getRotString(self):
        """Return a string for the grid (-box) name to identify the rotation
        of the box. In particular this string contains the values to be inserted
        into the Voxler transform filter to rotate corresponding gridboxes into
        the global coordinate system.
        """
        k, alpha = mat_toRodrigues(self.rotmat)
        # round rotation vector accuracy: 2-3 digits
        kInt = [int(round(x*100)) for x in k]
        # try to reduce number of digits
        for i in range(2):
            if all( ((x % 10)==0) for x in kInt ):
                kInt = [int(round(x/10)) for x in kInt]
            else:
                break
        # round angle (degrees)
        aInt = int(round(alpha))
        return "R%d%+d%+dR%d" % (kInt[0], kInt[1], kInt[2], aInt)

    @staticmethod
    def parseRotString(rotString):
        """Extract vector of rotation axis and rotation angle (in degrees) from
        a descriptive string similar to 'R10-34+100R12'.

        @kwarg rotString: A string describing the rotation of the box. Format:
          Rxxx+yyy-zzzRaa_a, R-xx+yy+zzR-aa. The R's stay verbatim, xxx, yyy,
          zzz give the vector of the rotaion axis; they are integer numbers
          with arbitrary lengths; xxx can contain a minus sign; + and - in
          front of yyy and zzz according to sign of the component. aaa after
          the second R give the angle; optional decimal point replaced by
          underscore '_'; optional minus sign.

        @returns: (rotVector, angle)-tuple; rotVector is *not* normalized to
          unit length, angle is in degrees.
        """
        rex = re.compile(
            r"R(?P<x>-?\d+)(?P<y>[-+]\d+)(?P<z>[-+]\d+)"
            r"R(?P<deg>-?\d+_?\d*)")
        res = rex.match(rotString)
        if not res:
            raise ValueError(
                "could not parse rotString argument: %s" % rotString)
        rotVector = map(int, res.groups()[:3])
        deg = float(res.group("deg").replace("_", "."))
        return rotVector, deg

    def getPoint(self, idx):
        """
        @param idx: Might be the integer index in the order of
           self.getPointsIter(). I.e.
            >>> self.getPoint(i)==list(self.getPointsIter())[i]
            True

           Or it might be an i,j,k tuple giving the grid index. I.e.
            >>> idx = list(self.getGridIdxIter())[i]
            >>> self.getPoint(idx)==list(self.getPointsIter())[i]
            True

           Grid indexes are zero based and the first changes fastest.
           In the first of the above cases --when idx is a single number--
           this index might be negative with the common semantics.
           I.e. self.getPoint(-1) returns the point in the corner opposite
           to the origin.
        """
        try:
            i, j, k = idx
        except TypeError:
            # accept negative indexes
            while idx<0:
                idx += len(self)
            k = int(idx / self.strides[2])
            idx -= k*self.strides[2]
            j = int(idx / self.strides[1])
            i = idx - j*self.strides[1]

        g = self.gridmat  # abbreviation
        return [self.origin[0] + g[0][0]*i + g[0][1]*j + g[0][2]*k,
                self.origin[1] + g[1][0]*i + g[1][1]*j + g[1][2]*k,
                self.origin[2] + g[2][0]*i + g[2][1]*j + g[2][2]*k]

    def getPointsIter(self):
        """
        Returns an iterator yielding all grid points one after the other. Each
        grid point is a float valued coordinate triple.
        """
        # some abbreviations
        o = self.origin
        g = self.gridmat
        ga = [gg[0] for gg in g]

        for k in xrange(self.gridPtNb[2]):
            for j in xrange(self.gridPtNb[1]):
                # starting pt at i=0
                pt = [o[0] + g[0][1]*j + g[0][2]*k,
                      o[1] + g[1][1]*j + g[1][2]*k,
                      o[2] + g[2][1]*j + g[2][2]*k]
                yield pt
                for i in xrange(self.gridPtNb[0]-1):
                    pt = vector_plus(pt, ga)  # we need to make a new pt
                    yield pt

    def getIdsInBox(self, box):
        """Not implemented. Just a placeholder to mask
        L{MeshStructuredPoints.getIdsInBox}()

        Returns an iterable (currently implemented: a list) of point ids for
        all points inside the given box. Points on the border included.

        This function might not be needed for the MeshStructuredPointsRot
        class. The only use so far is in L{Mesh.getElemCoordsForGrid()}
        and here we need MeshStructuredPoints.getIdsInBox.
        """
        raise NotImplementedError()

    def getPtIdsWeights(self, points, nbPoints=None):
        """Not implemented. Just a placeholder to mask
        L{MeshStructuredPoints.getPtIdsWeights()}
        """
        raise NotImplementedError()

    def getBoundingBox(self):
        """Returns a L{BoundingBox<bae.misc_01.BoundingBox>} that contains
        this rotated grid box (all corner points).
        """
        corners = (
            vector_sum(self.origin, *[
                [(n-1)*ijk*self.gridmat[ii][idx] for ii in range(3)]
                for idx, (n, ijk) in enumerate(izip(self.gridPtNb, (i,j,k)))])
            for k in (0,1) for j in (0,1) for i in (0,1))
        # ### DEBUG
        # corners = list(corners)
        # print "####### corners:"
        # print "\n".join(str(x) for x in corners)
        # ### DEBUG
        bb = BoundingBox()
        bb.update(corners)
        return bb


#---------------------------------------------------------
# mesh class for unstructured points as "grid"-instance
#---------------------------------------------------------
class MeshUnstructuredPoints(MeshBaseType, list):
    """Class to hold unstructured points as "grid"-instance.

    This is a subclass of list and must be initialized with a list of point
    coordinates
    """

    # not needed as long as nothing else happens in the constructor...
    # def __init__(self, pointsCoords):
    #     """Constructor.
    #     """
    #     list.__init__(self, pointsCoords)

    def __eq__(self, other):
        """Compare to other object of same class.

        @Note: This is rather expensive as we are calculating relative
        differences of the points in self and other.
        """
        if type(self)!=type(other):
            msg("MeshUnstructuredPoints.__eq__: different type of self (%s)"
                " and other (%s)." % (type(self), type(other)))
            return False

        charLen = max(length(x) for x in self)
        if any((length(vector(selfPt, otherPt))/charLen > 1E-6)
               for selfPt, otherPt in izip(self, other)):
            msg("MeshUnstructuredPoints differ: Different points."
                "\nself: %s\nother: %s" % (self[:100], other[:100]),
                debugLevel=10)
            return False

        # if all checks succeed then they are equal
        return True

    def __ne__(self, other):
        return not(self.__eq__(other))

    def getPointsIter(self):
        """
        Returns an iterator yielding all grid points one after the other. Each
        grid point is a float valued coordinate triple.
        """
        return iter(self)

    # alias needed for compatibility with L{MeshStructuredPoints} classes
    getPoint = list.__getitem__

    def getIdsInBox(self, box):
        """Returns an iterable (currently implemented: a list) of point ids for
        all points inside the given box. Points on the border included.

        @param box: [[min x, min y, min z], [max x, max y, max z]].

        @Note: If the points of self are changed in any way between subssequent
        calls to this function then the instance attribute self.kdtree must
        be removed, i.e. "del self.kdtree".
        """

        # if not done already initialize kdtree for self
        try:
            kdt = self.kdtree
        except AttributeError:
            msg("Initializing kd-tree.")
            kdt = KDTree(self)
            self.kdtree = kdt

        # find grid points in this bounding box: with KDTree
        return kdt.boxSearch(box)

    def getBoundingBox(self):
        """Returns a L{BoundingBox<bae.misc_01.BoundingBox>} with self.origin
        being one corner and the last point of the point grid being the
        opposite.
        """
        try:
            bb = self.boundingBox
        except AttributeError:
            bb = BoundingBox()
            bb.update(self)
            self.boundingBox = bb
        return bb


#---------------------------------------------------------
# interpolation class for structured points interpolated from a field defined
# on a tet mesh
#---------------------------------------------------------
class InterpolMeshBaseClass(object):
    """For interpolating field values from an unstructured (FEM-) mesh to
    certain points the derived children of this base class deliver
    interpolation data, in particular element and node numbers and weighting
    factors.

    When interpolating field values from a mesh to certain points/vertices
    you need to determine in which element the point sits and where exactly
    in the element it sits. The interpolation classes purpose is to
    generate/calculate this data and where this is costly also store and
    retrieve this interpolation data for and on subsequent uses. The actual
    interpolation is then done by Field.interpolate functions of the various
    flavours e.g. L{nodal<bae.field_01.FieldPosMeshNode.interpolate>} or
    L{integration point based<bae.field_01.FieldPosMeshElemIP.interpolate>}
    fields.

    This base class mainly provides the means to save and restore expensive
    interpolation data.

    @ivar grid: A MeshStructuredPoints or MeshStructuredPointsRot or
      MeshUnstructuredPoints object identifying the output points for the
      interpolation.
    @ivar mesh: A (tet-) Mesh object. Contains all node and element
      definitions involved in the interpolation.
    @ivar elemCoords: A list of tuples (element number, element coordinates)
      for the points in self.grid.

      elemCoords may contain tuples with None as element number, as
      Mesh.getElemCoords() provides for points outside the mesh.
    @cvar pickleFileFormat: This string is stored as first item in the
      pickle file to check compatibility between the pickle file and the
      version of the InterpolMeshToPoints class currently in use. This
      string has to be modified if the format of the pickle file is changed.
      Otherwise it must not change.

    @Note: Currently one interpolation object can only handle elements with
      all the same number of nodes and integration points. It's not possible
      to have tets and cohesives in one interpolation for instance.
    """

    # current format version:
    #  - grid: identifies the target points for the interpolation
    #  - mesh: node coords and element connectivity as bae.mesh_01.Mesh object
    #  - elemCoords: list of element numbers and coordinates
    # pickleFileFormat = "Interpolation Data Version 1.0: grid,mesh,elemCoords"
    pickleFileFormat = "Interpolation Data v2.0: binary grid,mesh,elemCoords"

    @staticmethod
    def _gridToFile(grid, pickleFile, pickleFileFormat):
        """Write grid data to a pickle file.
        """
        if isinstance(grid, MeshStructuredPointsRot):
            # STRUCTURED_POINTS_ROT
            # DIMENSIONS nx ny nz
            # ORIGIN x y z
            # SPACING sx sy sz
            # ROTVEC rx ry rz
            # ROTANGLE deg
            pickleFile.write("STRUCTURED_POINTS_ROT\n")
            pickleFile.write("DIMENSIONS %d %d %d\n" % tuple(grid.gridPtNb))
            pickleFile.write("ORIGIN %r %r %r\n" %tuple(grid.originBackRotated))
            pickleFile.write("SPACING %r %r %r\n" % tuple(grid.spacing))
            rotvec, deg = mat_toRodrigues(grid.rotmat)
            pickleFile.write("ROTVEC %r %r %r\n" % tuple(rotvec))
            pickleFile.write("ROTANGLE %r\n" % deg)
        elif isinstance(grid, MeshStructuredPoints):
            # STRUCTURED_POINTS
            # DIMENSIONS nx ny nz
            # ORIGIN x y z
            # SPACING sx sy sz
            pickleFile.write("STRUCTURED_POINTS\n")
            pickleFile.write("DIMENSIONS %d %d %d\n" % tuple(grid.gridPtNb))
            pickleFile.write("ORIGIN %r %r %r\n" % tuple(grid.origin))
            pickleFile.write("SPACING %r %r %r\n" % tuple(grid.spacing))
        elif isinstance(grid, MeshUnstructuredPoints):
            # UNSTRUCTURED_POINTS
            # POINTS n
            pickleFile.write("UNSTRUCTURED_POINTS\n")
            pickleFile.write("POINTS %d\n" % len(grid))
            a = array("f", (x for pt in grid for x in pt))
            a.tofile(pickleFile)

    @staticmethod
    def _gridFromFile(pickleFile, pickleFileFormat):
        """Read grid data and return corresponding MeshStructuredPoints-
        or MeshUnstructuredPoints- object.
        """
        if (pickleFileFormat
            =="Interpolation Data Version 1.0: grid,mesh,elemCoords"):
            return pickle.load(pickleFile)
        elif (pickleFileFormat
              =="Interpolation Data v2.0: binary grid,mesh,elemCoords"):
            gridType = pickleFile.readline().rstrip().upper()
            msg("Found grid of type %s..." % gridType, debugLevel=1)

            if gridType=="STRUCTURED_POINTS_ROT":
                # STRUCTURED_POINTS_ROT
                # DIMENSIONS nx ny nz
                # ORIGIN x y z
                # SPACING sx sy sz
                # ROTVEC rx ry rz
                # ROTANGLE deg
                gridPtNb=None; originBackRotated=None; spacing=None
                rotvec=None; deg=None
                for _ in range(5):
                    line = pickleFile.readline().rstrip().upper().split(None)
                    if line[0]=="DIMENSIONS":
                        gridPtNb = map(int, line[1:4])
                    elif line[0]=="ORIGIN":
                        originBackRotated = map(float, line[1:4])
                    elif line[0]=="SPACING":
                        spacing = map(float, line[1:4])
                    elif line[0]=="ROTVEC":
                        rotvec = map(float, line[1:4])
                    elif line[0]=="ROTANGLE":
                        deg = float(line[1])
                    else:
                        msg("WARNING: Unrecognized option %s for grid type %s."
                            % (line[0], gridType))
                if any(x is None for x in (
                        gridPtNb, originBackRotated, spacing, rotvec, deg)):
                    raise BadPickleFile(
                        "Options missing for grid type %s." % gridType)
                else:
                    return MeshStructuredPointsRot(
                        gridPtNb=gridPtNb,
                        originBackRotated=originBackRotated,
                        spacing=spacing,
                        rotVector=rotvec,
                        deg=deg,)

            elif gridType=="STRUCTURED_POINTS":
                # STRUCTURED_POINTS
                # DIMENSIONS nx ny nz
                # ORIGIN x y z
                # SPACING sx sy sz
                gridPtNb=None; origin=None; spacing=None
                for _ in range(3):
                    line = pickleFile.readline().rstrip().upper().split(None)
                    if line[0]=="DIMENSIONS":
                        gridPtNb = map(int, line[1:4])
                    elif line[0]=="ORIGIN":
                        origin = map(float, line[1:4])
                    elif line[0]=="SPACING":
                        spacing = map(float, line[1:4])
                    else:
                        msg("WARNING: Unrecognized option %s for grid type %s."
                            % (line[0], gridType))
                if (gridPtNb and origin and spacing):
                    return MeshStructuredPoints(
                        gridPtNb=gridPtNb, origin=origin, spacing=spacing)
                else:
                    raise BadPickleFile(
                        "Options missing for grid type %s." % gridType)

            elif gridType=="UNSTRUCTURED_POINTS":
                # UNSTRUCTURED_POINTS
                # POINTS n
                line = pickleFile.readline().rstrip().upper().split(None)
                if line[0] != "POINTS":
                    raise BadPickleFile(
                        "POINTS option missing for grid type %s." % gridType)
                else:
                    nbPoints = int(line[1])
                    arrPoints = array("f")
                    arrPoints.fromfile(pickleFile, nbPoints*3)
                    # CODE DESCRIPTION (speed beats readability here):
                    #  - for the izip(*[iter]*N)-phrase see ...
                    #    ... https://docs.python.org/2/library/...
                    #    ... itertools.html#itertools.izip
                    return MeshUnstructuredPoints(
                        list(x) for x in izip(*[iter(arrPoints)]*3))
            else:
                raise BadPickleFile("Unrecognized grid type %s." % gridType)
        else:
            raise BadPickleFile("Unrecognized pickleFileFormat %s"
                                % repr(pickleFileFormat[:80]))

    @staticmethod
    def _meshToFile(mesh, pickleFile, pickleFileFormat):
        if (pickleFileFormat
            =="Interpolation Data Version 1.0: grid,mesh,elemCoords"):
            pass
        elif (pickleFileFormat
              =="Interpolation Data v2.0: binary grid,mesh,elemCoords"):
            mesh.write(pickleFile)

    @staticmethod
    def _meshFromFile(pickleFile, pickleFileFormat):
        if (pickleFileFormat
            =="Interpolation Data Version 1.0: grid,mesh,elemCoords"):
            return pickle.load(pickleFile)
        elif (pickleFileFormat
              =="Interpolation Data v2.0: binary grid,mesh,elemCoords"):
            mesh = Mesh()
            mesh.read(pickleFile)
            return mesh
        else:
            raise BadPickleFile("Unrecognized pickleFileFormat %s"
                                % repr(pickleFileFormat[:80]))

    @staticmethod
    def _tetElemCoordsToFile(elemCoords, file_, pickleFileFormat):
        """service function to (space-, time-) efficiently write a list of
        (element number, element coordinates).

        Replacement for pickle.dump(elemCoords, file_).
        """
        # write nb of elemCoord-tuples (nb of points) and nb of coords per point
        nbCoords = set(
            len(c)
            for e, c in elemCoords
            if e is not None)
        if len(nbCoords)==1:
            # get that single length of all coords tuples
            nbCoords = nbCoords.pop()
            msg("_tetElemCoordsToFile: nbCoords=%d" % nbCoords, debugLevel=10)
        elif not nbCoords:
            nbCoords = 0
            msg("WARNING: Could not determine nb of coords per element."
                " len(elemCoords)=%d" % len(elemCoords))
        else:
            # coords tuples don't have all the same length
            raise NotImplementedError(
                "Element coords tuples don't have all the same length."
                " ...found: %s." % sorted(nbCoords))
            nbCoords = max(nbCoords)
            # when storing element coords in the file we'd have to fill up the
            # elemCoords-tuple with zeros to make them all of length nbCoords
            # not implemented, yet....

        a = array("I", (len(elemCoords), nbCoords))
        a.tofile(file_)

        # CODE DESCRIPTION (speed beats readability here):
        # The strange logical expressions (like "e or 0") have the effect of
        # conditional execution:
        # Points not in the mesh are represented by a (None, [])-item. In the
        # resulting array and binary file they are represented by element nb 0
        # and coordinates 0,0,0,0

        # write element numbers to file
        a = array("I", ((e or 0) for e,s in elemCoords))
        a.tofile(file_)

        # write element coords to file
        zeros = [0.0]*nbCoords
        a = array("f", (x
                        for e,s in elemCoords
                        for x in (e and s or zeros)))
        a.tofile(file_)

    @staticmethod
    def _tetElemCoordsFromFile(file_, pickleFileFormat):
        """service function to read the list of (element number, element
        coordinates) from file

        Replacement for pickle.load(file_).
        """
        arr = array("I")
        arr.fromfile(file_, 2)
        nbElems, nbCoords = arr
        msg("_tetElemCoordsFromFile: nbElems=%d, nbCoords=%d"
            % (nbElems, nbCoords), debugLevel=10)
        if nbCoords < 1 or nbCoords > 10:
            raise NotImplementedError(
                "Expected between one and ten element coordinates but got %d."
                % nbCoords)

        arrElem = array("I")
        arrElem.fromfile(file_, nbElems)
        msg("_tetElemCoordsFromFile: read %d element numbers." % len(arrElem),
            debugLevel=10)

        arrCoords = array("f")
        arrCoords.fromfile(file_, nbElems*nbCoords)
        msg("_tetElemCoordsFromFile: read %d element coords. First values:\n%s"
            % (len(arrCoords), arrCoords[:20]), debugLevel=10)

        # now assemble elemCoords as list of
        # (element number, element coordinates)-tuples
        # CODE DESCRIPTION (speed beats readability here):
        #  - for the izip(*[iter]*N)-phrase see ...
        #    ... https://docs.python.org/2/library/itertools.html#itertools.izip
        #  - the logical expressions have the effect of conditional execution:
        #    element nb 0 in the file results in (None, [])-item
        iterList = [arrElem,] + [iter(arrCoords)]*nbCoords
        elemCoords = [((x[0] or None), (x[0] and x[1:] or []))
                      for x in izip(*iterList)]
        msg("_tetElemCoordsFromFile: assembled elemCoords. First values:\n%s"
            % elemCoords[:20], debugLevel=10)

        return elemCoords

    def __len__(self):
        """
        Return the number of grid points in the interpolation:
         >>> interpolation = InterpolMeshToPoints(pickleName="mypic.interpol")
         >>> print len(interpolation)
        """
        try:
            return len(self.elemCoords)
        except AttributeError:
            return 0

    def getElemShapeFuncs(self):
        """Return an iterable of element shape function values for
        interpolation of nodal values. The items correspond to the
        points in self.grid and are tuples like this:
        (element number, [SF1, SF2, ..., SFN]). (E.g. for 10-node tets: N=10)

        >>> for elem, weights in interpolation.getElemShapeFuncs():
        >>>     print elem, weights   # generally a bad idea....
        ... 154023, [0.12,0.003,...,0.04]
        ... (something like that about 2mio times...)

        @Note: Always use this function to get those values instead of
          using the underlying instance-attributes directly. Later versions
          may change the implementation details!
        """
        return self.mesh.getElemShapeFuncs(self.elemCoords)

    def getElemIPInterpolFuncs(self):
        """Return an iterable of weighting factors for interpolation of
        values defined at integration points. The items correspond to
        the points in self.grid and are tuples like this:
        (element number, [N1, N2, ..., Nn]). (E.g. for 10-node tets: n=4)

        >>> for elem, weights in interpolation.getElemShapeFuncs():
        >>>     print elem, weights   # generally a bad idea....
        ... 154023, [0.24,0.03,0.32,0.41]
        ... (something like that about 2mio times...)

        @Note: Always use this function to get those values instead of
          using the underlying instance-attributes directly. Later versions
          may change the implementation details!
        """
        return self.mesh.getElemIPInterpolFuncs(self.elemCoords)

    def writeRemapParam(self, fileNameBase=None, notDefinedKey=None):
        """Calculates nodal weights and Gauss point weights and either writes
        them to two files or returns them as two strings.


        Output type
        ===========

        If fileNameBase is given then the function
        writes two files named fileNameBase+".nodalweights" and
        fileNameBase+".gaussptsweights". Those are the files that the
        initial version of remap requires for its interpolation.

        If fileNameBase is not given (None) then return a tuple of
        (nodalweights, gaussptsweights) --two strings that contain what would
        otherwise have been in the said files. In binary form.


        What to do with points outside the mesh
        =======================================

        If notDefinedKey is not given (None) then for each undefined point
        write existing dummy node or element numbers and set all weighting
        points to zero. This does not require special cleverness of the
        interpolation function. But it will always return zero for undefined
        values and is problematic for the material number field.

        If notDefinedKey is an integer than set node numbers and element
        numbers for undefined points to this value. It must not be a
        valid element or node number.


        @Note: remap also needs a grid data string: nx,ny,nz,x,y,z,sx,sy,sz.
        Create it like so (grid = self.grid):
         >>> # grid data string: nx,ny,nz,x,y,z,sx,sy,sz
         >>> gridNumbers = list()
         >>> gridNumbers.extend(grid.gridPtNb)
         >>> gridNumbers.extend(grid.origin)
         >>> gridNumbers.extend(grid.spacing)
         >>> remapGridParam =  ",".join(str(x) for x in gridNumbers)
         >>> del gridNumbers
        """

        msg("Preparing nodal and Gauss-pt weight data.")
        nbGridPts = len(self.grid)
        elNodes = self.mesh.elNodes  # just an abbreviation
        if not len(elNodes):
            raise ValueError(
                "writeRemapParam fails with no elements in the mesh.")

        # determine nb of nodes per element
        nodeNbs = set(len(elNodes[el])
                      for el, coords in self.elemCoords if el is not None)
        if len(nodeNbs) < 1:
            msg("WARNING: None of the requested output points lies inside an"
                " element of the mesh. Assuming 10 nodes on each element."
                "\nPLEASE CHECK!")
            nodesPerElem = 10
        elif len(nodeNbs) > 1:
            raise RuntimeError(
                "ERROR: We have different types of elements with different"
                " numbers of nodes (%s) in the interpolation."
                % sorted(nodeNbs))
        else:
            nodesPerElem = nodeNbs.pop()
            msg("Elements have %d nodes" % nodesPerElem, debugLevel=10)
        del nodeNbs

        # determine nb of GAUSS pts per element
        ipNbs = set(len(self.mesh.elShapeGaussPoints[self.mesh.elShape[el]])
                    for el, coords in self.elemCoords if el is not None)
        if len(ipNbs) < 1:
            msg("WARNING: None of the requested output points lies inside an"
                " element of the mesh. Assuming 4 GAUSS pts on each element."
                "\nPLEASE CHECK!")
            gaussptsPerElem = 4
        elif len(ipNbs) > 1:
            raise RuntimeError(
                "ERROR: We have different types of elements with different"
                " numbers of GAUSS pts (%s) in the interpolation."
                % sorted(ipNbs))
        else:
            gaussptsPerElem = ipNbs.pop()
            msg("Elements have %d Gauss points" % gaussptsPerElem,
                debugLevel=10)
        del ipNbs

        # binary data format
        nodalConv = Struct("if"*nodesPerElem)
        gaussptsConv = Struct("i%df" % gaussptsPerElem)

        # for points outside the mesh
        if notDefinedKey is None:
            dummyElem = iter(elNodes).next()   # just any valid element number
            dummyNode = elNodes[dummyElem][0]  # just any valid node number
            dummyNodeWeights = nodalConv.pack(*([dummyNode, 0.0]*nodesPerElem))
            dummyElemWeights = gaussptsConv.pack(
                dummyElem, *([0.0,]*gaussptsPerElem))
        else:
            dummyElem = notDefinedKey
            dummyNode = notDefinedKey
            dummyNodeWeights = nodalConv.pack(
                dummyNode, 1.0, *([dummyNode, 0.0]*(nodesPerElem-1)))
            dummyElemWeights = gaussptsConv.pack(
                dummyElem, 1.0, *([0.0,]*(gaussptsPerElem-1)))
        msg("notDefinedKey: %s, %d bytes dummyNodeWeights, %d bytes"
            " dummyElemWeights"
            % (notDefinedKey, len(dummyNodeWeights), len(dummyElemWeights)),
            debugLevel=10)

        # write ... .nodalweights
        if fileNameBase:
            outputName = "%s.nodalweights" % fileNameBase
            output = open(outputName, "wb")
        else:
            output = StringIO()
        ticker = MsgTicker(
            "Creating nodal weights data, grid pt %%d/%d" % len(self))
        for elem, weights in self.getElemShapeFuncs():
            if elem is None:
                output.write(dummyNodeWeights)
            else:
                ##### DEBUG
                # if not(isinstance(elem, (int, long))) or elem<=0:
                #     raise Exception("elem is <%s> %s" % (type(elem), elem))
                # if any(not(isinstance(w, float)) or w<-1 or w >1
                #        for w in weights):
                #     raise Exception("weights is %s" % weights)
                # if any(not(isinstance(n, int)) or n<=0
                #        for n in elNodes[elem]):
                #     raise Exception("elNodes[elem] is %s" % elNodes[elem])
                ##### DEBUG
                output.write(nodalConv.pack(
                    *[x for xx in izip(elNodes[elem], weights) for x in xx]))
            ticker.tick()
        del ticker
        fileSize = output.tell()
        msg("Created %d bytes of nodal weights data." % fileSize)
        if (fileSize != nbGridPts*nodesPerElem*8):
            msg("WARNING! Expected size of %d bytes but got %d."
                % (nbGridPts*nodesPerElem*8, fileSize))
        if fileNameBase:
            msg("Wrote nodal weights to %s" % outputName)
        else:
            nodalweights = output.getvalue()
        del output

        # write ... .gaussptsweights
        if fileNameBase:
            outputName = "%s.gaussptsweights" % fileNameBase
            output = open(outputName, "wb")
        else:
            output = StringIO()
        ticker = MsgTicker(
            "Creating GAUSS pts data, grid pt %%d/%d" % len(self))
        for elem, weights in self.getElemIPInterpolFuncs():
            if not elem:
                output.write(dummyElemWeights)
            else:
                output.write(gaussptsConv.pack(elem, *weights))
            ticker.tick()
        del ticker
        fileSize = output.tell()
        msg("Created %d bytes of GAUSS pts weights data." % fileSize)
        if (fileSize != nbGridPts*(1+gaussptsPerElem)*4):
            msg("WARNING! Expected size of %d bytes but got %d."
                % (nbGridPts*(1+gaussptsPerElem)*4, fileSize))
        if fileNameBase:
            msg("Wrote GAUSS pts weights to %s" % outputName)
            result = None
        else:
            gaussPtsweights = output.getvalue()
            result = (nodalweights, gaussPtsweights)
        del output
        return result


class InterpolMeshToPoints(InterpolMeshBaseClass):
    """For interpolating field values from an unstructured (FEM-) mesh to a
    regular points grid (L{MeshStructuredPoints} or L{MeshStructuredPointsRot})
    or just some arbitrary points (L{MeshUnstructuredPoints}) this class
    delivers interpolation data, in particular element and node numbers and
    weighting factors.

    When interpolating field values from a mesh to certain points/vertices
    you need to determine in which element the point sits and where exactly
    in the element it sits. The __init__ function of this class either
    calculates this data or loads it from a file generated previously and
    then stores it conveniently and efficiently for subsequent use if
    applicable.

    The actual interpolation is then done by Field.interpolate functions
    of the various flavours e.g.
    L{nodal<bae.field_01.FieldPosMeshNode.interpolate>} or
    L{integration point based<bae.field_01.FieldPosMeshElemIP.interpolate>}
    fields.

    Example:
     >>> points = [[0.1,0.4,5.0], [1.2,3.5,5.3], ...]
     >>> grid = MeshUnstructuredPoints(points)
     >>> if elems is None:
     >>>     getMeshCallBack = odb.getAbqModel()
     >>>     gmCallBackArgs = dict(recognizedBlocks="NODE,ELEMENT")
     >>> else:
     >>>     def getMeshCallBack():
     >>>         mesh = odb.getAbqModel(recognizedBlocks="NODE,ELEMENT,ELSET")
     >>>         elset = mesh.getUnionSet("elset", elems)
     >>>         mesh.initializePointSearch(
     >>>             box=grid.getBoundingBox(), elems=elset)
     >>>         return mesh
     >>>     gmCallBackArgs = dict()
     >>> interpolation = InterpolMeshToPoints(
     >>>     grid=grid,
     >>>     getMeshCallBack=getMeshCallBack,
     >>>     gmCallBackArgs=gmCallBackArgs,
     >>>     pickleName="MyProject_G01_checkPoints"
     >>>     )
     >>> field = odb.getFieldFromOdb(fieldName, odbFrame)
     >>> FieldClass = Field.classFromPosType(
     >>>     "%s_Pts"%fieldName, 'structPt', field.dataType)
     >>> ptValues = FieldClass(field.interpolate(interpolation))

    @ivar grid: A MeshStructuredPoints or MeshStructuredPointsRot or
      MeshUnstructuredPoints object identifying the output points for the
      interpolation.
    @ivar mesh: A (tet-) Mesh object. Contains all node and element
      definitions involved in the interpolation.
    @ivar elemCoords: A list of tuples (element number, element coordinates)
      for the points in self.grid.

      elemCoords may contain tuples with None as element number, as
      Mesh.getElemCoords() provides for points outside the mesh.
    @cvar pickleFileFormat: This string is stored as first item in the
      pickle file to check compatibility between the pickle file and the
      version of the InterpolMeshToPoints class currently in use. This
      string has to be modified if the format of the pickle file is changed.
      Otherwise it must not change.

    @Note: Currently one interpolation object can only handle elements with
      all the same number of nodes and integration points. It's not possible
      to have tets and cohesives in one interpolation for instance.
    """

    def __init__(self, grid=None,
                 mesh=None, getMeshCallBack=None, gmCallBackArgs={},
                 pickleName=None):
        """Initialize an interpolation of a regular grid of output points to
        a tet mesh. Calculate the element
        coordinates for each output grid point and possibly store them for
        subsequent use in a pickle file. Or load this data from the so created
        pickle file.

        Specify a pickleName argument for speed up on subsequent runs.

        There are two costly operations involved in the initialization of an
        interpolation: getting the mesh data (nodes and elements) and finding
        element coordinates for all grid points.

        There are different situations requiring a corresponding way to
        initialize this interpolation:

         1. The mesh (node coords and element connectivity) is not yet
            present in memory, there might or might not already be a pickle
            file.
            So the mesh argument can't be supplied as argument to __init__().
            In this case specify grid, getMeshCallBack and pickleName. If the
            pickle file does exist then it is read and the grid parameters
            given by the grid argument are checked against the contents of
            the pickle file. Mesh and interpolation data are being taken from
            the pickle file, getMeshCallBack() is not called in this case. If
            the pickle file does not exist yet, then getMeshCallBack() will
            be called, the element coordinates for the grid points will be
            calculated and the pickle file will be created for subsequent
            runs.
         2. The mesh data is already present in memory and will be supplied
            to __init__() as mesh argument. There might or might not already
            be a pickle file.
            In this case specify grid, mesh and pickleName. If the pickle
            file does exist then it is read and the grid parameters given by
            the grid argument are checked against the contents of the pickle
            file. Mesh and interpolation data are being taken from the pickle
            file, the mesh argument is being ignored in this case. If the
            pickle file does not exist yet, then the element coordinates for
            the grid points will be calculated and the pickle file will be
            created for subsequent runs.
         3. If you know that the pickle file is present you can skip the
            mesh or getMeshCallBack arguments. You may specify the grid
            argument to have it checked against the stored parameter, but
            this is optional. This is useful for the coupling procedure where
            a different script (prepareCouplingGrid.py) does the
            initialization.
         4. If you don't care about speeding up subsequent calls you may leave
            out the pickleName argument. In this case a pickle file is neither
            being read nor being written. You have to specify grid and mesh or
            getMeshCallBack

        @param grid: A L{MeshStructuredPoints} or L{MeshStructuredPointsRot}
          or L{MeshUnstructuredPoints} object identifying the output points
          for the interpolation. This is optional if pickleName is specified
          and the grid can be loaded from the corresponding file.

        @param mesh: A (tet-) L{Mesh} object. Ignored if it can be loaded from
          the specified pickleName. If it's expensive to create consider
          to specify getMeshCallBack and gmCallBackArgs instead for deferred
          creation.

        @param getMeshCallBack: This function will be called if a mesh object
          is not specified and can't be loaded from the pickle file. In case
          the creation of the mesh object is expensive this deferred creation
          may save you from doing it.

        @param gmCallBackArgs: For deferred creation of the mesh: A dict
          containing kwargs to be supplied to getMeshCallBack.

        @param pickleName: base name for a file containing the data of this
          interpolation. It will be saved and used on subsequent invocations.
          If you don't specify this argument the interpolation data will
          be calculated new and no pickle file will be generated for subsequent
          invocations.

          ".interpol" will be added as file name extension. The file name
          should identify the tet mesh and the set of points. This is to ensure
          that certain interpolation data is not accidently used for the wrong
          mesh or points.
          E.g. "RWAY2013_G01_gridBox1_h10.interpol" for RWAY mesh G01 and a
          point set labeled "gridBox1_h10".

          If the mesh or the requested points have changed then remove this
          file.

        @Note: Always use keyword arguments. Might add a single positional
          argument later. This would then be used as initializer. E.g. for
          copy constructor.
        """

        # see if we can find a pickle-file
        pickleFile = None
        if pickleName:
            pickleName = "%s.interpol" % pickleName
            try:
                pickleFile = open(pickleName, "rb")
            except IOError:
                msg("Did not find interpolation pickle file %s." % pickleName)

        # actually read the interpolation data from the pickle file
        if (pickleFile):
            msg("Reading interpolation pickle file %s." % pickleName)

            try:
                pickleFileFormat = pickleFile.readline().rstrip()
                msg("Found pickleFileFormat %s." % repr(pickleFileFormat),
                    debugLevel=1)

                # self.grid
                self.grid = self._gridFromFile(pickleFile, pickleFileFormat)

                # verify that the stored points grid is identical to the saved
                if grid and self.grid != grid:
                    text = ("ERROR: Specified grid data does not conform to"
                            " the grid data stored in %s." % pickleName)
                    msg(text)
                    # additional diagnostic information
                    # trigger the MeshStructuredPoints.__cmp__() method
                    log.setDebugLevel(10)  # to get all messages of __cmp__()
                    (self.grid == grid)
                    msg("Possible reasons:")
                    msg(" *** Not changing the name of the box but some of its"
                        " parameters. So we're at the wrong .interpol file.")
                    msg(" *** .interpol file has been created by an older"
                        " version of the bae-module. If point intervals don't"
                        " fit exactly into the box the upper corner of the box"
                        " has to be adjusted. The corresponding  procedure"
                        " has been changed in Oct 2015. Now we round the number"
                        " of intervals to the nearest integer. You may adjust"
                        " the gridbox to match what's stored in the .interpol"
                        " file.")
                    raise ValueError(text)
                msg("Interpolation points grid contains %d points."
                    % len(self.grid))

                # self.mesh
                self.mesh = self._meshFromFile(pickleFile, pickleFileFormat)
                msg("Interpolation mesh part contains %d elements and %d nodes."
                    % (len(self.mesh.elNodes), len(self.mesh.nodeCoords)))

                # self.elemCoords
                self.elemCoords = self._tetElemCoordsFromFile(
                    pickleFile, pickleFileFormat)
                msg("Loaded point positions for %d points."
                    % len(self.elemCoords))
                assert len(self.elemCoords)==len(self.grid)

                pickleFile.close()

            except BadPickleFile, exc:
                pickleFile.close()   # just for safety...
                pickleFile = None
                msg("Could not read interpolation pickle file %s:" % pickleName)
                msg(exc.args[0])
                import os
                bakPickleName = pickleName+".old"
                os.rename(pickleName, bakPickleName)
                msg("Renamed the old %s to %s. Will ignore it from now on."
                    % (pickleName, bakPickleName))

        # actually calculate position of the grid points in the mesh
        # (either we did not find a pickle file at all or problems occured
        # like unrecognized format or unexpected data in the file...)
        if not(pickleFile):
            if grid:
                self.grid = grid
            else:
                # this is an error, grid must be specified at this stage!
                if pickleName:
                    text = ("Grid not specified and file %s can't be loaded"
                            % pickleName)
                else:
                    text = "Neither grid nor pickleName have been specified."
                raise ValueError(text)

            if not mesh:
                msg("Getting the mesh.")
                mesh = getMeshCallBack(**gmCallBackArgs)
                msg("Got mesh: %s" % mesh)

            # list of (element number, elem coord)-tuples for each pt in
            # grid
            self.elemCoords = mesh.getElemCoordsForGrid(grid)
            msg("Calculated element coordinates for %d points, %d points are"
                " within the mesh." % (len(grid), len([
                    None for elem, coords in self.elemCoords
                    if elem is not None])))

            # mesh that only contains elements and nodes needed for the
            # interpolation
            allElems = set(elem for elem, coords in self.elemCoords)
            # we only want nodes and elements, so make sure we use the method
            # of mesh_01.Mesh
            self.mesh = Mesh.partMeshFromElements(mesh, allElems)
            msg("Created a submesh with %d elements and %d nodes."
                % (len(self.mesh.elNodes), len(self.mesh.nodeCoords)))
            del mesh

        # write position of the grid points in the mesh to pickle file
        if pickleName and (pickleFile is None):
            msg("Writing interpolation data pickle file %s ..."
                % pickleName)

            # check output directory
            import os
            pickleFileDir = os.path.dirname(pickleName)
            if (pickleFileDir and pickleFileDir != "."
                and not os.path.isdir(pickleFileDir)):
                msg("Directory %s does not yet exist. Creating it now."
                    % pickleFileDir)
                os.makedirs(pickleFileDir)

            pickleFile = open(pickleName, "wb")
            pickleFile.write("%s\n" % self.pickleFileFormat)
            self._gridToFile(self.grid, pickleFile, self.pickleFileFormat)
            msg("Saved grid info. Writing submesh...")
            self._meshToFile(self.mesh, pickleFile, self.pickleFileFormat)
            msg("Saved submesh. Writing point positions...")
            self._tetElemCoordsToFile(
                self.elemCoords, pickleFile, self.pickleFileFormat)
            msg("Saved point positions for %d points. Closing file ..."
                % len(self.elemCoords))
            del pickleFile
            msg("Finished writing mesh information pickle file %s."
                % pickleName)


class InterpolMeshToCentroids(InterpolMeshBaseClass):
    """Interpolation data for retrieving values at element centroids.

    For interpolating field values on an unstructured (FEM-) mesh to
    element centroids this class delivers interpolation data -weight factors-
    in a way consistent with the other interpolation classes.

    The actual interpolation is then done by Field.interpolate functions
    of the various flavours e.g.
    L{nodal<bae.field_01.FieldPosMeshNode.interpolate>} or
    L{integration point based<bae.field_01.FieldPosMeshElemIP.interpolate>}
    fields.

    @ivar grid: A MeshUnstructuredPoints object identifying the output points
      --the centroids' coordinates-- for the interpolation.
    @ivar mesh: A (tet-) Mesh object. Contains all node and element
      definitions involved in the interpolation.
    @ivar elemCoords: A list of tuples (element number, element coordinates)
      for the points in self.grid.
    @ivar elemList: list of element numbers of the points in grid.
      The order of entries in self.grid and self.elemList matches.

    @Note: Currently one interpolation object can only handle elements with
      all the same number of nodes and integration points. It's not possible
      to have tets and cohesives in one interpolation for instance.
    """

    def __init__(self, model, elems, boundingBox=None):
        """
        @param model: A L{bae.abq_model_02.Model} or a L{bae.mesh_01.Mesh}
           object.
        @param elems: Set of element labels.
           If model is a bae.abq_model_02.Model object then elems can also
           be an elset name or a list of elset names: anything that
           L{bae.abq_model_02.Model.getUnionSet}() accepts as items argument.
        @param boundingBox: L{BoundingBox<bae.misc_01.BoundingBox>}-object (or
           [[xmin,ymin,zmin],[xmax,ymax,zmax]]). Consider only elements whose
           centroid lies within this box. If None then all elements will be
           considered.

        @raises RuntimeError: if no elements are left after filtering.
        """

        msg("Calculating element centroids and weight factors.")
        if not elems:  # empty sets, strings or lists are ignored as well
            centroids = Mesh.getElCentroids(model)
            shapeElems = sorted(
                (shape, els)
                for shape, els in model.shapeEl.iteritems())
        else:
            if (isinstance(elems, set)
                    and isinstance(iter(elems).next(), int)):
                elset = elems
            elif hasattr(model, "getUnionSet"):
                elset = model.getUnionSet("elset", elems)
            else:
                raise ValueError(
                    "InterpolMeshToCentroids.__init__: model (%s) and elems"
                    " (%s) arguments don't fit." % (type(model), type(elems)))
            if not elset:
                raise ValueError(
                    "InterpolMeshToCentroids.__init__: elems argument doesn't"
                    " yield any elements. (Wrong elset name?)\n...elems: %s"
                    % str(elems)[:100])
            centroids = Mesh.getElCentroids(model, elset)
            # sorted list of (shape, list of elements)-tuples
            # elements only from elems-parameter
            shapeElems = sorted(
                (shape, els)
                for shape, els in (
                    (sh2, el2.intersection(elset))
                    for sh2, el2 in model.shapeEl.iteritems())
                if els)

        # filter boundingBox
        if boundingBox:
            centroids = dict(
                (el, xyz)
                for el, xyz in centroids.iteritems()
                if  boundingBox[0][0]<xyz[0]<boundingBox[1][0]
                and boundingBox[0][1]<xyz[1]<boundingBox[1][1]
                and boundingBox[0][2]<xyz[2]<boundingBox[1][2])
            shapeElems = [
                (shape, sorted(els.intersection(centroids)))
                for shape, els in shapeElems if els]
        else:
            # just fix the order of elements
            shapeElems = [(shape, sorted(els)) for shape, els in shapeElems]

        # check that we still have elements to process
        if not shapeElems:
            raise RuntimeError("ERROR: no elements left after filtering in"
                               " InterpolMeshToCentroids.__init__().")

        # list of element numbers
        self.elemList = [e for shape, els in shapeElems for e in els]
        self.grid = MeshUnstructuredPoints(centroids[e] for e in self.elemList)
        self.mesh = model
        self.elemCoords = [
            (elem, coords)
            for coords, ee in (
                    (model.shapeCentrCoords[shape], eee)
                    for shape, eee in shapeElems)
            for elem in ee]
        msg("Finished calculating element centroids and weight factors"
            " for %d elements." % len(self.elemCoords))
        return


class InterpolMeshToFaceCentroids(InterpolMeshBaseClass):
    """Interpolation data for retrieving values at element face centroids.

    For interpolating field values on an unstructured (FEM-) mesh to
    element face centroids this class delivers interpolation data, in
    particular element and node numbers and weighting factors in a way
    consistent with the other interpolation classes.

    The actual interpolation is then done by Field.interpolate functions
    of the various flavours e.g.
    L{nodal<bae.field_01.FieldPosMeshNode.interpolate>} or
    L{integration point based<bae.field_01.FieldPosMeshElemIP.interpolate>}
    fields.

    @ivar grid: A L{MeshUnstructuredPoints} object identifying the output points
      --the face centroids' coordinates-- for the interpolation.
    @ivar mesh: A (tet-) L{Mesh} object. Contains all node and element
      definitions involved in the interpolation.
    @ivar elemCoords: A list of tuples (element number, element coordinates)
      for the points in self.grid.
    @ivar elemFaceIdList: A list of (element number, faceId) tuples identifying
      the element and face corresponding to each output point.

    @Note: Currently one interpolation object can only handle elements with
      all the same number of nodes and integration points. It's not possible
      to have tets and cohesives in one interpolation for instance.
    """

    def __init__(self, *args, **kwargs):
        """
        @kwarg model: a L{bae.abq_model_02.Model}; requires faceEl argument,
           precludes surface argument. DEPRECATED, better supply surface.
        @kwarg faceEl: a "faceEl"-dictionary:
           {"S1": set((4346, 4768, 6821, ...)), ...}; requires model argument,
           precludes surface argument. DEPRECATED, better supply surface.
        @kwarg surface: a L{bae.surface_03.ElemFaceSurface} object; precludes
           model and faceEl arguments.
        @Note: model and faceEl arguments are provided for compatibility
           reasons only.
        """

        if (args and hasattr(args[0], "nodeCoords")
            and hasattr(args[1], "__getitem__")):
            # arguments are model and faceEl
            model = args[0]
            faceEl = args[1]
        elif args and hasattr(args[0], "mesh") and hasattr(args[0], "faceEl"):
            # argument is a ElemFaceSurface object
            surface = args[0]
            model = surface.mesh
            faceEl = surface.faceEl
        elif "surface" in kwargs:
            surface = kwargs["surface"]
            model = surface.mesh
            faceEl = surface.faceEl
        elif "model" in kwargs and "faceEl" in kwargs:
            model = kwargs["model"]
            faceEl = kwargs["faceEl"]
        else:
            raise ValueError(
                "Could not gobble up proper arguments for"
                " InterpolMeshToFaceCentroids. args: %s; kwargs: %s."
                % (map(type, args), sorted(kwargs)))

        # can't import at the module level because of cyclic dependence,
        # i.e. bae.fields imports from mesh_01 as well
        from bae.field_01 import FieldTypeVector

        # group elements ("elems") according to face-elsets ("S1", "S2", ...)
        # and element shape
        msg("Grouping according to element shape and face.", debugLevel=1)
        faceShapeElems = [
            [faceId2, shape2, elems]
            for (faceId2, shape2, elems) in (
                    (faceId, shape, sorted(faceEl[faceId]
                                           .intersection(elemsOfShape)))
                    for faceId in sorted(faceEl)
                    for shape, elemsOfShape in model.shapeEl.iteritems())
            if elems]
        msg("InterpolMeshToFaceCentroids: %d face-shape-groups"
            % len(faceShapeElems), debugLevel=10)

        # calculate face centroids
        msg("Calculating face centroids.", debugLevel=1)
        self.grid = MeshUnstructuredPoints()
        for faceId, shape, elems in faceShapeElems:
            faceCentroidWeights = model.elShapeFaceCentrWeights[shape][faceId]
            self.grid.extend(
                FieldTypeVector.weightedSum(
                    (model.nodeCoords[node] for node in model.elNodes[el]),
                    faceCentroidWeights)
                for el in elems)
        msg("InterpolMeshToFaceCentroids: %d face centroids" % len(self.grid),
            debugLevel=10)

        # store element coordinates for face centroids
        msg("Calculating face centroid weight factors.", debugLevel=1)
        self.elemCoords = list()
        self.elemFaceIdList = list()
        for faceId, shape, elems in faceShapeElems:
            faceCentrCoords = model.elShapeFaceCentrCoords[shape][faceId]
            self.elemCoords.extend((el, faceCentrCoords) for el in elems)
            self.elemFaceIdList.extend((el, faceId) for el in elems)
        msg("InterpolMeshToFaceCentroids: %d elemCoords, %d elemFaceIds"
            % (len(self.elemCoords), len(self.elemFaceIdList)), debugLevel=10)

        # list of element numbers
        self.mesh = model
        msg("Finished calculating face centroids and weight factors"
            " for %d faces." % len(self.elemCoords), debugLevel=1)
        return


class InterpolMeshFacesToSurface(object):
    """For interpolating field values from an L{ElemFaceSurface} in a
    unstructured (FEM-) L{mesh<Mesh>} to a L{TriangleSurface}-mesh this class
    delivers interpolation data.

    The actual interpolation is then done by Field.interpolate functions
    of the various flavours e.g.
    L{nodal<bae.field_01.FieldPosMeshNode.interpolate>} or
    L{integration point based<bae.field_01.FieldPosMeshElemIP.interpolate>}
    fields.

    @ivar sourceTopo: A L{ElemFaceSurface} object. The L{Field} object that
      uses this InterpolMeshFacesToSurface object for its interpolate method
      must be defined on sourceTopo.mesh.
    @ivar targetTopo: A L{TriangleSurface} object. The L{Field}.interpolate
      method that uses this InterpolMeshFacesToSurface object must return a
      L{Field} object that is defined on this (triangle-) mesh.
    @ivar volRef: A reference from the targetTopo elements to the original
      element and face in sourceTopo: a dict {tri elem label: (vol elem label,
      faceId, Gauss-pt-interpol-data)}.
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
    """

    def __init__(self, sourceTopo):
        """
        @param sourceTopo: A L{ElemFaceSurface} object. The L{Field} object
           that uses this InterpolMeshFacesToSurface object for its interpolate
           method must be defined on sourceTopo.mesh.
        """
        self.sourceTopo = sourceTopo
        (self.targetTopo, self.volRef) = sourceTopo.getTriangleSurface(
            quadTetTo4Tris=True, returnVolRef=True)
