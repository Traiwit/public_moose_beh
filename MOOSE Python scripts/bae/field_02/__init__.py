# -*- coding: utf-8 -*-
"""
Structure to handle time varying field data, i.e. scalars, vectors, tensors
(and other) that change in space and time.

The central object is a L{FieldsCollection<bae.field_02.FieldsCollection>}
which has
 - topo: e.g. a L{PointCloud}, L{StructuredGrid}, ...
 - fields: a dict {field name: ndarray}; ndarray.shape is (N,M, ...),
   N = len(topo), M = nb of time points, and '...' can be: nothing for a scalar
   field, 3 for a vector and 3,3 for a tensor field. '...'==6 can also be
   possible for a vector notation of a symmetrical tensor.
   Higher order tensores should be possible too, but consider memory footprint.
   fields are NxM(x...), Nx1(x...), or 1xM(x...)
 - times: a dict { field name: ndarray } of space-less persistent fields.
   ndarrays in times are 1xM(x...)
 - space: a dict { field name: ndarray } of time-less persistent fields.
   ndarrays in space are Nx1(x...)

Note that the dimensions of items in fields, times, space and topo must match:
Each field must have the same dimension N in space or 1 (spaceless field, constant
in space). And each field must have the same dimension M in time or 1 (timeless
field, constant in time, e.g. material code).

The 'persistent' fields in self.times and self.space are always copied over
to the result when a sub-field is created by self.__getitem__. Those fields
are meant to characterize the space resp. time points of the fields.

A FieldsCollection object is equally suitable (should be) for the special /
corner cases:
 - When there's only one field: See Examples "Create and use scalar, vector and
   tensor fields" below. See the L{FieldsCollection.array} attribute.
 - For only one time step (timeless field)
 - A time series without space dimension (spaceless field)


Examples:
=========

Working with FieldsCollection
-----------------------------

A FieldsCollection container is essentially a dict {field name: field / numpy
array} with common topo, space and times attributes for all fields.

Initializing values on creation:
Note that the dimensions of field values, topo, space and times must match:
Each field must have the same dimension N in space or 1 (spaceless field,
constant in space, i.e. a time series). And each field must have the same
dimension M in time or 1 (timeless field, constant in time, e.g. material
code).

 >>> flds = FieldsCollection(
 ...     topo=PointCloud(points),
 ...     times=frames,  # field names can not be called "topo" or "times"
 ...     U = fieldA,    # a numpy array as component 'U'
 ...     T = np.random.randn(nPts, nTimes, 3, 3),  # a numpy array as 'T'
 ...     )

You can add fields to the FieldsCollection and remove fields from it. E.g.
adding single fields:
 >>> for ii in range(3):
 >>>     name = 'V_%d' % ii
 >>>     flds[name] = np.random.randn(nPts, nTimes)

Adding contents of a dictionary or other FieldsCollection:
 >>> flds.update({'R', np.random.randn(nPts, nTimes)})

Removing fields:
 >>> del flds['R']

Accessing single fields. FieldsCollection behaves like a dict of numpy-arrays
if given a single string as key (or index if you prefer).
 >>> for fieldName in flds:
 >>>     print '%s: max=%f' % (fieldName, flds[fieldName].max())

You can pass a list of strings to get a new FieldsCollection object with just
those fields:
 >>> subF = flds[['U', 'T',]]
 >>> print "T_max=", subF["T"].max()

Numerical or array indexes return a new FieldsCollection object with all fields
sliced accordingly. Note: the second index is the time index.
 >>> refFlds = flds[:,0]  # get only first time point
 >>> U_rel = flds["U"] - refFlds["U"]

Selecting fields with a list of field names and array-slicing can be combined.
Creating subsets: The first index identifies the field by name.
 >>> subF = flds[['U', 'T',], 3:10, 4]

file I/O:  (PVF not yet implemented)
 >>> csvPath = 'testField.csv'
 >>> flds.writePvfCsv(csvPath)
 >>> ff = FieldsCollection().readPvfCsv(csvPath)

reshaping Field-objects stored in the FieldsCollection
 >>> fldS = flds.toScalars('U')
 >>> fldV = flds.toVector(['U_1', 'U_2', 'U_3'], 'V')

A FieldsCollection can be created from scratch and populated later:
 >>> new = FieldsCollection()
 >>> new.topo = flds.topo
 >>> new["Vx"] = np.zeros((len(flds.topo),2))
 >>> new["Vx"][:,0] = flds["V_1"][:,0]
 >>> new["Vx"][:,1] = flds["V_1"][:,1]
 >>> # store time frame names
 >>> new.times["names"] = np.array(["Y2010_M02", "Y2010_M03", ])


Create and use scalar, vector and tensor fields
-----------------------------------------------

Create a scalar field on N points at M time points:
 >>> from bae.field_02 import PointCloud, FieldsCollection
 >>> sFieldA = FieldsCollection(
 >>>     topo=PointCloud([[0.0,0.0,0.0], [0.2,1.5,0.7], ...N values...]),
 >>>     times=[(2,i) for i in range(M)],
 >>>     A=np.array([[1.0, 2.0, 1.5,...N values...],
 ...                 [...],
 ...                  ...M lists...]))

Vector math operations on vector fields:
 >>> npArr = np.norm(vFieldU["U"])  # numpy array
 >>> np.log(npArr, out=npArr)

...or if vFieldU contains only one field:
 >>> npArr = np.norm(vFieldU.array)  # numpy array

Results of math operations can be stored in the fields collection:
 >>> sFieldA["logU"] = npArr       # store as field "logU"

Combining components or extracting scalar values from vectors  (not yet implemented)
 >>> plunge, bearing = vFieldU.vectorToAngles("U") # returns scalars
 >>> u1, u2, u3 = vFieldA.toScalars()           #  components

Addressing single values and slicing:
 >>> aa = sFieldA[2:23]          # subset of some points
 >>> bb = sFieldA[:, 3:5]        # subset of some times
 >>> cc = vFieldU[:, :, 1]       # get a component, will return a scalar field
 >>> dd = vFieldU[2:80, :-2, 1]  # arbitrary index expressions
 >>> ee = sFieldA[sFieldA<0.0]   # boolean array indexing is also possible


Processing VTK files
--------------------
 >>> from bae.field_02 import FieldsCollection
 >>> flds = FieldsCollection.fromVtk("mydata.vtk")
 >>> u_1 = flds["U_1"]  # ... a numpy array, shape Nx1
 >>> u_2 = flds["U_2"]
 >>> u_hor = np.sqrt(u_1**2 + u_2**2)
 >>> flds["u_hor"] = u_hor   # add new field to container, topo from container
 >>> flds.toVtk("result.vtk")  # export as vtk
 >>> new = FieldsCollection()
 >>> new.topo = flds.topo
 >>> new["u_hor"] = u_hor
 >>> new.toVtk("singleField.vtk")  # export as vtk


Create a random VectorField using numpy.random
----------------------------------------------
 >>> nPts = 100
 >>> nTimes  = 7
 >>>
 >>> x = np.random.randn(nPts)
 >>> y = np.random.randn(nPts) + 10
 >>> z = np.random.randn(nPts) - 10
 >>> coords = np.vstack((x, y, z)).T
 >>> val = np.random.randn(nPts, nTimes, 3)
 >>>
 >>> # create the VectorField
 >>> from bae.field_02 import PointCloud, FieldsCollection
 >>>
 >>> # create time frames
 >>> frameIds = np.array([(2,i) for i in range(nTimes)])
 >>>
 >>> # create a topo object, some unstructured points
 >>> points = PointCloud(coords)
 >>>
 >>> # create a vector field object, field name "A"
 >>> vFieldA = FieldsCollection(topo=points, A=val)
 >>> vFieldA.times["frameIds"] = frameIds
"""

__version__ = "0.7"

_version_history_ = """
0.1 TR New
0.2 TR completely reworked - e.g. splitted fields into submodules
    fields+topo+frames
0.3 TR once more completely reworked
    - fields.Fields --> fieldscollection.FieldsCollection
    - simplified StaticEntities and their subclassing (PointCloud,
      StructuredGrid)
    - simplified init for topo-classes
    - removed some 'getIndexFrom...'-fuctions from topo-classes
    - some bug-fixes
    - added very premature version for UnstructuredMesh-topo
    find previous version in:
    mySvnRespos/python/tags/fieldData_00_preUnstructuredMesh
0.4 TR added fixed attribute and addFixed function
0.5 TR fixed multidim indexing of field objects, added iterVariables and
    iterFixed to StaticEntities
0.6 GP renamed to field_02, but still in alpha stage, not usable!
0.7 GP incompatible change: flds["T"] delivers np-array
"""

__todo__ = """
- topo-classes
  . mandatory method: getSubTopo taking a serial index or slice or the like
    for selecting from the space dimension of a FieldsCollection. Service
    function for FieldsCollection.__getitem__()
    (other name options discarded: getSubset, getPart)
  . optional additional __getitem__ for specialized purposes.
    returns the same type (always?) but only selected points,
  . getPoints() returns a nx3 array of coords
  . __len__
  . StructuredGrid not sparse
  . add class StructuredGridSparse (or other name?),
    refer to old StructurderGrid implementation
  . both structured grids with optional rotation matrix
  . new mixin in plot.py: FieldsCollectionPlots. Add to base classes of
    FieldsCollection
  . see topo/unstructuredmesh/__init__.py

. add explanation somewhere, why we don't follow the following suggestion:
  Because single indexes, index ranges, index lists are treated inconsistently
  in numpy and we prefer consistency. Plus: what if length of space index
  doesn't match length of time index?
. suggest change: make __getitem__([i0:iN,j0:jN]) return [i0,j0], [i1,j1], ...
  rather than the current set [i0:iN,:][:,j0:jN]
  . purpose: select a path in time and space
  . this is a rare use case but so is selecting space and time at the same time
  . semantics to get [i0:iN,:][:,j0:jN] by exactly this expression is simple
    and straightforward. -> only small benefit to do it in one operation by
    the expression [i0:iN,j0:jN]
  . on the flip side it's almost impossible to select a space-time-path if
    the expression [i0:iN,j0:jN] is meant to return [i0:iN,:][:,j0:jN]

- ??? add FieldsCollectionCore.addTimes(fld, persistent=False) and
  addSpace(...) (?)
  At least if arg persistent==False: fld (on input) has only space *or* time
  dimension and will now be reshaped N -> Nx1 or 1xN
- supply reshape functions: Abaqus-Voigt-vector to 3x3 tensor in
  fieldscollection
- from fieldData_00.Fields copy vectorToAngles to FieldsCollectionCore

- ??? checkShape(newShape) function (to be called from FC.__setitem__() et al. )
  checks newShape-argument against self.getShape()

"""

future_ideas = """
fieldsCollectionTimeIter class: iterator over time steps that holds data
for only one time step at a time.
- For transforming / converting data from the odb and other sources step by
  step.
- Pro: less danger of memory overflow, smaller memory footprint, less damage
  on interrupt (if the procedure crashes at step 30 then step 1 to 29 are
  already done)
- __init__ takes list of time steps (and topo?)
- stores current time step and can provide further data for the current time
  step, like frame name, Abaqus Transfer Number, ...
- different sub-classes with different data sources:
  . general / base class: next method returns a fieldsCollection
  . from odb: reads the data from successive odb frames
  . from vtk; accordingly from a set of vtk files
  . from csvs
  . from FieldsCollection, yielding all or some of its time steps. This can
    be used for time shifts (for dU3 or U_ref).
  . from several fieldsCollectionTimeIter, "merger".
  . from another fieldsCollectionTimeIter and composer functions
    (should this and the previous go in one class? i.e. simple merger be a
    special "composer")
- have output classes that take those iterators and create corresponding
  output, see vtk_02.FramesFieldsOnUnstructuredPoints. Maybe a subclass
  of fieldsCollectionTimeIter?
- what if topo changes over time? Consider subsidence on evolving pit surface.

possible restructuring for field_03: subclassing ndarray:
Derive a single fields object from ndarray
 - no name stored, needs to be done externally, can be an optional attribute
 - have a module global dict for topos. topo.__getitem__ selects partitions of
   topo and stores the new partition and stores the relation old topo --
   operation -- new topo.
 - possibly also store the space and times dicts for each topo-times
   configuration
 - possibly this is way too complicated, if we can't do it simpler it's
   rubbish.
"""


from .fieldscollection import FieldsCollectionCore
from .fileIO.readerWriter import FieldsCollectionReaderWriter
# from .plot import FieldsCollectionPlots
from .interpolation import FieldsCollectionInterpolation


class FieldsCollection(
        FieldsCollectionCore,
        FieldsCollectionReaderWriter,
        # FieldsCollectionPlots,
        FieldsCollectionInterpolation,
        ):
    pass


from .topo.pointcloud import PointCloud
from .topo.structuredgrid import StructuredGrid
from .topo.unstructuredmesh.meshbase import \
    HomoMesh, HeteroMesh, \
    MeshNodes, MeshElems, MeshGaussPts
from .topo.unstructuredmesh.trimesh import TriMesh

# from .plot import (
#     hist, scatter3D, plotOverFrames, quiver3D,
#     plungeBearingScatter, plungeBearingDensity,)
