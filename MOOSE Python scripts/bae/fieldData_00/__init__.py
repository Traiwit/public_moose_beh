# -*- coding: utf-8 -*-
"""
Structure to handle time varying field data, i.e. scalars, vectors, tensors
(and other) that change in space and time.

General structure of a Field object (L{ScalarField}, L{VectorField},
or L{TensorField})::
 Field
    |- name
    |- topo     <- indexable Topo-object of size nPoints
    |- frames   <- indexable Frame-object of size nFrames
    |- values   <- numpyArray of shape:
                         scalar field: (nPoints, nFrames)
                         vector field: (nPoints, nFrames, 3)
                         tensor field: (nPoints, nFrames, 3, 3)

Structure of L{FieldsCollection} (subclass of dict)::
 FieldsCollection
    |- topo     <- topo and frames attributes are referenced by all stored...
    |- frames   <- ...fields and the FieldsCollection object
    |- keys     <- field names
    |- values   <- Field objects


Examples:
=========

Create and use Scalar-, Vector- and TensorField(s)
--------------------------------------------------

Create a VectorFieldObject:
 >>> from bae.fieldData_00 import PointCloud, Frames, VectorField
 >>> vFieldA = VectorField(
 >>>     name='A',
 >>>     values=np.array([1.0, 2.0, 1.5,...]),
 >>>     topo=PointCloud([[0.0,0.0,0.0], [0.2,1.5,0.7], ...]),
 >>>     frames=Frames([(2,i) for i in range(7)]))

Vector math operations on the values attribute:
 >>> npArr = vFieldA.values  # numpy array
 >>> np.log(npArr, out=npArr)

Some math operations can be applied directly on the field object:
 >>> sFieldA = vFieldA.norm()                   # sFieldA is a ScalarField
 >>> plunge, bearing = vFieldA.vectorToAngles() # returns ScalarFields
 >>> u1, u2, u3 = vFieldA.toScalars()           # ScalarFields of components

Addressing single values and slicing:
 >>> aa = vFieldA[2:23]          # subset of some points
 >>> bb = vFieldA[:, 3:5]        # subset of some Frames
 >>> cc = vFieldA[:, :, 1]       # get a component, will return a ScalarField
 >>> dd = vFieldA[2:80, :-2, 1]  # arbitrary index expressions
 >>> ee = vFieldA[sFieldA<0.0]   # boolean array indexing is also possible

Index finding (see L{frames} and L{topo}-Objects for more methods)
 >>> iFr = vFieldA.frames.idxFromIds([(2,3), (2,5)])
 >>> subFr = vFieldA[:,iFr]
 >>> iPts = vField.topo.idxFromPath( vField.topo.boundingBox() )
 >>> diagPts = vField[iPts]

Working with FieldsCollection
-----------------------------

A FieldsCollection container is essentially a dict {field name: L{Field} object}
with a common topology and time points attribute for all fields.

Initializing values on creation:
Note that the dimensions of field values, topo and frame must match.
If Field objects are being used on initialization (or added later) then the
topo objects must be equal.
 >>> flds = FieldsCollection(
 ...     topo=points,
 ...     frames=frames,  # field names can not be called "topo" or "frames"
 ...     U = fieldA,  # a Field-object as component 'U'
 ...     T = np.random.randn(nPts, nFrs, 3, 3),  # a numpy array as 'T'
 ...     )

You can add fields to the FieldsCollection and remove fields from it. E.g.
adding single fields:
 >>> for ii in range(3):
 >>>     name = 'V_%d' % ii
 >>>     flds[name] = np.random.randn(nPts, nFrs)

Adding contents of a dictionary or other FieldsCollection:
 >>> flds.update({'R', np.random.randn(nPts, nFrs)})

Removing fields:
 >>> del flds['R']

Accessing single fields. (FieldsCollection is derived from dict.)
 >>> for fieldName in sorted(flds):
 >>>     field = fld[fieldName]
 >>>     print '%s: max=%f' % (fieldName, field.values.max())

Creating subsets: The first index identifies the field by name.
Note that f[['U',]] returns a new L{FieldsCollection} object while f['U']
returns the stored Field itself:
 >>> subF = fld[['U', 'T',], 3:10, 4]

file I/O:
 >>> csvPath = 'testField.csv'
 >>> flds.writePvfCsv(csvPath)
 >>> ff = FieldsCollection().readPvfCsv(csvPath)

reshaping Field-objects stored in the FieldsCollection
 >>> flds.reshape.toScalars('T', inplace=True)
 >>> flds.reshape.toVector(['V_1', 'V_2', 'V_3'], 'V')

Processing VTK files
--------------------
 >>> from bae.fieldData_00 import FieldsCollection, Frames
 >>> flds = FieldsCollection.fromVtks("mydata.vtk")
 >>> u_1 = flds["U_1"].values  # ... a numpy array, shape Nx1
 >>> u_2 = flds["U_2"].values
 >>> u_hor = np.sqrt(u_1**2 + u_2**2)
 >>> flds["u_hor"] = u_hor   # add new field to container, topo from container
 >>> flds.toVtks("result.vtk")  # export as vtk
 >>> new = FieldsCollection()
 >>> new.topo = flds.topo
 >>> new.frames = Frames(names=["Y2010_M02"])
 >>> new["u_hor"] = u_hor
 >>> new.toVtks("singleField.vtk")  # export as vtk

Create a random VectorField using numpy.random
----------------------------------------------
 >>> nPts = 100
 >>> nFrs  = 7
 >>>
 >>> x = np.random.randn(nPts)
 >>> y = np.random.randn(nPts) + 10
 >>> z = np.random.randn(nPts) - 10
 >>> coords = np.vstack((x, y, z)).T
 >>> val = np.random.randn(nPts, nFrs, 3)
 >>>
 >>> # create the VectorField
 >>> from bae.fieldData_00 import PointCloud, Frames, VectorField
 >>>
 >>> # create a FrameObject
 >>> frameIds = [(2,i) for i in range(nFrs)]
 >>> frames = Frames(ids=frameIds)
 >>>
 >>> # create a topo object, some unstructured points
 >>> points = PointCloud(coords)
 >>>
 >>> # create a VectorFieldObject
 >>> vFieldA = VectorField(name='A', values=val, topo=points, frames=frames)
"""

__version__ = "0.5"

_version_history_ = """
0.1 TR New
0.2 TR completely reworked - e.g. splitted fields into submodules fields+topo+frames
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
"""

__todo__ = """
- rename fieldData_00 to fields_02?

- rename Frames-Class and frames-Attribute to Time/time?
- rename topo-Attributes to space?

- reshape functions Voigt-vector to 3x3 tensor in fieldscollection
"""

from fields import (
    ScalarField, VectorField, TensorField,)

from fieldscollection import FieldsCollection

from topo import (PointCloud, StructuredGrid)

from frames import (Frames, )

###ok: I'd suggest to have a submodule fieldData_01.plots
from commonfieldplots import (
    hist, scatter3D, plotOverFrames, quiver3D,
    plungeBearingScatter, plungeBearingDensity,)
