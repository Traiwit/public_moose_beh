"""Package to access the Abaqus odb from Python.

Delivers L{OdbReader}, L{OdbReaderInterpolateToPoints} classes.

The package is intended to be imported like a module:
 >>> from bae.odb_03 import OdbReaderInterpolateToPoints

Usage example: simple odbToVtkGridPoints script:
 >>> from bae.odb_03 import OdbReaderInterpolateToPoints
 >>> from bae.mesh_01 import getMesh
 >>> from bae.vtk_02 import FieldsOnStructPointGrid
 >>>
 >>> # open the odb
 >>> odb = OdbReaderInterpolateToPoints(odbPath, version=abaqusVersion)
 >>>
 >>> # initialize the output grid and interpolation
 >>> grid = getMesh(**gridData)
 >>> odb.initOutputPoints(interpolPickleName=interpolPickleName, **gridData)
 >>>
 >>> frameList = [(2,i) for i in range(48, 76)]
 >>> fieldNames = ["SDV3", "S_MIN_PRINCIPAL", "U_3"]
 >>> for stepNb, frameNb in frameList:
 >>>     output = FieldsOnStructPointGrid(mesh=grid)
 >>>     outputFileName = "output_F%d%03d.vtk" % (stepNb, frameNb)
 >>>     for fieldName in fieldNames:
 >>>         msg("Reading field %s, step %d, frame %d."
 >>>             % (fieldName, stepNb, frameNb))
 >>>         field = odb.getFieldFromOdb(stepNb, frameNb, fieldName)
 >>>         output.updateFields(field)
 >>>     output.toVtk(
 >>>         outputFileName,
 >>>         description=("Abq-Grid-Data step %d, frame %d"
 >>>                      % (stepNb, frameNb)))
"""

__version__ = "3.14"
_version_history_ = r"""
Versions:

2.0 GP added version number
3.01 changed from module to package (folder with different module-files) 
3.02 changed: odbReader_ext.OdbReader.storeRemapWeights() requires nb of grid
    points as first argument
3.03 GP added: postData property;
    fixed: initOutputPoints() can be called multiple times for different grids
    on the same odb.
3.04 GP changed: make initOutputPoints more flexible on the grid argument,
    make odb_03 module import the version from odbReader
3.05 GP added getOdbPostData
3.06 GP added: try to load default binary extension module if one for the
    specified Abaqus version does not exist.
3.07 GP changed from multiprocessing.SyncManager to async_mananger-module
3.08 GP new clean design with odbAccessServer: separate Abq-C-Api executable
3.09 GP added initOutputAtCentroids using InterpolMeshToCentroids
        and initOutputAtFaceCentroids using InterpolMeshToFaceCentroids
3.10 GP added interpolationType argument for OdbReaderInterpolateToPoints
     ... .getFieldFromOdb() for piecewise constant interpolation.
     E.g. for the status field.
3.11 GP added setFilterElset, setFilterNset
3.12 GP changed initOutputAtFaceCentroids arguments, prefered now:
    ElemFaceSurface
3.13 GP changed OdbReader.getMatFieldFromOdb possibly modifies a given
    matNameToMatNumber dict in place, converting all mat names to upper case.
    Concentrated versioning information here (removed it from odbReader.py)
3.14 GP added odbOpen, odbClose functions to odbReader. You can now create
    an OdbReaderInterpolateToPoints-instance without opening the odb to be
    able to already prepare the interpolation.
"""

# programmer's note on the design of the odb_03 package:
#  - appears like a module
#  - extension module for (most) communication tasks with the Abaqus process
#  - different python modules for different classes possible

todo = """
- getFrameNbs (or similar) ... read number of frames for each odb step from
  the odb.

- incompatible change for odb_04:
  . rename odbPath argument of OdbReaderInterpolateToPoints to odb.
  . incorporate OdbReaderBase into OdbReader. I don't see the reason for this
    split anymore. Have one module odbReader with OdbReader class and another
    "interpolate" or so with the OdbReaderInterpolateToPoints class.

- incompatible change for odb_04 or be3.odb_01: move interpolation into
  field_02 package. => interpolator unstructmesh to structuredgrid, to
  arbitrary points, to elemFaceSurface -> trianglesurf

- for be3.odb_01:
  . provide OdbReader.readFields(frameList, fieldNames) returning a
    FieldsCollection object. What to do with material assignment field?
  . ... And getTopo returning an UnstructuredMesh. To be called by readFields
    as well. Silently store this for later requests? Possibly...
"""

#------------------------------------------------------------------------------

defaultNotDefinedKey = -1
defaultNotDefinedValue = {
    "scalar": 0.0,
    "vector": [0.0, 0.0, 0.0],
    "tensor": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
}

def getOdbPostData(odbPath):
    """Get the postData-dictionary associated with the specified odb.
    This dictionary is being fed from the corresponding ..._postData.py
    database.

    Usually contains some seqElsets... dictionaries like seqElsetsExcav in
    the exampe below. And a frameNames dictionary.

    The data will be read each time this function is called. It might be
    most efficient to store the data for subsequent access.

    Example: Suppose in frame 3 of step 2 labelled "YR2012_M02" the elsets
    "Pit_2012_02" and "Dev_2012_02" will be excavated and in frame 4
    "YR2012_M03" only "Dev_2012_03":
     >>> postData = getOdbPostData("/datapth/myodb.odb")
     >>> frameNames = postData["frameNames"]
     >>> print frameNames[2][3:5]
     ... ['YR2012_M02','YR2012_M03']
     >>> seqSets = odb.postData["seqElsetsExcav"]
     >>> print seqSets[2][3:5]
     ... [ ('Pit_2012_02','Dev_2012_02'),(Dev_2012_03,) ]

    @param odbPath: file name of the odb (complete path)

    @raises IOError: if the XXX_postData.py file can not be accessed. XXX is
    the name of the odb without the a trailing ".odb"-suffix.
    """

    # deduct postData-filename from odb name
    if odbPath[-4:].lower() == ".odb":
        postDataFileName = odbPath[:-4]+"_postData.py"
    else:
        postDataFileName = odbPath+"_postData.py"

    # load/parse postData from external file
    postData = dict()
    execfile(postDataFileName, postData)
    return postData

from abaqusVersion import getAbqVersionStr, getAbqExe
from odbReader import OdbReader, OdbReaderInterpolateToPoints
