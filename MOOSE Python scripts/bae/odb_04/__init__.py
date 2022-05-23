"""Package to access the Abaqus odb from Python.

Delivers L{OdbReader}, L{OdbReaderInterpolateToPoints} classes.

The package is intended to be imported like a module:
 >>> from bae.odb_04 import OdbReaderInterpolateToPoints

Usage example: simple odbToVtkGridPoints script:
 >>> from bae.odb_04 import OdbReaderInterpolateToPoints
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

Usage example: read some integration point values from the odb:
 >>> from bae.odb_04 import OdbReader
 >>>
 >>> # open the odb
 >>> odb = OdbReader(odbPath, version=abaqusVersion)
 >>> odb.openOdb()
 >>>
 >>> # get mesh data from the odb as a L{Model<bae.abq_model_02.Model>} object
 >>> mesh = odb.getAbqModel()
 >>>
 >>> # limit the following request to the specified elset
 >>> elsetName = odb.createFilterElset(elementLabels)
 >>>
 >>> # get field data from the odb (of type L{bae.field_01.Field}
 >>> fld_S1 = getFieldFromOdb("Step-2", 86, "S_MIN_PRINCIPAL")
 >>>
 >>> # clean-up and print results
 >>> del odb  # also closes the odb
 >>> print "element\tvalues for all integration points"
 >>> for elem in elementLabels:
 >>>     print "%d\t%s" % (elem, fld_S1[elem])


"""

__version__ = "4.01"
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
4.01 GP new version based on v3.14
    - OdbReaderInterpolateToPoints now *contains* an OdbReader-object
      following the "Composition over inheritance" principle
    - changed postData object from dict to Container
      (postData.sequence["frameNames"] ==> postData.sequence.frameNames)
    - odb steps are addressed by name not by number anymore.
"""

# programmer's note on the design of the odb_04 package:
#  - appears like a module
#  - extension module for (most) communication tasks with the Abaqus process
#  - different python modules for different classes possible

todo = """
- incompatible change for odb_04 or be3.odb_01: move interpolation into
  field_02 package. => interpolator unstructmesh to structuredgrid, to
  arbitrary points, to elemFaceSurface -> trianglesurf

- for be3.odb_01:
  . provide OdbReader.readFields(frameList, fieldNames) returning a
    FieldsCollection object. What to do with material assignment field?
  . ... And getTopo returning an UnstructuredMesh. To be called by readFields
    as well. Silently store this for later requests? Possibly...
"""

import os
import cPickle as pickle
from bae.log_01 import msg
from bae.misc_01 import Container

#------------------------------------------------------------------------------

defaultNotDefinedKey = -1
defaultNotDefinedValue = {
    "scalar": 0.0,
    "vector": [0.0, 0.0, 0.0],
    "tensor": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
}


def getOdbPostData(odbPath):
    """Get the postData-object associated with the specified odb.
    It is being fed from the corresponding ..._postData database / directory.

    Strip the ".odb" file name extension from the path of the corresponding
    odb and append "_postData" to get the name of the postData-directory.
    Each python file (with a ".py" file name extension) in this directory will
    result in a corresponding attribute of the postData object.

    Example: Suppose alongside "MyModel_R01_G01_Q02_M02.odb" you have the
    folder "MyModel_R01_G01_Q02_M02" containing the file "sequence.py".

    Content of MyModel_R01_G01_Q02_M02/sequence.py:
     >>> frameNames = {1: ["Start", "Y2012_Q4", "Y2013_Q1"]}
     >>> seqElsetsExcav = {1: [(),('Pit_2012_02','Dev_2012_02'),(Dev_2012_03,)]}

    Then we can access this data like so:
     >>> postData = getOdbPostData("MyModel_R01_G01_Q02_M02.odb")
     >>> frameNames = postData.sequence.frameNames
     >>> print frameNames[1][1:]
     ... ["Y2012_Q4", "Y2013_Q1"]
     >>> seqSets = postData.sequence.seqElsetsExcav
     >>> print seqSets[1][1:]
     ... [ ('Pit_2012_02','Dev_2012_02'),(Dev_2012_03,) ]

    Old-style postData in a single file XXX_postData.py alongside the odb
    will be recognized as well


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

    @param odbPath: file name of the odb (complete path) with or without the
    ".odb" file name extension.

    @returns: postData object of type L{Container}. Typical attributes:
      sequence, sequence.frameNames, sequence.seqElsetsExcav
      material, material.matNumberList, material.matProps
      varNames, varNames.varNameDamage, varNames.varNameState, ...

    @raises IOError: if the XXX_postData.py file can not be accessed. XXX is
    the name of the odb without the a trailing ".odb"-suffix.
    """

    # ignore possible odb file extension
    if odbPath[-4:].lower() == ".odb":
        odbPath = odbPath[:-4]

    # initialize result
    postData = Container()

    # try new format version: postData directory
    postDataDirName = odbPath+"_postData"
    try:
        postDataFiles = os.listdir(postDataDirName)
    except OSError:
        postDataFiles = None

    if postDataFiles is not None:
        postData.version = "1.01"
        for fn in postDataFiles:
            attrName, ext = os.path.splitext(fn)
            fullpath = os.path.join(postDataDirName, fn)
            if ext==".py":
                postDataDict = dict()
                execfile(fullpath, postDataDict)
                data = Container(**postDataDict)
            elif ext==".pickle":
                data = pickle.load(open(fullpath, "rb"))
            else:
                # silently ignoring other file types
                msg("When reading postData ignoring %s/%s"
                    % (postDataFiles, fn), debugLevel=1)
                continue
            if hasattr(postData, attrName):
                raise ValueError(
                    "Redundant post-data: When trying to read %s there already"
                    " is a %s attribute loaded earlier!"
                    % (fullpath, attrName))
            setattr(postData, attrName, data)

        # done, finished reading all files in postData directory
        return postData

    # try old format version "0.01": postData.py file -> postData.sequence...
    postDataFileName = odbPath+"_postData.py"
    try:
        pdf = open(postDataFileName, "r")
    except IOError:
        pdf = None

    if pdf is not None:
        # load/parse postData from external file
        postDataDict = dict()
        code = compile(pdf.read(), postDataFileName, "exec")
        exec(code, postDataDict)
        postData = Container()
        postData.version = "0.01"
        postData.sequence = Container(**postDataDict)

    # in any case return a container (be it empty...)
    return postData


from .abaqusVersion import getAbqVersionStr, getAbqExe
from .odbReader import OdbReader
from .odbReaderInterpolate import OdbReaderInterpolateToPoints
