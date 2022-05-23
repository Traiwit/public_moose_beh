"""Read interpolated field data from the specified odb and export in different
file formats. Exports regular grid data to legacy vtk files.


Usage:
======

For points on a regular grid:
 >>> from bae.utils.odbToPointsData_02 import odbToVtkGridBox
 >>> odbToVtkGridBox(...)

Or for unstructured points
 >>> from bae.utils.odbToPointsData_02 import odbToCsvPoints
 >>> odbToCsvPoints(...)
"""

__version__ = "2.35"

_version_history_ = """
1.0 GP: made it a module in bae/utils, derived from:
    3_post/pointsDataFromOdb/odbToVTKGridBox.py version 2.0
1.1 GP: split off getRemapGridParamString()
1.2 GP incompatible change: renamed module
1.3 GP renamed odbToVTKGridBox_wMaxFlds to odbToVTKGridBox. It can now run
  without fieldNamesMaximize options. Added relative fields.
1.4 GP incompatibel change: renamed getMeshCallback to getMeshCallBack
       added matRegionsOption
1.5 GP compileRemap(): don't copy abaqus_v6.env,
       added odbElCentroidsToCsvPoints (from odbCohsvToCsvPoints.py v1.02)
       pipe output from subprocesses to log channel(s)
1.6 AF,GP: made the material name filter optional (and switched off by default)
1.7 GP added rotated grid option
1.8 GP added skipExisting option: skip re-creating existing result files
1.9 GP incompatible change: renamed odbToVTKGridBox to odbToVtkGridBox
2.0 GP switched to odb_03
       incompatible change: if fieldNamesMaximize is given then the first frame
       in frameList is **the last frame to be ignored**.
       dropped skipExisting option.
2.1 GP changed order of arguments, added odbToCSVPoints
2.2 GP added outputFileNameComposer and interpolationDataDir arguments to
       odbToVtkGridBox
2.3 GP incompatible changes to odbToVtkGridBox:
        - removed outputFileNameComposer argument and added frameIdToName
          instead.
        - renamed argument boxName to gridName
       added relativeToReferenceVtkPath to odbToVtkGridBox
       added interpolationDataDir argument to odbToCSVPoints
2.4 GP incompatible change: renamed odbToCSVPoints to odbToCsvPoints
2.5 GP added "CsvPerFieldRowFrColPtName" as outputType argument to
       odbToCsvPoints
2.6 GP added odbToDriveClosureCsv
2.7 GP added list of frames to vtk file description file
2.8 GP incompatible change: renamed frameIdToName argument to stepFrameToName
       to avoid name clash with vtk_02 functions
2.9 GP added matNamesList arguments to odbToVtkGridBox and odbToCsvPoints
2.10 GP added fieldMaximizeAfter argument to odbToVtkGridBox and odbToCsvPoints
       and outputFields to odbToVtkGridBox only
2.11 GP simplify renaming fields with outputFields argument
2.12 GP added Converter_PlasticStrainTensor class for plastic strain tensor
2.13 GP added odbToDriveClosureLinesCsv as a variant of odbToDriveClosureCsv
2.14 GP added per point reference frames to odbToDriveClosureLinesCsv,
        incompatible interface change
2.15 GP added stepFrameToName accepts "fromPostData"
2.16 TR calling getMatFieldFromOdb with noMatCode=0: odbToVtkGridBox and
        odbToCSVPoints will set material(index) to 0 when called with
        matRegionsOption = True --> easier indexing for postprocessing
2.17 GP added relativeToFrame="frameNameFromPostData" option to
        odbToDriveClosureLinesCsv
2.18 GP added outputFields to odbToCsvPoints
2.19 GP major restructuring: added OdbToPointsData and derived classes
2.20 GP modify fieldnames argument: items can be
        (odbfieldname, interpolationtype)-tuples
2.21 GP added odbFaceCentroidsToCsvPoints
        and odbVariableFaceCentroidsToCsvPoints
2.22 GP added getExtraFields and corrected Converter_CombinedVtk
2.23 GP added Converter_StandardVtk, Converter_StressEigenDecomp
2.24 GP incompatible change for the getExtraFields argument: getExtraFields
        is now required to take an arbitrary number of keyword arguments.
2.25 GP add getU_1 to Converter_StandardVtk so that it's got the same interface
        as Converter_CombinedVtk, create the interpolationDataDir if necessary
2.26 GP added Converter_CombinedVtk_FS4
2.27 GP changed Converter_CombinedVtk similar to the FS4-variant.
2.28 GP added read FS4/Cavesim results from zip archive with same name as vtk
2.29 GP changed odbToVtkGridBox: removed default for relativeToFrame.
        Specifying both relativeToFrame and relativeToReferenceVtkPath is now
        considered an error.
2.30 GP changed Converter_CombinedVtk_FS4 vtk state: points in pit get
        status=cave if FS4 has particles in the cell
2.31 TR added Converter_StandardAndFOS
2.32 GP added Converter_CombinedVtk_FS4_multiBox
2.33 TR Converter_StandardAndFOS now can handle (var)fill material and anisotropy
        added mergeGetExtraFields to Converter (base class)
2.34 GP added skipExisting and raiseMissingData options to odbToVtkGridBox
2.35 GP corrected logP in converter-class
"""

todo = """
- getSurfaceAir option for odbVariableFaceCentroidsToCsvPoints

- add converter function for eigenvalue decomposition
- add x, y, z, name as possible output fields (constant fields)
  for odbToCsvPoints
- use OdbToPointsData-class for odbToDriveClosureLinesCsv and streamline the
  latter.
- add the layerToSetname-functionality in 3_Post/legacy/driveClosure/...
  ...closureLines_variant/getDriveClosure.py to odbToDriveClosureLinesCsv
- calculate fos for fill as well (would require to get varfill data from
  from vusubs_param.py and to assign varfilled elsets, sequence information
  can be derived from / toggled by state-field directly)

ideas for version ..._03:
=========================

1. composer function with "restype" attr

>>> def f(x):
...    print x
>>> f.restype = "scalar"

>>> f(4)
4
>>> f.restype
'scalar'

1a. even better with composer decorator:
def composer(dataType, fieldName=None):
    '''
    @param dataType:
    @param fieldName: optional
    '''
    def decorator(func):
        func.dataType = dataType
        func.fieldName = fieldName
        return func
    return decorator

@composer("scalar")
def composerLogP(cls, dat, grid):
    return cls(((x>0 and log10(x*1000+1)) or 0.0) for x in dat["SDV1"])

1b. maybe better: don't pass the cls argument to the composer and let the
    composer return an iterator (or iterable, like a list). The actual field
    object would then be created by the framework and the class would not be
    exposed to the user and the composer function. Example:
@composer("scalar")
def composerLogP(dat, grid):
    return (((x>0 and log10(x*1000+1)) or 0.0) for x in dat["SDV1"])

or equivalent (for more complicated cases):
@composer("scalar")
def composerLogP(dat, grid):
    for x in dat["SDV1"]:
        if x>0:
            yield log10(x*1000+1)
        else:
            yield 0.0


2. outputFields argument be a list of items of optional type:

2a. just string: read field of that name from odb and write to results
2b. (outName, odbFldName)-tuple of two strings: read odbFldName from odb and
    just rename to outName
2c. (None, odbFldName): read odbFldName from odb without direct output for later
    processing in composer functions, see below.
2d. (outName, composer-func): as in odbToPointsData_02: call composer function
    that creates a field of name outName
2e. (None, composer-func with fieldName): generate a field by means of the
    composer function. The field will be named according to the composer-
    attribute fieldName. This field is then accessible for other composer
    functions.
    That's to make multi-step composer functions more convenient. E.g. the
    eigenvalue calculation can be done with such an entry and then the result
    fields pull their data from this field created by this entry. Another
    use case would be data from external files (e.g. Cavesim results).

optional third True/False field, True means const field. A missing second item
in the tuple (option 2a) can be set to None.

maximize and relative-to-frame (or -to-external-file) as additional optional
items...


3. Material option automatic by requesting a field "Material" from the odb

3b. add automatic constant fields only if requested by listing it as odbFldName,
    i.e. an outputFields item "ptCoords", ("xyz", "ptCoords") or
    (None, "ptCoords"), the latter if needed in composer functions:
    ptCoords, ptName (also if outputType requires it), elsetName (only for
    odbElCentroidsToCsvPoints), element-nb.
    Replaces the structDataFieldNames-argument of legacy
    odbElCentroidsToCsvPoints()

4. getExtraFields:
  - should not return an iterator but store the fields by itself. I.e. get
    self (the OdbToVtkGRidBox object) and call its storeField method as
    often as needed. Because it's easier to have multiple getExtraFields
    functions from base classes or so. And calling it can just throw an
    exception easily differntiable from a possible error while storing...
  - this option is connected to the new suggests for outputFields.
    Maybe getEtraFields is not needed anymore or outputFields must be modified
"""

import os
import csv
import re
from math import log10
from itertools import izip, islice
from collections import defaultdict
import zipfile
import time
import traceback

from bae.future_01 import OrderedDict
from bae.misc_01 import RowIterFromCsvFile, CheckFileFinished
from bae.mesh_01 import getMesh, MeshUnstructuredPoints, Mesh
from bae.field_01 import createFieldClass, createFieldObject
from bae.vecmath_01 import vector, vector_plus, vector_sum, vector_scale, \
    vector_modif_scale, centre, mat_multvec, dot, norm, length
from bae.odb_03 import OdbReaderInterpolateToPoints, getOdbPostData
from bae.vtk_02 import FieldsOnStructPointGrid, \
    FramesFieldsOnUnstructuredPoints
from bae.surface_03 import ElemFaceSurface
from bae.makeHistory_03 import sequenceElsetsIter
from bae.log_01 import msg, MsgTicker

from bae.material_01 import MatPropsCollection
from bae.material_01.material import FactorOfSafetyIso, FactorOfSafetyAniso

try:
    import numpy as np
except ImportError:
    np = None


#------------------------------------------------------------------------------
# OdbToPointsData classes
#------------------------------------------------------------------------------
class OdbToPointsData(object):
    """Base class for L{odbToVtkGridBox}, L{odbToCsvPoints} and the like.

    Subclass this class supplying and overloading some functions then
    instantiate the new class and call it's run() method.

    These functions must or can be defined in the subclass:
     - L{__init__}: process arguments, invokes self.L{storeArgs}
     - L{initOutput}: prepare everything for the output of results
     - L{storeConstField}: store constant field
     - L{initRefData}: (optional) get reference data for relative fields,
       update self.refData.
     - L{processConstantFields}: (optional) process constant fields before we
       come to the frames-loop.
     - L{initOutputForFrame}: (optional) prepare for new frame data before
       looping over fields. Intended to initialize output and will only be
       called for frames for which output is requested.
     - L{storeField}: (optional) store the given field for output
     - L{storeFrame}: (optional) store frame data with all fields to output
     - L{closeOutput}: (optional) finalize output

    Usage:
     >>> class MyOdbToSomethingClass(OdbToPointsData):
     >>>     def __init__(self, odbPath, fieldNames, ...lots of arguments ...):
     >>>         ... some checks ...
     >>>         self.storeArgs(locals(), "odbPath", "fieldNames", ...)
     >>>     ... more overloaded functions
     >>> def myOdbToSomething(odbPath, fieldNames, ...lots of arguments ...):
     >>>     reader = MyOdbToSomethingClass(odbPath, fieldNames, ...args...)
     >>>     reader.run()
     >>>
     >>> myOdbToSomething(config.odbPath, ["U_1", "U_2"], ....)

    @ivar fieldNames: list of field name strings
    @ivar fieldInterpolType: {field name: type}, type being "default" or
       "const" for piecewise constant
    """
    def __init__(
            self, fieldNames,
            fieldNamesMaximize=None, fieldMaximizeAfter=None,
            relativeFieldNames=None, relativeToFrame=(2,0),
            ):
        """General and very basic processing of arguments:
        Some checks and storing as instance variables.

        Processes the fieldNames list with optional (name, interpol-type)-tuple
        items. Thereby creates the list of field names self.fieldNames and
        the dict self.fieldInterpolType {field name: interpolation type}

        Checks (but does not store!) the fieldNamesMaximize arguments:
        items must also be listed in the fieldNames arguments.
        Same for relativeFieldNames.

        Checks that fieldMaximizeAfter is a tuple of two ints. (Also not stored
        here!)
        """

        # for checking create a set of all fieldname - strings
        self.fieldNames = list()
        self.fieldInterpolType = dict()
        for fieldName in fieldNames:
            if isinstance(fieldName, (list, tuple)):
                self.fieldNames.append(fieldName[0])
                self.fieldInterpolType[fieldName[0]] = fieldName[1]
            else:
                self.fieldNames.append(fieldName)
                self.fieldInterpolType[fieldName] = "default"
        fieldNamesSet = set(self.fieldNames)
        msg("OdbToPointsData.__init__: got fieldNames: %s, %d times non-default"
            " interpolation-type."
            % (self.fieldNames,
               sum(1 for x in self.fieldInterpolType.itervalues()
                   if x != "default")), debugLevel=10)

        # do all checks, collect errors in valueErrors-list
        # ...and only report at the end
        valueErrors = list()

        # check arguments fieldNamesMaximize and fieldNames
        if fieldNamesMaximize:
            missingFields = [fn for fn in fieldNamesMaximize
                             if fn not in fieldNamesSet]
            if missingFields:
                valueErrors.append(
                    "ERROR: %d fields in argument fieldNamesMaximize are not"
                    " contained in argument fieldNames. That does not make"
                    " sense." % (len(missingFields)))
            del missingFields

            # check fieldMaximizeAfter argument
            if fieldMaximizeAfter:
                try:
                    if not(isinstance(fieldMaximizeAfter[0], int)
                           and isinstance(fieldMaximizeAfter[1], int)):
                        raise TypeError()
                except (TypeError, IndexError):
                    valueErrors.append(
                        "Incorrect value for fieldMaximizeAfter: %s"
                        % fieldMaximizeAfter)

        # check argument relativeFieldNames
        if relativeFieldNames:
            missingFields = [fn for fn in relativeFieldNames
                             if fn not in fieldNamesSet]
            if missingFields:
                valueErrors.append(
                    "ERROR: %d fields in argument relativeFieldNames are not"
                    " contained in argument fieldNames. That does not make"
                    " sense." % (len(missingFields)))
            del missingFields

        if valueErrors:
            raise ValueError("\n".join(valueErrors))
        return

    # service functions
    def storeArgs(self, locals_, *args):
        """Stores arguments listed as arguments in instance variables.

        @Note: Don't store fieldNames if you already called
        OdbToPointsData.__init__() before. The __init__-function prepares
        self.fieldNames for you already which you would overwrite with
        self.storeArgs(..., "fieldName", ...)
        """
        for arg in args:
            setattr(self, arg, locals_[arg])

    def prepareStepFrameToName(self):
        """Adapt the stepFrameToName argument in order to be able to use it
        for the frameIdToName argument of
        L{FramesFieldsOnUnstructuredPoints.initOutput}().
        """
        # convert stepFrameToName function:
        # The frameIdToName argument to vtk_02.FramesFieldsOnUnstructuredPoints.
        # ....initOutput --if its a function-- takes a single frameId argument:
        # a (stepNb, frameNb) - tuple.
        # On the other hand all bae.utils/odbToPointsData - methods
        # expect for their stepFrameToName argument a function taking two
        # arguments: stepNb, frameNb
        if callable(self.stepFrameToName):
            # strings don't need conversion, but functions do
            stepFrameToNameTwoArgs = self.stepFrameToName
            def stepFrameToName(frameId):
                return stepFrameToNameTwoArgs(*frameId)

        # if argument stepFrameToName=="fromPostData": define special function
        elif self.stepFrameToName=="fromPostData":
            # special case, frameNames from odb.postData
            if "frameNames" not in self.odb.postData:
                raise ValueError(
                    "frameNames not defined in associated _postData.py file for"
                    " odbPath=%s." % self.odbPath)
            def stepFrameToName(frameId):
                frameName = (
                    self.odb.postData["frameNames"][frameId[0]][frameId[1]])
                msg("frameName from postData: step %d, frame %d, frameName %s"
                    % (frameId[0], frameId[1], frameName), debugLevel=10)
                return frameName

        # no conversion otherwise
        else:
            stepFrameToName = self.stepFrameToName
        return stepFrameToName

    def frameStepperWithMaxFields(self):
        if not self.frameList:
            # no frames...
            return
        if self.frameList[0][0] != self.frameList[-1][0]:
            raise NotImplementedError(
                "Different steps with fieldNamesMaximize option not implemented"
                " (yet)!\nframeList: %s" % self.frameList)
        frameListAll = [
            (self.frameList[0][0], frameNb)
            for frameNb in range(self.frameList[0][1]+1,
                                 self.frameList[-1][1]+1)]
        frameList = set(self.frameList[1:])
        for stepNb, frameNb in frameListAll:
            yield (stepNb, frameNb, (stepNb, frameNb) in frameList)

    def frameStepperStandard(self):
        for stepNb, frameNb in self.frameList:
            yield (stepNb, frameNb, True)

    def openOdb(self):
        """prepare: open odb
        Initialize instance variable self.odb.
        """
        msg("Opening the odb %s." % self.odbPath)
        self.odb = OdbReaderInterpolateToPoints(
            self.odbPath, version=self.abaqusVersion)

    def initOutput(self):
        """Initialize output.
         1. check output directory
         2. preprocess self.outputFields:
            If outputFields is given then convert it into list of tuples:
             - if composer --third item of an outputFields-tuple-- is callable
               then create (field-class, composer function) tuples
             - if composer is string create (output field name, odb field name)
               tuples
             - if there is a fourth component "isConst" in the outputFields-
               tuple then this output field will be stored in
               self.outputFieldsConst instead of self.outputFields.

        Optionally overload this function to conduct further initalization.

        @Note: For a working subclass either initOutput or
        L{initOutputForFrame} have to create an attribute self.output that
        contains a mesh attribute. This attribute self.output.mesh might be
        supplied to the getExtraFields function.
        """
        # check output directory
        if self.outputDirName and not os.path.isdir(self.outputDirName):
            os.makedirs(self.outputDirName)
            msg("Created output dir %s" % self.outputDirName)

        # for outputFields option:
        # if outputFields argument is given convert into list of tuples
        # if composer is callable create (field-class, composer function) tuples
        # if composer is string create (output field name, odb field name)
        # tuples
        if self.outputFields:
            msg("Before converting outputFields:\n%s" % ",\n".join(
                isinstance(item, tuple) and (
                    "- fld name %s from %s"
                    % (item[0], (callable(item[2]) and "<func>" or item[2])))
                or ("- fld %s" % item)
                for item in self.outputFields), debugLevel=10)

            def convertOutputFldItem(item):
                # outputFields are (fieldName, fieldType, composer,  ...
                # ... optional isConst)-tuples
                if isinstance(item, basestring):
                    # if item is just a string
                    return (item, item)
                elif callable(item[2]):
                    return (
                        createFieldClass(
                            item[0], self.outputFieldType, item[1]), item[2])
                else:
                    return (item[0], item[2])

            self.outputFieldsConst = [
                convertOutputFldItem(item)
                for item in self.outputFields
                if len(item)>3 and item[3]]
            self.outputFields = [
                convertOutputFldItem(item)
                for item in self.outputFields
                if len(item)<=3 or not item[3]]
            msg("After converting outputFields:\n%s" % ",\n".join(
                "- output fld %s from odb fld %s"
                % (item[0], (callable(item[1]) and "<func>" or item[1]))
                for item in self.outputFields), debugLevel=10)

    def getMatRegions(self):
        """Read mat regions from the odb, write mat-names-csv file and return
        mat-number-field.
        """

        # create matNameToMatNumber dict from matNamesList
        if self.matNamesList:
            matNamesList = self.matNamesList
            matNameToMatNumber = dict(
                (n,i+1) for i,n in enumerate(matNamesList) if n)
        else:
            matNameToMatNumber = None
            matNamesList = list()  # will be extended with names from the odb

        # get the material data from the odb
        msg("Reading material properties from the odb.")
        (extraMatNamesList, matNumberField) = self.odb.getMatFieldFromOdb(
            nameFilter=self.matRegionsNameFilter,
            firstMatCode=len(matNamesList)+1,
            matNameToMatNumber=matNameToMatNumber,
            noMatCode=0)
        # converting to float
        matNumberField = type(matNumberField)(
            float(x) for x in matNumberField)

        # update matNamesList with data from the odb
        if matNamesList and extraMatNamesList:
            msg("WARNING: Found %d materials not listed in the specified"
                " matNamesList argument." % len(extraMatNamesList))
            msg("It's suggested to check spelling and the result.")
        matNamesList.extend(extraMatNamesList)

        # store possibly updated material names list in instance attribute
        self.matNamesList = matNamesList

        return matNumberField

    def writeMatNamesListToFile(self):
        """Write the instance attribute self.matNamesList to a csv file.
        This function should only be called after (the last invocation of)
        L{getMatRegions}.
        """
        outputPathMatList = os.path.join(
            self.outputDirName, self.fileNamePrefix+"_matList.csv")

        # write material number to names list
        msg("Writing material list to the csv file %s"
            % outputPathMatList)
        self.odb.writeMatNamesToCsv(outputPathMatList, self.matNamesList)

    def storeConstField(self, field):
        """Store constant field, e.g. the matNumberField
        """
        raise NotImplementedError(
            "storeConstField is not defined for class %s." % self.__class__)

    def initRefData(self):
        """read (and store) reference data for relative fields
        Updates instance variable self.refData.
        """
        stepNb, frameNb = self.relativeToFrame
        for fieldName in self.relativeFieldNames:
            msg("Reading reference field %s from step %d, frame %d."
                % (fieldName, stepNb, frameNb))
            self.refData[fieldName] = self.odb.getFieldFromOdb(
                stepNb, frameNb, fieldName,
                interpolationType=self.fieldInterpolType[fieldName])

    def processConstantFields(self):
        """Finally process constant fields before we come to the frames-loop.
        Overload this function in subclasses that need it.
        """
        return

    def initOutputForFrame(self, stepNb, frameNb):
        """Prepare for new frame data before looping over fields to be read
        from the odb. Intended to initialize output and will only be
        called for frames for which output is requested.
        Overload this function in subclasses that need it.

        @Note: For a working subclass either L{initOutput} or
        initOutputForFrame have to create an attribute self.output that
        contains a mesh attribute. This attribute self.output.mesh might be
        supplied to the getExtraFields function.
        """
        return

    def storeField(self, field):
        """Store the given field for output. Will be called once for each field
        and multiple times for each frame.
        Overload this function in subclasses.
        """
        raise NotImplementedError(
            "storeField is not defined for class %s." % self.__class__)

    def storeFrame(self, stepNb, frameNb):
        """Store frame data to output. This will be called once for each frame
        after all fields have been processed by self.storeField().

        It also needs to process self.outputFields by means of
        L{self.processOutputFields}()
        """
        raise NotImplementedError(
            "storeFrame is not defined for class %s." % self.__class__)

    def processOutputFields(self, odbFieldContainer):
        """Process self.outputFields for the current frame (fields from odb
        already read --passed in though odbFieldContainer argument).

        Generator function that yields field objects to be written to the
        output according to self.outputFields.

        @param odbFieldContainer: a dict containing the fields from the odb.
        """
        for outputClassOrName, composer in self.outputFields:
            if callable(composer):
                msg("Processing callable composer for field %s."
                    % outputClassOrName.fieldName, debugLevel=10)
                field = composer(
                    outputClassOrName, odbFieldContainer, self.output.mesh)
            else:
                msg("Processing non-composer field %s from odb-field %s."
                    % (outputClassOrName, composer), debugLevel=10)
                field = odbFieldContainer[composer]
                # make a copy of the type with changed fieldName
                FieldT = field.cloneType(fieldName=outputClassOrName)
                field = FieldT(field)
            yield field

    def closeOutput(self):
        """Finalize output. Will be called after all frames have been stored
        by means of L{storeFrame}.

        Overload this function in subclasses that need it.
        """
        return

    def run(self):
        """Main function that does the actual job. Calls all the other
        functions, loops over frames and fields.

        @returns: 0 on success, 1 or other on error. This value can be used
        as return code of the script in compliance with common practice.
        """
        retCode = 0  # ok, no error
        self.openOdb()
        self.initOutput()

        # mat regions preparation
        if self.matRegionsOption:
            self.storeConstField(self.getMatRegions())
            # also immidiately write the names to a csv file
            self.writeMatNamesListToFile()

        self.processConstantFields()

        # for fieldNamesMaximize option:
        # initialize frameListAll and fieldsMaximize-dict
        if self.fieldNamesMaximize:
            if self.fieldMaximizeAfter:
                # apply fieldMaximizeAfter to frameList
                self.frameList.insert(0, self.fieldMaximizeAfter)
            frameIter = self.frameStepperWithMaxFields()
            fieldsMaximize = dict()
            fieldNamesMaximize = set(self.fieldNamesMaximize)
        else:
            frameIter = self.frameStepperStandard()

            # we need to distinguish between fieldNamesMaximize not being
            # specified initially (=None) and fieldNamesMaximize becoming empty
            # later in the process. In the latter case it will be an empty set.
            fieldNamesMaximize = None

        # process relativeFieldNames: initialize refData
        self.refData = dict()
        if self.relativeFieldNames:
            self.initRefData()

        # get the field data from the odb frame by frame
        msg("Now starting to iterate over output frames.")
        stopNow = False
        for stepNb, frameNb, isOutputFrame in frameIter:

            if isOutputFrame:
                self.initOutputForFrame(stepNb, frameNb)
                currentFieldNames = self.fieldNames
            else:
                currentFieldNames = self.fieldNamesMaximize

            for fieldName in currentFieldNames:
                interpolationType = self.fieldInterpolType[fieldName]
                iptText = ""
                if interpolationType=="const":
                    iptText = " interpolating piecewise constant"
                msg("Reading field %s%s, step %d, frame %d."
                    % (fieldName, iptText, stepNb, frameNb))
                try:
                    field = self.odb.getFieldFromOdb(
                        stepNb, frameNb, fieldName,
                        interpolationType=interpolationType)
                except ValueError as exc:
                    if self.raiseMissingData:
                        raise
                    else:
                        msg("ERROR: Field %s for odb frame (%s, %d) not"
                            " available. Ending this loop now.\nError: %s"
                            % (fieldName, stepNb, frameNb, exc))
                        stopNow = True
                        retCode = 1
                        break

                # subtract reference field values (e.g. displacements)
                try:
                    refField = self.refData[fieldName]
                except KeyError:
                    pass
                else:
                    if field.dataType == "scalar":
                        field = type(field)(
                            (a-b) for a, b in izip(field, refField))
                    else:
                        field = type(field)(
                            [(a[i]-b[i]) for i in range(field.nbComponents)]
                            for a, b in izip(field, refField))

                # compute the current maximum of fields in fieldNamesMaximize
                if fieldNamesMaximize and fieldName in fieldNamesMaximize:
                    if field.dataType != "scalar":
                        msg("ERROR: Can only find maximum of scalar field."
                            " Field %s of dataType %s will not be maximized."
                            " Instead the current value will be stored in the"
                            " results file." % (fieldName, field.dataType))
                        fieldNamesMaximize.remove(fieldName)
                    else:
                        try:
                            oldMax = fieldsMaximize[fieldName]
                        except KeyError:
                            # no previous max value. Take current field as max.
                            pass
                        else:
                            msg("Updating max value for %s..." % fieldName)
                            field = type(oldMax)(
                                max(a, b) for a, b in izip(oldMax, field))
                        # store current maximum value in dict fieldsMaximize
                        fieldsMaximize[fieldName] = field

                # possibly also store the field for output
                if isOutputFrame:
                    self.storeField(field)
                    msg("Stored field %s in output." % fieldName, debugLevel=10)

            # after for fieldName ... loop
            if stopNow:
                break

            if isOutputFrame and self.getExtraFields:
                # fields external sources, e.g. from Cavesim vtk
                getExtraFieldsArgs = dict(
                    stepNb=stepNb, frameNb=frameNb,
                    points=self.output.mesh, odb=self.odb)
                try:
                    for fld in self.getExtraFields(**getExtraFieldsArgs):
                        self.storeField(fld)
                except IOError as exc:
                    if self.raiseMissingData:
                        raise
                    else:
                        msg("ERROR: Extra fields for odb frame (%s, %d) not"
                            " available. Ending this loop now.\nError: %s"
                            % (stepNb, frameNb, exc))
                        retCode = 2
                        break

                del getExtraFieldsArgs

            # self.storeFrame evaluates self.outputFields
            # ... and possibly writes files
            if isOutputFrame:
                self.storeFrame(stepNb, frameNb)

                if fieldNamesMaximize:
                    msg("Clearing %d fields from max-container"
                        % len(fieldsMaximize), debugLevel=10)
                    fieldsMaximize.clear()

        # finalize output
        self.closeOutput()
        return retCode


#------------------------------------------------------------------------------
# odbToVtkGridBox
#------------------------------------------------------------------------------
class OdbToVtkGridBox(OdbToPointsData):
    """Read interpolated field data from the specified odb at points on the
    specified regular grid and write the data to legacy vtk files.

    Don't use this class directly, its interface may change.
    Use L{odbToVtkGridBox} instead.
    """
    outputFieldType = "structPt"

    def __init__(
            self,
            odbPath, fieldNames, frameList,
            gridName, gridData,
            outputDirName, fileNamePrefix,
            projectPrefix, meshVersion,
            fieldNamesMaximize=None, fieldMaximizeAfter=None,
            relativeFieldNames=None, relativeToFrame=None,
            relativeToReferenceVtkPath=None,
            matRegionsOption=False, matRegionsNameFilter=None,
            matNamesList=None,
            outputFields=None,
            getExtraFields=None,
            stepFrameToName="F%d%03d",
            skipExisting=False,
            raiseMissingData=True,
            abaqusVersion=None,
            outputFormat="binary",
            interpolationDataDir=".",
            ):
        """Constructor does only basic argument checking and stores all
        arguments in instance variables
        """
        if relativeFieldNames:
            if relativeToFrame and relativeToReferenceVtkPath:
                raise ValueError(
                    "Both arguments relativeToFrame and"
                    " relativeToReferenceVtkPath are given but they are"
                    " mutually exclusive.")
            if relativeToFrame is None and relativeToReferenceVtkPath is None:
                raise ValueError(
                    "relativeFieldNames are given but neither relativeToFrame"
                    " nor relativeToReferenceVtkPath are given. We need"
                    " exactly one of them.")

        OdbToPointsData.__init__(
            self, fieldNames=fieldNames,
            fieldNamesMaximize=fieldNamesMaximize,
            fieldMaximizeAfter=fieldMaximizeAfter,
            relativeFieldNames=relativeFieldNames,
            relativeToFrame=relativeToFrame)

        # store arguments as instance variables
        # Note: fieldNames already processed and stored by
        # OdbToPointsData.__init__()

        # make a copy: fieldMaximizeAfter-option might alter the list...
        self.frameList = list(frameList)

        self.storeArgs(
            locals(),
            "odbPath",
            "gridName", "gridData",
            "outputDirName", "fileNamePrefix",
            "projectPrefix", "meshVersion",
            "fieldNamesMaximize", "fieldMaximizeAfter",
            "relativeFieldNames", "relativeToFrame",
            "relativeToReferenceVtkPath",
            "matRegionsOption", "matRegionsNameFilter",
            "matNamesList", "getExtraFields", "outputFields",
            "stepFrameToName",
            "skipExisting",
            "raiseMissingData",
            "abaqusVersion",
            "outputFormat",
            "interpolationDataDir")

    def initOutput(self):
        """Prepare everything for the output of results
         1. call L{OdbToPointsData.initOutput}
         2. compose output file names
         3. preprocess stepFrameToName
         4. initialize output grid and output containers
         5. initialize PtName field
        """
        # store outputFieldNames from self.outputFields before it's being
        # processed by OdbToPointsData.initOutput()
        if self.outputFields:
            outputFieldNames = [item[0] for item in self.outputFields]
        else:
            # in case of outputFields not given: copy(!) fieldNames
            outputFieldNames = list(self.fieldNames)
            if self.matRegionsOption:
                outputFieldNames.append("matcode")

        # check output directory and process outputFields
        OdbToPointsData.initOutput(self)

        # pickle data file name
        # from interpolationDataDir and others
        interpolPickleName = os.path.join(
            self.interpolationDataDir, "_".join((
                self.projectPrefix, self.meshVersion, self.gridName)))

        # check interpolationDataDir directory
        if self.interpolationDataDir and not os.path.isdir(
                self.interpolationDataDir):
            os.makedirs(self.interpolationDataDir)
            msg("Created directory for interpolation data %s"
                % self.interpolationDataDir)

        # output file names
        # from outputDirName and fileNamePrefix
        self.outputFileNameTemplate = os.path.join(
            self.outputDirName, "%s_%%s.vtk" % self.fileNamePrefix)

        # check argument stepFrameToName (preprocess if string)
        if isinstance(self.stepFrameToName, basestring):
            if self.stepFrameToName=="fromPostData":
                # special case, frameNames from odb.postData
                if "frameNames" not in self.odb.postData:
                    raise ValueError(
                        "frameNames not defined in associated _postData.py file"
                        " for odbPath=%s." % self.odbPath)
                def stepFrameToName(stepNb,frameNb):
                    frameName = self.odb.postData["frameNames"][stepNb][frameNb]
                    msg("frameName from postData: step %d, frame %d,"
                        " frameName %s"
                        % (stepNb, frameNb, frameName), debugLevel=10)
                    return frameName
            else:
                # stepFrameToName is format string
                stepFrameToNameStr = self.stepFrameToName
                def stepFrameToName(stepNb,frameNb):
                    return stepFrameToNameStr % (stepNb, frameNb)
            self.stepFrameToName = stepFrameToName

        # write csv file describing the resulting vtk files
        outputFileNameVtkDescription = os.path.join(
            self.outputDirName, self.fileNamePrefix+"_vtk_description.csv")
        if self.skipExisting and os.path.isfile(outputFileNameVtkDescription):
            msg("WARNING: vtk-description csv file %s already exists. Skipping"
                " creation of a new one. No check is performed if it's"
                " complete." % outputFileNameVtkDescription)
        else:
            output = csv.writer(open(outputFileNameVtkDescription, "wb"))

            # description of output fields
            output.writerow(["index", "component"])
            for i, fieldName in enumerate(outputFieldNames):
                output.writerow([i+1, fieldName])

            # description of output frames
            try:
                frameNames = self.odb.postData["frameNames"]
            except IOError:
                msg("WARNING: Did not find ..._postData.py file. The"
                    " vtk_description.csv file will not contain frame-date"
                    " information.")
            else:
                output.writerow([])
                output.writerow(["odb frame", "date"])
                for stepNb, frameNb in self.frameList:
                    frameLabel = self.stepFrameToName(stepNb, frameNb)
                    try:
                        frameName = frameNames[stepNb][frameNb]
                    except IndexError:
                        msg("ERROR: Could not find step %d frame %d in"
                            " postdata frameNames list. Check postData and"
                            " frameLists!"
                            % (stepNb, frameNb))
                        continue
                    output.writerow([frameLabel, frameName])

            # close description file
            del output
            msg("Wrote vtk file description file %s"
                % outputFileNameVtkDescription)

        # initialize the output grid and interpolation
        self.grid = getMesh(**self.gridData)
        self.odb.initOutputPoints(
            interpolPickleName=interpolPickleName, **self.gridData)

        # initialize Container for mat field
        self.constFieldsContainer = OrderedDict()

        # check if resultsfiles already exist
        if self.skipExisting:
            lastExisting = None
            for i, (stepNb, frameNb) in enumerate(self.frameList):
                # compose outputFileName
                frameName = self.stepFrameToName(stepNb, frameNb)
                outputFileName = self.outputFileNameTemplate % frameName
                if os.path.isfile(outputFileName):
                    lastExisting = (i, frameName)
            if lastExisting is not None:
                msg("WARNING: Found result files, last for frame %s. Will skip"
                    " everything up to that frame. NO CHECK IS PERFORMED IF"
                    " OLD DATA IS COMPLETE!"
                    % (lastExisting[1]))
                self.fieldMaximizeAfter = self.frameList[lastExisting[0]]
                self.frameList = self.frameList[lastExisting[0]+1:]
                msg("New fieldMaximizeAfter: %s. New frameList: %s"
                    % (self.fieldMaximizeAfter, self.frameList),
                    debugLevel=5)
            else:
                msg("No files to skip like %s, frameList: %s"
                    % (self.outputFileNameTemplate % "XXXXX", self.frameList),
                    debugLevel=5)

    def storeConstField(self, field):
        """Store constant field (i.e. "material") in self.constFieldsContainer.
        """
        self.constFieldsContainer[field.fieldName] = field

    def initRefData(self):
        """read (and store) reference data for relative fields
        Updates instance variable self.refData.

        This variant recognizes and processes the relativeToReferenceVtkPath
        argument.
        """
        if self.relativeToReferenceVtkPath:
            # read reference data from vtk file
            msg("Reading reference field data from vtk file %s."
                % self.relativeToReferenceVtkPath)
            vtk = FieldsOnStructPointGrid().fromVtk(
                self.relativeToReferenceVtkPath)

            # checks and diagnostic output
            if self.grid != vtk.mesh:
                msg("ERROR: Output grid <%s> not identical to the grid from"
                    " the reference data vtk file %s."
                    % (self.gridName, self.relativeToReferenceVtkPath))
            msg("Read reference fields from vtk file: %s"
                % ", ".join(fieldName for fieldName in vtk.data))
            if any(fieldName not in vtk.data
                   for fieldName in self.relativeFieldNames):
                raise ValueError(
                    "Didn't find the required fields in the reference data vtk"
                    " file %s. You specified relativeFieldNames=%s. In the vtk"
                    " file there are the following fields: %s."
                    % (self.relativeToReferenceVtkPath, self.relativeFieldNames,
                       vtk.data.keys()))

            # store the reference field data
            self.refData = vtk.data
        else:
            # read reference data from odb
            OdbToPointsData.initRefData(self)

    def initOutputForFrame(self, stepNb, frameNb):
        """Prepare for new frame data before looping over fields to be read
        from the odb.

        Initialize self.output
        """
        self.output = FieldsOnStructPointGrid(mesh=self.grid)

    def storeField(self, field):
        """Store the given field for output. Will be called once for each field
        and multiple times for each frame.
        """
        self.output.updateFields(field)

    def storeFrame(self, stepNb, frameNb):
        """Store frame data to output. This will be called once for each frame
        after all fields have been processed by self.storeField().
         1. Possibly adds material field to output
         2. Processes self.outputFields by means of
            L{OdbToPointsData.processOutputFields}().
        """
        # Possibly add material field to output
        if self.matRegionsOption:
            # add material number field
            matNumberField = self.constFieldsContainer["Material"]
            self.output.updateFields(matNumberField)

        # compose output fields if outputFields argument is given
        if self.outputFields:
            newOutput = FieldsOnStructPointGrid(mesh=self.grid)
            newOutput.updateFields(
                *self.processOutputFields(self.output.data))
            self.output = newOutput
            del newOutput
            msg("Converted output according to outputFields parameter.",
                debugLevel=10)

        # compose outputFileName
        frameName = self.stepFrameToName(stepNb, frameNb)
        outputFileName = self.outputFileNameTemplate % frameName

        # write data to vtk-file
        msg("Writing data for step %d, frame %d to %s"
            % (stepNb, frameNb, outputFileName))
        self.output.toVtk(
            outputFileName,
            description=("Abq-Grid-Data step %d, frame %d"
                         % (stepNb, frameNb)),
            outputFormat=self.outputFormat)
        del self.output

#------------------------------------------------------------------------------
def odbToVtkGridBox(
        odbPath, fieldNames, frameList,
        gridName, gridData,
        outputDirName, fileNamePrefix,
        projectPrefix, meshVersion,
        fieldNamesMaximize=None, fieldMaximizeAfter=None,
        relativeFieldNames=None, relativeToFrame=None,
        relativeToReferenceVtkPath=None,
        matRegionsOption=False, matRegionsNameFilter=None,
        matNamesList=None,
        outputFields=None,
        getExtraFields=None,
        stepFrameToName="F%d%03d",
        skipExisting=False,
        raiseMissingData=True,
        abaqusVersion=None,
        outputFormat="binary",
        interpolationDataDir=".",
        ):
    """Read interpolated field data from the specified odb at points on the
    specified regular grid and write the data to legacy vtk files.

    B{Option fieldNamesMaximize:} For selected fields take the maximum value of
    all intermediate time frames. Data for fields in fieldNamesMaximize for a
    particular frame listed in frameList will be read from one frame after the
    previous entry in frameList up the mentioned frame itself. The
    maximum for all those frames will be computed and stored in an output file
    for the mentioned frame in frameList together with all other requested
    fields from this frame. Then again only fields in fieldNamesMaximize will
    be read from the next frame (i.e. frameList[1][1]+1) up to the next frame
    listed in frameList yielding output for this latter frame (i.e. the third
    in the frameList).

    If fieldMaximizeAfter is given as well then reading of fields listed in
    fieldNamesMaximize starts at the frame following after fieldMaximizeAfter.
    All other fields are read for the first frame in frameList and the first
    output will be generated for this first frame in frameList as well.

    B{For example} fieldNames=["LogP", "RER"], fieldNamesMaximize=["RER"],
    fieldMaximizeAfter=(2,4), frameList=[(2,7), (2,10), (2,13)]:
    RER will be read for frame 5,6,7 and LogP will be read for frame 7 only.
    The first output file will be for frame 7 including LogP for frame 7 and
    the maximum of RER for frames 5,6,7. Then RER will be read for frame 8,9
    and 10 and LogP for frame 10 only. Second output file for frame 10 contains
    LogP from frame 10 and for RER the maximum of frame 8, 9, 10. Then from
    frames 11 and 12 only LogP will be read. From frame 13 RER and LogP will
    be read. Next output for frame 13 contains LogP (fr 13) and maximum RER
    from frame 11 to 13.

    If fieldMaximizeAfter is not given (or None) then the first item of
    frameList will be used in place of fieldMaximizeAfter. There will be no
    output for this first frame in this case.

    B{For example} fieldNames=["LogP", "RER"], fieldNamesMaximize=["RER"],
    fieldMaximizeAfter=None, frameList=[(2,7), (2,10), (2,13)]:
    RER will be read for frame 8,9,10 and LogP will be read for frame 10 only.
    The first output file will be for frame 10 including LogP for frame 10 and
    the maximum of RER for frames 8,9,10. Then RER will be read for frame 11,12
    and 13 and LogP for frame 13 only. Second output file for frame 13 contains
    LogP from frame 13 and for RER the maximum of frame 11, 12, 13.

    Note that if fieldNamesMaximize is specified then all time frames up to and
    including fieldMaximizeAfter or --if not specified-- the first in frameList
    will not be considered in any way. There will be no output file for the
    first frame in the frameList if fieldMaximizeAfter is not specified.

    B{Option relativeFieldNames:} For certain fields (e.g. displacement) output
    values relative to a reference frame / reference field.
    Specify the fields to report relative to a certain reference value in the
    argument relativeFieldNames. Additionally you need either argument
    relativeToFrame or relativeToReferenceVtkPath. These two arguments are
    mutually exclusive.

    B{Option matRegionsOption:} Add material region number field to every output
    vtk file. The material code field is named "Material".

    If you want the material region number only then choose a single
    arbitrary frame for the frameList parameter and make fieldNames an empty
    list and of course matRegionsOption=True.

    B{Option outputFields:} If the argument outputFields is given then only
    fields in this argument will appear in the output files. The fields listed
    in the argument fieldNames will be read from the odb and can be used to
    compute output fields.

    Please use keyword arguments at least for all optional arguments as their
    order might change in future (sub-)versions.

    @param odbPath: Complete path to the odb including .odb extension. Can be
      relative or absolute.

    @param fieldNames: List of field names to be exported to the output vtk.
      Must include all fields in argument fieldNamesMaximize.
      Example: ["UR_1","UR_2","UR_3","SDV4","S_MIN_PRINCIPAL","SDV12"]

      Instead of a field name the items of fieldNames can be
      (field name, interpolation type)- tuples. Interpolation type can be
      "default" or "const" for piecewise constant interpolation.
      Example: ["SDV4", ("SDV8", "const"), "SDV12"]
    @param fieldNamesMaximize: Optional list of field names to take the maximum
      value of intermediate frames of. For maxRER. Example: ["SDV12"]
    @param fieldMaximizeAfter: If given then the frameList argument again lists
      all frames where output is generated for, including the very first.
      Specify a (stepNr, frameNr)-tuple for fieldMaximizeAfter.
      For the first output frame (the first item in the frameList-argument) data
      from frames *after* this --i.e. starting with (stepNr, frameNr+1)-- will
      be considered.
    @param relativeFieldNames: Optional. From fields listed in this argument
      the value of the time frame relativeToFrame is being subtracted before
      export. E.g. for relative displacements (relative to the equilibrium
      step) let relativeFieldNames=["U_1", "U_2", "U_3"] and
      relativeToFrame=(2.0).
    @param relativeToFrame: A (stepNb, frameNb)-tuple specifying the reference
      time frame for the relativeFieldNames option. Reference field data will
      be read from the current odb from this time frame. This argument is
      mutually exclusive to relativeToReferenceVtkPath.
    @param relativeToReferenceVtkPath: Complete path to a vtk file that
      contains all fields listed in relativeFieldNames. This argument is
      mutually exclusive to relativeToFrame. Field names and output grid must
      match!
    @param matRegionsOption: Boolean flag: Write MATERIAL REGIONS data as well?
      If True then a field named "Material" will be generated.
    @param matRegionsNameFilter: This is for material name conversion. For
      example if there are region variants that shall be treated as one. I.e.
      if SEDIMENT, SEDIMENT_CAVE and SEDIMENT_UCUT shall all be treated as
      SEDIMENT.

      Must be a function otherwise it'll be ignored. It will be called as a
      conversion function taking the original material name from the odb as
      first argument and the section type (e.g. "solid", "beam", "shell") as
      second argument and returning the material name it should be treated as.

      Example: To take everything up to the first "_" or "-" use as
      matRegionsNameFilter:
       >>> lambda matname,secType: matname.split("_")[0].split("-")[0]

      That's equivalent to:
       >>> def convert(matname,secType):
       >>>     return matname.split("_")[0].split("-")[0]
       >>> odbToVtkGridBox(..., matRegionsNameFilter=convert, ...)

    @param matNamesList: Optional list determinining the material name to
      material number assignment. Material numbers in the output will be
      assigned in the order determined by the list starting with 1.

      The specified material names must match the final material names as
      returned from the odb and possibly modified by the nameFilter function
      if applicable.

      Typically these are the uppercase material names as specified
      in the Abaqus input file of the job.

      Example:
       >>> matNamesList=["HOST", "ULTRAMAFIC", ...],

      Note: The resulting ..._matList.csv file will contain all material names
      listed with this argument plus any additional material name that may
      appear in the odb appended to the end.

    @param getExtraFields: Optional generator function to supply more fields
      from sources other than odb fields. To be used in conjuction with the
      outputFields argument. Will be called with some keyword arguments.

      Currently implemented arguments are stepNb, frameNb, points, odb.
      points is a L{MeshStructuredPoints} object identifying the topology
      for the fields to be supplied. odb is the
      L{OdbReaderInterpolateToPoints
      <bae.odb_03.odbReader.OdbReaderInterpolateToPoints>} objects that is used
      to get all the data from the odb. This can be used to get all sorts of
      data other than field data from the odb. The odb attribute contains
      topological data of the points as well: through the
      L{bae.mesh_01.InterpolMeshToPoints} object odb.interpolation.

      B{Important:} The argument list in the definition of getExtraFields must
      always end with I{**dummy} because in future versions of this module more
      arguments might be passed to getExtraFields. See L{Converter_CombinedVtk}
      for an example.

      The generator getExtraFields yields L{Field} objects that can then be
      used in the outputFields argument. See L{Converter_CombinedVtk} for an
      example.

      B{Note:} In a derived subclass calling the super-class getExtraFields
      use keyword arguments --at least for all but the first three arguments!
      And consider the following:
      getExtraFields is a generator function. Its result acts as an iterator.
      Therefore after calling the super-class getExtraFields method you have
      to iterate over its result and pass those items to the caller through
      the yield command. You can't just return or yield the result of the
      super-class' getExtraFields method.

    @param outputFields: a list of (fieldName, fieldType, composer)-tuples.
      If given then output will only contain fields listed here.

      FieldName will be the name in the vtk file (visible for example in
      Paraview but not in voxler), fieldType can be "scalar", "vector",
      "tensor" or "str". fieldType can be None or arbitrary if composer is
      a string. composer describes how to generate the result field.
      composer can be a string or a function.

      If composer is a string then it will be taken as field name. The field
      with that name is taken from the odb without modification. fieldType will
      be ignored. This is to simply take fields from the odb as they are and
      possibly rename them.

      If composer is a function then it will be called once for each output
      frame with the following arguments: fieldClass, fieldFromOdbDict, grid.

      fieldClass is the suitable type for the result field derived from
      L{bae.field_01.Field}.

      fieldFromOdbDict is a dictionary {fieldName : data} containing the
      fields as they would occur in the output if the outputFields argument
      would have been ommited (values read from the odb as listed in the
      fieldNames argument; arguments fieldNamesMaximize, relativeFieldNames,
      matRegionsOption are recongized as well and have their usual effect,
      in this case on the data in fieldFromOdbDict)

      grid is a L{bae.mesh_01.MeshStructuredPoints} object.

      Example:
       >>> def composerLogP(cls, dat, mesh):
       >>>     return cls(log10(x*1000+1) for x in dat["SDV3"])
       >>> def composerS1(cls, dat, mesh):
       >>>     return cls(-x for x in dat["S_MIN_PRINCIPAL"])
       >>> ...
       >>> odbToVtkGridBox(...
       >>>   fieldNames=["SDV3", "S_MIN_PRINCIPAL"],
       >>>   matRegionsOption=True,
       >>>   outputFields=[
       >>>     ("logP", "scalar", composerLogP),
       >>>     ("RER", "scalar", "SDV6"),
       >>>     ("S1", "scalar", composerS1),
       >>>     ("xyz", "vector", lambda cls,dat,msh: cls(msh.getPointsIter())),
       >>>     ("Material", "scalar", "Material"),
       >>>     ],
       >>>   ...)

      Note how the "Material" field (i.e. the matcode) is transferred to the
      output without change. And SDV6 is only renamed to RER.

      Note: To pass values between different calls to the composer function(s)
      you can use global variables. Or give your composer an additional default
      argument with a list (or other container object) as default value.

      Note: If the Parameter outputFields is given then you have to have
      "Material" listed here otherwise this field will not be in the output
      even if you specified the material option.

    @param gridData: A dict with keys firstPoint, lastPoint and spacing and
      optionally rotString.
      All three values being lists of three floats. This will be passed on as
      keyword-arguments to the L{bae.mesh_01.MeshStructuredPoints} constructor.

    @param gridName: Used to identify the regular grid, in particular the
      corresponding interpolation data.

    @param frameList: List of output frames as (step number, frame number)
      tuples. Example: [ (2, 24), (2, 27) ]

      In case fieldNamesMaximize is given and fieldMaximizeAfter is not given
      or None then this list is treated slightly different: The first frame in
      this list is the last frame not to be considered at all. No output is
      generated for this first frame! For subsequent frames output will be
      generated.
      Example: [ (2, 20), (2, 24), (2, 28) ] will generate output for frame
      24 and frame 28. Not for frame 20! The results for frame 24 will include
      maximum RER from frame 21 to frame 24 (both included). For frame 28 we'll
      have maximum RER from frame 25 to frame 28 (both included). Provided
      that RER (i.e. SDV12 for example) is listed in fieldNamesMaximize.

    @param stepFrameToName: Defines how to build frame names in the vtk output
      file names from a (stepNb, frameNb)-tuple.
      Might be a function (stepNb, frameNb)-tuple --> frameName or a template
      string (like the default "F%d%03d").
      Can also be the string "fromPostData" in this case the frameName attribute
      from the associated _postData.py file is used.

      For restart analyses use something like this (We only have step 2 and 3,
      the latter sometimes referred to as 1):
       >>> def stepFrameToName(stepNb, frameNb):
       >>>     if stepNb==2:
       >>>         return "F2%03d" % frameNb
       >>>     else:
       >>>         return "F3%03d" % frameNb
       >>>
       >>> odbToVtkGridBox( ..., stepFrameToName=stepFrameToName, ...)

    @param skipExisting: Set to True if you want to add later frames to
      already existing files. All output files appording to the specified
      frameList are checked if they exist already. The last of those is
      considered and the processing of new files starts right after this one.
      The last frame for which an output file already exists is considered
      as value for fieldMaximizeAfter.

    @param raiseMissingData: If True (default) then raise an exception if
      data could not be found in the odb or external file. Otherwise just
      write an error message and finish.

    @param outputDirName: Vtk files will end up in this folder. Can be
      a relative or an absolute path.
    @param fileNamePrefix: Base filename for output files.
      Usually projectPrefix_runVersion_seqVersion_gridName
      Example: "Perse2015_R01_Q01_boxD3_h5"
    @param projectPrefix: Needed for the interpolation data filename.
    @param meshVersion: Needed for the interpolation data filename.

    @param outputFormat: "binary" or "ascii"
    @param abaqusVersion: Version string of Abaqus, something like "6.13-2".
       Note: separators must be exactly as in the example: dot and dash.
       For new style use something like "abq2018", ignoring the hot-fix
       postfix e.g. "hf4".
       If not specified defaults to what the 'abaqus' shell command
       invokes.
    @param interpolationDataDir: .interpol files will be in that subfolder

    @returns: 0 on success, 1 or other on error. This value can be used
       as return code of the script in compliance with common practice.
    """
    runner = OdbToVtkGridBox(
        odbPath, fieldNames, frameList,
        gridName, gridData,
        outputDirName, fileNamePrefix,
        projectPrefix, meshVersion,
        fieldNamesMaximize=fieldNamesMaximize,
        fieldMaximizeAfter=fieldMaximizeAfter,
        relativeFieldNames=relativeFieldNames,
        relativeToFrame=relativeToFrame,
        relativeToReferenceVtkPath=relativeToReferenceVtkPath,
        matRegionsOption=matRegionsOption,
        matRegionsNameFilter=matRegionsNameFilter,
        matNamesList=matNamesList,
        outputFields=outputFields,
        getExtraFields=getExtraFields,
        stepFrameToName=stepFrameToName,
        skipExisting=skipExisting,
        raiseMissingData=raiseMissingData,
        abaqusVersion=abaqusVersion,
        outputFormat=outputFormat,
        interpolationDataDir=interpolationDataDir,
        )
    retCode = runner.run()
    return retCode


#------------------------------------------------------------------------------
# odbToCsvPoints
#------------------------------------------------------------------------------
class OdbToCsvPoints(OdbToPointsData):
    """Read interpolated field data from the specified odb at specified points
    and write the data to csv files.

    Don't use this class directly, its interface may change.
    Use L{odbToCsvPoints} instead.
    """
    outputFieldType = "point"

    def __init__(
            self,
            odbPath, fieldNames, frameList,
            gridName, gridPoints,
            outputDirName, fileNamePrefix,
            projectPrefix, meshVersion,
            fieldNamesMaximize=None, fieldMaximizeAfter=None,
            relativeFieldNames=None, relativeToFrame=(2,0),
            matRegionsOption=False, matRegionsNameFilter=None,
            matNamesList=None,
            outputFields=None,
            getExtraFields=None,
            csvOutputType="CsvPerFrameRowPtColFld",
            filterEmpty=False,
            stepFrameToName="F%d%03d",
            raiseMissingData=True,
            abaqusVersion=None,
            interpolationDataDir=".",
            ):
        """Constructor does only basic argument checking and stores all
        arguments in instance variables
        """
        OdbToPointsData.__init__(
            self, fieldNames=fieldNames,
            fieldNamesMaximize=fieldNamesMaximize,
            fieldMaximizeAfter=fieldMaximizeAfter,
            relativeFieldNames=relativeFieldNames,
            relativeToFrame=relativeToFrame)

        # store arguments as instance variables
        # Note: fieldNames already processed and stored by
        # OdbToPointsData.__init__()

        # make a copy: fieldMaximizeAfter-option might alter the list...
        self.frameList = list(frameList)

        self.storeArgs(
            locals(),
            "odbPath",
            "gridName", "gridPoints",
            "outputDirName", "fileNamePrefix",
            "projectPrefix", "meshVersion",
            "fieldNamesMaximize", "fieldMaximizeAfter",
            "relativeFieldNames", "relativeToFrame",
            "matRegionsOption", "matRegionsNameFilter",
            "matNamesList",
            "getExtraFields",
            "outputFields",
            "csvOutputType",
            "filterEmpty",
            "stepFrameToName",
            "raiseMissingData",
            "abaqusVersion",
            "interpolationDataDir",
            )
        try:
            self.gridPoints = self.gridPoints.tolist()
        except AttributeError:
            pass


    def initOutput(self):
        """Prepare everything for the output of results
         1. call L{OdbToPointsData.initOutput}
         2. compose output file names
         3. preprocess stepFrameToName
         4. initialize output grid and output containers
         5. initialize PtName field
        """
        # check output directory and process outputFields
        OdbToPointsData.initOutput(self)

        # pickle data file name
        # from interpolationDataDir and others
        interpolPickleName = os.path.join(
            self.interpolationDataDir, "_".join((
                self.projectPrefix, self.meshVersion, self.gridName)))

        # check interpolationDataDir directory
        if self.interpolationDataDir and not os.path.isdir(
                self.interpolationDataDir):
            os.makedirs(self.interpolationDataDir)
            msg("Created dir for interpolation data %s"
                % self.interpolationDataDir)

        # output file names
        # from outputDirName and fileNamePrefix
        if self.csvOutputType.startswith("OneCsv"):
            outputFileNameTemplate = os.path.join(
                self.outputDirName, self.fileNamePrefix+".csv")
        else:
            outputFileNameTemplate = os.path.join(
                self.outputDirName, self.fileNamePrefix+"_%s.csv")

        # initialize the output grid
        if "PtName" in self.csvOutputType:
            self.gridPoints, ptNames = izip(*(
                (x[:3], x[3]) for x in self.gridPoints))
        else:
            ptNames = None
        points = MeshUnstructuredPoints(self.gridPoints)
        self.odb.initOutputPoints(
            grid=points, interpolPickleName=interpolPickleName)

        # initialize the output
        self.output = FramesFieldsOnUnstructuredPoints(points=points)
        self.output.initOutput(
            outputFileNameTemplate,
            outputType=self.csvOutputType, filterEmpty=self.filterEmpty,
            frameIdToName=self.prepareStepFrameToName())

        # initialize Container for constant fields
        self.constFieldsContainer = OrderedDict()

        # store ptNames field as constant field
        if ptNames:
            self.constFieldsContainer["ptName"] = createFieldObject(
                "ptName", "point", "str", initArgs=[ptNames,])

    def storeConstField(self, field):
        """Store constant field in self.constFieldsContainer
        """
        self.constFieldsContainer[field.fieldName] = field

    def processConstantFields(self):
        """Finally process constant fields before we come to the frames-loop.
         1. for outputFields option: store constant fields
         2. store constant field data in the output object
        """
        # for outputFields option:
        # store constant fields
        if self.outputFields and self.outputFieldsConst:
            newFldContainer = OrderedDict()
            for outputClassOrName, composer in self.outputFieldsConst:
                if callable(composer):
                    newFldContainer[outputClassOrName.fieldName] = composer(
                        outputClassOrName, self.constFieldsContainer,
                        self.output.mesh)
                else:
                    field = self.constFieldsContainer[composer]
                    # make a copy before changing the fieldName
                    field = type(field)(field)
                    field.fieldName = outputClassOrName
                    newFldContainer[outputClassOrName] = field
            self.constFieldsContainer = newFldContainer
            del newFldContainer
            msg("Converted constant field output according to outputFields"
                " parameter.", debugLevel=10)

        # actually store constant field data in the output object
        if self.constFieldsContainer:
            self.output.updateFields(*self.constFieldsContainer.itervalues())
            self.output.storeConstantData()
        return

    def initOutputForFrame(self, stepNb, frameNb):
        """Prepare for new frame data before looping over fields to be read
        from the odb.

        Initialize self.fieldsContainer
        """
        self.fieldsContainer = OrderedDict()

    def storeField(self, field):
        """Store the given field for output. Will be called once for each field
        and multiple times for each frame.
        """
        self.fieldsContainer[field.fieldName] = field

    def storeFrame(self, stepNb, frameNb):
        """Store frame data to output. This will be called once for each frame
        after all fields have been processed by self.storeField().

        Processes self.outputFields by means of
        L{OdbToPointsData.processOutputFields}().
        """
        # compose output fields if outputFields argument is given
        if self.outputFields:
            # join in constFieldsContainer
            if hasattr(self, "constFieldsContainer"):
                fieldsContainer = dict(self.fieldsContainer)
                fieldsContainer.update(self.constFieldsContainer)
            else:
                fieldsContainer = self.fieldsContainer

            # process self.outputFields
            newFldContainer = OrderedDict(
                (field.fieldName, field)
                for field in self.processOutputFields(fieldsContainer))
            self.fieldsContainer = newFldContainer
            del newFldContainer

        self.output.updateFields(*self.fieldsContainer.itervalues())
        del self.fieldsContainer
        self.output.storeFrame((stepNb, frameNb))

    def closeOutput(self):
        """Finalize output. Will be called after all frames have been stored
        by means of L{storeFrame}
        """
        self.output.closeCsvOutput()

#------------------------------------------------------------------------------
def odbToCsvPoints(
        odbPath, fieldNames, frameList,
        gridName, gridPoints,
        outputDirName, fileNamePrefix,
        projectPrefix, meshVersion,
        fieldNamesMaximize=None, fieldMaximizeAfter=None,
        relativeFieldNames=None, relativeToFrame=(2,0),
        matRegionsOption=False, matRegionsNameFilter=None,
        matNamesList=None,
        outputFields=None,
        getExtraFields=None,
        csvOutputType="CsvPerFrameRowPtColFld",
        filterEmpty=False,
        stepFrameToName="F%d%03d",
        raiseMissingData=True,
        abaqusVersion=None,
        interpolationDataDir=".",
        ):
    """Read interpolated field data from the specified odb at specified points
    and write the data to csv files.

    B{Option fieldNamesMaximize:} For selected fields take the maximum value of
    all intermediate time frames. Data for fields in fieldNamesMaximize for a
    particular frame listed in frameList will be read from one frame after the
    previous entry in frameList up the mentioned frame itself. The
    maximum for all those frames will be computed and stored in an output file
    for the mentioned frame in frameList together with all other requested
    fields from this frame. Then again only fields in fieldNamesMaximize will
    be read from the next frame (i.e. frameList[1][1]+1) up to the next frame
    listed in frameList yielding output for this latter frame (i.e. the third
    in the frameList).

    If fieldMaximizeAfter is given as well then reading of fields listed in
    fieldNamesMaximize starts at the frame following after fieldMaximizeAfter.
    All other fields are read for the first frame in frameList and the first
    output will be generated for this first frame in frameList as well.

    B{For example} fieldNames=["LogP", "RER"], fieldNamesMaximize=["RER"],
    fieldMaximizeAfter=(2,4), frameList=[(2,7), (2,10), (2,13)]:
    RER will be read for frame 5,6,7 and LogP will be read for frame 7 only.
    The first output file will be for frame 7 including LogP for frame 7 and
    the maximum of RER for frames 5,6,7. Then RER will be read for frame 8,9
    and 10 and LogP for frame 10 only. Second output file for frame 10 contains
    LogP from frame 10 and for RER the maximum of frame 8, 9, 10. Then from
    frames 11 and 12 only LogP will be read. From frame 13 RER and LogP will
    be read. Next output for frame 13 contains LogP (fr 13) and maximum RER
    from frame 11 to 13.

    If fieldMaximizeAfter is not given (or None) then the first item of
    frameList will be used in place of fieldMaximizeAfter. There will be no
    output for this first frame in this case.

    B{For example} fieldNames=["LogP", "RER"], fieldNamesMaximize=["RER"],
    fieldMaximizeAfter=None, frameList=[(2,7), (2,10), (2,13)]:
    RER will be read for frame 8,9,10 and LogP will be read for frame 10 only.
    The first output file will be for frame 10 including LogP for frame 10 and
    the maximum of RER for frames 8,9,10. Then RER will be read for frame 11,12
    and 13 and LogP for frame 13 only. Second output file for frame 13 contains
    LogP from frame 13 and for RER the maximum of frame 11, 12, 13.

    Note that if fieldNamesMaximize is specified then all time frames up to and
    including fieldMaximizeAfter or --if not specified-- the first in frameList
    will not be considered in any way. There will be no output file for the
    first frame in the frameList if fieldMaximizeAfter is not specified.

    B{Option relativeFieldNames:} For certain fields (e.g. displacement) output
    values relative to a reference frame / reference field.

    B{Option matRegionsOption:} Add material region number field to output.
    The material code field is named "Material".

    If you want the material region number only then choose a single
    arbitrary frame for the frameList parameter and make fieldNames an empty
    list and of course matRegionsOption=True.

    B{Option outputFields:} If the argument outputFields is given then only
    fields in this argument will appear in the output files. The fields listed
    in the argument fieldNames will be read from the odb and can be used to
    compute output fields.

    Please use keyword arguments at least for all optional arguments as their
    order might change in future (sub-)versions.

    @param odbPath: Complete path to the odb including .odb extension. Can be
      relative or absolute.

    @param fieldNames: list of field names (in capital letters)

      If you want a field like SDV4 just specify it as such: "SDV4", "POR",
      "STATUS", "U", "S"

      If you want a particular component of a vector or tensor specify the odb
      field name followed by an underscore and the component label: "U_1",
      "S_33"

      If you want a certain invariant of a vector or tensor specify the odb
      field name followed by an underscore and the invariant label in capital
      letters: "U_MAGNITUDE", "S_MIN_PRINCIPAL", ...

      Example:
       >>> fieldNames = ["U_1", "U_2", "U_3",
       >>>               "SDV4",   # LOG10(PEEQ*1000+1)
       >>>               "S_MIN_PRINCIPAL", "S_MAX_PRINCIPAL",
       >>>               "SDV12",   # max energy release rate per seqPeriod
       >>>               ]

      Instead of a field name items can be a (field name, interpolation type)-
      tuple instead. Interpolation type can be "default" or "const" for
      piecewise constant interpolation.
      Example: ["SDV4", ("SDV8", "const"), "SDV12"]

    @param fieldNamesMaximize: Optional list of field names to take the maximum
      value of intermediate frames of. For maxRER. Example: ["SDV12"]
    @param fieldMaximizeAfter: If given then the frameList argument again lists
      all frames where output is generated for, including the very first.
      Specify a (stepNr, frameNr)-tuple for fieldMaximizeAfter.
      For the first output frame (the first item in the frameList-argument) data
      from frames *after* this --i.e. starting with (stepNr, frameNr+1)-- will
      be considered.

    @param relativeFieldNames: optional. If you want fields
      -i.e. displacements- exported relative to a certain reference frame
      then specify a list of field names. E.g. ["U"]
    @param relativeToFrame: reference frame as
      (step number, frame number)-tuple. Defaults to (2, 0)

    @param matRegionsOption: write MATERIAL REGIONS data as well?
      For matRegion only output have frameList or fieldNames be empty lists.
    @param matRegionsNameFilter: This is for material name conversion. For
      example if there are region variants that shall be treated as one. I.e.
      if SEDIMENT, SEDIMENT_CAVE and SEDIMENT_UCUT shall all be treated as
      SEDIMENT.

      If this is a function it will be called as a conversion function
      taking the original material name from the odb as first argument and
      the section type as second argument and returning the material name it
      should be treated as.

      Example: To take everything up to the first "_" or "-" use as
      matRegionsNameFilter:
       >>> lambda matname,secType: matname.split("_")[0].split("-")[0]

      That's equivalent to:
       >>> def convert(matname,secType):
       >>>     return matname.split("_")[0].split("-")[0]
       >>> odbToVtkGridBox(..., matRegionsNameFilter=convert, ...)

    @param matNamesList: Optional list determinining the material name to
      material number assignment. Material numbers in the output will be
      assigned in the order determined by the list starting with 1.

      The specified material names must match the final material names as
      returned from the odb and possibly modified by the nameFilter function
      if applicable.

      Typically these are the uppercase material names as specified
      in the Abaqus input file of the job.

      Example:
       >>> matNamesList=["HOST", "ULTRAMAFIC", ...],

      Note: The resulting ..._matList.csv file will contain all material names
      listed with this argument plus any additional material name that may
      appear in the odb appended to the end.

    @param getExtraFields: Optional generator function to supply more fields
      from sources other than odb fields. To be used in conjuction with the
      outputFields argument. Will be called with some keyword arguments.

      Currently implemented arguments are stepNb, frameNb, points, odb.
      points is a L{MeshUnstructuredPoints} object identifying the topology
      for the fields to be supplied. odb is the
      L{OdbReaderInterpolateToPoints
      <bae.odb_03.odbReader.OdbReaderInterpolateToPoints>} objects that is used
      to get all the data from the odb. This can be used to get all sorts of
      data other than field data from the odb. The odb attribute contains
      topological data of the points as well: through the
      L{bae.mesh_01.InterpolMeshToPoints} object odb.interpolation.

      B{Important:} The argument list in the definition of getExtraFields must
      always end with I{**dummy} because in future versions of this module more
      arguments might be passed to getExtraFields. See L{Converter_CombinedVtk}
      for an example.

      The generator getExtraFields yields L{Field} objects that can then be
      used in the outputFields argument. See L{Converter_CombinedVtk} for an
      example.

      B{Note:} In a derived subclass calling the super-class getExtraFields
      use keyword arguments --at least for all but the first three arguments!
      And consider the following:
      getExtraFields is a generator function. Its result acts as an iterator.
      Therefore after calling the super-class getExtraFields method you have
      to iterate over its result and pass those items to the caller through
      the yield command. You can't just return or yield the result of the
      super-class' getExtraFields method.

    @param outputFields: a list of (fieldName, fieldType, composer, isConst)
      -tuples. If given then output will only contain fields listed here.

      FieldName will be the name in the vtk file (visible for example in
      Paraview but not in voxler), fieldType can be "scalar", "vector",
      "tensor" or "str". fieldType can be None or arbitrary if composer is
      a string. composer describes how to generate the result field.
      composer can be a string or a function.

      If composer is a string then it will be taken as field name. The field
      with that name is taken from the odb without modification. fieldType will
      be ignored. This is to simply take fields from the odb as they are and
      possibly rename them.

      If composer is a function then it will be called once for each output
      frame with the following arguments: fieldClass, fieldFromOdbDict, grid.

      fieldClass is the suitable type for the result field derived from
      L{bae.field_01.Field}.

      fieldFromOdbDict is a dictionary {fieldName : data} containing the
      fields as they would occur in the output if the outputFields argument
      would have been ommited: Values read from the odb as listed in the
      fieldNames argument. Arguments fieldNamesMaximize and relativeFieldNames
      are recongized as well and have their usual effect, in this case the
      corresponding fields in fieldFromOdbDict are maximized or taken relative.
      The matRegionsOption argument results in a field named "Material" being
      available in the fieldFromOdbDict-dictionary.

      grid is a L{bae.mesh_01.MeshUnstructuredPoints} object. (gridPoints
      argument)

      The optional (default=False) fourth item isConst can be used to identify
      this field as being constant over all time frames. See
      L{FramesFieldsOnUnstructuredPoints} for what happens with constant fields
      depending on the csvOutputType argument.

      Example:
       >>> def composerLogP(cls, dat, mesh):
       >>>     return cls(log10(x*1000+1) for x in dat["SDV3"])
       >>> def composerS1(cls, dat, mesh):
       >>>     return cls(-x for x in dat["S_MIN_PRINCIPAL"])
       >>> ...
       >>> odbToCsvPoints(...
       >>>   fieldNames=["SDV3", "S_MIN_PRINCIPAL"],
       >>>   matRegionsOption=True,
       >>>   outputFields=[
       >>>     ("logP", "scalar", composerLogP),
       >>>     ("RER", "scalar", "SDV6"),
       >>>     ("S1", "scalar", composerS1),
       >>>     ("xyz", "vector", lambda cls,dat,msh: cls(msh.getPointsIter())),
       >>>     ("Material", "scalar", "Material"),
       >>>     ],
       >>>   ...)

      Note how the "Material" field (i.e. the matcode) is transferred to the
      output without change. And SDV6 is only renamed to RER.

      Note: To pass values between different calls to the composer function(s)
      you can use global variables. Or give your composer an additional default
      argument with a list (or other container object) as default value.

      Note: If the Parameter outputFields is given then you have to have
      "Material" and "ptName" fields listed here otherwise those fields will
      not be in the output even if you specified the material option. For some
      output types the ptName field is required.

    @param gridPoints: Which points do you want results for? A list of [x,y,z]-
      lists.
      Or a list of [x,y,z,name]-lists in case of csvOutputType ==
      "CsvPerFieldRowFrColPtName". The name item will then be used for the name
      column in the output.

    @param gridName: Unique name for the gridPoints to distinguish the
      interpolation pickle file.

    @param frameList: frameList is a list of (step number, frame number list)-
      tuples. Example: [ (2, 24), (2, 27) ]

      In case fieldNamesMaximize is given and fieldMaximizeAfter is not given
      or None then this list is treated slightly different: The first frame in
      this list is the last frame not to be considered at all. No output is
      generated for this first frame! For subsequent frames output will be
      generated.
      Example: [ (2, 20), (2, 24), (2, 28) ] will generate output for frame
      24 and frame 28. Not for frame 20! The results for frame 24 will include
      maximum RER from frame 21 to frame 24 (both included). For frame 28 we'll
      have maximum RER from frame 25 to frame 28 (both included). Provided
      that RER (i.e. SDV12 for example) is listed in fieldNamesMaximize.

    @param filterEmpty: Don't write data line for points that have no data
      (i.e. lie outside the mesh or don't have the requested output quantity
      in any of the requested frames)
      Note: ignored for csvOutputType "OneCsvRowFrColPtField" and
      "OneCsvRowPtFldColFr"

    @param stepFrameToName: Defines how to build frame names in the csv file
      names or in the csv header line from a (stepNb, frameNb)-tuple.
      Might be a function (stepNb, frameNb)-tuple --> frameName or a template
      string (like the default "F%d%03d").
      Can also be the string "fromPostData" in this case the frameName attribute
      from the associated _postData.py file is used.

      For restart analyses use something like this (We only have step 2 and 3,
      the latter sometimes referred to as 1):
       >>> def stepFrameToName(stepNb, frameNb):
       >>>     if stepNb==2:
       >>>         return "F2%03d" % frameNb
       >>>     else:
       >>>         return "F3%03d" % frameNb
       >>>
       >>> odbToVtkGridBox( ..., stepFrameToName=stepFrameToName, ...)

      For old-style frame only output you may use this lambda expression:
       >>> ... stepFrameToName=lambda s,f: ('%03d' % f), ...

    @param outputDirName: Directory name for output files. This directory will
      be automatically created if it does not yet exist. It's strongly
      suggested to include the gridbox name in this outputDirName.
    @param fileNamePrefix: Base filename for output files.
      Usually projectPrefix_runVersion_seqVersion_gridName
      Example: "Perse2015_R01_Q01_boxD3_h5"
    @param projectPrefix: Needed for the interpolation data filename.
    @param meshVersion: Needed for the interpolation data filename.

    @param csvOutputType: Do you want one csv file per frame or all in one,
      and in which order? The following names are composed of three parts:
      First "OneCsv" (for all) or "CsvPerFrame" (many csv files, one for each
      frame) or "CsvPerField" (you guess it). Secondly "RowPt" or "RowFr":
      we'll have one row for each point or one row for each frame. The third
      and last part describes what we have the columns for.

       - "OneCsvRowPtColFldFr": one file for the whole data.
         Columns are: x,y,z, constant fields if given..., then
         field 1 frame 1, field 1 frame 2, field 1 frame 3, ...,
         field 2 frame 1, field 2 frame 2, field 2 frame 3, ...
         ... more fields ...
         There is exactly one row per point.
       - "CsvPerFrameRowPtColFld": one file per frame.
         Columns are: x,y,z, then constant fields, then remaining fields
         There is one row per data point and one file per odb frame.
         This format somehow corresponds to our usual vtk-gridbox data.
       - "CsvPerFieldRowPtColFr": one file per field.
         Columns are: x,y,z, then one column per frame
         There is exactly one row per point and one file per output field.
         This is the traditional style of pointsDataCohsv.py. It's
         identical to the "OneCsvRowPtColFldFr"-type in case of only one
         field being asked for.
       - "CsvPerFieldRowFrColPtName": one file per field.
         Columns are: frame name then one column for each point having the
         name in the column header. One row per frame.
         The gridPoints parameter must contain a fourth column stating the
         name for name column in the output.
       - "OneCsvRowFrColPtField": one file for the whole data.
         Columns are: frame, x,y,z, field1, field2,... x,y,z, field1,
         field2,... I.e. one row per frame and different points in columns.
         For each point we have the points coordinates first and then the
         field data.
       - "OneCsvRowPtFldColFr": one file for the whole data.
         Columns are: x,y,z, var, frame1, frame2, frame3, ...
         The var-column contains the field name. Order of rows is all points
         field 1 then all points field 2 and so on. Vector and tensor
         components come on succesive rows for the same point. I.e.::
           x,  y,  z,  var,  F2001, F2002, ...
           x1, y1, z1, logP, 0.0,   0.0,  ...
           x2, y2, z2, logP, 0.0,   0.1,  ...
           x3, y3, z3, logP, 0.0,   0.2,  ...
           x1, y1, z1, U_1,  0.0,  -0.1,  ...
           x1, y1, z1, U_2,  0.0,   0.2,  ...
           x1, y1, z1, U_3,  0.0,  -0.3,  ...
           x2, y2, z2, U_1,  0.0,  -0.1,  ...
           x2, y2, z2, U_2,  0.0,   0.2,  ...
           x2, y2, z2, U_3,  0.0,  -0.3,  ...
           ...

    @param raiseMissingData: If True (default) then raise an exception if
      data could not be found in the odb or external file. Otherwise just
      write an error message and finish.

    @param abaqusVersion: Version string of Abaqus, something like "6.13-2".
       Note: separators must be exactly as in the example: dot and dash.
       For new style use something like "abq2018", ignoring the hot-fix
       postfix e.g. "hf4".
       If not specified defaults to what the 'abaqus' shell command
       invokes.
    @param interpolationDataDir: .interpol files will be in that subfolder

    @returns: 0 on success, 1 or other on error. This value can be used
       as return code of the script in compliance with common practice.
    """

    runner = OdbToCsvPoints(
        odbPath, fieldNames, frameList,
        gridName, gridPoints,
        outputDirName, fileNamePrefix,
        projectPrefix, meshVersion,
        fieldNamesMaximize=fieldNamesMaximize,
        fieldMaximizeAfter=fieldMaximizeAfter,
        relativeFieldNames=relativeFieldNames,
        relativeToFrame=relativeToFrame,
        matRegionsOption=matRegionsOption,
        matRegionsNameFilter=matRegionsNameFilter,
        matNamesList=matNamesList,
        outputFields=outputFields,
        getExtraFields=getExtraFields,
        csvOutputType=csvOutputType,
        filterEmpty=filterEmpty,
        stepFrameToName=stepFrameToName,
        raiseMissingData=raiseMissingData,
        abaqusVersion=abaqusVersion,
        interpolationDataDir=interpolationDataDir,
        )
    retCode = runner.run()
    return retCode


#------------------------------------------------------------------------------
# odbElCentroidsToCsvPoints
#------------------------------------------------------------------------------
class OdbElCentroidsToCsvPoints(OdbToCsvPoints):
    """Read field data from the specified odb at element centroids
    and write the data to csv files.

    Don't use this class directly, its interface may change.
    Use L{odbElCentroidsToCsvPoints} instead.
    """
    def __init__(
            self,
            odbPath, fieldNames, frameList,
            boundingBox=None, elems='ALL_COHS',
            outputDirName=".", fileNamePrefix="PRJXXXX_RXX_COHS",
            fieldNamesMaximize=None, fieldMaximizeAfter=None,
            relativeFieldNames=None, relativeToFrame=(2,0),
            matRegionsOption=False, matRegionsNameFilter=None,
            matNamesList=None,
            outputFields=None,
            getExtraFields=None,
            csvOutputType="CsvPerFrameRowPtColFld",
            filterEmpty=False,
            stepFrameToName="F%d%03d",
            raiseMissingData=True,
            abaqusVersion=None,
            ):
        """Constructor does only basic argument checking and stores all
        arguments in instance variables
        """
        OdbToPointsData.__init__(
            self, fieldNames=fieldNames,
            fieldNamesMaximize=fieldNamesMaximize,
            fieldMaximizeAfter=fieldMaximizeAfter,
            relativeFieldNames=relativeFieldNames,
            relativeToFrame=relativeToFrame)

        # store arguments as instance variables
        # Note: fieldNames already processed and stored by
        # OdbToPointsData.__init__()

        # make a copy: fieldMaximizeAfter-option might alter the list...
        self.frameList = list(frameList)

        self.storeArgs(
            locals(),
            "odbPath",
            "boundingBox", "elems",
            "outputDirName", "fileNamePrefix",
            "fieldNamesMaximize", "fieldMaximizeAfter",
            "relativeFieldNames", "relativeToFrame",
            "matRegionsOption", "matRegionsNameFilter",
            "matNamesList",
            "getExtraFields",
            "outputFields",
            "csvOutputType",
            "filterEmpty",
            "stepFrameToName",
            "raiseMissingData",
            "abaqusVersion",
            )

    def initOutput(self):
        """Prepare everything for the output of results
         1. call L{OdbToPointsData.initOutput}
         2. compose output file names
         3. preprocess stepFrameToName
         4. initialize output grid and output containers
         5. initialize self.elsetNamesField field if "elsetName" in fieldNames
        """
        # check output directory and process outputFields
        OdbToPointsData.initOutput(self)

        # output file names from outputDirName and fileNamePrefix
        if self.csvOutputType.startswith("OneCsv"):
            outputFileNameTemplate = os.path.join(
                self.outputDirName, self.fileNamePrefix+".csv")
        else:
            outputFileNameTemplate = os.path.join(
                self.outputDirName, self.fileNamePrefix+"_%s.csv")

        # initialize centroids as output points
        recognizedBlocks = set(("NODE", "ELEMENT"))
        if isinstance(self.elems, (basestring, list)):
            recognizedBlocks.add("ELSET")
        msg("Reading model data from the odb.")
        model = self.odb.getAbqModel(recognizedBlocks=recognizedBlocks)
        msg("Read model: %s." % model)
        self.odb.initOutputAtCentroids(
            model=model, elems=self.elems, boundingBox=self.boundingBox)

        # initialize the output
        self.output = FramesFieldsOnUnstructuredPoints(
            points=self.odb.interpolation.grid)
        self.output.initOutput(
            outputFileNameTemplate,
            outputType=self.csvOutputType, filterEmpty=self.filterEmpty,
            frameIdToName=self.prepareStepFrameToName())

        # initialize Container for constant fields
        self.constFieldsContainer = OrderedDict()

        # create and store elsetName field
        if ("elsetName" in self.fieldNames
            and isinstance(self.elems, list)):
            elemToElset = dict()
            for eln in self.elems:
                elemToElset.update((e,eln) for e in model.elset[eln])
            elsetNamesField = createFieldObject(
                "elsetName", "point", "str", initArgs=[[
                    elemToElset[el]
                    for el in self.odb.interpolation.elemList
                ],])
            self.constFieldsContainer[elsetNamesField.fieldName] \
                = elsetNamesField

        # create and store element field
        if "element" in self.fieldNames:
            self.constFieldsContainer["element"] = createFieldObject(
                "element", "point", "scalar",
                initArgs=[self.odb.interpolation.elemList,])

#------------------------------------------------------------------------------
def odbElCentroidsToCsvPoints(
        odbPath, fieldNames, frameList,
        boundingBox=None, elems='ALL_COHS',
        outputDirName=".", fileNamePrefix="PRJXXXX_RXX_COHS",
        fieldNamesMaximize=None, fieldMaximizeAfter=None,
        relativeFieldNames=None, relativeToFrame=(2,0),
        matRegionsOption=False, matRegionsNameFilter=None,
        matNamesList=None,
        outputFields=None,
        getExtraFields=None,
        csvOutputType="CsvPerFrameRowPtColFld",
        filterEmpty=False,
        stepFrameToName="F%d%03d",
        raiseMissingData=True,
        abaqusVersion=None,
        projectPrefix="NotUsed", meshVersion="NotUsed",
        ):
    """Read field data from the specified odb at element centroids
    and write the data to csv files.

    Usage:
     >>> odbPath = config.odbPath
     >>> outputDirName = "CSV/Cadia2016033_R01_COHSV"
     >>> fileNamePrefix = "Cadia2016033_R01_COHSV"
     >>> fieldNames = ["PST", "RER"]
     >>> frameList = [(2, i) for i in range(10, 20)]
     >>>
     >>> odbElCentroidsToCsvPoints(
     >>>     odbPath,
     >>>     fieldNames=fieldNames, frameList=frameList,
     >>>     boundingBox=None, elems="COHSV",
     >>>     outputDirName=outputDirName, fileNamePrefix=fileNamePrefix,
     >>>     csvOutputType="CsvPerFieldRowPtColFr",
     >>>     )

    B{Option fieldNamesMaximize:} For selected fields take the maximum value of
    all intermediate time frames. Data for fields in fieldNamesMaximize for a
    particular frame listed in frameList will be read from one frame after the
    previous entry in frameList up the mentioned frame itself. The
    maximum for all those frames will be computed and stored in an output file
    for the mentioned frame in frameList together with all other requested
    fields from this frame. Then again only fields in fieldNamesMaximize will
    be read from the next frame (i.e. frameList[1][1]+1) up to the next frame
    listed in frameList yielding output for this latter frame (i.e. the third
    in the frameList).

    If fieldMaximizeAfter is given as well then reading of fields listed in
    fieldNamesMaximize starts at the frame following after fieldMaximizeAfter.
    All other fields are read for the first frame in frameList and the first
    output will be generated for this first frame in frameList as well.

    B{For example} fieldNames=["LogP", "RER"], fieldNamesMaximize=["RER"],
    fieldMaximizeAfter=(2,4), frameList=[(2,7), (2,10), (2,13)]:
    RER will be read for frame 5,6,7 and LogP will be read for frame 7 only.
    The first output file will be for frame 7 including LogP for frame 7 and
    the maximum of RER for frames 5,6,7. Then RER will be read for frame 8,9
    and 10 and LogP for frame 10 only. Second output file for frame 10 contains
    LogP from frame 10 and for RER the maximum of frame 8, 9, 10. Then from
    frames 11 and 12 only LogP will be read. From frame 13 RER and LogP will
    be read. Next output for frame 13 contains LogP (fr 13) and maximum RER
    from frame 11 to 13.

    If fieldMaximizeAfter is not given (or None) then the first item of
    frameList will be used in place of fieldMaximizeAfter. There will be no
    output for this first frame in this case.

    B{For example} fieldNames=["LogP", "RER"], fieldNamesMaximize=["RER"],
    fieldMaximizeAfter=None, frameList=[(2,7), (2,10), (2,13)]:
    RER will be read for frame 8,9,10 and LogP will be read for frame 10 only.
    The first output file will be for frame 10 including LogP for frame 10 and
    the maximum of RER for frames 8,9,10. Then RER will be read for frame 11,12
    and 13 and LogP for frame 13 only. Second output file for frame 13 contains
    LogP from frame 13 and for RER the maximum of frame 11, 12, 13.

    Note that if fieldNamesMaximize is specified then all time frames up to and
    including fieldMaximizeAfter or --if not specified-- the first in frameList
    will not be considered in any way. There will be no output file for the
    first frame in the frameList if fieldMaximizeAfter is not specified.

    B{Option relativeFieldNames:} For certain fields (e.g. displacement) output
    values relative to a reference frame / reference field.

    B{Option matRegionsOption:} Add material region number field to output.
    The material code field is named "Material".

    If you want the material region number only then choose a single
    arbitrary frame for the frameList parameter and make fieldNames an empty
    list and of course matRegionsOption=True.

    B{Option outputFields:} If the argument outputFields is given then only
    fields in this argument will appear in the output files. The fields listed
    in the argument fieldNames will be read from the odb and can be used to
    compute output fields.

    Please use keyword arguments at least for all optional arguments as their
    order might change in future (sub-)versions.

    @param odbPath: Complete path to the odb including .odb extension. Can be
      relative or absolute.
    @param outputDirName: Vtk files will end up in this folder. Can be
      a relative or an absolute path.
    @param fileNamePrefix: For result files and preliminary data files.
      Example: "Perse2015_R01_boxD3_h5"
    @param boundingBox: L{BoundingBox<bae.misc_01.BoundingBox>}-object (or
      [[xmin,ymin,zmin],[xmax,ymax,zmax]]). Consider only elements whose
      centroid lies within this box. If None then all elements will be
      considered.
    @param elems: Name of the elset to be exported. Defaults to
      "COHSV". If it's a list then consider all elements from all elsets in
      this list. A *set* of element labels is accepted as well.

    @param fieldNames: List of field names (in capital letters) to be exported
      to the output csv file.
      Must include all fields in argument fieldNamesMaximize.
      Example: ["UR_1","UR_2","UR_3","SDV4","S_MIN_PRINCIPAL","SDV12"]

      Note that "S_MIN_PRINCIPAL" is renamed to "S" in the output and
      "SDV_XXXXX" is renamed to "XXXXX".

      Add an item "elsetName" to get a field of elset names (from those given
      in elems). The elsets should not overlap otherwise the result is
      arbitrary in overlaps.

      Add an item "element" to get the element label (number).

    @param fieldNamesMaximize: Optional list of field names to take the maximum
      value of intermediate frames of. For maxRER. Example: ["SDV12"]
    @param fieldMaximizeAfter: If given then the frameList argument again lists
      all frames where output is generated for, including the very first.
      Specify a (stepNr, frameNr)-tuple for fieldMaximizeAfter.
      For the first output frame (the first item in the frameList-argument) data
      from frames *after* this --i.e. starting with (stepNr, frameNr+1)-- will
      be considered.

    @param relativeFieldNames: optional. If you want fields
      -i.e. displacements- exported relative to a certain reference frame
      then specify a list of field names. E.g. ["U"]
    @param relativeToFrame: reference frame as
      (step number, frame number)-tuple. Defaults to (2, 0)

    @param matRegionsOption: write MATERIAL REGIONS data as well?
      For matRegion only output have frameList or fieldNames be empty lists.
    @param matRegionsNameFilter: This is for material name conversion. For
      example if there are region variants that shall be treated as one. I.e.
      if SEDIMENT, SEDIMENT_CAVE and SEDIMENT_UCUT shall all be treated as
      SEDIMENT.

      If this is a function it will be called as a conversion function
      taking the original material name from the odb as first argument and
      the section type as second argument and returning the material name it
      should be treated as.

      Example: To take everything up to the first "_" or "-" use as
      matRegionsNameFilter:
       >>> lambda matname,secType: matname.split("_")[0].split("-")[0]

      That's equivalent to:
       >>> def convert(matname,secType):
       >>>     return matname.split("_")[0].split("-")[0]
       >>> odbToVtkGridBox(..., matRegionsNameFilter=convert, ...)

    @param matNamesList: Optional list determinining the material name to
      material number assignment. Material numbers in the output will be
      assigned in the order determined by the list starting with 1.

      The specified material names must match the final material names as
      returned from the odb and possibly modified by the nameFilter function
      if applicable.

      Typically these are the uppercase material names as specified
      in the Abaqus input file of the job.

      Example:
       >>> matNamesList=["HOST", "ULTRAMAFIC", ...],

      Note: The resulting ..._matList.csv file will contain all material names
      listed with this argument plus any additional material name that may
      appear in the odb appended to the end.

    @param getExtraFields: Optional generator function to supply more fields
      from sources other than odb fields. To be used in conjuction with the
      outputFields argument. Will be called with some keyword arguments.

      Currently implemented arguments are stepNb, frameNb, points, odb.
      points is a L{MeshUnstructuredPoints} object identifying the topology
      for the fields to be supplied. odb is the
      L{OdbReaderInterpolateToPoints
      <bae.odb_03.odbReader.OdbReaderInterpolateToPoints>} objects that is used
      to get all the data from the odb. This can be used to get all sorts of
      data other than field data from the odb. The odb attribute contains
      topological data of the points as well: through the
      L{bae.mesh_01.InterpolMeshToCentroids} object odb.interpolation and
      odb.elemList.

      B{Important:} The argument list in the definition of getExtraFields must
      always end with I{**dummy} because in future versions of this module more
      arguments might be passed to getExtraFields. See L{Converter_CombinedVtk}
      for an example.

      The generator getExtraFields yields L{Field} objects that can then be
      used in the outputFields argument. See L{Converter_CombinedVtk} for an
      example.

      B{Note:} In a derived subclass calling the super-class getExtraFields
      use keyword arguments --at least for all but the first three arguments!
      And consider the following:
      getExtraFields is a generator function. Its result acts as an iterator.
      Therefore after calling the super-class getExtraFields method you have
      to iterate over its result and pass those items to the caller through
      the yield command. You can't just return or yield the result of the
      super-class' getExtraFields method.

    @param outputFields: a list of (fieldName, fieldType, composer, isConst)
      -tuples. If given then output will only contain fields listed here.

      FieldName will be the name in the vtk file (visible for example in
      Paraview but not in voxler), fieldType can be "scalar", "vector",
      "tensor" or "str". fieldType can be None or arbitrary if composer is
      a string. composer describes how to generate the result field.
      composer can be a string or a function.

      If composer is a string then it will be taken as field name. The field
      with that name is taken from the odb without modification. fieldType will
      be ignored. This is to simply take fields from the odb as they are and
      possibly rename them.

      If composer is a function then it will be called once for each output
      frame with the following arguments: fieldClass, fieldFromOdbDict, grid.

      fieldClass is the suitable type for the result field derived from
      L{bae.field_01.Field}.

      fieldFromOdbDict is a dictionary {fieldName : data} containing the
      fields as they would occur in the output if the outputFields argument
      would have been ommited: Values read from the odb as listed in the
      fieldNames argument. Arguments fieldNamesMaximize and relativeFieldNames
      are recongized as well and have their usual effect, in this case the
      corresponding fields in fieldFromOdbDict are maximized or taken relative.
      The matRegionsOption argument results in a field named "Material" being
      available in the fieldFromOdbDict-dictionary.

      grid is a L{bae.mesh_01.MeshUnstructuredPoints} object. (gridPoints
      argument)

      The optional (default=False) fourth item isConst can be used to identify
      this field as being constant over all time frames. See
      L{FramesFieldsOnUnstructuredPoints} for what happens with constant fields
      depending on the csvOutputType argument.

      Example:
       >>> def composerLogP(cls, dat, mesh):
       >>>     return cls(log10(x*1000+1) for x in dat["SDV3"])
       >>> def composerS1(cls, dat, mesh):
       >>>     return cls(-x for x in dat["S_MIN_PRINCIPAL"])
       >>> ...
       >>> odbToCsvPoints(...
       >>>   fieldNames=["SDV3", "S_MIN_PRINCIPAL"],
       >>>   matRegionsOption=True,
       >>>   outputFields=[
       >>>     ("logP", "scalar", composerLogP),
       >>>     ("RER", "scalar", "SDV6"),
       >>>     ("S1", "scalar", composerS1),
       >>>     ("xyz", "vector",
       >>>      lambda cls,dat,msh: cls(msh.getPointsIter()), True),
       >>>     ("Material", "scalar", "Material", True),
       >>>     ],
       >>>   ...)

      Note how the "Material" field (i.e. the matcode) is transferred to the
      output without change. And SDV6 is only renamed to RER.

      Note: To pass values between different calls to the composer function(s)
      you can use global variables. Or give your composer an additional default
      argument with a list (or other container object) as default value.

      Note: If the Parameter outputFields is given then you have to have
      "Material" and "ptName" fields listed here otherwise those fields will
      not be in the output even if you specified the material option. For some
      output types the ptName field is required.

    @param frameList: List of output frames as (step number, frame number)
      tuples. Example: [ (2, 24), (2, 27) ]

      In case fieldNamesMaximize is given this list is treated slightly
      different: The first frame in this list is the first frame to
      consider at all. No output is generated for this first frame! For
      subsequent frames output will be generated.
      Example: [ (2, 20), (2, 24), (2, 27) ] will generate output for frame
      24 and frame 27. Not for frame 20! The results for frame 24 will include
      maximum RER from frame 20 to frame 24 (both included). For frame 27 we'll
      have maximum RER from frame 25 to frame 27 (both included). Provided
      that RER (i.e. SDV12 for example) is listed in fieldNamesMaximize.

    @param csvOutputType:
      See L{bae.vtk_02.FramesFieldsOnUnstructuredPoints.initOutput}
       - "OneCsvRowPtColFldFr": one file for the whole data.
         Columns are: x,y,z, constant fields if given..., then
         field 1 frame 1, field 1 frame 2, field 1 frame 3, ...,
         field 2 frame 1, field 2 frame 2, field 2 frame 3, ...
         ... more fields ...
         There is exactly one row per point.
       - "CsvPerFrameRowPtColFld": one file per frame.
         Columns are: x,y,z, then constant fields, then remaining fields
         There is one row per data point and one file per odb frame.
         This format somehow corresponds to our usual vtk-gridbox data.
       - "CsvPerFieldRowPtColFr": one file per field.
         Columns are: x,y,z, then one column per frame
         There is exactly one row per point and one file per output field.
         This is the traditional style of pointsDataCohsv.py. It's
         identical to the "OneCsvRowPtColFldFr"-type in case of only one
         field being asked for.
       - "OneCsvRowFrColPtField": one file for the whole data.
         Columns are: frame, x,y,z, field1, field2,... x,y,z, field1,
         field2,... I.e. one row per frame and different points in columns.
         For each point we have the points coordinates first and then the
         field data.

    @param stepFrameToName: Defines how to build frame names in the csv file
      names or in the csv header line from a (stepNb, frameNb)-tuple.
      Might be a function (stepNb, frameNb)-tuple --> frameName or a template
      string. Default if not given, i.e. None: "F%d%03d".

      For restart analyses use something like this (We only have step 2 and 3,
      the latter sometimes referred to as 1):
       >>> def stepFrameToName(stepNb, frameNb):
       >>>     if stepNb==2:
       >>>         return "F2%03d" % frameNb
       >>>     else:
       >>>         return "F3%03d" % frameNb
       >>>
       >>> odbToVtkGridBox( ..., stepFrameToName=stepFrameToName, ...)

    @param raiseMissingData: If True (default) then raise an exception if
      data could not be found in the odb or external file. Otherwise just
      write an error message and finish.

    @param abaqusVersion: Version string of Abaqus, something like "6.13-2".
       Note: separators must be exactly as in the example: dot and dash.
       For new style use something like "abq2018", ignoring the hot-fix
       postfix e.g. "hf4".
       If not specified defaults to what the 'abaqus' shell command
       invokes.

    @returns: 0 on success, 1 or other on error. This value can be used
       as return code of the script in compliance with common practice.
    """
    runner = OdbElCentroidsToCsvPoints(
        odbPath, fieldNames, frameList,
        boundingBox=boundingBox, elems=elems,
        outputDirName=outputDirName, fileNamePrefix=fileNamePrefix,
        fieldNamesMaximize=fieldNamesMaximize,
        fieldMaximizeAfter=fieldMaximizeAfter,
        relativeFieldNames=relativeFieldNames,
        relativeToFrame=relativeToFrame,
        matRegionsOption=matRegionsOption,
        matRegionsNameFilter=matRegionsNameFilter,
        matNamesList=matNamesList,
        outputFields=outputFields,
        getExtraFields=getExtraFields,
        csvOutputType=csvOutputType,
        filterEmpty=filterEmpty,
        stepFrameToName=stepFrameToName,
        raiseMissingData=raiseMissingData,
        abaqusVersion=abaqusVersion,
        )
    retCode = runner.run()
    return retCode


#------------------------------------------------------------------------------
# odbFaceCentroidsToCsvPoints
#------------------------------------------------------------------------------
class OdbFaceCentroidsToCsvPoints(OdbToPointsData):
    """Read field data from the specified odb at element face centroids
    of a specified surface and write the data to csv files.

    See L{OdbVariableFaceCentroidsToCsvPoints} if the surface changes over time.

    Don't use this class directly, its interface may change.
    Use L{odbFaceCentroidsToCsvPoints} instead.
    """
    outputFieldType = "point"

    def __init__(
            self,
            odbPath, fieldNames, frameList, surface,
            outputDirName=".", fileNamePrefix="PRJXXXX_RXX_COHS",
            fieldNamesMaximize=None, fieldMaximizeAfter=None,
            relativeFieldNames=None, relativeToFrame=(2,0),
            matRegionsOption=False, matRegionsNameFilter=None,
            matNamesList=None,
            outputFields=None,
            getExtraFields=None,
            csvOutputType="CsvPerFrameRowPtColFld",
            filterEmpty=False,
            stepFrameToName="F%d%03d",
            raiseMissingData=True,
            abaqusVersion=None,
            ):
        """Constructor does only basic argument checking and stores all
        arguments in instance variables
        """
        OdbToPointsData.__init__(
            self, fieldNames=fieldNames,
            fieldNamesMaximize=fieldNamesMaximize,
            fieldMaximizeAfter=fieldMaximizeAfter,
            relativeFieldNames=relativeFieldNames,
            relativeToFrame=relativeToFrame)

        # convert single surface name to list with this name as single item
        if isinstance(surface, basestring):
            surface = [surface,]

        # store arguments as instance variables
        # Note: fieldNames already processed and stored by
        # OdbToPointsData.__init__()

        # make a copy: fieldMaximizeAfter-option might alter the list...
        self.frameList = list(frameList)

        self.storeArgs(
            locals(),
            "odbPath", "surface",
            "outputDirName", "fileNamePrefix",
            "fieldNamesMaximize", "fieldMaximizeAfter",
            "relativeFieldNames", "relativeToFrame",
            "matRegionsOption", "matRegionsNameFilter",
            "matNamesList",
            "getExtraFields",
            "outputFields",
            "csvOutputType",
            "filterEmpty",
            "stepFrameToName",
            "raiseMissingData",
            "abaqusVersion",
            )

    def initOutput(self):
        """Prepare everything for the output of results
         1. call L{OdbToPointsData.initOutput}
         2. compose output file names
         3. preprocess stepFrameToName
         4. initialize output grid and output containers
         5. initialize self.surfaceNamesField field if "surfaceName"
            in fieldNames
        """
        # check output directory and process outputFields
        OdbToPointsData.initOutput(self)

        # output file names from outputDirName and fileNamePrefix
        if self.csvOutputType.startswith("OneCsv"):
            outputFileNameTemplate = os.path.join(
                self.outputDirName, self.fileNamePrefix+".csv")
        else:
            outputFileNameTemplate = os.path.join(
                self.outputDirName, self.fileNamePrefix+"_%s.csv")

        # initialize centroids as output points
        recognizedBlocks = set(("NODE", "ELEMENT"))
        if isinstance(self.surface, list):
            recognizedBlocks.update(("ELSET", "SURFACE"))
        msg("Reading model data from the odb.")
        model = self.odb.getAbqModel(recognizedBlocks=recognizedBlocks)
        msg("Read model: %s." % model)

        msg("Preparing surface(s) and output locations...")
        faceEl = defaultdict(set)
        for surfName in self.surface:
            surf = ElemFaceSurface().updateFromModel(model, surfName)
            for faceId, elems in surf.faceEl.iteritems():
                faceEl[faceId].update(elems)
        faceEl.default_factory = None  # convert to ordinary dict
        self.odb.initOutputAtFaceCentroids(model=model, faceEl=faceEl)
        # Maybe the previous chunk could be replaced by...
        # surf = ElemFaceSurface()
        # for surfName in self.surface:
        #     surf.updateFromModel(model, surfName)
        # self.odb.initOutputAtFaceCentroids(model=model, faceEl=surf.faceEl)

        # initialize the output
        self.output = FramesFieldsOnUnstructuredPoints(
            points=self.odb.interpolation.grid)
        self.output.initOutput(
            outputFileNameTemplate,
            outputType=self.csvOutputType, filterEmpty=self.filterEmpty,
            frameIdToName=self.prepareStepFrameToName())

        # initialize Container for constant fields
        self.constFieldsContainer = OrderedDict()

        # create and store surfaceName field
        if ("surfaceName" in self.fieldNames
            and isinstance(self.surface, list)):
            msg("Preparing surfaceName field...")
            elemFaceIdToSurfName = dict()
            for surfName in self.surface:
                surf = ElemFaceSurface().updateFromModel(model, surfName)
                elemFaceIdToSurfName.update(
                    ((e, faceId), surfName)
                    for faceId, elems in surf.faceEl.iteritems()
                    for e in elems)
            surfaceNamesField = createFieldObject(
                "surfaceName", "point", "str", initArgs=[[
                    elemFaceIdToSurfName[elemFace]
                    for elemFace in self.odb.interpolation.elemFaceIdList],])
            self.constFieldsContainer[surfaceNamesField.fieldName] \
                = surfaceNamesField

#------------------------------------------------------------------------------
def odbFaceCentroidsToCsvPoints(
        odbPath, fieldNames, frameList, surface,
        outputDirName=".", fileNamePrefix="PRJXXXX_RXX_COHS",
        fieldNamesMaximize=None, fieldMaximizeAfter=None,
        relativeFieldNames=None, relativeToFrame=(2,0),
        matRegionsOption=False, matRegionsNameFilter=None,
        matNamesList=None,
        outputFields=None,
        getExtraFields=None,
        csvOutputType="CsvPerFrameRowPtColFld",
        filterEmpty=False,
        stepFrameToName="F%d%03d",
        raiseMissingData=True,
        abaqusVersion=None,
        ):
    """Read field data from the specified odb at element face centroids
    of a specified surface and write the data to csv files.

    Usage:
     >>> odbPath = config.odbPath
     >>> outputDirName = "CSV/Cadia2016033_R01_subsidence"
     >>> fileNamePrefix = "Cadia2016033_R01_subsidence"
     >>> fieldNames = ["U_3", "logP"]
     >>> frameList = [(2, i) for i in range(10, 20)]
     >>>
     >>> odbFaceCentroidsToCsvPoints(
     >>>     odbPath, fieldNames, frameList, "TOP_SURF",
     >>>     outputDirName=outputDirName, fileNamePrefix=fileNamePrefix,
     >>>     csvOutputType="CsvPerFieldRowPtColFr",
     >>>     )

    ...or if the surface is not defined in the odb but in some external Abaqus
    input file:
     >>> model = Model().read("mysurf.inp")
     >>> surf = ElemFaceSurface().updateFromModel(model, "MYSURF")
     >>> odbFaceCentroidsToCsvPoints(
     >>>     odbPath, fieldNames, frameList,
     >>>     outputDirName, fileNamePrefix,
     >>>     surface=surf.faceEl,
     >>>     csvOutputType="CsvPerFieldRowPtColFr",
     >>>     )

    B{Option fieldNamesMaximize:} Maximize selected fields between output
    frames. Data will be read from the specified odb from the first frame
    specified in frameList to its last including all intermediate frames not
    just those contained in frameList. For frames not specified in frameList
    read only fields in fieldNamesMaximize. For frames listed in frameList
    -except for the first- generate output vtks. Take the maximum value since
    the previous output frame -or the first frame in frameList- for fields
    listed in fieldNamesMaximize.

    B{Option relativeFieldNames:} For certain fields (e.g. displacement) output
    values relative to a reference frame / reference field.

    B{Option matRegionsOption:} Add material region number field to every output
    vtk file. If you want the material region number only then choose a single
    arbitrary frame for the frameList parameter and make fieldNames an empty
    list and of course matRegionsOption=True.

    B{Option outputFields:} If the argument outputFields is given then only
    fields in this argument will appear in the output files. The fields listed
    in the argument fieldNames will be read from the odb and can be used to
    compute output fields.

    Please use keyword arguments at least for all optional arguments as their
    order might change in future (sub-)versions.

    @param odbPath: Complete path to the odb including .odb extension. Can be
      relative or absolute.
    @param outputDirName: Vtk files will end up in this folder. Can be
      a relative or an absolute path.
    @param fileNamePrefix: For result files and preliminary data files.
      Example: "Perse2015_R01_boxD3_h5"
    @param surface: A surface name defined in the odb or a list of such names.
      or a dict: {"S1": set((4346, 4768, 6821, ...)), ...}

    @param fieldNames: List of field names to be exported to the output vtk.
      Must include all fields in argument fieldNamesMaximize.
      Example: ["UR_1","UR_2","UR_3","SDV4","S_MIN_PRINCIPAL","SDV12"]

      Note that "S_MIN_PRINCIPAL" is renamed to "S" in the output and
      "SDV_XXXXX" is renamed to "XXXXX".

      Add an item "surfaceName" to get a field of surface names (from those
      given in the surface argument).
    @param fieldNamesMaximize: Optional list of field names to take the maximum
      value of intermediate frames of. For maxRER. Example: ["SDV12"]
    @param fieldMaximizeAfter: If given then the frameList argument again lists
      all frames where output is generated for, including the very first.
      Specify a (stepNr, frameNr)-tuple for fieldMaximizeAfter.
      For the first output frame (the first item in the frameList-argument) data
      from frames *after* this --i.e. starting with (stepNr, frameNr+1)-- will
      be considered.
    @param relativeFieldNames: Optional. From fields listed in this argument
      the value of the time frame relativeToFrame is being subtracted before
      export. E.g. for relative displacements (relative to the equilibrium
      step) let relativeFieldNames=["U_1", "U_2", "U_3"] and
      relativeToFrame=(2.0).
    @param relativeToFrame: Reference time frame for the relativeFieldNames
      option. A (stepNb, frameNb)-tuple.
    @param matRegionsOption: Boolean flag: Write MATERIAL REGIONS data as well?
      If True then a field named "Material" will be generated.
    @param matRegionsNameFilter: This is for material name conversion. For
      example if there are region variants that shall be treated as one. I.e.
      if SEDIMENT, SEDIMENT_CAVE and SEDIMENT_UCUT shall all be treated as
      SEDIMENT. For further details see the corresponding parameter of
      L{odbToVtkGridBox}.
    @param matNamesList: Optional list determinining the material name to
      material number assignment. Material numbers in the output will be
      assigned in the order determined by the list starting with 1.
      For further details see the corresponding parameter of L{odbToVtkGridBox}.

    @param getExtraFields: Optional generator function to supply more fields
      from sources other than odb fields. To be used in conjuction with the
      outputFields argument. Will be called with some keyword arguments.

      Currently implemented arguments are stepNb, frameNb, points, odb.
      points is a L{MeshUnstructuredPoints} object identifying the topology
      for the fields to be supplied. odb is the
      L{OdbReaderInterpolateToPoints
      <bae.odb_03.odbReader.OdbReaderInterpolateToPoints>} objects that is used
      to get all the data from the odb. This can be used to get all sorts of
      data other than field data from the odb. The odb attribute contains
      topological data of the points as well: through the
      L{bae.mesh_01.InterpolMeshToFaceCentroids} object odb.interpolation.

      B{Important:} The argument list in the definition of getExtraFields must
      always end with I{**dummy} because in future versions of this module more
      arguments might be passed to getExtraFields. See L{Converter_CombinedVtk}
      for an example.

      The generator getExtraFields yields L{Field} objects that can then be
      used in the outputFields argument. See L{Converter_CombinedVtk} for an
      example.

      B{Note:} In a derived subclass calling the super-class getExtraFields
      use keyword arguments --at least for all but the first three arguments!
      And consider the following:
      getExtraFields is a generator function. Its result acts as an iterator.
      Therefore after calling the super-class getExtraFields method you have
      to iterate over its result and pass those items to the caller through
      the yield command. You can't just return or yield the result of the
      super-class' getExtraFields method.

    @param outputFields: a list of (fieldName, fieldType, composer)-tuples.
      If given then output will only contain fields listed here.
      For further details see the corresponding parameter of L{odbToVtkGridBox}.

    @param frameList: List of output frames as (step number, frame number)
      tuples. Example: [ (2, 24), (2, 27) ]

      In case fieldNamesMaximize is given this list is treated slightly
      different: The first frame in this list is the first frame to
      consider at all. No output is generated for this first frame! For
      subsequent frames output will be generated.
      Example: [ (2, 20), (2, 24), (2, 27) ] will generate output for frame
      24 and frame 27. Not for frame 20! The results for frame 24 will include
      maximum RER from frame 20 to frame 24 (both included). For frame 27 we'll
      have maximum RER from frame 25 to frame 27 (both included). Provided
      that RER (i.e. SDV12 for example) is listed in fieldNamesMaximize.

    @param csvOutputType:
      See L{bae.vtk_02.FramesFieldsOnUnstructuredPoints.initOutput}
       - "OneCsvRowPtColFldFr": one file for the whole data.
         Columns are: x,y,z, constant fields if given..., then
         field 1 frame 1, field 1 frame 2, field 1 frame 3, ...,
         field 2 frame 1, field 2 frame 2, field 2 frame 3, ...
         ... more fields ...
         There is exactly one row per point.
       - "CsvPerFrameRowPtColFld": one file per frame.
         Columns are: x,y,z, then constant fields, then remaining fields
         There is one row per data point and one file per odb frame.
         This format somehow corresponds to our usual vtk-gridbox data.
       - "CsvPerFieldRowPtColFr": one file per field.
         Columns are: x,y,z, then one column per frame
         There is exactly one row per point and one file per output field.
         This is the traditional style of pointsDataCohsv.py. It's
         identical to the "OneCsvRowPtColFldFr"-type in case of only one
         field being asked for.
       - "OneCsvRowFrColPtField": one file for the whole data.
         Columns are: frame, x,y,z, field1, field2,... x,y,z, field1,
         field2,... I.e. one row per frame and different points in columns.
         For each point we have the points coordinates first and then the
         field data.

    @param stepFrameToName: Defines how to build frame names in the csv file
      names or in the csv header line from a (stepNb, frameNb)-tuple.
      Might be a function (stepNb, frameNb)-tuple --> frameName or a template
      string. Default if not given, i.e. None: "F%d%03d".

      For restart analyses use something like this (We only have step 2 and 3,
      the latter sometimes referred to as 1):
       >>> def stepFrameToName(stepNb, frameNb):
       >>>     if stepNb==2:
       >>>         return "F2%03d" % frameNb
       >>>     else:
       >>>         return "F3%03d" % frameNb
       >>>
       >>> odbToVtkGridBox( ..., stepFrameToName=stepFrameToName, ...)

    @param raiseMissingData: If True (default) then raise an exception if
      data could not be found in the odb or external file. Otherwise just
      write an error message and finish.

    @param abaqusVersion: Version string of Abaqus, something like "6.13-2".
       Note: separators must be exactly as in the example: dot and dash.
       For new style use something like "abq2018", ignoring the hot-fix
       postfix e.g. "hf4".
       If not specified defaults to what the 'abaqus' shell command
       invokes.

    @returns: 0 on success, 1 or other on error. This value can be used
       as return code of the script in compliance with common practice.
    """
    runner = OdbFaceCentroidsToCsvPoints(
        odbPath, fieldNames, frameList,
        surface=surface,
        outputDirName=outputDirName, fileNamePrefix=fileNamePrefix,
        fieldNamesMaximize=fieldNamesMaximize,
        fieldMaximizeAfter=fieldMaximizeAfter,
        relativeFieldNames=relativeFieldNames,
        relativeToFrame=relativeToFrame,
        matRegionsOption=matRegionsOption,
        matRegionsNameFilter=matRegionsNameFilter,
        matNamesList=matNamesList,
        outputFields=outputFields,
        getExtraFields=getExtraFields,
        csvOutputType=csvOutputType,
        filterEmpty=filterEmpty,
        stepFrameToName=stepFrameToName,
        raiseMissingData=raiseMissingData,
        abaqusVersion=abaqusVersion,
        )
    retCode = runner.run()
    return retCode


#------------------------------------------------------------------------------
# odbVariableFaceCentroidsToCsvPoints
#------------------------------------------------------------------------------
class OdbVariableFaceCentroidsToCsvPoints(OdbToPointsData):
    """Read field data from the specified odb at element face centroids
    of a specified surface that differs from one time frame to the next.
    Write the data to a set of csv files, one per time frame.

    See L{OdbFaceCentroidsToCsvPoints} if the surface does not change over time.

    Don't use this class directly, its interface may change.
    Use L{odbFaceCentroidsToCsvPoints} instead.
    """
    outputFieldType = "point"

    def __init__(
            self,
            odbPath, fieldNames, frameList,
            initialSurface="top", seqSetsNameExcav=None, seqSetsNameFill=None,
            outputDirName=".", fileNamePrefix="PRJXXXX_RXX_topSurf",

            # maximize option not implemented:
            # fieldNamesMaximize=None, fieldMaximizeAfter=None,

            relativeFieldNames=None, relativeToFrame=(2,0),

            matRegionsOption=False, matRegionsNameFilter=None,
            matNamesList=None,

            outputFields=None,
            getExtraFields=None,

            # not applicable because the points might change:
            # csvOutputType="CsvPerFrameRowPtColFld",
            # filterEmpty=False,

            stepFrameToName="F%d%03d",
            raiseMissingData=True,
            abaqusVersion=None,
            ):
        """Constructor does only basic argument checking and stores all
        arguments in instance variables
        """
        OdbToPointsData.__init__(self, fieldNames=fieldNames)

        # store arguments as instance variables
        # Note: fieldNames already processed and stored by
        # OdbToPointsData.__init__()
        self.storeArgs(
            locals(),
            "odbPath", "frameList",
            "initialSurface", "seqSetsNameExcav", "seqSetsNameFill",
            "outputDirName", "fileNamePrefix",
            # "fieldNamesMaximize", "fieldMaximizeAfter",
            "relativeFieldNames", "relativeToFrame",
            "matRegionsOption", "matRegionsNameFilter",
            "matNamesList",
            "getExtraFields",
            "outputFields",
            "stepFrameToName",
            "raiseMissingData",
            "abaqusVersion",
            )

    def initOutput(self):
        """Prepare everything for the output of results
         1. call L{OdbToPointsData.initOutput}
         2. compose output file names
         3. preprocess stepFrameToName
         4. initialize output grid and output containers
         5. initialize current surface self.surf
        """
        # check output directory and process outputFields
        OdbToPointsData.initOutput(self)

        # output file names from outputDirName and fileNamePrefix
        self.outputFileNameTemplate = os.path.join(
            self.outputDirName, self.fileNamePrefix+"_%s.csv")

        # get the mesh and elsets
        msg("Reading model data from the odb.")
        self.model = self.odb.getAbqModel(recognizedBlocks="NODE,ELEMENT,ELSET")
        msg("Read model: %s." % self.model)

        # initialize the initial surface
        # note: self.model is a quadratic mesh as abq_model_02.Model object
        # self.surf.mesh is a linear mesh (with consistent element numbers)
        # of type mesh_01.Mesh
        if self.initialSurface == "top":
            msg("Preparing initial surface as top surface of the model...")
            self.surf = self.model.getTopSurface()
        elif self.initialSurface == "none":
            msg("Starting with a nonexisting initial surface...")
            ticker = MsgTicker("Stitching the mesh: ...%s")
            from bae.utils.splitting_01 import rejoinSplitMesh
            linMesh = Mesh()
            linMesh.nodeCoords = self.model.nodeCoords
            linMesh.updateElems(dict(
                (el, self.model.elNodes[el][:4])
                for el in self.model.shapeEl["TET"]), "TET_L")
            ticker.msg("copied all tet elements.")
            linMesh.updateElems(dict(
                (el, self.model.elNodes[el])
                for el in self.model.shapeEl["WEDGE"]), "WEDGE_L")
            ticker.msg("copied all cohesive elements.")
            rejoinSplitMesh(linMesh)
            ticker.msg("stitched the cohesive gaps.")
            self.surf = ElemFaceSurface(model=linMesh)
            ticker.msg("ready.")
        else:
            raise NotImplementedError(
                "Unrecognized value for the argument initialSurface: '%s'."
                " Expecting 'top' or 'none'."
                % self.initialSurface)

        # convert mesh to Abaqus model and then also store elsets
        self.surf.mesh = self.surf.mesh.asAbqModel()
        self.surf.mesh.elset.update(self.model.elset)

        # initialize elsetsExcav and elsetsFill iterators yielding lists of
        # elset names from seqSetsNameExcav and seqSetsNameFill
        if self.seqSetsNameExcav:
            if self.seqSetsNameExcav not in self.odb.postData:
                raise KeyError(
                    "Sequence elsets '%s' not defined in associated postData.py"
                    " file for odbPath=%s."
                    % (self.seqSetsNameExcav, self.odbPath))
            elsetsExcav = sequenceElsetsIter(
                self.odb.postData, frameList=self.frameList,
                seqElsets=self.seqSetsNameExcav, incremental=True)
        else:
            elsetsExcav = None
        if self.seqSetsNameFill:
            if self.seqSetsNameFill not in self.odb.postData:
                raise KeyError(
                    "Sequence elsets '%s' not defined in associated postData.py"
                    " file for odbPath=%s."
                    % (self.seqSetsNameFill, self.odbPath))
            elsetsFill = sequenceElsetsIter(
                self.odb.postData, frameList=self.frameList,
                seqElsets=self.seqSetsNameFill, incremental=True)
        else:
            elsetsFill = None

        # initialize self.seqSurfIter
        self.seqSurfIter = self.surf.updateFromSequence(
            model=self.surf.mesh,
            elsetsExcav=elsetsExcav, elsetsFill=elsetsFill)

    def initOutputForFrame(self, stepNb, frameNb):
        """prepare everything for the current time step

        Initialize self.fieldsContainer

        @param stepNb: not used for this variant
        @param frameNb: not used for this variant
        """

        # update surface, init output points
        surf = self.seqSurfIter.next()  # surf == self.surf, btw.
        msg("...updated surface: %s" % surf, debugLevel=10)

        # temporarily take the quad mesh for the interpolation
        # otherwise there'll be problems interpolating to the face centres
        linMesh = self.surf.mesh
        self.surf.mesh = self.model
        self.odb.initOutputAtFaceCentroids(surface=surf)
        self.surf.mesh = linMesh

        # initialize the output
        self.output = FramesFieldsOnUnstructuredPoints(
            points=self.odb.interpolation.grid)
        self.output.initOutput(
            self.outputFileNameTemplate,
            outputType="CsvPerFrameRowPtColFld",
            frameIdToName=self.prepareStepFrameToName())

        # initialize temporal data container
        self.fieldsContainer = OrderedDict()

    def storeField(self, field):
        """Store the given field for output. Will be called once for each field
        and multiple times for each frame.
        """
        self.fieldsContainer[field.fieldName] = field

    def storeFrame(self, stepNb, frameNb):
        """Store frame data to output. This will be called once for each frame
        after all fields have been processed by self.storeField().

        Processes self.outputFields by means of
        L{OdbToPointsData.processOutputFields}().
        """
        # compose output fields if outputFields argument is given
        if self.outputFields:
            # process self.outputFields
            newFldContainer = OrderedDict(
                (field.fieldName, field)
                for field in self.processOutputFields(self.fieldsContainer))
            self.fieldsContainer = newFldContainer
            del newFldContainer

        self.output.updateFields(*self.fieldsContainer.itervalues())
        del self.fieldsContainer
        self.output.storeFrame((stepNb, frameNb))
        self.output.closeCsvOutput()

    def run(self):
        """Main function that does the actual job. Calls all the other
        functions, loops over frames and fields.

        We need a special version here for our points (i.e. the topgraphy /
        the grid) is changing from step to step

        @returns: 0 on success, 1 or other on error. This value can be used
        as return code of the script in compliance with common practice.
        """

        retCode = 0  # ok, no error

        if self.relativeFieldNames:
            relativeFieldNames = set(self.relativeFieldNames)
            refStep, refFrame = self.relativeToFrame
        else:
            relativeFieldNames = None

        self.openOdb()
        self.initOutput()

        # get the field data from the odb frame by frame
        msg("Now starting to iterate over output frames.")
        stopNow = False
        for stepNb, frameNb in self.frameList:

            self.initOutputForFrame(stepNb, frameNb)

            # mat regions preparation
            if self.matRegionsOption:
                self.storeField(self.getMatRegions())

            # other fields from odb
            for fieldName in self.fieldNames:
                interpolationType = self.fieldInterpolType[fieldName]
                iptText = ""
                if interpolationType=="const":
                    iptText = " interpolating piecewise constant"

                msg("Reading field %s%s, step %d, frame %d."
                    % (fieldName, iptText, stepNb, frameNb))
                try:
                    field = self.odb.getFieldFromOdb(
                        stepNb, frameNb, fieldName,
                        interpolationType=interpolationType)
                except ValueError as exc:
                    if self.raiseMissingData:
                        raise
                    else:
                        msg("ERROR: Field %s for odb frame (%s, %d) not"
                            " available. Ending this loop now.\nError: %s"
                            % (fieldName, stepNb, frameNb, exc))
                        stopNow = True
                        retCode = 1
                        break

                if relativeFieldNames and fieldName in relativeFieldNames:
                    msg("Reading reference field %s%s, from step %d, frame %d."
                        % (fieldName, iptText, refStep, refFrame))
                    refField = self.odb.getFieldFromOdb(
                        refStep, refFrame, fieldName,
                        interpolationType=interpolationType)
                    if field.dataType == "scalar":
                        field = type(field)(
                            (a-b) for a, b in izip(field, refField))
                    else:
                        field = type(field)(
                            [(a[i]-b[i]) for i in range(field.nbComponents)]
                            for a, b in izip(field, refField))

                # store the field for output
                self.storeField(field)
                msg("Stored field %s in output." % fieldName, debugLevel=10)

            # after for fieldName ... loop
            if stopNow:
                break

            if self.getExtraFields:
                # fields from external sources, e.g. from Cavesim vtk
                getExtraFieldsArgs = dict(
                    stepNb=stepNb, frameNb=frameNb,
                    points=self.output.mesh, odb=self.odb)
                for fld in self.getExtraFields(**getExtraFieldsArgs):
                    self.storeField(fld)
                del getExtraFieldsArgs

            # self.storeFrame evaluates self.outputFields and writes files
            self.storeFrame(stepNb, frameNb)

        # finally write material names to csv file
        if self.matRegionsOption:
            self.writeMatNamesListToFile()

        return retCode

#------------------------------------------------------------------------------
def odbVariableFaceCentroidsToCsvPoints(
        odbPath, fieldNames, frameList,
        initialSurface="top", seqSetsNameExcav=None, seqSetsNameFill=None,
        outputDirName=".", fileNamePrefix="PRJXXXX_RXX_topSurf",

        # maximize option not implemented:
        # fieldNamesMaximize=None, fieldMaximizeAfter=None,

        relativeFieldNames=None, relativeToFrame=(2,0),

        matRegionsOption=False, matRegionsNameFilter=None,
        matNamesList=None,

        outputFields=None,
        getExtraFields=None,

        stepFrameToName="F%d%03d",
        raiseMissingData=True,
        abaqusVersion=None,
        ):
    """Read field data from the specified odb at element face centroids
    of a surface changing from time step to time step.

    Write the data to a set of csv files, one per time frame.

    See L{odbFaceCentroidsToCsvPoints} if the surface does not change over time.

    Usage:
     >>> odbPath = config.odbPath
     >>> outputDirName = "CSV/Cadia2016033_R01_subsidence"
     >>> fileNamePrefix = "Cadia2016033_R01_subsidence"
     >>> fieldNames = ["U_3", "logP"]
     >>> frameList = [(2, i) for i in range(10, 20)]
     >>>
     >>> odbVariableFaceCentroidsToCsvPoints(
     >>>     odbPath, fieldNames=fieldNames, frameList=frameList,
     >>>     initialSurface="top", seqSetsNameExcav="seqElsetsPit",
     >>>     outputDirName=outputDirName, fileNamePrefix=fileNamePrefix,
     >>>     )

    ...or for an underground excavation:
     >>> odbVariableFaceCentroidsToCsvPoints(
     >>>     odbPath, fieldNames=fieldNames, frameList=frameList,
     >>>     outputDirName=outputDirName, fileNamePrefix=fileNamePrefix,
     >>>     initialSurface="none", seqSetsNameExcav="seqElsetsCrusher",
     >>>     )

    B{Option matRegionsOption:} Add material region number field to every output
    vtk file. If you want the material region number only then choose a single
    arbitrary frame for the frameList parameter and make fieldNames an empty
    list and of course matRegionsOption=True.

    B{Option outputFields:} If the argument outputFields is given then only
    fields in this argument will appear in the output files. The fields listed
    in the argument fieldNames will be read from the odb and can be used to
    compute output fields.

    Please use keyword arguments at least for all optional arguments as their
    order might change in future (sub-)versions.

    @param odbPath: Complete path to the odb including .odb extension. Can be
      relative or absolute.
    @param outputDirName: Vtk files will end up in this folder. Can be
      a relative or an absolute path.
    @param fileNamePrefix: For result files and preliminary data files.
      Example: "Perse2015_R01_boxD3_h5"

    @param fieldNames: List of field names to be exported to the output vtk.
      Must include all fields in argument fieldNamesMaximize.
      Example: ["UR_1","UR_2","UR_3","SDV4","S_MIN_PRINCIPAL","SDV12"]

      Note that "S_MIN_PRINCIPAL" is renamed to "S" in the output and
      "SDV_XXXXX" is renamed to "XXXXX".

      Add an item "surfaceName" to get a field of surface names (from those
      given in the surface argument).
    @param matRegionsOption: Boolean flag: Write MATERIAL REGIONS data as well?
      If True then a field named "Material" will be generated.
    @param matRegionsNameFilter: This is for material name conversion. For
      example if there are region variants that shall be treated as one. I.e.
      if SEDIMENT, SEDIMENT_CAVE and SEDIMENT_UCUT shall all be treated as
      SEDIMENT. For further details see the corresponding parameter of
      L{odbToVtkGridBox}.
    @param matNamesList: Optional list determinining the material name to
      material number assignment. Material numbers in the output will be
      assigned in the order determined by the list starting with 1.
      For further details see the corresponding parameter of L{odbToVtkGridBox}.

    @param getExtraFields: Optional generator function to supply more fields
      from sources other than odb fields. To be used in conjuction with the
      outputFields argument. Will be called with some keyword arguments.

      Currently implemented arguments are stepNb, frameNb, points, odb.
      points is a L{MeshUnstructuredPoints} object identifying the topology
      for the fields to be supplied. odb is the
      L{OdbReaderInterpolateToPoints
      <bae.odb_03.odbReader.OdbReaderInterpolateToPoints>} objects that is used
      to get all the data from the odb. This can be used to get all sorts of
      data other than field data from the odb. The odb attribute contains
      topological data of the points as well: through the
      L{bae.mesh_01.InterpolMeshToFaceCentroids} object odb.interpolation.

      B{Important:} The argument list in the definition of getExtraFields must
      always end with I{**dummy} because in future versions of this module more
      arguments might be passed to getExtraFields. See L{Converter_CombinedVtk}
      for an example.

      The generator getExtraFields yields L{Field} objects that can then be
      used in the outputFields argument. See L{Converter_CombinedVtk} for an
      example.

      B{Note:} In a derived subclass calling the super-class getExtraFields
      use keyword arguments --at least for all but the first three arguments!
      And consider the following:
      getExtraFields is a generator function. Its result acts as an iterator.
      Therefore after calling the super-class getExtraFields method you have
      to iterate over its result and pass those items to the caller through
      the yield command. You can't just return or yield the result of the
      super-class' getExtraFields method.

    @param outputFields: a list of (fieldName, fieldType, composer)-tuples.
      If given then output will only contain fields listed here.
      For further details see the corresponding parameter of L{odbToVtkGridBox}.

    @param relativeFieldNames: optional. If you want fields
      -i.e. displacements- exported relative to a certain reference frame
      then specify a list of field names. E.g. ["U"]
    @param relativeToFrame: reference frame as
      (step number, frame number)-tuple. Defaults to (2, 0)

    @param frameList: List of output frames as (step number, frame number)
      tuples. Example: [ (2, 24), (2, 27) ]

      In case fieldNamesMaximize is given this list is treated slightly
      different: The first frame in this list is the first frame to
      consider at all. No output is generated for this first frame! For
      subsequent frames output will be generated.
      Example: [ (2, 20), (2, 24), (2, 27) ] will generate output for frame
      24 and frame 27. Not for frame 20! The results for frame 24 will include
      maximum RER from frame 20 to frame 24 (both included). For frame 27 we'll
      have maximum RER from frame 25 to frame 27 (both included). Provided
      that RER (i.e. SDV12 for example) is listed in fieldNamesMaximize.

    @param initialSurface: Can be "top" for the top surface of the model.
      Can be "none" for no surface initially, e.g. to develop excavation voids.

    @param seqSetsNameExcav: A string specifying the sequence elsets object in
      the postData file for regions to be extracted step by step.
      This string must match the variable name in the postData file. I.e.
      to reference the sequence from the example postdata file sketched
      below you would specify seqElsets="seqElsets".
      See L{bae.makeHistory_03.sequenceElsetsIter} for further explanation.

    @param seqSetsNameFill: Same as seqSetsNameExcav for regions to be added
      step by step. Like dumps if you're after the top surface or refill if
      you look at an underground excavation.

    @param stepFrameToName: Defines how to build frame names in the csv file
      names or in the csv header line from a (stepNb, frameNb)-tuple.
      Might be a function (stepNb, frameNb)-tuple --> frameName or a template
      string. Default if not given, i.e. None: "F%d%03d".

      For restart analyses use something like this (We only have step 2 and 3,
      the latter sometimes referred to as 1):
       >>> def stepFrameToName(stepNb, frameNb):
       >>>     if stepNb==2:
       >>>         return "F2%03d" % frameNb
       >>>     else:
       >>>         return "F3%03d" % frameNb
       >>>
       >>> odbToVtkGridBox( ..., stepFrameToName=stepFrameToName, ...)

    @param raiseMissingData: If True (default) then raise an exception if
      data could not be found in the odb or external file. Otherwise just
      write an error message and finish.

    @param abaqusVersion: Version string of Abaqus, something like "6.13-2".
       Note: separators must be exactly as in the example: dot and dash.
       For new style use something like "abq2018", ignoring the hot-fix
       postfix e.g. "hf4".
       If not specified defaults to what the 'abaqus' shell command
       invokes.

    @returns: 0 on success, 1 or other on error. This value can be used
       as return code of the script in compliance with common practice.
    """
    runner = OdbVariableFaceCentroidsToCsvPoints(
        odbPath, fieldNames, frameList,
        initialSurface=initialSurface,
        seqSetsNameExcav=seqSetsNameExcav, seqSetsNameFill=seqSetsNameFill,
        outputDirName=outputDirName, fileNamePrefix=fileNamePrefix,
        relativeFieldNames=relativeFieldNames, relativeToFrame=relativeToFrame,
        matRegionsOption=matRegionsOption,
        matRegionsNameFilter=matRegionsNameFilter,
        matNamesList=matNamesList,
        outputFields=outputFields,
        getExtraFields=getExtraFields,
        stepFrameToName=stepFrameToName,
        raiseMissingData=raiseMissingData,
        abaqusVersion=abaqusVersion,
        )
    retCode = runner.run()
    return retCode


# ### CONSTRUCTION SITE AHEAD ###
# #------------------------------------------------------------------------------
# # odbToSurfaceCsv
# #------------------------------------------------------------------------------
# class OdbToSurfaceCsv(object):  # or derive from (OdbToPointsData)?

#     def __init__(
#             self,
#             odbPath, fieldNames, frameList,
#             gridName, gridData,

#             outputDirName=".", fileNamePrefix="PRJXXXX_RXX_topSurf",

#             # maximize option not implemented:
#             # fieldNamesMaximize=None, fieldMaximizeAfter=None,

#             relativeFieldNames=None,
#             relativeToFrame=(2,0), relativeToOdbPath=None,

#             outputFields=None,
#             getExtraFields=None,

#             # not applicable because the points might change:
#             # csvOutputType="CsvPerFrameRowPtColFld",
#             # filterEmpty=False,

#             stepFrameToName="F%d%03d",
#             raiseMissingData=True,
#             varNameState = None,
#             abaqusVersion=None):
#         """
#         @param varNameState: Usually not required. If not specified then this
#           will be determined from the postData file. For cases this is not
#           possible you can specify something like "SDV8" here.
#         """

#         self.odb = OdbReader(odbPath)
#         if varNameState:
#             self.varNameState = varNameState
#         else:
#             self.varNameState = self.odb.postData["varNameState"]

#     @staticmethod
#     def stateFilterFuncDefault2021(gpVals):
#         "Return False for elements to be removed."
#         return all(x in (1, 4) for x in gpVals)

#     def extractFreeSurface(
#             self, odbStep, odbFrame,
#             stateFilterFunc=stateFilterFuncDefault2021):

#         self.odb.openOdb()

#         # get mesh if not yet there
#         if not hasattr(self, mesh):
#             self.mesh = self.odb.getAbqModel()

#         # read status and determine elements to be removed
#         stateFld = self.odb.getFieldFromOdb(
#             odbStep, odbFrame, self.varNameState)
#         removeElems = [
#             e for e, gpVals in stateFld
#             if stateFilterFunc]

#         # update mesh: remove elements based on state
#         self.mesh.removeElems(removeElems)
#         topSurf = self.mesh.getTopSurface(surfZTol=0.0)

#         return topSurf

#     def run(self):

#         pass


# def odbToSurfaceCsv(
#         odbPath, fieldNames, frameList,
#         gridName, gridData,

#         outputDirName=".", fileNamePrefix="PRJXXXX_RXX_topSurf",

#         # maximize option not implemented:
#         # fieldNamesMaximize=None, fieldMaximizeAfter=None,

#         relativeFieldNames=None,
#         relativeToFrame=(2,0), relativeToOdbPath=None,
            
#         outputFields=None,
#         getExtraFields=None,

#         # not applicable because the points might change:
#         # csvOutputType="CsvPerFrameRowPtColFld",
#         # filterEmpty=False,

#         stepFrameToName="F%d%03d",
#         raiseMissingData=True,
#         abaqusVersion=None):
#     """Read field data from the specified odb at the current free top surface
#     of the model.
#     Write the data to a set of csv files, one per time frame.
#     """
#     pass


# ### CONSTRUCTION SITE END ###

#------------------------------------------------------------------------------
# odbToDriveClosureCsv
#------------------------------------------------------------------------------

def odbToDriveClosureCsv(
        odbPath, fieldNames, frameList,

        gridName, gridDataPath,
        outputDirName, fileNamePrefix,
        projectPrefix, meshVersion,

        displacementFieldName="U", stressField="S",
        relativeFieldNames=None, relativeToFrame=(2,0),

        csvOutputType="OneCsvRowPtColFldFr",
        filterEmpty=False,
        stepFrameToName="F%d%03d",
        abaqusVersion=None,
        interpolationDataDir=".",
        ):
    """Read interpolated field data from the specified odb at specified points
    and write the data to csv files.

    @param odbPath: Complete path to the odb including .odb extension. Can be
      relative or absolute.

    @param fieldNames: A list of one or more of the following (case sensitive):
      "closure_horiz", "closure_vert",
      "Snormal_left", "Snormal_right", "Snormal_top", "Snormal_floor",
      "Sshear_left", "Sshear_right", "Sshear_top", "Sshear_floor",
      "point_left", "point_right", "point_top", "point_floor"

      Note: The order of the requested fields might change according to the
      order given above.

    @param relativeFieldNames: optional. If you want fields
      -i.e. displacements- exported relative to a certain reference frame
      then specify a list of field names. E.g. ["U"]
    @param relativeToFrame: reference frame as
      (step number, frame number)-tuple. Defaults to (2, 0)

    @param displacementFieldName: field name to be used as displacement vector
      for calculating the closure.
    @param stressField: field name to be used as stress tensor for calculating
      the stress output values.

    @param frameList: frameList is a list of (step number, frame number list)-
      tuples.

      The step number is the number in the step name of the step in the odb.
      I.e. step number 3 identifies odb step "Step-3". In the ordinary two
      steps analysis the possible values are 1 and 2. Note that if the first
      step in the odb is called "Step-4" the corresponding step number is 4.

      The frame number list contains frame numbers as in the odb. Instead of
      a frame number list there can be None which means: yield all frames for
      this step. Negative values can be used to count from the end like usual
      list indices.

      Frame numbers that do not exist in the odb are ignored without any
      warning.

    @param gridName: Unique name for the gridPoints to distinguish the
      interpolation pickle file.

    @param gridDataPath: path to the csv file containing the grid points and
      point labels. It's got to be a csv file with following columns:
      x, y, z, name. No header lines.

      The BE-Rhino-tools function exportPointNameLayerCsv can be used with
      the layer serving as name: The function creates such csv files
      with an additional fourth column for the object name -usually empty,
      which we will ignore here- and the layer name in the fifth column.
      The layer might contain a parent layer prefix for this reason anything
      up to and including a double colon "::" will be discarded.

      Note: The points must be ordered in groups of four per station. I.e. the
      function takes the first four points and works out which of them is
      left, right, above and below the drive assuming they belong to one
      cross section of the drive. Then it takes the next four doing the same.
      The name (Rhino layer) must be the same for those four points belonging
      to one station.
      It seems that Rhino exports the points in the order they have been
      created (possibly reversed).
      Left vs right is determined such that "left" is further West (neg. x)
      and "right" is further East (pos. x). Once "left" and "right" is
      established for the first cross section following cross sections
      take the previous direction from the "left" point to the "right" as
      guide. This is to prevent jumps in case the drive is oriented approx.
      East-West or for curved drives.

    @param stepFrameToName: Defines how to build frame names in the csv file
      names or in the csv header line from a (stepNb, frameNb)-tuple.
      Might be a function (stepNb, frameNb)-tuple --> frameName or a template
      string (like the default "F%d%03d").
      Can also be the string "fromPostData" in this case the frameName attribute
      from the associated _postData.py file is used.

      For restart analyses use something like this (We only have step 2 and 3,
      the latter sometimes referred to as 1):
       >>> def stepFrameToName(stepNb, frameNb):
       >>>     if stepNb==2:
       >>>         return "F2%03d" % frameNb
       >>>     else:
       >>>         return "F3%03d" % frameNb
       >>>
       >>> odbToVtkGridBox( ..., stepFrameToName=stepFrameToName, ...)

    @param outputDirName: Directory name for output files. This directory will
      be automatically created if it does not yet exist. It's strongly
      suggested to include the gridbox name in this outputDirName.

    @param fileNamePrefix: Base filename for output files.
      Usually projectPrefix_runVersion_seqVersion_gridName

    @param csvOutputType: Do you want one csv file per frame or all in one,
      and in which order? The following names are composed of three parts:
      First "OneCsv" (for all) or "CsvPerFrame" (many csv files, one for each
      frame) or "CsvPerField" (you guess it). Secondly "RowPt" or "RowFr":
      we'll have one row for each point or one row for each frame. The third
      and last part describes what we have the columns for.

       - "OneCsvRowPtColFldFr": one file for the whole data.
         Columns are: x,y,z, constant fields if given..., then
         field 1 frame 1, field 1 frame 2, field 1 frame 3, ...,
         field 2 frame 1, field 2 frame 2, field 2 frame 3, ...
         ... more fields ...
         There is exactly one row per point.
       - "CsvPerFrameRowPtColFld": one file per frame.
         Columns are: x,y,z, then constant fields, then remaining fields
         There is one row per data point and one file per odb frame.
         This format somehow corresponds to our usual vtk-gridbox data.
       - "CsvPerFieldRowPtColFr": one file per field.
         Columns are: x,y,z, then one column per frame
         There is exactly one row per point and one file per output field.
         This is the traditional style of pointsDataCohsv.py. It's
         identical to the "OneCsvRowPtColFldFr"-type in case of only one
         field being asked for.
       - "CsvPerFieldRowFrColPtName": one file per field.
         Columns are: frame name then one column for each point having the
         name in the column header. One row per frame.
         The gridPoints parameter must contain a fourth column stating the
         name for name column in the output.
       - "OneCsvRowFrColPtField": one file for the whole data.
         Columns are: frame, x,y,z, field1, field2,... x,y,z, field1,
         field2,... I.e. one row per frame and different points in columns.
         For each point we have the points coordinates first and then the
         field data.

    @param interpolationDataDir: .interpol files will be in that subfolder

    @param abaqusVersion: Version string of Abaqus, something like "6.13-2".
       Note: separators must be exactly as in the example: dot and dash.
       For new style use something like "abq2018", ignoring the hot-fix
       postfix e.g. "hf4".
       If not specified defaults to what the 'abaqus' shell command
       invokes.

    @Note: Possibly add later: option to group values: all, average, min, max
       for points with same name
    """

    # load csv : xyz, layer name
    msg("Reading grid points from %s" % gridDataPath)
    gridDataFile = csv.reader(open(gridDataPath, "rb"))
    gridPoints, names = zip(*(
        (map(float, row[:3]), row[3:])
        for row in gridDataFile))
    del gridDataFile, gridDataPath

    if len(names[0])>1 and names[0][1]:
        # we have a non-empty fifth column -> use that
        # ... and remove "parentlayer::"-prefix from layer name
        names = [x[1].rsplit("::", 1)[-1] for x in names]
    else:
        # else just take the fourth column as it is
        # (assuming that Rhino layer names only occur in the fifth column)
        names = [x[0] for x in names]
    msg("... read %d points." % len(gridPoints))

    # identify drive positions
    # prepare point list with points in order left, right, top, bottom
    lineCnt = 0
    leftToRight = [1,0,0]
    gridPointsIter = iter(gridPoints)
    namesIter = iter(names)
    gridPoints = list()  # new list to be
    names = list()
    centroids = list()
    normals = list()  # normals in order left, right, top, bottom

    while 1:

        # get next four points (assumed to belong to one cross section)
        try:
            points = [gridPointsIter.next() for i in range(4)]
        except StopIteration:
            break

        # identify and check names
        lineCnt += 4
        name = [namesIter.next() for i in range(4)]
        if any((name[0]!=x) for x in name[1:]):
            msg("ERROR: Four successive points have different names:"
                " %s, starting on line nb %d."
                "\nNote: This might indicate that the order of points is not"
                " correct. We need four points (left, right, top, bottom) for"
                " one drive cross section then four for the next cross section"
                " and so on. Anyway: For now taking the first name."
                % (",".join(name), lineCnt-3))
        names.append(name[0])

        # sort points according to z-coord
        points.sort(key=lambda x:x[2])
        # check z-coords: |left-right| / top-bottom << 1
        check = float(points[2][2]-points[1][2]) / (points[3][2]-points[0][2])
        if check<0 or check>0.25:
            msg("WARNING: Expecting top point to be considerably higher and"
                " bottom point to be considerably lower than left and right"
                " points. This is not the case for the quadruple starting at"
                " line %d. Maybe the order of points is not correct. We need"
                " four points (left, right, top, bottom) for one drive cross"
                " section then four for the next cross section and so on."
                % lineCnt-3)

        # sort out left vs right
        dir12 = norm(vector(points[1], points[2]))
        if dot(dir12, leftToRight)<0:
            points[1], points[2] = points[2], points[1]
            vector_modif_scale(dir12, -1)
        leftToRight = dir12

        # store points (left, right, top, bottom)
        gridPoints.extend(points[i] for i in (1, 2, 3, 0))
        # store centroid and four normals
        centroids.append(vector_scale(vector_sum(*points), 0.25))
        dirUpDown = norm(vector(points[3], points[0]))
        normals.extend([
            vector_scale(dir12, -1), dir12,
            vector_scale(dirUpDown, -1), dirUpDown])
    msg("Sorted %d points into left, right, top, bottom order."
        % len(gridPoints))

    #--- get data from odb and store in FramesFieldsOnUnstructuredPoints
    msg("Opening the odb %s." % odbPath)
    odb = OdbReaderInterpolateToPoints(odbPath, version=abaqusVersion)

    # check output directory
    if outputDirName and not os.path.isdir(outputDirName):
        os.makedirs(outputDirName)
        msg("Created output dir %s" % outputDirName)

    # pickle data file name
    # from interpolationDataDir and others
    interpolPickleName = os.path.join(
        interpolationDataDir, "_".join((projectPrefix, meshVersion, gridName)))

    # check interpolationDataDir directory
    if interpolationDataDir and not os.path.isdir(interpolationDataDir):
        os.makedirs(interpolationDataDir)
        msg("Created folder for interpolation data %s" % interpolationDataDir)

    # default file names from outputDirName and fileNamePrefix
    if csvOutputType.startswith("OneCsv"):
        outputFileNameTemplate = os.path.join(
            outputDirName, fileNamePrefix+".csv")
    else:
        outputFileNameTemplate = os.path.join(
            outputDirName, fileNamePrefix+"_%s.csv")

    # initialize the odb reader's output grid
    odb.initOutputPoints(
        grid=MeshUnstructuredPoints(gridPoints),
        interpolPickleName=interpolPickleName)

    # read (and store) reference data for relative fields
    refData = dict()
    if relativeFieldNames:
        stepNb, frameNb = relativeToFrame
        for fieldName in relativeFieldNames:
            msg("Reading reference field %s from step %d, frame %d."
                % (fieldName, stepNb, frameNb))
            refData[fieldName] = odb.getFieldFromOdb(
                stepNb, frameNb, fieldName)
        del stepNb, frameNb

    # convert stepFrameToName function:
    # The frameIdToName argument to vtk_02.FramesFieldsOnUnstructuredPoints.
    # ....initOutput --if its a function-- takes a single frameId argument:
    # a (stepNb, frameNb) - tuple.
    # On the other hand all bae.utils/odbToPointsData - methods
    # expect for their stepFrameToName argument a function taking two arguments:
    # stepNb, frameNb
    if callable(stepFrameToName):
        # string don't need conversion, but functions do
        stepFrameToNameTwoArgs = stepFrameToName
        def stepFrameToName(frameId):
            return stepFrameToNameTwoArgs(*frameId)

    # if argument stepFrameToName=="fromPostData": define special function
    elif stepFrameToName=="fromPostData":
        # special case, frameNames from odb.postData
        if "frameNames" not in odb.postData:
            raise ValueError(
                "frameNames not defined in associated _postData.py file for"
                " odbPath=%s." % odbPath)
        def stepFrameToName(frameId):
            frameName = odb.postData["frameNames"][frameId[0]][frameId[1]]
            msg("frameName from postData: step %d, frame %d, frameName %s"
                % (frameId[0], frameId[1], frameName), debugLevel=10)
            return frameName

    # initialize the output
    output = FramesFieldsOnUnstructuredPoints(points=centroids)
    output.initOutput(
        outputFileNameTemplate,
        outputType=csvOutputType, filterEmpty=filterEmpty,
        frameIdToName=stepFrameToName)
    output.updateFields(
        createFieldObject("ptName", "point", "str", initArgs=[names,]))
    output.storeConstantData()

    # store point coordinates if requested (in the requested order)
    sectionPointNames = ["left", "right", "top", "floor"]
    nameToSectPt = dict(
        ("point_%s" % x, i) for i, x in enumerate(sectionPointNames))
    sectionPts = [
        (x, nameToSectPt[x]) for x in fieldNames if x in nameToSectPt]
    for ptName, i in sectionPts:
        output.updateFields(
            createFieldObject(ptName, "point", "vector",
                              initArgs=[gridPoints[i::4],]))
    if sectionPts:
        output.storeConstantData()

    # sort requested fields (retain order)
    closureFieldsDict = {"closure_horiz": (0,1), "closure_vert": (2,3)}
    closureFields = [
        (x, closureFieldsDict[x], createFieldClass(x, "point", "scalar"))
        for x in fieldNames if x in closureFieldsDict]
    stressFieldsDict = dict(
        ("Snormal_%s" % x, ("normal", i))
        for i, x in enumerate(sectionPointNames))
    stressFieldsDict.update(
        ("Sshear_%s" % x, ("shear", i))
        for i, x in enumerate(sectionPointNames))
    stressFields = [
        (x, stressFieldsDict[x], createFieldClass(x, "point", "scalar"))
        for x in fieldNames if x in stressFieldsDict]

    # fieldsToRead  ... from odb
    fieldsToRead = [displacementFieldName]  # we always want closure, right?
    if stressFields:
        # we also want stress
        fieldsToRead.append(stressField)

    # get the field data from the odb frame by frame
    msg("Now starting to iterate over output frames.")
    for stepNb, frameNb in frameList:

        fieldData = dict()
        for fieldName in fieldsToRead:
            msg("Reading field %s from the odb, step %d, frame %d."
                % (fieldName, stepNb, frameNb))
            field = odb.getFieldFromOdb(
                stepNb, frameNb, fieldName)

            # subtract reference field values (e.g. displacements)
            try:
                refField = refData[fieldName]
            except KeyError:
                pass
            else:
                if field.dataType == "scalar":
                    field = type(field)(
                        (a-b) for a, b in izip(field, refField))
                else:
                    field = type(field)(
                        [(a[i]-b[i]) for i in range(field.nbComponents)]
                        for a, b in izip(field, refField))
            # store
            fieldData[fieldName] = field

        # closure
        fieldU = fieldData[displacementFieldName]
        for fieldName, (ptId1, ptId2), FieldClass in closureFields:
            # closure = dot( (u[ptId2]-u[ptId1]) , normal[ptId1] )
            # (ptId1, ptId2) is (left, right) or (top, floor)
            field = FieldClass(
                dot(vector(u1, u2), normal1)
                for u1, u2, normal1 in izip(
                        islice(fieldU, ptId1, None, 4),
                        islice(fieldU, ptId2, None, 4),
                        islice(normals, ptId1, None, 4)))
            output.updateFields(field)
        del fieldU, fieldName, FieldClass

        # stress
        if stressFields:
            # calculate all relevant values (no matter if actually needed)
            stress = {"normal": [], "shear": []}
            for S, normal in izip(fieldData[stressField], normals):
                # S components are S_11, S_22, S_33, S_12, S_13, S_23
                s = [[S[0],S[3],S[4]],
                     [S[3],S[1],S[5]],
                     [S[4],S[5],S[2]]]
                svec = mat_multvec(s, normal)
                sn = -dot(svec, normal)
                stress["normal"].append(sn)
                stress["shear"].append(
                    length(vector(svec, vector_scale(normal, -sn))))

            # now store what you need
            for fieldName, (type_, i), FieldClass in stressFields:
                field = FieldClass(islice(stress[type_], i, None, 4))
                output.updateFields(field)

        # finally store all the data
        output.storeFrame((stepNb, frameNb))

    # finalize output
    output.closeCsvOutput()


#------------------------------------------------------------------------------
# odbToDriveClosureLinesCsv
#------------------------------------------------------------------------------

def odbToDriveClosureLinesCsv(
        odbPath, fieldNames, frameList,

        gridName, gridDataPath,
        outputDirName, fileNamePrefix,
        projectPrefix, meshVersion,

        displacementFieldName="U", stressFieldName="S",
        relativeToFrame=None,

        csvOutputType="OneCsvRowPtColFldFr",

        outputDirNameGroups=None,
        csvOutputTypeGroups=None,

        filterEmpty=False,
        stepFrameToName="F%d%03d",
        abaqusVersion=None,
        interpolationDataDir=".",
        ):
    """Read interpolated field data from the specified odb at the end points of
    specified lines and write closure and support load data to csv file(s).

    There will be two (sets of) output files: One with data for each individual
    station (line). And a second with data averaged per layer, i.e data for
    each group.

    @param odbPath: Complete path to the odb including .odb extension. Can be
      relative or absolute.

    @param fieldNames: A list of one or more (case sensitive) field names.

      Field names relevant for the per-station output are: "closure",
      "point_start", "point_end", "Snormal" or "Sshear" (for averaged values)
      and combinations of "Snormal" or "Sshear" with one of the postfixes
      "_start", "_end".

      Field names relevant for the per-group output are: "closure",
      "Snormal" or "Sshear" (for averaged values) or combinations of
      "closure", "Snormal" or "Sshear" with one of the postfixes
      "_min", "_max".

      For example: ["closure", "Snormal_mean", "point_start", "point_end"]

      Note: The order of the requested fields might change: point_start and
      point_end are constant values and are being treated differently.

    @param relativeToFrame: If you want closure relative to a certain reference
      frame then simply specify the reference frame as
      (step number, frame number)-tuple.

      For proper handling of sequenced drive sections there are two additional
      options available, a simple and a flexible one. For both options the
      content of the layer column of the closure-stations csv file is used to
      defer the reference frame of the particular closure station. This string
      will be referred to as layerName in the following description.
      Additionally the frame names as specified in the ..._postData.py file
      can be used through the L{postData<bae.odb_03.OdbReader.postData>}
      attribute of the current odb. Those strings will be referrend to as
      frameNames in the following description.

      The simple option is to specify the string "frameNameFromPostData".
      In this case we try to match frame names  from the ..._postData.py file
      to the end of a given layerName.

      Note that with this option the frame names will be replaced by whatever
      the stepFrameToName argument specifies in the header of the resulting
      output csv file.

      Example: Suppose frameNames "Y2017_06", "Y2017_12", "Y2018_06" refer to
      frames (2,1), (2,2) and (2,3) according to the ..._postData.py file.
      Having the sequenced closure station lines in Rhino on layers "Y2017_06",
      "Y2017_12", "Y2018_06" we export those with the
      "export bolt lines..."-button. Layer names being composed as
      "objectName_layerName" we end up with layerNames like
      "horiz_XCUT_Y2017_06", "vert_DP_Y2017_12" in the csv file specified by
      the gridDataPath argument. The reference frames for those two layer names
      would be (2,1) and (2,2). Suppose the stepFrameToName argument is
      "F%d%03d" then the column headers will be "horiz_XCUT_F2001",
      "vert_DP_F2002".

      For the flexible option a function can be passed to this argument that
      returns the reference frame number for a particular closure station.
      This function gets the content of the layer column of the
      closure-stations csv file and the
      L{postData<bae.odb_03.OdbReader.postData>} attribute of the current odb.
      It is expected to return the associated (stepNb, frameNb)-tuple.

      IMPORTANT NOTE: To maintain compatibility with possible future versions
      of this module the function should take an arbitrary-length tuple of
      positional arguments!
      I.e. use this phrase: C{def relativeToFrame(*args)}.

      Simple example, frame number as last part of the layer name "horz_F02":
       >>> def relativeToFrame(*args):
       >>>     layerName = args[0]
       >>>     frameNb = int(layerName.rsplit("_F", 1)[-1])
       >>>     return (2, frameNb)
       >>>
       >>> # ....
       >>> odbToDriveClosureLinesCsv(
       >>>     ...,
       >>>     relativeToFrame=relativeToFrame,
       >>>     ...)

      Example, layer names contain frame name:
       >>> dateToRefFrame = dict()
       >>> def relativeToFrame(*args):
       >>>     global dateToRefFrame
       >>>     # Note: Possibly more than two args in future versions
       >>>     layerName, postData = args[:2]
       >>>
       >>>     # initialize dateToRefFrame only once at first invocation
       >>>     if not dateToRefFrame:
       >>>         frameNames = postData["frameNames"]
       >>>         dateToRefFrame.update(
       >>>             (name, (stepNb, frameNb))
       >>>             for stepNb, frameData in frameNames.iteritems()
       >>>             for frameNb, name in enumerate(frameData))
       >>>
       >>>     # slice off the frame name / date from layerName
       >>>     name = layerName.split("_", 1)[-1]
       >>>
       >>>     # look up right frame for elset
       >>>     frame = dateToRefFrame[name]
       >>>     return frame
       >>>
       >>> # ....
       >>> odbToDriveClosureLinesCsv(
       >>>     ...,
       >>>     relativeToFrame=relativeToFrame,
       >>>     ...)

      Example, layer names contain seq elset names:
       >>> postDataAttrSeqElsets = "seqElsetsGS"
       >>> seqElsetToRefFrame = dict()
       >>> def relativeToFrame(*args):
       >>>     global seqElsetToRefFrame, postDataAttrSeqElsets
       >>>     # Note: Possibly more than two args in future versions
       >>>     layerName, postData = args[:2]
       >>>
       >>>     # initialize seqElsetToRefFrame only once at first invocation
       >>>     if not seqElsetToRefFrame:
       >>>         seqElsets = postData[postDataAttrSeqElsets]
       >>>         seqElsetToRefFrame.update(
       >>>             (el, (stepNb, frameNb))
       >>>             for stepNb, frameData in seqElsets.iteritems()
       >>>             for frameNb, elsets in enumerate(frameData)
       >>>             for el in elsets)
       >>>
       >>>     # slice off the name of the sequence elset from layerName
       >>>     elset = layerName.split("_", 1)[-1]
       >>>
       >>>     # look up right frame for elset
       >>>     frame = seqElsetToRefFrame[elset]
       >>>     return frame
       >>>
       >>> # ....
       >>> odbToDriveClosureLinesCsv(
       >>>     ...,
       >>>     relativeToFrame=relativeToFrame,
       >>>     ...)

    @param displacementFieldName: field name to be used as displacement vector
      for calculating the closure.
    @param stressFieldName: field name to be used as stress tensor for
      calculating the stress output values.

    @param frameList: frameList is a list of (step number, frame number list)-
      tuples.

      The step number is the number in the step name of the step in the odb.
      I.e. step number 3 identifies odb step "Step-3". In the ordinary two
      steps analysis the possible values are 1 and 2. Note that if the first
      step in the odb is called "Step-4" the corresponding step number is 4.

      The frame number list contains frame numbers as in the odb. Instead of
      a frame number list there can be None which means: yield all frames for
      this step. Negative values can be used to count from the end like usual
      list indices.

      Frame numbers that do not exist in the odb are ignored without any
      warning.

    @param gridName: Unique name for the gridPoints to distinguish the
      interpolation pickle file.

    @param gridDataPath: Path to the csv file containing the closure stations.
      It's got to be a csv file as you get from the BE-Rhino-tools function
      "Export bolt lines as csv file".
      This csv file starts with the following header line:
      layer,x0,y0,z0,x1,y1,z1,dx,dy,dz

      Columns to follow contain the corresponding data.

    @param stepFrameToName: Defines how to build frame names in the csv file
      names or in the csv header line from a (stepNb, frameNb)-tuple.
      Might be a function (stepNb, frameNb)-tuple --> frameName or a template
      string (like the default "F%d%03d").
      Can also be the string "fromPostData" in this case the frameName attribute
      from the associated _postData.py file is used.

      For restart analyses use something like this (We only have step 2 and 3,
      the latter sometimes referred to as 1):
       >>> def stepFrameToName(stepNb, frameNb):
       >>>     if stepNb==2:
       >>>         return "F2%03d" % frameNb
       >>>     else:
       >>>         return "F3%03d" % frameNb
       >>>
       >>> odbToVtkGridBox( ..., stepFrameToName=stepFrameToName, ...)

    @param outputDirName: Directory name for output files for
      all-stations-data. This directory will be automatically created if it
      does not yet exist.
      It's strongly suggested to include the gridbox name in this outputDirName.

    @param outputDirNameGroups: Directory name for output files for
      data grouped by layer. This directory will be automatically created if it
      does not yet exist.
      It's strongly suggested to include the gridbox name in this outputDirName.
      Defaults to outputDirName.

    @param fileNamePrefix: Base filename for output files.
      Usually projectPrefix_runVersion_seqVersion_gridName. "_stations.csv" and
      "_groups.csv" will be appended for the final file name.
      Or "_stations_<frameId>.csv" and "_groups_<frameId>.csv" in case of
      parameter L{csvOutputType} referring to "per-frame"-files. See also
      parameter L{stepFrameToName}.

    @param csvOutputType: For the all-stations-data output:
      Do you want one csv file per frame or all in one, and in which order?
      The following names are composed of three parts:
      First "OneCsv" (for all) or "CsvPerFrame" (many csv files, one for each
      frame) or "CsvPerField" (you guess it). Secondly "RowPt" or "RowFr":
      we'll have one row for each point or one row for each frame. The third
      and last part describes what we have the columns for.

       - "OneCsvRowPtColFldFr": one file for the whole data.
         Columns are: x,y,z, constant fields if given..., then
         field 1 frame 1, field 1 frame 2, field 1 frame 3, ...,
         field 2 frame 1, field 2 frame 2, field 2 frame 3, ...
         ... more fields ...
         There is exactly one row per point.
       - "CsvPerFrameRowPtColFld": one file per frame.
         Columns are: x,y,z, then constant fields, then remaining fields
         There is one row per data point and one file per odb frame.
         This format somehow corresponds to our usual vtk-gridbox data.
       - "CsvPerFieldRowPtColFr": one file per field.
         Columns are: x,y,z, then one column per frame
         There is exactly one row per point and one file per output field.
         This is the traditional style of pointsDataCohsv.py. It's
         identical to the "OneCsvRowPtColFldFr"-type in case of only one
         field being asked for.
       - "CsvPerFieldRowFrColPtName": one file per field.
         Columns are: frame name then one column for each point having the
         name in the column header. One row per frame.
         The gridPoints parameter must contain a fourth column stating the
         name for name column in the output.
       - "OneCsvRowFrColPtField": one file for the whole data.
         Columns are: frame, x,y,z, field1, field2,... x,y,z, field1,
         field2,... I.e. one row per frame and different points in columns.
         For each point we have the points coordinates first and then the
         field data.

    @param csvOutputTypeGroups: For the data grouped by layer:
      Do you want one csv file per frame or all in one, and in which order?
      See parameter L{csvOutputType} for possible values.
      Defaults to csvOutputType.

    @param interpolationDataDir: .interpol files will be in that subfolder

    @param abaqusVersion: Version string of Abaqus, something like "6.13-2".
       Note: separators must be exactly as in the example: dot and dash.
       For new style use something like "abq2018", ignoring the hot-fix
       postfix e.g. "hf4".
       If not specified defaults to what the 'abaqus' shell command
       invokes.

    @Note: Possibly add later: option to group values: all, average, min, max
       for points with same name
    """

    # initialize the output points
    msg("Reading closure stations from %s" % gridDataPath)
    tab = RowIterFromCsvFile(
        open(gridDataPath, "rb"),
        columns=[
            "layer",
            (["x0", "y0", "z0"], float, "pt0"),
            (["x1", "y1", "z1"], float, "pt1"),
            (["dx", "dy", "dz"], float, "direct"),
            ], )
    gridPoints = []  # 2*N items: [x,y,z] lists of floats
    normals = []     # 2*N items: unit vectors along lines
    names = []       # N items: strings
    centroids = []   # N items: [x,y,z]-coords of station centroid
    nameIds = defaultdict(list)  # {name: list of ids i in [0,N]}
    for idx, line in enumerate(tab.getContainerIter()):
        gridPoints.extend((line.pt0, line.pt1))
        direct = norm(line.direct)
        normals.extend((vector_scale(direct,-1), direct))
        names.append(line.layer)
        centroids.append(vector_scale(vector_plus(line.pt0, line.pt1), 0.5))
        nameIds[line.layer].append(idx)
    del tab, direct, idx
    # convert nameIds to a list of (name, indexlist)-tuples
    nameIds = [(n, nameIds[n]) for n in sorted(nameIds)]
    msg("Found %d stations (%d points) in %d groups."
        % (len(centroids), len(gridPoints), len(nameIds)))

    #--- get data from odb and store in FramesFieldsOnUnstructuredPoints
    msg("Opening the odb %s." % odbPath)
    odb = OdbReaderInterpolateToPoints(odbPath, version=abaqusVersion)

    # check output directory
    if outputDirName and not os.path.isdir(outputDirName):
        os.makedirs(outputDirName)
    if outputDirNameGroups and not os.path.isdir(outputDirNameGroups):
        os.makedirs(outputDirNameGroups)
    if not outputDirNameGroups:
        outputDirNameGroups = outputDirName

    # pickle data file name
    # from interpolationDataDir and others
    interpolPickleName = os.path.join(
        interpolationDataDir, "_".join((projectPrefix, meshVersion, gridName)))

    # check interpolationDataDir directory
    if interpolationDataDir and not os.path.isdir(interpolationDataDir):
        os.makedirs(interpolationDataDir)
        msg("Created folder for interpolation data %s" % interpolationDataDir)

    # default file names from outputDirName and fileNamePrefix
    if csvOutputType.startswith("OneCsv"):
        outputAllFileNameTemplate = os.path.join(
            outputDirName, fileNamePrefix+"_stations.csv")
    else:
        outputAllFileNameTemplate = os.path.join(
            outputDirName, fileNamePrefix+"_stations_%s.csv")

    if not csvOutputTypeGroups:
        csvOutputTypeGroups = csvOutputType
    if csvOutputTypeGroups.startswith("OneCsv"):
        outputCumFileNameTemplate = os.path.join(
            outputDirNameGroups, fileNamePrefix+"_groups.csv")
    else:
        outputCumFileNameTemplate = os.path.join(
            outputDirNameGroups, fileNamePrefix+"_groups_%s.csv")

    # initialize the odb reader's output grid
    odb.initOutputPoints(
        grid=MeshUnstructuredPoints(gridPoints),
        interpolPickleName=interpolPickleName)

    # initialize flags about what we need
    withClosure = any(x.startswith("closure") for x in fieldNames)
    withStress = any(x.startswith("Snormal") or x.startswith("Sshear")
                     for x in fieldNames)

    # convert stepFrameToName function:
    # The frameIdToName argument to vtk_02.FramesFieldsOnUnstructuredPoints.
    # ....initOutput --if its a function-- takes a single frameId argument:
    # a (stepNb, frameNb) - tuple.
    # On the other hand all bae.utils.odbToPointsData - methods
    # expect for their stepFrameToName argument a function taking two arguments:
    # stepNb, frameNb
    # and also define convStepFrameToName that is always a function
    if callable(stepFrameToName):
        # strings don't need conversion, but functions do
        stepFrameToNameTwoArgs = stepFrameToName
        def stepFrameToName(frameId):
            return stepFrameToNameTwoArgs(*frameId)

        # define converter function convStepFrameToName
        convStepFrameToName = stepFrameToName

    # if argument stepFrameToName=="fromPostData": define special function
    elif stepFrameToName=="fromPostData":
        # special case, frameNames from odb.postData
        if "frameNames" not in odb.postData:
            raise ValueError(
                "frameNames not defined in associated _postData.py file for"
                " odbPath=%s." % odbPath)
        def stepFrameToName(frameId):
            frameName = odb.postData["frameNames"][frameId[0]][frameId[1]]
            msg("frameName from postData: step %d, frame %d, frameName %s"
                % (frameId[0], frameId[1], frameName), debugLevel=10)
            return frameName

        # define converter function convStepFrameToName
        convStepFrameToName = stepFrameToName

    else:
        # define converter function convStepFrameToName from string
        def convStepFrameToName(x):
            return stepFrameToName % tuple(x)

    # prepare nameToFrame
    # for special relativeToFrame-option "frameNameFromPostData"
    if withClosure and relativeToFrame=="frameNameFromPostData":
        frameNames = odb.postData["frameNames"]
        dateToRefFrame = dict(
            (date, (stepNb, frameNb))
            for stepNb, frameData in frameNames.iteritems()
            for frameNb, date in enumerate(frameData))
        nameToFrame = dict()
        # to rename layernames consistently with frame names: layerNameOldToNew
        layerNameOldToNew = dict()
        for layerName in set(names):

            # slice off the frame name (usually a date string) from the end of
            # the layerName
            # We might find no, one or even more "dates" to fit to the end of
            # the layerName...
            foundNames = [date for date in dateToRefFrame
                          if layerName.endswith("_%s" % date)]
            if not foundNames:
                sortedDateFrames = sorted(dateToRefFrame.iteritems(),
                                          key=lambda x: x[1])
                firstFrame = sortedDateFrames[0][1]
                msg("WARNING: Found no matching frame name in"
                    " postData['frameNames'] for the name '%s' from the"
                    " 'layer'-column of the closure-stations csv file. Assuming"
                    " frame %s %s as reference frame. Following frameNames have"
                    " been found:\n%s" % (
                        layerName, sortedDateFrames[0][0], firstFrame,
                        "\n".join(x[0] for x in sortedDateFrames)))
                return firstFrame
                ##### THIS IS WRONG!!!! date = firstFrame might be correct
                #### BUT: more documentation is necessary!!!!!
                #### What's the layer-column used for, how is it interpreted???
            elif len(foundNames)==1:
                date = foundNames[0]
            else:
                date = sorted(foundNames, key=lambda x: len(x))[-1]
                msg("WARNING: Found more than one matching frame name in"
                    " postData['frameNames'] for the name '%s' from the"
                    " 'layer'-column of the closure-stations csv file."
                    " Selecting %s from the following list of possible matches:"
                    "\n%s" % (
                        layerName, date, foundNames))

            # look up right frame for elset
            refFrame = dateToRefFrame[date]
            nameToFrame[layerName] = refFrame

            layerNameOldToNew[layerName] = (
                layerName[:-len(date)] + convStepFrameToName(refFrame))

    # prepare nameToFrame for relativeToFrame-function
    elif withClosure and callable(relativeToFrame):
        nameToFrame = dict((name, relativeToFrame(name, odb.postData))
                           for name in set(names))

    # prepare refFrameField
    if withClosure and (
            callable(relativeToFrame)
            or relativeToFrame=="frameNameFromPostData"):

        # list of reference frames for each grid point
        # names has only N elements, there are 2*N points
        refFrameField = [nameToFrame[name] for name in names for i in (0,1)]
        msg("Identified reference frame for relative displacements for all %d"
            " closure measurement points." % len(refFrameField), debugLevel=5)

    # modify layer names (names of closure stations)
    # for special relativeToFrame-option "frameNameFromPostData"
    if (withClosure and relativeToFrame=="frameNameFromPostData"):
        # apply layerNameOldToNew to names
        names = [layerNameOldToNew[name] for name in names]
        # apply layerNameOldToNew nameIds: list of (name, indexlist)-tuples
        nameIds = [(layerNameOldToNew[n], ids) for n, ids in nameIds]

        msg("Converted frame part of %d different layer names."
            % len(layerNameOldToNew), debugLevel=5)

    # initialize outputAll ... data for each individual station
    outputAll = FramesFieldsOnUnstructuredPoints(points=centroids)
    outputAll.initOutput(
        outputAllFileNameTemplate,
        outputType=csvOutputType, filterEmpty=filterEmpty,
        frameIdToName=stepFrameToName)

    # store point name field
    outputAll.updateFields(
        createFieldObject("ptName", "point", "str", initArgs=[names,]))

    # initialize outputCum ... cumulated output like average, min, max
    centroidsCum = [centre(*(centroids[i] for i in ids)) for n, ids in nameIds]
    outputCum = FramesFieldsOnUnstructuredPoints(points=centroidsCum)
    outputCum.initOutput(
        outputCumFileNameTemplate,
        outputType=csvOutputTypeGroups, filterEmpty=filterEmpty,
        frameIdToName=stepFrameToName)

    # store point name field
    outputCum.updateFields(
        createFieldObject("ptName", "point", "str", initArgs=[
            (n for n, ids in nameIds),]))

    # store point coordinates if requested (in the requested order)
    rex = re.compile(r"^point_(start|end)$")
    for fn in fieldNames:
        res = rex.match(fn)
        if not res:
            continue
        type_ = res.group(1)
        if type_=="start":
            outputAll.updateFields(createFieldObject(
                fn, "point", "vector", initArgs=[gridPoints[0::2],]))
        elif type_=="end":
            outputAll.updateFields(createFieldObject(
                fn, "point", "vector", initArgs=[gridPoints[1::2],]))

    # store point names and (if applicable) point_XXX fields as constant data
    outputAll.storeConstantData()
    outputCum.storeConstantData()

    # initialize reference displacements
    frameToFieldU = dict()
    if withClosure and (
            callable(relativeToFrame)
            or relativeToFrame=="frameNameFromPostData"):
        # relative displacements relative to different frames per
        # closure-station, according to layerName
        # creates refFrameField, refFieldU
        allRefFrames = sorted(set(nameToFrame.itervalues()))
        msg("Reading displacements for %d reference frames."
            % len(allRefFrames))
        for stepNb, frameNb in allRefFrames:
            msg("Reading field %s from the odb, step %d, frame %d."
                % (displacementFieldName, stepNb, frameNb))
            frameToFieldU[(stepNb, frameNb)] = odb.getFieldFromOdb(
                stepNb, frameNb, displacementFieldName)
        del allRefFrames

        # reference displacements potentially from different frames
        refFieldU = [
            frameToFieldU[refFr][i]
            for i, refFr in enumerate(refFrameField)]

    elif withClosure and relativeToFrame:
        # relative displacements relative to the same frame for all points
        # creates refFieldU
        stepNb, frameNb = relativeToFrame
        refFieldU = odb.getFieldFromOdb(
            stepNb, frameNb, displacementFieldName)
        frameToFieldU[relativeToFrame] = refFieldU
    else:
        # no relative displacements
        # creates dummy refFieldU
        refFieldU = None

    # get the field data from the odb frame by frame
    msg("Now starting to iterate over output frames.")
    for stepNb, frameNb in frameList:

        # precompute closure
        if withClosure:
            try:
                fieldU = frameToFieldU[(stepNb, frameNb)]
            except KeyError:
                msg("Reading field %s from the odb, step %d, frame %d."
                    % (displacementFieldName, stepNb, frameNb))
                fieldU = odb.getFieldFromOdb(
                    stepNb, frameNb, displacementFieldName)
            else:
                msg("Skip reading field %s from the odb, step %d, frame %d."
                    " Already got it as reference data."
                    % (displacementFieldName, stepNb, frameNb))

            # calculate relative displacements
            if relativeToFrame:
                if callable(relativeToFrame):
                    fieldU = type(fieldU)(
                        (((stepNb, frameNb)>refFr) and vector(u0,u1)
                         or [0.0, 0.0, 0.0])
                        for u0, u1, refFr in izip(
                                refFieldU, fieldU, refFrameField))
                elif (stepNb, frameNb) > relativeToFrame:
                    fieldU = type(fieldU)(
                        vector(u0,u1)
                        for u0, u1 in izip(refFieldU, fieldU))
                else:
                    fieldU = type(fieldU)(([0.0, 0.0, 0.0],)*len(fieldU))

            # compute closure
            iterU = iter(fieldU)
            iterN = iter(normals)
            fieldCl = createFieldObject(
                "closure", "point", "scalar", initArgs=[
                    (dot(vector(u1,u2), n1)
                     for u1, u2, n1, n2 in izip(iterU, iterU, iterN, iterN)),])

        # precompute stress
        if withStress:
            # calculate all relevant values (no matter if actually needed)
            msg("Reading field %s from the odb, step %d, frame %d."
                % (stressFieldName, stepNb, frameNb))
            fieldS = odb.getFieldFromOdb(
                stepNb, frameNb, stressFieldName)
            Snormal = []
            Sshear = []
            for S, normal in izip(fieldS, normals):
                # S components are S_11, S_22, S_33, S_12, S_13, S_23
                s = [[S[0],S[3],S[4]],
                     [S[3],S[1],S[5]],
                     [S[4],S[5],S[2]]]
                svec = mat_multvec(s, normal)
                sn = -dot(svec, normal)
                Snormal.append(sn)
                Sshear.append(length(vector(svec, vector_scale(normal, -sn))))

        for fieldName in fieldNames:

            # closure
            if fieldName=="closure":
                outputAll.updateFields(fieldCl)
                outputCum.updateFields(createFieldObject(
                    fieldName, "point", "scalar", initArgs=[
                        ((sum(fieldCl[i] for i in ids)/len(ids))
                         for n, ids in nameIds),]))
            elif fieldName=="closure_min":
                outputCum.updateFields(createFieldObject(
                    fieldName, "point", "scalar", initArgs=[
                        (min(fieldCl[i] for i in ids)
                         for n, ids in nameIds),]))
            elif fieldName=="closure_max":
                outputCum.updateFields(createFieldObject(
                    fieldName, "point", "scalar", initArgs=[
                        (max(fieldCl[i] for i in ids)
                         for n, ids in nameIds),]))

            # stress
            iterS = None  # initialize as flag. None means: not a stress field
            if fieldName.startswith("Snormal"):
                iterS = iter(Snormal)
            elif fieldName.startswith("Sshear"):
                iterS = iter(Sshear)
            if iterS is not None:
                try:
                    fn1, fn2 = fieldName.split("_", 1)
                except ValueError:
                    # no "_" in fieldName, i.e. fieldName=="Snormal", "Sshear"
                    # average stress requested
                    fldS = createFieldObject(
                        fieldName, "point", "scalar", initArgs=[
                            (0.5*(s0+s1) for s0,s1 in izip(*(iterS, iterS))),])
                    outputAll.updateFields(fldS)
                    outputCum.updateFields(createFieldObject(
                        fieldName, "point", "scalar", initArgs=[
                            ((sum(fldS[i] for i in ids)/len(ids))
                             for n, ids in nameIds),]))
                else:
                    # found "_" in fieldName, fn2 is suffix "start" or "end"
                    if fn2=="start":
                        fldS = createFieldObject(
                            fieldName, "point", "scalar", initArgs=[
                                (s0 for s0,s1 in izip(*(iterS, iterS))),])
                    elif fn2=="end":
                        fldS = createFieldObject(
                            fieldName, "point", "scalar", initArgs=[
                                (s1 for s0,s1 in izip(*(iterS, iterS))),])
                    outputAll.updateFields(fldS)

        # finally store all the data
        outputAll.storeFrame((stepNb, frameNb))
        outputCum.storeFrame((stepNb, frameNb))

    # finalize output
    outputAll.closeCsvOutput()
    outputCum.closeCsvOutput()


#------------------------------------------------------------------------------
# converter classes for odbToVtkGridBox and others
#------------------------------------------------------------------------------

class Converter(object):
    @staticmethod
    def mergeGetExtraFields(converters):
        """Wraps getExtraFields-funtions/generators of multiple converters.
        This can be usefull when you want to use more than one converter during 
        the call of odbToVTKGridBox (and others).

        Example:
            >>> converterFOS = Converter_StandardAndFOS( ... )
            >>> converterCombinedVTKs = Converter_CombinedVtk_FS4( ... )
            >>> # do all your preparation here
            >>> mergedGetExtraFields = Converter.mergeGetExtraFields(
            ...                         [converterFOS, converterCombinedVTKs])
            >>> odbToVtkGridBox(...,
            >>>                 getExtraFields=mergedGetExtraFields,
            >>>                 ...)

        @param converters: iterable of Converter-instances 
        """
        def getExtraFields(*args, **kwargs):
            for converter in converters:
                if not hasattr(converter, 'getExtraFields'):
                    continue
                for res in converter.getExtraFields(*args, **kwargs):
                    yield res

        return getExtraFields


class Converter_PlasticStrainTensor(Converter):
    """
    Usage:
     >>> import os
     >>> from bae.abq_model_02 import Model
     >>> from bae.utils.odbToPointsData_02 import \\
     >>>     odbToVtkGridBox, \\
     >>>     Converter_PlasticStrainTensor as Converter
     >>>
     >>> #-----
     >>> boxName = "box2_mine_h06"
     >>> gridData = ...
     >>> ...
     >>>
     >>> #-----
     >>>
     >>> matNameToElProps = Converter.getMatNameToElPropsFromModel(
     >>>     Model().read("MyProject_material_M01.inp"),
     >>>     iPropsG=10, iPropsK=11)
     >>>
     >>> stateSDV = "SDV8"
     >>> def maskFieldFunc(data, mesh):
     >>>     return [(state<1.5) for state in data[stateSDV]]
     >>>
     >>> converter = Converter(
     >>>     matNameToElProps=matNameToElProps,
     >>>     stateSDV=stateSDV, maskFieldFunc=maskFieldFunc)
     >>>
     >>> odbToVtkGridBox(
     >>>     odbPath, frameList=frameList,
     >>>     gridName=boxName, gridData=gridData,
     >>>     outputDirName=outputDirName, fileNamePrefix=fileNamePrefix,
     >>>     projectPrefix=..., meshVersion=...,
     >>>     **converter.fieldArguments)
     >>>
     >>> # remove the matList.csv file
     >>> os.remove(os.path.join(outputDirName, fileNamePrefix+"_matList.csv"))
    """

    def __init__(
            self,
            matNameToElProps,
            fieldName="EP",
            stateSDV="SDV8",
            maskFieldFunc=(
                lambda data, mesh: [(state<1.5) for state in data["SDV8"]]),
            valueNotDefined=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]):
        """
        @param matNameToElProps: a dict {matName : (G, K)} giving shear and
           bulk modulus. See L{getMatNameToElPropsFromModel}.
        @param fieldName: Name of the resulting field
        @param stateSDV: Name of the state field in the odb.
        @param fieldName: Name of the result field
        @param maskFieldFunc: A function generating a mask field: True where
           the plastic strain shall be calculated and False where we want the
           constant valueNotDefined as output
        @param valueNotDefined: Take this value for points that are masked out
           by maskFieldFunc and for which we don't have results in the odb
           (for example above the surface).
        """
        self.fieldName=fieldName
        self.stateSDV = stateSDV
        self.maskFieldFunc = maskFieldFunc
        self.valueNotDefined = valueNotDefined
        self.fieldPos = None

        # dummy class with arbitrary position
        # will be corrected from first call to self.compose1()
        self.FieldClass = createFieldClass(
            self.fieldName, "structPt", "tensor")

        # prepare material props
        self.matNamesList = sorted(matNameToElProps)
        self.matNbToElProps = dict(
            (i+1, matNameToElProps[matName])
            for i, matName in enumerate(self.matNamesList))

    @staticmethod
    def getMatNameToElPropsFromModel(model, iPropsG=10, iPropsK=11):
        """Return the matNameToElProps dictionary required for the
        Converter_PlasticStrainTensor.__init__ function.

        @param model: a L{bae.abq_model_02.Model} object containing material
            data. E.g.
             >>> model = Model().read("Cadia2016033_material_M21.inp")

        @param iPropsG: zero based index of the shear modulus in the data lines
            of the Abaqus *USER MATERIAL key card
        @param iPropsK: zero based index of the bulk modulus
        """

        matNameToElProps = dict()  # {matName: (G, K)}
        for matName, mat in model.material.iteritems():
            try:
                G = mat['Umat'][iPropsG]
                K = mat['Umat'][iPropsK]
            except KeyError:
                raise KeyError(
                    "Material %s is no user material, other types not"
                    " implemented yet. (Note: this might be very easy, look at"
                    " the comments in odbToPointsData_02."
                    "Converter_PlasticStrainTensor.)" % matName)
                # For materials using the Abaqus key card *ELASTIC query
                # mat["EMod"] and mat["Nue"]. Then use these formulas:
                # G = 0.5 * EMod / (1+Nue)
                # K = EMod / (3 - 6*Nue)
            except IndexError:
                raise IndexError(
                    "Wrong values for parameters iPropsG or iPropsK. Material"
                    " %s has only %d properties." % len(mat['Umat']))
            # store G and K in matNbToElProps
            matNameToElProps[matName] = (G, K)
        msg("Got elastic properties for %d materials. G in [%g, %g], K in"
            " [%g, %g]" % (len(matNameToElProps),
                           min(x[0] for x in matNameToElProps.itervalues()),
                           max(x[0] for x in matNameToElProps.itervalues()),
                           min(x[1] for x in matNameToElProps.itervalues()),
                           max(x[1] for x in matNameToElProps.itervalues()),
                           ))
        return matNameToElProps

    @staticmethod
    def calc_E_plast(E, S, G, K):
        """return e_plast = E - E_elast,
        E_elast = compliance matrix * Stress"""
        a = (3.0*K+G)/(9.0*K*G)  # = 1/3G + 1/9K
        b = (G-1.5*K)/(9.0*K*G)  # = 1/9K - 1/6G
        c = 0.5/G
        return [
            E[0] - a*S[0] + b*S[1] + b*S[2],
            E[1] - b*S[0] + a*S[1] + b*S[2],
            E[2] - b*S[0] + b*S[1] + a*S[2],
            E[3] - c*S[3],
            E[4] - c*S[4],
            E[5] - c*S[5]]

    def composer1(self, cls, dat, mesh):
        "computes plastic strain tensor"

        # get some structural info from result class on very first invocation
        # and prepare class for intermediate storage and componentsDict:
        if self.fieldPos is None:
            self.fieldPos = cls.position
            self.FieldClass = createFieldClass(
                self.fieldName, self.fieldPos, "tensor")
            # self.componentsDict = {"EP_11":0, "EP_22":1, ...}
            self.componentsDict = dict(
                (x,i) for i, x in enumerate(
                    self.FieldClass.getComponentNamesList()))

        # create maskfield and matPropsField
        maskField = self.maskFieldFunc(dat, mesh)
        matPropsField = [self.matNbToElProps.get(matNb, None)
                         for matNb in dat["Material"]]

        # initialize store
        self.store = self.FieldClass(
            (mask and matGK and self.calc_E_plast(E, S, matGK[0], matGK[1])
             or self.valueNotDefined)
            for E, S, mask, matGK in izip(
                dat["E"], dat["S"], maskField, matPropsField))

        # return correct component
        i = self.componentsDict[cls.fieldName]
        return cls(vec[i] for vec in self.store)

    def composer2(self, cls, dat, mesh):
        "return previously computed strain tensor"
        # return correct component
        i = self.componentsDict[cls.fieldName]
        return cls(vec[i] for vec in self.store)

    @property
    def fieldArguments(self):
        # a list of ("EP_11", "scalar", self.composerX) - tuples
        outputFields=[
            (x, "scalar", ((i==0) and self.composer1 or self.composer2))
            for i, x in enumerate(self.FieldClass.getComponentNamesList())]

        return dict(
            fieldNames=["E", "S", self.stateSDV],
            relativeFieldNames=["E", "S"],
            matRegionsOption=True,  # creates a field Material
            matNamesList=self.matNamesList,
            outputFields=outputFields,
            )


class Converter_StressEigenDecomp(Converter):
    """
    Usage:
     >>> import os
     >>> from bae.abq_model_02 import Model
     >>> from bae.utils.odbToPointsData_02 import \\
     >>>     odbToVtkGridBox, \\
     >>>     Converter_StressEigenDecomp as Converter
     >>>
     >>> #-----
     >>> boxName = "box2_mine_h06"
     >>> gridData = ...
     >>> ...
     >>>
     >>> #-----
     >>>
     >>> converter = Converter(
     >>>     outputFields=[
     >>>         "S1", "Bearing1", "Plunge1", "S3", "Bearing3", "Plunge3"])
     >>>
     >>> odbToVtkGridBox(
     >>>     odbPath, frameList=frameList,
     >>>     gridName=boxName, gridData=gridData,
     >>>     outputDirName=outputDirName, fileNamePrefix=fileNamePrefix,
     >>>     projectPrefix=..., meshVersion=...,
     >>>     **converter.fieldArguments)
     >>>
     >>> # remove the matList.csv file
     >>> os.remove(os.path.join(outputDirName, fileNamePrefix+"_matList.csv"))


    Note:

    This requires numpy >= 1.8.
    See /mnt/boab2/it/install/PYTHON/python-numpy-upgrade4_SLES11SP3.txt
    how to install the required version on nuku SLES11-SP3 servers.
    """

    def __init__(self, outputFields, stressFieldName="S",
                 miningConvention=True):
        """
        @param stressFieldName: optionally override input data odb-field name

        @param outputFields: A list of the desired output components. Which can
           be any combination of the eigenvalues "S1", "S2", "S3", normalized
           eigenvector components "v1x", "v1y", "v1z", "v2x", "v2y", ..., "v3z",
           "Bearing1", "Plunge1", "Bearing2", "Plunge2", "Bearing3", "Plunge3",
           as well as any of the cartesian components "S_11", "S_22", "S_33",
           "S_12", "S13", "S_23". If stressFieldName is given different from
           the default "S" then the real component names of the odb field must
           be used instead of "S_11", ...

        @param miningConvention: if True, the sign and order of eigenValues
           will be inverted
        """
        self.stressFieldName = stressFieldName
        self.outputFields = outputFields
        self.miningConvention = miningConvention

    def composer1(self, cls, dat, mesh):
        "computes eigenvalues and vectors"

        import numpy as np
        if map(int, np.__version__.split(".")) <= [1,8]:
            raise NotImplementedError(
                "This composer requires numpy v1.8 or greater. Found %s."
                % np.__version__)

        self.store = dict()

        msg("Storing cartesian components.")
        stressField = dat[self.stressFieldName]
        Svec = np.array(stressField)
        compStrToIndex = dict()
        for i, name in enumerate(stressField.getComponentNamesList()):
            self.store[name] = Svec[:,i]
            compStrToIndex[name.rsplit("_",1)[-1]] = i

        msg("Preparing stress tensors.")
        N = Svec.shape[0]
        Smat = np.zeros((N, 3,3), float)
        for j, k, compStr in [
                (0,0, "11"), (1,1, "22"), (2,2, "33"),
                (0,1, "12"), (1,0, "12"), (0,2, "13"),
                (2,0, "13"), (1,2, "23"), (2,1, "23")]:
            iVec = compStrToIndex[compStr]
            Smat[:, j, k] = Svec[:, iVec]

        msg("Computing eigenvalues and -vectors.")
        # With using eigh (for symetric matrices) the EigenVectors v will
        # be normalized and the order will be ascending
        w, v = np.linalg.eigh(Smat)

        if self.miningConvention:
            w *= -1     # invert sign for eigenValues but preserve order
        else:
            w = w[:,::-1]   # preserve sign but reversed order
            v = v[:,:,::-1]

        msg("Storing eigenvalues and vectors.")
        self.store["S1"]  = w[:,0]
        self.store["v1x"] = v[:,0,0]
        self.store["v1y"] = v[:,1,0]
        self.store["v1z"] = v[:,2,0]
        self.store["S2"]  = w[:,1]
        self.store["v2x"] = v[:,0,1]
        self.store["v2y"] = v[:,1,1]
        self.store["v2z"] = v[:,2,1]
        self.store["S3"]  = w[:,2]
        self.store["v3x"] = v[:,0,2]
        self.store["v3y"] = v[:,1,2]
        self.store["v3z"] = v[:,2,2]

        # calculate Plunge and Bearing if requested
        needsPlg = any(['plunge' in fld.lower() for fld in self.outputFields])
        needsBrg = any(['bearing' in fld.lower() for fld in self.outputFields])
        if needsPlg or needsBrg:
            msg("Calculating and storing bearing and plunge.")
            for i in range(3):

                # Bearing will be in the x-y plane...
                bearing = np.arctan2(v[:,1,i], v[:,0,i])

                # Plunge is the angle between the line and the x-y plane
                plunge = np.arcsin(-v[:,2,i])

                # Convert back to azimuths in degrees..
                plunge, bearing = np.degrees(plunge), np.degrees(bearing)
                bearing = 90 - bearing
                bearing[bearing < 0] += 360

                # If the plunge angle is upwards, get the opposite end of the
                # line
                upwards = plunge < 0
                plunge[upwards] *= -1
                bearing[upwards] -= 180
                bearing[upwards & (bearing < 0)] += 360

                # storing results
                self.store["Bearing%d" % (i+1)] = bearing
                self.store["Plunge%d" % (i+1)] = plunge

        # return correct component, convert to field_01-cls
        return cls(self.store[cls.fieldName])

    def composer2(self, cls, dat, mesh):
        "return previously computed result component, convert to field_01-cls"

        return cls(self.store[cls.fieldName])

    @property
    def fieldArguments(self):
        """This property yields the fieldNames and outputFields arguments
        suitable for L{odbToVtkGridBox} et al. if you want only the fields
        generated by this composer.
        See the outputFields argument of the
        L{constructor<Converter_StressEigenDecomp.__init__>}.
        """
        # A list of ("S1", "scalar", self.composerX) - tuples
        # The first item --the field name-- is taken from self.outputFields.
        outputFields=[
            (x, "scalar", ((i==0) and self.composer1 or self.composer2))
            for i, x in enumerate(self.outputFields)]

        return dict(
            fieldNames=[self.stressFieldName],
            outputFields=outputFields,
            )


class Converter_StandardVtk(Converter):
    """Converter to create standard state and logP values from Abaqus output.

    For cave coupling converter see L{Converter_CombinedVtk_FS4} and
    L{Converter_CombinedVtk}.

    Usage:
     >>> import os.path
     >>> from bae.utils.odbToPointsData_02 import \\
     >>>     odbToVtkGridBox, Converter_StandardVtk
     >>>
     >>> # from outputVarNames.py
     >>> varNameDamage = "SDV1"
     >>> varNameState = "SDV8"
     >>>
     >>> converter = Converter_StandardVtk(
     >>>     varNameU="U_%d",
     >>>     varNameDamage=varNameDamage, varNameState=varNameState)
     >>>
     >>> odbToVtkGridBox(
     >>>     odbPath, frameList=frameList,
     >>>     gridName=boxName, gridData=gridData,
     >>>     outputDirName=outputDirName, fileNamePrefix=fileNamePrefix,
     >>>     projectPrefix=..., meshVersion=...,
     >>>     fieldNames=[
     >>>         converter.varNameU%1,converter.varNameU%2,converter.varNameU%3,
     >>>         "S_MIN_PRINCIPAL",
     >>>         varNameDamage,
     >>>         (varNameState, "const"),],
     >>>     outputFields=[
     >>>         converter.getU_1, converter.getU_2, converter.getU_3,
     >>>         ("logP", "scalar", converter.getLogP),
     >>>         ("S1", "scalar", converter.getS1),
     >>>         ("state", "scalar", converter.getState),
     >>>     ])
    """

    # The functions L{odbToVtkGridBox} and the like accept a parameter
    # getExtraFields. Thereby a method to read extra data e.g. from Cavesim
    # can be supplied.
    # This converter doesn't need extra fields. Passing this None-valued
    # attribute of the converter to odbToVtkGridBox makes it skip this
    # functionality.
    # This attribute is provided for convenience. So the call to
    # odbToVtkGridBox can have the same set of arguments whether
    # Converter_StandardVtk or Converter_CombinedVtk is used.
    getExtraFields = None

    def __init__(self, varNameDamage="SDV1", varNameState="SDV8",
                 varNameU="U_%d"):
        """
        @param varNameDamage: variable name for plastic strain, e.g. "SDV1"
        @param varNameState: variable name for state, e.g. "SDV8"
        @param varNameU: "U_%d" or "UR_%d" or something like that
        """
        self.varNameState = varNameState
        self.varNameDamage = varNameDamage

        self.varNameU = varNameU
        self.getU_1 = varNameU % 1
        self.getU_2 = varNameU % 2
        self.getU_3 = varNameU % 3

    def getLogP(self, cls, dat, mesh):
        "provides the logP field"
        return cls((
            # 2.0 if cave (state==5) ... even if damage<=0
            2.0 if 4.5<state<5.5 else

            # 0.0 if degrading (state==2) or void (state==3) or deleted ...
            #        ... (state==0) or air (state==6) or negative damage
            0.0 if 1.5<state<3.5 or state<=0.5 or state>=5.5 or damage<=0.0 else

            # logP-formula if rock (state==1) or fill (state==4)
            log10(damage*1000+1))

            # provide the series of values...
            for damage, state in izip(
                    *map(dat.get, [self.varNameDamage, self.varNameState])))

    def getS1(self, cls, dat, mesh):
        return cls(-x for x in dat["S_MIN_PRINCIPAL"])

    def getS3(self, cls, dat, mesh):
        return cls(-x for x in dat["S_MAX_PRINCIPAL"])

    def getState(self, cls, dat, mesh):
        """Convert state into:
         0 : air (deleted elements, outside the mesh)
         1 : void (Abq-state 2, 3, 6)
         2 : rock (Abq-state 1)
         3 : fill (Abq-state 4)
         4 : CAVE (Abq-state 5)
        """
        conv = {
            # note: results must be floats, otherwise Voxler creates two
            # input objects
            0: 0.0,
            1: 2.0,
            2: 1.0,
            3: 1.0,
            4: 3.0,
            5: 4.0,
            6: 1.0}
        return cls(conv.get(int(round(x)),-1.0) for x in dat[self.varNameState])


class Converter_CombinedVtk(Converter_StandardVtk):
    """Converter to create "combined Abaqus/Cavesim vtks".

    For Abaqus-only result converter see L{Converter_StandardVtk}.

    Usage:
     >>> import os.path
     >>> from bae.utils.odbToPointsData_02 import \\
     >>>     odbToVtkGridBox, Converter_CombinedVtk
     >>>
     >>> # from outputVarNames.py
     >>> varNameDamage = "SDV1"
     >>> varNameState = "SDV8"
     >>>
     >>> converter = Converter_CombinedVtk(
     >>>     csOutputDir="../../2_RUN/R02b_RST02/CavesimOutput",
     >>>     stepFrameToCsStep={1: [-1,]*8 + range(1,73)},
     >>>     varNameU="U_%d",
     >>>     varNameDamage=varNameDamage, varNameState=varNameState)
     >>>
     >>> odbToVtkGridBox(
     >>>     odbPath, frameList=frameList,
     >>>     gridName=boxName, gridData=gridData,
     >>>     outputDirName=outputDirName, fileNamePrefix=fileNamePrefix,
     >>>     projectPrefix=..., meshVersion=...,
     >>>     fieldNames=[
     >>>         converter.varNameU%1,converter.varNameU%2,converter.varNameU%3,
     >>>         "S_MIN_PRINCIPAL", varNameDamage,
     >>>         (varNameState, "const")],
     >>>     getExtraFields=converter.getExtraFields,
     >>>     outputFields=[
     >>>         ("U_1", "scalar", converter.getU_1),
     >>>         ("U_2", "scalar", converter.getU_2),
     >>>         ("U_3", "scalar", converter.getU_3),
     >>>         ("logP", "scalar", converter.getLogP),
     >>>         ("S1", "scalar", converter.getS1),
     >>>         ("state", "scalar", converter.getState),
     >>>     ])
    """

    def __init__(self, stepFrameToCsStep, csOutputDir="./CavesimOutput",
                 varNameU="U_%d", varNameDamage="SDV1", varNameState="SDV8"):
        """
        @param csOutputDir: Folder of the Cavesim result files. This folder
            must contain the folders FLOW_VTK and VELOCITY_VTK.
        @param stepFrameToCsStep: provides the Cavesim/FS4 step for a given odb
            step and frame number: stepFrameToCsStep[stepNb][frameNb] must
            yield the right Cavesim step number --an integer >0. For frames
            without Cavesim step --i.e. not coupled steps, e.g. before
            cave coupling starts-- this should yield -1.
        @param varNameU: "U_%d" or "UR_%d" or something like that
        @param varNameDamage: variable name for plastic strain, e.g. "SDV1"
        @param varNameState: variable name for state, e.g. "SDV8"
        """
        self.stepFrameToCsStep = stepFrameToCsStep
        self.csOutputDir = csOutputDir
        self.varNameU = varNameU
        self.varNameState = varNameState
        self.varNameDamage = varNameDamage

        # flag indicating if Cavesim results are there
        # Defaults to True: if getExtraFields is not being called then an error
        # occurs
        self.csData = True

    @staticmethod
    def readVTKfromFileOrZipArch(vtkPath, pauseTime=5):
        """Try to find the cavesim output. Either as plain vtk file vtkPath
        or in a zip-archive. In the latter case the vtk file is temporarily
        extracted from the zip-archive.

         - Read the vtk file from vtkPath if present.
         - Check if there is a zip file with the same name, i.e. for
           FS4out_step_064.vtk try FS4out_step_064.vtk.zip and
           FS4out_step_064.zip.
         - Assume that the folder directly holding the vtk file has been zipped,
           i.e. for vtkPath=".../CavesimOutput/VELOCITY_VTK/CAVESIM_VEL_10.vtk"
           try: in archive ".../CavesimOutput/VELOCITY_VTK.zip" members
           "CAVESIM_VEL_10.vtk" or "VELOCITY_VTK/CAVESIM_VEL_10.vtk"
         - Assume that a folder-structure has been zipped and the zip-file has
           the name of the folder plus extension ".zip", i.e.
           for vtkPath=".../OUTPUT.zip/VELOCITY_VTK/CAVESIM_VEL_10.vtk"
           try: in archive ".../OUTPUT.zip" members
           "VELOCITY_VTK/CAVESIM_VEL_10.vtk" or
           "OUTPUT/VELOCITY_VTK/CAVESIM_VEL_10.vtk"
         - Assume that a folder-structure has been zipped and the zip-file has
           the name of the folder plus a prefix like "NCA_Perse_" plus
           extension ".zip", i.e. for
           vtkPath=".../NCA_Perse_OUTPUT.zip/VELOCITY_VTK/CAVESIM_VEL_10.vtk"
           try: in archive ".../NCA_Perse_OUTPUT.zip" members
           "VELOCITY_VTK/CAVESIM_VEL_10.vtk" or
           "OUTPUT/VELOCITY_VTK/CAVESIM_VEL_10.vtk"

        @param vtkPath: vtk-file name
        @param pauseTime: seconds to wait after zip file has changed before
           trying to read again
        @returns: vtk object
        """
        failedPaths = list()

        # simply try the path as vtk file
        if os.path.isfile(vtkPath):
            vtkData = FieldsOnStructPointGrid().fromVtk(vtkPath)
            return vtkData
        else:
            failedPaths.append(vtkPath)

        # Failed to locate the file directly, now trying some heuristics.
        # This while clause will either be left with a break after defining
        # the variables archivePath and extractPathCandidates or raise an error
        while 1:

            # if there is a zip file with the given file name
            # i.e. for FS4out_step_064.vtk try FS4out_step_064.vtk.zip
            archivePath = vtkPath + ".zip"
            if os.path.isfile(archivePath):
                extractPathCandidates = [os.path.basename(vtkPath)]
                break
            else:
                failedPaths.append(archivePath)
                del archivePath

            # if there is a zip file with the given file name ex .vtk
            # i.e. for FS4out_step_064.vtk try FS4out_step_064.zip
            archivePath = os.path.splitext(vtkPath)[0] + ".zip"
            if os.path.isfile(archivePath):
                extractPathCandidates = [os.path.basename(vtkPath)]
                break
            else:
                failedPaths.append(archivePath)
                del archivePath

            # try to detect zipped folders...
            # firstly prepare list of folders in vtkPath
            fileName = os.path.basename(vtkPath)

            dirNames = []
            path = vtkPath
            while 1:
                path, folder = os.path.split(path)
                if folder != "":
                    dirNames.append(folder)
                else:
                    if path != "":
                        dirNames.append(path)
                    break
            dirNames.reverse()
            del path, folder

            # if zip file in path, i.e. ...
            # for vtkPath=".../OUTPUT.zip/VELOCITY_VTK/CAVESIM_VEL_10.vtk"
            # try archive ".../OUTPUT.zip"
            for depth in reversed(range(len(dirNames))):
                dirName = dirNames[depth]
                base, ext = os.path.splitext(dirName)
                if ext.lower()==".zip" and not dirName.startswith("."):
                    archivePath = os.path.join(*dirNames[:depth+1])
                    remainingPath = os.path.join(*dirNames[depth+1:])
                    break

            if ext.lower()==".zip":
                if os.path.isfile(archivePath):
                    extractPathCandidates = [
                        # try VELOCITY_VTK/CAVESIM_VEL_10.vtk
                        remainingPath,
                        # try OUTPUT/VELOCITY_VTK/CAVESIM_VEL_10.vtk
                        os.path.join(base, remainingPath),
                        # for archivePath=".../NCA_Perse_OUTPUT.zip"
                        # try OUTPUT/VELOCITY_VTK/CAVESIM_VEL_10.vtk
                        os.path.join(base.rsplit("_", 1)[-1], remainingPath),
                        ]
                    break
                else:
                    failedPaths.append(archivePath)

            # try folder name as zip-file path:
            # i.e. for vtkPath=".../CsOutput/VELOCITY_VTK/CAVESIM_VEL_10.vtk"
            # try archive ".../CsOutput/VELOCITY_VTK.zip"
            # or for vtkPath=".../OUTPUT/VELOCITY_VTK/CAVESIM_VEL_10.vtk"
            # try archive ".../OUTPUT/VELOCITY_VTK.zip" then ".../OUTPUT.zip"
            for depth in reversed(range(len(dirNames))):
                if dirNames[depth].startswith("."):
                    break
                archivePath = os.path.join(*dirNames[:depth+1]) + ".zip"
                if os.path.isfile(archivePath):
                    remainingPath = os.path.join(*dirNames[depth+1:])
                    extractPathCandidates = [
                        # try CAVESIM_VEL_10.vtk
                        os.path.join(remainingPath, fileName),
                        # try VELOCITY_VTK/CAVESIM_VEL_10.vtk
                        os.path.join(dirNames[depth], remainingPath, fileName),]
                    break
                else:
                    failedPaths.append(archivePath)

            # failed to locate vtk file or zip archive
            raise IOError(
                "Can't find vtkPath or corresponding zip-archive.\nTried: %s"
                % (', '.join(failedPaths)))

        # process the archive we just found
        archive = zipfile.ZipFile(archivePath, 'r')
        stored = archive.namelist()

        extractPath = None
        for extractPath in extractPathCandidates:
            if extractPath in stored:
                break
            else:
                extractPath = None
        if not extractPath:
            raise IOError(
                "Can't identify the vtk %s in zip-archive %s. Tried: %s"
                % (os.path.basename(vtkPath), archivePath,
                   ", ".join(extractPathCandidates)))

        msg("Reading %s from zip archive %s" % (extractPath, archivePath),
            debugLevel=1)

        tmpVtkfile = None
        for cnt in range(20):  # try 20 times to extract from zip file

            checker = CheckFileFinished(archivePath)
            try:
                tmpVtkfile = archive.extract(extractPath)
            except KeyError, exc:
                raise KeyError("Path %s vanished from zip archive %s.\n%s"
                               % (extractPath, archivePath, exc.args[0]))
            except zipfile.BadZipfile:
                if checker.finished():
                    raise
                msg("WARNING: The zip archive %s seems to have changed"
                    " while trying to extract %s. Trying again after"
                    " waiting for %s sec."
                    % (archivePath, extractPath, pauseTime))
                time.sleep(pauseTime)
            else:
                break

        # finished extracting
        if tmpVtkfile is None:
            raise IOError(
                "Can't find the vtk %s in zip-archive %s after %d attempts"
                % (os.path.basename(vtkPath), archivePath), cnt)

        # getting the data and cleaning up
        vtkData = FieldsOnStructPointGrid().fromVtk(tmpVtkfile)
        os.remove(tmpVtkfile)

        return vtkData

    def getExtraFields(self, stepNb, frameNb, points, **dummy):
        """Method to read extra data from Cavesim output vtk files.
        This method must be supplied to the parameter getExtraFields of
        functions L{odbToVtkGridBox} and the like.
        """
        try:
            csStep = self.stepFrameToCsStep[stepNb][frameNb]
        except (KeyError, IndexError):
            msg("ERROR: odb steps provided in stepFrameToCsStep: %s"
                "\nRequested step %d, frame %d."
                % (sorted(self.stepFrameToCsStep), stepNb, frameNb))
            raise
        if csStep<0:
            # no Cavesim output for current step/frame
            self.csData = False
            return

        # read the Cavesim vtk file
        fname = os.path.join(
            self.csOutputDir, "VELOCITY_VTK", "CAVESIM_VEL_%d.vtk" % csStep)
        msg("Reading %s" % fname)
        vtkVel = self.readVTKfromFileOrZipArch(fname)

        # Cavesim data available
        self.csData = True

        # interpolate Cavesim displacement fields
        for csFieldName, fieldName in [
                ("TOTAL_X_MAG", "UCS_1"),
                ("TOTAL_Y_MAG", "UCS_2"),
                ("TOTAL_Z_MAG", "UCS_3")]:
            msg("Interpolating field %s from Cavesim field %s"
                % (fieldName, csFieldName))

            # convert displacement components: x>9998 => x=0
            fld = vtkVel.data[csFieldName]
            fld = type(fld)(((x <= 9998) and x) or 0.0 for x in fld)

            # interpolate to output points
            fld = createFieldObject(fieldName, "point", "scalar", initArgs=[
                fld.interpolateToPoints(vtkVel.mesh, points),])

            yield fld

    def getU_1(self, cls, dat, mesh):
        "provides the combined U_1 field"
        if not self.csData:
            return cls(dat[self.varNameU % 1])

        return cls((
            # if not in coupling-box: take Abq-values
            # if 0.5 < state_Abq < 4.5: take Abq-values
            u_Abq if (u_CS is None) or (0.5<state_Abq<4.5)

            # if 4.5<abq_state (abq-cave,-air) or 0.5<abq_state (abq-deleted):
            # take FS4-values
            else u_CS if (state_Abq>=4.5 or state_Abq<=0.5)

            # else: return 0.0
            else 0.0)

            # provide the series of values...
            for u_CS, u_Abq, state_Abq in izip(
                    *map(dat.get, ["UCS_1",self.varNameU%1,self.varNameState])))

    def getU_2(self, cls, dat, mesh):
        "provides the combined U_2 field"
        if not self.csData:
            return cls(dat[self.varNameU % 2])

        return cls((
            # if not in coupling-box: take Abq-values
            # if 0.5 < state_Abq < 4.5: take Abq-values
            u_Abq if (u_CS is None) or (0.5<state_Abq<4.5)

            # if 4.5<abq_state (abq-cave,-air) or 0.5<abq_state (abq-deleted):
            # take FS4-values
            else u_CS if (state_Abq>=4.5 or state_Abq<=0.5)

            # else: return 0.0
            else 0.0)

            # provide the series of values...
            for u_CS, u_Abq, state_Abq in izip(
                    *map(dat.get, ["UCS_2",self.varNameU%2,self.varNameState])))

    def getU_3(self, cls, dat, mesh):
        "provides the combined U_3 field"
        if not self.csData:
            return cls(dat[self.varNameU % 3])

        return cls((
            # if not in coupling-box: take Abq-values
            # if 0.5 < state_Abq < 4.5: take Abq-values
            u_Abq if (u_CS is None) or (0.5<state_Abq<4.5)

            # if 4.5<abq_state (abq-cave,-air) or 0.5<abq_state (abq-deleted):
            # take FS4-values
            else u_CS if (state_Abq>=4.5 or state_Abq<=0.5)

            # else: return 0.0
            else 0.0)

            # provide the series of values...
            for u_CS, u_Abq, state_Abq in izip(
                    *map(dat.get, ["UCS_3",self.varNameU%3,self.varNameState])))

    def getState(self, cls, dat, mesh):
        """
        Creates those values::
          0 : AIR - In Abaqus: outside model or deleted element;
                    in Cavesim / FS4: no particles
          1 : void
          2 : rock
          3 : fill
          4 : CAVE
        """
        abqOnlyStateFld = Converter_StandardVtk.getState(self, cls, dat, mesh)
        if not self.csData:
            return abqOnlyStateFld

        return cls((
            # Correct the state==0 if there is (vertical) movement in the
            # Cavesim data
            4.0 if (state_Abq==0.0) and (u_CS!=0.0)

            # ... else take state from Abaqus
            else state_Abq)

            # provide the series of values...
            for u_CS, state_Abq in izip(
                    dat["UCS_3"], abqOnlyStateFld))


class Converter_CombinedVtk_FS4(Converter_CombinedVtk):
    """Converter to create "combined Abaqus/FS4 vtks".

    For Abaqus-only result converter see L{Converter_StandardVtk}, for the
    Cavesim pendant see L{Converter_CombinedVtk}.

    Usage:
     >>> import os.path
     >>> from bae.utils.odbToPointsData_02 import \\
     >>>     odbToVtkGridBox, Converter_CombinedVtk_FS4
     >>>
     >>> # from outputVarNames.py
     >>> varNameDamage = "SDV1"
     >>> varNameState = "SDV8"
     >>>
     >>> converter = Converter_CombinedVtk_FS4(
     >>>     csOutputTemplate="../../2_RUN/R02b_RST02/FS4Out/FS4out_%03d.vtk",
     >>>     stepFrameToCsStep={1: [-1,]*8 + range(1,73)},
     >>>     varNameU="U_%d",
     >>>     varNameDamage=varNameDamage, varNameState=varNameState)
     >>>
     >>> odbToVtkGridBox(
     >>>     odbPath, frameList=frameList,
     >>>     gridName=boxName, gridData=gridData,
     >>>     outputDirName=outputDirName, fileNamePrefix=fileNamePrefix,
     >>>     projectPrefix=..., meshVersion=...,
     >>>     fieldNames=[
     >>>         converter.varNameU%1,converter.varNameU%2,converter.varNameU%3,
     >>>         "S_MIN_PRINCIPAL", varNameDamage,
     >>>         (varNameState, "const")],
     >>>     getExtraFields=converter.getExtraFields,
     >>>     outputFields=[
     >>>         ("U_1", "scalar", converter.getU_1),
     >>>         ("U_2", "scalar", converter.getU_2),
     >>>         ("U_3", "scalar", converter.getU_3),
     >>>         ("logP", "scalar", converter.getLogP),
     >>>         ("S1", "scalar", converter.getS1),
     >>>         ("state", "scalar", converter.getState),
     >>>     ])

    If you need special values it's suggested to derive a subclass that
    overrides the corresponding getXXX-function. For example:
     >>> from bae.utils.odbToPointsData_02 import \\
     >>>     odbToVtkGridBox, Converter_CombinedVtk_FS4
     >>>
     >>> class MyConverter(Converter_CombinedVtk_FS4):
     >>>     def getState(self, cls, dat, mesh):
     >>>         normalState = super(MyConverter, self).getState(
     >>>             self, cls, dat, mesh)
     >>>         if not self.csData:
     >>>             return normalState
     >>>
     >>>         return cls(1.0 if ucs<0 else ns
     >>>                    for ns, ucs in izip(normalState, dat["UCS_MAG"]))
     >>>
     >>> converter = MyConverter(
     >>>     ...arguments as for Converter_CombinedVtk_FS4...)
     >>>
     >>> odbToVtkGridBox(
     >>>     odbPath, frameList=frameList, ...
     >>>     ...arguments as usual...)
    """

    def __init__(self, stepFrameToCsStep,
                 csOutputTemplate="./FlowsimOutput/FS4out_step_%03d.vtk",
                 varNameU="U_%d", varNameDamage="SDV1", varNameState="SDV8",
                 varNameUFS4="UPC_%s"):
        """
        @param csOutputTemplate: File path template including a "%03d"
            placeholder to be replaced by the Cavesim/FS4 step number. These
            vtk files must contain the fields identified by the varNameUFS4
            argument.
        @param stepFrameToCsStep: provides the Cavesim/FS4 step for a given odb
            step and frame number: stepFrameToCsStep[stepNb][frameNb] must
            yield the right Cavesim step number --an integer >0. For frames
            without Cavesim step --i.e. not coupled steps, e.g. before
            cave coupling starts-- this should yield -1.
        @param varNameU: "U_%d" or "UR_%d" or something like that
        @param varNameDamage: variable name for plastic strain, e.g. "SDV1"
        @param varNameState: variable name for state, e.g. "SDV8"
        @param varNameUFS4: variable name for the FS4 displacements to be
            combined with the Abaqus data. Needs a "%s" component placeholder
            that will be replaced by 1, 2, 3, "MAG".
        """
        self.stepFrameToCsStep = stepFrameToCsStep
        self.csOutputTemplate = csOutputTemplate
        self.varNameU = varNameU
        self.varNameState = varNameState
        self.varNameDamage = varNameDamage

        # flag indicating if Cavesim results are there
        # Defaults to True: if getExtraFields is not being called then an error
        # occurs
        self.csData = True
        self.varNameUFS4 = varNameUFS4

    def getExtraFields(self, stepNb, frameNb, points, **dummy):
        """Method to read extra data from FS4 output vtk files.
        This method must be supplied to the parameter getExtraFields of
        functions L{odbToVtkGridBox} and the like.
        """

        try:
            csStep = self.stepFrameToCsStep[stepNb][frameNb]
        except (KeyError, IndexError):
            msg("ERROR: odb steps provided in stepFrameToCsStep: %s"
                "\nRequested step %d, frame %d."
                % (sorted(self.stepFrameToCsStep), stepNb, frameNb))
            raise
        if csStep<0:
            # no Cavesim/FS4 output for current step/frame
            self.csData = False
            return

        # read the Cavesim/FS4 vtks
        fname = self.csOutputTemplate % csStep
        msg("Reading %s" % fname)
        vtkFlow = self.readVTKfromFileOrZipArch(fname)

        # Cavesim data available
        self.csData = True

        # interpolate fields and "yield"
        # density (0: air, -1: intact rock, >0 : number of particles in a cell)
        # for fieldName in ["du1", "du2", "du3", "dumag", "density"]:
        for x in [1, 2, 3, "MAG"]:
            csFieldName = self.varNameUFS4 % x
            fieldName = "UCS_%s" % x
            msg("Interpolating field %s from FS4-field %s"
                % (fieldName, csFieldName))

            # convert displacement components: x>9998 => x=0
            try:
                fld = vtkFlow.data[csFieldName]
            except KeyError:
                # ignore component if not available
                continue

            # interpolate to output points
            fld = createFieldObject(
                fieldName, "point", "scalar", initArgs=[
                    fld.interpolateToPoints(vtkFlow.mesh, points),])
            yield fld

    def getState(self, cls, dat, mesh):
        """
        Creates those values::
          0 : AIR - In Abaqus: outside model or deleted element;
                    in Cavesim / FS4: no particles
          1 : void
          2 : rock
          3 : fill
          4 : CAVE
        """
        abqOnlyStateFld = Converter_StandardVtk.getState(self, cls, dat, mesh)
        if not self.csData:
            return abqOnlyStateFld

        return cls((
            # Correct the state==0 if there is movement in the FS4 data
            4.0 if (state_Abq==0.0) and (u_CS>0.0)

            # ... else take state from Abaqus
            else state_Abq)

            # provide the series of values...
            for u_CS, state_Abq in izip(
                    dat["UCS_MAG"], abqOnlyStateFld))


class Converter_CombinedVtk_FS4_multiBox(Converter_CombinedVtk_FS4):
    """Converter to create "combined Abaqus/FS4 vtks" from multiple coupling
    boxes.

    Works exactly like L{Converter_CombinedVtk_FS4} (the base class) except for
    different arguments to the constructor.
    """
    def __init__(
            self, caveCouplingVtkNb, stepFrameToCsStep, csOutputTemplate,
            varNameDamage, varNameState, varNameU="U_%d",
            varNameUFS4="UPC_%s"):
        """
        @param caveCouplingVtkNb: number of coupling boxes
        @param csOutputTemplate: A list of file path templates including a
            "%03d" placeholder to be replaced by the Cavesim/FS4 step number.
            One for each coupling box. These vtk files must contain the fields
            identified by the varNameUFS4 argument.
        @param stepFrameToCsStep: A list of relations: Abaqus odb-frame to FS4
            mining step number (==Abaqus Transfer Number --ATN). One item
            for each coupling box. Each provides the Cavesim/FS4 step for a
            given odb step and frame number:
            stepFrameToCsStep[iBox][stepNb][frameNb] must yield the right
            Cavesim step number --an integer >0 for the coupling box with
            index iBox. For frames without Cavesim step --i.e. not coupled
            steps, e.g. before cave coupling starts-- this should yield -1.
        @param varNameU: "U_%d" or "UR_%d" or something like that
        @param varNameDamage: variable name for plastic strain, e.g. "SDV1"
        @param varNameState: variable name for state, e.g. "SDV8"
        @param varNameUFS4: variable name for the FS4 displacements to be
            combined with the Abaqus data. Needs a "%s" component placeholder
            that will be replaced by 1, 2, 3, "MAG".
        """

        # check special arguments
        if (not isinstance(stepFrameToCsStep, (list, tuple))
            or len(stepFrameToCsStep) != caveCouplingVtkNb):
            raise ValueError(
                "ERROR: Argument stepFrameToCsStep must be a list with %d"
                " dictionaries, corresponding to the given caveCouplingVtkNb."
                % (caveCouplingVtkNb, caveCouplingVtkNb))
        if (not isinstance(csOutputTemplate, (list, tuple))
            or len(csOutputTemplate) != caveCouplingVtkNb):
            raise ValueError(
                "ERROR: Argument csOutputTemplate must be a list with %d"
                " dictionaries, corresponding to the given caveCouplingVtkNb."
                % (caveCouplingVtkNb, caveCouplingVtkNb))

        self.caveCouplingVtkNb = caveCouplingVtkNb
        self.stepFrameToCsStep = stepFrameToCsStep
        self.csOutputTemplate = csOutputTemplate
        self.varNameU = varNameU
        self.varNameState = varNameState
        self.varNameDamage = varNameDamage
        self.varNameUFS4 = varNameUFS4

        # flag indicating if Cavesim results are there
        # Defaults to True: if getExtraFields is not being called then an error
        # occurs
        self.csData = True

    def getExtraFields(self, stepNb, frameNb, points, **dummy):
        """Method to read extra data from FS4 output vtk files.
        This method must be supplied to the parameter getExtraFields of
        functions L{odbToVtkGridBox} and the like.
        """

        vtkDataList = list()  # list of vtk data potentially for multiple boxes
        self.csData = False
        fileNames = list()  # for diagnostic output
        for iBox in range(self.caveCouplingVtkNb):

            try:
                csStep = self.stepFrameToCsStep[iBox][stepNb][frameNb]
            except (KeyError, IndexError):
                msg("ERROR: odb steps provided in stepFrameToCsStep[%d]: %s"
                    "\nRequested step %d, frame %d."
                    % (iBox, sorted(self.stepFrameToCsStep), stepNb, frameNb))
                raise

            if csStep < 0:
                # Cavesim/FS4 output not available for current step/frame
                # for this box
                continue

            self.csData = True

            # read the Cavesim/FS4 vtks
            fname = self.csOutputTemplate[iBox] % csStep
            msg("Reading %s" % fname)
            vtkDataList.append(self.readVTKfromFileOrZipArch(fname))
            fileNames.append(fname)

        # stop here if no coupling box has data for this mining step
        if not self.csData:
            return

        # interpolate fields and "yield"
        # density (0: air, -1: intact rock, >0 : number of particles in a cell)
        # for fieldName in ["du1", "du2", "du3", "dumag", "density"]:
        for x in [1, 2, 3, "MAG"]:
            csFieldName = self.varNameUFS4 % x
            fieldName = "UCS_%s" % x
            msg("Interpolating field %s from FS4-field %s"
                % (fieldName, csFieldName))

            # iterate over coupling boxes
            for iBox, vtkData in enumerate(vtkDataList):

                if csFieldName not in vtkData.data:
                    msg("WARNING: field <%s> not available in %s. (Maybe we"
                        " don't need it?)" % (csFieldName, fileNames[0]))
                    # ignore component if not available
                    continue

                newFld = vtkData.data[csFieldName]

                # interpolate to output points
                if iBox==0:
                    fld = newFld.interpolateToPoints(vtkData.mesh, points)
                else:
                    # combine: overwrite if possible
                    fld = ((x if y is None else y) for x, y in izip(
                        fld, newFld.interpolateToPoints(vtkData.mesh, points)))

            # finallize: store in correct field object
            fld = createFieldObject(
                fieldName, "point", "scalar", initArgs=[fld,])
            yield fld


class Converter_CombinedVtk_Old(Converter_StandardVtk):
    """Converter to create "combined Abaqus/Cavesim vtks".

    DEPRECATED version. Use L{Converter_CombinedVtk} instead or
    L{Converter_CombinedVtk_FS4}.

    This old version here uses both the VELOCITY_VTK and FLOW_VTK Cavesim
    output files. It was designed to cope with the artifacts due to linear
    interpolation of the Abaqus-state-SDV. Nowadays we use piecewise constant
    interpolation.

    For Abaqus-only result converter see L{Converter_StandardVtk}.

    Usage:
     >>> import os.path
     >>> from bae.utils.odbToPointsData_02 import \\
     >>>     odbToVtkGridBox, Converter_CombinedVtk
     >>>
     >>> # from outputVarNames.py
     >>> varNameDamage = "SDV1"
     >>> varNameState = "SDV8"
     >>>
     >>> converter = Converter_CombinedVtk(
     >>>     csOutputDir="../../2_RUN/R02b_RST02/CavesimOutput",
     >>>     stepFrameToCsStep={1: [-1,]*8 + range(1,73)},
     >>>     varNameU="U_%d",
     >>>     varNameDamage=varNameDamage, varNameState=varNameState)
     >>>
     >>> odbToVtkGridBox(
     >>>     odbPath, frameList=frameList,
     >>>     gridName=boxName, gridData=gridData,
     >>>     outputDirName=outputDirName, fileNamePrefix=fileNamePrefix,
     >>>     projectPrefix=..., meshVersion=...,
     >>>     fieldNames=[
     >>>         converter.varNameU%1,converter.varNameU%2,converter.varNameU%3,
     >>>         "S_MIN_PRINCIPAL", varNameDamage,
     >>>         (varNameState, "const")],
     >>>     getExtraFields=converter.getExtraFields,
     >>>     outputFields=[
     >>>         ("U_1", "scalar", converter.getU_1),
     >>>         ("U_2", "scalar", converter.getU_2),
     >>>         ("U_3", "scalar", converter.getU_3),
     >>>         ("logP", "scalar", converter.getLogP),
     >>>         ("S1", "scalar", converter.getS1),
     >>>         ("state", "scalar", converter.getState),
     >>>     ])
    """

    def __init__(self, stepFrameToCsStep, csOutputDir="./CavesimOutput",
                 varNameU="U_%d", varNameDamage="SDV1", varNameState="SDV8"):
        """
        @param csOutputDir: Folder of the Cavesim result files. This folder
            must contain the folders FLOW_VTK and VELOCITY_VTK.
        @param stepFrameToCsStep: provides the Cavesim/FS4 step for a given odb
            step and frame number: stepFrameToCsStep[stepNb][frameNb] must
            yield the right Cavesim step number --an integer >0. For frames
            without Cavesim step --i.e. not coupled steps, e.g. before
            cave coupling starts-- this should yield -1.
        @param varNameU: "U_%d" or "UR_%d" or something like that
        @param varNameDamage: variable name for plastic strain, e.g. "SDV1"
        @param varNameState: variable name for state, e.g. "SDV8"
        """
        self.stepFrameToCsStep = stepFrameToCsStep
        self.csOutputDir = csOutputDir
        self.csVelThreshold = None
        self.varNameU = varNameU
        self.varNameState = varNameState
        self.varNameDamage = varNameDamage
        self.csData = False  # flag indicating if Cavesim results are there

    @staticmethod
    def _readVTKfromFileOrZipArch(vtkPath, pauseTime=5):
        """
        Reads the vtk file from path if present. If not: try to unzip the data
        from zip-archive temporarily.

        @param vtkPath: vtk-file name
        @param pauseTime: seconds to wait after zip file has changed before
           trying to read again
        @returns: vtk object
        """
        if os.path.isfile(vtkPath):
            vtkData = FieldsOnStructPointGrid().fromVtk(vtkPath)
        elif os.path.isfile(os.path.dirname(vtkPath) + '.zip'):
            archivePath = os.path.dirname(vtkPath) + '.zip'
            archive = zipfile.ZipFile(archivePath, 'r')
            stored = archive.namelist()

            extractPath = None
            for extractPath in [
                    os.path.basename(vtkPath),
                    os.path.join(*vtkPath.split('/')[-2:])]:
                if extractPath in stored:
                    break
                else:
                    extractPath = None
            if not extractPath:
                raise IOError(
                    "Can't identify the vtk %s in zip-archive %s"
                    % (os.path.basename(vtkPath), archivePath))

            msg("Reading %s from zip archive %s" % (extractPath, archivePath),
                debugLevel=1)

            tmpVtkfile = None
            for cnt in range(20):  # try 20 times to extract from zip file

                checker = CheckFileFinished(archivePath)
                try:
                    tmpVtkfile = archive.extract(extractPath)
                except KeyError, exc:
                    raise KeyError("Path %s vanished from zip archive %s.\n%s"
                                   % (extractPath, archivePath, exc.args[0]))
                except zipfile.BadZipfile:
                    if checker.finished():
                        raise
                    msg("WARNING: The zip archive %s seems to have changed"
                        " while trying to extract %s. Trying again after"
                        " waiting for %s sec."
                        % (archivePath, extractPath, pauseTime))
                    time.sleep(pauseTime)
                else:
                    break

            # finished extracting
            if tmpVtkfile is None:
                raise IOError(
                    "Can't find the vtk %s in zip-archive %s after %d attempts"
                    % (os.path.basename(vtkPath), archivePath), cnt)

            # getting the data and cleaning up
            vtkData = FieldsOnStructPointGrid().fromVtk(tmpVtkfile)
            os.remove(tmpVtkfile)
        else:
            raise IOError("Can't find vtkPath %s or zip-archive named %s" %
                          (vtkPath, os.path.dirname(vtkPath) + '.zip'))
        return vtkData

    def _getExtraFields_Prepare(self, stepNb, frameNb, points):
        """Part of the method getExtraFields has moved to this hidden method
        in order to avoid duplicate code: This part of the method would be
        identical in the derived class L{Converter_CombinedOutputOnSurface}.
        """
        try:
            csStep = self.stepFrameToCsStep[stepNb][frameNb]
        except (KeyError, IndexError):
            msg("ERROR: odb steps provided in stepFrameToCsStep: %s"
                "\nRequested step %d, frame %d."
                % (sorted(self.stepFrameToCsStep), stepNb, frameNb))
            raise
        if csStep<0:
            # no Cavesim output for current step/frame
            self.csData = False
            return None, None, None

        # read the Cavesim vtks
        fname = os.path.join(
            self.csOutputDir, "VELOCITY_VTK", "CAVESIM_VEL_%d.vtk" % csStep)
        msg("Reading %s" % fname)
        vtkVel = self._readVTKfromFileOrZipArch(fname)

        fname = os.path.join(
            self.csOutputDir, "FLOW_VTK", "CAVESIM_FLOW_%d.vtk" % csStep)
        msg("Reading %s" % fname)
        vtkFlow = self._readVTKfromFileOrZipArch(fname)

        # check that all Cavesim data is on the same grid
        if vtkVel.mesh != vtkFlow.mesh:
            raise ValueError(
                "Grids of the Cavesim output files differ.\nVEL: %s\nFLOW: %s"
                % (vtkVel.mesh, vtkFlow.mesh))

        # Cavesim data available
        self.csData = True

        # prepare Cavesim displacement fields
        dispFields = list()
        for fieldName in ["TOTAL_X_MAG", "TOTAL_Y_MAG", "TOTAL_Z_MAG"]:
            msg("Interpolating field %s" % fieldName)

            # convert displacement components: x>9998 => x=0
            fld = vtkVel.data[fieldName]
            fld = type(fld)(((x <= 9998) and x) or 0.0 for x in fld)

            # interpolate to output points
            fld = createFieldObject(fieldName, "point", "scalar", initArgs=[
                fld.interpolateToPoints(vtkVel.mesh, points),])

            dispFields.append(fld)

        # calculate threshold
        self.csVelThreshold = min(
            x for x in vtkVel.data["TOTAL_MAG"] if x > 0) / 2

        fieldName = "UMAG_CS"
        msg("Calculating %s" % fieldName)
        fld = createFieldObject(fieldName, "point", "scalar", initArgs=[
            # find None values and replace by 0.0
            (((x is not None) and (x**2+y**2+z**2)**0.5) or 0.0
             for x,y,z in izip(*dispFields)),])
        dispFields.append(fld)

        return vtkVel, vtkFlow, dispFields

    def getExtraFields(self, stepNb, frameNb, points, **dummy):
        """Method to read extra data from Cavesim output vtk files.
        This method must be supplied to the parameter getExtraFields of
        functions L{odbToVtkGridBox} and the like.
        """
        # part of the code has been out-sourced to a common hidden method
        vtkVel, vtkFlow, dispFields = self._getExtraFields_Prepare(
            stepNb, frameNb, points)
        if vtkVel is None:
            return
        for fld in dispFields:
            yield fld

        # there is only 0,1,2 values, air: <1, cave: >=1
        fieldName = "CellParticleState"
        msg("Interpolating field %s" % fieldName)
        fld = vtkFlow.data[fieldName]
        fld = createFieldObject(fieldName, "point", "scalar", initArgs=[(
            # convert None to -1 because storeField converts Nones to 0.0
            ((x is None) and -1.0 or x)
            for x in fld.interpolateToPoints(vtkFlow.mesh, points)),])
        yield fld

    def getU_1(self, cls, dat, mesh):
        """Here is the original expression from Behrooz:
         >>> if uMag_CS>threshold:
         >>>     newField_U1.append(u_CS)
         >>> elif (state_CS<0.9 and state_CS is not None) or state_Abq == 0.0:
         >>>     newField_U1.append(0.0)
         >>> else:
         >>>     newField_U1.append(u_Abq)
        """
        if not self.csData:
            return cls(dat[self.varNameU % 1])

        threshold = self.csVelThreshold
        return cls((
            # how to read this:
            # (cond1 and val1) or (cond2 and val2) or val3
            # translates to:
            # if cond1 then val1 elif cond2 then val2 else val3
            (uMag_CS>threshold and u_CS)  # this is cond1 and val1
            or (uMag_CS<=threshold and (state_CS<0 or state_CS>=0.9)
                and state_Abq!=0.0 and u_Abq)  # this is cond2 and val2
            or 0.0)  # this is val3
            for uMag_CS, u_CS, state_CS, u_Abq, state_Abq in izip(
                    *map(dat.get, [
                        "UMAG_CS", "TOTAL_X_MAG", "CellParticleState",
                        self.varNameU%1, self.varNameState])))

    def getU_2(self, cls, dat, mesh):
        """Same as U_1
        """
        if not self.csData:
            return cls(dat[self.varNameU % 2])

        threshold = self.csVelThreshold
        return cls((
            # how to read this:
            # (cond1 and val1) or (cond2 and val2) or val3
            # translates to:
            # if cond1 then val1 elif cond2 then val2 else val3
            (uMag_CS>threshold and u_CS)  # this is cond1 and val1
            or (uMag_CS<=threshold and (state_CS<0 or state_CS>=0.9)
                and state_Abq!=0.0 and u_Abq)  # this is cond2 and val2
            or 0.0)  # this is val3
            for uMag_CS, u_CS, state_CS, u_Abq, state_Abq in izip(
                    *map(dat.get, [
                        "UMAG_CS", "TOTAL_Y_MAG", "CellParticleState",
                        self.varNameU % 2, self.varNameState])))

    def getU_3(self, cls, dat, mesh):
        """Same as U_1
        """
        if not self.csData:
            return cls(dat[self.varNameU % 3])

        threshold = self.csVelThreshold
        return cls((
            # how to read this:
            # (cond1 and val1) or (cond2 and val2) or val3
            # translates to:
            # if cond1 then val1 elif cond2 then val2 else val3
            (uMag_CS>threshold and u_CS)  # this is cond1 and val1
            or (uMag_CS<=threshold and (state_CS<0 or state_CS>=0.9)
                and state_Abq!=0.0 and u_Abq)  # this is cond2 and val2
            or 0.0)  # this is val3
            for uMag_CS, u_CS, state_CS, u_Abq, state_Abq in izip(
                    *map(dat.get, [
                        "UMAG_CS", "TOTAL_Z_MAG", "CellParticleState",
                        self.varNameU % 3, self.varNameState])))

    def getLogP(self, cls, dat, mesh):
        """Here is the original expression from Behrooz:
         >>> if uMag_CS>threshold:
         >>>     newField_logP.append(2.0)
         >>> elif (state_CS<0.9 and state_CS is not None) or state_Abq == 0.0:
         >>>     newField_logP.append(0.0)
         >>> else:
         >>>     newField_logP.append(logP_Abq)

        ... and the formula for logP is:
         >>> if PST>0 then log10(PST*1000+1) else 0.0
        """
        if not self.csData:
            return Converter_StandardVtk.getLogP(self, cls, dat, mesh)

        threshold = self.csVelThreshold
        return cls((
            # how to read this:
            # (cond1 and val1) or (cond2 and val2) or val3
            # translates to:
            # if cond1 then val1 elif cond2 then val2 else val3
            (uMag_CS>threshold and 2.0)  # this is cond1 and val1
            or ((state_CS<0 or state_CS>=0.9) and state_Abq!=0.0
                and damage_Abq>0 and log10(damage_Abq*1000+1))  # cond2 and val2
            or 0.0)  # this is val3
            for uMag_CS, state_CS, damage_Abq, state_Abq in izip(
                    *map(dat.get, [
                        "UMAG_CS", "CellParticleState",
                        self.varNameDamage, self.varNameState])))

    def statusFunc(self, uMag_CS, state_CS, state_Abq):
        """This is the original expression from Behrooz.
        ("state_CS>=0" was "state_CS is not None" originally)
        """
        threshold = self.csVelThreshold
        if (state_CS<0.99 and state_CS>=0):
            return 0.0
        elif (state_Abq < 0.5 or state_Abq > 5.01):
            return 0.0
        elif uMag_CS>threshold or state_Abq>=4.5:
            return 4.0
        elif (state_Abq>=2.5 and state_Abq<3.5):
            return 1.0
        elif (state_Abq>=0.5 and state_Abq<1.5):
            return 2.0
        elif (state_Abq>=3.5 and state_Abq<4.5):
            return 3.0
        else:
            return 2.0

    def getState(self, cls, dat, mesh):
        """
        Creates those values::
          0 : AIR - where CaveSIM is air (component_2<0.9) and Abq is air
              (component_7=0)
          1 : void - this will be fixed after correct interpolation technique
              is implemented in creating VTKs
          2 : rock
          3 : fill - this will be fixed after correct interpolation technique
              is implemented in creating VTKs
          4 : CAVE - where (uMag_CS>threshold)
        """
        if not self.csData:
            return Converter_StandardVtk.getState(self, cls, dat, mesh)

        return cls(
            self.statusFunc(uMag_CS, state_CS, state_Abq)
            for uMag_CS, state_CS, state_Abq in izip(*map(dat.get, [
                    "UMAG_CS", "CellParticleState", self.varNameState])))


class Converter_CombinedOutputOnSurface(Converter_CombinedVtk_Old):
    """Variant of L{Converter_CombinedVtk} suitable for generating output
    at the top surface, like surface subsidence.

    Points sitting in a cell that is air (CellParticleState==0) will be
    supplied with the value of the first non-air cell below.

    All other formulas and processes are the same as for
    L{Converter_CombinedVtk}.

    Usage:
     >>> import os.path
     >>> from bae.utils.odbToPointsData_02 import \\
     >>>     odbVariableFaceCentroidsToCsvPoints, \\
     >>>     Converter_CombinedOutputOnSurface
     >>>
     >>> # from outputVarNames.py
     >>> varNameDamage = "SDV1"
     >>> varNameState = "SDV8"
     >>>
     >>> converter = Converter_CombinedOutputOnSurface(
     >>>     stepFrameToCsStep={1: [-1,]*8 + range(1,73)},,
     >>>     csOutputDir="../../2_RUN/R02b_RST02/CavesimOutput",
     >>>     varNameU="U_%d",
     >>>     varNameDamage=varNameDamage, varNameState=varNameState)
     >>>
     >>> odbVariableFaceCentroidsToCsvPoints(
     >>>     odbPath, frameList=frameList,
     >>>     initialSurface="top", seqSetsNameExcav="seqElsetsPit",
     >>>     outputDirName=outputDirName, fileNamePrefix=fileNamePrefix,
     >>>     fieldNames=["U_1", "U_2", "U_3", varNameDamage,
     >>>                 (varNameState, "const")],
     >>>     getExtraFields=converter.getExtraFields,
     >>>     outputFields=[
     >>>         ("U_1", "scalar", converter.getU_1),
     >>>         ("U_2", "scalar", converter.getU_2),
     >>>         ("U_3", "scalar", converter.getU_3),
     >>>         ("logP", "scalar", converter.getLogP),
     >>>     ])
    """

    @staticmethod
    def interpolateOnSurface(oldField, filterField, oldMesh, points):
        """Interpolate values from Cavesim (oldData on oldMesh) onto given
        points in a piecewise constant manner. If the corresponding grid point
        (in Cavesim output) doesn't have valid data then take that of the next
        point vertically below that has proper data. This is a generator
        function suitable as initializer for a list/field-object

        @param oldField: Field(pos="structPt") object holding the data to be
           interpolated
        @param filterField: Field is True in cells that contains valid data.
           Otherwise try the cell below...
        @param oldMesh: MeshStructuredPoints object for oldField (and
           filterField)
        @param points: iterable containing the point coordinates where
           interpolated values are seeked for.
        """

        stride = oldMesh.strides[2]  # abbreviation

        # list of point coordinates/L{MeshUnstructuredPoints}
        ptIdsWeights = oldMesh.getPtIdsWeights(points)
        for ptIds, weights in ptIdsWeights:
            if ptIds is None:
                yield None
            else:

                # index of the largest element coordinate / closest corner
                # Note: in pure python this method is faster than:
                # max(enumerate(coords), key=lambda x:x[1])[0]
                # max(ids, key=coords.__getitem__)
                # max((c,i) for i,c in enumerate(coords))[1]
                # max((c,i) for c,i in izip(coords, ids))[1]
                idx = weights.index(max(weights))

                # if not filterField go down one cell
                while not filterField[idx] and idx>=stride:
                    idx -= stride
                yield oldField[idx]

    def getExtraFields(self, stepNb, frameNb, points, **dummy):
        """Method to read extra data from Cavesim output vtk files.
        This method must be supplied to the parameter getExtraFields of
        functions L{odbToVtkGridBox} and the like.
        """
        # part of the code has been out-sourced to a common hidden method
        vtkVel, vtkFlow, dispFields = (
            Converter_CombinedVtk._getExtraFields_Prepare(
                self, stepNb, frameNb, points))
        if vtkVel is None:
            return
        for fld in dispFields:
            yield fld

        # prepare filter field: True if not air
        msg("Preparing filter field.")
        fld = vtkFlow.data["CellParticleState"]  # air: 0
        airFilter = type(fld)(x>0.5 for x in fld)

        # there is only 0,1,2 values, air: <1, cave: >=1
        fieldName = "CellParticleState"
        msg("Interpolating field %s" % fieldName)
        fld = vtkFlow.data[fieldName]
        fld = createFieldObject(fieldName, "point", "scalar", initArgs=[
            self.interpolateOnSurface(fld, airFilter, vtkFlow.mesh, points),])
        yield fld


class Converter_StandardAndFOS(Converter_StandardVtk):
    """Variant of L{Converter_StandardVtk} suitable for generating output
    for local factor of safety and damage variable (unified pst dependent
    strength ratio).

    The different types of FOS will be calculated with respect to local
    material type, damage (pst), (fill-)state and anisotropic orientation.

    Special settings/arguments for FOS calculation can be set via
    param-dictionarys (ivars, see)

    Usage:
        >>> from bae.utils.odbToPointsData_02 import odbToVtkGridBox
        >>> from bae.utils.odbToPointsData_02 import Converter_StandardAndFOS
        >>> from bae.misc_01 import Config
        >>> config = Config("configBE.py")
        >>> matPropsInps = [os.path.join(os.path.dirname(config.odbPath),
        ...                 '%(projectPrefix)s_MATPROPS_%(matVersion)s.inp'
        ...                  % vars(config))]
        >>> converter = Converter_StandardAndFOS(matPropsInps,
        ...                 isotropic=False,
        ...                 varNameCohsNormals=config.varNameCohsNormals,
        ...                 varNameDamage=config.vumatVerDamage,
        ...                 varNameState=config.vumatVerSTATE,
        ...                 varNameTimeFill=config.varNameTimeFill
        ...                 )
        >>> fieldNames = converter.requiredFieldNames
        >>> # something special for FOSMohrCoulomb here:
        >>> #     - use a fixed frictionAngle=35deg
        >>> #     - force cohesion to pass UCS
        >>> #     - compare against intact rock (>> pst = 0)
        >>> converter.paramMC['phi'] = 35.0
        >>> converter.paramMC['constrainCohesionToUCS'] = True
        >>> converter.paramMC['pst'] = 0.0
        >>> outputFields = [
        ...     ("PST", "scalar", config.vumatVerDamage),
        ...     ("STATE", "scalar", converter.getState),
        ...     ("DAMAGEVAR","scalar", converter.damageVariable),
        ...     ("MATERIAL", "scalar", "Material"),
        ...     ("FOS_MC_IntactCohs_phi35", "scalar",converter.fosMohrCoulomb),
        ...     ("FOS_MISES", "scalar", converter.fosMises),
        ...     ("FOS_P", "scalar", converter.fosHydrostatic),
        ...     ]
        >>> boxName = 'myFOSBox'
        >>> fileNamePrefix = "_".join((config.projectPrefix, config.runVersion,
        ...                            config.seqVersion, boxName))
        >>> odbToVtkGridBox(
        ...     odbPath=config.odbPath,
        ...     frameList=config.myFrameList,
        ...     fieldNames=fieldNames,
        ...     outputFields=outputFields,
        ...     matRegionsOption=matRegionsOption,
        ...     getExtraFields=converter.getExtraFields,
        ...     gridName=boxName,
        ...     gridData=config.gridDataDict[boxName],
        ...     outputDirName='./VTK',
        ...     fileNamePrefix=fileNamePrefix,
        ...     projectPrefix=config.projectPrefix,
        ...     meshVersion=config.meshVersion,
        ...     abaqusVersion="abq6132",
        ...     outputFormat="binary",
        ...     interpolationDataDir='./',
        ...    )
    """

    #: Parameters for MohrCoulomb approximation of FOS. See
    #: L{bae.material_01.material.FactorOfSafetyIso.fosMohrCoulomb<fosMohrCoulomb>}
    paramMC = dict(
        pst=None,  # if not None: sets a fixed value for pst-evaluation
        c=None,    # if not None: sets a fixed cohesion
        phi=None,  # if not None: sets a fixed friction angle

        # if True and c is None: forces fit to match UCS
        constrainCohesionToUCS=False,
        fittingRange=[0,None],        # restricts fitting range of p
        measureFromBisection=True,    # measures from S1=S3 instead of S1=0
        )

    #: Parameters for HoekBrown approximation of FOS. See
    #: L{bae.material_01.material.FactorOfSafetyIso.fosHoekBrown<fosHoekBrown>}
    paramHB = dict(
        pst=None,  # if not None: sets a fixed value for pst-evaluation
        measureFromBisection=True,   # measures from S1=S3 instead of S1=0
        )

    #: Parameters for Mises approximation of FOS.
    paramMises = dict(
        pst=None,  # if not None: sets a fixed value for pst-evaluation
        )

    #: Parameters for Hydrostatic pressure approximation of FOS.
    paramHydrostatic = dict(
        pst=None,  # if not None: sets a fixed value for pst-evaluation
        )

    #: infinite fos values will set to this value
    infValue = 10

    #: UCSi threshold value that assumes matierial to be elastic --> skippes fos
    UCSiElastic = 500.E6

    #: offset value for (var)FillTime
    varFillOffset = 1000

    def __init__(self, matInputFiles, **kwargs):
        """
        @param matInputFiles: single path of matInputdeck or list of those
        @kwarg varNameDamage: variable name for plastic strain, e.g. "SDV1"
        @kwarg varNameState: variable name for state, e.g. "SDV8"
        @kwarg varNameU: "U_%d" or "UR_%d" or something like that
        @kwarg fixedFillValue: default is None. If not None, filled cells will
            get this fos-value. No fos-calculation to fill/varFill will be
            applied.
        @kwarg initSDVsVTKPathDummy: default is None. For not-restart jobs the
            fill time and thus the varFill(offset) usually gets stored in the
            first (equilibrium) step and initSDVsVTKPathDummy is not required.
            For restarts you'll have to run odbToVtkGridBox_getFillTime.py first
            and set initSDVsVTKPathDummy accordingly.
        @kwarg varNameTimeFill: variable name for fillTime. Not required if
            initSDVsVTKPathDummy is not None or fixedFillValue is set.
        @kwarg vusubsParamPyPath: default is None. If not specified, the
            vusubs_param.py is assumed to be at odbPath/source/vusubs_param.py.
        @kwarg isotropic: default is True
        @kwarg varNameCohsNormals: variable name for (first) anistropic
            angle(s). Can be None/unset if isotropic is True or
            initSDVsVTKPathDummy is set.
        """
        ## init some attributes
        self.rockLookup = None
        self.fillLookup = None
        self.matField = None
        self.fillField = None
        self.anisoAlpha = None
        self.anisoBeta = None
        self.odbReader = None

        self._currentStepFrame = (None, None)
        self._lastStepFrame = (None, None)
        self._lastS1, self._lastS2, self._lastS3 = None, None, None
        self._lastPst, self._lastState = None, None
        self._lastMesh = None

        ## parse kwargs
        # required inputs
        self.matInputFiles = matInputFiles

        # we will pop the variables from kwargs to bypass the rest to
        # StandardConverter.__init__ later
        self.isotropic = kwargs.pop('isotropic', True)
        self.fixedFillValue = kwargs.pop('fixedFillValue', None)
        self.varNameTimeFill = kwargs.pop('varNameTimeFill', None)
        self.varNameCohsNormals = kwargs.pop('varNameCohsNormals', None)
        self.vusubsParamPyPath = kwargs.pop('vusubsParamPyPath', None)
        self.initSDVsVTKPathDummy = kwargs.pop('initSDVsVTKPathDummy', None)

        ## check inputs
        if self.initSDVsVTKPathDummy:
            raise NotImplementedError(
                'Getting FillTime and anisotropic angles from'
                ' existing vtk is not implemented/tested yet')

        if not self.isotropic and self.varNameCohsNormals is None:
            raise ValueError('varNameCohsNormals (e.g. "SDV15") required for'
                             ' anisotropic fos-calculation')

        if self.fixedFillValue is None and self.varNameTimeFill is None:
            raise ValueError(
                'Set a fixedFillValue to skip fos calculation for fill'
                ' or define the varNameTimeFill (e.g. "SDV13")')

        super(Converter_StandardAndFOS, self).__init__(**kwargs)


    @property
    def requiredFieldNames(self):
        """List of Fielnames required to calculate FOS-data.
        If material is isotropic (or anisoS, anisoN = 1, 1) only eigen stresses
        are required. An anisotropic material requires the complete stress
        tensor.

        Note that stress values should be 'interpolated' by using the piecwise
        constant method (takes nearest integration point) to avoid interpolation
        artefacts like negative FOS values.
        """
        if self.isotropic:
            fieldNames = [
                self.varNameDamage, self.varNameState,
                "S_MIN_PRINCIPAL", "S_MID_PRINCIPAL", "S_MAX_PRINCIPAL"]
        else:
            fieldNames = [
                self.varNameDamage, self.varNameState,
                "S_11", "S_22", "S_33", "S_12", "S_23", "S_13"]

        fieldNames = [(f, 'const') if f.startswith('S_') else f
                      for f in fieldNames]
        return fieldNames


    @property
    def standardFOSOutputFields(self):
        return [("REL_STRENGTHFRACTURE", "scalar", self.getDamageVariable),
                ("FOS_vonMISES", "scalar", self.getFosMises),
                ("FOS_HYDROSTATIC", "scalar", self.fosHydrostatic),]


    @property
    def fosType(self):
        if self.isotropic:
            return FactorOfSafetyIso
        else:
            return FactorOfSafetyAniso

    def createLookups(self):
        """Creates the
            - matCode to material lookup dict (self.rockLookup),
            - matCode field on grid (self.matField)
            - *varFillCode to fillMaterial lookup dict (self.fillLookup)
            - *varFill field on grid (self.fillField)
            - *anisotropy-normal-field (angles, self.anisoNormalsField)
        Data tagged with * will only be created if required.

        Usually, this function is just called once during the first
        call of L{getFosInput}.

        The rock and fill materials will be stored as
        L{PST-samples<bae.material_01.material.FactorOfSafetyIso>}
        in lookup dicts.
        """
        # init lookup-dicts holding matKey: material
        self.rockLookup, self.fillLookup = dict(), dict()

        # init matCode-field
        matNamesList, matField = self.odbReader.getMatFieldFromOdb()
        self.matField = np.round(matField).astype(int)

        # reading matprops from abaqus input
        mpc = MatPropsCollection('_tmp').importAbaqusInput(self.matInputFiles)
        rock = mpc.rock
        fill = mpc.fill  # 'default'-fill if varFill-offset is 0
        for rockName in rock:
            if rockName.upper() == rockName:
                continue
            rock.renameItem(rockName, rockName.upper())
            fill.renameItem(rockName, rockName.upper())

        # create rockMaterial lookup dict
        for ii, name in enumerate(matNamesList, start=1):
            msg('    adding %s to rockLookup' % name)
            fos = self.fosType(rock[name])
            setattr(fos, 'name', name)
            self.rockLookup[ii] = fos

            # assign defaultFills to fillLookup
            fos = self.fosType(fill[name])
            setattr(fos, 'name', name)
            self.fillLookup[ii] = fos

        # create (var)FillMaterial lookup dict
        if self.fixedFillValue is not None:
            # fill lookups not required
            return

        if self.vusubsParamPyPath is None:
            # assume path is odbPath + source/vusubs_param.py
            vusubs = os.path.join(os.path.dirname(self.odbReader.odbPath),
                                  'source', 'vusubs_param.py')
        else:
            vusubs = self.vusubsParamPyPath

        mpc = MatPropsCollection('_tmp')
        mpc.fill.importVusubsVarFill(vusubs)
        #note that varfill stored in mpc has the same order as found in
        #vusubs_param.py
        for ii, (name, varFillMat) in enumerate(mpc.fill.iteritems(), start=1):
            msg('    adding varFill %s to fillLookup' % name)
            fos = self.fosType(varFillMat)
            setattr(fos, 'name', name)
            self.fillLookup[ii*self.varFillOffset] = fos

        # derive varFillNumber from TimeFill-SDV
        if self.initSDVsVTKPathDummy:
            # getting TimeFill-SDV from existing VTK (required for restart)
            from bae.field_02 import FieldsCollection
            vtkPath = self.initSDVsVTKPathDummy % (self.gridName)
            flds = FieldsCollection.fromVtk(vtkPath)
            fillTime = flds['TimeFill']
        else:
            # getting TimeFill-SDV from equilibrium step
            msg('getting TimeFill from odb')
            step, frame = 1, 0
            fillTime = self.odbReader.getFieldFromOdb(
                step, frame, self.varNameTimeFill, interpolationType='const')
            fillTime = np.array(fillTime)

        self.fillField = self.matField.copy()   #first: default fill everywhere
        mask = (fillTime >= self.varFillOffset) #second: overwrite with varfill
        varFills = (fillTime[mask] // self.varFillOffset)*self.varFillOffset
        self.fillField[mask] = varFills.astype(int)

        if self.isotropic:
            return

        # get normal angles for anisotropy
        if self.initSDVsVTKPathDummy:
            # getting TimeFill-SDV from existing VTK (required for restart)
            from bae.field_02 import FieldsCollection
            vtkPath = self.initSDVsVTKPathDummy % (self.gridName)
            flds = FieldsCollection.fromVtk(vtkPath)
            self.anisoAlpha = flds['ANISO_ALPHA']
            self.anisoBeta = flds['ANISO_BETA']
        else:
            # getting TimeFill-SDV from equilibrium step
            step, frame = 1, 0
            alphaSDV = self.varNameCohsNormals
            try:
                betaSDV = 'SDV%d' % (int(alphaSDV.replace('SDV','')) + 1)
            except Exception as e:
                e.message += (' Expeceted something like "SDV$$" for'
                              ' varNameCohsNormals, but is %s'
                              % alphaSDV)
                msg('ERROR: %s' % e.message)
                raise e

            msg('getting anisotropic directions/angles from odb')
            self.anisoAlpha = np.array(self.odbReader.getFieldFromOdb(
                step, frame, alphaSDV, interpolationType='const'))
            self.anisoBeta = np.array(self.odbReader.getFieldFromOdb(
                step, frame, betaSDV, interpolationType='const'))


    def getExtraFields(self, stepNb, frameNb, points, odb, **dummy):
        """Enables access to the odbReader which is required to set up the
        material lookup (see L{createLookups}).
        The odbReader is stored in class attribute 'odbReader'.
        """
        self.odbReader = odb
        self._currentStepFrame = (stepNb, frameNb)
        return
        yield


    def getFosInput(self, dat, mesh):
        """Returns numpy arrays of principle stresses, pst, state and matcode.
        Make sure all fields are specified as fieldNames for odbToPointsData_02.
        """
        if self._lastMesh is not mesh:
            # mesh has changed --> new lookup required
            self.createLookups()

        if (self._lastMesh is mesh
                and self._lastStepFrame == self._currentStepFrame):
            # to avoid overhead in if getFosInput gets called multiple times
            # in same step.
            return (
                self._lastS1.copy(), self._lastS2.copy(), self._lastS3.copy(),
                self._lastPst.copy(), self._lastState.copy())

        pst = np.array(dat[self.varNameDamage]).astype(float)
        state = np.round(dat[self.varNameState]).astype(int)

        if self.isotropic:
            S3 = np.array(dat["S_MIN_PRINCIPAL"]).astype(float)
            S2 = np.array(dat["S_MID_PRINCIPAL"]).astype(float)
            S1 = np.array(dat["S_MAX_PRINCIPAL"]).astype(float)
        else:
            from bae.material_01.yieldsurfaces import rotVoigtTransversal
            # setup full stress tensor and apply the anisotropic
            # distorsion for rock material

            msg("ATTENTION! Anisotropy is only applied to rock material")
            S = np.empty((pst.shape[0], 3, 3))

            isRock = (state==1)

            # voxels that are not in rock material
            # only assing lower triangle of S as required for eigenvalsh
            S[:,0,0] = np.array(dat["S_11"])
            S[:,1,1] = np.array(dat["S_22"])
            S[:,2,2] = np.array(dat["S_33"])
            S[:,1,0] = np.array(dat["S_12"])
            S[:,2,1] = np.array(dat["S_23"])
            S[:,2,0] = np.array(dat["S_13"])

            # loop rock materials and apply anisotropic distorsion of S
            for matIdx in np.unique(self.matField[isRock]):
                material = self.rockLookup[matIdx]
                anisoS, anisoN = material.anisoS, material.anisoN
                ii = isRock & (self.matField == matIdx)
                if not ii.any():
                    continue

                msg(('Transforming StressTensor for %d sampling points'
                     ' with material %s') %(ii.sum(), material.name))
                try:
                    # get normal-vector angles
                    a,b = self.anisoAlpha[ii], self.anisoBeta[ii]

                    # apply rotation to stressTensor
                    V = [S[ii,0,0], S[ii,1,1], S[ii,2,2],
                         S[ii,1,0], S[ii,2,1], S[ii,2,0]]

                    v11, v22, v33, v12, v23, v31 = rotVoigtTransversal(V, a, b)

                    # scale rotated stess with anisotropic parameters
                    # only assing lower triangle of S as required for eigenvalsh
                    S[ii,0,0] = anisoN * v11
                    S[ii,1,1] = anisoN * v22
                    S[ii,2,2] = 1 * v33
                    S[ii,1,0] = anisoN * v12
                    S[ii,2,1] = anisoS * v23
                    S[ii,2,0] = anisoS * v31

                except Exception as e:
                    msg('ERROR: %s, %s' % (type(e), e.message))
                    raise

            try:
                # eigvalsh calculates eigenVals from Lower triangle (UPLO='L')
                # eigenValues already sorted sorted
                msg("getting eigen values of anisotropic deformed stress"
                    " tensor at %d sampling points" % len(isRock))
                E = np.linalg.eigvalsh(S, UPLO='L')
                S3, S2, S1 = E.T  # E.T. home phone

            except Exception as e:
                msg('ERROR during eigenvaluecalculation: %s, %s'
                    % (type(e), e.message))
                raise

            self._lastStepFrame = self._currentStepFrame
            self._lastS1, self._lastS2, self._lastS3 = S1, S2, S3
            self._lastPst, self._lastState = pst, state
            self._lastMesh = mesh

        return S1, S2, S3, pst, state

    def damageVariable(self, cls, dat, mesh):
        """Creates a field of pst-dependent damage variable of rock mass which
        is defined as::
           damVar = (UCS(PST) - UCS_Res) / (UCS_Peak - UCS_Res)

        which ranges from 1 (intact strength) to 0 (resudial strength). See
        L{bae.material_01.DefaultPst.getDamageVariable} for more details. Fill
        will be flagged with value -1.

        @note: poorly tested yet
        """
        if self._lastMesh is not mesh:
            # mesh has changed --> new lookup required
            self.createLookups()
            self._lastMesh = mesh

        msg('## calculating damage variable')

        # getting required fields
        pst = np.array(dat[self.varNameDamage]).astype(float)
        state = np.round(dat[self.varNameState]).astype(int)

        # initialize output field
        damVar = np.full_like(pst, -1)
        isRock = (state==1)
        isFill = (state==4)
        damVar[isRock | isFill] = np.nan  # make sure that we got everything

        for matIdx in np.unique(self.matField[isRock]):
            material = self.rockLookup[matIdx]
            ii = isRock & (self.matField == matIdx)
            if not ii.any():
                continue
            damVar[ii] = material.getDamageVariable(pst[ii])

        if self.fixedFillValue is not None:
            damVar[isFill] = -.5
            return cls(damVar)

        for fillIdx in np.unique(self.fillField[isFill]):
            material = self.fillLookup[fillIdx]
            ii = isFill & (self.fillField == fillIdx)
            if not ii.any():
                continue
            damVar[ii] = material.getDamageVariable(pst[ii])

        return cls(damVar)


    def _fosMain(self, dat, mesh, pst, fosFun, miningNotation):
        """FOS-'workflow' is the same for all implemented variants of
        FOS. Just the actual calculation, done by fosFun, differs.

        @param dat: odb-fields as provided for each converter function
        @param pst: plastic strain. Can be a fixed value or a (odb)field
        @param fosFun: wrapped fos-function that takes
            material, S1, S2, S3, pst
        @param miningNotation: forces reversed sorted negative S1, S2, S3
            -values to be passed to fosFun
        """

        S1, S2, S3, pStrain, state = self.getFosInput(dat, mesh)

        if pst is None:
            pst = pStrain  # evaluate pst dependent yieldSurfaces

        if miningNotation:
            S1, S2, S3 = -S3, -S2, -S1

        fos = np.zeros_like(S1)
        isRock = (state==1)
        isFill = (state==4)
        fos[isRock | isFill] = np.nan  # make sure that we got everything

        ## fos for rock material
        for matIdx in np.unique(self.matField[isRock]):
            material = self.rockLookup[matIdx]
            ii = isRock & (self.matField == matIdx)

            if not ii.any():
                continue

            if material.UCSi < self.UCSiElastic:
                msg('    - calculating fos for %d positions in rock %s'
                    % (ii.sum(), material.name))
                try:                # to get pstField
                    _pst = pst[ii]  # --> pst is a field
                except TypeError:   # fails for a scalar float value
                    _pst = pst      # --> pst is a fixed value

                matFos = fosFun(material, S1[ii], S2[ii], S3[ii], _pst)
            else:  # is elastic rock --> set fos = infValue
                msg('    - rock %s is elastic (UCSi>=%.1f)'
                    % (material.name, self.UCSiElastic))
                matFos = self.infValue

            fos[ii] = matFos

        ## fos for fill material
        if self.fixedFillValue is not None:
            # sets a fixed value to sampling points in fill
            fos[isFill] = self.fixedFillValue
            return fos

        for fillIdx in np.unique(self.fillField[isFill]):
            material = self.fillLookup[fillIdx]
            ii = isFill & (self.fillField == fillIdx)
            if not ii.any():
                continue

            if material.UCSi < self.UCSiElastic:
                msg('    - calculating fos for %d positions in fill %s'
                    % (ii.sum(), material.name))
                try:                # to get pstField
                    _pst = pst[ii]  # --> pst is a field
                except TypeError:   # fails for a scalar float value
                    _pst = pst      # --> pst is a fixed value

                matFos = fosFun(material, S1[ii], S2[ii], S3[ii], _pst)
            else:  # is elastic fill --> set fos = infValue
                msg('    - fill %s is elastic (UCSi>=%.1f)'
                    % (material.name, self.UCSiElastic))
                matFos = self.infValue

            fos[ii] = matFos
        return fos


    def fosMises(self, cls, dat, mesh):
        """FOS that compares the vonMises stress of local stress values
        and their related vonMises stress on the (current) yieldSurface.

        This fos focuses on the influence of shear deformation.

        See L{fosMises<bae.material_01.material.FactorOfSafetyIso.fosMises>}
        for more information. Define further parameters via instance variable
        L{paramMises}.

        @note: This method is time consuming because the vonMises stress
            on yield surface has to be solved iteratively. If you have many
            points and many materials this may become quite painfull.
        """
        msg('## calculating fosMises')
        ## get parameters
        pst = self.paramMises['pst']

        def fosFun(mat, s1, s2, s3, p):
            return mat.fosMises(s1, s2, s3, p, infValue=self.infValue)
        try:
            fos = self._fosMain(dat, mesh, pst, fosFun, False)
        except Exception as e:  # write error and traceback to log
            msg('Error: %s' % type(e))
            msg(traceback.format_exc())
            msg(e.message)
            raise
        return cls(fos)


    def fosHydrostatic(self, cls, dat, mesh):
        """FOS that compares the hydrostatic pressure of local stress values
        and their related pressure on the (current) yieldSurface.

        This fos focuses on the influence of (loose of) confinement.

        See L{fosHydrostatic<bae.material_01.material.FactorOfSafetyIso.fosHydrostatic>}
        for more information. Define further parameters via instance variable
        L{paramHydrostatic}.
        """
        msg('## calculating fosHydrostatic')
        ## get parameters
        pst = self.paramHydrostatic['pst']

        def fosFun(mat, s1, s2, s3, p):
            return mat.fosHydrostatic(s1, s2, s3, p, infValue=self.infValue)

        try:
            fos = self._fosMain(dat, mesh, pst, fosFun, False)
        except Exception as e:  # write error and traceback to log
            msg('Error: %s' % type(e))
            msg(traceback.format_exc())
            msg(e.message)
            raise
        return cls(fos)


    def fosHoekBrown(self, cls, dat, mesh):
        """FOS from HoekBrown approximation. This will call
        L{fosHoekBrown<bae.material_01.FactorOfSafetyIso.fosHoekBrown>}. Define
        further parameters via instance variable L{paramHB}.
        """
        msg('## calculating fosHoekBrown')
        ## get parameters
        pst = self.paramHB['pst']
        measureFromBisection = self.paramHB['measureFromBisection']

        def fosFun(mat, s1, s2, s3, p):
            return mat.fosHoekBrown(s1, s3, p,
                                    measureFromBisection=measureFromBisection,
                                    infValue=self.infValue)

        try:
            fos = self._fosMain(dat, mesh, pst, fosFun, True)
        except Exception as e:  # write error and traceback to log
            msg('Error: %s' % type(e))
            msg(traceback.format_exc())
            msg(e.message)
            raise
        return cls(fos)


    def fosMohrCoulomb(self, cls, dat, mesh):
        """FOS from MohrCoulomb approximation. This will call
        L{fosMohrCoulomb<bae.material_01.FactorOfSafetyIso.fosMohrCoulomb>}.
        Define further parameters via instance variable L{paramMC}.

        Example for outputFields used in odbToVtk and others
         >>> converter = Converter_StandardAndFOS('matprops.inp')
         >>> converter.paramMC['phi'] = 35.0
         >>> converter.paramMC['constrainCohesionToUCS'] = True
         >>> outputFields = [
         ...    ("S3", "scalar", converter.getS3),
         ...    ("FOS_MC35", "scalar", converter.fosMohrCoulomb),
         ...    ]

        @note: If you don't specify phi and c (or constrain c to UCS) this
            method might be very slow because when evaluated on many points
            and materials. In a first step the true yieldSurface has to be
            evaluated and in a second step this will be fitted. This is fast
            for a few points but may be horrible for 25Mio samples.
        """
        msg('## calculating fosMohrCoulomb')

        ## get parameters
        pst = self.paramMC['pst']
        c   = self.paramMC['c']
        phi = self.paramMC['phi']
        constrainCohesionToUCS = self.paramMC['constrainCohesionToUCS']
        fittingRange = self.paramMC['fittingRange']
        measureFromBisection = self.paramMC['measureFromBisection']

        def fosFun(mat, s1, s2, s3, p):
            return mat.fosMohrCoulomb(
                s1, s3, p,
                pstMax=None, infValue=self.infValue,
                miningConvention=True, c=c, phi=phi,
                constrainCohesionToUCS=constrainCohesionToUCS,
                fittingRange=fittingRange,
                measureFromBisection=measureFromBisection)

        try:
            fos = self._fosMain(dat, mesh, pst, fosFun, True)
        except Exception as e:  # write error and traceback to log
            msg('Error: %s' % type(e))
            msg(traceback.format_exc())
            msg(e.message)
            raise
        return cls(fos)
