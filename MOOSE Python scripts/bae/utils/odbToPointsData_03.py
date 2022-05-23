"""Read interpolated field data from the specified odb and export in different
file formats. Exports regular grid data to legacy vtk files.


Usage:
======

For points on a regular grid:
 >>> from bae.utils.odbToPointsData_03 import odbToVtkGridBox
 >>> odbToVtkGridBox(...)

Or for unstructured points
 >>> from bae.utils.odbToPointsData_03 import odbToCsvPoints
 >>> odbToCsvPoints(...)
"""

__version__ = "3.01"

_version_history_ = """
3.01 GP new based on 2.35
  - switched from odb_03 to odb_04
  - odbStep identified by name (string)
  - inputFieldsPerFrame argument with interpolation, maximize and
    relative options
  - split off odbToPointsConverter_03.py
  - framelist states only output frames, fieldMaximizeAfter mandatory
"""

todo = """
todo for version ..._03:
========================

- test reading materialdata from postData

- clean up stepframe-to-name:
  . use modern sting.format() syntax so we can use postData.sequence.frameNames
  . or remove string-template option alltogether?
  . always supply framename from postdata
- take matNamesList from postdata as a list of (number, name)-tuples
- clean way to get constant and time varying input from other sources as well:
  . FS4
  . from precomputed vtks (like reference displacements)
  . new time constant fields from composer or so...

- move getExtrafields as attribute into inputFieldsPerFrame argument
- new inputFieldsConst argument like inputFieldsPerFrame with getExtraFields-function
- rename getExtraFields


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

ideas for version ..._04:
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
    exception easily differentiable from a possible error while storing...
  - this option is connected to the new suggests for outputFields.
    Maybe getEtraFields is not needed anymore or outputFields must be modified
"""

import os
import csv
import re
from itertools import izip
from collections import defaultdict

from bae.future_01 import OrderedDict
from bae.misc_01 import RowIterFromCsvFile
from bae.mesh_01 import getMesh, MeshUnstructuredPoints, Mesh
from bae.field_01 import createFieldClass, createFieldObject
from bae.vecmath_01 import vector, vector_plus, vector_scale, \
    centre, mat_multvec, dot, norm, length
from bae.odb_04 import OdbReaderInterpolateToPoints, OdbReader, getOdbPostData
from bae.vtk_02 import FieldsOnStructPointGrid, \
    FramesFieldsOnUnstructuredPoints
from bae.surface_03 import ElemFaceSurface
from bae.makeHistory_03 import sequenceElsetsIter
from bae.log_01 import msg, MsgTicker


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
     >>>     def __init__(self, odbPath, inputFieldsPerFrame,...
     >>>                  ...lots of arguments ...):
     >>>         ... some checks ...
     >>>         self.storeArgs(locals(), "odbPath", ...)
     >>>     ... more overloaded functions
     >>> def myOdbToSomething(odbPath, inputFieldsPerFrame,...
     >>>                      ...lots of arguments ...):
     >>>     reader = MyOdbToSomethingClass(
     >>>         odbPath, inputFieldsPerFrame, ...args...)
     >>>     reader.run()
     >>>
     >>> myOdbToSomething(config.odbPath, ["U_1", "U_2"], ....)

    @ivar fieldNames: list of field name strings of fields to be read from the
       odb
    @ivar fieldInterpolType: {field name: type}, type being "default" or
       "const" for piecewise constant
    @ivar fieldNamesMaximize: set of field name strings containing all fields
       to be maximized per output frame over intermediate frames
    @ivar relativeFieldNames: list of field name strings containing all fields
       to be taken relative to a reference value.
    @ivar outputTopo: grid/points/topo for the output data. Typically a
       L{MeshStructuredPoints}- or L{MeshUnstructuredPoints}-instance created
       by the L{initOutput}() function or the L{initOutputForFrame}() function
    """
    def __init__(
            self, inputFieldsPerFrame,
            fieldMaximizeAfter=None,
            relativeToFrame=("Step-2",0),
            ):
        """General and very basic processing of arguments:
        Some checks and storing as instance variables.

        Processes the inputFieldsPerFrame list with (name, interpol-type,
        maximize-flag, relative-flag)-tuple items. Thereby creates the list of
        field names self.fieldNames, the dict self.fieldInterpolType {field
        name: interpolation type}, the set self.fieldNamesMaximize and the list
        relativeFieldNames.

        Checks that fieldMaximizeAfter is a (string,int)-tuple. (Also not
        stored here!)
        """

        # for checking create a set of all fieldname - strings
        self.fieldNames = list()
        self.fieldInterpolType = dict()
        self.fieldNamesMaximize = set()
        self.relativeFieldNames = list()
        for field in inputFieldsPerFrame:
            # for simplicity of further processing: make sure field is a tuple
            if isinstance(field, basestring):
                field = (field, )

            # store data in different containers
            fieldName = field[0]
            self.fieldNames.append(fieldName)  # field name

            try:  # interpolation type
                self.fieldInterpolType[fieldName] = field[1]
            except IndexError:
                self.fieldInterpolType[fieldName] = "default"

            try:  # maximize flag
                if field[2]:
                    self.fieldNamesMaximize.add(fieldName)
            except IndexError:
                pass

            try:  # relativeField
                if field[3]:
                    self.relativeFieldNames.append(fieldName)
            except IndexError:
                pass

        fieldNamesSet = set(self.fieldNames)
        msg("OdbToPointsData.__init__: got fieldNames: %s, %d times non-default"
            " interpolation-type."
            % (self.fieldNames,
               sum(1 for x in self.fieldInterpolType.itervalues()
                   if x != "default")), debugLevel=10)

        # check fieldMaximizeAfter argument
        if self.fieldNamesMaximize and fieldMaximizeAfter:
            try:
                if not(isinstance(fieldMaximizeAfter[0], basestring)
                       and isinstance(fieldMaximizeAfter[1], int)):
                    raise TypeError()
            except (TypeError, IndexError):
                raise ValueError(
                    "Incorrect value for fieldMaximizeAfter: %s"
                    % fieldMaximizeAfter)
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
        # a (odbStep, odbFrame) - tuple.
        # On the other hand all bae.utils/odbToPointsData - methods
        # expect for their stepFrameToName argument a function taking two
        # arguments: odbStep, odbFrame
        if callable(self.stepFrameToName):
            # strings don't need conversion, but functions do
            stepFrameToNameTwoArgs = self.stepFrameToName
            def stepFrameToName(frameId):
                return stepFrameToNameTwoArgs(*frameId)

        # if argument stepFrameToName=="fromPostData": define special function
        elif self.stepFrameToName=="fromPostData":
            # special case, frameNames from postData
            def stepFrameToName(frameId):
                try:
                    frameNames = self.odb.odbReader.postData.sequence.frameNames
                except AttributeError:
                    raise ValueError(
                        "frameNames not defined in associated XXX_postData"
                        " database for odbPath=%s." % self.odbPath)
                frameName = frameNames[frameId[0]][frameId[1]]
                msg("frameName from postData: step %s, frame %d, frameName %s"
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
        if (self.frameList[0][0] != self.frameList[-1][0]
            or self.frameList[0][0] != self.fieldMaximizeAfter[0]):
            raise NotImplementedError(
                "Different steps with fieldNamesMaximize option not implemented"
                " (yet)!\nframeList: %s" % self.frameList)
        frameListAll = [
            (self.frameList[0][0], odbFrame)
            for odbFrame in range(self.fieldMaximizeAfter+1,
                                  self.frameList[-1][1]+1)]
        frameList = set(self.frameList)
        for odbStep, odbFrame in frameListAll:
            yield (odbStep, odbFrame, (odbStep, odbFrame) in frameList)

    def frameStepperStandard(self):
        for odbStep, odbFrame in self.frameList:
            yield (odbStep, odbFrame, True)

    def openOdb(self):
        """prepare: open odb
        Initialize instance variable self.odb.
        """
        msg("Opening the odb %s." % self.odbPath)
        odb = OdbReader(self.odbPath, version=self.abaqusVersion)
        odb.openOdb()
        self.odb = OdbReaderInterpolateToPoints(odb)

    def initOutput(self):
        """Initialize output.
         1. check output directory
         2. preprocess self.outputFields and self.outputFieldsConst:
            If outputFields is given then convert it into list of tuples:
             - if composer --third item of an outputFields-tuple-- is callable
               then create (field-class, composer function) tuples
             - if composer is string create (output field name, odb field name)
               tuples
            If outputFieldsConst is given then also convert it in the same way.

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

        # for outputFields and outputFieldsConst arguments:
        # if outputFields argument is given convert into list of tuples
        # if composer is callable create (field-class, composer function) tuples
        # if composer is string create (output field name, odb field name)
        # tuples
        def convertOutputFldItem(item):
            # outputFields are (fieldName, fieldType, composer)-tuples
            if isinstance(item, basestring):
                # if item is just a string
                return (item, item)
            elif callable(item[2]):
                return (
                    createFieldClass(
                        item[0], self.outputFieldType, item[1]), item[2])
            else:
                return (item[0], item[2])

        if self.outputFields:
            try:
                diagTxt = ",\n".join(
                    ("- fld name %s from %s" % (
                        item[0], ("<func>" if callable(item[2]) else item[2]))
                     if isinstance(item, tuple)
                     else ("- fld %s" % item))
                    for item in self.outputFields)
            except IndexError:
                raise ValueError(
                    "outputFields argument accepts strings or tuples of three"
                    " (!) values. Example: ('matCode', 'scalar', 'Material')")
            msg("Before converting outputFields:\n%s" % diagTxt, debugLevel=10)
            self.outputFields = [
                convertOutputFldItem(item)
                for item in self.outputFields]
            msg("After converting outputFields:\n%s" % ",\n".join(
                "- output fld %s from odb fld %s"
                % (item[0], (callable(item[1]) and "<func>" or item[1]))
                for item in self.outputFields), debugLevel=10)
        else:
            msg("OutputFields parameter not given, will export all fields in"
                " fieldNames parameter without modification.")
            self.outputFields = [(x, x) for x in self.fieldNames]

        if self.outputFieldsConst:
            try:
                diagTxt = ",\n".join(
                    ("- fld name %s from %s" % (
                        item[0], ("<func>" if callable(item[2]) else item[2]))
                     if isinstance(item, tuple)
                     else ("- fld %s" % item))
                    for item in self.outputFieldsConst)
            except IndexError:
                raise ValueError(
                    "outputFieldsConst argument accepts strings or tuples of"
                    " three (!) values. Example: ('matCode', 'scalar',"
                    " 'Material')")
            msg("Before converting outputFieldsConst:\n%s" % diagTxt,
                debugLevel=10)
            self.outputFieldsConst = [
                convertOutputFldItem(item)
                for item in self.outputFieldsConst]
            msg("After converting outputFieldsConst:\n%s" % ",\n".join(
                "- output fld %s from odb fld %s"
                % (item[0], (callable(item[1]) and "<func>" or item[1]))
                for item in self.outputFieldsConst), debugLevel=10)
        else:
            self.outputFieldsConst = []
            if self.matRegionsOption:
                msg("MatRegionsOption is set and will induce a time constant"
                    " vtk output file.")
                self.outputFieldsConst.append(("Material", "Material"))

    def getMatRegions(self):
        """Read mat regions from the odb, write mat-names-csv file and return
        mat-number-field.

        @returns: returns a field of a type that you get from this call:
           createFieldClass("Material", "structPt", "scalar")
        """

        # create matNameToMatNumber dict from postdata
        if hasattr(self, "matNumberList"):
            matNameToMatNumber = dict((n,i) for i,n in self.matNumberList)
            firstMatCode = max(i for i,n in self.matNumberList) + 1
        else:
            try:
                self.matNumberList = (
                    self.odb.odbReader.postData.material.matNumberList)
            except AttributeError:
                msg("WARNING: No material number list in the post-data"
                    " database. Will make up new numbers.")
                matNameToMatNumber = None
                self.matNumberList = list()
                firstMatCode = 1
            else:
                matNameToMatNumber = dict((n,i) for i,n in self.matNumberList)
                firstMatCode = max(i for i,n in self.matNumberList) + 1

        # get the material data from the odb
        msg("Reading material properties from the odb.")
        (extraMatNamesList, matNumberField) = self.odb.getMatFieldFromOdb(
            nameFilter=self.matRegionsNameFilter,
            firstMatCode=firstMatCode,
            matNameToMatNumber=matNameToMatNumber,
            noMatCode=0)
        # converting to float
        matNumberField = type(matNumberField)(
            float(x) for x in matNumberField)

        # update matNumberList with data from the odb
        # store possibly updated material names list in instance attribute
        if self.matNumberList and extraMatNamesList:
            msg("WARNING: Found %d materials not listed in the specified"
                " matNamesList argument." % len(extraMatNamesList))
            msg("It's suggested to check spelling and the result.")
        self.matNumberList.extend(
            (i+firstMatCode, n) for i, n in enumerate(extraMatNamesList))

        return matNumberField

    def writeMatNamesListToFile(self):
        """Write the instance attribute self.matNumberList to a csv file.
        This function should only be called after (the last invocation of)
        L{getMatRegions}.
        """
        import csv
        outputPathMatList = os.path.join(
            self.outputDirName, self.fileNamePrefix+"_matList.csv")

        # write material number to names list
        msg("Writing material list to the csv file %s"
            % outputPathMatList)
        output = csv.writer(open(outputPathMatList, "wb"))
        output.writerow(["material number", "material name"])
        output.writerows(self.matNumberList)

    def storeConstField(self, field):
        """Store constant field, e.g. the matNumberField
        """
        raise NotImplementedError(
            "storeConstField is not defined for class %s." % self.__class__)

    def initRefData(self):
        """read (and store) reference data for relative fields
        Updates instance variable self.refData.
        """
        odbStep, odbFrame = self.relativeToFrame
        for fieldName in self.relativeFieldNames:
            msg("Reading reference field %s from step %s, frame %d."
                % (fieldName, odbStep, odbFrame))
            self.refData[fieldName] = self.odb.getFieldFromOdb(
                odbStep, odbFrame, fieldName,
                interpolationType=self.fieldInterpolType[fieldName])

    def processConstantFields(self):
        """Finally process constant fields before we come to the frames-loop.
        Overload this function in subclasses that need it.
        """
        return

    def initOutputForFrame(self, odbStep, odbFrame):
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

    def storeFrame(self, odbStep, odbFrame):
        """Store frame data to output. This will be called once for each frame
        after all fields have been processed by self.storeField().

        It also needs to process self.outputFields by means of
        L{self.processOutputFields}()
        """
        raise NotImplementedError(
            "storeFrame is not defined for class %s." % self.__class__)

    def processOutputFields(self, odbFieldContainer, outputFields):
        """Process self.outputFields for the current frame (fields from odb
        already read --passed in though odbFieldContainer argument).

        Generator function that yields field objects to be written to the
        output according to self.outputFields.

        @param odbFieldContainer: a dict containing the fields from the odb.
        @param outputFields: a list of output field descriptions. Pass
          either self.outputFields or self.outputFieldsConst.
        """
        outputGrid = self.odb.interpolation.grid
        for outputClassOrName, composer in outputFields:
            if callable(composer):
                msg("Processing callable composer for field %s."
                    % outputClassOrName.fieldName, debugLevel=10)
                field = composer(
                    cls=outputClassOrName,
                    data=odbFieldContainer,
                    topo=outputGrid,
                    postData=self.odb.odbReader.postData,
                    ctrl=self)
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
        for odbStep, odbFrame, isOutputFrame in frameIter:

            if isOutputFrame:
                self.initOutputForFrame(odbStep, odbFrame)
                currentFieldNames = self.fieldNames
            else:
                currentFieldNames = self.fieldNamesMaximize

            for fieldName in currentFieldNames:
                interpolationType = self.fieldInterpolType[fieldName]
                iptText = ""
                if interpolationType=="const":
                    iptText = " interpolating piecewise constant"
                msg("Reading field %s%s, step %s, frame %d."
                    % (fieldName, iptText, odbStep, odbFrame))
                try:
                    field = self.odb.getFieldFromOdb(
                        odbStep, odbFrame, fieldName,
                        interpolationType=interpolationType)
                except ValueError as exc:
                    if self.raiseMissingData:
                        raise
                    else:
                        msg("ERROR: Field %s for odb frame (%s, %d) not"
                            " available. Ending this loop now.\nError: %s"
                            % (fieldName, odbStep, odbFrame, exc))
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
                try:
                    self.getExtraFields(
                        ctrl=self, odbStep=odbStep, odbFrame=odbFrame,
                        topo=self.outputTopo, odb=self.odb)
                except IOError as exc:
                    if self.raiseMissingData:
                        raise
                    else:
                        msg("ERROR: Extra fields for odb frame (%s, %d) not"
                            " available. Ending this loop now.\nError: %s"
                            % (odbStep, odbFrame, exc))
                        retCode = 2
                        break

            # self.storeFrame evaluates self.outputFields
            # ... and possibly writes files
            if isOutputFrame:
                self.storeFrame(odbStep, odbFrame)

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
            odbPath,
            inputFieldsPerFrame,
            frameList,
            gridName, gridData,
            outputDirName, fileNamePrefix,
            projectPrefix, meshVersion,
            fieldMaximizeAfter=None,
            relativeToFrame=None, relativeToReferenceVtkPath=None,
            matRegionsOption=False, matRegionsNameFilter=None,
            outputFields=None,
            outputFieldsConst=None,
            getExtraFields=None,
            stepFrameToName=lambda s,f: "F%s%03d" % (s[-1],f),
            skipExisting=False,
            raiseMissingData=True,
            abaqusVersion=None,
            outputFormat="binary",
            interpolationDataDir=".",
            ):
        """Constructor does only basic argument checking and stores all
        arguments in instance variables.
        """
        OdbToPointsData.__init__(
            self, inputFieldsPerFrame=inputFieldsPerFrame,
            fieldMaximizeAfter=fieldMaximizeAfter,
            relativeToFrame=relativeToFrame)

        if self.relativeFieldNames:
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

        # store arguments as instance variables
        # Note: fieldNames already processed and stored by
        # OdbToPointsData.__init__()

        self.storeArgs(
            locals(),
            "odbPath", "frameList",
            "gridName", "gridData",
            "outputDirName", "fileNamePrefix",
            "projectPrefix", "meshVersion",
            "fieldMaximizeAfter",
            "relativeToFrame", "relativeToReferenceVtkPath",
            "matRegionsOption", "matRegionsNameFilter",
            "getExtraFields",
            "outputFields", "outputFieldsConst",
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
                def stepFrameToName(odbStep,odbFrame):
                    try:
                        frameNames = (
                            self.odb.odbReader.postData.sequence.frameNames)
                    except AttributeError:
                        raise ValueError(
                            "frameNames not defined in associated XXX_postData"
                            " database for odbPath=%s." % self.odbPath)
                    frameName = frameNames[odbStep][odbFrame]
                    msg("frameName from postData: step %s, frame %d,"
                        " frameName %s"
                        % (odbStep, odbFrame, frameName), debugLevel=10)
                    return frameName
            else:
                # stepFrameToName is format string
                stepFrameToNameStr = self.stepFrameToName
                def stepFrameToName(odbStep,odbFrame):
                    return stepFrameToNameStr % (odbStep, odbFrame)
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

            if self.outputFieldsConst:
                # description of time-constant output fields
                output.writerow(["components in time-constant vtk files",])
                output.writerow(["index", "component"])
                for i, field in enumerate(self.outputFieldsConst):
                    if isinstance(field[0], basestring):
                        fieldName = field[0]
                    else:
                        fieldName = field[0].fieldName
                    output.writerow([i+1, fieldName])
                output.writerow([])  # empty line

            # description of output fields
            output.writerow(["components in vtk files for time frames",])
            output.writerow(["index", "component"])
            if self.outputFields:
                for i, field in enumerate(self.outputFields):
                    if isinstance(field[0], basestring):
                        fieldName = field[0]
                    else:
                        fieldName = field[0].fieldName
                    output.writerow([i+1, fieldName])
            else:
                # in case of outputFields not given: take fieldNames
                for i, fieldName in enumerate(self.fieldNames):
                    output.writerow([i+1, fieldName])

            # description of output frames
            try:
                frameNames = self.odb.odbReader.postData.sequence.frameNames
            except AttributeError:
                msg("WARNING: Did not find frame names in XXX_postData. The"
                    " vtk_description.csv file will not contain frame-date"
                    " information.")
            else:
                output.writerow([])
                output.writerow(["odb frame", "date"])
                for odbStep, odbFrame in self.frameList:
                    frameLabel = self.stepFrameToName(odbStep, odbFrame)
                    try:
                        frameName = frameNames[odbStep][odbFrame]
                    except IndexError:
                        msg("ERROR: Could not find step %s frame %d in"
                            " postdata frameNames list. Check postData and"
                            " frameLists!"
                            % (odbStep, odbFrame))
                        continue
                    output.writerow([frameLabel, frameName])

            # close description file
            del output
            msg("Wrote vtk file description file %s"
                % outputFileNameVtkDescription)

        # initialize the output grid and interpolation
        self.outputTopo = getMesh(**self.gridData)
        self.odb.initOutputPoints(
            interpolPickleName=interpolPickleName, **self.gridData)

        # initialize Container for time-const fields
        self.constFieldsContainer = OrderedDict()

        # check if resultsfiles already exist
        if self.skipExisting:
            lastExisting = None
            for i, (odbStep, odbFrame) in enumerate(self.frameList):
                # compose outputFileName
                frameName = self.stepFrameToName(odbStep, odbFrame)
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

    def processConstantFields(self):
        """Finally process constant fields before we come to the frames-loop.
        For odbToVtkGridBox write time-const vtk file.
        """
        if not self.outputFieldsConst:
            # nothing to write...
            msg("No constant fields to write.", debugLevel=10)
            return

        output = FieldsOnStructPointGrid(mesh=self.outputTopo)
        output.updateFields(
            *self.processOutputFields(
                self.constFieldsContainer, self.outputFieldsConst))

        # compose outputFileName
        frameName = "constant"
        outputFileName = self.outputFileNameTemplate % frameName

        # write data to vtk-file
        msg("Writing time constant data to %s" % outputFileName)
        output.toVtk(
            outputFileName,
            description="Abq-Grid-Data time-constant data",
            outputFormat=self.outputFormat)

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
            if self.outputTopo != vtk.mesh:
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

    def initOutputForFrame(self, odbStep, odbFrame):
        """Prepare for new frame data before looping over fields to be read
        from the odb.

        Initialize self.framesFieldsContainer
        """
        # initialize Container for time-varying (per frame) fields
        self.framesFieldsContainer = OrderedDict()

    def storeField(self, field):
        """Store the given field for output. Will be called once for each field
        and multiple times for each frame.
        """
        self.framesFieldsContainer[field.fieldName] = field

    def storeFrame(self, odbStep, odbFrame):
        """Store frame data to output.

        This will be called once for each frame after all fields have been
        processed by self.storeField(). Processes self.outputFields by means
        of L{OdbToPointsData.processOutputFields}().
        """
        # make available time-constant fields to frame data as well
        self.framesFieldsContainer.update(self.constFieldsContainer)

        output = FieldsOnStructPointGrid(mesh=self.outputTopo)
        output.updateFields(
            *self.processOutputFields(self.framesFieldsContainer,
                                      self.outputFields))

        # compose outputFileName
        frameName = self.stepFrameToName(odbStep, odbFrame)
        outputFileName = self.outputFileNameTemplate % frameName

        # write data to vtk-file
        msg("Writing data for step %s, frame %d to %s"
            % (odbStep, odbFrame, outputFileName))
        output.toVtk(
            outputFileName,
            description=("Abq-Grid-Data step %s, frame %d"
                         % (odbStep, odbFrame)),
            outputFormat=self.outputFormat)

#------------------------------------------------------------------------------
def odbToVtkGridBox(
        odbPath, inputFieldsPerFrame, frameList,
        gridName, gridData,
        outputDirName, fileNamePrefix,
        projectPrefix, meshVersion,
        fieldMaximizeAfter=None,
        relativeToFrame=None, relativeToReferenceVtkPath=None,
        matRegionsOption=False, matRegionsNameFilter=None,
        outputFields=None, outputFieldsConst=None,
        getExtraFields=None,
        stepFrameToName=lambda s,f: "F%s%03d" % (s[-1],f),
        skipExisting=False,
        raiseMissingData=True,
        abaqusVersion=None,
        outputFormat="binary",
        interpolationDataDir=".",
        ):
    """Read interpolated field data from the specified odb at points on the
    specified regular grid and write the data to legacy vtk files (that Voxler
    can read).


    Field Input from the ODB
    ========================

    The argument inputFieldsPerFrame identifies the fields to be read
    from the odb for each odb-frame. It's a list with each item being either
    just the field name --like "SDV1", "U_3"-- or a tuple with this field name
    followed by up to three optional "pre-processing" switches.

    The option switches are:
     1. interpolation type: "default" or "const",
     2. maximize flag: True or False (default),
     3. relative flag: True or False (default).
    Default values at the end of this tuple can be omitted.

    Simple field name strings and tuples can be mixed.

    Example:
     >>> odbToVtkGridBox( ...,
     >>>     inputFieldsPerFrame=[
     >>>         # name, interpolation, maximize, relative
     >>>         "S_MIN_PRINCIPAL",
     >>>         ("SDV8, "const"),   # status piecewise const interpolation
     >>>         ("U_3", "default", False, True),  # U_3 relative to BOM
     >>>         "SDV1",
     >>>         ("SDV6", "default", True),], ...)

    B{Option maximize fields:}
    For selected fields take the maximum value of all intermediate time
    frames. Only applicable if there are gaps in the list of requested output
    frames.

    Reading of fields to be maximized starts at the frame following after the
    frame given by the parameter fieldMaximizeAfter. Data for those fields
    will be read for each frame up to the first mentioned frame in frameList.
    The maximum for all those frames will be computed and stored in the output
    file together with all other (ordinary) fields for this first frame in
    frameList.

    Then again only fields that are requested to be maximized will be read
    from the next frame (i.e. frameList[0][1]+1) up to the next frame listed
    in frameList, maximized and combined with the other fields from only the
    frame in frameList.

    All other fields are read for the first frame in frameList and the first
    output will be generated for this first frame in frameList as well.

    B{For example} fieldNames=["LOGP", ("RER", "default", True)],
    fieldMaximizeAfter=("Step-2",4), frameList=[("Step-2",7), ("Step-2",10),
    ("Step-2",13)]: RER will be read for frame 5,6,7 and LogP will be read
    for frame 7 only. The first output file will be for frame 7 including
    LogP for frame 7 and the maximum of RER for frames 5,6,7. Then RER will
    be read for frame 8,9 and 10 and LogP for frame 10 only. Second output
    file for frame 10 contains LogP from frame 10 and for RER the maximum of
    frame 8, 9, 10. Then from frames 11 and 12 only LogP will be read. From
    frame 13 RER and LogP will be read. Next output for frame 13 contains LogP
    (frame 13) and maximum RER from frame 11 to 13.

    B{Option relative fields:}  CONSTRUCTION SITE
    For certain fields (e.g. displacement) values are taken relative to a
    reference frame / reference field. In addition to the corresponding option
    in the inputFieldsPerFrame argument you need either argument
    relativeToFrame or relativeToReferenceVtkPath. These two arguments are
    mutually exclusive.

    B{Option matRegionsOption:} If this optional argument is given and True
    then read the material region number field from the odb. The material code
    field is named "Material".


    Output Fields
    =============

    The argument outputFieldsConst describes the fields to be written to a
    single time-constant output file.

    The argument outputFields describes the fields to be written to per-frame
    output files.

    Those two arguments can optionally further describe the output fields.
    If none of these arguments are given then all input fields from the odb
    will be exported without modification. In this case and if
    matRegionsOption=True then a material code field will additionally be
    written to a time-constant output file.

    Everything described below for the outputFields argument applies to the
    outputFieldsConst argument as well.

    OutputFields is a list, each item describes one field in the desired
    output. Each item can be a simple field name string or a tuple of three
    items. The tuple would contain the field name, a field type ("scalar",
    "vector", "tensor", "str") and a composer.

    For simple string items in the outputFields argument the corresponding
    input field from the odb is passed on to the output without modification.

    If the item in outputFields is a tuple and its last item --composer-- is
    a string then the field from the odb is simply renamed on output. In this
    case the first item in the tuple states the field name in the output file
    and the third item in the tuple states the odb name as in the argument
    inputFieldsPerFrame. The second item is ignored.

    If the item in outputFields is a tuple and its last item --composer-- is
    callable --a function-- then this function will be called for each frame
    in frameList to supply the desired output. The function will be passed
    the following keyword arguments: cls, data, topo, postData, ctrl.
     - cls is the suitable type (class) for the resulting output field
       --derived from L{bae.field_01.Field},
     - data is a dictionary {fieldName: data} containing all input fields from
       the odb (and possibly from the getExtraFields function) for the current
       frame (or all constant fields when supplying outputFieldsConst),
     - topo is the L{MeshStructuredPoints} object describing the output points
       grid,
     - postData is the postData container accompanying the odb, and
     - ctrl is the L{OdbToVtkGridPoints}-object for even more sophisticated
       demands.

    There should always be a trailing dummy keyword-argument container in the
    function declaration to catch additional arguments that might be
    implemented in future versions of this module. This argument ("**_" in the
    example below) also enables you to only explicitly list those arguments
    from the above list that your implementation actually needs.

    Example:
     >>> def composerLogP(cls, data, **_):
     >>>     return cls(log10(x*1000+1) for x in data["SDV3"])
     >>> def composerS1(cls, data, **_):
     >>>     return cls(-x for x in data["S_MIN_PRINCIPAL"])
     >>> def composerCoord(cls, topo, **_):
     >>>     return cls(topo.getPointsIter())
     >>> ...
     >>> odbToVtkGridBox(...
     >>>   inputFieldsPerFrame=["SDV3", "S_MIN_PRINCIPAL"],
     >>>   matRegionsOption=True,
     >>>   outputFields=[
     >>>     ("logP", "scalar", composerLogP),
     >>>     ("RER", "scalar", "SDV6"),
     >>>     ("S1", "scalar", composerS1),
     >>>     ("xyz", "vector", composerCoord),
     >>>     ("Material", "scalar", "Material"),
     >>>     ],
     >>>   ...)

    Note how the "Material" field (i.e. the matcode) is transferred to the
    output without change. And SDV6 is only renamed to RER.

    @param odbPath: Complete path to the odb including .odb extension. Can be
      relative or absolute.

    @param inputFieldsPerFrame: List of field names to be read from the
      odb. CONSTRUCTION SITE:

      Example: ["UR_1","UR_2","UR_3","SDV4","S_MIN_PRINCIPAL","SDV12"]

      Instead of a field name the items of fieldNames can be
      (field name, interpolation type, maximize, relative)- tuples.
      Interpolation type can be "default" or "const" for piecewise constant
      interpolation.
      Example: ["SDV4", ("SDV8", "const"), "SDV12"]

      For further explanation see section "Field Input from the ODB" above.

    @param fieldMaximizeAfter: a (odbStep, odbFrame)-tuple, odbStep being a
      string. For the first output frame (the first item in the
      frameList-argument) data from frames *after* this --i.e. starting with
      (stepNr, frameNr+1)-- will be considered.
    @param relativeToFrame: A (odbStep, odbFrame)-tuple specifying the reference
      time frame for the relativeFieldNames option. Reference field data will
      be read from the current odb from this time frame. This argument is
      mutually exclusive to relativeToReferenceVtkPath.
    @param relativeToReferenceVtkPath: Complete path to a vtk file that
      contains all fields listed in relativeFieldNames. This argument is
      mutually exclusive to relativeToFrame. Field names and output grid must
      match!
    @param matRegionsOption: Boolean flag: Read MATERIAL REGIONS data as well?
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

    @param getExtraFields: Optional function to supply more fields from sources
      other than odb fields. Will be called with some keyword arguments.

      Currently implemented arguments are ctrl, odbStep, odbFrame, topo, odb.

      ctrl is the L{OdbToVtkGridPoints}-object running the procedure. Call
      ctrl.L{storeField()<OdbToVtkGridBox.storeField>} or
      ctrl.L{storeConstField<OdbToVtkGridBox.storeConstField>} in the supplied
      function to actually store a new field with data from external sources.

      topo is a L{MeshStructuredPoints} object identifying the topology
      for the fields to be supplied.

      odb is the L{OdbReaderInterpolateToPoints
      <bae.odb_03.odbReader.OdbReaderInterpolateToPoints>} objects that is used
      to get all the data from the odb. This can be used to get all sorts of
      data other than field data from the odb. Topological data of the output
      points is stored in odb.interpolation which is a
      L{bae.mesh_01.InterpolMeshToPoints} object. The postData database is
      accessible through odb.odbReader.postData.

      B{Important:} The argument list in the definition of getExtraFields must
      always end with I{**_} because in future versions of this module more
      arguments might be passed to getExtraFields. See L{Converter_CombinedVtk}
      for an example.

    @param outputFields: a list of (fieldName, fieldType, composer)-tuples
      describing the fields to be written to per-frame output files.

      Output will only contain fields listed here. Note in
      particular that setting matRegionsOption=True is not sufficient to get a
      material code field: A "Material" entry in outputFields or
      outputFieldsConst is required as well.

      FieldName will be the name in the vtk file (visible for example in
      Paraview but not in Voxler), fieldType can be "scalar", "vector",
      "tensor" or "str". fieldType can be None or arbitrary if composer is
      a string. composer describes how to generate the result field.
      composer can be a string or a function.

      If composer is a string then it will be taken as field name. The field
      with that name is taken from the odb without modification. fieldType will
      be ignored. This is to simply take fields from the odb as they are and
      possibly rename them.

      If composer is a function then it will be called once for each output
      frame with the following keyword arguments: cls, data, topo, postData,
      ctrl.

      For further explanation see section "Output Fields" above.

    @param outputFieldsConst: a list of (fieldName, fieldType, composer)-tuples
      describing the fields to be written to a single time-constant output
      file.

      Time constant output will only contain fields listed here. Note in
      particular that setting matRegionsOption=True is not sufficient to get a
      material code field: A "Material" entry in outputFieldsConst is required
      as well.

      For further explanation see section "Output Fields" above.

    @param gridData: A dict with keys firstPoint, lastPoint and spacing and
      optionally rotString.
      All three values being lists of three floats. This will be passed on as
      keyword-arguments to the L{bae.mesh_01.MeshStructuredPoints} constructor.

    @param gridName: Used to identify the regular grid, in particular the
      corresponding interpolation data.

    @param frameList: List of output frames as (step number, frame number)
      tuples. Example: [ ("Step-2", 24), ("Step-2", 27) ]

    @param stepFrameToName: Can be a function taking a
      (odbStep, odbFrame)-tuple and delivering the frameName. Or can be the
      string "fromPostData". In this case the frameName attribute from the
      associated _postData database is used.

      Example:
       >>> def stepFrameToName(odbStep, odbFrame):
       >>>     return "F%s%03d" % (odbStep[-1], odbFrame)
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
        odbPath, inputFieldsPerFrame, frameList,
        gridName, gridData,
        outputDirName, fileNamePrefix,
        projectPrefix, meshVersion,
        fieldMaximizeAfter=fieldMaximizeAfter,
        relativeToFrame=relativeToFrame,
        relativeToReferenceVtkPath=relativeToReferenceVtkPath,
        matRegionsOption=matRegionsOption,
        matRegionsNameFilter=matRegionsNameFilter,
        outputFields=outputFields,
        outputFieldsConst=outputFieldsConst,
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
            odbPath, inputFieldsPerFrame, frameList,
            gridName, gridPoints,
            outputDirName, fileNamePrefix,
            projectPrefix, meshVersion,
            fieldMaximizeAfter=None,
            relativeToFrame=(2,0),
            matRegionsOption=False, matRegionsNameFilter=None,
            outputFields=None, outputFieldsConst=None,
            getExtraFields=None,
            csvOutputType="CsvPerFrameRowPtColFld",
            filterEmpty=False,
            stepFrameToName=lambda s,f: "F%s%03d" % (s[-1],f),
            raiseMissingData=True,
            abaqusVersion=None,
            interpolationDataDir=".",
            ):
        """Constructor does only basic argument checking and stores all
        arguments in instance variables
        """
        OdbToPointsData.__init__(
            self, inputFieldsPerFrame=inputFieldsPerFrame,
            fieldMaximizeAfter=fieldMaximizeAfter,
            relativeToFrame=relativeToFrame)

        # store arguments as instance variables
        # Note: fieldNames already processed and stored by
        # OdbToPointsData.__init__()
        self.storeArgs(
            locals(),
            "odbPath", "frameList",
            "gridName", "gridPoints",
            "outputDirName", "fileNamePrefix",
            "projectPrefix", "meshVersion",
            "fieldMaximizeAfter",
            "relativeToFrame",
            "matRegionsOption", "matRegionsNameFilter",
            "getExtraFields",
            "outputFields", "outputFieldsConst",
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
        self.outputTopo = MeshUnstructuredPoints(self.gridPoints)
        del self.gridPoints  # not needed anymore, duplicate of self.outputTopo
        self.odb.initOutputPoints(
            grid=self.outputTopo, interpolPickleName=interpolPickleName)

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

    def initOutputForFrame(self, odbStep, odbFrame):
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

    def storeFrame(self, odbStep, odbFrame):
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
                for field in self.processOutputFields(
                        fieldsContainer, self.outputFields))
            self.fieldsContainer = newFldContainer
            del newFldContainer

        self.output.updateFields(*self.fieldsContainer.itervalues())
        del self.fieldsContainer
        self.output.storeFrame((odbStep, odbFrame))

    def closeOutput(self):
        """Finalize output. Will be called after all frames have been stored
        by means of L{storeFrame}
        """
        self.output.closeCsvOutput()

#------------------------------------------------------------------------------
def odbToCsvPoints(
        odbPath, inputFieldsPerFrame, frameList,
        gridName, gridPoints,
        outputDirName, fileNamePrefix,
        projectPrefix, meshVersion,
        fieldMaximizeAfter=None,
        relativeToFrame=(2,0),
        matRegionsOption=False, matRegionsNameFilter=None,
        outputFields=None, outputFieldsConst=None,
        getExtraFields=None,
        csvOutputType="CsvPerFrameRowPtColFld",
        filterEmpty=False,
        stepFrameToName=lambda s,f: "F%s%03d" % (s[-1],f),
        raiseMissingData=True,
        abaqusVersion=None,
        interpolationDataDir=".",
        ):
    """Read interpolated field data from the specified odb at specified points
    and write the data to csv files.


####### CONSTRUCTION SITE AHEAD  ##########################

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

    @param getExtraFields: Optional generator function to supply more fields
      from sources other than odb fields. To be used in conjuction with the
      outputFields argument. Will be called with some keyword arguments.

      Currently implemented arguments are odbStep, odbFrame, points, odb.
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
      names or in the csv header line from a (odbStep, odbFrame)-tuple.
      Might be a function (odbStep, odbFrame)-tuple --> frameName or a template
      string (like the default "F%d%03d").
      Can also be the string "fromPostData" in this case the frameName attribute
      from the associated _postData.py file is used.

      For restart analyses use something like this (We only have step 2 and 3,
      the latter sometimes referred to as 1):
       >>> def stepFrameToName(odbStep, odbFrame):
       >>>     if odbStep==2:
       >>>         return "F2%03d" % odbFrame
       >>>     else:
       >>>         return "F3%03d" % odbFrame
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
        odbPath, inputFieldsPerFrame, frameList,
        gridName, gridPoints,
        outputDirName, fileNamePrefix,
        projectPrefix, meshVersion,
        fieldMaximizeAfter=fieldMaximizeAfter,
        relativeToFrame=relativeToFrame,
        matRegionsOption=matRegionsOption,
        matRegionsNameFilter=matRegionsNameFilter,
        outputFields=outputFields,
        outputFieldsConst=outputFieldsConst,
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
            odbPath, inputFieldsPerFrame, frameList,
            boundingBox=None, elems='ALL_COHS',
            outputDirName=".", fileNamePrefix="PRJXXXX_RXX_COHS",
            fieldMaximizeAfter=None,
            relativeToFrame=(2,0),
            matRegionsOption=False, matRegionsNameFilter=None,
            outputFields=None, outputFieldsConst=None,
            getExtraFields=None,
            csvOutputType="CsvPerFrameRowPtColFld",
            filterEmpty=False,
            stepFrameToName=lambda s,f: "F%s%03d" % (s[-1],f),
            raiseMissingData=True,
            abaqusVersion=None,
            ):
        """Constructor does only basic argument checking and stores all
        arguments in instance variables
        """
        OdbToPointsData.__init__(
            self, inputFieldsPerFrame=inputFieldsPerFrame,
            fieldMaximizeAfter=fieldMaximizeAfter,
            relativeToFrame=relativeToFrame)

        # store arguments as instance variables
        # Note: fieldNames already processed and stored by
        # OdbToPointsData.__init__()
        self.storeArgs(
            locals(),
            "odbPath", "frameList",
            "boundingBox", "elems",
            "outputDirName", "fileNamePrefix",
            "fieldMaximizeAfter",
            "relativeToFrame",
            "matRegionsOption", "matRegionsNameFilter",
            "getExtraFields",
            "outputFields", "outputFieldsConst",
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
        model = self.odb.odbReader.getAbqModel(
            recognizedBlocks=recognizedBlocks)
        msg("Read model: %s." % model)
        self.odb.initOutputAtCentroids(
            model=model, elems=self.elems, boundingBox=self.boundingBox)

        # store as standard-attribute outputTopo
        self.outputTopo = self.odb.interpolation.grid

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
        odbPath, inputFieldsPerFrame, frameList,
        boundingBox=None, elems='ALL_COHS',
        outputDirName=".", fileNamePrefix="PRJXXXX_RXX_COHS",
        fieldMaximizeAfter=None,
        relativeToFrame=(2,0),
        matRegionsOption=False, matRegionsNameFilter=None,
        outputFields=None, outputFieldsConst=None,
        getExtraFields=None,
        csvOutputType="CsvPerFrameRowPtColFld",
        filterEmpty=False,
        stepFrameToName=lambda s,f: "F%s%03d" % (s[-1],f),
        raiseMissingData=True,
        abaqusVersion=None,
        projectPrefix="NotUsed", meshVersion="NotUsed",
        ):
    """Read field data from the specified odb at element centroids
    and write the data to csv files.


####### CONSTRUCTION SITE AHEAD  ##########################

    Usage:
     >>> odbPath = config.odbPath
     >>> outputDirName = "CSV/Cadia2016033_R01_COHSV"
     >>> fileNamePrefix = "Cadia2016033_R01_COHSV"
     >>> fieldNames = ["PST", "RER"]
     >>> frameList = [(2, i) for i in range(10, 20)]
     >>>
     >>> odbElCentroidsToCsvPoints(
     >>>     odbPath,
     >>>     inputFieldsPerFrame=inputFieldsPerFrame, frameList=frameList,
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

    @param getExtraFields: Optional generator function to supply more fields
      from sources other than odb fields. To be used in conjuction with the
      outputFields argument. Will be called with some keyword arguments.

      Currently implemented arguments are odbStep, odbFrame, points, odb.
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
      names or in the csv header line from a (odbStep, odbFrame)-tuple.
      Might be a function (odbStep, odbFrame)-tuple --> frameName or a template
      string. Default if not given, i.e. None: "F%d%03d".

      For restart analyses use something like this (We only have step 2 and 3,
      the latter sometimes referred to as 1):
       >>> def stepFrameToName(odbStep, odbFrame):
       >>>     if odbStep==2:
       >>>         return "F2%03d" % odbFrame
       >>>     else:
       >>>         return "F3%03d" % odbFrame
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
        odbPath, inputFieldsPerFrame, frameList,
        boundingBox=boundingBox, elems=elems,
        outputDirName=outputDirName, fileNamePrefix=fileNamePrefix,
        fieldMaximizeAfter=fieldMaximizeAfter,
        relativeToFrame=relativeToFrame,
        matRegionsOption=matRegionsOption,
        matRegionsNameFilter=matRegionsNameFilter,
        outputFields=outputFields,
        outputFieldsConst=outputFieldsConst,
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
            odbPath, inputFieldsPerFrame, frameList, surface,
            outputDirName=".", fileNamePrefix="PRJXXXX_RXX_COHS",
            fieldMaximizeAfter=None,
            relativeToFrame=(2,0),
            matRegionsOption=False, matRegionsNameFilter=None,
            outputFields=None, outputFieldsConst=None,
            getExtraFields=None,
            csvOutputType="CsvPerFrameRowPtColFld",
            filterEmpty=False,
            stepFrameToName=lambda s,f: "F%s%03d" % (s[-1],f),
            raiseMissingData=True,
            abaqusVersion=None,
            ):
        """Constructor does only basic argument checking and stores all
        arguments in instance variables
        """
        OdbToPointsData.__init__(
            self, inputFieldsPerFrame=inputFieldsPerFrame,
            fieldMaximizeAfter=fieldMaximizeAfter,
            relativeToFrame=relativeToFrame)

        # convert single surface name to list with this name as single item
        if isinstance(surface, basestring):
            surface = [surface,]

        # store arguments as instance variables
        # Note: fieldNames already processed and stored by
        # OdbToPointsData.__init__()
        self.storeArgs(
            locals(),
            "odbPath", "frameList", "surface",
            "outputDirName", "fileNamePrefix",
            "fieldMaximizeAfter",
            "relativeToFrame",
            "matRegionsOption", "matRegionsNameFilter",
            "getExtraFields",
            "outputFields", "outputFieldsConst",
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
        model = self.odb.odbReader.getAbqModel(
            recognizedBlocks=recognizedBlocks)
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

        # store as standard-attribute outputTopo
        self.outputTopo = self.odb.interpolation.grid

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
        odbPath, inputFieldsPerFrame, frameList, surface,
        outputDirName=".", fileNamePrefix="PRJXXXX_RXX_COHS",
        fieldMaximizeAfter=None,
        relativeToFrame=(2,0),
        matRegionsOption=False, matRegionsNameFilter=None,
        outputFields=None, outputFieldsConst=None,
        getExtraFields=None,
        csvOutputType="CsvPerFrameRowPtColFld",
        filterEmpty=False,
        stepFrameToName=lambda s,f: "F%s%03d" % (s[-1],f),
        raiseMissingData=True,
        abaqusVersion=None,
        ):
    """Read field data from the specified odb at element face centroids
    of a specified surface and write the data to csv files.


####### CONSTRUCTION SITE AHEAD  ##########################

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
      option. A (odbStep, odbFrame)-tuple.
    @param matRegionsOption: Boolean flag: Write MATERIAL REGIONS data as well?
      If True then a field named "Material" will be generated.
    @param matRegionsNameFilter: This is for material name conversion. For
      example if there are region variants that shall be treated as one. I.e.
      if SEDIMENT, SEDIMENT_CAVE and SEDIMENT_UCUT shall all be treated as
      SEDIMENT. For further details see the corresponding parameter of
      L{odbToVtkGridBox}.

    @param getExtraFields: Optional generator function to supply more fields
      from sources other than odb fields. To be used in conjuction with the
      outputFields argument. Will be called with some keyword arguments.

      Currently implemented arguments are odbStep, odbFrame, points, odb.
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
      names or in the csv header line from a (odbStep, odbFrame)-tuple.
      Might be a function (odbStep, odbFrame)-tuple --> frameName or a template
      string. Default if not given, i.e. None: "F%d%03d".

      For restart analyses use something like this (We only have step 2 and 3,
      the latter sometimes referred to as 1):
       >>> def stepFrameToName(odbStep, odbFrame):
       >>>     if odbStep==2:
       >>>         return "F2%03d" % odbFrame
       >>>     else:
       >>>         return "F3%03d" % odbFrame
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
        odbPath, inputFieldsPerFrame, frameList,
        surface=surface,
        outputDirName=outputDirName, fileNamePrefix=fileNamePrefix,
        fieldMaximizeAfter=fieldMaximizeAfter,
        relativeToFrame=relativeToFrame,
        matRegionsOption=matRegionsOption,
        matRegionsNameFilter=matRegionsNameFilter,
        outputFields=outputFields,
        outputFieldsConst=outputFieldsConst,
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
            odbPath, inputFieldsPerFrame, frameList,
            initialSurface="top", seqSetsNameExcav=None, seqSetsNameFill=None,
            outputDirName=".", fileNamePrefix="PRJXXXX_RXX_topSurf",

            # maximize option not implemented:
            # fieldMaximizeAfter=None,

            relativeToFrame=(2,0),

            matRegionsOption=False, matRegionsNameFilter=None,

            outputFields=None, outputFieldsConst=None,
            getExtraFields=None,

            # not applicable because the points might change:
            # csvOutputType="CsvPerFrameRowPtColFld",
            # filterEmpty=False,

            stepFrameToName=lambda s,f: "F%s%03d" % (s[-1],f),
            raiseMissingData=True,
            abaqusVersion=None,
            ):
        """Constructor does only basic argument checking and stores all
        arguments in instance variables
        """
        OdbToPointsData.__init__(
            self, inputFieldsPerFrame=inputFieldsPerFrame)

        # store arguments as instance variables
        # Note: fieldNames already processed and stored by
        # OdbToPointsData.__init__()
        self.storeArgs(
            locals(),
            "odbPath", "frameList",
            "initialSurface", "seqSetsNameExcav", "seqSetsNameFill",
            "outputDirName", "fileNamePrefix",
            # "fieldMaximizeAfter",
            "relativeToFrame",
            "matRegionsOption", "matRegionsNameFilter",
            "getExtraFields",
            "outputFields", "outputFieldsConst",
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
        moveMatFldFlag = self.matRegionsOption and not self.outputFieldsConst
        OdbToPointsData.initOutput(self)
        if moveMatFldFlag:
            self.outputFieldsConst.remove(("Material", "Material"))
            self.outputFields.insert(0, ("Material", "Material"))

        # output file names from outputDirName and fileNamePrefix
        self.outputFileNameTemplate = os.path.join(
            self.outputDirName, self.fileNamePrefix+"_%s.csv")

        # get the mesh and elsets
        msg("Reading model data from the odb.")
        self.model = self.odb.odbReader.getAbqModel(
            recognizedBlocks="NODE,ELEMENT,ELSET")
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
            try:
                getattr(self.odb.odbReader.postData.sequence,
                        self.seqSetsNameExcav)
            except AttributeError:
                raise KeyError(
                    "Sequence elsets '%s' not defined in associated"
                    " XXX_postData database for odbPath=%s."
                    % (self.seqSetsNameExcav, self.odbPath))
            elsetsExcav = sequenceElsetsIter(
                vars(self.odb.odbReader.postData.sequence),
                frameList=self.frameList,
                seqElsets=self.seqSetsNameExcav, incremental=True)
        else:
            elsetsExcav = None
        if self.seqSetsNameFill:
            try:
                getattr(self.odb.odbReader.postData.sequence,
                        self.seqSetsNameFill)
            except AtributeError:
                raise KeyError(
                    "Sequence elsets '%s' not defined in associated"
                    " XXX_postData database for odbPath=%s."
                    % (self.seqSetsNameFill, self.odbPath))
            elsetsFill = sequenceElsetsIter(
                vars(self.odb.odbReader.postData.sequence),
                frameList=self.frameList,
                seqElsets=self.seqSetsNameFill, incremental=True)
        else:
            elsetsFill = None

        # initialize self.seqSurfIter
        self.seqSurfIter = self.surf.updateFromSequence(
            model=self.surf.mesh,
            elsetsExcav=elsetsExcav, elsetsFill=elsetsFill)

    def initOutputForFrame(self, odbStep, odbFrame):
        """prepare everything for the current time step

        Initialize self.fieldsContainer

        @param odbStep: not used for this variant
        @param odbFrame: not used for this variant
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

        # store as standard-attribute outputTopo
        self.outputTopo = self.odb.interpolation.grid

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

    def storeFrame(self, odbStep, odbFrame):
        """Store frame data to output. This will be called once for each frame
        after all fields have been processed by self.storeField().

        Processes self.outputFields by means of
        L{OdbToPointsData.processOutputFields}().
        """
        # process output fields (possibly apply composer functions)
        # and store in output
        self.output.updateFields(
            *self.processOutputFields(
                self.fieldsContainer, self.outputFields))
        del self.fieldsContainer
        self.output.storeFrame((odbStep, odbFrame))
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
        for odbStep, odbFrame in self.frameList:

            self.initOutputForFrame(odbStep, odbFrame)

            # mat regions preparation
            if self.matRegionsOption:
                self.storeField(self.getMatRegions())

            # other fields from odb
            for fieldName in self.fieldNames:
                interpolationType = self.fieldInterpolType[fieldName]
                iptText = ""
                if interpolationType=="const":
                    iptText = " interpolating piecewise constant"

                msg("Reading field %s%s, step %s, frame %d."
                    % (fieldName, iptText, odbStep, odbFrame))
                try:
                    field = self.odb.getFieldFromOdb(
                        odbStep, odbFrame, fieldName,
                        interpolationType=interpolationType)
                except ValueError as exc:
                    if self.raiseMissingData:
                        raise
                    else:
                        msg("ERROR: Field %s for odb frame (%s, %d) not"
                            " available. Ending this loop now.\nError: %s"
                            % (fieldName, odbStep, odbFrame, exc))
                        stopNow = True
                        retCode = 1
                        break

                if relativeFieldNames and fieldName in relativeFieldNames:
                    msg("Reading reference field %s%s, from step %s, frame %d."
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
                self.getExtraFields(
                    ctrl=self, odbStep=odbStep, odbFrame=odbFrame,
                    topo=self.outputTopo, odb=self.odb)

            # self.storeFrame evaluates self.outputFields and writes files
            self.storeFrame(odbStep, odbFrame)

        # finally write material names to csv file
        if self.matRegionsOption:
            self.writeMatNamesListToFile()

        return retCode

#------------------------------------------------------------------------------
def odbVariableFaceCentroidsToCsvPoints(
        odbPath, inputFieldsPerFrame, frameList,
        initialSurface="top", seqSetsNameExcav=None, seqSetsNameFill=None,
        outputDirName=".", fileNamePrefix="PRJXXXX_RXX_topSurf",

        # maximize option not implemented:
        # fieldMaximizeAfter=None,

        relativeToFrame=(2,0),

        matRegionsOption=False, matRegionsNameFilter=None,

        outputFields=None, outputFieldsConst=None,
        getExtraFields=None,

        stepFrameToName=lambda s,f: "F%s%03d" % (s[-1],f),
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
     >>> inpFieldNames = ["U_3", "logP"]
     >>> frameList = [("Step-2", i) for i in range(10, 20)]
     >>>
     >>> odbVariableFaceCentroidsToCsvPoints(
     >>>     odbPath,
     >>>     inputFieldsPerFrame=inpFieldNames,
     >>>     frameList=frameList,
     >>>     initialSurface="top", seqSetsNameExcav="seqElsetsPit",
     >>>     outputDirName=outputDirName, fileNamePrefix=fileNamePrefix,
     >>>     )

    ...or for an underground excavation:
     >>> odbVariableFaceCentroidsToCsvPoints(
     >>>     odbPath,
     >>>     inputFieldsPerFrame=inpFieldNames,
     >>>     frameList=frameList,
     >>>     outputDirName=outputDirName, fileNamePrefix=fileNamePrefix,
     >>>     initialSurface="none", seqSetsNameExcav="seqElsetsCrusher",
     >>>     )

    B{Option matRegionsOption:} Add material region number field to every output
    vtk file. If you want the material region number only then choose a single
    arbitrary frame for the frameList parameter and make fieldNames an empty
    list and of course matRegionsOption=True.


####### CONSTRUCTION SITE AHEAD  ##########################

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

    @param getExtraFields: Optional generator function to supply more fields
      from sources other than odb fields. To be used in conjuction with the
      outputFields argument. Will be called with some keyword arguments.

      Currently implemented arguments are odbStep, odbFrame, points, odb.
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
      names or in the csv header line from a (odbStep, odbFrame)-tuple.
      Might be a function (odbStep, odbFrame)-tuple --> frameName or a template
      string. Default if not given, i.e. None: "F%d%03d".

      For restart analyses use something like this (We only have step 2 and 3,
      the latter sometimes referred to as 1):
       >>> def stepFrameToName(odbStep, odbFrame):
       >>>     if odbStep==2:
       >>>         return "F2%03d" % odbFrame
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
        odbPath, inputFieldsPerFrame, frameList,
        initialSurface=initialSurface,
        seqSetsNameExcav=seqSetsNameExcav, seqSetsNameFill=seqSetsNameFill,
        outputDirName=outputDirName, fileNamePrefix=fileNamePrefix,
        relativeToFrame=relativeToFrame,
        matRegionsOption=matRegionsOption,
        matRegionsNameFilter=matRegionsNameFilter,
        outputFields=outputFields,
        outputFieldsConst=outputFieldsConst,
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

#             stepFrameToName=lambda s,f: "F%s%03d" % (s[-1],f),
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
#             self.mesh = self.odb.odbReader.getAbqModel()

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

#         stepFrameToName=lambda s,f: "F%s%03d" % (s[-1],f),
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

        displacementFieldName="U", stressFieldName="S",
        relativeToFrame=None,

        csvOutputType="OneCsvRowPtColFldFr",

        outputDirNameGroups=None,
        csvOutputTypeGroups=None,

        filterEmpty=False,
        stepFrameToName=lambda s,f: "F%s%03d" % (s[-1],f),
        abaqusVersion=None,
        interpolationDataDir=".",
        ):
    """Not fully converted to odb_04 and new postData structure, yet!!!


    Read interpolated field data from the specified odb at the end points of
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
      It is expected to return the associated (odbStep, odbFrame)-tuple.

      IMPORTANT NOTE: To maintain compatibility with possible future versions
      of this module the function should take an arbitrary-length tuple of
      positional arguments!
      I.e. use this phrase: C{def relativeToFrame(*args)}.

      Simple example, frame number as last part of the layer name "horz_F02":
       >>> def relativeToFrame(*args):
       >>>     layerName = args[0]
       >>>     odbFrame = int(layerName.rsplit("_F", 1)[-1])
       >>>     return (2, odbFrame)
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
       >>>             (name, (odbStep, odbFrame))
       >>>             for odbStep, frameData in frameNames.iteritems()
       >>>             for odbFrame, name in enumerate(frameData))
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
       >>>             (el, (odbStep, odbFrame))
       >>>             for odbStep, frameData in seqElsets.iteritems()
       >>>             for odbFrame, elsets in enumerate(frameData)
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
      names or in the csv header line from a (odbStep, odbFrame)-tuple.
      Might be a function (odbStep, odbFrame)-tuple --> frameName or a template
      string (like the default "F%d%03d").
      Can also be the string "fromPostData" in this case the frameName attribute
      from the associated _postData.py file is used.

      For restart analyses use something like this (We only have step 2 and 3,
      the latter sometimes referred to as 1):
       >>> def stepFrameToName(odbStep, odbFrame):
       >>>     if odbStep==2:
       >>>         return "F2%03d" % odbFrame
       >>>     else:
       >>>         return "F3%03d" % odbFrame
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
    # a (odbStep, odbFrame) - tuple.
    # On the other hand all bae.utils.odbToPointsData - methods
    # expect for their stepFrameToName argument a function taking two arguments:
    # odbStep, odbFrame
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
        def stepFrameToName(frameId):
            try:
                frameNames = odb.odbReader.postData.sequence.frameNames
            except AttributeError:
                    raise ValueError(
                        "frameNames not defined in associated XXX_postData"
                        " database for odbPath=%s." % self.odbPath)
            frameName = frameNames[frameId[0]][frameId[1]]
            msg("frameName from postData: step %s, frame %d, frameName %s"
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
        frameNames = odb.odbReader.postData.sequence.frameNames
        dateToRefFrame = dict(
            (date, (odbStep, odbFrame))
            for odbStep, frameData in frameNames.iteritems()
            for odbFrame, date in enumerate(frameData))
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
                    " postData.sequence.frameNames for the name '%s' from the"
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
                    " postData.sequence.frameNames for the name '%s' from the"
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
        nameToFrame = dict((name, relativeToFrame(name, odb.odbReader.postData))
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
        for odbStep, odbFrame in allRefFrames:
            msg("Reading field %s from the odb, step %s, frame %d."
                % (displacementFieldName, odbStep, odbFrame))
            frameToFieldU[(odbStep, odbFrame)] = odb.getFieldFromOdb(
                odbStep, odbFrame, displacementFieldName)
        del allRefFrames

        # reference displacements potentially from different frames
        refFieldU = [
            frameToFieldU[refFr][i]
            for i, refFr in enumerate(refFrameField)]

    elif withClosure and relativeToFrame:
        # relative displacements relative to the same frame for all points
        # creates refFieldU
        odbStep, odbFrame = relativeToFrame
        refFieldU = odb.getFieldFromOdb(
            odbStep, odbFrame, displacementFieldName)
        frameToFieldU[relativeToFrame] = refFieldU
    else:
        # no relative displacements
        # creates dummy refFieldU
        refFieldU = None

    # get the field data from the odb frame by frame
    msg("Now starting to iterate over output frames.")
    for odbStep, odbFrame in frameList:

        # precompute closure
        if withClosure:
            try:
                fieldU = frameToFieldU[(odbStep, odbFrame)]
            except KeyError:
                msg("Reading field %s from the odb, step %s, frame %d."
                    % (displacementFieldName, odbStep, odbFrame))
                fieldU = odb.getFieldFromOdb(
                    odbStep, odbFrame, displacementFieldName)
            else:
                msg("Skip reading field %s from the odb, step %s, frame %d."
                    " Already got it as reference data."
                    % (displacementFieldName, odbStep, odbFrame))

            # calculate relative displacements
            if relativeToFrame:
                if callable(relativeToFrame):
                    fieldU = type(fieldU)(
                        (((odbStep, odbFrame)>refFr) and vector(u0,u1)
                         or [0.0, 0.0, 0.0])
                        for u0, u1, refFr in izip(
                                refFieldU, fieldU, refFrameField))
                elif (odbStep, odbFrame) > relativeToFrame:
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
            msg("Reading field %s from the odb, step %s, frame %d."
                % (stressFieldName, odbStep, odbFrame))
            fieldS = odb.getFieldFromOdb(
                odbStep, odbFrame, stressFieldName)
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
        outputAll.storeFrame((odbStep, odbFrame))
        outputCum.storeFrame((odbStep, odbFrame))

    # finalize output
    outputAll.closeCsvOutput()
    outputCum.closeCsvOutput()
