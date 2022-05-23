r"""Module to handle odbs

first of all: read data easily

All scripts using this module must be run with "abaqus python myscript.py" or
"abaqus cae noGUI=myscript.py".
"""

__version__ = "2.24"
_version_history_ = r"""\
Versions:

1.01 new     - versions not documented so far
1.02 GP added: seperate getOdbFieldClass function, OdbReader now accepts an
     already open odb as input
1.03 GP added: OdbReader.getOdbPartInstance()
1.04 GP fixed: getFrameFieldIterator(frameList=range(1000), ...) looked for
     frame nb 999 for the field information
1.05 GP added: OdbReader.getFrameIterator2()
1.06 GP added: OdbReader.getFrame2()
1.07 GP modif: incompatible change to OdbReader.getFrame2() and
     OdbReader.getFrameIterator2(): stepNo==3 now refers to Step-3 regardless
     whether it's the third step in the odb or not.
1.08 GP modif: removed some apparently not needed imports, now use bae.log_01
        added: OdbReaderInterpolateToGrid
1.09 GP added: OdbReaderInterpolateToGrid.getMatFieldFromOdb()
1.10 GP added: getFrameIterator3()
1.11 GP added: sensible error message if position="elemSP" but field not
     defined at section points.
1.12 GP added: accepts MeshUnstructuredPoints as grid as well
1.13 AF added: added optional kwarg meshinfopickleName for initOutputGrid() in
               class OdbReaderInterpolateToGrid
1.14 GP added: new name and postData properties for the OdbReader object
1.15 GP changed: renamed _OdbField to OdbField, only to make everything more
     transparent
1.16 GP added: OdbReader.getAbqModel() now reads elsets
2.17 GP changed: OdbReaderInterpolateToGrid.initOutputGrid() now uses
     mesh_01.InterpolMeshToPoints and writes a
     different pickle file called "%s_.interpol"%boxName also
     containing essential grid properties.
     changed version number to 2.xxx because it's odb_02
     added OdbReader.getFieldFromOdb() and OdbReader.getMatFieldFromOdb()
2.18 GP updated getOdbFieldClass() to correctly process named SDV fields like
     "SDV_PST" and so on.
2.19 GP unused material types are now being filtered out by
     OdbReaderInterpolateToGrid.getMatFieldFromOdb(...)
2.20 GP OdbReaderInterpolateToGrid.getMatFieldFromOdb(...) now only considers
     elements (and material codes) of section type "solid"; filter for
     unused material types is now optional (and off by default)
2.21 GP adapted to incompatible change in mesh_01: renamed getMeshCallback
     to getMeshCallBack
2.22 GP incompatible change: stepNb arguments now refers to the sequential
     number of steps in the odb. Earlier it was used to construct a name
     "Step-xx" and use this as key in the odb.steps-repository. I.e. in a
     restart analysis step number now start with 1 whereas before they started
     for example with 3 if the base analysis contained 2 steps.
2.23 GP replace call to MeshStructuredPoints with mesh_01.getMesh so that
     rotated and aligned grids can be handled transparently.
2.24 GP added getOdbFieldClassFromOdbField
"""

try:
    import odbAccess
except ImportError:
    raise ImportError(
        "No module named odbAccess. Please run with <abaqus python>.")

from warnings import warn
import re

from bae.future_01 import defaultdict
from bae.misc_01 import Container
from bae.log_01 import log, msg, MsgTicker
import bae.abq_model_02 as abq_model
from bae.field_01 import Field
from bae.mesh_01 import getMesh, MeshUnstructuredPoints, \
    InterpolMeshToPoints
# , \
#     FieldPosMeshNode, FieldPosMeshElement, FieldPosMeshElemIP, \
#     FieldTypeScalar, FieldTypeVector, FieldTypeTensor

_debug = False
# _debug = True

class OdbReadWarn(UserWarning):
    pass

class OdbFieldNotFound(ValueError):
    pass


#---------------------------------------------------------------------
#--- read data from the odb
#... wrapper class for reading from the odb 

class OdbReader(object):
    r"""Class to read data from the odb

    You may access the Abaqus odb object as self.odb (as long as the
    self.close() function has not been called).

    @ivar name: either a copy of the odb argument to self.__init__() if it's
       been a string. It may be retrieved from Abq/cae's odb repository.
       Otherwise you might specify it manually to activate funcionality
       like the .postData property.
    """

    def __init__(self, odb, logfile="default_logfile"):
        """Constructor

        @param odb: file name of the odb or the already opened odb object.

        @param logfile: DEPRECATED! Use bae.log_01 instead!
        Already opened log file for diagnostic output, specify
        logfile=None to suppress diagnostic output, default is stdout
        """

        if isinstance(odb, basestring):
            self.odb = odbAccess.openOdb(path=odb, readOnly=True)
            self.name = odb
        elif str(type(odb))=="<type 'Odb'>":
            self.odb = odb
            self.name = "unnamed_odb"  # dummy in case we don't know

            # Abq/cae has a session object containing an odbs-dict
            import sys
            if sys.argv[0].endswith("ABQvwrK.exe"):
                from abaqus import session
                for name, obj in session.odbs.items():
                    if obj is odb:
                        self.name = name
                        break
        else:
            raise ValueError(
                "odb argument must either be a string specifying the file name"
                " or an open odb. Instead it's of type %s." % type(odb))

        self._postData = None

        if logfile is None:
            log.quiet()
        elif logfile != "default_logfile":
            log.toFile(logfile)

        # initialize default start frame
        self.setFrameZero(0,0)

    @property
    def postData(self):
        """A dictionary being fed from the corresponding ..._postData.py
        database.

        Usually contains some seqElsets... dictionaries like seqElsetsExcav in
        the exampe below. And a frameNames dictionary.

        Example: Suppose in frame 3 of step 2 labelled "YR2012_M02" the elsets
        "Pit_2012_02" and "Dev_2012_02" will be excavated and in frame 4
        "YR2012_M03" only "Dev_2012_03":
         >>> odb = OdbReader("/datapth/myodb.odb")
         >>> frameNames = odb.postData["frameNames"]
         >>> print frameNames[2][3:5]
         ... ['YR2012_M02','YR2012_M03']
         >>> seqSets = odb.postData["seqElsetsExcav"]
         >>> print seqSets[2][3:5]
         ... [ ('Pit_2012_02','Dev_2012_02'),(Dev_2012_03,) ]

        The data will only be read when the portData property is being accessed
        for the first time.

        @raises IOError: If the XXX_postData.py file can not be accessed. XXX is
          name of the odb stripped of a trailing ".odb"
        """
        if self._postData is None:
            # deduct postData-filename from odb name
            if self.name[-4:].lower() == ".odb":
                postDataFileName = self.name[:-4]+"_postData.py"
            else:
                postDataFileName = self.name+"_postData.py"
            # load/parse postData from external file
            self._postData = dict()
            execfile(postDataFileName, self._postData)

        return self._postData

    def close(self):
        self.odb.close()

    def setFrameZero(self, step, frame):
        r"""Specify the odb step and frame that will be treated as frame with
        index zero in self.frameData. frameData frames are consecutively
        numbered regardless of how many steps they span.

        Only consider the last step of the analysis:
         >>> self.setFrameZero(-1,0)

        step specifies the odb step and might be a string or a zero base index
        (-1 == last)

        frame is the zero based odb frame number (-1 == last)

        Creates self.frameNoToStepFrame a list of tuples of odb step names and
        frame numbers starting with the specified frame. This list can be used
        to find the odb step and frame for a certain consecutive frame number
        used in self.frameData.

        [(odb step name, odb frame number), (... for second frame ...), ...]

        See also: getFrameIterator(), getNumberOfFrames()
        """
        allStepNames = self.odb.steps.keys()

        if type(step)==int:
            firstStepNo = step
            while firstStepNo<0:
                firstStepNo += len(allStepNames)
            if firstStepNo>=len(allStepNames):
                raise ValueError("There is no %dth step in the odb."
                                 % firstStepNo)
        elif type(step)==str:
            try:
                firstStepNo = allStepNames.index(step)
            except ValueError:
                raise ValueError("Step %s not in the odb."%step)
        else:
            raise ValueError("Don't know, what to do with step %s of type %s."
                             % (str(step), str(type(step))))
        if type(frame)!=int:
            raise ValueError("Don't know, what to do with frame %s of type %s."
                             % (str(frame), str(type(frame))))
        firstFrameNo = frame

        self.frameNoToStepFrame = list()
        for stepNo in range(firstStepNo, len(allStepNames)):
            nb_frames = len(self.odb.steps[allStepNames[stepNo]].frames)

            # self.nbFrame=nb_frames
            if firstFrameNo!=None:
                # first step: start at firstFrameNo
                frameNo = firstFrameNo
                while frameNo<0:
                    frameNo += nb_frames
                firstFrameNo = None
            else:
                # one of the subsequent steps: start at second frame
                frameNo = 1

            for frameNo in range(frameNo, nb_frames):
                self.frameNoToStepFrame.append(
                    (allStepNames[stepNo], frameNo))


    def getFrame(self, frameId):
        r"""Get the OdbFrame object associated with the given frame number
        based on the setting by setFrameZero().

        @param frameId: frame index. Index zero identifies the frame specified
          by self.setFrameZero(). Negative values can be used to count from the
          end like usual list indices.
        """
        stepName, frameNo = self.frameNoToStepFrame[frameId]
        return self.odb.steps[stepName].frames[frameNo]


    def getFrameIterator(self, frameIds=None):
        r"""Iterator over all frames, yields a tuple (frame number, frame
        object) in each iteration. The frame number is based on the setting by
        setFrameZero() and the frame object is an OdbFrame object.

        Usage:
         >>> for frameNo, frame in self.getFrameIterator():
         ...     print ("picture %d shows odb frame %d, simulation time %g"
         ...            % (frameNo, frame.incrementNumber, frame.frameValue))

        using this iterator might be faster than the following alternative:
         >>> for stepName, frameNo in self.frameNoToStepFrame:
         ...     frame = self.odb.steps[stepName].frames[frameNo]

        @param frameIds: list of frame indexes. Index zero identifies the frame
          specified by self.setFrameZero(). If None or not specified: yield all
          frames up to the end. Negative values can be used to count from the
          end like usual list indices. They are converted to their corresponding
          positive value on output. Frame indexes that do not exist in the odb,
          more precisely in self.frameNoToStepFrame, are ignored but a warning
          is written to the log file. If frameIds[-1]>=900 then it is assumed
          that all frames to the end were requested and the warning mentioned
          above is suppressed.
        """
        frameIds = self.getFrameNumberList(frameIds)
        oldstepName = None
        for frameId in frameIds:
            stepName, frameNo = self.frameNoToStepFrame[frameId]
            if stepName != oldstepName:
                currentStep = self.odb.steps[stepName]
                oldstepName = stepName
            yield frameId, currentStep.frames[frameNo]


    def getFrame2(self, stepNo, frameNo):
        r"""Get the OdbFrame object associated with the given step and frame
        numbers.

        @param stepNo: step number. The first step in the odb has got stepNo=1
        @param frameNo: frame number in the step. Negative values can be used
          to count from the end like usual list indices.
        """
        stepName = self.odb.steps.keys()[stepNo-1]
        return self.odb.steps[stepName].frames[frameNo]


    def getFrameIterator2(self, framesList=None):
        r"""Iterator over frames identified by step nb and frame nb, yields a
        (step number, frame number, OdbFrame object) - tuple in each iteration.

        Usage:
         >>> framesList = [
         ...    (2, [10, 11, 12]),
         ...    (3, [1, 2, 3, 4, 5]),
         ...    ]
         >>> for stepNb, frameNb, frame in self.getFrameIterator2(framesList):
         ...     print ("odb frame %d, simulation time %g"
         ...            % (frame.incrementNumber, frame.frameValue))

        @param framesList: list (or other iterable) of (step number, frame
          number list) - tuples. If None or not specified: yield all frames up
          to the end. 

          The step number: In the ordinary two steps analysis the possible
          values are 1 and 2. Note that if the first step in the odb is called
          step-4 the corresponding step number is 1.

          The frame number list contains frame numbers as in the odb. Instead
          of a frame number list there can be None which means: yield all
          frames for this step. Negative values can be used to count from the
          end like usual list indices. Frame numbers that do not exist in the
          odb are ignored without any warning.

        @Returns: a (stepNo, frameNb, odbFrame) tuple on each iteration.
        """

        if framesList is None:
            framesList = [[s, None] for s in self.odb.steps.keys()]

        for stepNo, frameIdList in framesList:
            stepName = self.odb.steps.keys()[stepNo-1]
            odbStep = self.odb.steps[stepName]
            if frameIdList is None:
                frameIdList = range(len(odbStep.frames))
            else:
                nbFrs = len(odbStep.frames)
                try:
                    frameIdList = [i for i in frameIdList if -nbFrs<=i<nbFrs]
                except TypeError:
                    raise ValueError(
                        "framesList must be a list of (step number, frame"
                        " number list) - tuples. For step number %d we got"
                        " %s as frame number list which is not iterable."
                        " You may want to look at the documentation or"
                        " choose from getFrameIterator(), getFrameIterator2()"
                        " or getFrameIterator3()."
                        % (stepNo, frameIdList))
            for frameNb in frameIdList:
                odbFrame = odbStep.frames[frameNb]
                yield (stepNo, frameNb, odbFrame)


    def getFrameIterator3(self, framesList):
        r"""Iterator over frames identified by a list of (step nb, frame nb)
        tuples. Yields a (step number, frame number, OdbFrame object) - tuple
        in each iteration.

        Usage:
         >>> framesList = [
         ...    (2,10), (2,11), (2,12),
         ...    (3,1), (3,2), (3,3), (3,4), (3,5),
         ...    ]
         >>> for stepNb, frameNb, frame in self.getFrameIterator3(framesList):
         ...     print ("odb frame %d, simulation time %g"
         ...            % (frame.incrementNumber, frame.frameValue))

        @param framesList: list (or other iterable) of (step number, frame
          number) - tuples.

          The step number is corresponding to the sequence of steps in the odb
          starting with 1. In the ordinary two steps analysis the possible
          values are 1 and 2. Note that if the first step in the odb is called
          step-4 the corresponding step number is 1.

          The frame number is as in the odb. Negative values can be used to
          count from the end like usual list indices. Frame numbers that do not
          exist in the odb are ignored without any warning.

        @Returns: a (stepNo, frameNb, odbFrame) tuple on each iteration.
        """

        oldStepNo = None
        for stepNo, frameNo in framesList:

            # find odb step
            if stepNo != oldStepNo:
                try:
                    stepName = self.odb.steps.keys()[stepNo-1]
                    odbStep = self.odb.steps[stepName]
                except IndexError:
                    continue
                oldStepNo = stepNo
                nbFrs = len(odbStep.frames)  # last frame: frameNo = nbFrs-1

            if frameNo>=0 and frameNo<nbFrs:
                odbFrame = odbStep.frames[frameNo]
                yield (stepNo, frameNo, odbFrame)


    def getNumberOfFrames(self):
        r"""Returns the number of frames in the odb starting at the one
        specified by self.setFrameZero(). This is the number of frames
        self.getFrameIterator() yields WHEN CALLED WITHOUT ARGUMENTS.
        """
        return len(self.frameNoToStepFrame)


    def getFrameNumberList(self, frameIds=None):
        r"""Return a list of frame indexes that are actually in the odb. This
        method just filters out frameIds to which the corresponding frame
        actually does not exist in the odb.

        @param frameIds: list of frame indexes. Index zero identifies the frame
          specified by self.setFrameZero(). If None or not specified: yield all
          frames up to the end. Negative values can be used to count from the
          end like usual list indices. They are converted to their corresponding
          positive value on output. Frame indexes that do not exist in the odb,
          more precisely in self.frameNoToStepFrame, are ignored but a warning
          is written to the log file. If frameIds[-1]>=900 then it is assumed
          that all frames to the end were requested and the warning mentioned
          above is suppressed.
        """
        if frameIds is None:
            # just return all possible frame numbers
            frameNumberList = range(len(self.frameNoToStepFrame))
        else:
            # convert neg. numbers and check if valid
            nbFrames = len(self.frameNoToStepFrame)
            frameNumberList = list()
            for frameId in frameIds:
                # convert negative values, needed for the output
                while frameId<0:
                    frameId+=len(self.frameNoToStepFrame)
                # look up the corresponding step, frame in the odb
                if frameId < nbFrames:
                    frameNumberList.append(frameId)
                else:
                    if frameIds[-1]<900:
                        msg(
                            "WARNING: OdbReader.getFrameNumberList() got an"
                            " invalid frame id %d in its argument frameIds."
                            " Must be below %d."
                            % (frameId, len(self.frameNoToStepFrame)))
                    continue
        return frameNumberList

    def getOdbFieldClass(self, fieldName, odbFrame=None):
        """
        Create a new class based on (being child of) L{bae.field_01.Field} and
        L{OdbField} with all the info how to get the data from the odb.

        Returns a class, not an object! You can create multiple objects of this,
        one for each frame for example.  You will then typically call their
        L{readFromFrame<OdbField.readFromFrame>} or
        L{readFromOdbField<OdbField.readFromOdbField>} methods to load actual
        data from the odb.

        The new class has the following class variables
         - fieldName: The fieldName argument.
         - odbFieldDataAttr: The name of the data attribute of the
        FieldValue-object that holds the actual data. It's "data" or
        "dataDouble".
         - odbFieldName: "S_MISES" -> "S", "U_1" -> "U", "SDV4" -> "SDV4"
         - odbFieldComponent: "S_MISES"->None, "SDV4"->None, "U_1"->0,
        "U_2"->1
         - odbFieldInvariant: "S_MISES"->"mises",
        "S_MIN_PRINCIPAL"->"minPrincipal", "SDV4"->None, "U_1"->None

        Usage:
         >>> odb = OdbReader("mydata.odb")
         >>> FieldU = odb.getOdbFieldClass("U")
         >>> frameFieldU = dict()
         >>> for frameNo, frame in odb.getFrameIterator([10, 12, -1]):
         ...     try:
         ...         odbField = frame.fieldOutputs[FieldU.odbFieldName]
         ...     except KeyError:
         ...         print "WARNING: No field U in frame %d." % frameNo
         ...         continue
         ...     newfield = FieldU()
         ...     newfield.readFromOdbField(odbField)
         ...     queryNodeNb = 10
         ...     print ("Displacement at node %d = %s"
         ...            % (queryNodeNb, newfield[queryNodeNb]))
         ...     frameFieldU[frameNo] = newfield

        @param fieldName: examples: "S", "V", "S_MISES", "S_MIN_PRINCIPAL",
        "U_1", "S_33", "SDV4"

        @param odbFrame: OdbFrame object in the odb where the field information
        like scalar/tensor, nodal/element should be looked up. Defaults to the
        last frame in the odb.

        @returns: A new class. It is a subclass of OdbField and two of the
        Field... classes.

        @Raises OdbFieldNotFound: if there is no field in the odb matching
        fieldName.

        @Note: This might not work for vector or tensor components in other
        than 3D solid elements because of the assumption in the
        odbComponentToArrayIndex dictionary! In this case use
        getScalarField(componentLabel="S11") to get the component.
        """

        if odbFrame is None:
            odbFrame = self.getFrame(-1)

        return getOdbFieldClass(fieldName, odbFrame)

    def getFieldIterator(self, frame, fieldClassList, frameName):
        r"""Iterator over the fields in an odb frame.

        @param frame: OdbFrame object
        @param fieldClassList: A list of OdbField classes (not objects!)
        @param frameName: A string identifying the frame for diagnostic output.

        @returns: An iterator that yields (fieldClass, FieldValue) tuples.
          fieldClass is one of the items in fieldClassList and FieldValue is
          the corresponding FieldValue object from the odb.
        """
        odbFields = frame.fieldOutputs

        for fieldClass in fieldClassList:
            # odbField
            try:
                odbField = odbFields[fieldClass.odbFieldName]
            except KeyError:
                if fieldClass.odbFieldName != fieldClass.fieldName:
                    fieldStr = "%s (for %s)" % (fieldClass.odbFieldName,
                                                fieldClass.fieldName)
                else:
                    fieldStr = fieldClass.fieldName
                warn(OdbReadWarn("No field %s in frame %s."
                                 % (fieldStr, frameName)))
                continue
            yield (fieldClass, odbField)


    def getFrameFieldIterator(self, frameIds=None, fieldList=None):
        r"""Iterator over all frames, yields a tuple (frame number, fields)
        in each iteration. The frame number is based on the setting by
        setFrameZero(). fields is an iterator of the kind returned by
        self.getFieldIterator().

        Usage:
         >>> frameVarDict = dict()
         >>> for frameNo, fields in self.getFrameFieldIterator(
         ...     [1,10, 100], ["U", SDV4]):
         ...     print ("data for frame %d:" % frameNo)
         ...     newFrame = dict()
         ...     for FieldClass, odbField in fields:
         ...         newfield = FieldClass().readFromOdbField(odbField)
         ...         newFrame[FieldClass.fieldName] = newfield
         ...     frameVarDict[frameNo] = newFrame
         >>> print ("displacement in frame 10: %s" % (frameVarDict[10]["U"]))

        @param frameIds: list of frame indexes. Index zero identifies the frame
          specified by self.setFrameZero(). If None or not specified: yield all
          frames up to the end. Negative values can be used to count from the
          end like usual list indices. They are converted to their corresponding
          positive value on output.
        @param fieldList: A field name or a list of field names,
          e.g. ["S", "U_3", "SDV4", "U_MAGNITUDE"]
        """
        if isinstance(fieldList, basestring):
            fieldList=[fieldList, ]

        # indicate first iteration in getFrameIterator - loop
        fieldClasses = None

        for frameNo, frame in self.getFrameIterator(frameIds=frameIds):

            if fieldClasses is None:
                fieldClasses = [
                    self.getOdbFieldClass(fieldName, odbFrame=frame)
                    for fieldName in fieldList]

            yield frameNo, self.getFieldIterator(
                frame, fieldClasses, str(frameNo))

    def getFieldFromOdb(self, fieldName, odbFrame):
        """read the data for the given field from the odb.

        @Raises OdbFieldNotFound: if there is no field in the odb matching
        fieldName.
        """

        FieldClass = self.getOdbFieldClass(fieldName, odbFrame)

        msg("Reading values for field %s from the odb."%fieldName,debugLevel=1)
        field = FieldClass()

        # get the right odb field instance (right field, right frame)
        try:
            odbField = field.odbFieldFromFrame(odbFrame)
        except KeyError:
            raise Exception("ERROR: Field %s not found but getFieldClasses()"
                            " succeded. Weird..." % fieldName)

        # read values from odb
        field.readFromOdbField(odbField)
        msg("Read %d values." % len(field), debugLevel=1)
        return field

    def getMatFieldFromOdb(self, nameFilter=None, firstMatCode=1, noMatCode=-1,
                           secTypeFilter=None):
        """Read material section name from the odb.

        @param nameFilter: This is for material name conversion. For example
          if there are region variants that shall be treated as one. I.e. if
          SEDIMENT, SEDIMENT_CAVE and SEDIMENT_UCUT shall all be treated as
          SEDIMENT.

          If this is a function it will be called as a conversion function
          taking the original material name from the odb as only
          argument and returning the material name it should be treated as.

        @param firstMatCode: Material code number of the first material in
          the list of material names. Or if you prefer: Offset of the material
          code numbers in the matNumberField field and the index of the
          corresponding material in the matNamesList.

          I.e. if firstMatCode=1 then all points with assigned material number
          of 1 are of the first material type in the list matNamesList which is
          matNamesList[0]. All points having material number 5 are of material
          type matNamesList[4].

          Specify firstMatCode=0 if you want the material numbers in
          matNumberField to directly specify the index in matNumberField. The
          default value of 1 seems to be more appropriate for the usual use
          case of exporting matNumberField to a vtk file and matNamesList to a
          seperate csv file.

        @param noMatCode: Material number to represent points that are not
          being assinged any material to.
          Currently this does not work because each and any element is
          associated to a section in Abaqus.

        @param secTypeFilter: Iterable (i.e. list, tuple) of section types
          to be recognized. Elements of any other section type will be
          silently ignored. Can contain "solid", "Unknown type", ...(?).
          Tet elements are of type "solid", cohesive elements are of type
          "Unknown type" (reason not known).

          If None or an empty list then all elements will be considered.

          If you need other types like "shell", "beam" then please check if
          those types work as expected. You get a printout of section types
          in debugging information to be switched on by setting:
          from bae.log_01 import log; log.setDebugLevel(10)

        @returns: (matNamesList, matNumberField)-tuple
          matNamesList is the list of material names.
          matNumberField is the field (i.e. a {element label: material number}
          -dict) of material numbers
        """

        # extract material region info
        reSectName = re.compile(r"(\S*)\s*<\s*(.*?)\s*>")

        MatNumberFieldClass = Field.classFromPosType(
            "Material", "element", "scalar")
        matNumberField = MatNumberFieldClass()

        matNameToMatNumber = dict()
        filteredMatNameList = []
        secTypeFilter = set(secTypeFilter)
        filteredSecTypes = set()

        instance = self.getOdbPartInstance()

        tmt = ("... extracting material region info: element %%d/%d\n"
               % len(instance.elements) )
        ticker = MsgTicker(tmt)

        for e in instance.elements:
            ticker.tick()

            # get section name from the odb
            secName = e.sectionCategory.name

            # split the section string into section type and material name
            res = reSectName.search(secName)
            try:
                secType, matName = res.groups()
            except AttributeError:
                secType = "Unknown type"
                matName = secName

            if secTypeFilter:
                if secType in filteredSecTypes:
                    continue
                if secType not in secTypeFilter:
                    msg("Found element(s) of section type %s. They will be"
                        " ignored." % secType, debugLevel=1)
                    filteredSecTypes.add(secType)
                    continue

            # find old or create new material number
            try:
                matNumber = matNameToMatNumber[matName]
            except KeyError:

                msg("New material found: %s, section type: %s"
                    % (matName, secType), debugLevel=1)
                # replace some parts of the material name as specified
                if callable(nameFilter):
                    filteredMatName = nameFilter(matName)
                else:
                    filteredMatName = matName

                # find material number for filtered name
                try:
                    matNumber = matNameToMatNumber[filteredMatName]
                except KeyError:
                    matNumber = len(filteredMatNameList) + firstMatCode
                    filteredMatNameList.append(filteredMatName)
                    matNameToMatNumber[filteredMatName] = matNumber
                    msg("... stored as %s with new number %d."
                        %(filteredMatName,matNumber), debugLevel=1)
                else:
                    msg("... stored as %s."%filteredMatName, debugLevel=1)

                # store the correct matNumber for this region name
                matNameToMatNumber[matName] = matNumber

            # store the correct matNumber for the current element
            matNumberField[e.label] = matNumber

        del ticker
        msg("Identified %d different material codes for %d elements."
            % (len(filteredMatNameList), len(matNumberField)))
        msg("... %d elements have been ignored. (Because of the secTypeFilter"
            " presumably.)" % (len(instance.elements)-len(matNumberField)),
            debugLevel=1)
        assert secTypeFilter or len(matNumberField)==len(instance.elements)

        return (filteredMatNameList, matNumberField)

    @staticmethod
    def writeMatNamesToCsv(outputFileName, matNamesList, firstMatCode=1):
        """Convenience function to write the material number - material regions
        relations to a csv file.
        """
        import csv
        output = csv.writer(open(outputFileName, "wb"))
        matList = [ [cnt+firstMatCode, matName]
                    for cnt, matName in enumerate(matNamesList) ]
        output.writerow(["material number", "material name"])
        output.writerows(matList)

    def getOdbPartInstance(self):
        r"""Returns the first instance in the rootAssembly object of the odb.
        This contains the model definition, i.e. the mesh, material props and
        so on.

         >>> odb = OdbReader("mydata.odb")
         >>> instance = odb.getOdbPartInstance()
         >>> # read node info
         >>> for node in instance.nodes:
         >>>     output.write("%d:%s\n" % (node.label, list(node.coordinates)))
         >>> # read element info
         >>> for elem in instance.elements:
         >>>     output.write("%d:%s\n" % (elem.label, elem.type))

        @Note: There is a getAbqModel() - method to read the model data. It's
        probably more convenient if it serves your needs.
        """
        return self.odb.rootAssembly.instances.values()[0]

    def getAbqModel(self, logfile="UseTheSameLogFile", recognizedBlocks=None):
        r"""get an abq_model_02.Model-object from the odb model data,

        reads nodes and elements and elsets

        assumes there is only one part instance in the model

        ideas for the future:
         - add the possibility to selectively read only elements, nodes (see
           recognizedBlocks argument of abq_model_02.Model.read() / write())
         - read nsets, surfaces, all the rest

        @param recognizedBlocks: If specified only get the specified property
        from the odb. Implemented so far: "NODE", "ELEMENT", "ELSET". Can be
        specified as tuple of string or comma separated string.

        @Note: parameter logfile is being ignored. It's there for compatibility
        reasons only.
        @Note: elsets are being read from the first part instance and from the
        root assembly as well. If the same elset name occurs twice the elset
        in the root assembly is ignored.
        @Note: Look in abq_model_02.Model.read how recognizedBlocks is treated
        there and adapt this here accordingly.
        """
        model = abq_model.Model()

        # mesh data from the first instance
        part_instance = self.getOdbPartInstance()

        # find dimensionality of our problem
        if part_instance.embeddedSpace == odbAccess.TWO_D_PLANAR:
            self.dim = 2
        elif part_instance.embeddedSpace == odbAccess.THREE_D:
            self.dim = 3
        else:
            msg = "Can't handle dimensionality %s."%part_instance.embeddedSpace
            msg(msg)
            raise Exception(msg)

        # extract undeformed node coordinates
        if recognizedBlocks is None or "NODE" in recognizedBlocks:
            msgTemplate = ("Read %%d/%d nodes from the odb."
                           % len(part_instance.nodes))
            ticker = MsgTicker(msgTemplate, printLast=True)
            for node in part_instance.nodes:
                model.nodeCoords[node.label] = list(node.coordinates)
                ticker.tick()
            del ticker

        # extract element/cell connectivity
        if recognizedBlocks is None or "ELEMENT" in recognizedBlocks:
            msgTemplate = ("Read %%d/%d elements from the odb."
                           % len(part_instance.elements))
            ticker = MsgTicker(msgTemplate, printLast=True)
            for elt in part_instance.elements:
                model.elNodes[elt.label] = list(elt.connectivity)
                model.elType[elt.label] = str(elt.type)
                model.typeEl[elt.type].add(elt.label)
                ticker.tick()
            del ticker

        # extract element sets
        if recognizedBlocks is None or "ELSET" in recognizedBlocks:
            elsetDictAndName = [
                (elset, e)
                for elset in [self.odb.rootAssembly.elementSets,
                              part_instance.elementSets]
                for e in elset.keys()]
            ticker = MsgTicker("Read %%d/%d elsets from the odb."
                               % len(elsetDictAndName), printLast=True)
            for elset, elsetName in elsetDictAndName:
                elset2 = elset[elsetName].elements
                if type(elset2).__name__ == 'OdbMeshElementArray':
                    pass
                elif type(elset2).__name__ == 'tuple':
                    elset2 = elset2[0]
                else:
                    raise Exception("Element set of type %s not implemented so"
                                    " far" % type(elset2).__name__)
                model.elset[elsetName] = set( e.label for e in elset2 )
                ticker.tick()
            del ticker

        return model

    def getFields(self, frameIdx=-1):
        """which fields are there in the odb in the specified frame
        (does not read the data itself, use self.readData() for that

        frameIdx is the consecutive frame number used for self.frameData based
        on the frame zero specified with self.setFrameZero().

        yields a dict {fieldname: additional data}, where additional data is a
        Container object with the following attributes:
         - description (a string)
         - dataType ("scalar" / "vector" / "tensor")
         - componentLabels (list of strings)
         - validInvariants (list of strings)
        """
        stepName, frameNo = self.frameNoToStepFrame[frameIdx]
        frame = self.odb.steps[stepName].frames[frameNo]
        fields = frame.fieldOutputs

        fieldDict = dict([
            (str(f.name), Container(
                description = str(f.description),
                dataType = _odbFieldTypeToFieldType[f.type],
                componentLabels = list(f.componentLabels),
                validInvariants = [str(x) for x in f.validInvariants]
            ))
            for f in fields.values()])
        return fieldDict

    def odbSetFromElements(self, elset):
        """
        Create an odb set containing the given elements. This can be used to
        filter the output.

        Usage example needed!

        @param elset: elememt numbers
        @returns: the odb set
        """
        # future: also accept elset names... if isinstance(elset, basestring):

        elset = sorted(elset)

        # mesh data is in the first instance
        part_instance = self.getOdbPartInstance()

        odbItemNumberLastUsed = 1
        while 1:
            setName = "gen_elset_%d" % odbItemNumberLastUsed
            if not(part_instance.elementSets.has_key(setName)):
                break
            odbItemNumberLastUsed += 1

        # create the new odb elset
        return part_instance.ElementSetFromElementLabels(
            name=setName, elementLabels=elset)

    def odbSetFromNodes(self, nset):
        """
        Create an odb set containing the given nodes. This can be used to
        filter the output.

        Usage example needed!

        @param nset: node numbers
        @returns: the odb set
        """
        # future: also accept nset names... if isinstance(nset, basestring):

        nset = list(nset)
        nset.sort()

        # mesh data is in the first instance
        part_instance = self.getOdbPartInstance()

        odbItemNumberLastUsed = 1
        while 1:
            setName = "gen_nset_%d" % odbItemNumberLastUsed
            if not(part_instance.nodeSets.has_key(setName)):
                break
            odbItemNumberLastUsed += 1

        # create the new odb nset
        return part_instance.NodeSetFromNodeLabels(
            name=setName, nodeLabels=nset)


#---------------------------------------------------------------------
#--- odb field class
#... add-on for bae.field_01 to get field info from the odb and then
#... to also read values from the odb

class OdbField(object):
    r"""Additional base class for Field classes that should read their data from
    the odb. Has additional attributes and methods to load the field data from
    the odb.

    @cvar fieldName: identifies the field (The fieldName argument to
       L{getOdbFieldClass}.)
    @cvar odbFieldDataAttr: The name of the data attribute of the
       FieldValue-object that holds the actual data. It's "data" or
       "dataDouble".
    @cvar odbFieldName: Field identifier in the odb. E.g. "S_MISES" -> "S",
       "U_1" -> "U", "SDV4" -> "SDV4"
    @cvar odbFieldComponent: "S_MISES"->None, "SDV4"->None, "U_1"->0,
       "U_2"->1
    @cvar odbFieldInvariant: "S_MISES"->"mises",
       "S_MIN_PRINCIPAL"->"minPrincipal", "SDV4"->None, "U_1"->None

    @Note: You will usually use dynamically created classes (types) based on
    L{bae.field_01.Field} and OdbField.
    Call L{OdbReader.getOdbFieldClass} or L{getOdbFieldClass} to get such a
    class (type).
    Objects (instances) of those types can not easily be pickled using the
    ordinary pickle.dump() method. Use the pickle_dump() and pickle_load()
    methods instead.
    """
    def odbFieldFromFrame(self, odbFrame):
        r"""Return the odb field from the given frame of the odb that
        corresponds to this Field instance.

        Usage:
         >>> MyField = odb.getOdbFieldClass("U")
         >>> filterSet = odb.odbSetFromElements([10624,2627,8912,...])
         >>> for frameNb, odbFrame in odb.getFrameIterator():
         >>>     field = MyField()
         >>>     odbField = field.odbFieldFromFrame(odbFrame)
         >>>     partData = odbField.getSubset(region=filterSet)
         >>>     field.readFromOdbField(partData)

        @raise KeyError: if the corresponding odb field is not contained in the
        specified odb frame.
        """
        return odbFrame.fieldOutputs[self.odbFieldName]

    def readFromFrame(self, odbFrame):
        r"""To read the field values from the given frame of the odb.
        Usage:
         >>> field = dict()
         >>> MyField = odb.getOdbFieldClass("S_MISES")
         >>> for frameNb, odbFrame in odb.getFrameIterator():
         >>>     field[frameNb] = MyField().readFromFrame(odbFrame)

        This is exactly equivalent to:
         >>> ...
         >>> for frameNb, odbFrame in odb.getFrameIterator():
         >>>     odbField = odb.odbFieldFromFrame(odbFrame)
         >>>     field[frameNb] = MyField().readFromOdbField(odbField)
        """
        return self.readFromOdbField(self.odbFieldFromFrame(odbFrame))

    def readFromOdbField(self, odbField):
        r"""To read the field values from the odb.

        Usage:
         >>> MyField = odb.getOdbFieldClass("U")
         >>> filterSet = odb.odbSetFromElements([10624,2627,8912,...])
         >>> for frameNb, odbFrame in odb.getFrameIterator():
         >>>     field = MyField()
         >>>     odbField = field.odbFieldFromFrame(odbFrame)
         >>>     partData = odbField.getSubset(region=filterSet)
         >>>     field.readFromOdbField(partData)

        @param odbField: FieldOutput object in the odb. This may be created by
        self.odb.getSubset() or it may also be odbFrame.fieldOutputs[name]
        @returns: self
        """
        if self.position=="elemIP":
            # value at integration points
            newValues = defaultdict(list)
            if self.odbFieldInvariant is not None:
                for value in odbField.values:
                    newValues[getattr(value, self.odbFieldLabelAttr)].append(
                        (value.integrationPoint,
                         getattr(value, self.odbFieldInvariant)))
            elif self.odbFieldComponent is not None:
                for value in odbField.values:
                    newValues[getattr(value, self.odbFieldLabelAttr)].append(
                        (value.integrationPoint,
                         getattr(value, self.odbFieldDataAttr)[
                             self.odbFieldComponent]))
            elif self.dataType=="scalar":
                # scalar field
                for value in odbField.values:
                    newValues[getattr(value, self.odbFieldLabelAttr)].append(
                        (value.integrationPoint,
                         getattr(value, self.odbFieldDataAttr)))
            else:
                # vector or tensor field
                for value in odbField.values:
                    newValues[getattr(value, self.odbFieldLabelAttr)].append(
                        (value.integrationPoint,
                         list(getattr(value, self.odbFieldDataAttr))))

            # sort integration points in ascending order, extract the values
            self.update((label, [val for ip, val in
                                 sorted(valueList, key=lambda tup: tup[0])])
                        for label, valueList in newValues.iteritems())

        elif self.position=="elemSP":
            # values at integration points on beam or shell elements with
            # different section points
            newValues = defaultdict(list)
            try:
                if self.odbFieldInvariant is not None:
                    for value in odbField.values:
                        newValues[getattr(value, self.odbFieldLabelAttr)]. \
                            append((value.integrationPoint,
                                    value.sectionPoint.number,
                                    getattr(value, self.odbFieldInvariant)))
                elif self.odbFieldComponent is not None:
                    for value in odbField.values:
                        newValues[getattr(value, self.odbFieldLabelAttr)]. \
                            append((value.integrationPoint,
                                    value.sectionPoint.number,
                                    getattr(value, self.odbFieldDataAttr)
                                    [self.odbFieldComponent]))
                elif self.dataType=="scalar":
                    # scalar field
                    for value in odbField.values:
                        newValues[getattr(value, self.odbFieldLabelAttr)]. \
                            append((value.integrationPoint,
                                    value.sectionPoint.number,
                                    getattr(value, self.odbFieldDataAttr)))
                else:
                    # vector or tensor field
                    for value in odbField.values:
                        newValues[getattr(value, self.odbFieldLabelAttr)]. \
                            append((
                                value.integrationPoint,
                                value.sectionPoint.number,
                                list(getattr(value, self.odbFieldDataAttr))))
            except AttributeError, exc:
                if value.sectionPoint is None:
                    text = ("ERROR: Field %s is not defined at section points"
                            " (for element %d)." % (
                                self.fieldName,
                                getattr(value, self.odbFieldLabelAttr)))
                    msg(text)
                    raise ValueError(text)
                raise

            # sort integration points in ascending order, extract the values
            self.update((label, [ tup2[-1] for tup2 in sorted(
                valueList, key=lambda tup1: (tup1[0],tup1[1])) ] )
                        for label, valueList in newValues.iteritems())

        else:
            # nodal or element value
            if self.odbFieldInvariant is not None:
                self.update(((getattr(value, self.odbFieldLabelAttr),
                              getattr(value, self.odbFieldInvariant))
                             for value in odbField.values))
            elif self.odbFieldComponent is not None:
                self.update([(getattr(value, self.odbFieldLabelAttr),
                              getattr(value, self.odbFieldDataAttr)[
                                  self.odbFieldComponent])
                             for value in odbField.values])
            elif self.dataType=="scalar":
                # scalar field
                self.update([(getattr(value, self.odbFieldLabelAttr),
                              getattr(value, self.odbFieldDataAttr))
                             for value in odbField.values])
            else:
                # vector or tensor field
                self.update([(getattr(value, self.odbFieldLabelAttr),
                              list(getattr(value, self.odbFieldDataAttr)))
                             for value in odbField.values])
        return self


#---------------------------------------------------------------------
#--- module data and function
#... for reading field information from the odb
#... add-on for bae.field_01 to get field info from the odb and then
#... to also read values from the odb

# dict {odb field type: (Field.position, odb label attribute)}
_odbPositionToFieldPosLabelAttr = {
    odbAccess.NODAL: ("node", "nodeLabel"),
    odbAccess.WHOLE_ELEMENT: ("element", "elementLabel"),
    odbAccess.INTEGRATION_POINT: ("elemIP", "elementLabel"),
}

# dict {odb field type: Field.dataType}
_odbFieldTypeToFieldType = {
    odbAccess.SCALAR: "scalar",
    odbAccess.VECTOR: "vector",
    odbAccess.TENSOR_3D_SURFACE: "vector",
    odbAccess.TENSOR_3D_FULL: "tensor",
    }

# dict {number part of the component label: index in the data array}
_odbComponentToArrayIndex = {
    "1": 0, "2": 1, "3": 2,
    # tensor components in the odb: xx, yy, zz, xy, xz, yz
    "11": 0, "22": 1, "33": 2,
    "12": 3, "13": 4, "23": 5 }

def getOdbFieldClass(fieldName, odbFrame):
    r"""
    Create a new class based on (being child of) L{bae.field_01.Field} and
    L{OdbField} with all the info how to get the data from the odb.

    Returns a class, not an object! You can create multiple objects of this,
    one for each frame for example. You will then typically call their
    L{readFromFrame<OdbField.readFromFrame>} or
    L{readFromOdbField<OdbField.readFromOdbField>} methods to load actual data
    from the odb.

    The new class has the following class variables
     - fieldName: The fieldName argument.
     - odbFieldDataAttr: The name of the data attribute of the
    FieldValue-object that holds the actual data. It's "data" or
    "dataDouble".
     - odbFieldName: "S_MISES" -> "S", "U_1" -> "U", "SDV4" -> "SDV4"
     - odbFieldComponent: "S_MISES"->None, "SDV4"->None, "U_1"->0,
    "U_2"->1
     - odbFieldInvariant: "S_MISES"->"mises",
    "S_MIN_PRINCIPAL"->"minPrincipal", "SDV4"->None, "U_1"->None

    Usage:

     >>> odb = OdbReader("mydata.odb")
     >>> lastframe = odb.odb.steps["Step-2"].frame[10]
     >>> FieldU = getOdbFieldClass("U", lastframe)
     >>> frameFieldU = dict()
     >>> for frameNo, frame in odb.getFrameIterator([10, 12, -1]):
     ...     try:
     ...         odbField = frame.fieldOutputs[FieldU.odbFieldName]
     ...     except KeyError:
     ...         print "WARNING: No field U in frame %d." % frameNo
     ...         continue
     ...     newfield = FieldU()
     ...     newfield.readFromOdbField(odbField)
     ...     queryNodeNb = 10
     ...     print ("Displacement at node %d = %s"
     ...            % (queryNodeNb, newfield[queryNodeNb]))
     ...     frameFieldU[frameNo] = newfield

    In case there are values at different section points, i.e. integration
    point values on beam or shell elements, manually switch the position
    attribute of the newly created class from "elemIP" to "elemSP":
     >>> FieldS = getOdbFieldClass("S", lastframe)
     >>> FieldS.position = "elemSP"
     >>> currentStressField = FieldS()
     >>> # this can now be used to read data from the odb:
     >>> currentStressField.readFromFrame( odb.getFrame(currentFrameNo) )

    Integration point values are ordered according to ascending integration
    point index.

    Values at different section points are ordered according to ascending
    (integration pt index, section point number)-tuples, *only* if you changed
    the position attribute to "elemSP"!

    I.e. suppose you ordered stress values at section points 1, 3, 7, 11, 15
    for a beam element with 2 integration points you'll get a list of stress
    tensor values (i.e. lists by themselves) in the following order:
    [ (ip 1 sp 1), (ip 1 sp 3), ..., (ip 1 sp 15), (ip 2 sp 1), (ip 2 sp 3),
    ..., (ip 2 sp 15) ] each "(ip x sp y)" standing for one stress tensor.

    Querying values for elements that don't have section points (i.e. tets)
    with the position attribute set to "elemSP" will stop the program with an
    error! This inconvenience is due to speed considerations.

    @param fieldName: examples: "S", "V", "S_MISES", "S_MIN_PRINCIPAL",
    "U_1", "S_33", "SDV4"

    @param odbFrame: OdbFrame object in the odb where the field information
    like scalar/tensor, nodal/element should be looked up.

    @returns: A new class. It is a subclass of OdbField and two of the
    Field... classes.

    @Raises OdbFieldNotFound: if there is no field in the odb matching
    fieldName.

    @Note: This might not work for vector or tensor components in other
    than 3D solid elements because of the assumption in the
    odbComponentToArrayIndex dictionary! In this case use
    getScalarField(componentLabel="S11") to get the component.

    @Note: possible improvements: accept None as odbFrame. In this case
    get the data from a yet-to-create dictionary of the odb-module instead
    of the odb itself.
    """

    # we will collect additional class variables for the new class
    classValues = dict()

    #--- field, invariant and component can and will be derived from the name
    invariant = None
    component = None

    if fieldName.startswith("SDV"):
        field = fieldName
    elif '_' in fieldName:
        field, refinement = fieldName.split('_', 1)
        if refinement.isdigit():
            # component of a vector or tensor
            try:
                component = _odbComponentToArrayIndex[refinement]
            except KeyError:
                raise ValueError(
                    "Invalid component %s, must be one of %s"
                    % (refinement, ",".join(_odbComponentToArrayIndex)))
        else:
            # invariant (like magnitude or so)
            invariant = refinement.lower()
            try:
                while 1:
                    i = invariant.index('_')
                    invariant = (invariant[:i] + invariant[i+1].upper()
                                 + invariant[i+2:])
            except ValueError:
                pass
    else:
        field = fieldName

    classValues["odbFieldName"] = field
    classValues["odbFieldComponent"] = component
    classValues["odbFieldInvariant"] = invariant

# at the moment there is one read function that treates all cases:
# OdbField.readFromOdbField(). It could be more efficient to define different
# ones and assign them similar to that:
#     if invariant!=None:
#         classValues["readFromOdbField"] = OdbField.readInvariant
#     elif component!=None:
#         classValues["readFromOdbField"] = OdbField.readComponent
#     else:
#         classValues["readFromOdbField"] = OdbField.readField

    #--- type, precision and position will be looked up in the odb
    odbFields = odbFrame.fieldOutputs
    try:
        odbField = odbFields[field]
    except KeyError:
        raise OdbFieldNotFound(
            "No field %s (%s) in the specified frame."
            % (field, fieldName))

    # odbType
    if component is None and invariant is None:
        odbType = odbField.type
    else:
        odbType = odbAccess.SCALAR
    try:
        fieldType = _odbFieldTypeToFieldType[odbType]
    except KeyError:
        raise ValueError("odb field type %s not implemented yet." % odbType)

    # remove SDV_ from fieldName from now on
    if fieldName.startswith("SDV_"):
        fieldName = fieldName[4:]

    # get precision from the corresponding attribute of the first value
    precision = odbField.values[0].precision
    if precision==odbAccess.SINGLE_PRECISION:
        classValues["odbFieldDataAttr"] = "data"
    elif precision==odbAccess.DOUBLE_PRECISION:
        classValues["odbFieldDataAttr"] = "dataDouble"
    else:
        raise ValueError("Precision %s not implemented." % precision)

    # odbPosition, self.position
    positions = set([f_loc.position
                     for f_loc in odbField.locations])
    odbPosition = positions.pop()
    if len(positions)>0:
        warn(OdbReadWarn("Multiple positions for field %s. All but"
                         " position %s will be ignored."
                         % (field, odbPosition)))
    fieldPos, labelAttr = _odbPositionToFieldPosLabelAttr[odbPosition]
    classValues["odbFieldLabelAttr"] = labelAttr

    #--- now create the new class
    NewClass = Field.classFromPosType(
        fieldName, fieldPos, fieldType, addBaseClasses=(OdbField,))
    # append the class attributes collected in classValues
    for attr, val in classValues.iteritems():
        setattr(NewClass, attr, val)

    return NewClass


def getOdbFieldClassFromOdbField(
        odbField, component=None, invariant=None, fieldName=None):
    """Create a new class based on (being child of) L{bae.field_01.Field} and
    L{OdbField} with all the info how to get the data from the odb.

    Returns a class, not an object!

    See also L{getOdbFieldClass}, this function only differes in that it takes
    a FieldOutput object from the abaqus API.

    @param odbField: FieldOutput object
    @param component: array index if you want a component of a vector or
        tensor, see L{_odbComponentToArrayIndex}
    @param invariant: abaqus API attribute name like "minPrincipal", "mises",
        "magnitude"
    @param fieldName: will be copied to the corresponding L{OdbField} objects
        fieldName attribute. Defaults to odbField.name
    """

    # we will collect additional class variables for the new class
    classValues = dict()
    classValues["odbFieldName"] = odbField.name
    classValues["odbFieldComponent"] = component
    classValues["odbFieldInvariant"] = invariant

    # odbType
    if component is None and invariant is None:
        odbType = odbField.type
    else:
        odbType = odbAccess.SCALAR
    try:
        fieldType = _odbFieldTypeToFieldType[odbType]
    except KeyError:
        raise ValueError("odb field type %s not implemented yet." % odbType)

    if fieldName is None:
        fieldName = odbField.name

    # remove SDV_ from fieldName from now on
    if fieldName.startswith("SDV_"):
        fieldName = fieldName[4:]

    # get precision from the corresponding attribute of the first value
    precision = odbField.values[0].precision
    if precision==odbAccess.SINGLE_PRECISION:
        classValues["odbFieldDataAttr"] = "data"
    elif precision==odbAccess.DOUBLE_PRECISION:
        classValues["odbFieldDataAttr"] = "dataDouble"
    else:
        raise ValueError("Precision %s not implemented." % precision)

    # odbPosition, self.position
    positions = set([f_loc.position
                     for f_loc in odbField.locations])
    odbPosition = positions.pop()
    if len(positions)>0:
        warn(OdbReadWarn("Multiple positions for field %s. All but"
                         " position %s will be ignored."
                         % (fieldName, odbPosition)))
    fieldPos, labelAttr = _odbPositionToFieldPosLabelAttr[odbPosition]
    classValues["odbFieldLabelAttr"] = labelAttr

    #--- now create the new class
    NewClass = Field.classFromPosType(
        fieldName, fieldPos, fieldType, addBaseClasses=(OdbField,))
    # append the class attributes collected in classValues
    for attr, val in classValues.iteritems():
        setattr(NewClass, attr, val)

    return NewClass


#---------------------------------------------------------------------
#--- OdbReader variant
#... intended to read and interpolate data to a regular grid

class OdbReaderInterpolateToGrid(OdbReader):
    r"""Read data from the odb and interpolate to a regular grid or to a "grid"
    of unstructured points.
    This object stores all the administrative data and provides the means
    to retrieve the data from the odb.

     - initialize interpolation (also create filtersets)
     - initialize fields (output types)
     - read(field, odbFrame): Field object with regular grid reference

    As for the parent class OdbReader, you may access the Abaqus odb object as
    self.odb (as long as the self.close() function has not been called).

    Example 1; Read some Fields from one frame:
     >>> from bae.odb_02 import OdbReaderInterpolateToGrid
     >>> from bae.vtk_02 import FieldsOnStructPointGrid
     >>> from bae.log_01 import msg
     >>> 
     >>> #--- parameter section
     >>> 
     >>> odbName = "/work/abaqus/myproj/myproj.odb"
     >>> gridData = dict(
     >>>     boxName="CadiaM21_bigbox100m",
     >>>     firstPoint=[14000, 21000, 4700],
     >>>     lastPoint=[15000, 22000, 5700],
     >>>     spacing=100)
     >>> odbStepFrame = [2, 10]
     >>> fieldNames = ["SDV17", "SDV18", "SDV19", "SDV20"]
     >>> 
     >>> #--- end parameter section
     >>> 
     >>> 
     >>> odb = OdbReaderInterpolateToGrid(odbName)
     >>> 
     >>> # initialize the output grid
     >>> grid = odb.initOutputGrid(**gridData)
     >>> 
     >>> odbFrame = odb.getFrame2(*odbStepFrame)
     >>> output = FieldsOnStructPointGrid(mesh=grid)
     >>> for fieldName in fieldNames:
     >>>     msg("reading field %s" % fieldName)
     >>>     field = odb.getFieldFromOdb(fieldName=fieldName, odbFrame=odbFrame)
     >>>     output.updateFields(field)
     >>> 
     >>> outputFileName = "%s_frame10.vtk"%gridData["boxName"]
     >>> output.toVtk(outputFileName,
     ...              description="Data from Cadia M24 step-2, frame 10")
     >>> msg("Wrote results to %s" % outputFileName)

    Example 2; Read data from some frames, modify output fields:
     >>> from bae.odb_02 import OdbReaderInterpolateToGrid as Odb
     >>> from bae.vtk_02 import FieldsOnStructPointGrid as Vtk
     >>> from bae.field_01 import Field
     >>> from itertool import izip
     >>> odb = Odb("myjob.odb")
     >>> 
     >>> # initialize the output grid
     >>> grid = odb.initOutputGrid(boxName="myjob_smallbox2m",
     ...                           firstPoint=[1000, 200, -500],
     ...                           lastPoint=[1200, 350, -400],
     ...                           spacing=2)
     >>> 
     >>> # get base displacement
     >>> baseStep=2; baseFrame=0
     >>> baseU3 = odb.getFieldFromOdb(
     ...     fieldName="U_3", odbFrame=odb.getFrame2(baseStep, baseFrame))
     >>> 
     >>> # Type/class for objects that will later hold newly generated data
     >>> UU3Type = Field.classFromPosType("UU_3", "structPt", "scalar")
     >>> 
     >>> # main loop: read and interpolate and export data for some frames
     >>> framesList = [ [2, [10, 11, 12]], ]
     >>> for stepNb, frameNb, odbFrame in self.getFrameIterator2(framesList):
     >>>     fieldU3 = odb.getFieldFromOdb(fieldName="U_3", odbFrame)
     >>>     fieldUU3 = UU3Type((a-b) for a, b in izip(fieldU3, baseU3))
     >>>     output = Vtk(mesh=grid)
     >>>     output.updateFields(fieldU3, fieldUU3)
     >>>     output.toVtk("myjob_Uvert_F%02d.vtk"%frameNb,
     ...                  description="Uvert for frame %d"%frameNb)

    Example 3; Read field and material data
     >>> from bae.odb_02 import OdbReaderInterpolateToGrid
     >>> from bae.vtk_02 import FieldsOnStructPointGrid
     >>> 
     >>> odb = Odb("myjob.odb")
     >>> 
     >>> # initialize the output grid
     >>> grid = odb.initOutputGrid(boxName="myjob_smallbox2m",
     ...                           firstPoint=[1000, 200, -500],
     ...                           lastPoint=[1200, 350, -400],
     ...                           spacing=2)
     >>> 
     >>> # initalize output data buffer
     >>> output = FieldsOnStructPointGrid(mesh=grid)
     >>> 
     >>> # get field data from frame
     >>> odbFrame = odb.getFrame2(2, 10)
     >>> for fieldName in ["SDV17", "SDV18", "SDV19", "SDV20"]:
     >>>     msg("reading field %s" % fieldName)
     >>>     field = odb.getFieldFromOdb(fieldName=fieldName, odbFrame=odbFrame)
     >>>     output.updateFields(field)
     >>> 
     >>> # get the material data
     >>> (matNamesList, matNumberField) = odb.getMatFieldFromOdb()
     >>> output.updateFields(matNumberField)
     >>> 
     >>> # write data to vtk-file
     >>> outputFileName = "%s_frame10.vtk"%gridData["boxName"]
     >>> output.toVtk(outputFileName,
     ...              description="Data from step-2, frame 10")
     >>> 
     >>> # write material number to names list
     >>> outputFileNameMatList = "%s_matList.csv"%gridData["boxName"]
     >>> odb.writeMatNamesToCsv(outputFileNameMatList, matNamesList)

    Example 4; Read some Fields from one frame:
     >>> from bae.odb_02 import OdbReaderInterpolateToGrid
     >>> import csv
     >>> from bae.log_01 import msg
     >>> 
     >>> odbName = "/work/abaqus/myproj/myproj.odb"
     >>> pointCoords = [
     ...     map(float, row[:3])
     ...     for row in csv.reader(open("MyPoints.csv","rb"))]
     >>> odbStepFrame = [2, 10]
     >>> fieldNames = ["SDV17", "SDV18", "SDV19", "SDV20"]
     >>> 
     >>> odb = OdbReaderInterpolateToGrid(odbName)
     >>> odb.initOutputPoints(pointCoords, "MyPoints")
     >>> 
     >>> odbFrame = odb.getFrame2(*odbStepFrame)
     >>> for fieldName in fieldNames:
     >>>     msg("reading field %s" % fieldName)
     >>>     field = odb.getFieldFromOdb(fieldName=fieldName, odbFrame=odbFrame)
     >>>     doSomething(field)

    @ivar grid: A L{MeshStructuredPoints<bae.mesh_01.MeshStructuredPoints>}
      or L{MeshUnstructuredPoints<bae.mesh_01.MeshUnstructuredPoints>}
      object representing the output grid
    @ivar fieldNameTypes: A dictionary holding already created subclasses of
      L{bae.field_01.Field} that will be used for the resulting field by
      L{getFieldFromOdb}:
      {fieldName (like 'U_3', 'S_MISES') : corresponding subclass of Field}
    """

    def __init__(self, odb):
        r"""
        @param odb: file name of the odb or the already opened odb object.
        """
        OdbReader.__init__(self, odb)
        self.fieldNameTypes = dict()

    def initOutputGrid(self, **kwargs):
        """Initialize a regular grid of output points. Calculate the element
        coordinates for each output grid point and possibly store them for
        subsequent use in a pickle file. Or load this data from the so created
        pickle file.

        You may specify either of the following argument combinations:
         - firstPoint, lastPoint and spacing
         - gridPtNbs, origin and spacing

        Specify a interpolPickleName argument for speed up on subsequent runs.

        Creates the following instance variables:
         - grid: A mesh_01.MeshStructuredPoints object representing the output
              grid
         - mesh: mesh object of only those elements that include one of the
              grid points
         - elemCoords: list of (element number, elem coord)-tuples for each
              grid point
         - filterElset: odbset containing all elements of which output data
              will be read from the odb.
         - filterNset: odbset containing all nodes of which output data
              will be read from the odb.

        @kwarg gridPtNbs: A triple of integers specifying the number of grid
          points in x, y, z direction
        @kwarg origin: first point (min x, y, z) (synonym to firstPoint)
        @kwarg box: A bounding box [[xmin,ymin,zmin], [xmax, ymax, zmax]].
          You need to supply spacing as well.
        @kwarg spacing: A triple of doubles: Point spacing in x, y, z direction
          or just one number the same spacing in each direction.
        @kwarg firstPoint: first point (min x, y, z) (synonym to origin)
        @kwarg lastPoint: last point (max x, y, z)

        @kwarg interpolPickleName: A file containing the
          position of the grid points in the mesh (i.e. element numbers and
          element coordinates) will be saved and used on subsequent
          invocations. Its filename will be assembled using this string.
          See L{bae.mesh_01.InterpolMeshToPoints} for details.

          If the mesh or the requested points have changed remove
          this file. If you don't specify this argument the point position will
          be calculated and no pickle file will be generated for subsequent
          invocations.
        @kwarg boxName: synonym for interpolPickleName if interpolPickleName
          is not given. interpolPickleName takes precedence.
        @kwarg meshinfopickleName: deprecated synonym for interpolPickleName

        @returns: the grid, a @L{bae.mesh_01.MeshStructuredPoints} object
        """

        # get pickelName from meshinfopickleName or boxname
        pickleName = None
        if "interpolPickleName" in kwargs:
            pickleName = kwargs["interpolPickleName"]
            # note: it's save to delete from the kwargs dict because even if
            # __init__() is called with a dictionary argument, i.e. as
            # __init__(self, **myparams) kwargs will be a copy of (not a
            # reference to) myparams. Automatically...
            del kwargs["interpolPickleName"]
        elif "boxName" in kwargs:
            pickleName = kwargs["boxName"]
        if "boxName" in kwargs:
            del kwargs["boxName"]

        # for compatibility also process deprecated arguments:
        elif "meshinfopickleName" in kwargs:
            pickleName = kwargs["meshinfopickleName"]
            del kwargs["meshinfopickleName"]

        if pickleName is None:
            msg("Not using using pickle file for the interpolation.",
                debugLevel=20)
        else:
            msg("Using pickle file %s" % pickleName, debugLevel=20)
        grid = getMesh(**kwargs)

        self.interpolation = InterpolMeshToPoints(
            grid=grid,
            getMeshCallBack=self.getAbqModel,
            gmCallBackArgs=dict(recognizedBlocks="NODE,ELEMENT"),
            pickleName=pickleName
            )

        # for backwards compatibility:
        self.elemCoords = self.interpolation.elemCoords

        # create a filterset for the output data
        msg("Creating filterset for the output data.")
        self.filterElset = self.odbSetFromElements(
            self.interpolation.mesh.elNodes.iterkeys())
        self.filterNset = self.odbSetFromNodes(
            self.interpolation.mesh.nodeCoords.iterkeys())
        msg("len(filterElset): %d"%len(self.filterElset.elements),debugLevel=10)
        msg("len(filterNset): %d"%len(self.filterNset.nodes), debugLevel=10)

        return grid

    def initOutputPoints(self, pointsCoords,
                         interpolPickleName=None, pointsName=None,
                         elems=None):
        """Initialize output points. Calculate the element coordinates for each
        output point and possibly store them for subsequent use in a pickle
        file. Or load this data from the so created pickle file.

        @param pointsCoords: list of point coordinates ([x,y,z]-lists)
        @param interpolPickleName: A file "%s.interpol"%pointsName containing
          the position of the output points in the mesh (i.e. element numbers
          and element coordinates) will be saved and used on subsequent
          invocations. If the mesh or the requested points have changed remove
          this file. If you don't specify this argument the point position will
          be calculated and no pickle file will be generated for subsequent
          invocations.
        @param pointsName: synonym for interpolPickleName
        @param elems: Elements to be considered. Points outside those elements
          will not return a value like points outside of the mesh. If None
          (default) consider all elements of the mesh. Otherwise must be
          anything that Model.getUnionSet("elset",...) accepts as input.
        @returns: self.grid
        """

        pickleName = None
        # get pickleName from meshinfopickleName or boxname
        if "interpolPickleName" is not None:
            pickleName = interpolPickleName
            del interpolPickleName
        elif "boxName" is not None:
            pickleName = pointsName
            del pointsName

        # create grid object
        msg("Creating unstructured points 'grid' with %d points."
            % len(pointsCoords))
        grid = MeshUnstructuredPoints(pointsCoords)
        if not len(grid):
            raise Exception("No points in the grid.")

        if elems is None:
            getMeshCallBack = self.getAbqModel
            gmCallBackArgs = dict(recognizedBlocks="NODE,ELEMENT")
        else:
            def getMeshCallBack():
                model = self.getAbqModel(recognizedBlocks="NODE,ELEMENT,ELSET")
                elset = model.getUnionSet("elset", elems)
                msg("filterset(s) %s containing %d elements"
                    % (elems, len(elset)), debugLevel=10)
                model.initializePointSearch(
                    box=grid.getBoundingBox(), elems=elset)
                return model
            gmCallBackArgs = dict()

        self.interpolation = InterpolMeshToPoints(
            grid=grid,
            getMeshCallBack=getMeshCallBack,
            gmCallBackArgs=gmCallBackArgs,
            pickleName=pickleName
            )

        # for backwards compatibility:
        self.elemCoords = self.interpolation.elemCoords

        # create a filterset for the output data
        msg("Creating filterset for the output data.")
        self.filterElset = self.odbSetFromElements(
            self.interpolation.mesh.elNodes.iterkeys())
        self.filterNset = self.odbSetFromNodes(
            self.interpolation.mesh.nodeCoords.iterkeys())
        msg("len(filterElset): %d"%len(self.filterElset.elements),debugLevel=10)
        msg("len(filterNset): %d"%len(self.filterNset.nodes), debugLevel=10)

        return grid

    def getFieldClasses(self, fieldName, odbFrame):
        """Return proper classes to hold the raw odb data and the interpolated
        data for the specified field name.

        Merely for internal use.

        @Raises OdbFieldNotFound: if there is no field in the odb matching
        fieldName.
        """

        try:
            FieldClassInp, FieldClassOut = self.fieldNameTypes[fieldName]
        except KeyError:
            # create the classes
            FieldClassInp = getOdbFieldClass(fieldName, odbFrame)
            FieldClassOut = Field.classFromPosType(
                fieldName, 'structPt', FieldClassInp.dataType)
            self.fieldNameTypes[fieldName] = (FieldClassInp, FieldClassOut)

        return (FieldClassInp, FieldClassOut)

    def getFieldFromOdb(self, fieldName, odbFrame):
        """read the data for the given field from the odb and interpolate to
        the grid.

        @Raises OdbFieldNotFound: if there is no field in the odb matching
        fieldName.
        """

        FieldClassInp, FieldClassOut = self.getFieldClasses(fieldName,odbFrame)

        msg("Reading values for field %s from the odb."%fieldName,debugLevel=1)
        field = FieldClassInp()

        # get the right odb field instance (right field, right frame)
        try:
            odbField = field.odbFieldFromFrame(odbFrame)
        except KeyError:
            raise Exception("ERROR: Field %s not found but getFieldClasses()"
                            " succeded. Weird..." % fieldName)

        # get the subset containing only values needed for interpolation
        if field.position=="node":
            partData = odbField.getSubset(region=self.filterNset)
        else:
            partData = odbField.getSubset(region=self.filterElset)
        msg("len(odbField.values): %d" % len(odbField.values), debugLevel=10)
        msg("len(partData.values): %d" % len(partData.values), debugLevel=10)

        # read values from odb
        field.readFromOdbField(partData)
        msg("Read %d values. Interpolating data." % len(field), debugLevel=1)

        # interpolate to the grid
        ptValues = FieldClassOut(field.interpolate(self.interpolation))
        return ptValues

    def getMatFieldFromOdb(self, nameFilter=None, firstMatCode=1, noMatCode=-1,
                           filterUnused=False):
        """Read material section name from the specified odb at grid point
        positions.

        @param nameFilter: This is for material name conversion. For example
          if there are region variants that shall be treated as one. I.e. if
          SEDIMENT, SEDIMENT_CAVE and SEDIMENT_UCUT shall all be treated as
          SEDIMENT.

          If this is a function it will be called as a conversion function
          taking the original material name from the odb as only
          argument and returning the material name it should be treated as.

        @param firstMatCode: Material code number of the first material in
          the list of material names. Or if you prefer: Offset of the material
          code numbers in the matNumberField field and the index of the
          corresponding material in the matNamesList.

          I.e. if firstMatCode=1 then all points with assigned material number
          of 1 are of the first material type in the list matNamesList which is
          matNamesList[0]. All points having material number 5 are of material
          type matNamesList[4].

          Specify firstMatCode=0 if you want the material numbers in
          matNumberField to directly specify the index in matNumberField. The
          default value of 1 seems to be more appropriate for the usual use
          case of exporting matNumberField to a vtk file and matNamesList to a
          seperate csv file.

        @param noMatCode: Material number to represent points that are not
          being assinged any material to (because they sit outside the mesh).

        @param filterUnused: If True then remove material codes from the list
          of material names that are not present/used in/by the given grid.

        @returns: (matNamesList, matNumberField)-tuple
          matNamesList is the list of material names.
          matNumberField is the field (i.e. a list) of material
        """
        (matNameList, matNumberElemField) = OdbReader.getMatFieldFromOdb(
            self, nameFilter, firstMatCode, noMatCode, secTypeFilter=["solid"])

        msg("Interpolating field of material numbers to grid points")
        MatNumberFieldClass = Field.classFromPosType(
            "Material", "structPt", "scalar")
        matNumberField = MatNumberFieldClass(
            matNumberElemField.interpolate(self.interpolation,
                                           externalPtValue=noMatCode))

        if filterUnused:
            msg("Filtering material types keeping only required types.")
            usedMatNumbers = set(matNumberField)

            nbMatNames = len(matNameList)
            if noMatCode in usedMatNumbers:
                nbMatNames += 1

            if len(usedMatNumbers) != nbMatNames:
                msg("Filter diagnostic:...\nusedMatNumbers: %s\n"
                    "matNameList: %s"
                    % (str(sorted(usedMatNumbers)),
                       ",".join(matNameList)), debugLevel=1)
                oldToNew = dict()
                oldToNew[noMatCode] = noMatCode
                newMatNameList = list()
                for i, matName in enumerate(matNameList):
                    if (i+firstMatCode) in usedMatNumbers:
                        oldToNew[i+firstMatCode] = (
                            len(newMatNameList)+firstMatCode)
                        newMatNameList.append(matName)
                msg("Filter diagnostic:... oldToNew: %s"
                    % oldToNew, debugLevel=1)
                matNumberField = MatNumberFieldClass(
                    oldToNew[matCode] for matCode in matNumberField)
                msg("Originally we had %d material types. Only %d of them are"
                    " actually being used. I removed the other: %s."
                    % (len(matNameList), len(newMatNameList), ", ".join(
                        set(matNameList).difference(newMatNameList))))
                matNameList = newMatNameList

        msg("Identified %d different material codes for %d grid points."
            % (len(matNameList), len(self.interpolation.elemCoords)))
        return (matNameList, matNumberField)
