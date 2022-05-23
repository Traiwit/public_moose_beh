r"""Proxy script to be run with abaqus python. This script will be launched by
the bae.odb_03.OdbReader class with two file descriptors as argument: the input
(for commands) and output (for results from the odb):

abaqus python odb_proxy_03.py fd_inp fd_out
"""

__version__ = "1.01"
_version_history_ = r"""\
Versions:

1.01 new     - based on odb_02 version 1.08
"""

import os, sys
from struct import Struct
from array import array

# provides pickle (=cPickle.dumps) and unpickle (=cPickle.loads)
from cPickle import loads as unpickle, dumps as pickle_dumps
pickle = lambda obj: pickle_dumps(obj, 2)

from collections import defaultdict

import odbAccess

class OdbReader(object):

    def __init__(self):
        self.lenStruct = Struct("i")
        self._fd_inp = int(sys.argv[1])
        self._fd_out = int(sys.argv[2])

        self.odbObjContainer = dict()

    def sendStr(self, string):
        os.write(self._fd_out, self.lenStruct.pack(len(string)))
        os.write(self._fd_out, string)

    def receiveStr(self):
        length = self.lenStruct.unpack(
            os.read(self._fd_inp, self.lenStruct.size))[0]
        string = ""
        while length>0:
            addstr = os.read(self._fd_inp, length)
            length -= len(addstr)
            string += addstr
        return string

    def mainLoop(self):
        """Receive messages and run the corresponding commands.
        """
        while 1:

            cmd = self.receiveStr()
            print "proxy received: %s" % cmd #####
            if cmd=="terminate":
                break

            kwargs = unpickle(self.receiveStr())
            print "proxy received cmd kwargs: %s" % kwargs #####

            print "proxy running cmd." #####
            result = getattr(self, cmd)(**kwargs)
            self.sendStr(pickle(result))

        # cleanup: close odb
        self.odb.close()
        print "proxy now stoppping." #####
        return

    def newOdbObject(self, odbObj):
        """store the given odb object in the container and return it's id
        (key in the container self.odbObjContainer)
        """
        try:
            odbId = max(self.odbObjContainer) + 1
        except ValueError:
            # first id if odbObjContainer empty
            odbId = 1
        self.odbObjContainer[odbId] = odbObj
        return odbId

    def removeOdbObject(self, odbId):
        del self.OdbObjContainer[odbId]

    #-------------------------------------------------------------

    def openOdb(self, odb):
        self.odb = odbAccess.openOdb(path=odb, readOnly=True)

    def _getOdbPartInstance(self):
        r"""Returns the first instance in the rootAssembly object of the odb.
        This contains the model definition, i.e. the mesh, material props and
        so on.

         >>> odb = OdbReader("mydata.odb")
         >>> instance = odb._getOdbPartInstance()
         >>> # read node info
         >>> for node in instance.nodes:
         >>>     output.write("%d:%s\n" % (node.label, list(node.coordinates)))
         >>> # read element info
         >>> for elem in instance.elements:
         >>>     output.write("%d:%s\n" % (elem.label, elem.type))

        @Note: There is a getMesh() - method to read the model data. It's
        probably more convenient if it serves your needs.
        """
        return self.odb.rootAssembly.instances.values()[0]

    def getMesh(self):
        r"""reads nodes and elements from the odb model data,

        assumes there is only one part instance in the model

        ideas for an additional future getAbqModel method:
         - add the possibility to selectively read only elements, nodes (see
           recognizedBlocks argument of abq_model_02.Model.read() / write())
         - read elsets, nsets, surfaces, all the rest
        """

        # mesh data from the first instance
        part_instance = self._getOdbPartInstance()

        # find dimensionality of our problem
        if part_instance.embeddedSpace != odbAccess.THREE_D:
            msg = "Can't handle dimensionality %s."%part_instance.embeddedSpace
            raise Exception(msg)

        # check coordinate values type code (of the last node)
        assert part_instance.nodes[0].coordinates.typecode()=="d", \
            "Coordinates must be doubles."

        # extract undeformed node coordinates
        packLabel = Struct("L").pack
        nodes = "".join(packLabel(node.label)+node.coordinates.tostring()
                        for node in part_instance.nodes)
        self.sendStr(nodes)

        # extract element/cell connectivity
        elTypes = _CountingDict()
        connectivity = "".join(
            array("L", [elt.label, elTypes[elt.type], len(elt.connectivity)
                        ] + list(elt.connectivity) ).tostring()
            for elt in part_instance.elements )

        self.sendStr(pickle(elTypes.sortedkeys()))
        self.sendStr(connectivity)
        return

    def getStepsList(self):
        """Returns a list of (step name, number of frames) tuples for all steps
        in the odb.
        """
        result = [ (stepname, len(step.frames))
                   for stepname, step in self.odb.steps.items()]
        return result

    def getFrame(self, stepName, frameNo):
        return self.newOdbObject(self.odb.steps[stepName].frames[frameNo])

    def initFrameIter(self, stepName, frameIdList):
        """call this to initialize iteration over specified frames.

        Usage:
         >>> self.simpleCmd("initFrameIter", stepName=step, frameIdList=frames)
         >>> while 1:
         >>>     frNo, odbFr = self.simpleCmd("frameIterNext")
         >>>     if frNo is None:
         >>>         break
         >>>     dosomething(frNo, odbFr)
        """
        odbStep = self.odb.steps[stepName]
        if frameIdList is None:
            frameIdList = range(len(odbStep.frames))
        else:
            nbFrs = len(odbStep.frames)
            frameIdList = [i for i in frameIdList if -nbFrs<=i<nbFrs]
        frameIdList.reverse()
        self.currentOdbStep = (odbStep, frameIdList)

    def frameIterNext(self):
        """Next frame of the current step as initialized with
        self.initFrameIter() or None.

        See L{OdbReader.initFrameIter} for details.
        """
        try:
            frameNo = self.currentOdbStep[1].pop()
        except IndexError:
            # no more frame
            del self.currentOdbStep
            return (None, None)
        except AttributeError:
            # initFrameIter not called yet
            return None
        else:
            return (frameNo, self.newOdbObject(odbStep.frames[frameNo]))
        























    #---------------------------------------------------------------------
    #{ get odb field info
    #... for reading field information from the odb
    #... add-on for bae.field_01 to get field info from the odb and then
    #... to also read values from the odb

    # dict {odb field type: (Field.position, odb label attribute)}
    _odbPositionToFieldPosLabelAttr = {
        odbAccess.NODAL : ("node", "nodeLabel"),
        odbAccess.WHOLE_ELEMENT : ("element", "elementLabel"),
        odbAccess.INTEGRATION_POINT : ("elemIP", "elementLabel"),
    }

    # dict {odb field type: Field.dataType}
    _odbFieldTypeToFieldType = {
        odbAccess.SCALAR : "scalar",
        odbAccess.VECTOR : "vector",
        odbAccess.TENSOR_3D_SURFACE : "vector",
        odbAccess.TENSOR_3D_FULL : "tensor",
        }

    # dict {number part of the component label: index in the data array}
    _odbComponentToArrayIndex = {
        "1" : 0, "2" : 1, "3" : 2,
        # tensor components in the odb: xx, yy, zz, xy, xz, yz
        "11" : 0, "22" : 1, "33" : 2,
        "12" : 3, "13" : 4, "23" : 5 }

    def getOdbFieldInfo(self, fieldName, odbFrame):
        r"""
        Create a new child class of OdbField with all the info how to look up
        the given field in the odb.

        Returns a class, not an object! You can create multiple objects of this,
        one for each frame for example.

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

        @returns: A new class. It is a subclass of _OdbField and two of the
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

        if '_' in fieldName:
            field, refinement = fieldName.split('_', 1)
            if refinement.isdigit():
                # component of a vector or tensor
                try:
                    component = self._odbComponentToArrayIndex[refinement]
                except KeyError:
                    raise ValueError(
                        "Invalid component %s, must be one of %s"
                        % (refinement, ",".join(self._odbComponentToArrayIndex)))
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
    # _OdbField.readFromOdbField(). It could be more efficient to define different ones and
    # assign them similar to that:
    #     if invariant!=None:
    #         classValues["readFromOdbField"] = _OdbField.readInvariant
    #     elif component!=None:
    #         classValues["readFromOdbField"] = _OdbField.readComponent
    #     else:
    #         classValues["readFromOdbField"] = _OdbField.readField

        #--- type, precision and position will be looked up in the odb
        odbFields = self.odbObjContainer[odbFrame].fieldOutputs
        try:
            odbField = odbFields[field]
        except KeyError:
            ########## don't store anything, return Errormarker
            ### old: to be removed
            raise OdbFieldNotFound(
                "No field %s (%s) in the specified frame."
                % (field, fieldName))

        # odbType
        if component==None and invariant==None:
            odbType = odbField.type
        else:
            odbType = odbAccess.SCALAR
        try:
            fieldType = _odbFieldTypeToFieldType[odbType]
        except KeyError:
            raise ValueError("odb field type %s not implemented yet." % odbType)

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
        fieldPos, labelAttr = self._odbPositionToFieldPosLabelAttr[odbPosition]
        classValues["odbFieldLabelAttr"] = labelAttr

        #--- now create the new class
        NewClass = Field.classFromPosType(
            fieldName, fieldPos, fieldType, addBaseClasses=(_OdbField,))
        # append the class attributes collected in classValues
        for attr, val in classValues.iteritems():
            setattr(NewClass, attr, val)

        return NewClass

    #} end of stuff to get odb field info



#-------------------------------------------------------------


class OdbSet(object):
    pass
class OdbElset(OdbSet):
    pass
class OdbNset(OdbSet):
    pass

class OdbField(object):

    def getFieldFromOdb(self, fieldName, odbFrame):
        """read the data for the given field from the odb and interpolate to
        the grid.

        @Raises OdbFieldNotFound: if there is no field in the odb matching
        fieldName.
        """
        pass

    def getSubset(self, region):
        """

        @returns: OdbField object
        """
        pass


#-------------------------------------------------------------

class _CountingDict(defaultdict):
    """This is a dictionary that stores an automatically increasing counter to
    each newly inserted key. Since being derived from defaultdict insertion is
    performed by querying the value for a key that is not yet present.

    I.e.:
     >>> d = _CountingDict()
     >>> print d[5], d[4], d[1], d[5], d[1], d[2]
     0, 1, 2, 0, 2, 3
     >>> print d.sortedkeys()
     [5, 4, 1, 2]
    """
    def __init__(self, *args):
        dict.__init__(self, *args)
    def __missing__(self, key):
        new = len(self)
        self[key]=new
        return new
    def sortedkeys(self):
        return [y[0] for y in sorted(self.iteritems(), key=lambda x:x[1])]

#-------------------------------------------------------------

odb = OdbReader()
odb.mainLoop()


exit()
#-------------------------------------------------------------


from warnings import warn

from bae.future_01 import defaultdict
from bae.misc_01 import Container
import bae.abq_model_02 as abq_model
from bae.field_01 import Field
from bae.mesh_01 import MeshStructuredPoints


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
    """

    def __init__(self, odb):
        """Constructor

        @param odb: file name of the odb or the already opened odb object.
        """

        if isinstance(odb, basestring):
            self.odb = odbAccess.openOdb(path=odb, readOnly=True)
        elif str(type(odb))=="<type 'Odb'>":
            self.odb = odb
        else:
            raise ValueError(
                "odb argument must either be a string specifying the file name"
                " or an open odb. Instead it's of type %s." % type(odb))


    def getFrame(self, stepNo, frameNo):
        r"""Get the OdbFrame object associated with the given step and frame
        numbers.

        @param stepNo: step number. Identifies the step in the odb by its name
          which must be "Step-%s"%stepNo. Usually the first step in the odb
          has got stepNo=1
        @param frameNo: frame number in the step. Negative values can be used
          to count from the end like usual list indices.
        """
        stepName = "Step-%s" % stepNo
        return self.odb.steps[stepName].frames[frameNo]


    def getFrameIterator(self, framesList=None):
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

          The step number is the number in the step name of the step in the
          odb. I.e. step number 3 identifies odb step "Step-3". In the ordinary
          two steps analysis the possible values are 1 and 2. Note that if the
          first step in the odb is called step-4 the corresponding step number
          is 4.

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
            odbStep = self.odb.steps["Step-%s" % stepNo]
            if frameIdList is None:
                frameIdList = range(len(odbStep.frames))
            else:
                nbFrs = len(odbStep.frames)
                frameIdList = [i for i in frameIdList if -nbFrs<=i<nbFrs]
            for frameNb in frameIdList:
                odbFrame = odbStep.frames[frameNb]
                yield (stepNo, frameNb, odbFrame)


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
        if frameIds==None:
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
        Create a new child class of OdbField with all the info how to look up
        the given field in the odb.

        Returns a class, not an object! You can create multiple objects of this,
        one for each frame for example.

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

        @returns: A new class. It is a subclass of _OdbField and two of the
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

    def getAbqModel(self, logfile="UseTheSameLogFile", recognizedBlocks = None):
        r"""get an abq_model_02.Model-object from the odb model data,

        reads nodes and elements

        assumes there is only one part instance in the model

        ideas for the future:
         - add the possibility to selectively read only elements, nodes (see
           recognizedBlocks argument of abq_model_02.Model.read() / write())
         - read elsets, nsets, surfaces, all the rest

        @param recognizedBlocks: Is being ignored at the moment. For future use.

        @Note: parameter logfile is being ignored. It's there for compatibility
        reasons only.
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
        msgTemplate = ("Read %s/%d nodes from the odb."
                       % ("%d", len(part_instance.nodes)))
        ticker = MsgTicker(msgTemplate, printLast=True)
        for cnt, node in enumerate(part_instance.nodes):
            model.nodeCoords[node.label] = list(node.coordinates)
            ticker.msg(cnt+1)
        del ticker

        # extract element/cell connectivity
        msgTemplate = ("Read %s/%d elements from the odb."
                       % ("%d", len(part_instance.elements)))
        ticker = MsgTicker(msgTemplate, printLast=True)
        for cnt, elt in enumerate(part_instance.elements):
            model.elNodes[elt.label] = list(elt.connectivity)
            model.elType[elt.label] = str(elt.type)
            model.typeEl[elt.type].add(elt.label)
            ticker.msg(cnt+1)
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

        elset = list(elset)
        elset.sort()

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

class _OdbField(object):
    r"""Additional base class for Field classes that should read their data from
    the odb. Has additional attributes and methods to load the field data from
    the odb.

    @Note: You will usually use dynamically created classes (types) based on
    class Field and _OdbField. The class (type) will be supplied by
    OdbReader.getOdbFieldClass(). Objects (instances) of those types can not
    easily be pickled using the ordinary pickle.dump() method. Use the
    pickle_dump() and pickle_load() methods instead.
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
            if self.odbFieldInvariant!=None:
                for value in odbField.values:
                    newValues[getattr(value, self.odbFieldLabelAttr)].append(
                        (value.integrationPoint,
                         getattr(value, self.odbFieldInvariant)))
            elif self.odbFieldComponent!=None:
                for value in odbField.values:
                    newValues[getattr(value, self.odbFieldLabelAttr)].append(
                        (value.integrationPoint,
                         getattr(value, self.odbFieldDataAttr)[self.odbFieldComponent]))
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
            self.update((label, [val for ip, val in sorted(valueList, key=lambda tup: tup[0])])
                        for label, valueList in newValues.iteritems())

        elif self.position=="elemSP":
            # values at integration points on beam or shell elements with
            # different section points
            newValues = defaultdict(list)
            if self.odbFieldInvariant!=None:
                for value in odbField.values:
                    newValues[getattr(value, self.odbFieldLabelAttr)].append(
                        (value.integrationPoint, value.sectionPoint.number,
                         getattr(value, self.odbFieldInvariant)))
            elif self.odbFieldComponent!=None:
                for value in odbField.values:
                    newValues[getattr(value, self.odbFieldLabelAttr)].append(
                        (value.integrationPoint, value.sectionPoint.number,
                         getattr(value, self.odbFieldDataAttr)[self.odbFieldComponent]))
            elif self.dataType=="scalar":
                # scalar field
                for value in odbField.values:
                    newValues[getattr(value, self.odbFieldLabelAttr)].append(
                        (value.integrationPoint, value.sectionPoint.number,
                         getattr(value, self.odbFieldDataAttr)))
            else:
                # vector or tensor field
                for value in odbField.values:
                    newValues[getattr(value, self.odbFieldLabelAttr)].append(
                        (value.integrationPoint, value.sectionPoint.number,
                         list(getattr(value, self.odbFieldDataAttr))))

            # sort integration points in ascending order, extract the values
            self.update((label, [ tup2[-1] for tup2 in sorted(
                            valueList, key=lambda tup1: (tup1[0],tup1[1])) ] )
                        for label, valueList in newValues.iteritems())

        else:
            # nodal or element value
            if self.odbFieldInvariant!=None:
                self.update(((getattr(value, self.odbFieldLabelAttr),
                              getattr(value, self.odbFieldInvariant))
                             for value in odbField.values))
            elif self.odbFieldComponent!=None:
                self.update([(getattr(value, self.odbFieldLabelAttr),
                              getattr(value, self.odbFieldDataAttr)[self.odbFieldComponent])
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
    odbAccess.NODAL : ("node", "nodeLabel"),
    odbAccess.WHOLE_ELEMENT : ("element", "elementLabel"),
    odbAccess.INTEGRATION_POINT : ("elemIP", "elementLabel"),
}

# dict {odb field type: Field.dataType}
_odbFieldTypeToFieldType = {
    odbAccess.SCALAR : "scalar",
    odbAccess.VECTOR : "vector",
    odbAccess.TENSOR_3D_SURFACE : "vector",
    odbAccess.TENSOR_3D_FULL : "tensor",
    }

# dict {number part of the component label: index in the data array}
_odbComponentToArrayIndex = {
    "1" : 0, "2" : 1, "3" : 2,
    # tensor components in the odb: xx, yy, zz, xy, xz, yz
    "11" : 0, "22" : 1, "33" : 2,
    "12" : 3, "13" : 4, "23" : 5 }
    
def getOdbFieldClass(fieldName, odbFrame):
    r"""
    Create a new child class of OdbField with all the info how to look up
    the given field in the odb.

    Returns a class, not an object! You can create multiple objects of this,
    one for each frame for example.

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

    @returns: A new class. It is a subclass of _OdbField and two of the
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

    if '_' in fieldName:
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
# _OdbField.readFromOdbField(). It could be more efficient to define different ones and
# assign them similar to that:
#     if invariant!=None:
#         classValues["readFromOdbField"] = _OdbField.readInvariant
#     elif component!=None:
#         classValues["readFromOdbField"] = _OdbField.readComponent
#     else:
#         classValues["readFromOdbField"] = _OdbField.readField

    #--- type, precision and position will be looked up in the odb
    odbFields = odbFrame.fieldOutputs
    try:
        odbField = odbFields[field]
    except KeyError:
        raise OdbFieldNotFound(
            "No field %s (%s) in the specified frame."
            % (field, fieldName))

    # odbType
    if component==None and invariant==None:
        odbType = odbField.type
    else:
        odbType = odbAccess.SCALAR
    try:
        fieldType = self._odbFieldTypeToFieldType[odbType]
    except KeyError:
        raise ValueError("odb field type %s not implemented yet." % odbType)

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
        fieldName, fieldPos, fieldType, addBaseClasses=(_OdbField,))
    # append the class attributes collected in classValues
    for attr, val in classValues.iteritems():
        setattr(NewClass, attr, val)

    return NewClass


#---------------------------------------------------------------------
#--- OdbReader variant
#... intended to read and interpolate data to a regular grid

class OdbReaderInterpolateToGrid(OdbReader):
    r"""Read data from the odb and interpolate to a regular grid.
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
     ... gridData = dict(
     ...     boxName="CadiaM21_bigbox100m",
     ...     firstPoint=[14000, 21000, 4700],
     >>>     lastPoint=[15000, 22000, 5700],
     >>>     spacing=100)
     >>> odbStepFrame = [2, 10]
     >>> fieldNames = ["SDV17", "SDV18", "SDV19", "SDV20"]
     ... 
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
     ... 
     ... outputFileName = "%s_frame00.vtk"%gridData["boxName"]
     ... output.toVtk(outputFileName,
     >>>              description="Data from Cadia M24 step-1, frame 0")
     >>> msg("Wrote results to %s" % outputFileName)

    Example 2; Read data from some frames, modify output fields:
     >>> from bae.odb_02 import OdbReaderInterpolateToGrid as Odb
     >>> from bae.vtk_02 import FieldsOnStructPointGrid as Vtk
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

    @ivar grid: A L{MeshStructuredPoints<bae.mesh_01.MeshStructuredPoints>}
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
        """
        You may specify either of the following argument combinations:
         - firstPoint, lastPoint and spacing
         - gridPtNbs, origin and spacing

        Specify a boxName argument for speed up on subsequent runs.

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
        @kwarg spacing: A triple of doubles: Point spacing in x, y, z direction
          or just one number the same spacing in each direction.
        @kwarg firstPoint: first point (min x, y, z) (synonym to origin)
        @kwarg lastPoint: last point (max x, y, z)

        @kwarg boxName: A file "%s_coordinates.pickle"%boxName containing the
          position of the grid points in the mesh (i.e. element numbers and
          element coordinates) will be saved and used on subsequent
          invocations. If the mesh or the requested points have changed remove
          this file. If you don't specify this argument the point position will
          be calculated and no pickle file will be generated for subsequent
          invocations.

        @returns: self.grid
        """

        def getOrigin(kwargs):
            try:
                return kwargs["origin"]
            except KeyError:
                pass

            try:
                return kwargs["firstPoint"]
            except KeyError:
                raise KeyError("Neither origin nor firstPoint in keyw args.")

        def getGridPtNb(kwargs, spacing):

            # try "lastPoint"
            try:
                return [int((x1-x0)/dx)+1 for x0, x1, dx
                            in zip(origin, kwargs["lastPoint"], spacing)]
            except KeyError:
                pass

            # try "gridPtNbs"
            try: 
                return kwargs["gridPtNbs"]
            except KeyError:
                raise KeyError("Neither lastPoint nor gridPtNbs in keyw args.")

        try:
            # "origin" or "firstPoint" must be there in any case
            origin = getOrigin(kwargs)

            # "spacing" must be there
            spacing = kwargs["spacing"]
            if isinstance(spacing, (int, float)):
                spacing = [spacing,spacing,spacing]

            # "lastPoint" or "gridPtNbs" must be there
            gridPtNb = getGridPtNb(kwargs, spacing)

        except KeyError:
            raise ValueError(
                "Unsufficient arguments to define the output grid.")

        # create grid object
        msg("Creating grid points list.")
        self.grid = MeshStructuredPoints(
            gridPtNb=gridPtNb, origin=origin, spacing=spacing)
        gridPoints = list(self.grid.getPointsIter())
        msg("... finished grid with %d points" % len(gridPoints))
        if not len(gridPoints):
            raise Exception("No points in the grid.")

        # calculate or look up position of the grid points in the mesh
        boxName = kwargs.get("boxName", None)
        if boxName:
            meshinfopickleName = "%s_coordinates.pickle" % boxName
            try:
                meshPickleFile = open(meshinfopickleName, "rb")
            except IOError:
                # get meshinfo from the odb
                msg("Did not find mesh information pickle file %s."
                    % meshinfopickleName)
                getMeshinfoFromPickle = False
                writeMeshinfoPickle = True
            else:
                getMeshinfoFromPickle = True
                writeMeshinfoPickle = False

        else:
            getMeshinfoFromPickle = False
            writeMeshinfoPickle = False

        # actually read the meshinfo from the pickle file
        if getMeshinfoFromPickle:
            msg("Reading mesh information pickle file %s."%meshinfopickleName)
            msg("Loading point positions for interpolation.")
            self.mesh = pickle.load(meshPickleFile)
            msg("Loaded mesh with %d elements and %d nodes."
                % (len(self.mesh.elNodes), len(self.mesh.nodeCoords)))
            self.elemCoords = pickle.load(meshPickleFile)
            msg("Loaded point positions for %d points." % len(self.elemCoords))
            del meshPickleFile

        # actually calculate position of the grid points in the mesh
        else:
            msg("Reading mesh from the odb.")
            model = self.getAbqModel(recognizedBlocks = "NODE,ELEMENT")
            msg("Finished reading mesh information from the odb.")

            # initializing point search only within the bounding box of points
            model.initializePointSearch(box=self.grid.getBoundingBox())

            # list of (element number, elem coord)-tuples for each pt in gridPoints
            self.elemCoords = model.getElemCoords(gridPoints)
            msg("Calculated element coordinates for %d points, %d points are"
                " within the mesh." % (len(gridPoints), len([
                            None for elem, coords in self.elemCoords
                            if elem is not None])))
            allElems = set(elem for elem, coords in self.elemCoords)
            self.mesh = model.partMeshFromElements(allElems)
            msg("Created a submesh with %d elements and %d nodes."
                % (len(self.mesh.elNodes), len(self.mesh.nodeCoords)))
            del model

        assert len(gridPoints)==len(self.elemCoords)

        # write position of the grid points in the mesh to pickle file
        if writeMeshinfoPickle:
            msg("Writing mesh information pickle file %s."
                % meshinfopickleName)
            meshPickleFile = open(meshinfopickleName, "wb")
            pickle.dump(self.mesh, meshPickleFile)
            msg("Saved submesh.")
            pickle.dump(self.elemCoords, meshPickleFile)
            msg("Saved point positions for %d points." % len(self.elemCoords))
            del meshPickleFile

        # create a filterset for the output data
        msg("Creating filterset for the output data.")
        self.filterElset = self.odbSetFromElements(self.mesh.elNodes.iterkeys())
        self.filterNset = self.odbSetFromNodes(self.mesh.nodeCoords.iterkeys())
        msg("len(filterElset): %d"%len(self.filterElset.elements), debugLevel=10)
        msg("len(filterNset): %d"%len(self.filterNset.nodes), debugLevel=10)

        return self.grid

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
        ptValues = FieldClassOut(field.interpolateToPoints(
                self.mesh, self.elemCoords))
        return ptValues
