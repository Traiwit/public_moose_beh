"""Module to handle odbs

first of all: read data easily

All scripts using this module must be run with "abaqus python myscript.py" or
"abaqus cae noGUI=myscript.py".
"""

import odbAccess
from warnings import warn
from sys import stdout
import copy

from bae.future_01 import *
from bae.misc_01 import MsgTicker, QuietLog, Container, BoundingBox
import bae.abq_model_02 as abq_model
from bae.volume_02 import VolumeFromMesh

_debug = False
# _debug = True

class OdbReadWarn(UserWarning):
    pass

def _getOdbLookUpInfo(fieldName):
    """derive some info, how to look up the given field in the odb
    for internal use only

    returns a tuple containing:
     - FieldOutput.name,
     - a boolean: True if field shall be taken from an invariant
     - component name or attribute name of invariant,
    examples:
     - "S_MISES"         -> "S", True, "mises"
     - "S_MIN_PRINCIPAL" -> "S", True, "minPrincipal"
     - "U_1"             -> "U", False, "1"
     - "SDV4"            -> "SDV4", False, None
    """
    if '_' in fieldName:
        field, component = fieldName.split('_', 1)
        isInvariant = not(component.isdigit())
        if isInvariant:
            component = component.lower()
            try:
                while 1:
                    i = component.index('_')
                    component = (component[:i] + component[i+1].upper()
                                 + component[i+2:])
            except ValueError:
                pass
        else:
            component = field + component
    else:
        field = fieldName
        isInvariant = False
        component = None
    
#     print 'field: ', field, type(field)
    return field, isInvariant, component



class MapOdbType(dict):
    """dict {odb field type : "scalar" / "vector"/ "tensor"}

    there is one instance defined in the module: mapOdbType
    """
    def __init__(self):
        dict.__init__(self, {
                odbAccess.SCALAR : "scalar",
                odbAccess.VECTOR : "vector",
                odbAccess.TENSOR_3D_SURFACE : "vector",
                odbAccess.TENSOR_3D_FULL : "tensor",
                })
        def __getitem__(self, key):
            try:
                return dict.__getitem__(self, key)
            except KeyError:
                raise ValueError("odb field type %s not implemented yet."%key)
mapOdbType = MapOdbType()


class OutputField(dict):
    """
    Class for those objects that hold the data for one field at one frame.

    They are basically dicts {node/element number: value} with the extra
    attributes position, type and componentLabels

     >>> field = OutputFieldXXXX()
    
    field.position :
       - "node"    : field[node number] = nodal value
       - "element" : field[element n.] = element value
                     e.g. mean of integration point values
       - "elist"   : field[element n.] = list of values
                     one value for each integration point

    field.type specifies the type of the values:
       - "scalar" : a single float
       - "vector" : array of three floats
       - "tensor" : array of six floats, components: xx, yy, zz, xy, xz, yz

    field.componentLabels contains the data from the corresponding attribute of
    FieldOutput.values[...] in the odb. It's a list specifying the components
    in case field.type is "vector" or "tensor".

    @Bug:
    Only checks the precision attribute of the first value. If those attributes
    differ from value to value we are in trouble. (That's for speed.)
    """

    __slots__ = ["position", "type", "componentLabels"]  # additional data

    def __init__(self):
        pass
    


#     def updateFromOdb(self, fieldName, treatIP="mean", considerElems=None):
#         """
#         @param treatIP: only relevant for data at element integration points
#         like stress, strain. May be:
#          - 'list': save a list of all individual values for all IP.
#            self.position="elist" in this case.
#          - 'x2node': extrapolate integration point values to the nodes
#            averaging over the contribution from different elements
#            self.position="node" (not implemented yet)
#          - 'mean' take the arithmetic mean of all IP values for this element
#            self.position="element" (this is the default)

#         @param considerElems: a set (or so) specifying all elements to consider
#           if "elem not in considerElems" then elem is ignored. Used for
#           treatIP = 'x2node'
#         """


class OdbReader(object):
    """Class to read and store data from the odb

    member attributes filled by self.readModel():
     - nodeCoords : {node number: [coordinates (undeformed)]}
     - elNodes    : {element number: [node numbers]}
     - elType     : {element number: element type string}
     - typeEl     : {element type string: set of element numbers}

    frameData member attribute (a dict {frame:{fieldname:OutputField object}})
    is filled by self.readData():

    >>> field = self.frameData[frameidx][fieldname]  # OutputField object

    field.position :
       - "node"    : field[node number] = nodal value (float or array of
         floats)
       - "element" : field[element n.] = element value (float or array of
         floats), e.g. mean of integration point values
       - "elist"   : field[element n.] = list of values, one value
                     one float for each integration point
    values in the list above are floats in case of scalar quantities and arrays
    (as in the module array, they are no lists) of floats in case of vector or
    tensor quantities.

    You may as well access the Abaqus odb object as self.odb (as long as the
    self.close() function has not been called).
    """

    def __init__(self, odbName, logfile = stdout):
        """Constructor

        odbName is the file name of the odb to process.

        logfile ... already opened log file for diagnostic output
        logfile = None ... to suppress diagnostic output
        logfile not given ... write diagnostic output to stdout
        """

        self.odb = odbAccess.openOdb(path=odbName, readOnly=True)

        if logfile == None:
            self.logfile = QuietLog()
        else:
            self.logfile = logfile

        # permanent member attributes - model data
        self.nodeCoords = dict()
        self.elNodes = dict()
        self.elType = dict()
        self.typeEl = defaultdict(set)

        #-- copied from abq_model.Model, not used so far
        # self.elset = internal.ElsetsDict()
        # self.nset = internal.NsetsDict()
        # self.boundary = internal.BoundaryDict()
        # self.nodeBC = internal.NodeBoundaryDict(self)
        # self.surface = internal.SurfacesDict()

        # permanent member attributes - output data
        self.frameData = dict()

        # initialize default start frame
        self.setFrameZero(0,0)

    def close(self):
        self.odb.close()


    def setFrameZero(self, step, frame):
        """specify the odb step and frame that will be treated as frame with
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
                raise ValueError("There is no %dth step in the odb." % firstStepNo)
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

#             self.nbFrame=nb_frames
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


    def getFrameIterator(self, frameIds=None):
        """iterator over all frames, basically for internal use
        usage:
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
          positive value on output.
        """

        if frameIds==None:
            frameIds = range(len(self.frameNoToStepFrame))

        oldstepName = None
        for frameId in frameIds:
            # convert negativ values, needed for the output
            while frameId<0:
                frameId+=len(self.frameNoToStepFrame)
            # look up the corresponding step, frame in the odb
            try:
                stepName, frameNo = self.frameNoToStepFrame[frameId]
            except IndexError:
                self.logfile.write(
                    "WARNING: OdbReader.getFrameIterator() got an invalid"
                    " frame id %d in its argument frameIds. Must be below %d.\n"
                    % (frameId, len(self.frameNoToStepFrame)))
                raise StopIteration()
            if stepName != oldstepName:
                currentStep = self.odb.steps[stepName]
                oldstepName = stepName
            yield frameId, currentStep.frames[frameNo]


    def getNumberOfFrames(self):
        """returns the number of frames self.getFrameIterator() yields
        """
        return len(self.frameNoToStepFrame)


    def asAbqModel(self, logfile="UseTheSameLogFile"):
        """get an abq_model.Model-object from the odb model data
        you have to invoke self.readModel() first."""
        if logfile == "UseTheSameLogFile":
            logfile=self.logfile
        model = abq_model.Model(logfile=logfile)

        # copy the data
        model.nodeCoords = abq_model.internal.NodesDict(self.nodeCoords)
        model.elNodes = abq_model.internal.ElementNodesDict(self.elNodes)
        model.elType = dict(self.elType)
        model.typeEl = copy.deepcopy(self.typeEl)
        
        return model

    def getAbqModel(self, logfile="UseTheSameLogFile"):
        """Deprecated alias for asAbqModel(), only kept for compatibility.
        
        Don't use this anymore. Use asAbqModel()!
        """
        return self.asAbqModel(logfile=logfile)


    def readModel(self):
        """read nodes and elements

        assumes there is only one part instance in the model

        overwrites self.nodeCoords

        ideas for the future:
         - add the possibility to selectively read only elements, nodes (see
           recognizedBlocks argument of abq_model_02.Model.read() / write())
         - read elsets, nsets, surfaces, all the rest
        always use same names as abq_model_02
        """

        # mesh data from the first instance
        part_instance = self.odb.rootAssembly.instances.values()[0]

        # find dimensionality of our problem
        if part_instance.embeddedSpace == odbAccess.TWO_D_PLANAR:
            self.dim = 2
        elif part_instance.embeddedSpace == odbAccess.THREE_D:
            self.dim = 3
        else:
            msg = "Can't handle dimensionality %s."%part_instance.embeddedSpace
            self.logfile.write(msg)
            raise Exception(msg)

        # extract undeformed node coordinates
        self.nodeCoords = dict()
        msgTemplate = ("Read %s/%d nodes from the odb.\n"
                       % ("%d", len(part_instance.nodes)))
        ticker = MsgTicker(self.logfile, msgTemplate, printLast=True)
        for cnt, node in enumerate(part_instance.nodes):
            self.nodeCoords[node.label] = node.coordinates
            ticker.msg(cnt+1)
        del ticker

        # extract element/cell connectivity
        self.elNodes = dict()
        self.elType = dict()
        self.typeEl = defaultdict(set)
        msgTemplate = ("Read %s/%d elements from the odb.\n"
                       % ("%d", len(part_instance.elements)))
        ticker = MsgTicker(self.logfile, msgTemplate, printLast=True)
        for cnt, elt in enumerate(part_instance.elements):
            self.elNodes[elt.label] = list(elt.connectivity)
            self.elType[elt.label] = str(elt.type)
            self.typeEl[elt.type].add(elt.label)
            ticker.msg(cnt+1)
        del ticker


    def getFields(self, frameIdx=-1):
        """which fields are there in the odb in the specified frame
        (does not read the data itself, use self.readData() for that

        frameIdx is the consecutive frame number used for self.frameData based
        on the frame zero specified with self.setFrameZero().

        yields a dict {fieldname: additional data}, where additional data is a
        Container object with the following attributes:
         - description (a string)
         - type ("scalar" / "vector" / "tensor")
         - componentLabels (list of strings)
         - validInvariants (list of strings)
        """
        stepName, frameNo = self.frameNoToStepFrame[frameIdx]
        frame = self.odb.steps[stepName].frames[frameNo]
        fields = frame.fieldOutputs

        fieldDict = dict([
                (str(f.name), Container(
                        description = str(f.description),
                        type = mapOdbType[f.type],
                        componentLabels = list(f.componentLabels),
                        validInvariants = [str(x) for x in f.validInvariants]
                        ))
                for f in fields.values()])
        return fieldDict


    def filterElements(self, elset):
        """
        specify elements that are to be recognized by readData exclusively.

        @param elset: If None take all elements. Otherwise only those contained
          in this argument.
        """
        # future: also accept elset names... if isinstance(elset, basestring):

        if elset==None:
            self.filterSet = None
        else:
            elset = list(elset)
            elset.sort()

            # mesh data is in the first instance
            part_instance = self.odb.rootAssembly.instances.values()[0]

            # create the new odb elset
            self.filterSet = part_instance.ElementSetFromElementLabels(
                name="currentFilterSet", elementLabels=elset)


    def filterNodes(self, nset):
        """
        specify nodes that are to be recognized by readData exclusively.

        @param nset: If None take all nodes and elements. Otherwise only those
          contained in this argument.
        """
        if nset==None:
            self.filterSet = None
        else:
            nset = list(nset)
            nset.sort()

            # mesh data is in the first instance
            part_instance = self.odb.rootAssembly.instances.values()[0]

            # create the new odb nset
            self.filterSet = part_instance.NodeSetFromNodeLabels(
                name="currentFilterSet", nodeLabels=nset)


    def readData(self, fields, frameIds=None):
        """adds data to the frameData-attribute

        @param fields: is a list of field names like "U", "S_MISES", "S_11" or
        just a single field name. It may also be a list of tuples [field name,
        integration point treatment]. The components of such a tuple are:
         - field name : E.g. "SDV3", "PEEQ"
         - integration point treatment: 
           - None if not an integration point value
           - "mean" to take the mean value
           - "list" yields a list of the values at all integration points

        @param frameIds: may be
         - None or not specified: read all frames starting with the one
           specified by self.setFrameZero()
         - an int: read frame with that consecutive frame number (as used for
           self.frameData based on the frame zero specified with
           self.setFrameZero() )
         - a list of consecutive frame numbers to read
        
        @note:
        If the field name doesn't contain a "_" the value as in the
        FieldValue.data attribute will be taken. E.g. "U", "SDV4", "S"
        If the fieldname contains a "_" with a number following, the
        corresponting component of the field will be taken. E.g. "U_1", "S_11"
        If the fieldname contains a "_" with something else following, this
        something else will be treated as an invariant identifier and the value
        of this invariant will be taken. E.g. "S_MISES", "S_PRESS",
        "U_MAGNITUDE", "LE_MIN_PRINCIPAL"

        """

        if _debug:
            print '\n1) In readData...'

        # interpret frameIds argument
        framesText = ""
        if frameIds==None:
            nbOfFrames = self.getNumberOfFrames()
            framesText = "all "
        else:
            if type(frameIds)==int:
                frameIds = [frameIds]

            frameIds = [f for f in frameIds if f<self.getNumberOfFrames()]
            nbOfFrames = len(frameIds)

        if nbOfFrames==0:
            self.logfile.write(
                "WARNING: The request to read data contains no valid frame id."
                " Nothing read.\n")
            return

        # interpret fields argument in case it's just a single name
        if isinstance(fields, basestring):
            fields = [fields]

        # just the message...
        if nbOfFrames==1:
            framesText = "frame %d" % frameIds[0]
        else:
            framesText += "%d frames" % nbOfFrames
        self.logfile.write("Will read %s with %d field(s)...\n"
                           %(framesText, len(fields)))
        del framesText

        # get the iterator for the main loop
        frameIterator = self.getFrameIterator(frameIds)

        # make fields a list of tuples (field name, treatIP option)
        newFields = list()
        for field in fields:
            if type(field) == str:
                newFields.append((field, "mean"))
            else:
                newFields.append(field)
        fields = newFields

        # here the actual data is read
        if _debug:
            print '\n2) Start reading data...'

        msgTemplate = ("Read output data for frame %d (%d/"
                       +str(nbOfFrames)
                       +"), currently field %s.\n")
        ticker = MsgTicker(self.logfile, msgTemplate=msgTemplate)
        for frameCnt, (frameNo, frame) in enumerate(frameIterator):
            odbFields = frame.fieldOutputs
            if _debug:
                print '   Copy Abaqus frame.fieldOutputs repository'
                print 'type(frame): ', type(frame)
            for fieldName, treatIP in fields:
                # Creation of object (dictionary):
                field = OutputField()
                if _debug:
                    print 'Updating from odb...'

                odbFieldName, isInvariant, component = _getOdbLookUpInfo(fieldName)

                if _debug:
                    print 'odbFieldName: ', odbFieldName
                    print 'isInvariant: ', isInvariant
                    print 'component: ', component, type(component)

                    print '\n\nodbFieldName: ', odbFieldName
                # odbField
                try:
                    odbField = odbFields[odbFieldName]
                except KeyError:
                    warn(OdbReadWarn(
                            "No field %s (%s) in frame %d."
                            % (odbFieldName, fieldName, frameNo)))
                    return

                # odbType, self.type
                if component==None:
                    odbType = odbField.type
                else:
                    odbType = odbAccess.SCALAR
                field.type = mapOdbType[odbType]

                # set field.componentLabels
                if (component==None and odbType!=odbAccess.SCALAR):
                    field.componentLabels = list(odbField.componentLabels)
                else:
                    field.componentLabels = None

                # componentIndex
                if component!=None and not(isInvariant):
                    try:
                        if _debug:
                            print 'odbField.componentLabels: ', odbField.componentLabels
                        componentIndex = list(odbField.componentLabels).index(component)
                        component = None
                    except ValueError:
                        raise ValueError("No component %s in field %s."%(component, fieldName))
                else:
                    componentIndex = None

                # component (if None so far)
                if component==None:
                    # check precision attribute of the first value
                    precision = odbField.values[0].precision
                    if precision==odbAccess.SINGLE_PRECISION:
                        component = "data"
                    elif precision==odbAccess.DOUBLE_PRECISION:
                        component = "dataDouble"
                    else:
                        raise ValueError("Precision %s not implemented."
                                         %precision)

                # odbPosition, self.position
                positions = set([f_loc.position
                                 for f_loc in odbField.locations])
                if len(positions)==1:
                    odbPosition = positions.pop()
                else:
                    odbPosition = positions.pop()
                    warn(OdbReadWarn(
                            "Multiple positions for field %s. All but"
                            " position %s will be ignored."
                            % (odbFieldName, odbPosition)))

                # filter filterSet if applicable
                if self.filterSet!=None:
                    odbField = odbField.getSubset(region=self.filterSet)

                # get the actual data from the odb
                if odbPosition==odbAccess.NODAL:
                    field.position = "node"
                    if componentIndex==None:
                        field.update([(value.nodeLabel, getattr(value, component))
                                     for value in odbField.values])
                    else:
                        field.update([(value.nodeLabel,
                                      getattr(value, component)[componentIndex])
                                     for value in odbField.values])
                elif odbPosition==odbAccess.WHOLE_ELEMENT:
                    field.position = "element"
                    assert componentIndex==None
                    field.update([(value.elementLabel, getattr(value, component))
                                 for value in odbField.values])
                elif (odbPosition==odbAccess.INTEGRATION_POINT
                      and treatIP=="list"):
                    field.position = "elist"
                    newValues = defaultdict(list)
                    for value in odbField.values:
                        newValues[value.elementLabel].append(getattr(value, component))
                    field.update(newValues)
                    del newValues
                elif (odbPosition==odbAccess.INTEGRATION_POINT
                      and treatIP=="x2node"):
                    field.position = "node"
                    raise ValueError("Integration point treatment %s not"
                                     " implemented yet."%treatIP)
                elif (odbPosition==odbAccess.INTEGRATION_POINT
                      and treatIP=="mean"):
                    field.position = "element"
                    newValues = defaultdict(list)
                    for value in odbField.values:
                        newValues[value.elementLabel].append(getattr(value, component))
                    if field.type == "scalar":
                        for label, values in newValues.iteritems():
                            field[label] = sum(values)/len(values)
                    else:
                        for label, values in newValues.iteritems():
                            valuesIter = iter(values)
                            avgVal = list(valuesIter.next())
                            for newVal in valuesIter:
                                for i, item in enumerate(newVal):
                                    avgVal[i]+=item
                            for i in range(len(avgVal)):
                                avgVal[i] /= len(values)
                            field[label] = avgVal
                    del newValues
                elif odbPosition==odbAccess.INTEGRATION_POINT:
                    raise ValueError("Integration point treatment %s not"
                                     " implemented yet."%treatIP)
                else:
                    raise ValueError("odbPosition %s not implemented yet."%odbPosition)

                # store the field just read in self.frameData[frameNo]
                try:
                    thisFrameData = self.frameData[frameNo]
                except KeyError:
                    thisFrameData = dict()
                    self.frameData[frameNo] = thisFrameData
                thisFrameData[fieldName] = field
                ticker.msg(frameNo, frameCnt+1, fieldName)
        ##-- end of loop
        del ticker


    def getFieldValueAtPoints(self, pointCoords, fields, frameIds=None):
        """
        Get field values at specified arbitrary points (not nodes).

        This function does not use the cae path feature but does a linear
        interpolation of nodal values.

        @param pointCoords: a list of point coordinates
        @param fields: see self.readData()
        @param frameIds: may be
         - None or not specified: read all frames starting with the one
           specified by self.setFrameZero()
         - an int: read frame with that consecutive frame number (as used for
           self.frameData based on the frame zero specified with
           self.setFrameZero() )
         - a list of consecutive frame numbers to read

        @returns: A dictionary: {frame:{fieldname:list of values}} The list of
        values contains one value for each given point.

        @Note: So far this only works for nodal fields like U.

        @Note: The function does not check, if the necessary data from the odb
        has already been read previously. It just reads everything it needs
        (again).
        """

        #--- find elements and their nodes to each point

        # make sure the model data has been read
        if len(self.nodeCoords)==0:
            # assume that the model data has not been read so far
            self.readModel()

        # find min, max for all point coordinates
        box = BoundingBox()
        box.update(pointCoords)
        self.logfile.write("Requested values for %d points in the following"
                           " box:\n%s\n" % (len(pointCoords), box))

        # create a volume instance to perform the element lookup
        model = self.asAbqModel()
        vol = VolumeFromMesh(model)
        vol.initializePointSearchInBox(box=box)

        # findig elements and shape function values for each point
        # elemShapeFunc = [vol.pointIsInside(point, returnLinShapeFunc=True)
        #                  for point in pointCoords]
        elemShapeFunc = list()
        msgTemplate = ("Finding elements for points (%s/%d done)\n"
                       % ("%d", len(pointCoords)))
        ticker = MsgTicker(self.logfile, msgTemplate=msgTemplate)
        for cnt, point in enumerate(pointCoords):
            elemShapeFunc.append(
                vol.pointIsInside(point, returnLinShapeFunc=True))
            ticker.msg(cnt+1)
        del ticker

        # nodes of those elements
        allElems = set((pt[0] for pt in elemShapeFunc if pt != False))
        allNodes = model.getNodesFromElset(elements=allElems)

        #--- read nodal data from the odb
        self.filterNodes(allNodes)
        self.readData(fields=fields, frameIds=frameIds)
        frameIds = self.frameData.keys()
        frameIds.sort()
        fields = self.frameData[frameIds[0]].keys()

        #--- interpolate to points and store results
        result = dict()
        ticker = MsgTicker(self.logfile, msgTemplate=(
                "Interpolating nodal values to points. In frame %s (%s/%d)"
                " %s/%d fields done so far.\n"
                % ("%d", "%d", len(frameIds), "%d", len(fields))))
        
        for cntFr, frameNo in enumerate(frameIds):
            fieldDict = self.frameData[frameNo]
            frameData = dict()
            for cntFi, (fieldName, field) in enumerate(fieldDict.iteritems()):
                assert field.position == "node", (
                    "Only nodal values can be read by getFieldValueAtPoints()"
                    " so far. %s is not suitable (%s)."
                    % (fieldName, field.position))
#                 pointsField = list()
#                 for elem, shapeFunc in elemShapeFunc:
#                     nodes = model.elNodes[elem]
                    
#                     if field.type == "scalar":
#                         ptValue = sum(weight*field[node]
#                                       for weight, node in zip(shapeFunc,nodes))
#                     else:
#                         ptValue = [sum(weight*component
#                                        for weight, node in zip(shapeFunc,nodes))
#                                    for component in field[node]]
#                     pointsField.append(ptValue)

                if field.type == "scalar":
                    pointsField = []
                    for res in elemShapeFunc:
                        if res != False:
                            elem, shapeFunc = res
                            pointsField.append(sum(
                                weight*field[node]
                                for weight,node
                                in zip(shapeFunc, model.elNodes[elem])))
                        else:
                            pointsField.append(0.0)
                else: # type = vector or tensor
                    pointsField = []
                    for res in elemShapeFunc:
                        if res != False:
                            elem, shapeFunc = res
                            weightValList = [
                                (weight, field[node])
                                for weight, node
                                in zip(shapeFunc, model.elNodes[elem])]
                            pointsField.append([
                                    sum(weight*val[i] for weight, val in weightValList)
                                    for i in xrange(len(weightValList[0][1]))])
                        else:
                            if field.type == "vector":
                                pointsField.append([0.0, 0.0, 0.0])
                            else: # field.type == "tensor"
                                pointsField.append([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

                frameData[fieldName] = pointsField
                ticker.msg(frameNo, cntFr+1, cntFi+1)

            result[frameNo] = frameData
        del ticker

        #--- fini
        return result
