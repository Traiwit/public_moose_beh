r"""Module to read field values and other data from an Abaqus output database
(odb).

Provides classes L{OdbReader} for reading nodal or Gauss point data as is and
L{OdbReaderInterpolateToPoints} for reading values interpolated to arbitrary
positions.
"""

### IMPORTANT NOTE:
# The versioning information is stored in bae/odb_04/__init__.py

import os
import re
from struct import Struct, error as StructModuleError
from subprocess import Popen, PIPE
from array import array

import bae.abq_model_02 as abq_model
from bae.field_01 import createFieldClass
from bae.log_01 import msg

from bae.odb_04 import getAbqVersionStr, getAbqExe, getOdbPostData
from bae.odb_04.odbReader_ext import Communicator


#------------------------------------------------------------------------
# odb reader classes

class RemoteError(Exception):
    def __str__(self):
        return ('\n' + '-'*75 + '\n' + str(self.args[0]) + '-'*75)

class OdbNotOpenError(Exception):
    pass

_iStruct = Struct("i")
_fStruct = Struct("f")
_typecodeDict = {int: "i", long: "i", float: "f"}


class OdbReader(object):
    """Class for getting data from the odb.

    Normal operation
    ================

    Example, get min prinipal stress from some specified elements
     >>> odb = OdbReader(odbPath, version=abaqusVersion)
     >>> odb.openOdb()
     >>> mesh = odb.getAbqModel()
     >>> elsetName = odb.createFilterElset(elementLabels)
     >>> fld_S1 = getFieldFromOdb("Step-2", 86, "S_MIN_PRINCIPAL")
     >>> del odb  # also closes the odb
     >>> print "element\tvalues for all integration points"
     >>> for elem in elementLabels:
     >>>     print "%d\t%s" % (elem, fld_S1[elem])


    Implementation details, odbAccessServer
    =======================================

    An odbAccessServer is running in a separate child process that opens the
    odb and does all the communication with the Abaqus API.
    The OdbReader object provides the means for communicating with the
    odbAccessServer to read data from the odb. It then processes this data.

    Simple communication through L{_callmethod}, example:
    -----------------------------------------------------
     >>> elsetName = self._callmethod(
     >>>     "createFilterElset", sorted(filterElset), resultTypeCodes="s")

    Complex communication:
    ----------------------
    ... using L{_startRequest}, L{_finalizeRequest}, L{_pipeWrite},
    L{_pipeFlush}, and L{_pipeRead}. Example:
     >>> # call odbAccessServer-function getMatNames with one argument
     >>> self._startRequest("getMatNames", firstMatCode)
     >>>
     >>> # pass additional data (here for demostration)
     >>> self._pipeWrite(1, "section")
     >>> self._pipeFlush()  # always required after last _pipeWrite
     >>>
     >>> # get results from the pipe (just an example...)
     >>> resultData = self._pipeRead("isss")
     >>>
     >>> # all data received, close transaction
     >>> self._finalizeRequest()

    Communication between _startRequest and _finalizeRequest can as well be
    delegated to a separate task or object like the class
    L{bae.odb_03.odbReader_ext.Communicator}.
    """

    _odbComponentToArrayIndex = {
        "1": 0, "2": 1, "3": 2,
        # tensor components in the odb: xx, yy, zz, xy, xz, yz
        "11": 0, "22": 1, "33": 2,
        "12": 3, "13": 4, "23": 5}

    _odbInvariantNameToNb = {
        "magnitude": 1,
        "mises": 2,
        "tresca": 3,
        "press": 4,
        "inv3": 5,
        "max_principal": 6,
        "mid_principal": 7,
        "min_principal": 8,
        "max_in_plane_principal": 9,
        "min_in_plane_principal": 10,
        "out_of_plane_principal": 11,
        }

    #{ housekeeping
    def __init__(self, odbPath, version=None):
        """
        @param odbPath: file name of the odb.
        @param version: Version string of Abaqus, something like "6.13-2".
           Note: separators must be exactly as in the example: dot and dash.
           For new style use something like "abq2018", ignoring the hot-fix
           postfix e.g. "hf4".
           If not specified defaults to what the 'abaqus' shell command
           invokes.
        """
        msg('OdbReaderBase.__init__: PID %d; PPID: %d'
            % (os.getpid(), os.getppid()), debugLevel=10)

        self.currentCommand = None

        #--- start subprocess with odb C-API: odbAccessServer

        # determine command to execute
        # It's in the same directory as this module.
        exeName = "odbAccessServer_%s" % getAbqVersionStr(version)
        exeDir = os.path.abspath(os.path.dirname(__file__))
        cmd = [getAbqExe(version), os.path.join(exeDir, exeName)]

        # start odbAccessServer
        msg('OdbReaderBase.__init__: starting %s' % cmd,
            debugLevel=10)
        try:
            self.odbAccessServer = Popen(cmd, stdin=PIPE, stdout=PIPE)
        except OSError as exc:
            msg("ERROR: OdbReaderBase.__init__ received an OSError when"
                " launching the odbAccessServer:\n%s\n" % exc.child_traceback)
            raise
        msg('OdbReaderBase.__init__: started odbAccessServer.',
            debugLevel=10)

        # check if server is (still) running
        # this check seems to not work because the child process needs more
        # time to decide whether it's running fine or not. The check would
        # require "import time; time.sleep(0.5)" in advance.
        # ...but we're not that patient...
        if not self.odbAccessServer.poll() is None:
            errMsg = (
                "ERROR: OdbReaderBase.__init__: server died immediately."
                " Its stdout:\n%s\n" % self.odbAccessServer.stdout.read())
            msg(errMsg)
            raise RemoteError(errMsg)

        # call init on server
        # Why don't we do it on odbAccessServer-startup or later in openOdb?
        # 1. We want all odb related stuff in odbAccessServer/odbFunctions.cc.
        # 2. We may eventually need some procedure separate from openOdb.
        # 3. This is a communication test as well: Is the server still running
        #    and happy to process commands?
        try:
            self._callmethod("init")
        except RemoteError as exc:
            res = re.search(
                r'Command line option (".*") may not be used with "analysis"',
                exc.args[0])
            if res:
                errMsg = (
                    "ERROR: OdbReaderBase.__init__: Did not find the"
                    " odbAccessServer executable %s." % res.group(1))
                msg(errMsg)
                raise RemoteError(
                    "%s\n\n... Error message from odbAccessServer:\n%s"
                    % (errMsg, exc.args[0]))
        msg('OdbReaderBase.__init__: called init in odbAccessServer.',
            debugLevel=10)

        # initialize communicator that'll receive bulk data from server
        self.comm = Communicator(
            dataChannelIn=self.odbAccessServer.stdout.fileno(),
            dataChannelOut=self.odbAccessServer.stdin.fileno())

        self.odbPath = odbPath
        self.odbIsOpen = False

    def __del__(self):

        returncode = self.odbAccessServer.poll()
        if returncode is None:
            msg("OdbReaderBase.__del__: Sending #STOP to the odbAccesServer",
                debugLevel=10)
            self.closeOdb()
            self._pipeWrite("#STOP")
            self._pipeFlush()
            self.odbAccessServer.wait()
        else:
            msg("OdbReaderBase.__del__: odbAccesServer has already terminated.",
                debugLevel=10)
    #} end of housekeeping

    #{ communications
    def _pipeWrite(self, *args):
        """Write some values to the pipe.

        Usage:
         >>> self._pipeWrite(5, 7.314, "hallo")

        Arrays, lists and tuples are preceded by the number of items.
        Arrays are written with their array.tofile() method. Lists and tuples
        are converted to a corresponding array and then treated alike. Their
        items must all be of the same type either int or float.

        @param args: all arguments will be written to the pipe

        @note: L{_pipeFlush} must be called after the last _pipeWrite.
        """
        f = self.odbAccessServer.stdin
        for x in args:
            if isinstance(x, basestring):
                f.write(_iStruct.pack(len(x)))
                f.write(x)
            elif isinstance(x, int):
                f.write(_iStruct.pack(x))
            elif isinstance(x, float):
                f.write(_fStruct.pack(x))
            elif isinstance(x, array):
                f.write(_iStruct.pack(len(x)))
                x.tofile(f)
            elif isinstance(x, (tuple, list)):
                f.write(_iStruct.pack(len(x)))
                if len(x):
                    tc = _typecodeDict[type(x[0])]
                    array(tc, x).tofile(f)
            else:
                raise ValueError(
                    "_pipeWrite: Could not write <%s> to pipe." % type(x))

    def _pipeFlush(self):
        """This must be called after all arguments/data have been written
        to the pipe in order to flush the IO buffer. That's because
        self.comm doesn't use these IO buffers!
        """
        self.odbAccessServer.stdin.flush()

    def _pipeRead(self, typecodes):
        """Read some values from the pipe

        Usage:
         >>> (i, x, text) = self._pipeRead("ifs")
         >>>
         >>> # CAUTION: to read just one value do...
         >>> (value, ) = self._pipeRead("f")

        @param typecodes: a string of single characters, s for string,
           i for int, f for float
        """
        f = self.odbAccessServer.stdout
        result = list()
        for tc in typecodes:
            try:
                if tc=="s":
                    length = _iStruct.unpack(f.read(_iStruct.size))[0]
                    result.append(f.read(length))
                elif tc=="i":
                    result.append(_iStruct.unpack(f.read(_iStruct.size))[0])
                elif tc=="f":
                    result.append(_fStruct.unpack(f.read(_fStruct.size))[0])
                else:
                    raise ValueError(
                        "_pipeRead: Don't know what to read from pipe as type"
                        " <%s>." % tc)
            except StructModuleError as exc:
                if "unpack requires a string argument of length" in exc.args[0]:
                    errMsg = (
                        "ERROR: odbAccessServer seems to have stopped"
                        " communication.\n\n"
                        "Could not read data of expected length from it:\n%s"
                        % exc.args[0])
                    msg(errMsg)
                    raise RemoteError(errMsg)
                else:
                    raise

        return result

    def _pipeReadArray(self, typecode, length):
        """Read an array of values from the pipe

        Usage:
         >>> (nbElems,) = self._pipeRead("i")
         >>> elems = set(self._pipeReadArray("i", nbElems)

        @param typecode: to be passed on the the constructor of array.array:
           i for int, f for float
        @param length: of the data array to be read from the pipe.
           Must be provided!

        @Note: This does not read the length of the array, just the raw data.
        Typically the length is supplied through the pipe just ahead of the
        raw data. See the usage example above.
        """
        f = self.odbAccessServer.stdout
        arr = array(typecode)
        try:
            arr.fromfile(f, length)
        except EOFError as exc:
            errMsg = (
                "ERROR: odbAccessServer seems to have stopped"
                " communication. Or at least it didn't provide the expected"
                " number of %d values of type<%s>\n\n%s"
                % (length, typecode, exc.args[0]))
            msg(errMsg)
            raise RemoteError(errMsg)
        return arr

    def _startRequest(self, func, *args):
        """Call a method of the odbAccessServer. The odbAccessServer will be
        processing the request in a separate process asyncronously.

        Communication between the main process -that called _startRequest()-
        and the odbAccessServer is possible through self.odbAccessServer.stdin
        (to push data to the server) and self.odbAccessServer.stdout (to read
        data from the server). See also L{_pipeWrite}, L{_pipeFlush},
        L{_pipeRead}, and L{_pipeReadArray}.

        This function must eventually be followed by self.L{_finalizeRequest}().
        """
        self._pipeWrite(func, *args)
        self._pipeFlush()
        self.currentCommand = func
        return

    def _finalizeRequest(self, resultTypeCodes=""):
        """Receive the result from the function previously started by
        self._startRequest().

        This function must be called after a corresponding call to
        self._startRequest().

        @kwarg resultTypeCodes: A string defining the types of the return
           values. One character for each result item: i for int, f for float,
           s for string.

        @returns: Only if resultTypeCodes is given then return the result
           value(s). If len(resultTypeCodes)==1 then return the corrsponding
           value. If len(resultTypeCodes)>1 then return a list of values.
        """
        assert self.currentCommand
        if resultTypeCodes:
            result = self._pipeRead(resultTypeCodes)

            # make result a single value if only one value is expected
            if len(resultTypeCodes)==1:
                result = result[0]
        else:
            result = None
        (endMarker,) = self._pipeRead("s")
        if endMarker == "#ERROR":
            # error message is generated by the function in odbFuncion.cc
            (errmsg, ) = self._pipeRead("s")
            raise RemoteError(
                "Command <%s> returned an error:\n%s"
                % (self.currentCommand, errmsg))
        elif endMarker == "#EXCEPTION":
            # exception issued by the function in odbFuncion.cc
            # (through the use of pipeException())
            (excTypeStr, errmsg) = self._pipeRead("ss")
            ExcType = eval(excTypeStr)
            raise ExcType(errmsg)
        elif endMarker != "#RETURN":
            errmsg = (
                "Command <%s> did not finalise properly. Received <%s> instead"
                " of <#RETURN>" % (self.currentCommand, endMarker))
            if resultTypeCodes:
                errmsg += " Got this result: %s" % result
            raise RemoteError(errmsg)
        self.currentCommand = None
        return result

    def _callmethod(self, func, *args, **kwargs):
        """Call a method of the odbAccessServer and return the
        result.

        Exxample (in OdbReaderInterpolateToPoints.initOutputPoints):
         >>> elsetName = self._callmethod(
         >>>     "createFilterElset", sorted(filterElset), resultTypeCodes="s")

        @kwarg resultTypeCodes: A string defining the types of the return
           values. One character for each result item: i for int, f for float,
           s for string.

        @returns: Only if resultTypeCodes is given then return the result
           value(s). If len(resultTypeCodes)==1 then return the corrsponding
           value. If len(resultTypeCodes)>1 then return a list of values.
        """
        self._startRequest(func, *args)
        try:
            resultTypeCodes = kwargs["resultTypeCodes"]
        except KeyError:
            resultTypeCodes = ""
        return self._finalizeRequest(resultTypeCodes)
    #} end of communications

    #{ open/close
    def openOdb(self):
        """Open the odb.
        """
        try:
            self._callmethod("openOdb", self.odbPath)
        except RemoteError:
            msg("ERROR: Could not open odb %s" % self.odbPath)
            raise
        msg('OdbReaderBase.__init__: opened odb %s.'
            % self.odbPath, debugLevel=10)
        self.odbIsOpen = True
        return

    def closeOdb(self):
        if self.odbIsOpen:
            self._callmethod("closeOdb")
            self.odbIsOpen = False
        return
    #} end of open/close

    #{ read model data from odb
    def getAbqModel(self, recognizedBlocks=None):
        r"""get an L{abq_model_02.Model<bae.abq_model_02.Model>}-object from
        the odb model data,

        reads nodes, elements and elsets

        assumes there is only one part instance in the model

        ideas for the future:
         - read nsets, surfaces, all the rest

        @param recognizedBlocks: If specified only get the specified property
        from the odb. Not implemented so far. Reads "NODE", "ELEMENT", "ELSET".
        Can be specified as tuple of strings or comma separated string.

        @Note: elsets are being read from the first part instance and from the
        root assembly as well. If the same elset name occurs twice the elset
        in the root assembly is ignored.
        @Note: Look in abq_model_02.Model.read how recognizedBlocks is treated
        there and adapt this here accordingly.
        @Note: Can only handle 3D models: Node coords are always float-triples.
        """
        model = abq_model.Model()

        # preprocess the recognizedBlocks argument
        recognizedBlocks = model._checkRecognizedBlocks(recognizedBlocks)

        if not self.odbIsOpen:
            raise OdbNotOpenError()

        if "NODE" in recognizedBlocks:
            # get the nodes
            self._startRequest("getNodeCoordsFromOdb")
            self.comm.updateNodeCoordsFromPipe(model.nodeCoords)
            self._finalizeRequest()

        if "ELEMENT" in recognizedBlocks:
            # get the elements
            self._startRequest("getElementsFromOdb")
            self.comm.updateElNodesFromPipe(
                model.elNodes, model.elType, model.typeEl)
            self._finalizeRequest()

        if "ELSET" in recognizedBlocks:
            # get the elements
            self._startRequest("getElsetsFromOdb")
            while 1:
                (elsetName,) = self._pipeRead("s")
                if not elsetName:
                    break
                (nbElems,) = self._pipeRead("i")
                elset = set(self._pipeReadArray("i", nbElems))
                model.elset[elsetName] = elset
            self._finalizeRequest()

        msg("OdbReader.getAbqModel: Got model from the odb: %s" % model,
            debugLevel=10)

        return model

    def getMatFieldFromOdb(self,
                           nameFilter=None, firstMatCode=None,
                           matNameToMatNumber=None):
        """Read material section name from the odb at grid point
        positions.

        @param nameFilter: A callable or None for material name conversion.
          For example if there are region variants that shall be treated as
          one. I.e. if SEDIMENT, SEDIMENT_CAVE and SEDIMENT_UCUT shall all be
          treated as SEDIMENT.

          Must be a function otherwise it'll be ignored. It will be called as a
          conversion function taking the original material name from the odb as
          first argument and a dummy section type "Unknown type" as second
          argument. It's assumed to return the desired converted material name.
          (The section type was intended to read something like "solid",
          "beam", "shell" but this functionality is deprecated. The second
          argument is only left there for compatibility and should be ignored.)

          Example: To take everything up to the first "_" or "-" use as
          matRegionsNameFilter:
           >>> lambda matname,secType: matname.split("_")[0].split("-")[0]

          That's equivalent to:
           >>> def convert(matname,secType):
           >>>     return matname.split("_")[0].split("-")[0]
           >>> odbToVtkGridBox(..., matRegionsNameFilter=convert, ...)

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

          If not specified then it defaults to 1 or the next after what's in
          the matNameToMatNumber argument.

        @param matNameToMatNumber: Optional dictionary determinining the
          material name to material number assignment. Keys are the final
          material names as returned by the nameFilter function if applicable.

          Typically the keys are the uppercase material names as specified
          in the Abaqus input file of the job.

          Note: Should a material name appear that is not in the supplied
          dictionary yet it will be added.

          Note: If matNameToMatNumber is specified then the resulting
          matNamesList only contains new material names that were not
          specified in the matNameToMatNumber argument.

          Note: All keys (mat names) are silently converted to uppercase
          if needed. This is done in place, i.e. the dictionary provided
          by the calling function is modified as a side effect.

        @returns: (matNamesList, matNumberField)-tuple
          matNamesList is the list of material names.
          matNumberField is the field (i.e. a list) of material

        @Note: Only materials that are found at any of the given points will
          be listed in the resulting matNamesList. (See optional parameter
          matNameToMatNumber for an exception to this.) Material names from
          other regions of the model are not considered.
        """

        # determine or check firstMatCode
        if firstMatCode is None:
            # if firstMatCode not given then determine from matNameToMatNumber
            if matNameToMatNumber:
                firstMatCode = max(matNameToMatNumber.itervalues()) + 1
            else:
                firstMatCode = 1
        elif (matNameToMatNumber and
                firstMatCode<=max(matNameToMatNumber.itervalues())):
            # if firstMatCode is given and interferes with matNameToMatNumber
            # then issue a warning
            lastName, lastNumber = max(
                matNameToMatNumber.iteritems(), key=lambda x,y: y)
            msg("WARNING"
                " from OdbReader.getMatFieldFromOdb():"
                " firstMatCode %d <= last specified mat number from"
                " matNameToMatNumber %d (%s). Automatically added mat numbers"
                " starting at %d may conflict with the given"
                " matNameToMatNumber dict without notice."
                % (firstMatCode, lastNumber, lastName))
            del lastName, lastNumber

        # check that given mat names are upper case
        if matNameToMatNumber:
            for matName, matNb in matNameToMatNumber.items():
                nameUpper = matName.upper()
                if nameUpper != matName:
                    matNameToMatNumber[nameUpper] = matNameToMatNumber[matName]
                    del matNameToMatNumber[matName]

        if not self.odbIsOpen:
            raise OdbNotOpenError()

        # get the actual data from the odb
        msg("OdbReader.getMatFieldFromOdb: starting getMatNames on the"
            " odbAccessServer.", debugLevel=10)
        self._startRequest("getMatNames", firstMatCode, 0)
        matNumbersPrelim = dict()
        self.comm.updateFieldValuesConstInt(matNumbersPrelim)
        msg("OdbReader.getMatFieldFromOdb: Got matNumbersPrelim for %d pts."
            % len(matNumbersPrelim), debugLevel=10)

        # get matNamesList from the pipe
        (nbMatNames,) = self._pipeRead("i")
        msg("OdbReader.getMatFieldFromOdb: will read %d matNames from the"
            " pipe." % nbMatNames, debugLevel=10)
        matNamesList = list(self._pipeRead("s"*nbMatNames))
        self._finalizeRequest()
        msg("OdbReader.getMatFieldFromOdb: read %d mat numbers and %d mat"
            " names."
            % (len(matNumbersPrelim), len(matNamesList)), debugLevel=10)

        # extract material region info from section names and apply nameFilter
        nextNewNumber = firstMatCode
        oldToNewMatNumber = dict()
        newMatNamesList = list()
        if matNameToMatNumber is None:
            matNameToMatNumber = dict()
        for oldNumber, matName in enumerate(matNamesList, start=firstMatCode):
            if callable(nameFilter):
                matName = nameFilter(matName, "Unknown type")
            try:
                newNumber = matNameToMatNumber[matName]
            except KeyError:
                newNumber = nextNewNumber
                nextNewNumber += 1
                matNameToMatNumber[matName] = newNumber
                newMatNamesList.append(matName)
            oldToNewMatNumber[oldNumber] = newNumber
        matNamesList = newMatNamesList
        del newMatNamesList

        # create appropriate bae.field_01.Field object converting matNumbers
        MatNumberFieldClass = createFieldClass("Material", "element", "scalar")
        try:
            matNumberField = MatNumberFieldClass(
                (e, oldToNewMatNumber[x])
                for (e, x) in matNumbersPrelim.iteritems())
        except KeyError:
            raise KeyError(
                "ERROR when translating material numbers: Not all mat numbers"
                " %s have a corresponding entry in the list of mat names. This"
                " list has %d entries and starts with %d."
                % (sorted(set(matNumbersPrelim.itervalues())),
                   len(oldToNewMatNumber), firstMatCode))

        msg("Identified %d different new material codes for %d elements."
            " (Not counting the %d preassigned material codes.)"
            % (len(matNamesList), len(matNumberField),
               len(matNameToMatNumber)-len(matNamesList)))
        return (matNamesList, matNumberField)
    #} end of read model data from odb

    #{ read fielddata from odb
    def getFieldProps(self, fieldName):
        """Determine odb field name, invariant and component from the
        fieldName.

        @returns: a 6-tuple of the following values:
         - fieldName: may be stripped of a leading "SDV_"
         - position: i.e. node, elemIP, ...
         - dataType: i.e. vector, scalar, tensor
         - odbFieldName: string, key in the FieldValues repository of the odb
         - component: integer, -1 if not a vector or tensor-component
         - invariant: integer according to the L{_odbInvariantNameToNb} dict
        """
        invariant = 0
        component = -1

        if fieldName.startswith("SDV"):
            odbFieldName = fieldName
        elif '_' in fieldName:
            odbFieldName, refinement = fieldName.split('_', 1)
            if refinement.isdigit():
                # component of a vector or tensor
                try:
                    component = self._odbComponentToArrayIndex[refinement]
                except KeyError:
                    raise ValueError(
                        "Invalid component %s, must be one of %s"
                        % (refinement, ",".join(
                            self._odbComponentToArrayIndex)))
            else:
                # invariant (like magnitude or so)
                try:
                    invariant = self._odbInvariantNameToNb[refinement.lower()]
                except KeyError:
                    raise NotImplementedError(
                        "OdbReaderBase.getFieldProps: Invariant %s not known."
                        % refinement)
        else:
            odbFieldName = fieldName
        msg("OdbReaderBase.getFieldProps:"
            " Field %s: name %s, component %s, invariant %s"
            % (fieldName, odbFieldName, component, invariant), debugLevel=10)

        if not self.odbIsOpen:
            raise OdbNotOpenError()

        # get position (i.e. node, elemIP, ...; do we need it?) and type
        # (like vector, scalar, tensor) from the odb
        (position, dataType) = self._callmethod(
            "getFieldProps", odbFieldName, resultTypeCodes="ss")
        if invariant or component>=0:
            dataType = "scalar"
        msg("OdbReaderBase.getFieldProps: ... position %s, dataType %s."
            % (position, dataType), debugLevel=10)

        # remove SDV_ from fieldName from now on
        if fieldName.startswith("SDV_"):
            fieldName = fieldName[4:]

        return fieldName, position, dataType, odbFieldName, component, invariant

    def createFilterElset(self, elset=None):
        """Set a elset as filter for all subsequent read operations for
        element based fields.
        @param elset: set or other iterable of element labels (ints) Supply an
           empty list to disable the filter.
        @returns: the Abaqus-API internal set name. This set name can later be
           used for L{setFilterElset} to re-activate this filter set.
        """
        if elset:
            elset = sorted(elset)
        else:
            elset = []

        if not self.odbIsOpen:
            raise OdbNotOpenError()

        elsetName = self._callmethod(
            "createFilterElset", elset, resultTypeCodes="s")
        return elsetName

    def createFilterNset(self, nset=None):
        """Set a nset as filter for all subsequent read operations for
        nodal fields.
        @param nset: set or other iterable of node labels (ints). Supply an
           empty list to disable the filter.
        @returns: the Abaqus-API internal set name. This set name can later be
           used for L{setFilterNset} to re-activate this filter set.
        """
        if nset:
            nset = sorted(nset)
        else:
            nset = []

        if not self.odbIsOpen:
            raise OdbNotOpenError()

        nsetName = self._callmethod(
            "createFilterNset", sorted(nset), resultTypeCodes="s")
        return nsetName

    def setFilterElset(self, elsetName=""):
        """Set a elset as filter for all subsequent read operations for
        element based fields.
        @param elsetName: name of an elset already in the odb or created
          earlier. Supply an empty string to disable the filter.
        @returns: nb of elements in the filter elset; 0 if the filter is
          disabled; -1 if the specified elset name does not exist.
        """
        if not elsetName:
            elsetName = ""
        elif not isinstance(elsetName, basestring):
            raise ValueError("setFilterElset expects a string as elsetName"
                             " argument but got %s." % elsetName)

        if not self.odbIsOpen:
            raise OdbNotOpenError()

        nbItems = self._callmethod(
            "setFilterElset", elsetName, resultTypeCodes="i")
        return nbItems

    def setFilterNset(self, nsetName=""):
        """Set a nset as filter for all subsequent read operations for
        nodal fields.
        @param nsetName: name of an nset already in the odb or created
          earlier. Supply an empty string to disable the filter.
        @returns: nb of nodes in the filter nset; 0 if the filter is
          disabled; -1 if the specified nset name does not exist.
        """
        if not nsetName:
            nsetName = ""
        elif not isinstance(nsetName, basestring):
            raise ValueError("setFilterNset expects a string as nset name but"
                             " got %s." % nsetName)

        if not self.odbIsOpen:
            raise OdbNotOpenError()

        nbItems = self._callmethod(
            "setFilterNset", nsetName, resultTypeCodes="i")
        return nbItems

    def getOdbStepNames(self):
        """Return a list of odb step names.
        """

        if not self.odbIsOpen:
            raise OdbNotOpenError()

        self._startRequest("getOdbStepNames")

        (nbOdbSteps,) = self._pipeRead("i")
        msg("OdbReader.getOdbStepNames: will read %d step names from the"
            " pipe." % nbOdbSteps, debugLevel=10)
        odbStepNames = list(self._pipeRead("s"*nbOdbSteps))

        self._finalizeRequest()
        msg("OdbReader.getOdbStepNames: read %d step names. ending."
            % len(odbStepNames), debugLevel=10)

        return odbStepNames

    def getOdbFrames(self, odbStep):
        """Return the number of frames in the specified odb step.
        @param odbStep: (string) name of the odb step, e.g. "Step-1"
        @returns: the number of frames in this odb step, an integer
        """
        framesCnt = self._callmethod(
            "getOdbFrames", odbStep, resultTypeCodes="i")
        return framesCnt

    def getFieldFromOdb(self, odbStep, odbFrame, fieldName):
        """
        Example:
         >>> odb = OdbReader("myodb.odb")
         >>> odb.setFilterElset("ALL_TETS")
         >>> fld = odb.getFieldFromOdb(2, 11, "S")
         >>> print fld.intPtSectPt, fld.dataType, fld.position
         [(1,-1), (2,-1), (3,-1), (4,-1)] tensor elemIP
         >>> print fld[101234]  # get stress tensors for elem 101234
         [[1.434E6, 3.123E4, ... 6 components for Gauss pt 1...], [6 comps
         ... for Gauss pt 2], [...], [...]]
         >>> odb.setFilterElset("BEAMS")
         >>> fld = odb.getFieldFromOdb(2, 11, "S_11")
         >>> print fld.intPtSectPt, fld.dataType, fld.position
         [(1,3), (1,7), (1,11), (1,15)] scalar elemIP
         >>> print fld[9012]  # get stress component 11 for elem 9012
         [3.1415E6, 3.333E5, ... 4 components for 4 section points]

        @param odbStep: odb step name (a string) of the step in the odb.
           Something like "Step-3".
        @param odbFrame: integer frame number
        @param fieldName: examples: "S", "V", "S_MISES", "S_MIN_PRINCIPAL",
           "U_1", "S_33", "SDV4", "SDV_PST"

        @returns: A L{Field<bae.field_01.Field>} object.

           If possible the field has an attribute intPtSectPt: a list of (Gauss
           pt number, section point number) tuples. The items in this list
           correspond to the items in each value of the returned field, see
           examples provided.

           This attribute is missing if elements of different types with
           different numbers of integration points and section points are being
           requested.

           If you need this attribute and don't get it yet then use
           L{setFilterElset} or L{createFilterElset} to constrain the set of
           elements returned by a specific call to this function.

        @raises ValueError: if the specified frame is not in the odb.

        @Note: For tensor fields of structural elements the full tensor might
           not be available and cause an error. In this case take individual
           components one by one. The reason is that Abaqus has different
           types of tensors (TENSOR_3D_FULL, TENSOR_3D_SURFACE,
           TENSOR_3D_PLANAR) that are not fully implemented in this module.
        @Note: For beam elements (B31 in A/Explicit) it's suggested to only
           consider S_11. Gero doesn't understand the rest.
        @Note: For shell elements it might be necessary to query the
           components as S_11, S_22, S_33 in order to get S_11, S_22, S_12
           due to limitations of this module.
        """
        # Calls two different methods on the odbAccessServer to
        # accomplish the task of getting field data from the odb.

        # First it calls L{getFieldProps} to determine some properties of the
        # field. Then it initializes the resulting field values object (a
        # L{bae.field_01.Field} object) and then calls the getFieldData.

        if not self.odbIsOpen:
            raise OdbNotOpenError()

        # set frame, can raise ValueError
        assert isinstance(odbStep, basestring)
        assert isinstance(odbFrame, int)
        self._callmethod("setFrame", odbStep, odbFrame)

        # field, invariant and component will be derived from the name
        (fieldName, position, dataType, odbFieldName, component, invariant)\
            = self.getFieldProps(fieldName)

        # create appropriate bae.field_01.Field
        FieldClass = createFieldClass(fieldName, position, dataType)
        field = FieldClass()

        # get the actual data from the odb
        self._startRequest("getFieldData", odbFieldName, component, invariant)
        self.comm.updateFieldValues(position, field)
        self._finalizeRequest()
        msg("OdbReader.getFieldFromOdb: read %d field values. ending."
            % len(field), debugLevel=10)
        return field
    #} end of  read fielddata from odb

    #{ other
    @property
    def postData(self):
        """An object with attributes being fed from the corresponding
        ..._postData.py file (old) or the _postData folder (new).

        Usually contains the attribute "sequence" with some seqElsets...
        dictionaries like seqElsetsExcav in the exampe below. And a frameNames
        dictionary. And the attribute "material" with some material properties
        data.

        Example: Suppose in frame 3 of step 2 labelled "YR2012_M02" the elsets
        "Pit_2012_02" and "Dev_2012_02" will be excavated and in frame 4
        "YR2012_M03" only "Dev_2012_03":
         >>> odb = OdbReader("/datapth/myodb.odb")
         >>> frameNames = odb.postData.sequence.frameNames
         >>> print frameNames[2][3:5]
         ... ['YR2012_M02','YR2012_M03']
         >>> seqSets = odb.postData.sequence.seqElsetsExcav
         >>> print seqSets[2][3:5]
         ... [ ('Pit_2012_02','Dev_2012_02'),(Dev_2012_03,) ]

        The data will only be read when the portData property is being accessed
        for the first time.

        Raises IOError if the XXX_postData.py file can not be accessed. XXX is
        the name of the odb without the a trailing ".odb"-suffix.
        """
        if not hasattr(self, "_postData"):
            if self.odbPath:
                self._postData = getOdbPostData(self.odbPath)
            else:
                self._postData = dict()
                

        return self._postData
    #} end of other
