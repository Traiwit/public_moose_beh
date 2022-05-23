r"""Module provides the base class for data exchange with the odb.
"""

__version__ = "3.08"

_version_history_ = r"""\
Versions:

3.08 GP new, split off from odbReader.py v3.08
"""

import os
import re
from struct import Struct, error as StructModuleError
from subprocess import Popen, PIPE
from array import array

from bae.log_01 import msg

from bae.odb_03 import getAbqVersionStr, getAbqExe
from bae.odb_03.odbReader_ext import Communicator


class RemoteError(Exception):
    def __str__(self):
        return ('\n' + '-'*75 + '\n' + str(self.args[0]) + '-'*75)

_iStruct = Struct("i")
_fStruct = Struct("f")
_typecodeDict = {int: "i", long: "i", float: "f"}

class OdbReaderBase(object):
    """Base class for getting data from the odb.

    Launches a separate process with odbAccessServer and provides the means for
    communication with it.

    Derived classes provide convenient methods that access the odb through the
    odbAccessServer. This base class provides basic routines to communicate
    with the odbAccessServer.


    Simple communication through L{_callmethod}, example:
    =====================================================
     >>> elsetName = self._callmethod(
     >>>     "createFilterElset", sorted(filterElset), resultTypeCodes="s")


    Complex communication through L{_startRequest}, L{_finalizeRequest},
    L{_pipeWrite}, L{_pipeFlush}, and L{_pipeRead}, example:
    ========================================================
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
        "12": 3, "13": 4, "23": 5 }

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

    def __init__(self, odbPath, version=None, openOdb=True):
        """
        @param odbPath: file name of the odb.
        @param version: Version string of Abaqus, something like "6.13-2".
           Note: separators must be exactly as in the example: dot and dash.
           For new style use something like "abq2018", ignoring the hot-fix
           postfix e.g. "hf4".
           If not specified defaults to what the 'abaqus' shell command
           invokes.
        @param openOdb: If True (the default) then automatically open the odb.
           If set to False then the odb must be opened manually with
           L{openOdb}.
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
        except OSError, exc:
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
        except RemoteError, exc:
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

        #--- open odb
        if openOdb:
            self.openOdb()

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
            except StructModuleError, exc:
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
        except EOFError, exc:
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
