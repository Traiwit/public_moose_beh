r"""Module for easy event logging.

Usage:
======

In the main script do:
 >>> from log_01 import log, msg
 >>> msg("Hi there, let's go.")  # will write to stdout
 >>> log.close()   # might be needed to automatically write the end message

(msg being just a shortcut for log.msg, which is obviously the most important
method of log)

In modules to be used by the main script do
 >>> # module myworker_01.py
 >>> from log_01 import log, msg
 >>> def myworker():
 >>>    msg("Now this is a message from the worker.")

In the main script you may assign output channels. You may optionally assign
to only a specific range of debuglevel values:
 >>> from log_01 import log, msg
 >>> from myworker_01 import myworker
 >>> log.toFile("auto") # standard messages go to logfile with a generic name
 >>> log.addFile("myscript.someDebug.log", maxDebugLevel=3)
 >>> log.addFile("myscript.verbose.log", maxDebugLevel=None)
 >>> msg("normal output") # will go to all files
 >>> msg("debug 3", debugLevel=3) # will go to last two files
 >>> msg("debug 5", debugLevel=5) # will only go to last file
 >>> myworker()  # will also send output to all three files
 >>> log.close()   # might be needed to properly close the files

To have all output suppressed do:
 >>> from log_01 import log, msg
 >>> from myworker_01 import myworker
 >>> log.quiet()
 >>> msg("Anybody listening?")   # this message goes nowhere.
 >>> myworker()  # output of subroutines will also be suppressed


There is an initialization string that's being written when the first message
is issued. And an end string which is being written when the module is deleted
(or the LogManager instance is deleted or the log.close() method is being
called).

@Note: Ordinary messages for ordinary users have debugLevel=0. Messages with
debugLevel>0 are debugging messages that you don't want to see in normal
operation. The higher the debugLevel the more sophisticated the output and
the less people want to see it.
"""

__version__ = "1.12"

_version_history = r"""\
Versions:
=========

1.0 GP: first version
1.1 GP: added the "date" special for the LogManager.addFile() logfile argument.
1.2 GP: added log.setDebugLevel() method.
1.3 AF: the log.channels[] list of (write,flush,close,...)-tuples is changed
        to a list of _LogChannelFile-instances.
        also, some formatting enhancements in LogManager.msg().
1.4 GP: fixed:
   - default debugLevel=0, otherwise msg without debugLevel arg won't work
   - changed order of arguments in setDebugLevel
   - corrected docs
   - toFile("auto") did not work because of assignment to arg[0]
1.5 GP: removed writePattern (replaced by Arnd's linePrefixes)
        added some Rhino specific features, cleaned up a bit.
1.6 GP: clarified and fixed: msg with debugLevel>0 are debugging messages that
        are to be suppressed in normal cases.
1.7 GP added name and file attributes to _LogChannelFile
       added getLogFileForSubProcess
1.8 GP changed default mode for logfile="auto" and "date" to "a" (append).
1.9 GP rename _LogChannelFile to LogChannelFile making it visible in the docs
1.10 GP changed logfile="auto", "date" to write to the current working dir
1.11 GP added with log.sleeping(): capability
1.12 GP added LogManager.fmtTime and msgWrite()-function for raw output
"""

__bugs__ = r"""\
toFile("auto") option does not work with "abaqus view noGUI=..."
see an implementation in bae.cae_01
"""

import sys, os, time
import atexit
from getpass import getuser
from socket import gethostname
from contextlib import contextmanager

# Using abaqus2017/18 platform.system() results in a fatal (uncatchable) error.
# As a hofix ossysytemname now bypasses the platform property from sys-module
# which results in a 'similar' platform string (e.g. win64 instead of Windows)  
#from platform import system as ossystemname
def ossystemname():
    return sys.platform

try:
    import rhinoscriptsyntax as rhino
except ImportError:
    rhino = None
    


class LogChannelFile(object):
    """
    A file object for diagnostic output.
    Automatically writes some header and footer upon initialization and
    deletion.

    Some attributes of the corresponding file like object (i.e. the log file)
    are accessible through this LogChannelFile-object. This is done by simply
    assigning a reference of the file attribute to a corresponding attribute
    of self: E.g. self.write = logfile.write.

    Here is the list of those attributes: write, close, flush, name.

    @ivar file: is a reference to the underling file object itself
    @ivar linewidth: number of chars in one line, see corresponding argument
       of self.__init__
    @ivar doubleSep: separator line string of "=" characters
    @ivar singleSep: separator line string of "-" characters
    @ivar minDebugLevel: lower limit to output-messages, i.e. a message
       will be written to this particular channel only if its debugLevel
       argument is higher or equal to this value.
    @ivar maxDebugLevel: upper limit to output-messages, i.e. a message
       will be written to this particular channel only if its debugLevel
       argument is lower or equal to this value.
    @ivar skipInitStrings: if true then don't write header and footer strings.

    @Note: minDebugLevel and maxDebugLevel are only stored in this object.
       It's up to the LogManager to decide whether to write a message or skip
       it. Most likely this will be done based on the values of these
       attributes.

    @Note: Usually an instance of this class should write a final message when
       it's being deleted. However this is not guaranteed (and appears to
       usually not happen) if the object still exists when the python
       interpreter exits. So it is most save to call self.close(). (del self
       is not that sensible because there might still be another object holding
       a reference to self.
    """

    def __init__(self, logfile="stdout", mode="w",
                 linewidth=None, minDebugLevel=-999, maxDebugLevel=0):
        """
        Constructor, can be used like open()

        @param logfile: If this equals "stdout" then write to sys.stdout.
            If this argument is an already opened file object, write to this.
            Otherwise opens a file with the specified name.

        @param mode: like the corresponding argument to open(). Default "w",
            i.e. open for write / overwrite. Use "a" to append to an existing
            file.

        @param linewidth: number of chars in one line (when used with linebreak
           option in msg-call). Defaults to 79 on Windows and else 80. (It is
           useful to restrict it to 79 for output in Windows 'cmd'-window.)

        @param minDebugLevel: lower limit to output-messages, i.e. a message
           will be written to this particular channel only if its debugLevel
           argument is higher or equal to this value.
           Defaults to -999 (no sensible restriction)

        @param maxDebugLevel: upper limit to output-messages, i.e. a message
           will be written to this particular channel only if its debugLevel
           argument is lower or equal to this value.
           Defaults to 0. (Print only ordinary messages with debugLevel<=0,
           no debugging information with debugLevel>0.)

        @Note: The logfile argument has already been processed by the calling
           LogManager.addFile() function. The special values "auto" and "date"
           are treated at that place.

        @Note: The debug-level limits for a channel restrict output of
           msg() calls that have the debugLevel within these limits. Generally,
           the larger the debug number, the more debuggy the output information
           gets.
           see parameter debugLevel, or alias dbg, in L{LogManager.msg()}
           method for more details.
        """

        ### linewidth and separator-lines
        ###----------------------------------------
        if linewidth is None:
            if sys.platform=="win32":
                self.linewidth = 79
            else:
                self.linewidth = 80
        else:
            self.linewidth = linewidth
        self.doubleSep = "="*self.linewidth+"\n"
        self.singleSep = "-"*self.linewidth+"\n"

        ### debugging limits
        ###----------------------------------------
        self.minDebugLevel = minDebugLevel
        self.maxDebugLevel = maxDebugLevel

        ### setup logfile
        ###----------------------------------------
        # dummy close or flush method
        def doNothing():
            pass

        ### case 'stdout': use stdout for output
        if logfile=="stdout":
            logfile=sys.stdout
            self.close = doNothing

        ### case 'other string than stdout': create a new logfile with the
        ### given string as filename
        elif isinstance(logfile, basestring):
            logfile = open(logfile, mode)
            try:
                self.close = logfile.close
            except AttributeError:
                self.close = doNothing

        ### assume already opened file-like object if it's got a write method
        ### Note: existing files will not be closed
        elif hasattr(logfile, "write"):
            self.close = doNothing

        ### don't know what to do with this logfile type
        else:
            raise ValueError("logfile argument is neither a string nor a file.")

        ### flush method
        try:
            self.flush = logfile.flush
        except AttributeError:
            self.flush = doNothing

        ### name attribute
        try:
            self.name = logfile.name
        except AttributeError:
            pass

        ### file attribute
        self.file = logfile

        ### write method
        ### the write method for Rhino stdout is more complex
        ### It is defined separately as class methods, see below
        if not (rhino and logfile==sys.stdout):
            self.write = logfile.write

        ### some log channels don't want initialization strings
        ### (they just don't work right)
        if rhino and logfile==sys.stdout:
            self.skipInitStrings = True
        else:
            self.skipInitStrings = False


    if rhino:
        # special write method for Rhino writing to stdout
        # otherwise write is not defined here or it's overlaid by instance
        # variables defined in __init__()
        def write(self, txt):
            sys.stdout.write(txt)
            rhino.Prompt(txt.strip())

    ### ADDITIONALLY, however, I still might consider defining a msg-method
    ### within this class. (because of reformatting the text-message for
    ### different linewidths of various channel-instances)
    #def msg(...):
    #    pass
    ### NO, I DO NOT THINK THIS WILL BE REQUIRED, THOUGH.


class LogManager(object):
    """
    class to manage different log channels. Provides the msg function
    to submit messages.

    @ivar channels: list of L{LogChannelFile} instances
    @ivar fmtTime: Format of the time prefix for each message. Defaults to
       "%Hh%M.%S". See time.strftime for how to specify the format.

    @bug: scriptName variable in __init__ not suitable for use with abaqus
       cae/viewer. Need to detect this case and then get the correct name.
       Use cae_01.openDefaultLogFile() to extract the correct script name.
    """

    def __init__(self):

        ### default list of LogChannelFile instances
        self.channels = [
            LogChannelFile(),  # use all default values
            ]

        ### used for sleeping context manager and sleep/wakeup methods
        # "sleep Debug Level Offset"
        # normally = 0, will be set to 1 by self.sleep(). This offset will
        # be added to the message debugLevel before comparing to
        # LogChannelFile.maxDebugLevel
        self.sleepDLO = 0

        ### extract name of script
        ### Note: this might be not sophisitcated enough, check with versions
        ###       where main script is called
        ###       - python mainscript.py        -> ??
        ###       - python thisscript.py        -> ??
        ###       - abq python mainscript.py    -> ??
        ###       - abq python thisscript.py    -> ??
        ###       - abq vi noGUI=mainscript.py    -> ???

        if sys.argv[0].endswith("ABQvwrK.exe"):
            try:
                scriptNamePos = sys.argv.index('noGUI') + 1
            except ValueError:
                # there is no noGUI argument
                self.scriptName = "<unknown>"
            else:
                if len(sys.argv)<=scriptNamePos:
                    # noGUI has been the last argument
                    self.scriptName = "<unknown>"
                else:
                    # next argument is the script name
                    self.scriptName = sys.argv[scriptNamePos]
                    if self.scriptName=="-":
                        self.scriptName = "<unknown>"
        else:
            self.scriptName = sys.argv[0]

        self.scriptFullname = os.path.abspath(self.scriptName)
        tmpNames = self.scriptFullname.replace(os.sep,"/").rsplit("/",1)
        assert(len(tmpNames)==2)
        self.scriptPathname = tmpNames[0].replace("/",os.sep)+os.sep
        self.scriptPathCWD  = os.getcwd()+os.sep
        self.scriptFilename = tmpNames[1]
        self.scriptFileextn = self.scriptFilename.rsplit(".",1)[-1]
        del tmpNames

        self.scriptArgs = []
        if len(sys.argv)>1:
            self.scriptArgs = sys.argv[1:]

        ### provide a string for summarizing script-execution information
        initString = [
            "script '%s' executed"%(self.scriptFilename),
            "",
            "  by %s on %s (%s)"%(getuser(),gethostname(),ossystemname()),
            "  in %s"%(os.getcwd()+os.sep),
            "  on %s"%(time.strftime("%a, %Y-%m-%d at %H:%M:%S")),
            "",
            ]
        if self.scriptArgs:
            initString.extend([
                "and with these commandline arguments:",
                "",
                ]+[
                "  sys.argv[%d] = '%s'"%(i+1,str(strArg))
                for i,strArg in enumerate(self.scriptArgs)
                ]+[
                "",
                ])
        else:
            initString.extend([
                "with no additional commandline arguments.",
                ])
        self.writeInitString = "\n".join(initString)+"\n"

        # default time format in message prefix
        self.fmtTime = "%Hh%M.%S"

        # only write end string if at least one message has been issued
        self.writeFiniString = False

        # initialize linePrefixes
        self.linePrefixes = []
        self.linePrefix = ""

    def quiet(self):
        """
        Quiet out this log.

        This will just send every subsequent message into the black hole
        at the centre of the galaxy whose name was also swallowed by that very
        black hole.

        @Note: after usage of quiet(), the channels cannot be recovered.
        Use "with log.L{sleeping}():" otherwise. Or in special circumstances
        use L{sleep}(), L{wakeup}().
        """
        self.channels = []

    def sleep(self):
        """Generally: DON'T USE THIS METHOD. Use "with log.L{sleeping}():"
        instead.

        Quiet out this log (all channels). Recovery of channels is possible
        with L{wakeup}().

        See also L{sleeping} -generally the preferred method. sleep() and
        wakeup() may be used if they are not being called in the same
        subroutine. However the following example is an example of very bad
        coding style:
         >>> from bae.log_01 import log, msg
         >>>
         >>> def somethingWithWakeup()
         >>>     do_something()  # possible messages still swallowed
         >>>     log.wakeup()
         >>>     maybe_do_some_more()
         >>>     msg("This message should pop up.")
         >>>     return
         >>>
         >>> msg("This is important.")  # will be logged
         >>> log.sleep()
         >>> do_something_with_nasty_messages()  # messages suppressed
         >>> somethingWithWakeup()
         >>> msg("This message should be logged as well.")
        """
        self.sleepDLO = 1

    def wakeup(self):
        """
        Recovery of LogChannels that were put to sleep before. See L{sleep}().

        Also consider using "with log.L{sleeping}():"
        """
        self.sleepDLO = 0

    @contextmanager
    def sleeping(self):
        """Return a U{context manager
        <https://docs.python.org/2.7/library/stdtypes.html#typecontextmanager>}
        that makes the log channel suppress all output. to be used with the
        with statement:
         >>> from bae.log_01 import log, msg
         >>> log.addFile("auto")
         >>> msg("This will be logged.")
         >>> with log.sleeping():
         >>>     msg("This message will not be printed or logged anywhere.")
         >>>     do_something_with_a_nasty_message()
         >>>     msg("No one will see the nasty message (neither this one).")
         >>> msg("Here we are back again. This message will be logged.")

        See also L{quiet}() to permanently silence the log. And L{sleep}() and
        L{wakeup}().
        """
        self.sleep()
        try:
            yield
        finally:
            self.wakeup()

    def toFile(self, *args, **kwargs):
        """
        Switches output to the specified file. This will be the only output
        channel from now on, as long as you don't add more channels.

        Same arguments as L{addFile}.
        """
        self.quiet()
        self.addFile(*args, **kwargs)

    def addFile(self, *args, **kwargs):
        """
        Add the specified file as an additional output channel.

        The arguments logfile and mode might as well be specified as the first
        and second positional arguments. This is to make the interface similar
        to that of the builtin open().

        @kwarg logfile: There are three possible special strings:
          - "stdout" ... write to sys.stdout
          - "auto"   ... open and write to a file named as the current script
                         with extension '*.log' in the current working directory
          - "date"   ... as "auto" but additionally appends a date and time
                         stamp to the file name.
                         format: 'scriptFilename_YYYYMMDD-hhmm.log'

        If this argument is an already opened file object, write to this.
        Otherwise opens a file with the specified name.

        @kwarg mode: like the corresponding argument to open(). Default "w",
        i.e. open for write / overwrite. Use "a" to append to an existing file.
        Default mode for logfile="auto" and "date" is "a" (append)

        @kwarg maxDebugLevel: Optional. Only write messages with a debugLevel
           smaller than or equal to this. None means write all messages.
           Defaults to 0 (only ordinary diagnostic output, no debugging
           information).

        @kwarg minDebugLevel: Optional and usually not specified. Only write
           messages with a debugLevel higher than or equal to this. This might
           be used if only debugging messages should go to a dedicated file and
           not the ordinary diagnostic output. I doubt that it makes much sense
           to use this... None means not restriction. Defaults to None.

        @Note: Any more than two positional arguments are silently ignored!
        """

        ### process positional arguments
        ###-----------------------------
        if len(args)>0:
            kwargs["logfile"] = args[0]
        if len(args)>1:
            kwargs["mode"] = args[1]

        ### process argument logfile
        ###-----------------------------
        try:
            logfile = kwargs["logfile"]
        except KeyError:
            raise ValueError("No logfile argument given!")

        ### process various logfile-options
        ###---------------------------------
        if logfile=="auto" or logfile=="date":

            # remove .py extension
            scriptBaseName = self.scriptName.rsplit(".py",1)[0]

            # remove directory prefix
            scriptBaseName = os.path.basename(scriptBaseName)

            # possibly add date stamp
            if logfile=="date":
                scriptBaseName = "%s_%s"%(scriptBaseName,
                                          time.strftime("%Y%m%d-%H%M"))

            # add .log extension
            logfile = scriptBaseName+".log"

            # store logfile name in arguments
            kwargs["logfile"] = logfile

            # default mode for "auto" and "date" is "a" (append)
            if "mode" not in kwargs:
                kwargs["mode"] = "a"

        ### add a LogChannelFile instance
        ###----------------------------
        self.channels.append(LogChannelFile(**kwargs))

    def msgWrite(self, text, debugLevel=0):
        """Write the given text without any prefix or postfix (i.e. end of
        line) to the log file(s).
        This is essentially the same function as L{msg} but doesn't append the
        time stamp nor the end of line.

        @param text: message to be written, will be converted to string as
           str() does

        @param debugLevel: Only write this message to file(s) if this argument
           is smaller or equal than the maxDebugLevel argument of the
           corresponding file and if this argument is larger or equal to the
           minDebugLevel argument of the corresponding file. See
           L{LogManager.addFile} for details how to specify minDebugLevel and
           maxDebugLevel.

           Don't use debugLevels greater than 1Mio.
        """
        for ch in self.channels:
            if debugLevel>=ch.minDebugLevel and \
               debugLevel + self.sleepDLO <= ch.maxDebugLevel:
                ch.write(text)
                ch.flush()

    ### Arndt would like to include:
    ###   - addThisPrefix
    ###   - remLastPrefix
    ###   - linebreak
    ### Geros opinion:
    ###  - go ahead (but don't break it)
    ###  - For the addThisPrefix, remLastPrefix I suggest to add methods to
    ###    the LogManager class:
    ###     . log.addPrefix(" "*4)  # indent all following messages
    ###     . log.delPrefix()  # removes the last prefix
    ###     . delPrefix can have an optional number argument to remove more
    ###       than just the last prefix
    ###     . actually Arnd already added a linePrefix string and a
    ###       linePrefixes list as attributes to the log object.
    ###  - I think: linebreak is superceded by the end argument?
    ###  - fmtTime functionality is already implemented by the instance var
    ###    fmtTime
    def msg(self, text, end="\n", debugLevel=0):
        r"""Write messages to the log file(s).

        The command
         >>> msg(message)

        is basically doing the following (The actual implementation is
        different and additionally there is this debugLevel feature, see
        below):
         >>> for logfile in log.channels:
         >>>     logfile.write("%s : %s\n" %(time.strftime("%Hh%M.%S"),message))
         >>>     logfile.flush()

        @param text: message to be written, will be converted to string as
           str() does
        @param end: append this string to the text. Defaults to newline so text
           does not need to have newline included. If you don't want to finish
           the current line then specify end="".

           I.e. the following two commands result in two lines of output.
            >>> msg("This and that happens.")
            >>> msg("Now something else happens.")

           Whereas the following code will write the first message then do
           something and append the finishing message to the very same line.
            >>> msg("Something is going to happen ...", end="")
            >>> do_something()
            >>> msgWrite(" finished. Something happened.\n")

           In this case you normally don't get the end time. Whether this is a
           desirable feature remains an open question...

        @param debugLevel: Only write this message to file(s) if this argument
           is smaller or equal than the maxDebugLevel argument of the
           corresponding file and if this argument is larger or equal to the
           minDebugLevel argument of the corresponding file. See
           L{LogManager.addFile} for details how to specify minDebugLevel and
           maxDebugLevel.

           Don't use debugLevels greater than 1Mio.

        @Note: Ideas for future parameters (see comment in source code):
        addThisPrefix, remLastPrefix
        """

        # @param addThisPrefix: (future parameter yet to be implemented)
        #    a string that is appended to the indents that
        #    preceed the actual message

        # @param remLastPrefix: (future parameter yet to be implemented)
        #    a boolean wheather to remove the last prefix
        #    from the indents list after the message has been written

        # @param linebreak: (future parameter yet to be implemented)
        #    a boolean wheather to wrap the message in the next
        #    lines when exceeding the channels linewidth. Note the indent
        #    preceeding the message will be repeated with blanks for all next
        #    lines.

        if self.writeInitString is not False:
            # write initialization string
            for ch in self.channels:
                if (ch.minDebugLevel<=0 and ch.maxDebugLevel>=0
                        and not ch.skipInitStrings):
                    ch.write("%s%s%s" % (
                        ch.doubleSep, self.writeInitString, ch.doubleSep))
                    ch.flush()
            # write writeInitString only once
            self.writeInitString = False
            # only write end string if at least one message has been issued
            self.writeFiniString = True

        msgString = "%s : %s%s" % (time.strftime(self.fmtTime), text, end)
        msgWrite(msgString, debugLevel)

    def close(self):
        """Closes all output channels

        @Note: Usually an instance of this class should write a final message
           when it's being deleted. However this is not guaranteed (and appears
           to usually not happen) if the object still exists when the python
           interpreter exits. So it is most save to call self.close(). (del
           self is not that sensible because there might still be another
           object holding a reference self.
        """

        try:
            endString = ("script '%s' finished %s\non %s (with ...)\n"
                         %(self.scriptFilename,"",
                           time.strftime("%a, %Y-%m-%d at %H:%M:%S")
                         ))
        except AttributeError:
            # this usually happens when the interpreter exits before self.close
            # has been called and the order of object destruction is determined
            # by the python interpreter. The time module might not be there
            # anymore.
            endString = "Fini.\n"

        for ch in self.channels:
            if (ch.minDebugLevel<=0 and ch.maxDebugLevel>=0
                and self.writeFiniString and not ch.skipInitStrings):
                ch.write("\n%s%s%s" % (
                    ch.doubleSep, endString, ch.doubleSep))
                ch.flush()
            # close all channels
            ch.close()

        self.writeFiniString = False

    def __del__(self):
        self.close()

    def setWritePattern(self, pattern=None):
        """
        Method is **deprecated**. keep for compatibilty reasons only.
        Does not do anything!
        """
        pass

    def setDebugLevel(self, maxDebugLevel=None, minDebugLevel=None):
        """
        Change the min- and/or maxDebugLevel attributes of **all** log channels.

        @kwarg maxDebugLevel: Optional. Only write messages with a debugLevel
           smaller than or equal to this. None means write all messages.
           Defaults to 0 (only ordinary diagnostic output, no debugging
           information).

        @kwarg minDebugLevel: Optional and usually not specified. Only write
           messages with a debugLevel higher than or equal to this. This might
           be used if only debugging messages should go to a dedicated file and
           not the ordinary diagnostic output. I doubt that it makes much sense
           to use this... None means not restriction. Defaults to None.
        """
        for ch in self.channels:
            if minDebugLevel is not None:
                ch.minDebugLevel = minDebugLevel
            if maxDebugLevel is not None:
                ch.maxDebugLevel = maxDebugLevel

    def getLogFileForSubProcess(self):
        """Return a log file that a subprocess can send diagnostic output
        to. Choose the first channel that accepts messages with debugLevel=0.
        If no such channel can be found then return sys.stdout.

        Note that the output sent to this file will not get the time stamp
        prefix of messages passed through the msg() function.

        Example:
         >>> from subprocess import call, STDOUT
         >>> from bae.log_01 import log, msg
         >>> call("myprogram",
         >>>      stderr=STDOUT, stdout=log.getLogFileForSubProcess())
        """
        res = None
        for ch in self.channels:
            if ch.maxDebugLevel>=0 and ch.minDebugLevel<0:
                res = ch.file
                break
        if not res:
            res = sys.stdout
        return res


###----------------------------------------------------------------------------
### instantiate the global log and msg
###----------------------------------------------------------------------------
log = LogManager()
msg = log.msg
msgWrite = log.msgWrite

# this makes sure that the end string is being written before the file objects
# are being deleted
atexit.register(log.close)


class MsgTicker(object):
    r"""
    Takes frequent messages (e.g. from a loop) as input and passes them on
    to std LogManager instance for diagnostic output only at a specified
    intervall.

    Usage:
     >>> from bae.log_01 import log, msg, MsgTicker
     >>> jobList = [...]
     >>> msg("Starting to process %d jobs" % len(jobList))
     >>> msgTemplate = "processed %d/"+str(len(jobList))
     >>> ticker = MsgTicker(msgTemplate)
     >>> for cnt, job in enumerate(jobList):
     ...     job.process() # this does the real work
     ...     ticker.msg(cnt+1)
     >>> del ticker

    The del command is necessary to enforce the writing of the last message. If
    you leave this away, the last message will be created when the ticker
    object is destroyed and this could lead to very confusing output.

    You may also change the message template during the lifetime of the ticker
    object. Simply modify the self.msgTemplate attribute. The following example
    creates a special last message. Note the extra call to ticker.msg() with
    arguments appropriate for the new msgTemplate (none in this case).
     >>> msgTemplate = "processing job %d of "+str(len(jobList))
     >>> ticker = MsgTicker(msgTemplate)
     >>> for cnt, job in enumerate(jobList):
     ...     ticker.msg(cnt+1)
     ...     job.process() # this does the real work
     >>> ticker.msgTemplate = "Finished the last job."
     >>> ticker.msg()
     >>> del ticker

    @ivar msgTemplate: message template. May be changed at any time. Mind that
      the arguments of subsequent self.msg() calls must match the string
      conversion items in the new msgTemplate.

    @ivar end: append this string to each output. See the corresponding
      argument of self.__init__(). May be changed at any time.

    @ivar debugLevel: Only write the messages to file(s) with applicable
      maxDebugLevel and minDebugLevel values. See arguments of
      L{LogManager.msg} for details. May be changed at any time. Don't use
      debugLevels greater than 1Mio.

    @cvar defaultFirstPeriod: Default value for the corresponding argument to
      the constructor. Can be changed to save you from supplying your value
      as argument at each invocation.
    @cvar defaultPeriod: Default value for the corresponding argument to
      the constructor. Can be changed to save you from supplying your value
      as argument at each invocation.

    @note: Of course this object can only print something when its msg()
      method is called. If it takes longer between two calls to msg() than the
      specified period for output, output will be delayed accordingly.

    @note: If you need fancier functionality consider subclassing this. It's
      easy and elegant.

    @note: The ticker tries to output as close to the requested interval as
      possible instead of just "period" seconds after the last output: Suppose
      you specify firstPeriod=10 and period=10, so you expect output at roughly
      10s intervals. Now suppose the msg() function is called at times
      11, 19, 35, 41 then there will be output printed at 11, 35, 41. Note
      firstly that there is no output in the interval between 20 and 30. And
      secondly that there is output at 41 although the last output at 35
      occurred only less than 10s ago. That's because the ticker tries to
      write as close to (of course only after) the times 10, 20, 30, 40.
    """
    # you can change those default values once instead of supplying your value
    # as argument at each invocation
    # Change defaultFirstPeriod to Zero for the use in test functions to make
    # sure the message string is evaluated at least once. (Test functions
    # usually take so much time...)
    defaultFirstPeriod=5
    defaultPeriod=60

    def __init__(self, msgTemplate="%s", end="\n", debugLevel=0,
                 firstPeriod=None, period=None,
                 printLast=False):
        """
        @param msgTemplate: message template with %- conversion specifiers
          matching the arguments to the subsequent self.msg() function calls.

        @param end: append this string to the text. Defaults to newline so the
           msgTemplate does not need to have newline included. If you don't
           want to finish the current line then specify end="".

        @param debugLevel: Only write the messages to file(s) with applicable
           maxDebugLevel and minDebugLevel values. See arguments of
           L{LogManager.msg} for details. Don't use debugLevels greater than
           1Mio.

        @param firstPeriod: nb of seconds until the first message is passed.
          If not specified or None, MsgTicker.defaultFirstPeriod will be used.

        @param period: nb of seconds until each subsequent message is passed
          through.
          If not specified or None, MsgTicker.defaultPeriod will be used.
          If Zero each msg will be displayed.

        @param printLast: if True, (at least) the last message is printed in
          any case. This is done when either self.close() is called or the
          object (self) is deleted.

        @Note: If at least one message has already been printed (i.e. the job
          being reported by this MsgTicker lasted longer than firstPeriod),
          the very last message will be printed as well regardless of the
          printLast argument. Setting this argument to True will make sure that
          at least one message (the last) is being printed regardless how long
          the MsgTicker object lives, i.e. how fast the job is done. The
          default value False for the printLast argument will silence the
          MsgTicker if it doesn't last at least for firstPeriod. (More
          precisely: if the last call to its msg() function occurs within
          this period.)
        """
        self.msgTemplate = msgTemplate
        self.end = end
        self.debugLevel = debugLevel
        self.printLast = printLast

        if firstPeriod is None:
            firstPeriod = MsgTicker.defaultFirstPeriod
        self.nextOutputTime = time.time()+firstPeriod

        if period is None:
            self.period = MsgTicker.defaultPeriod
        else:
            self.period = period

        self.lastArgs = None
        self.cnt = 0  # for the self.tick() method

    def msg(self, *args, **kwargs):
        """
        Specify a set of state parameters that may be used to generate
        output if enough time has passed since the last output.

        There are two possible ways to call this function depending on the
        msgTemplate argument to __init__():

        If the msgTemplate argument includes %(name)-style format options,
        the arguments to this msg() function must be keyword arguments:

        >>> ticker=MsgTicker("finished %(cnt)d items. State: %(state)s")
        ...    ...
        ...    ticker.msg(cnt=1, state="go")

        If the msgTemplate argument includes format options without names,
        the arguments to this msg() function must be positional arguments:

        >>> ticker=MsgTicker("finished %d items. State: %s")
        ...    ...
        ...    ticker.msg(1, "go")
        """
        if len(args)==0:
            args = kwargs
        timetowait = self.nextOutputTime - time.time()
        if timetowait>0:
            self.lastArgs = args
        else:
            if self.period!=0:
                steps = int(-timetowait/self.period) + 1
                self.nextOutputTime += self.period*steps

            msg(self.msgTemplate % args,
                end=self.end, debugLevel=self.debugLevel)
            self.lastArgs = None
            self.printLast = True

    def tick(self):
        """
        Like self.msg() except it generates consecutive integers
        (1, 2, 3, 4, 5, ...) as arguments for the message template each call.

        Usage:
         >>> items = ...
         >>> ticker = MsgTicker("processed %%d/%d"%len(items))
         >>> for item in items:
         >>>     ticker.tick()
         >>>     dosomething(item)
         >>> del ticker

        would produce something like::
         processed 1/100
         processed 13/100
         ...
         processed 100/100

        @returns: the current integer counter
        """
        self.cnt += 1
        timetowait = self.nextOutputTime - time.time()
        if timetowait>0:
            self.lastArgs = self.cnt
            return self.cnt
        else:
            if self.period!=0:
                steps = int(-timetowait/self.period) + 1
                self.nextOutputTime += self.period*steps

            msg(self.msgTemplate % self.cnt,
                end=self.end, debugLevel=self.debugLevel)
            self.lastArgs = None
            self.printLast = True
        return self.cnt

    def periodPassed(self):
        """
        Like self.msg() but without doing anything with messages. Just
        returns True each time the specified period has passed once more.

        Intended mainly for use in specialized subclasses so you don't have
        to fiddle around with time comparison and so on.
        """
        if self.nextOutputTime <= time.time():
            self.nextOutputTime += self.period
            return True
        else:
            return False

    def close(self):
        """
        If the printLast argument to self.__init__() was True or if at least
        one message has been printed already, this method ensures that the last
        message is printed.

        @note: Does _not_ close the output stream. Does not prevent you from
        further writing to the MsgTicker object, though this does not make any
        sense. Better use del to remove the MsgTicker object completely.
        """
        if self.printLast and self.lastArgs is not None:
            msg(self.msgTemplate % self.lastArgs,
                end=self.end, debugLevel=self.debugLevel)
            self.lastArgs = None

    def __del__(self):
        self.close()
