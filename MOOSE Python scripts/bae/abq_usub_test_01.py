"""abq_usub_test_01.py

Classes for testing Abaqus User Subroutines (Fortran and C) in python code.
Compiling is done implicitly and automatically by this module.

Usage
=====
 >>> from bae.abq_usub_test_01 import \\
 >>>     AbqUsubTest, ModuleProcedure, ModuleParameter
 >>> mytest = AbqUsubTest("source")
 >>> mytest.openLogFile(6, "abqTest.log")
 >>> myfunc = ModuleProcedure(mytest, "test_", [
 >>>     ("myfloat", ct.c_double, True),
 >>>     ("string", ct.c_char*80),
 >>>     ("array", ct.c_double*10),    ])
 >>> myvars = ModuleParameter(mytest.lib, [ ... ])
 >>> res = myfunc(myfloat=1.25, string="gugu", array=range(5))
 >>> mytest.closeLogFile(6)
 >>> print "result:", res[0]
 >>> print "static vars:", myvars
"""

__version__ = "1.09"

_version_history_ = """
1.00 : GP new
1.01 : added parseFortranConstants
1.02 : added matrixToFortran, matrixFromFortran
1.03 : added timeStepper, checkVector
1.04 : added: stop at compile error
1.05 : added: patch option to AbqUsubTest for variants in source code.
1.06 : added AbqUsubTest.openLogFile, closeLogFile
1.07 : when parsing fortran constants skip subroutine and function body
1.08 : add abaqusMainLoop and more features to AbqUsubTest
       incompatible change to ModuleProcedure and ModuleParameter __init__()
       now taking ctype.CDLL object
1.09 : GP added AbqUsubTest.getVumatCoords, AbqUsubTest.getVufieldCoords
"""

import ctypes as ct
import os
import re
from random import uniform
from subprocess import Popen, PIPE
from itertools import izip, chain

from bae.future_01 import OrderedDict

from bae.log_01 import msg

class CompileError(Exception):
    pass

def getGccVersion():
    """Returns the g++ compiler version as tuple of three integers like (6,3,0)
    """
    try:
        verStr = Popen(["g++", "--version"], stdout=PIPE).communicate()[0]
    except OSError:
        msg("g++ compiler executable not found.", debugLevel=1)
        return (0,0,0)
    verStr = verStr.splitlines()[0]  # get the first line only
    res = re.search(r"(\d+)\.(\d+)\.(\d+)", verStr)
    ver = map(int, res.groups())
    return ver

class FortranCompilerEnvironment(object):
    """Base class for different Fortran compiler interfaces.
    Subclasses are for specific compilers, i.e. Intel or GNU
    """
    def __init__(self):
        self.env = dict(os.environ)

    def getCommand(self, outputName, includeDirs=[]):
        cmd = list(self.cmd)

        for x in includeDirs:
            cmd.extend(("-I", x))

        cmd.extend(("-o", outputName))
        return cmd

    def getSymbol(self, genericName):
        """Each compiler creates different symbol names (at least other than
        ifort). Those names can be determined using objdump -t libtest.so.

        This function delivers the individula variant of those names for this
        compiler in order to be able to retrieve the symbols in the compiled
        library.
        If adapting the caller to those symbol names is not an option then
        the symbols could also be changed afterwards with
        objcopy --redefine-syms=filename. The file specified hereby contains
        "old new" lines.
        """
        raise NotImplementedError(
            "The class %s requires a getSymbol() function" % self.__class__)


class FortranEnvIfort(FortranCompilerEnvironment):
    """Compiler interface for ifort
    """

    def __init__(self, extraOptionsFortran=[]):
        """

        @param extraOptionsFortran: list of additional compiler (or linker)
           options to pass on to the Fortran compiler.
           E.g. extraOptionsFortran=["-check", "all"]
        """
        FortranCompilerEnvironment.__init__(self)

        # add some environment varibles (taken from abaqus_v6.env, I think)
        self.env["FOR_IGNORE_EXCEPTIONS"] = "1"
        self.env["FOR_DISABLE_STACK_TRACE"] = "1"
        self.env["IPATH_NO_CPUAFFINITY"] = "1"

        # ifort command taken from compile scripts
        self.cmd = ['/opt/intel/Compiler/11.1/072/bin/intel64/ifort',
               '-shared', '-fPIC', '-auto', '-mP2OPT_hpo_vec_divbyzero=F',
               '-extend_source', '-WB', '-cxxlib', '-threads', '-Wl,-Bdynamic',
               "-Xlinker", "-rpath",
               "-Xlinker", "/opt/intel/Compiler/11.1/072/lib/intel64",
               '-i-dynamic', '-lifport', '-lifcoremt',
               ]
        self.cmd += extraOptionsFortran
    
    def checkExecutable(self):
        "Check if the compiler can be called at all"
        if not os.path.isfile(self.cmd[0]):
            raise OSError("Fortran compiler executable does not exist: %s"
                          % self.cmd[0])

    def getSymbol(self, genericName):
        """convert generic symbol name to compiler dependent symbol name.
        I.e. convert "seqperiodcounter.seqPeriodIsFirstInc"
        to "seqperiodcounter_mp_seqperiodisfirstinc_"

        Each compiler creates different symbol names (at least other than
        ifort). Those names can be determined using objdump -t libtest.so.

        This function delivers the individula variant of those names for this
        compiler in order to be able to retrieve the symbols in the compiled
        library.
        If adapting the caller to those symbol names is not an option then
        the symbols could also be changed afterwards with
        objcopy --redefine-syms=filename. The file specified hereby contains
        "old new" lines.
        """
        try:
            a, b = genericName.split(".")
        except ValueError:
            return genericName
        return "%s_mp_%s_" % (a.lower(), b.lower())


class FortranEnvGfortran(FortranCompilerEnvironment):
    """Compiler interface for gfortran
    """
    def __init__(self, extraOptionsFortran=[]):
        FortranCompilerEnvironment.__init__(self)

        self.cmd = [
            'gfortran',
            '-shared', '-fPIC', '-ffixed-line-length-132',
            '-fbounds-check',
            ]
        # gfortran creates different symbol names than ifort. This can be
        # changed afterwards with objcopy --redefine-syms=filename. The file
        # specified hereby contains "old new" lines. The names to be redefined
        # can be determined using objdump -t libtest.so

        
    def checkExecutable(self):
        "Check if the compiler can be called at all"
        return

    def getSymbol(self, genericName):
        """convert generic symbol name to compiler dependent symbol name.
        I.e. convert "seqperiodcounter.seqPeriodIsFirstInc"
        to "__seqperiodcounter_MOD_seqperiodislastinc"

        Each compiler creates different symbol names (at least other than
        ifort). Those names can be determined using objdump -t libtest.so.

        This function delivers the individula variant of those names for this
        compiler in order to be able to retrieve the symbols in the compiled
        library.
        If adapting the caller to those symbol names is not an option then
        the symbols could also be changed afterwards with
        objcopy --redefine-syms=filename. The file specified hereby contains
        "old new" lines.
        """
        try:
            a, b = genericName.split(".")
        except ValueError:
            return genericName
        return "__%s_MOD_%s" % (a.lower(), b.lower())


class AbqUsubTest(object):
    """Interface object for Fortran subroutine calls.

    Basic usage:
     >>> from bae.abq_usub_test_01 import \\
     >>>     AbqUsubTest, ModuleProcedure
     >>> mytest = AbqUsubTest("source")
     >>> myfunc = ModuleProcedure(mytest, "test_", [ ... ])
     >>> res = myfunc(...)
     >>> print "result:", res

    @cvar moduleVarsDescr: optional list of (name, symbol, type)-tuples used to
       initialize the instance variable L{moduleVars}. To be defined in a
       subclass.

       To find out about the name (symbol) of a particular subroutine or
       variable (on Linux) use objdump::
         objdump -t libtest/libtest.so | less

    @ivar libPath: (absolute) path to the dynamic/shared library
    @ivar lib: the dynamic library, a ctypes.CDLL object
    @ivar moduleVars: A L{ModuleParameter} object referring to static variables
       in the library. Can be used to read or write them. Additionally the
       initial values of all variables in this dictionary are stored in the
       instance variable L{initValues} and will be restored to those values
       before each test procedure.
    @ivar initValues: dictionary of moduleVars with their initial values after
       the module has been compiled.

    @ivar nBlock: Must be defined in derived subclass. Number of integration
       points per call to vumat subroutine.
    @ivar nstatev: Must be defined in derived subclass. Number of SDVs per
       integration point.
    @ivar nprops: Must be defined in derived subclass. Number of material props.
    @ivar nfieldv: Must be defined in derived subclass. Number of field
       variables.
    @ivar nBlockNodes: Must be defined in derived subclass. Number of nodes
       per call to vufield.
    """
    extraTestSourcesSubDir = "abq_usub_test_extra"
    patchId = 1

    def __init__(self, sourceDir, patch=None,
                 sourceFilesFortran=["vumat.f"], sourceFilesC=["lrm.c"],
                 sourceFilesCSubDir=".",
                 libDir=None, patchedSourceDir="patchedSource",
                 abqLogFileName="abqTest.log",
                 extraOptionsFortran=[],
                 fortranCompiler=None):
        """
        @param sourceDir: subdir containing all source files.
        @param patch: describes modifications to the source files prior to
           compiling, for example to set parameters in the code. This is for
           testing variations and different versions of modules. Can be a
           function or a filename or a sequence of replace instructions in the
           form of regular expressions.

           Instead of supplying this parameter to __init__ there can as well be
           an attribute defined in the derived sub-class. A function must be
           defined as staticmethod in this case. The argument to __init__
           takes precedence over the class attribute.

           If it's a function then this function will be called with
           patchedSourceDir as single argument and is expected to apply the
           required changes to file(s) in that folder.

           If it's not callable and not a string it's assumed to be a sequence
           of (filename, patches)-tuples. filename identifies the file in the
           source folder that is to be modified, e.g. "vusubs_param.f" or
           "modules/seqPeriodCounter.f". patches is a list of changes which
           contains tuples of arguments to the re.sub() function. I.e. pattern,
           replacement and optionally count.
            >>> patch = [ ("vusubs_param.f", [
            ...     (r"(?i)(.*jobDir\s*=\s*)(.*)", "\\1'%s/'" % os.getcwd()),
            ...     # ...further replacements for vusubs_param.f
            ...     ] ), ]

           You can use multiline regular expressions. A backup of the original
           state will be made, with the filename "gaga.orig.f" for the original
           "gaga.f".

           If it's a filename (i.e a string) this file is read as a patch file
           that describes changes to the source code.

           A patch file should be created by the following commands.
           First create a copy of the source subdirectory named patchedSource.
           Secondly apply all required changes to patchedSource. Then run the
           diff command. The final sed command is required to make patch act
           on the patchedSource directory rather than source::
            >>> $ cp -r source patchedSource
            >>> $ emacs patchedSource/vusubs_param.f
            >>> $ diff -ruN source/ patchedSource/ > myversion.patch
            >>> $ sed -i 's:source/:patchedSource/:g' myversion.patch

           To apply a patch manually (not needed for this operation):
            >>> $ cp -r source patchedSource
            >>> $ patch -s -p0 <myversion.patch

        @param sourceFilesFortran: Note that include files need not be listed
           here. Just the main Fortran source files.
        @param sourceFilesC: List of main C source files. This is relative to
           the folder sourceFilesCSubDir.
        @param sourceFilesCSubDir: optional directory name where the C source
           files are. The C compiler will run in that folder. To be specified
           relative to sourceDir.
           
        @param libDir: Directory to hold intermediate files and the final
           object-file of the library. Defaults to "libtest_%s" % patchId
           with patchId being 1, 2, 3,... for subsequent invocations.
        @param patchedSourceDir: temporary subfolder for the modified (patched)
           source files. (Only if you don't like "patchedSource" as directory
           name.)
        @param abqLogFileName: filename of a file to be opened for writing
           as Fortran unit 6 which is supposed to be the Abaqus log file.

        @param fortranCompiler: instance of subclass of
           L{FortranCompilerEnvironment}, default to ifort
        """

        if not os.path.isdir(sourceDir):
            raise OSError("Source directory %s does not exist." % sourceDir)

        # select default fortran compiler
        if fortranCompiler is None:
            fortranCompiler = FortranEnvIfort()

        extraTestSourcesDir = os.path.abspath(
            os.path.join(
                os.path.dirname(__file__), self.extraTestSourcesSubDir))
        if not os.path.isdir(extraTestSourcesDir):
            raise CompileError(
                "Extra test sources directory %s not found."
                % extraTestSourcesDir)

        # If argument patch or instance has a patch attribute
        # then create a copy of the source directory to then patch it
        if patch or hasattr(self, "patch"):
            # clear patched source directory
            if os.path.isdir(patchedSourceDir):
                msg("Clearing secondary source directory for the patch: %s"
                    % patchedSourceDir, debugLevel=1)
                rc = Popen(['rm', '-rf', patchedSourceDir]).wait()
                if rc != 0:
                    raise OSError("Could not remove %s" % patchedSourceDir)

            # copy source
            rc = Popen(['cp', '-R', sourceDir, patchedSourceDir]).wait()
            if rc != 0:
                raise OSError("Could not copy %s to %s"
                              % (sourceDir, patchedSourceDir))

            # copy abq_usub_test_extra into the patchedSourceDir
            rc = Popen([
                'cp', '-R', extraTestSourcesDir, patchedSourceDir]).wait()
            if rc != 0:
                raise OSError("Could not copy %s to %s"
                              % (extraTestSourcesDir, patchedSourceDir))
            # later use the copied version as include path
            extraTestSourcesDir = os.path.join(
                os.path.abspath(patchedSourceDir), self.extraTestSourcesSubDir)

            # now apply the actual patch
            # (argument to __init__ takes precedence over class attribute)
            if callable(patch):
                # call patch function
                patch(patchedSourceDir)

            elif isinstance(patch, basestring):
                # apply patch file with the patch command
                rc = Popen(['patch', '-s', '-p0', '--input=%s' % patch]).wait()
                if rc != 0:
                    raise CompileError("Patch %s could not be applied."
                                       % patch)

                # check if any reject file has been created by the patch
                output = Popen(["find", patchedSourceDir, "-name", "*.rej"],
                               stdout=PIPE).communicate()[0]
                if output:
                    raise CompileError(
                        "Problems occurred when applying patch %s. The"
                        " following file(s) contain(s) modifications that"
                        " could not be applied:\n%s" % (patch, output))

            # If patch not given as argument then look for attribute of
            # actual instance of subclass.
            elif hasattr(self, "patch"):
                # call patch function
                self.patch(patchedSourceDir)
                
            else:
                # process a list of (filename, list-of-changes) -tuple
                for filename, patches in patch:
                    msg("Patching %s (%d patches)." % (filename, len(patches)),
                        debugLevel=1)
                    patchFile(os.path.join(patchedSourceDir, filename), patches)

            # use patchedSourceDir as sourceDir
            sourceDir = patchedSourceDir

        # diagnostic output
        msg("Extra test sources directory: %s" % extraTestSourcesDir,
            debugLevel=1)

        # compiled lib will be found here...
        if libDir is None:
            # We need a unique library name for each initialization because
            # the dynamic library loader caches libraries.
            libDir="libtest_%s" % AbqUsubTest.patchId
            AbqUsubTest.patchId += 1
        libDir = os.path.join(os.getcwd(), libDir)
        if not os.path.isdir(libDir):
            msg("Creating library directory %s" % libDir, debugLevel=1)
            os.makedirs(libDir)

        #-- compile c sources
        if sourceFilesC:

            # determine depth of C subdir
            subDirDepth = sourceFilesCSubDir.rstrip("/").count("/")

            # and then the adapted output path for the C object file
            croutinesPath = ["../"]*subDirDepth + [libDir, 'croutines.o']
            croutinesPath = os.path.join(*croutinesPath)

            if getGccVersion()[0]>=6:
                # use g++ to compile c routines
                cmd = [
                    "g++", "-c", "-fPIC", "-fpermissive", "-mpc64",
                    "-DABQ_LINUX", "-DABQ_LNX86_64", "-DFOR_TRAIL",
                    "-DHAS_BOOL", "-DASSERT_ENABLED", "-D_BSD_TYPES",
                    "-D_DEFAULT_SOURCE", "-D_GNU_SOURCE", "-D_POSIX_SOURCE",
                    "-D_XOPEN_SOURCE_EXTENDED", "-D_XOPEN_SOURCE",
                    "-DHAVE_OPENGL", "-DHKS_OPEN_GL", "-DTYPENAME=typename",
                    "-DGL_GLEXT_PROTOTYPES", "-D_LARGEFILE64_SOURCE",
                    "-D_FILE_OFFSET_BITS=64",
                    "-O0",
                    '-o', croutinesPath]
            else:
                # use Intel compiler for c code
                # icpc command taken from compile scripts
                cmd = ['/opt/intel/Compiler/11.1/072/bin/intel64/icpc',
                       '-c', '-cxxlib', '-Kc++eh', '-fPIC', '-Krtti', '-Kc++',
                       '-pc64', '-restrict', '-DABQ_LINUX', '-DABQ_LNX86_64',
                       '-DFOR_TRAIL', '-DHAS_BOOL', '-DASSERT_ENABLED',
                       '-D_BSD_TYPES', '-D_BSD_SOURCE', '-D_GNU_SOURCE',
                       '-D_POSIX_SOURCE', '-D_XOPEN_SOURCE_EXTENDED',
                       '-D_XOPEN_SOURCE', '-DHAVE_OPENGL', '-DHKS_OPEN_GL',
                       '-DTYPENAME=typename', '-DGL_GLEXT_PROTOTYPES',
                       '-D_LARGEFILE64_SOURCE', '-D_FILE_OFFSET_BITS=64',
                       '-we1011', '-we120', '-we117', '-we556', '-we144',
                       '-we268', '-we1224', '-we167', '-we880', '-O0',
                       # this is a mistake in the compile scripts, has no effect
                       # '-I%I',
                       '-o', croutinesPath]

                # check if the Intel compiler is available
                if not os.path.isfile(cmd[0]):
                    raise OSError("C compiler executable does not exist: %s"
                                  % cmd[0])

            cmd += sourceFilesC
            msg("Running the C compiler with the following command:\n%s"
                % cmd, debugLevel=2)
            p = Popen(cmd, cwd=os.path.join(sourceDir, sourceFilesCSubDir),
                      stdout=PIPE, stderr=PIPE)
            (pOut, pErr) = p.communicate()
            rc = p.returncode
            del p
            if rc != 0:
                raise CompileError(
                    "The C compiler returned an error (rc=%d)."
                    "\n...command:\n%s"
                    "\n...compiler output:\n%s"
                    "\n...compiler error output:\n%s"
                    % (rc, " ".join(cmd), pOut, pErr))
            msg("Compiled C sources: %s" % ",".join(sourceFilesC), debugLevel=1)

        #-- compile and link the fortran code
        libPath = os.path.join(libDir, "libtest.so")  # needed again later...
        fortranCompiler.checkExecutable()

        env = fortranCompiler.env

        cmd = fortranCompiler.getCommand(
            libPath, includeDirs=[
                extraTestSourcesDir,
                os.getcwd(),  # include files in the current folder
                ])
        cmd += sourceFilesFortran
        cmd.append(os.path.join(extraTestSourcesDir, "abqdummy.f"))
        if sourceFilesC:
            cmd.append(croutinesPath)

        msg("Running the Fortran compiler with the following command:\n%s"
            % cmd, debugLevel=2)
        p = Popen(cmd, cwd=sourceDir, env=env, stdout=PIPE, stderr=PIPE)
        (pOut, pErr) = p.communicate()
        rc = p.returncode
        msg("Finished compiling, return code: %d" % rc, debugLevel=2)
        msg("...compiler output\n%s" % pOut, debugLevel=2)
        msg("...compiler error output\n%s" % pErr, debugLevel=2)
        del p
        if rc != 0:
            raise CompileError(
                "The Fortran compiler returned an error (rc=%d)."
                "\n...command:\n%s"
                "\n...compiler output:\n%s"
                "\n...compiler error output:\n%s"
                % (rc, " ".join(cmd), pOut, pErr))

        self.libPath = libPath
        self.lib = ct.CDLL(libPath)

        # initialize self.moduleVars ... static variables in the library
        if hasattr(self, "moduleVarsDescr"):
            # apply compiler dependent symbol name pattern to the given
            # generic symbol names in self.moduleVarsDescr which can be
            # defined as sub-class attribute.
            # I.e. convert "seqperiodcounter.seqPeriodIsFirstInc"
            # to "seqperiodcounter_mp_seqperiodisfirstinc_"
            moduleVarsDescr = [
                (name, fortranCompiler.getSymbol(symbol), type_)
                for (name, symbol, type_) in self.moduleVarsDescr]
            self.moduleVars = ModuleParameter(self.lib, moduleVarsDescr)
            self.initValues = self.moduleVars.save()

        # open Abaqus log file on Fortran unit 6
        if abqLogFileName:
            self.openLogFile(6, abqLogFileName)

    def __del__(self):
        """Close Abaqus log file on Fortran unit 6.
        """
        self.closeLogFile(6)

        # close the dynamic library
        # This seems to not work: Calling dlclose results in a segmentation
        # fault. Sometimes, e.g. on Debian Stretch. Not on Debian Wheezy,
        # though.
        # Make sure that different versions of the library have different paths
        # otherwise the dlopen-mechanism doesn't realize that this is intended
        # to be a different/modified library!
        # msg("Closing the dynamic library...", debugLevel=10)
        # handle = self.lib._handle  # obtain the SO handle
        # msg("... got lib._handle", debugLevel=10)
        # ct.cdll.LoadLibrary('libdl.so').dlclose(handle)
        # msg("... finished closing the library.", debugLevel=10)

    def initTest(self):
        """To be called once before each individual test to initialize all
        variables, arguments...
        Usually this method is being overloaded by subclasses to do further
        initializations like L{initVumat} and L{initVufield}.
        """
        self.moduleVars.set(**self.initValues)

    def initVumat(self):
        """Initialize vumat and its arguments. Suggested usage is in derived
        classes that define the following class variables. They are required
        for this function to work: nBlock, nfieldv, nstatev, nprops.
        """
        self.vumat = VumatPrototype(
            self.lib, NBLOCK=self.nBlock,
            NSTATEV=self.nstatev, NFIELDV=self.nfieldv, NPROPS=self.nprops)

    def initVufield(self):
        """Initialize vufield and its arguments. Suggested usage is in derived
        classes that define the following class variables. They are required
        for this function to work: nBlockNodes, nfieldv
        """
        # number of nodes: nBlockNodes
        # node numbers are 1, 2, 3, 4, ..., nBlockNodes
        self.vufield = ModuleProcedure(
            self.lib, "vufield_",   # module, symbol
            [   # rUserField(nBlock,nComp,nField); nComp=1
                ("RUSERFIELD",
                 ct.c_double*(self.nBlockNodes*self.nfieldv), True),

                ("NBLOCK", ct.c_int, False, self.nBlockNodes),
                ("NFIELD", ct.c_int, False, self.nfieldv),
                ("KFIELD", ct.c_int, False, -1),
                ("NCOMP", ct.c_int, False, 1),

                ("KSTEP", ct.c_int, False, 1),
                ("JFLAGS", ct.c_int*2, False, [1,0]),
                ("JNODEUID",
                 ct.c_int*self.nBlockNodes, False,range(1,self.nBlockNodes+1)),
                ("TIMESVEC", ct.c_double*4, False),

                ("COORDS", ct.c_double*(3*self.nBlockNodes), False),
                ("U", ct.c_double*(8*self.nBlockNodes), False),
                ("V", ct.c_double*(8*self.nBlockNodes), False),
                ("A", ct.c_double*(8*self.nBlockNodes), False)
            ])

    def openLogFile(self, unit, filename=None):
        """Open file for the given fortran unit.
        If filename not given write to a scratch file.
        """
        if not filename:
            self.lib.openlogfilescratch_(ct.byref(ct.c_int(unit)))
        else:
            self.lib.openlogfile_(
                ct.byref(ct.c_int(unit)),
                ct.c_char_p(filename), ct.c_size_t(len(filename)))
        return

    def closeLogFile(self, unit):
        """Close file for the given fortran unit.
        """
        self.lib.closelogfile_(ct.byref(ct.c_int(unit)))
        return

    def setVumatRandomCoords(self, zeroLevel=0.0, zRelMin=-1000.0, zRelMax=0.0,
                             x=5555.5, y=-3333.3):
        """initialize random integration point coordinates for vumat
        """
        self.zRelVec = [uniform(zRelMin, zRelMax) for _ in range(self.nBlock)]
        coords = [ [x, y, zRel+zeroLevel] for zRel in self.zRelVec ]
        self.vumat.setParam( COORDMP=matrixToFortran(coords) )

    def setVufieldRandomCoords(
            self, zeroLevel=0.0, zRelMin=-1000.0, zRelMax=0.0,
            x=5588.8, y=-3366.6):
        """initialize random node coordinates for vufield
        """
        self.zRelVecNodes = [uniform(zRelMin, zRelMax)
                             for _ in range(self.nBlockNodes)]
        coords = [ [x]*self.nBlockNodes, [y]*self.nBlockNodes,
                   [zRel+zeroLevel for zRel in self.zRelVec], ]
        self.vufield.setParam( COORDS=matrixToFortran(coords) )

    def getSdvInitVals(self, sdvInitVals=None):
        """Duplicate the given vector of sdvs (i.e. 20 floats) to an array
        of vectors, i.e. 100x20 values. Values are deep copied so you can later
        change single values for just one integration point.

        Possible use case:
         >>> sdvs = testLib.getSdvInitVals()
         >>> sdvs[0][iSdvState] = 4.0
         >>> testLib.vumat( ..., STATEOLD=matrixToFortran(sdvs), ...)

        ...Or:
         >>> sdvs = testLib.getSdvInitVals()
         >>> sdvs[0][iSdvState] = 4.0
         >>> testLib.vumat.setParam(STATEOLD=matrixToFortran(sdvs))

        @param sdvInitVals: optional vector of self.nstatev values. If not
          given then self.defaultSdvs is being used.

        @returns: list of list of values to be modified and passed on to
        a vumat-call later via matrixToFortran().
        """
        if sdvInitVals is None:
            sdvInitVals = self.defaultSdvs
        sdvs = [list(sdvInitVals) for _ in range(self.nBlock)]
        return sdvs

    def getVumatCoords(self):
        """Get current values of the integration point coordinates passed to
        the vumat subroutine.

        Values may have been initialized by L{setVumatRandomCoords}.

        @returns: a NBLOCK x 3 matrix. I.e. Get the y-coord for the 50th
           integration point by:
            >>> coords = self.testLib.getVumatCoords()
            >>> y50 = coords[49][1]
        """
        coords = matrixFromFortran(self.vumat.param["COORDMP"], self.nBlock)
        return coords

    def getVufieldCoords(self):
        """Get current values of the nodal coordinates passed to
        the vufield subroutine.

        Values may have been initialized by L{setVufieldRandomCoords}.

        @returns: a 3 x NBLOCK matrix. I.e. Get the y-coord for the 50th
           node by:
            >>> coords = self.testLib.getVufieldCoords()
            >>> y50 = coords[1][49]
        """
        coords = matrixFromFortran(self.vufield.param["COORD"], 3)
        return coords

    def setSDV(self, iIP, iSDV, value):
        """set single values in the initial sdvs (STATEOLD)

        @param iIP: index of the integration point (zero based!)
        @param iSDV: index of the SDV (zero based!)
        @param value: value to assign
        """
        self.vumat.param["STATEOLD"][iIP + iSDV*self.nBlock] = value

    def setStep2Time(self, steptime, dt=None):
        """set steptime and totaltime and seqPeriodIdx for/in step 2

        @param steptime: ... time in step-2
        @param dt: optional, if not given then assume seqPeriodLength/100
           Don't set seqPeriodIdx if the following call to vumat is supposed to
           be the first increment in the step indicated by dt==steptime. Leave
           seqPeriodIdx be the initial value in this case.
        """
        eqStepTime = self.vusubsParam["eqStepTime"]
        seqPeriodLength = self.vusubsParam["seqPeriodLength"]
        if dt is None:
            dt = min(steptime, 0.01 * seqPeriodLength)

        # don't set seqPeriodIdx in first call to vumat!
        if steptime != dt:
            self.moduleVars.set(
                seqPeriodIdx=int((steptime-dt)/seqPeriodLength)+1)

        self.vumat.setParam(
            STEPTIME=steptime,
            TOTALTIME=steptime+eqStepTime,
        )

    def abaqusMainLoop(
            self,
            dt, tEnd, checkTimes=None, stepStartTime=0.0, stepNb=1,
            rUserField=None, UMatrix=None,
            stress=None, sdvs=None,
            updateUMatrix=None,
            getDeltaEps=None,
            checkAfterVufield=None,
            checkAfterVumat=None,
            nodeWeights4Gpts=None,
            nbCallsPerIncr=1):
        """Run vufield and vumat for one Abaqus step (i.e. many increments,
        potentially many seqPeriods)

        @param checkTimes: If given then check-functions are called only at
          listed step times.
        @param stepStartTime: set to eqStepTime for step=2 (mining step)
        @param rUserField: initial values for vufield, a list of nbCallsPerIncr
          matrices with dimension nBlockNodes x nfieldv each.
        @param UMatrix:initial values for displacements, a list of
          nbCallsPerIncr matrices with dimension 8 x nBlockNodes each.
        @param stress, sdvs: initial values for vumat, will be modified. Both
          are lists of nbCallsPerIncr matrices. The matrices are nBlock x 6 and
          nBlock x nstatev respectively.
        @param nodeWeights4Gpts: for interpolation of nodal values to
          Gauss-point values. For each Gauss-point (max. number is self.nBlock)
          we have a list of (node-index, weight)-tuples. node-index is the
          zero-based index in ruserfield (max. self.nBlockNodes-1)

          Example: For the first Gauss pt take the value at node 12 (index 11)
          and for the second Gauss pt take the average of node 22 and 23:
           >>> nodeWeights4Gpts = [
           >>>     [(11, 1.0),],
           >>>     [(21, 0.5), (22, 0.5)],  ]

        @param updateUMatrix: a function to update displacements. It takes the
          arguments cntCallsPerIncr, steptime, dt, and the old UMatrix and
          delivers an updated version of UMatrix.
          cntCallsPerIncr is an int to track the multiple calls for each time
          increment. It will be zero for the first invocation at a particular
          steptime and one for the next invocation for the same steptime.
          UMatrix is the 8 x nBlockNodes matrix of current displacements for
          input to vufield. It will be passed in the value at the previous
          increment.
        @param getDeltaEps: a function that delivers the new strain increment.
          It takes the arguments cntCallsPerIncr, steptime, dt, and deltaEpsMat
          and delivers the new strain increments / deltaEpsMat.
          cntCallsPerIncr is an int to track the multiple calls for each time
          increment. It will be zero for the first invocation at a particular
          steptime and one for the next invocation for the same steptime.
          deltaEpsMat is a template nBlock x 6 matrix initialized to all zeros
          on input to getDeltaEps at each invocation.
        @param checkAfterVufield: function for checks, takes arguments
          cntCallsPerIncr, steptime, dt, ruserfield (as
          nBlockNodes x nfieldv matrix)
          Will be called right after the call to vufield if steptime has
          reached the next item of checkTimes.
        @param checkAfterVumat: function for checks, takes arguments
          cntCallsPerIncr, steptime, dt, stress (nBlock x 6 matrix),
          sdvs (nBlock x nstatev matrix)
          Will be called right after the call to vumat if steptime has
          reached the next item of checkTimes.
        @param nbCallsPerIncr: vumat and vufield will be called multiple times
          for the same increment (same steptime). Stress and state will be
          conserved and passed on separately for each of those lines of calls.

        @returns: Fields required for transition to the next Abaqus step:
          rUserField, UMatrix, stress, sdvs
        """

        # initial values for rUserField, later it's only modified by vufield
        # nBlockNodes x nFieldv matrix
        if rUserField is None:
            rUserField = [
                [0.0]*(self.nBlockNodes*self.nfieldv)
                for i in range(nbCallsPerIncr) ]
        else:
            if len(rUserField)!=nbCallsPerIncr:
                raise ValueError(
                    "We need nbCallsPerIncr sets of rUserField-matrices. You"
                    " supplied nbCallsPerIncr=%d and len(rUserField)=%d."
                    % (nbCallsPerIncr, len(rUserField)))

        # initial values for UMatrix; 8 x nBlockNodes matrix
        if UMatrix is None:
            UMatrix = [
                [[0.0]*self.nBlockNodes for i in range(8)]
                for i in range(nbCallsPerIncr) ]
        else:
            if len(UMatrix)!=nbCallsPerIncr:
                raise ValueError(
                    "We need nbCallsPerIncr sets of UMatrix-matrices. You"
                    " supplied nbCallsPerIncr=%d and len(UMatrix)=%d."
                    % (nbCallsPerIncr, len(UMatrix)))

        # initial values for stress; a nBlock x 6 matrix
        if stress is None:
            stress = [ ([ [0.0]*6 ]*self.nBlock)
                       for i in range(nbCallsPerIncr) ]
        else:
            if len(stress)!=nbCallsPerIncr:
                raise ValueError(
                    "We need nbCallsPerIncr sets of stress. You supplied"
                    " nbCallsPerIncr=%d and len(stress)=%d."
                    % (nbCallsPerIncr, len(stress)))

        # initial values for sdvs; a nBlock x nstatev matrix
        if sdvs is None:
            sdvs = [self.getSdvInitVals() for i in range(nbCallsPerIncr) ]
        else:
            if len(sdvs)!=nbCallsPerIncr:
                raise ValueError(
                    "We need nbCallsPerIncr sets of sdvs. You supplied"
                    " nbCallsPerIncr=%d and len(sdvs)=%d."
                    % (nbCallsPerIncr, len(sdvs)))

        # fieldOld, fieldNew: interpolated field values for vumat:
        # nBlock x nfieldv matrix
        fieldOld = [
            [0.0]*(self.nfieldv*self.nBlock)
            for i in range(nbCallsPerIncr) ]
        # just in case there is no nodeWeights4Gpts: initialize with zeros
        fieldNew = [0.0]*(self.nfieldv*self.nBlock)

        incrementNb = 0
        for (steptime, dt, checkFlag) in timeStepper(dt, tEnd, checkTimes):

            incrementNb += 1

            for cntCallsPerIncr in range(nbCallsPerIncr):
                if updateUMatrix is not None:
                    UMatrix[cntCallsPerIncr] = updateUMatrix(
                        cntCallsPerIncr, steptime, dt,
                        UMatrix[cntCallsPerIncr])
                    self.vufield.setParam(U=matrixToFortran(
                        UMatrix[cntCallsPerIncr]))

                # call vufield, result is a list with one item: rUserField
                res = self.vufield(
                    KSTEP=stepNb,
                    JFLAGS=[incrementNb,0],
                    TIMESVEC=[steptime, dt, steptime, stepStartTime+steptime],
                    RUSERFIELD=rUserField[cntCallsPerIncr])
                rUserField[cntCallsPerIncr][:] = res[0]  # copy result array
                rufld = rUserField[cntCallsPerIncr]  # abbreviation

                if checkFlag and checkAfterVufield is not None:
                    checkAfterVufield(
                        cntCallsPerIncr,
                        steptime, dt,
                        matrixFromFortran(rufld, self.nBlockNodes))

            for cntCallsPerIncr in range(nbCallsPerIncr):

                # interpolate field from nodes to Gauss-points
                if nodeWeights4Gpts is not None:
                    rufld = rUserField[cntCallsPerIncr]  # abbreviation
                    fieldNew = [0.0]*(self.nfieldv*self.nBlock)  # init zeros
                    for iComp in range(self.nfieldv):
                        i1 = iComp*self.nBlock
                        i2 = i1+len(nodeWeights4Gpts)
                        fieldNew[i1:i2] = [
                            sum(weight*rufld[iNode+self.nBlockNodes*iComp]
                                for iNode, weight in nodeWeights)
                            for nodeWeights in nodeWeights4Gpts]

                if getDeltaEps is not None:
                    deltaEpsMat = getDeltaEps(
                         cntCallsPerIncr, steptime, dt,
                         [[0.0]*6 for i in range(self.nBlock)])
                    self.vumat.setParam(STRAININC=matrixToFortran(deltaEpsMat))

                # call vumat
                # (fieldOld has been set at the end of previous iteration)
                stressnew, statenew = self.vumat(
                    STEPTIME=steptime,
                    TOTALTIME=stepStartTime+steptime,
                    DT=dt,
                    STRESSOLD=matrixToFortran(stress[cntCallsPerIncr]),
                    STATEOLD=matrixToFortran(sdvs[cntCallsPerIncr]),
                    FIELDOLD=fieldOld[cntCallsPerIncr],
                    FIELDNEW=fieldNew,
                    )
                stress[cntCallsPerIncr] = matrixFromFortran(stressnew, self.nBlock)
                sdvs[cntCallsPerIncr] = matrixFromFortran(statenew, self.nBlock)
                if checkFlag and checkAfterVumat is not None:
                    checkAfterVumat(
                        cntCallsPerIncr,
                        steptime, dt,
                        stress[cntCallsPerIncr], sdvs[cntCallsPerIncr])

                # copy fieldNew to fieldOld for next iteration
                fieldOld[cntCallsPerIncr] = fieldNew

        # for transition to the next Abaqus step return results from vumat
        return (rUserField, UMatrix, stress, sdvs)


class ModuleProcedure(object):
    """Object representing a Fortran procedure including a datastructure to
    hold all parameters needed for a function call.
     >>> from bae.abq_usub_test_01 import \\
     >>>     AbqUsubTest, ModuleProcedure
     >>> mytest = AbqUsubTest("source")
     >>> vumat = ModuleProcedure(
     >>>     mytest.lib, "vumat_", [
     >>>         ("myfloat", ct.c_double, False, None),
     >>>         ("number", ct.c_int, False, None),
     >>>         ("string", ct.c_char*80, False, None),
     >>>         ("array", ct.c_double*10, False, None),
     >>>         ("array2", ct.c_int*8, False, None),     ])
     >>> vumat.setParam(myfloat=1.25, string="gugu")
     >>> vumat.param["array"][0] = 5.0
     >>> vumat.setParam(array2=range(8), number=8)
     >>> vumat()

    There are function prototypes in ctypes so why do we need this class?
    Because our Fortran procedures get parameters passed in by reference. And
    this interface seems simpler to me.

    @ivar lib: a ctypes.CDLL library object like AbqUsubTest.lib passed in
        through __init__()
    @ivar symbol: The internal name of the procedure in the module/library.
    @ivar param: A dict referring to the arguments of the procedure. This can
        be used to query values and to set individual components of an array:
         >>> vumat.param["array"][0] = 5.0

        Setting single values or complete arrays must be done through the
        setParam() method.
    @ivar namesInp: List of just the parameter names.
    @ivar namesOut: List of parameter names with the outputFlag set. I.e. names
        of all output parameters.

    @Note: Arrays are always flat. Note that Fortran iterates multiple
       indexes (i,j,k) in the order i,j,k. E.g. the Fortran matrix
       COORDMP(10, 3) is stored and would be assigned in the following order
       ("row major order"):
       vumat["COORDMP"] = [ COORDMP_1_1, COORDMP_2_1, COORDMP_3_1, ...,
       COORDMP_9_1, COORDMP_10_1, COORDMP_1_2, COORDMP_2_2, ...,
       COORDMP_9_2, COORDMP_10_2, COORDMP_1_3, ..., COORDMP_10_3 ]

    """

    def __init__(self, lib, symbol, parameters=[]):
        """
        @param lib: a ctypes.CDLL library object like AbqUsubTest.lib
        @param symbol: The internal name of the procedure in the module.
        @param parameters: list of (name, type, outputFlag, default)
           -tuples describing the parameters of the procedure call. The last
           two items are optional.

           name is the key by which this particular parameter shall be
           accessible. type is a ctypes-type.

           outputFlag determines if this argument is considered a output
           quantity and will therefore be in the result-tuple of a procedure
           call or not. Even if this flag is False the corresponding value
           will be accessible after the procedure call like all other
           arguments to the procedure. (I.e. its value can be retrieved and it
           can therefore be used as an output quantity.)

           default can be None or otherwise it serves as
           initial (default) value.
        """
        self.lib = lib
        self.symbol = symbol
        self.param = dict()
        self.namesInp = list()
        self.namesOut = list()
        self.fortranArgs = list()
        stringLengths = list()

        for paraTuple in parameters:
            name, type_ = paraTuple[:2]

            try:
                default = paraTuple[3]
            except IndexError:
                default = None

            # initialize ctypes-variable, assign default if given
            if default is None:
                param = type_()
            elif isinstance(default, basestring):
                param = type_()
                param[:] = default.ljust(len(param))
            elif isinstance(default, (list, tuple)):
                param = type_(*default)
            else:
                param = type_(default)
            self.param[name] = param

            # update self.namesInp
            self.namesInp.append(name)

            # update self.namesOut
            try:
                outputFlag = paraTuple[2]
            except IndexError:
                pass
            else:
                if outputFlag:
                    self.namesOut.append(name)

            # update self.fortranArgs: argument(s) for the fortran procedure
            self.fortranArgs.append(ct.byref(param))
            if issubclass(type_, ct.Array) and issubclass(type_._type_, ct.c_char):
                stringLengths.append(len(param))

        # the lengths of string arguments (e.g. CHARACTER*80) are passed as
        # long int arguments after all other
        self.fortranArgs.extend(ct.c_size_t(x) for x in stringLengths)

    def __str__(self):
        """Nice printable list of parameters
        """
        valuesList = list()
        for name in self.namesInp:
            param = self.param[name]
            try:
                value = param.value
            except AttributeError:
                if isinstance(param, ct.Array):
                    value = list(param)
                else:
                    value = None
            valuesList.append((name, value))
        return "\n".join(("  %s: %s" % x) for x in valuesList)

    def setParam(self, *args, **kwargs):
        """Set paramater value(s) for next procedure call.
        Set many values in one go. Accepts positional and keyword arguments.

         >>> vumat = ModuleProcedure(
         >>>     mytest.lib, "vumat_", [
         >>>         ("number", ct.c_int, True, 10),
         >>>         # for Fortran: CHARACTER*80 string
         >>>         ("string", ct.c_char*80, False, None), ])
         >>> vumat.setParam(5, string="hallo")

        Arrays:

        To set the whole array use this method setParam(), possibly with the
        help of matrixToFortran():
         >>> N = ...
         >>> vumat = ModuleProcedure(
         >>>     mytest.lib, "vumat_", [
         >>>         # for Fortran: dimension myArray(N)
         >>>         ("myArray", ct.c_double*N, True, 10),
         >>>         # for Fortran: dimension coords(N,3)
         >>>         ("coords", ct.c_double*(N*3), True, 10), ])
         >>> vumat.setParam(myArray=[1.0, 10.4, 4.2e5, ...])
         >>> vumat.setParam(coords=matrixToFortran([
         >>>     [0.5, 1.0, 3.2],
         >>>     [1.5, 4.1, 1.1],
         >>>     ... ]))

        Setting individual components of an array can be done like so:
         >>> vumat.param["myArray"][0] = 2.0   # before was: 1.0
         >>> idxRow = 0; idxCol = 1
         >>> # before: coords[idxRow][idxCol] == 1.0
         >>> vumat.param["coords"][idxRow+N*idxCol] = 2.0
        """
        # programmers note:
        # - izip ignores surplus items if one of its arguments is longer than
        #   the other. I.e. the result is as long as the shortest of its
        #   arguments. If args is empty the izip yield nothing.
        # - chain chains one iterator after the other.
        # ... the result is that we first assign positional arguments to the
        # corresponding self.param - key-value-tuple according to self.namesInp
        # and then process the keyword arguments.
        for name, value in chain(
                izip(self.namesInp, args), kwargs.iteritems()):

            param = self.param[name]
            if isinstance(param, ct.Array):
                type_ = param._type_
                len_ = len(param)

                if issubclass(type_, ct.c_char):
                    # special treatment of character arrays:
                    # fill up to full array length with spaces
                    self.param[name][:] = value.ljust(len_)
                else:
                    # array of supposedly numeric values
                    # value is some kind of list or tuple and will be assigned
                    # assuming self.param[name] is an array
                    try:
                        self.param[name][:len_] = value
                    except TypeError, exc:
                        raise TypeError(
                            "%s\n... trying to assign %s to %s"
                            % (exc.args[0], value, name))

            else:
                # value is an atomic value and will be assigned directly
                self.param[name].value = value

    def __call__(self, *args, **kwargs):
        """Call the procedure with the current set of parameters.

        Optional positional and keyword arguments might be given to alter
        parameters beforehand.

        The result is a tuple of values for all output parameters as marked by
        the outputFlag in the parameters-argument to __init__().

         >>> vumat = ModuleProcedure(
         >>>     mytest.lib, "vumat_", [
         >>>         ("number", ct.c_int, True, 10),
         >>>         ("string", ct.c_char*80, False, None), ])
         >>> res = vumat(5, string="hallo")
         >>> print "f(number) = ", res[0]
         f(number) = 20
         >>> print "f(number) = ", vumat["number"]
         f(number) = c_int(20)
        """
        # process new parameters
        self.setParam(*args, **kwargs)
        
        # call the procedure
        procedure = getattr(self.lib, self.symbol)

        procedure(*self.fortranArgs)

        # collect output parameters into result value
        res = list()
        for pname in self.namesOut:
            try:
                res.append(self.param[pname].value)
            except AttributeError:
                res.append(self.param[pname])
        return tuple(res)

class VumatPrototype(ModuleProcedure):
    """Subclass for vumat calls. Initialize the parameters with
    some base parameters like NPROPS and NSTATEV. Before calling it you only
    need to set the remaining arguments like STEPTIME, DT, STRAININC (still
    a whole lot...)
     >>> from bae.abq_usub_test_01 import \\
     >>>     AbqUsubTest, VumatPrototype
     >>> mytest = AbqUsubTest("source")
     >>> vumat = VumatPrototype(mytest.lib, 100, 8, 0, 21)
     >>> stressnew, statenew = vumat(
     >>>     STEPTIME=1.0,
     >>>     TOTALTIME=2.0,
     >>>     COORDMP=[x for vec in zip(*coords) for x in vec],
     >>>     PROPS=(UCS, alpha, trallala, ....),
     >>>     STRAININC=[x for vec in zip(*deltaEps) for x in vec],
     >>>     ...,
     >>>     )

    This is assuming coords is a list of xyz-tuples and deltaEps is a list
    of strain increment 6-tuples, each of the tuples representing a tensor
    with following components: 11,22,33,12,23,31.
    """
    def __init__(self, lib, NBLOCK, NSTATEV, NFIELDV, NPROPS,
                 NDIR=3, NSHR=3, symbol="vumat_"):

        ModuleProcedure.__init__(
            self, lib, symbol,
            [
                ("NBLOCK",  ct.c_int, False, NBLOCK),
                ("NDIR",    ct.c_int, False, NDIR),    # usually: 3
                ("NSHR",    ct.c_int, False, NSHR),    # usually: 3
                ("NSTATEV", ct.c_int, False, NSTATEV),
                ("NFIELDV", ct.c_int, False, NFIELDV),
                ("NPROPS",  ct.c_int, False, NPROPS),

                ("LANNEAL",   ct.c_int),
                ("STEPTIME",  ct.c_double),
                ("TOTALTIME", ct.c_double),
                ("DT",        ct.c_double),
                ("CMNAME",    ct.c_char*80),     # CHARACTER*80 CMNAME

                ("COORDMP",       ct.c_double*(NBLOCK*3)),
                ("CHARLENGTH",    ct.c_double*NBLOCK),
                ("PROPS",         ct.c_double*NPROPS),
                ("DENSITY",       ct.c_double*NBLOCK),
                ("STRAININC",     ct.c_double*(NBLOCK*(NDIR+NSHR))),
                ("RELSPININC",    ct.c_double*(NBLOCK*NSHR)),
                ("TEMPOLD",       ct.c_double*NBLOCK),
                ("STRETCHOLD",    ct.c_double*(NBLOCK*(NDIR+NSHR))),
                ("DEFGRADOLD",    ct.c_double*(NBLOCK*(NDIR+NSHR+NSHR))),
                ("FIELDOLD",      ct.c_double*(NBLOCK*NFIELDV)),
                ("STRESSOLD",     ct.c_double*(NBLOCK*(NDIR+NSHR))),
                ("STATEOLD",      ct.c_double*(NBLOCK*NSTATEV)),
                ("ENERINTERNOLD", ct.c_double*NBLOCK),
                ("ENERINELASOLD", ct.c_double*NBLOCK),
                ("TEMPNEW",       ct.c_double*NBLOCK),
                ("STRETCHNEW",    ct.c_double*(NBLOCK*(NDIR+NSHR))),
                ("DEFGRADNEW",    ct.c_double*(NBLOCK*(NDIR+NSHR+NSHR))),
                ("FIELDNEW",      ct.c_double*(NBLOCK*NFIELDV)),
                ("STRESSNEW",     ct.c_double*(NBLOCK*(NDIR+NSHR)), True),
                ("STATENEW",      ct.c_double*(NBLOCK*NSTATEV), True),
                ("ENERINTERNNEW", ct.c_double*NBLOCK),
                ("ENERINELASNEW", ct.c_double*NBLOCK),

                ])


class ModuleParameter(object):
    """A datastructure to easily access parameters statically stored in a
    library.
    It behaves like a dict.
     >>> from bae.abq_usub_test_01 import \\
     >>>     AbqUsubTest, ModuleParameter
     >>> mytest = AbqUsubTest("source")
     >>> params = ModuleParameter(
     >>>     mytest.lib, fields=[
     >>>         ("myfloat", "mymodule_mp_myfloat_", ct.c_double),
     >>>         ("number", "mymodule_mp_number_", ct.c_int),
     >>>         ("string", "mymodule_mp_string_", ct.c_char*80),
     >>>         ("array", "mymodule_mp_array_", ct.c_double*10),
     >>>         ("array2", "mymodule_mp_array2_", ct.c_int*8),     ])
     >>> params["myfloat"] = 1.25
     >>> params["string"] = "gugu"
     >>> params["array"][0] = 5.0
     >>> params.set(array2=range(8), number=8)
     >>> print params

    @ivar fields: OrderedDict {name: (symbol, type)-tuple}. Will be
       initialized by the fields argument to __init__().
    @Note: Arrays are always flat. Note that Fortran iterates multiple
       indexes (i,j,k) in the order i,j,k. E.g. the Fortran matrix
       COORDMP(10, 3) is stored and would be assigned in the following order:
       params["COORDMP"] = [ COORDMP_1_1, COORDMP_2_1, COORDMP_3_1, ...,
       COORDMP_9_1, COORDMP_10_1, COORDMP_1_2, COORDMP_2_2, ...,
       COORDMP_9_2, COORDMP_10_2, COORDMP_1_3, ..., COORDMP_10_3 ]
    """

    def __init__(self, lib, fields=[]):
        """
        @param lib: a ctypes.CDLL library object like AbqUsubTest.lib
        @param fields: list of (name, symbol, type)-tuples identifying the
           module variables to be accessed. Used to initialize the instance
           variable L{ModuleParameter.fields}.
        """
        self.lib = lib
        self.fields = OrderedDict(
            (name, (symbol, type_))
            for (name, symbol, type_) in fields)

    def __getitem__(self, key):
        symbol, type_ = self.fields[key]
        return type_.in_dll(self.lib, symbol)

    def __setitem__(self, key, value):
        symbol, type_ = self.fields[key]
        try:
            # if value is a ctypes-type then use it's .value-attribute
            value = value.value
        except AttributeError:
            pass

        try:
            # if value is an atomic value then this fails with TypeError
            len_ = len(value)
        except TypeError:
            # value is an atomic value and will be assigned directly
            type_.in_dll(self.lib, symbol).value = value
        else:
            # value is some kind of list or tuple and will be assigned
            # assuming type_.in_dll() is an array as well
            type_.in_dll(self.lib, symbol)[:len_] = value

    def __str__(self):
        """Nice printable list of module variables
        """
        return "\n".join(
            "%s: %s" % (name, self[name])
            for name in self.fields)

    def set(self, *args, **kwargs):
        """Set many values in one go. Accepts positional and keyword arguments.

         >>> params = ModuleParameter(mytest.lib, [
         >>>     ("number", "mymodule_mp_number_", ct.c_int),
         >>>     ("string", "mymodule_mp_string_", ct.c_char*80), ])
         >>> params.set(5, string="hallo")
        """
        for name, value in izip(self.fields, args):
            self[name] = value
        for name, value in kwargs.iteritems():
            self[name] = value

    def save(self):
        """Return current contents of all module variables (as defined by
        self.fields) to a dict. This dict would be suitable to restore those
        values by passing them to self.set():
         >>> savedVals = params.save()
         >>> print params["string"]
         gugu
         >>> params["string"] = "hallo"
         >>> print params["string"]
         hallo
         >>> params.set(**savedVals)
         >>> print params["string"]
         gugu

        The values of the dict are independet copies -not references- of the
        module internal values. If the values in the module are changed, the
        values returned by this method are not.

        The C- or Fortran-types of the module are being converted to the
        corresponding Python types.
        """
        return dict(
            # ....value: convert to Python type and make a copy
            (name, self[name].value)
            for name in self.fields)


def parseFortranConstants(fileName, firstLineMarker=None, lastLineMarker=None):
    """Fortran constants (at least non-arrays / single values) seem to have no
    symbol in the compiled object file. In order to retreive some constants
    from the Fortran sources this function tries to parse the Fortran source
    code and supply a dictionary of values found thereby.

    Usage:
     >>> vumatParam = parseFortranConstants("source/vumat_param.f")
     >>> eqStepTime = vumatParam["EqStepTime"]

    @param fileName: ... of the Fortran source file

    @param firstLineMarker: Start parsing with the first line that contains
        this string. If not given start with the first line of the file.

    @param lastLineMarker: Stop parsing after the first line after the first
        parsed line that contains this string. If not given then parse all to
        the end of the file.

        Note: The lastLineMarker string will be ignored if it occurs on a
        continuation line.

    @returns: A dict {Fortran identifier: value}.

    @Note: Case matters for Python but not for Fortran. This function will fail
       if the same variable is spelled differently in the Fortran code. I.e.
       ramptime = 10; eqsteptime = RAMPTIME+5 might be ok for Fortran but not
       for this function.
    """

    def readAhead(file_):
        """Generator function yielding complete lines of code concatenating
        continuation lines and removing pure comment lines and empty lines.
        """
        line = ""
        parsing = firstLineMarker is None
        for newline in file_:

            parsing = parsing or (firstLineMarker in newline)
            if not parsing:
                continue

            # ignore empty and comment lines
            if not(newline) or newline[0].upper()=="C":
                continue
            newlinestripped = newline.strip()
            if not(newlinestripped) or newlinestripped[0]=="!":
                continue

            # remove trailing spaces
            newline = newline.rstrip()

            # continuation line
            if newline[5] != " ":

                # remove optional continuation mark at the end of last line
                if line[-1]==newline[5]:
                    line = line[:-1]

                line += newline[6:]
                continue

            # just ordinary line...
            # return the last (stored) line
            if line:
                yield line

            # store current line for next iteration
            line = newline

            # check for lastLineMarker
            if lastLineMarker and (lastLineMarker in newline):
                break

        # also return the last line
        if line:
            yield line

    rexExpo = re.compile(r"([.\d])[Dd](\d)")
    def parseAssignment(paraDict, string):
        """Parse a string like "name = value" and store it as
        paraDict["name"] = value.
        """
        # remove comments
        commentPos = string.find("!")
        if commentPos > -1:
            string = string[:commentPos]

        # Replace array brackets; Fortran arrays -> Python lists
        # If not an array then we have multiple assignments on one line
        # Ignoring the possiblity of multiple array assignments, though
        if "(/" in string:
            string = string.replace("(/", "[").replace("/)", "]")
        else:
            string = string.replace(",", ";")

        # replace .True. and .False.
        string = string.replace(".True.", "True").replace(".False.", "False")

        # replace double precision exponent "D" by e
        string = rexExpo.sub(r"\1E\2", string)

        # execute the string
        try:
            exec string in paraDict
        except Exception, exc:
            raise type(exc)(exc.args[0] + (
                "\nError occurred when parsing the following string:\n%s"
                % string))

    sixSpaces = " "*6
    rexF90Declaration = re.compile(
        r"^(?P<type>.*?)\s*::\s*(?P<assignment>.+=.+)$", re.IGNORECASE)
    rexEndSub = re.compile(r"(?i)^end\s*subroutine")
    rexEndFunc = re.compile(r"(?i)^end\s*function")
    inSubroutine = False
    inFunction = False

    paraDict = dict()
    for line in readAhead(open(fileName, "r")):
        if not line.replace("\t", sixSpaces)[:6] == sixSpaces:
            msg("WARNING: Ignoring illegal line in %s:\n%s"
                % (fileName, line))
            continue

        line = line.strip()
        lineUp = line.upper()

        if rexEndSub.search(line):
            inSubroutine = False
        if rexEndFunc.search(line):
            inFunction = False
        elif line.lower().startswith("subroutine"):
            inSubroutine = True
        elif "function" in line.lower():
            inFunction = True
        if inSubroutine or inFunction:
            continue

        if lineUp.startswith("PARAMETER "):
            parseAssignment(paraDict, line[10:])
            continue

        res = rexF90Declaration.match(line)
        if res:
            parseAssignment(paraDict, res.group("assignment"))
            continue

        msg("Ignoring line without parameter definition:\n%s"
            % line, debugLevel=5)

    del paraDict['__builtins__']
    return paraDict


class PatchStringNotFound(Exception):
    pass

def patchFile(filename, patches):
    """Apply a list of changes to a single file in one piece possibly using
    multiline regular expressions. A backup of the original state will be made.

    The list of changes (argument patches) contains tuples of arguments to the
    re.sub() function. I.e. pattern, replacement and optionally count.
    """
    rexList = [(re.compile(p[0]), list(p[1:])) for p in patches]
    for rex, subArgs in rexList:
        subArgs.insert(1, None)  # placeholder for string to operate on

    origName = "%s.orig%s" % os.path.splitext(filename)
    os.rename(filename, origName)

    data = open(origName, "r").read()

    for cnt, (rex, subArgs) in enumerate(rexList):
        subArgs[1] = data
        data, found = rex.subn(*subArgs)
        if not found:
            raise PatchStringNotFound(
                "ERROR: Did not find the %d-th patch in file %s"
                % (cnt+1, filename))

    open(filename, "w").write(data)


def matrixToFortran(mat):
    """Flatten a python matrix (list of lists) to a vector that will be
    interpreted as matrix by Fortran::

    result(i,j) <- mat[i-1][j-1]

    Note that fortran interates the first index fastest. This is called
    column major order.
    """
    return [x for vec in izip(*mat) for x in vec]

def matrixFromFortran(matFortran, N):
    """Create a python matrix (list of lists) from a Fortran matrix
    matFortran(N,*)::

    matFortran(i,j) -> result[i-1][j-1] , i=1..N

    @param N: number of columns (first index)
    """
    return [ matFortran[i::N] for i in range(N) ]

def timeStepper(dt, tEnd, checkTimes=None, tStart=0.0):
    """Generator function that yields (steptime, dt, checkFlag)-tuples.
    steptime is the steptime at the end of the increment as the vumat expects
    it. dt is generally the same as the given argument except at the end of
    the iteration it might be truncated to exactly match tEnd. checkFlag is
    True if a checkTime has been passed during this increment.

    Usage:
     >>> for (steptime, dt, checkFlag) in timeStepper(
     >>>         dt=0.05,   # in real life e.g. 5E-3
     >>>         tEnd=10.0,
     >>>         checkTimes=[0.1, 0.2, 0.3, 1.9, 2.0, 2.1,] ):
     >>>
     >>>     stressnew, statenew = self.vumat(
     >>>         STEPTIME=steptime,
     >>>         TOTALTIME=steptime+eqStepTime,
     >>>         DT=dt,
     >>>         STRAININC=deltaEps,
     >>>         STRESSOLD=matrixToFortran(stress),
     >>>         )
     >>>     stress = matrixFromFortran(stressnew, self.nBlock)
     >>>
     >>>     if checkFlag:
     >>>         self.assert_(...)

    @param tStart: a float
    @param tEnd: a float
    @param dt: a float
    @param checkTimes: a list of floats.
       If not specified then return checkFlag=True on each iteration.
       If it's an empty list then never check.
    """
    if checkTimes is not None:
        checkTimes = iter( list(checkTimes)+[2*tEnd] )
        try:
            nextCheckTime = checkTimes.next()
        except StopIteration:
            nextCheckTime = 9e99   # never check
    steptime = tStart
    while steptime < tEnd:
        steptime += dt
        if steptime > tEnd:
            dt = tEnd - (steptime-dt)
            steptime = tEnd

        # checks only at certain times
        if checkTimes is None:
            checkFlag = True
        else:
            if steptime >= nextCheckTime:
                while nextCheckTime<=steptime:
                    nextCheckTime = checkTimes.next()
                checkFlag = True
            else:
                checkFlag = False

        yield steptime, dt, checkFlag

def checkVector(vec, reference, minRefNorm=1E-10):
    """Calculate the relative error norm squared of vec compared to reference.

    @param minRefNorm: To get the relative error we have to divide by the
        squared norm of reference. In order to prevent division by zero you can
        specify a positive value to use instead should this squared norm be
        smaller.
    """
    absSqRef = max(minRefNorm, sum(x**2 for x in reference))
    errSq = sum((x-y)**2 for x,y in izip(vec, reference)) / absSqRef
    return errSq
