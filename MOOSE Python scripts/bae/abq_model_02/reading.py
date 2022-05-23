"""abq_model_02 - service classes for reading the Abaqus input file.
Those classes don't actually hold the bulk data they are just for interpreting
the input file. Not supposed to be used directly but only through the class
abq_model_02.Model.


Adding new structures to read from the input file:
==================================================

That's quite easy. For now see the other Block...classes and ask Gero.
... and please update the documentation.
"""

# no __version__ info here, it's in __init__.py

_todo_ = """
ideas:
======

- to speed up reading:
  . read data block until next command filtering comments with Python, then use
    pandas or numpy to read the bulk data.
  . use Cython or Pyrex to write a compiled extendion module for reading
  . both options must be optional, because RhinoPython will not work with
    neither of them.
- for further ideas do the following:
   >>> import cProfile
   >>> from bae.abq_model_01 import Model
   >>> m = Model()
   >>> cProfile.run('m.read("somebigmodel.inp")')
- Very frequent Abaqus cards like our gravity amplitude definitions and *DLOAD
  cards may suffer from expensive object initialization. To circumvent:
  In the Model-object at initialization, create one static Block-object
  for each Block-type. Blocks to be ignored are of type SkipBlock (or
  similar but specialized). Then those instances don't have to be created
  and destroyed while reading. Attributes like the options-warnflag may
  be changed after creating the Model-object before actually reading.
"""


import re
import warnings
from itertools import izip

from bae.log_01 import msg

from bae.misc_01 import isIronPython

from bae.future_01 import NumpyVersion

try:
    import numpy as np
except ImportError:
    np = None

try:
    import pandas as pd
    pandasVersion = NumpyVersion(pd.__version__)
except:
    pd = None
    pandasVersion = None


#DEBUG
# np = None
# pd = None

from StringIO import StringIO



from .container import \
    NodesDict, ElementNodesDict, Material, \
    ConnectorBehavior, Orientation, Amplitude, \
    initCondTypeMultLine, initCondTypeToClass

#==============================================================================


class ItemAlreadyThere(ValueError):
    pass
class ReadError(Exception):
    pass
class UnknownOption(ReadError):
    pass
class UnionSetInvalidItem(Exception):
    pass
class UnknownElType(Exception):
    pass

cmdOptionTypeToClass = dict()


######################################################
###    Command Options service classes
######################################################

class CmdOption(object):
    """Descriptor class that finds one command option in the ABAQUS input file,
    something like NSET=myset or GENERATE, and creates an appropriate
    attribute/member in the owner class instance.

    Uses a regular expression to find a certain kind of command option
    in the string argument passed to the member function find_in.

    If the member warnflag is set to True, each occurence of this option in the
    input file will trigger a warning in the log.
    """
    def __init__(self, keyword, nbSigLetters=-1, warnflag=False):
        """
        The first argument to the constructor "keyword" states the name of the
        option in the ABAQUS input file (e.g. generate) as well as the name of
        the attribute in the owner class instance. While the option name is case
        insensitive, the attribute name of the owner class is not, therefore
        it's suggested to use lowercase letters.

        The second argument nbSigLetters (an int) states how many letters of the
        option keyword are significant. Only those significant letters are
        checked by the regular expression.
        """
        self.attr_name = keyword
        self.nbSigLetters = nbSigLetters
        self.warnflag = warnflag

    def init_member_attr(self, instance):
        setattr(instance, self.attr_name, self.default)


class CmdOptionNoArgs(CmdOption):
    """Descriptor class that finds one command option in the ABAQUS input file,
    something like NSET=myset or GENERATE, and creates an appropriate
    attribute/member in the owner class instance.

    Uses a regular expression to find a certain kind of command option
    in the string argument passed to the member function find_in.

    If the member warnflag is set to True, each occurence of this option in the
    input file will trigger a warning in the log.

    This is for options without arguments e.g. GENERATE.
    The value of the corresponding attribute is True (found) or False
    """

    def __init__(self, keyword, nbSigLetters=-1, warnflag=False):
        r"""
        The first argument to the constructor "keyword" states the name of the
        option in the ABAQUS input file (e.g. generate) as well as the name of
        the attribute in the owner class instance. While the option name is case
        insensitive, the attribute name of the owner class is not, therefore
        it's suggested to use lowercase letters.

        The second argument nbSigLetters (an int) states how many letters of the
        option keyword are significant. Only those significant letters are
        checked by the regular expression.
        """
        CmdOption.__init__(self, keyword, nbSigLetters, warnflag)
        self.default = False

    def find_in(self, instance, arg, model):
        "set option if arg contains it. string values are converted uppercase"
        arg = arg.strip()
        if (self.attr_name.upper().startswith(arg.upper())
            and (self.nbSigLetters<0
                 or len(arg)>=self.nbSigLetters)):
            setattr(instance, self.attr_name, True)

            # issue a warning, if warnflag is True
            if self.warnflag:
                msg("WARNING: Option %s for command *%s occured on line %d."
                    % (self.attr_name.upper(), instance.name, model.linecnt))
            return True
        else:
            return False
cmdOptionTypeToClass["NoArgs"] = CmdOptionNoArgs


class CmdOptionWiArg(CmdOption):
    """Base class for command options with arguments.
    """
    generalRegExp = re.compile(r"(\S+)\s*=\s*(.*)", re.IGNORECASE)

    def __init__(self, keyword, nbSigLetters=-1, warnflag=False):
        r"""
        The first argument to the constructor "keyword" states the name of the
        option in the ABAQUS input file (e.g. elset) as well as the name of the
        attribute in the owner class instance. While the option name is case
        insensitive, the attribute name of the owner class is not, therefore
        it's suggested to use lowercase letters.

        The second argument nbSigLetters (an int) states how many letters of the
        option keyword are significant. Only those significant letters are
        checked by the regular expression.
        """
        CmdOption.__init__(self, keyword, nbSigLetters, warnflag)
        self.regex = re.compile(keyword + r"\s*=\s*(.*)", re.IGNORECASE)
        self.default = None

    def find_in(self, instance, arg, model, conversionFunc):
        "set option if arg contains it. string values are converted uppercase"
        set_fnd = self.regex.match(arg)
        if set_fnd:
            res = set_fnd.group(1)
        else:
            set_fnd = self.generalRegExp.match(arg)
            if set_fnd:
                keyFromFile, res = set_fnd.groups()
                set_fnd = (
                    self.attr_name.upper().startswith(keyFromFile.upper())
                    and self.nbSigLetters>0
                    and len(keyFromFile)>=self.nbSigLetters)

        if set_fnd:
            try:
                res = conversionFunc(res)
            except ValueError, exc:
                raise ReadError(
                    "The value %s for the option %s for command *%s on line %d"
                    " can not be converted to the correct type\n%s"
                    % (res, self.attr_name.upper(), instance.name,
                       model.linecnt, exc.args[0]))
            setattr(instance, self.attr_name, res)

            # issue a warning, if warnflag is True
            if self.warnflag:
                msg("WARNING: Option %s for command *%s occured on line %d."
                    % (self.attr_name.upper(), instance.name, model.linecnt))
        return bool(set_fnd)


class CmdOptionOneName(CmdOptionWiArg):
    """Descriptor class that finds one command option in the ABAQUS input file,
    something like NSET=myset, and creates an appropriate attribute/member in
    the owner class instance.

    Uses a regular expression to find a certain kind of command option
    in the string argument passed to the member function find_in.

    This is for options with one argument, e.g. ELSET=myset.
    Its value is None (keyword not found) or the options argument string

    If the member warnflag is set to True, each occurence of this option in the
    input file will trigger a warning in the log.
    """
    def find_in(self, instance, arg, model):
        return CmdOptionWiArg.find_in(
            self, instance, arg, model, str.upper)
cmdOptionTypeToClass["OneName"] = CmdOptionOneName


class CmdOptionOneInt(CmdOptionWiArg):
    """Descriptor class that finds one command option in the ABAQUS input file,
    something like NSET=myset, and creates an appropriate attribute/member in
    the owner class instance.

    Uses a regular expression to find a certain kind of command option
    in the string argument passed to the member function find_in.

    This is for options with one argument, e.g. DELETE=10.
    Its value is None (keyword not found) or the options argument converted to
    an integer.

    If the member warnflag is set to True, each occurence of this option in the
    input file will trigger a warning in the log.
    """
    def find_in(self, instance, arg, model):
        return CmdOptionWiArg.find_in(
            self, instance, arg, model, int)
cmdOptionTypeToClass["OneInt"] = CmdOptionOneInt


class CmdOptionOneFloat(CmdOptionWiArg):
    """Descriptor class that finds one command option in the ABAQUS input file,
    something like NSET=myset, and creates an appropriate attribute/member in
    the owner class instance.

    Uses a regular expression to find a certain kind of command option
    in the string argument passed to the member function find_in.

    This is for options with one argument, e.g. ALPHA=0.5.
    Its value is None (keyword not found) or the options argument string

    If the member warnflag is set to True, each occurence of this option in the
    input file will trigger a warning in the log.
    """
    def find_in(self, instance, arg, model):
        return CmdOptionWiArg.find_in(
            self, instance, arg, model, float)
cmdOptionTypeToClass["OneFloat"] = CmdOptionOneFloat


class OptionList(list):
    """List of L{CmdOption} - objects

    A list of -at least- (keyword, type)-tuples. keyword is the Abaqus keyword.
    Type is one of the keys in the global dict L{cmdOptionTypeToClass} with the
    corresponding L{CmdOption}-subclass as value.

    A third optional element of this tuple can be the number of significant
    letters of the option. An option identifier can be abbreviated to at least
    this number of letters. The default -1 indicates no abbreviation possible.

    A fourth optional element of this tuple can be the string "warning". If so
    the attribute warnflag of the correspondig CmdOption-object will be set to
    True and a warning will be issued, when this option occurs in the input
    file.
    """
    def __init__(self, key_type_list):
        for opt_dat in key_type_list:
            # mandatory items in the tuples
            key, typ = opt_dat[0:2]
            # optional item in the tuples: nb of significant letters
            try:
                nbSigLetters = opt_dat[2]
                assert type(nbSigLetters) is int, "Fatal Error in OptionList."
            except (IndexError, ValueError):
                nbSigLetters = -1
            # optional item in the tuples: warn flag
            try:
                warnflag = (opt_dat[3][:4].lower()=='warn')
            except IndexError:
                warnflag = False
            # add option to self
            cmdOptionClass = cmdOptionTypeToClass[typ]
            self.append(cmdOptionClass(key, nbSigLetters, warnflag))

    def getNameList(self):
        "Return a list of the names of the options in this list."
        return [opt.attr_name for opt in self]

blockdict = dict()


######################################################
###    Blocks: (abstract) super classes
######################################################

class Block(object):
    """base class for all data blocks

    Any child class should call Block.__init__

    Any child class that is to be used for reading must implement a
    readDataLines method that reads data lines from the input file.
    It returns the first line it can't parse as data line, usually containing
    the following command. Plus the number of lines read.

    If the child class is a subblock of a superblock, its readDataLines method
    must accept an additional argument superblock stating the superblock-object.

    Any block-class that needs special treatment for skipping this block (e.g.
    a *STEP-block does not end until its *END STEP counterpart), has to
    redefine the "SkipThisBlock" class member as a reference to a previously
    defined specialized class (possibly derived from SkipBlock).
    (example: class SkipMyBlock: ....

    then: class MyBlock(Block): SkipThisBlock = SkipMyBlock;def __init__ ....)
    The Block.SkipThisBlock class attribute is assigned after the definition of
    the class "SkipBlock".
    """

    re_key = re.compile(r"\*(.*?)\s*(,.*)*$")
    options = list()

    def __init__(self, name, line, model, options=True):
        """__init__ has to be called from any child instance
        if options-flag is True, initializes and reads the command options
        """
        self.name = name.upper()
        if options:
            for opt in self.options:
                opt.init_member_attr(self)
            self.read_options(line, model)

    @classmethod
    def read_key(self, line):
        """
        extract the abaqus command key word from line

        line has to be stripped of spaces from the right end already:
        ... line=line.rstrip() already performed

        returns None or False if line is empty
        """
        keymatch = self.re_key.match(line)
        if not(keymatch):
            key = "???"
            raise ReadError(
                "Invalid beginning of block:\n%s" % line)
        key = keymatch.group(1).upper().replace(" ", "")
        return key

    def read_options(self, line, model):
        """
        extracts command options from the first line of a ABAQUS command
        block and assigns to appropriate member attributes"""
        for arg in line.split(",")[1:]:
            arg = arg.strip()  # ignore white space
            if not arg:
                continue  # ignore empty options (double comma or so)
            found = False
            for opt in self.options:
                found = opt.find_in(self, arg, model)
                if found:
                    break
            # not a known option
            if not found:
                raise UnknownOption("Unknown option of the *%s-command:\n%s"
                                    % (self.name, arg))

    def readDataLines(self, inp, model):
        "read data lines until next block starts"
        raise NotImplementedError(
            "ERROR: readDataLines method not implemented for type %s"
            % type(self))
        # typical implementation follows:
        cnt = 0
        for cnt, line in self.filteredLines(inp, 1):
            if line[0]=="*":
                break
            data = line.split(",")
            #... process data lines
            line = None  # in case inp hits EOF: return empty line
        return line, cnt

    @staticmethod
    def createListFromGenerate(settype, line):
        """create list of nodes / elements for an *ELSET / *NSET command
        with the GENERATE option

        @param settype: either "ELSET" or "NSET"
        @param line: one the data line from the input file
        """
        data = line.split(",")
        if len(data) < 2:
            raise ReadError(
                "*%s,GENERATE has less than 2 values in data line."
                % settype)
        try:
            start = int(data[0])
            stop = int(data[1])
        except ValueError:
            raise ReadError(
                "*%s,GENERATE has non-integers as start ('%s') or stop"
                " ('%s') values in the data line."
                % (settype, data[0], data[1]))
        if len(data) > 2:
            try:
                inc = int(data[2])
            except ValueError:
                raise ReadError(
                    "*%s,GENERATE has the non-integer '%s' as increment"
                    " value in the data line." % (settype, data[2]))
        else:
            inc = 1
        return range(start, stop+1, inc)

    def getOptionsDict(self):
        "Return a dictionary {option name: value}"
        return dict([(optName, getattr(self, optName))
                     for optName in self.options.getNameList()])

    @staticmethod
    def filteredLines(inp, cnt):
        """return an iterator of (counter, lines) tuples from inp with comment
        and empty lines filtered out
        """
        return ((cnt+i, line) for i, line in enumerate(inp)
                if not(len(line)>=2 and line[0:2]=="**") and line.strip())


class SkipBlock(Block):
    """class for skipping a block that is to be skipped, because
    it is not in the list recognized_blocks or because it is a block of
    unknown type.

    If a block contains subblocks (e.g. *STEP lasts until *END STEP)
    this blockclass has to provide its own specialized SkipBlock class
    with a special readDataLines function. The skipblock-method of this block
    must then return a new object of this specialized SkipBlock class.
    """
    def __init__(self, line, model):
        Block.__init__(self, line.split(",")[0].strip().upper(),
                       line, model, options=False)

    def readDataLines(self, inp, model):
        "skip lines until next block starts"
        cnt = -1
        line = None
        for cnt, line in enumerate(inp):
            if line[0]=="*" and line[1]!="*":
                break
            line = None  # in case inp hits EOF: return empty line
        return line, cnt+1


# initialize default SkipThisBlock - class member as class SkipBlock
Block.SkipThisBlock = SkipBlock


class SuperBlock(Block):
    """base class for super blocks like *MATERIAL that have dependent subblocks
    following like *ELASTIC

    Any child class should call Block.__init__

    see description of class Block
    """
    validSubBlocks = dict()

    @staticmethod
    def SkipThisBlock(line, model):
        abq_key = Block.read_key(line)
        blockName = blockdict[abq_key].__name__
        raise Exception(
            "Program ERROR: In abq_model_02.reading there is missing a"
            " statement like the following: %s.SkipThisBlock = Skip%s."
            " And probably the definition of the class Skip%s as well.\n"
            "For the moment you might get away with including '%s' into"
            " the recognizedBlocks argument of Model.read() or not specifying"
            " this argument at all."
            % (blockName, blockName, blockName, abq_key))
# add lines like this for each SuperBlock-child class (here BlockXXXX)
# class SkipBlockXXXX(SkipSuperBlock):
#     validSubBlocks = BlockXXXX.validSubBlocks
# BlockXXXX.SkipThisBlock = SkipBlockXXXX


class SkipSuperBlock(Block):
    """class for skipping a super block and all its sub blocks
    """

    def __init__(self, line, model):
        Block.__init__(self, line.split(",")[0].strip().upper(),
                       line, model, options=False)

    def readDataLines(self, inp, model):
        "skip lines until next super-block starts"
        cnt = -1
        line = None
        for cnt, line in enumerate(inp):
            if line[0]=="*" and line[1]!="*":
                key = self.read_key(line)
                if key not in self.validSubBlocks:
                    break
            line = None  # in case inp hits EOF: return empty line
        return line, cnt+1

# add lines like this for each SuperBlock-child class (here BlockXXXX)
# class SkipBlockXXXX(SkipSuperBlock):
#     validSubBlocks = BlockXXXX.validSubBlocks
# BlockXXXX.SkipThisBlock = SkipBlockXXXX


######################################################
###    Blocks: general
######################################################

class BlockHeading(Block):
    "Class for reading the *HEADING - block"

    def __init__(self, line, model):
        Block.__init__(self, "HEADING", line, model, options=False)

    def readDataLines(self, inp, model):
        "read lines until next block starts"
        cnt = -1
        line = None
        for cnt, line in enumerate(inp):
            if line[0]=="*" and line[1]!="*":
                break
            line = None  # in case inp hits EOF: return empty line
        return line, cnt+1
blockdict["HEADING"] = BlockHeading


######################################################
###    Blocks: model definition - nodes
######################################################

class BlockNode(Block):
    """class for reading of a *NODE - block in the input file

    creates:
    model.nodeCoords (NodesDict nodenumber -> list of coords)
    possibly model.nset (NsetsDict nset-name -> set of nodenumbers)

    if optional argument onlySet==True: read only node numbers and create nset
    if optional argument noSets==True: don't import nsets
    """

    options = OptionList([("nset", "OneName"),])

    def __init__(self, line, model, onlySet=False, noSets=False):
        Block.__init__(self, "NODE", line, model)
        self.onlySet = onlySet
        if noSets and self.nset:
            self.nset = None

        # make sure nodeCoords-attribute exists
        # should not be necessary
        if not(hasattr(model, "nodeCoords")):
            model.nodeCoords = NodesDict()

    def readDataLinesContNoNset(self, inp, nodeCoords):
        """Read a contigous block of n,x,y,z data lines

        @param nodeCoords: list to be filled with tuples
           (node number, list-of-coords)
        """
        # No benefit from using pandas or numpy here:
        #    - pandas.read_csv would require an additional StringIO-object
        #    generation to return line and cnt correctly
        #    - np.fromstr seems to be slightly slower than mapping float
        for cnt, line in enumerate(inp):
            try:
                #mapping float to all and then int(first) is slightly faster
                #than to do this do indexed line.split()
                vals = map(float, line.split(','))
                nodeCoords[int(vals[0])] = vals[1:]
                line = ''
            except ValueError as e:
                if not line or not line.strip() or line[:2] == '**':
                    #skip comments and empty lines
                    continue
                elif line[0] == '*':
                    #new keyword
                    break
                else:
                    raise ReadError("Can't read/process NODE line %d: %s"
                                    % (cnt,str(e)))

        return cnt, line


    def readDataLinesContNsetCoords(self, inp, nset, nodeCoords):
        """Read a contigous block of n,x,y,z data lines

        @param nodeCoords: list to be filled with tuples
           (node number, list-of-coords)
        """
        cnt = 0
        while 1:

            # try to read a contiguous block of n,x,y,z data lines
            line = None
            cnt2 = -1
            try:
                for cnt2, line in enumerate(inp):
                    n, x, y, z = line.split(",")
                    n = int(n)
                    nodeCoords[n] = [float(x), float(y), float(z)]
                    nset.add(n)
                    line = None
            except ValueError:
                pass
            cnt += (cnt2+1)

            # stop if not comment and line not empty
            if not(line) or (line[0:2]!="**" and line.strip()):
                break

        # found something else: new command or dataline other than n,x,y,z
        return cnt, line

    def readDataLinesContNsetOnly(self, inp, nset):
        """Read a contigous block of n,x,y,z data lines

        @param nset: list to be filled with tuples
           (node number, list-of-coords)
        """
        cnt = 0
        while 1:

            # try to read a contiguous block of n,x,y,z data lines
            line = None
            cnt2 = -1
            try:
                for cnt2, line in enumerate(inp):
                    nset.add(int(line.split(",", 1)[0]))
                    line = None
            except (ValueError, IndexError):
                pass
            cnt += (cnt2+1)

            # stop if not comment and line not empty
            if not(line) or (line[0:2]!="**" and line.strip()):
                break

        # found something else: new command or dataline other than n,x,y,z
        return cnt, line

    def readDataLines(self, inp, model):
        "read *NODE data lines until next block starts"

        if self.nset:
            try:
                nset = model.nset[self.nset]
            except KeyError:
                nset = set()
                assignNset = True
            else:
                assignNset = False

            # read the data
            if self.onlySet:
                cnt, line = self.readDataLinesContNsetOnly(inp, nset)
            else:
                cnt, line = self.readDataLinesContNsetCoords(
                    inp, nset, model.nodeCoords)

            # only assing the nset if there are new nodes...
            if assignNset and nset:
                model.nset[self.nset] = nset
        else:
            # read the data
            cnt, line = self.readDataLinesContNoNset(inp, model.nodeCoords)

        # if line is still a data line:
        # ...read data lines with more flexibility
        while line and line[0]!="*":

            data = line.split(",")
            try:
                nodenum = int(data[0])
            except ValueError:
                raise ReadError("Invalid node number %s." % data[0])
            if self.nset:
                try:
                    model.nset[self.nset].add(nodenum)
                except KeyError:
                    model.nset[self.nset] = set((nodenum,))

            if not self.onlySet:
                if len(data) > 7:
                    errstr = (
                        "more than 6 node coordinates. Whole data line:\n%s"
                        % line)
                    raise ReadError(errstr)
                try:
                    coords = map(float,data[1:7])
                except ValueError:
                    raise ReadError(
                        "Invalid node coordinates %s." % str(data[1:7]))
                model.nodeCoords[nodenum] = coords

            # get next line
            while 1:
                try:
                    line = inp.next()
                except StopIteration:
                    line = None
                    break
                cnt += 1

                # ignore comments and empty lines:
                # stop reading if not comment and line not empty
                if line[0:2]!="**" and line.strip():
                    break

        # fini
        return line, cnt
blockdict["NODE"] = BlockNode


class BlockNset(Block):
    """class for reading a *NSET - block in the input file

    creates:
    model.nset (NsetsDict nset-name -> list of nodenumbers)
    """

    options = OptionList((
        ("nset", "OneName"),
        ("generate", "NoArgs"),
        ("instance", "OneName", -1, "warning"),
        ("internal", "NoArgs", -1, "warning"),
        ("elset", "OneName"),
    ))

    def __init__(self, line, model):
        Block.__init__(self, "NSET", line, model)

        # check nset-name
        if not self.nset:
            raise ReadError("*NSET command without NSET-option.")

        # elset Option
        # contains the elset name from *NSET, NSET=..., ELSET=...
        if self.elset:

            # 1. generate option not suitable in common with elset option
            if self.generate:
                msg("WARNING: *NSET command with ELSET and GENERATE option"
                    " on line %d, NSET %s. The GENERATE option and all"
                    " datalines will be ignored."
                    % (model.linecnt, self.nset))

            # get element numbers (-> elems) and...
            # 2. check: is the elset defined?
            try:
                elems = model.elset[self.elset]
            except KeyError:
                msg("WARNING: Elset %s not defined for creating nset %s by a"
                    " command like this: *NSET, NSET=%s, ELSET=%s on line %d."
                    " Ignoring this *NSET - card."
                    % (self.elset, self.nset, self.nset, self.elset,
                       model.linecnt))
                return

            # 3. check: do we have element definitions at all?
            if not model.elNodes:
                msg("WARNING: No element definitions available for creating"
                    " nset %s by a command like this: *NSET, NSET=%s, ELSET=%s"
                    " on line %d. Ignoring this *NSET - card."
                    % (self.nset, self.nset, self.elset, model.linecnt))
                return

            # 4. check: are all elements defined for the given elset?
            missingElems = elems.difference(model.elNodes)
            if missingElems:
                msg("WARNING: %d of %d elements in elset %s lack a definition"
                    " when creating nset %s by a command like this:"
                    " *NSET, NSET=%s, ELSET=%s on line %d. The resulting nset"
                    " might not be complete."
                    % (len(missingElems), len(elems), self.elset,
                       self.nset, self.nset, self.elset, model.linecnt))

            # cretate nset from elems
            model.nset[self.nset] = model.getConnectedNodes(elems)

    def readDataLines(self, inp, model):
        "read *NSET data lines until next block starts"
        cnt = 0
        line = None
        if self.elset:
            for cnt, line in self.filteredLines(inp, 1):
                if line[0]=="*":
                    break
                line = None
            if cnt>0 and not line:
                msg("WARNING: *NSET, NSET=%s, ELSET=%s appears to have a data"
                    " line on line %d. It will be ignored."
                    % (self.nset, self.elset, model.linecnt))

        else:
            try:
                nset = model.nset[self.nset]
            except KeyError:
                nset = set()
                assignNset = True
            else:
                assignNset = False

            if self.generate:
                for cnt, line in self.filteredLines(inp, 1):
                    if line[0]=="*":
                        break
                    nset.update(self.createListFromGenerate("NSET", line))
                    line = None  # in case inp hits EOF: return empty line

            else:
                for cnt, line in self.filteredLines(inp, 1):
                    if line[0]=="*":
                        break
                    data = line.split(",")

                    # try all items are numbers
                    try:
                        nset.update(int(d) for d in data)
                    except ValueError:
                        pass
                    else:
                        line = None  # in case inp hits EOF: return empty line
                        continue

                    # at least some items are not a number
                    for d in data:
                        try:
                            nset.add(int(d))
                        except ValueError:
                            d = d.strip().upper()
                            if d:
                                try:
                                    nset.update(model.nset[d])
                                except KeyError:
                                    msg("WARNING: NSET %s specified as part of"
                                        " NSET %s, but %s is not (yet) defined"
                                        " in the model." % (d, self.nset, d))
                        # item is an empty string (or similar)
                        # ...ignore this item
                        continue

                    # line processed, continue
                    line = None  # in case inp hits EOF: return empty line

            # store the new set if applicable
            if assignNset and nset:
                model.nset[self.nset] = nset
        return line, cnt
blockdict["NSET"] = BlockNset


######################################################
###    Blocks: model definition - elements
######################################################

class BlockElement(Block):
    """class for reading of an *ELEMENT block in the input file

    creates:
    model.elNodes (ElementNodesDict: element number -> list of node numbers)
    model.elType (dict: element number -> element type)
    model.typeEl (dict: element type -> set of element numbers)
    possibly model.elset (ElsetsDict: elset name -> list of element numbers)

    If you want to change this BlockElement-implementation: It is probably
    tempting to use model.update_elem()-method. However this would probably be
    slower than this.

    If optional argument onlySet==True: read only element numbers and create
    nset. If optional argument noSets==True: don't import elsets.
    """

    options = OptionList((
        ("type", "OneName"),
        ("elset", "OneName")))

    def __init__(self, line, model, onlySet=False, noSets=False):
        Block.__init__(self, "ELEMENT", line, model)
        self.onlySet = onlySet
        if noSets and self.elset:
            self.elset = None

        # make sure elNodes-, elType- and typeEl-attributes exist
        if not(hasattr(model, "elNodes")):
            model.elNodes = ElementNodesDict()
        if not(hasattr(model, "elType")):
            model.elType = dict()
        if not(hasattr(model, "typeEl")):
            model.typeEl = dict()

        self.data = list()  # temporary list of element numbers

    def readDataLines(self, inp, model):
        "read *ELEMENT data lines until next block starts"

        #1) buffering to StringIO
        strIO = StringIO()

        # try to read a contiguous block of data lines
        cnt = -1
        line = None
        for cnt, line in enumerate(inp):
            lineS = line.strip() #remove trailing linebreaks, spaces
            if lineS and not lineS[0] == '*':
                # not empty, no comment, no command => is dataline
                if lineS[-1] == ',':
                    # remove linebreak
                    strIO.write(lineS)
                else:
                    # resume linebreak
                    strIO.write(line)
            elif lineS and lineS[:2] != '**':
                #new keyword
                break
            line = None

        # end of contiguous block
        cnt += 1

        # check last line
        if strIO.getvalue()[-1] == ',':
            msg("WARNING: Element data not complete on line %d. Line ends with"
                " trailing comma. Will be ignored." % model.linecnt+cnt)
            return line, cnt


        #3) store data in model
        if strIO.tell():

            # reset strIO
            strIO.seek(0)

            # find elset
            if self.elset:
                try:
                    elset = model.elset[self.elset]
                except KeyError:
                    elset = set()
                    assignElset = True
                else:
                    assignElset = False

            # convert to numbers and update elNodes, elType, typeEl
            if self.onlySet:  # implies that self.elset is set!
                try:
                    data = (int(row.split(',', 1)[0]) for row in strIO)
                    elset.update(data)
                except ValueError:
                    badLine = strIO.getvalue()
                    badPosE = strIO.tell()
                    badPosA = badLine[:badPosE].rfind('\n')+1
                    badLine = badLine[badPosA:badPosE]
                    raise ReadError("Invalid element number %s." % badLine)

            else:
                try:
                    if pd:
                        # 0. pandas - fastest
                        intRows = pd.read_csv(strIO, header=None,
                                              dtype=np.uint32)
                        elems = intRows.iloc[:,0]
                        if pandasVersion >= '0.24.0':
                            nodes = intRows.iloc[:,1:].to_numpy().tolist()
                        elif pandasVersion >= '0.22.0':
                            nodes = intRows.iloc[:,1:].values.tolist()
                        else:
                            nodes = intRows.iloc[:,1:].as_matrix().tolist()
                        data = izip(elems, nodes)

                    elif np:
                        # 1. np.fromstring fastest for int
                        # no benefit using map here
                        intRows = np.array(
                            [np.fromstring(row, dtype=int, sep=',')
                             for row in strIO])
                        # 3. med fast - inbuild str-conversation
                        # intRows = np.array(
                        #         [np.array(row.split(','), dtype=int)
                        #         for row in strIO])
                        # 4. np.genfromtxt - very slow
                        # nNodes = len(strIO.next().split(','))
                        # strIO.seek(0)
                        # intRows = np.genfromtxt(strIO,
                        #                        dtype=[('v%d'%d,'<u4')
                        #                        for d in range(nNodes)],
                        #                        delimiter=',')
                        # intRows = intRows.view(np.uint32).reshape(
                        #        intRows.shape + (-1,))
                        # 5. np.loadtxt - slowest, by far
                        # intRows = np.loadtxt(strIO, delimiter=',', dtype=int)

                        elems = intRows[:,0]
                        nodes = intRows[:,1:].tolist()
                        data = izip(elems, nodes)

                    else:
                        intRows = (map(int, row.split(',')) for row in strIO)
                        data = [(e[0], e[1:]) for e in intRows]
                        elems = [e[0] for e in data]

                except Exception as e:
                    badLine = strIO.getvalue()
                    badPosE = strIO.tell()
                    badPosA = badLine[:badPosE].rfind('\n')+1
                    badLine = badLine[badPosA:badPosE]
                    raise ReadError(str(e) + "\n" +
                                    "Invalid element or node number: %s."
                                    % badLine)

                # store elements in model
                if isIronPython():
                    # Rhinos ironPython has a corrupt/buggy implementation of
                    # hashtables --> dict.update can be extremly slow if the
                    # dictionary allready has items and the new items have
                    # clumsy keys (e.g. is ranges from 0 to 1Mio and appended
                    # ranges from 10Mio to ...). This could slow down reading
                    # for hours. The workaround is to create a  dictionary with
                    # continuous keys (min(of all) to max(of all)), updating
                    # this one with the existing and the new data and delete
                    # all 'None'-set items. The elNodes dict will be replaced
                    # by this dictionary.
                    # In worst case this workaround slowes down reading for
                    # a few seconds, in best case it saves hours.
                    dataDict = dict(data)
                    msg('storing %d elements to model.elnodes' % len(data))
                    if len(model.elNodes):
                        maxId = max(max(dataDict.keys()), max(model.elNodes))
                        minId = min(min(dataDict.keys()), min(model.elNodes))
                        tmpDict = type(model.elNodes)().fromkeys(
                            range(minId, maxId))
                        dict.update(tmpDict, model.elNodes)
                        dict.update(tmpDict, data)
                        for k,v in tmpDict.items():
                            if v is None:
                                del tmpDict[k]
                        model.elNodes = tmpDict
                    else:
                        dict.update(model.elNodes, dataDict)
                    msg('... done')
                else:
                    dict.update(model.elNodes, data)

                model.elType.update(dict.fromkeys(elems, self.type))

                model.typeEl[self.type].update(elems)
                if self.elset:
                    elset.update(elems)

            # store elset
            if self.elset and assignElset and elset:
                model.elset[self.elset] = elset

        # to make sure that strIO will die
        strIO.close()
        del strIO

        # end
        return line, cnt
blockdict["ELEMENT"] = BlockElement


class BlockElset(Block):
    """class for reading an *ELSET - block

    creates:
    model.elset (ElsetsDict: elset name -> list of element numbers)
    """

    options = OptionList((
        ("elset", "OneName"),
        ("generate", "NoArgs"),
        ("instance", "OneName", -1, "warning"),
        ("internal", "NoArgs", -1, "warning")))

    def __init__(self, line, model):
        Block.__init__(self, "ELSET", line, model)

        # check: option ELSET mandatory
        if not self.elset:
            raise ReadError("*ELSET command without ELSET-option.")

    def readDataLines(self, inp, model):
        "read *ELSET data lines until next block starts"
        cnt = 0
        line = None
        try:
            elset = model.elset[self.elset]
        except KeyError:
            elset = set()
            assignElset = True
        else:
            assignElset = False

        if self.generate:
            for cnt, line in self.filteredLines(inp, 1):
                if line[0]=="*":
                    break
                elset.update(self.createListFromGenerate("ELSET", line))
                line = None  # in case inp hits EOF: return empty line

        else:
            for cnt, line in self.filteredLines(inp, 1):
                if line[0]=="*":
                    break
                data = line.split(",")

                # try all items are numbers
                try:
                    elset.update(int(d) for d in data)
                except ValueError:
                    pass
                else:
                    line = None  # in case inp hits EOF: return empty line
                    continue

                # at least some items are not a number
                for d in data:
                    try:
                        elset.add(int(d))
                    except ValueError:
                        d = d.strip().upper()
                        if d:
                            try:
                                elset.update(model.elset[d])
                            except KeyError:
                                msg("WARNING: ELSET %s specified as part of"
                                    " ELSET %s, but %s is not (yet) defined in"
                                    " the model." % (d, self.elset, d))
                    # item is an empty string (or similar)
                    # ...ignore this item
                    continue

                # line processed, continue
                line = None

        # store the new set if applicable
        if assignElset and elset:
            model.elset[self.elset] = elset
        return line, cnt
blockdict["ELSET"] = BlockElset


######################################################
###    Blocks: model definition - boundary conditions
######################################################

class BlockBoundary(Block):
    """Class for reading of a *BOUNDARY block in the input file

    This is not very useful yet, since it's not recognized whether the BC's
    are set within a certain step.

    updates:
     - model.boundary
     - model.node_bc

    @Note: OP-option as in *BOUNDARY, OP=NEW is being ignored so far, a warning
    is issued.
    """

    options = OptionList((
        ("type", "OneName"),
        ("amplitude", "OneName", 3),
        ("op", "OneName", -1, "warning"),
        ("submodel", "NoArgs", -1, "warning"),
        ("step", "OneInt", -1, "warning"),
    ))

    def __init__(self, line, model):
        Block.__init__(self, "BOUNDARY", line, model)

    def readDataLines(self, inp, model):
        "read *BOUNDARY data lines until next block starts"
        cnt = 0
        for cnt, line in self.filteredLines(inp, 1):
            if line[0]=="*":
                break
            data = line.split(",")

            # which nodes / nsets
            try:
                key = int(data[0])
            except ValueError:
                key = data[0].strip()

            # which flavour of BC definition?
            try:
                dof_1 = int(data[1])
            except ValueError:
                # zero-valued bc's using the "type" format (model data only)
                bc = data[1].strip().upper()

                # check: no option
                if not all([getattr(self, opt) is None
                            for opt in self.options.getNameList()]):
                    msg('WARNING: *BOUNDARY command with options and "type"'
                        ' format <%s> on line %d.'
                        'All command options will be ignored.'
                        % (bc, model.linecnt))
            else:
                # using the "direct" format
                bc = [dof_1, dof_1, 0.0]
                if len(data)>2:
                    bc[1] = int(data[2])
                if len(data)>3:
                    bc[2] = float(data[3])

            # store data in dictionary
            model.boundary.insert(
                key, bc, type=self.type, amplitude=self.amplitude)
            line = None  # in case inp hits EOF: return empty line
        return line, cnt
blockdict["BOUNDARY"] = BlockBoundary


######################################################
###    blocks: model definition - mpc
######################################################

class BlockMPC(Block):
    """class for reading of a *MPC block in the input file

    updates:
    model.mpc ... list of tuples, each tuple consisting of:
                  1. string defining the type (first item of the data line)
                  2. list of node numbers
    """

    options = OptionList((
        ("user", "NoArgs"),
        ("mode", "OneName")))

    def __init__(self, line, model):
        Block.__init__(self, "MPC", line, model)
        if self.mode:
            model.mpc.mode_option = self.mode.upper()

    def readDataLines(self, inp, model):
        "read *MPC data lines until next block starts"
        cnt = 0
        for cnt, line in self.filteredLines(inp, 1):
            if line[0]=="*":
                break
            data = line.split(",")
            points = list()
            for item in data:
                try:
                    item = int(item)
                except ValueError:
                    item = item.strip()
                points.append(item)

            if points[0] == 0:
                model.mpc[-1][1].extend(points[1:])
            else:
                model.mpc.append((points[0], points[1:]))

            line = None  # in case inp hits EOF: return empty line
        return line, cnt
blockdict["MPC"] = BlockMPC


######################################################
###    Blocks: model definition - surface
######################################################

class BlockSurface(Block):
    """class for reading of a *SURFACE block in the input file

    updates:
    model.surface which is a dict: surface name -> list of datalines of
    the *SURFACE command

    If type="ELEMENT" each dataline with two items is stored in a dict:
    {Face identifier (second value on data line) : elset name / element number
    (first value on data line)}
    If type != "ELEMENT" each dataline is stored as list of string values.
    """

    options = OptionList((
        ("name", "OneName"),
        ("type", "OneName")))

    def __init__(self, line, model):
        Block.__init__(self, "SURFACE", line, model)

        # default type
        if not(self.type):
            self.type = "ELEMENT"

    def readDataLines(self, inp, model):
        "read *SURFACE data lines until next block starts"
        cnt = 0
        for cnt, line in self.filteredLines(inp, 1):
            if line[0]=="*":
                break
            data = [x.strip() for x in line.split(",")]
            if self.type == "ELEMENT":
                if self.name not in model.surface:
                    model.surface[self.name] = dict()
                if len(data)==2:
                    model.surface[self.name][data[1]] = data[0]
                else:
                    model.surface["all exterior faces"] = data
            else:
                try:
                    model.surface[self.name].append(data)
                except KeyError:
                    model.surface[self.name] = [data,]

            line = None  # in case inp hits EOF: return empty line
        return line, cnt
blockdict["SURFACE"] = BlockSurface


######################################################
###    Blocks: model definition - materials
######################################################

class BlockMaterial(SuperBlock):
    """
    class for reading a *MATERIAL block in the input file

    creates:
    model.material: a L{container.MaterialDict} { material name : material }
    material is a L{Material} object, basically another dict
    { parameter name: value/list of values}
    """

    options = OptionList((
        ("name", "OneName"),))

    def __init__(self, line, model):
        Block.__init__(self, "MATERIAL", line, model)

        # check: option NAME mandatory
        if not self.name:
            raise ReadError(
                "*MATERIAL command without NAME-option.")

        if self.name in model.material:
            msg("WARNING: *MATERIAL, NAME=%s is being redefined on"
                " line %d. The first definition will be ignored."
                % (self.name, model.linecnt))
        self.material = Material(self.name)
        model.material[self.name] = self.material

    def readDataLines(self, inp, model):
        """no data line for *MATERIAL, ignore any data lines until next block
        starts"""
        cnt = 0
        err = False
        for cnt, line in self.filteredLines(inp, 1):
            if line[0]=="*":
                break
            err = True
            line = None  # in case inp hits EOF: return empty line
        if err:
            msg("WARNING: *MATERIAL, NAME=%s on line %d appears to have %d"
                " data line(s). It will be ignored."
                % (self.name, model.linecnt, cnt))
        return line, cnt
blockdict["MATERIAL"] = BlockMaterial


class SkipBlockMaterial(SkipSuperBlock):
    validSubBlocks = BlockMaterial.validSubBlocks
BlockMaterial.SkipThisBlock = SkipBlockMaterial


class BlockElastic(Block):
    """class for reading an *ELASTIC block in the input file

    creates:
    model.material[...]['EMod']
    model.material[...]['Nue']
    """

    options = OptionList((
        ("type", "OneName", -1),))

    def __init__(self, line, model, superblock):
        Block.__init__(self, "ELASTIC", line, model)
        self.superblock = superblock

        # check: TYPE=ISOTROPIC
        if self.type and (
                not("ISOTROPIC".startswith(self.type))
                or len(self.type)<2):
            raise ReadError(
                "ELASTIC-TYPE must be ISOTROPIC. Found %s"
                % self.type)

    def readDataLines(self, inp, model):
        "read *ELASTIC data lines until next block starts"
        cnt = 0
        err = False
        for cnt, line in self.filteredLines(inp, 1):
            if line[0]=="*":
                break

            # evaluate data line
            if cnt==1:
                data = line.split(",")
                self.superblock.material['EMod'] = float(data[0])
                self.superblock.material['Nue'] = float(data[1])
            else:
                err = True
            line = None  # in case inp hits EOF: return empty line

        if err:
            msg("WARNING: Temperature or field dependent data for *ELASTIC"
                " will be ignored for material %s on line %d."
                % (self.superblock.name, model.linecnt))
        return line, cnt
BlockMaterial.validSubBlocks["ELASTIC"] = BlockElastic


class BlockPlastic(Block):
    """class for reading an *PLASTIC block in the input file

    creates:
    model.material[...]['Plastic']
    """

    options = OptionList((
        ("hardening", "OneName", -1),))

    def __init__(self, line, model, superblock):
        Block.__init__(self, "PLASTIC", line, model)
        self.superblock = superblock

        # check: hardening=ISOTROPIC
        if self.hardening and (
                not("ISOTROPIC".startswith(self.hardening))
                or len(self.type)<2):
            raise ReadError(
                "PLASTIC-hardening must be ISOTROPIC. Found %s"
                % self.hardening)

        # initialize list to collect tuples of
        # (stress, strain, temp, field1, field2, ...)
        # store in material['Plastic']
        self.para = list()
        self.superblock.material['Plastic'] = self.para

    def readDataLines(self, inp, model):
        """read *PLASTIC data lines until next block starts

        Convert data lines to list of floats and store in self.para.
        """
        cnt = 0
        for cnt, line in self.filteredLines(inp, 1):
            if line[0]=="*":
                break
            # convert to list of floats and store in self.para
            self.para.append(map(float, line.split(",")))
            line = None  # in case inp hits EOF: return empty line
        return line, cnt
BlockMaterial.validSubBlocks["PLASTIC"] = BlockPlastic


class BlockDensity(Block):
    """class for reading a *DENSITY block in the input file

    creates:
    model.material[...]['Density']
    """

    def __init__(self, line, model, superblock):
        Block.__init__(self, "DENSITY", line, model)
        self.superblock = superblock

    def readDataLines(self, inp, model):
        "read data lines until next block starts"
        cnt = 0
        err = False
        for cnt, line in self.filteredLines(inp, 1):
            if line[0]=="*":
                break
            if cnt==1:
                data = line.split(",")
                self.superblock.material['Density'] = float(data[0])
            else:
                err = True
            line = None  # in case inp hits EOF: return empty line
        if err:
            msg("WARNING: Temperatur dependent data for *DENSITY will be"
                "ignored for material %s on line %d."
                % (self.superblock.name, model.linecnt))
        return line, cnt
BlockMaterial.validSubBlocks["DENSITY"] = BlockDensity


class BlockDepvar(Block):
    """class for reading a *DEPVAR block in the input file

    *DEPVAR,DELETE=<UmatDelete>
    <UmatDepvarNum>

    creates:
    model.material[...]['UmatDepvarNum']
    model.material[...]['UmatDelete']
    """

    options = OptionList((
        ("delete", "OneInt"),))

    def __init__(self, line, model, superblock):
        Block.__init__(self, "DEPVAR", line, model)
        self.superblock = superblock
        if hasattr(self, "delete") and self.delete is not None:
            self.superblock.material['UmatDelete'] = self.delete

    def readDataLines(self, inp, model):
        "read data lines until next block starts"
        cnt = 0
        err = False
        for cnt, line in self.filteredLines(inp, 1):
            if line[0]=="*":
                break
            if cnt==1:
                data = line.split(",")
                self.superblock.material['UmatDepvarNum'] = int(data[0])
            else:
                err = True
            line = None  # in case inp hits EOF: return empty line
        if err:
            msg("WARNING: More than one data line for *DEPVAR card for"
                " material %s on line %d. Ignoring all but the first."
                % (self.superblock.name, model.linecnt))
        return line, cnt
BlockMaterial.validSubBlocks["DEPVAR"] = BlockDepvar


class BlockDamping(Block):
    """class for reading a *DAMPING block in the input file

    *DAMPING, ALPHA=0.5

    creates:
    model.material[...]['Alpha']
    """

    options = OptionList((
        ("alpha", "OneFloat"),))

    def __init__(self, line, model, superblock):
        Block.__init__(self, "DAMPING", line, model)
        self.superblock = superblock
        if hasattr(self, "alpha") and self.alpha is not None:
            self.superblock.material['DampingAlpha'] = self.alpha
        else:
            msg("WARNING: No Parameter ALPHA for *DAMPING, material %s,"
                " line %d. *DAMPING command will be ignored."
                % (self.superblock.name, model.linecnt))

    def readDataLines(self, inp, model):
        """no data line for *DAMPING, ignore any data lines until next block
        starts"""
        cnt = 0
        err = False
        for cnt, line in self.filteredLines(inp, 1):
            if line[0]=="*":
                break
            err = True
            line = None  # in case inp hits EOF: return empty line
        if err:
            msg("WARNING: Datalines for *DAMPING, material %s, line %d,"
                " will be ignored."
                % (self.superblock.name, model.linecnt, cnt))
        return line, cnt

BlockMaterial.validSubBlocks["DAMPING"] = BlockDamping


class BlockUserDefinedField(Block):
    """class for reading a *USER DEFINED FIELD block in the input file

    creates:
    model.material[...]['UserDefinedField'] = True
    """

    def __init__(self, line, model, superblock):
        Block.__init__(self, "USERDEFINEDFIELD", line, model)
        self.superblock = superblock
        self.superblock.material['UserDefinedField'] = True

    def readDataLines(self, inp, model):
        """no data line for *USER DEFINED FIELD, ignore any data lines until
        next block starts"""
        cnt = 0
        err = False
        for cnt, line in self.filteredLines(inp, 1):
            if line[0]=="*":
                break
            err = True
            line = None  # in case inp hits EOF: return empty line
        if err:
            msg("WARNING: Datalines for *USER DEFINED FIELD, material %s, line %d,"
                " will be ignored."
                % (self.superblock.name, model.linecnt, cnt))
        return line, cnt

BlockMaterial.validSubBlocks["USERDEFINEDFIELD"] = BlockUserDefinedField


class BlockUserMaterial(Block):
    """class for reading an *User Material block in the input file

    *USER MATERIAL, CONSTANTS=29

    creates:
    model.material[...]['Umat'] = list of umat parameter
    """

    options = OptionList((
        ("constants", "OneInt"),
        ("type", "OneName")))

    def __init__(self, line, model, superblock):
        Block.__init__(self, "USERMATERIAL", line, model)
        self.superblock = superblock
        if not(hasattr(self, "constants") and self.constants is not None):
            msg("WARNING: No Parameter CONSTANTS for *USER MATERIAL,"
                " material %s, line %d."
                % (self.superblock.name, model.linecnt))
        self.para = list()
        self.superblock.material['Umat'] = self.para
        if hasattr(self, "type") and self.type is not None:
            self.superblock.material['UmatType'] = self.type

    def readDataLines(self, inp, model):
        "read *USER MATERIAL data lines until next block starts"
        cnt = 0
        for cnt, line in self.filteredLines(inp, 1):
            if line[0]=="*":
                break
            try:
                data = [float(d) for d in line.split(",")]
            except ValueError:
                raise ReadError(
                    "At least one parameter of User Material %s, line %d is"
                    " not a number."
                    % (self.superblock.name, model.linecnt))
            self.para.extend(data)
            line = None  # in case inp hits EOF: return empty line

        # check number of data items
        if len(self.para) != self.constants:
            msg("WARNING: Number of Umat parameters (=%d) differs from the"
                " value for the CONSTANTS option (=%d), material %s, line"
                " %d" % (len(self.para), self.constants,
                         self.superblock.name, model.linecnt))
        return line, cnt
BlockMaterial.validSubBlocks["USERMATERIAL"] = BlockUserMaterial


######################################################
###    Blocks: model definition - connector behavior
######################################################

class BlockConnectorBehavior(SuperBlock):
    """
    class for reading a *CONNECTOR BEHAVIOR block in the input file

    creates:
    model.connectorBehavior: a dict { name : parameter value dict } with
    parameter value dict being another dict { parameter name: value/list
    of values}
    """

    options = OptionList((
        ("name", "OneName"),))

    def __init__(self, line, model):
        Block.__init__(self, "CONNECTORBEHAVIOR", line, model)

        # check: option NAME mandatory
        if not self.name:
            raise ReadError("*CONNECTOR BEHAVIOR command without NAME-option.")

        # make sure connectorBehavior-attribute exists
        if self.name in model.connectorBehavior:
            msg("WARNING: *CONNECTOR BEHAVIOR, NAME=%s is being redefined on"
                " line %d. The first definition will be ignored."
                % (self.name, model.linecnt))
        self.behavior = ConnectorBehavior(self.name)
        model.connectorBehavior[self.name] = self.behavior

    def readDataLines(self, inp, model):
        """no data line for *CONNECTOR BEHAVIOR, ignore any data lines
        until next block starts
        """
        cnt = 0
        err = False
        for cnt, line in self.filteredLines(inp, 1):
            if line[0]=="*":
                break
            err = True
            line = None  # in case inp hits EOF: return empty line
        if err:
            msg("WARNING: *CONNECTOR BEHAVIOR, NAME=%s appears to have %d data"
                " line(s) on line %d. They will be ignored."
                % (self.name, cnt, model.linecnt))
        return line, cnt
blockdict["CONNECTORBEHAVIOR"] = BlockConnectorBehavior


class SkipBlockConnectorBehavior(SkipSuperBlock):
    validSubBlocks = BlockConnectorBehavior.validSubBlocks
BlockConnectorBehavior.SkipThisBlock = SkipBlockConnectorBehavior


class BlockConnectorElasticity(Block):
    """class for reading an *CONNECTOR ELASTICITY block in the input file

    creates:
    model.connectorBehavior[...]['elasticity']
    """

    options = OptionList((
        ("rigid", "NoArgs"),
        ("component", "OneInt"),))

    def __init__(self, line, model, superblock):
        Block.__init__(self, "CONNECTORELASTICITY", line, model)
        self.superblock = superblock

        try:
            elastList = self.superblock.behavior['elasticity']
        except KeyError:
            elastList = list()
            self.superblock.behavior['elasticity'] = elastList

        # check: option RIGID
        if self.rigid:
            elastList.append((self.component, "rigid"))
        else:
            raise ReadError(
                "*CONNECTOR ELASTICITY command without RIGID-option not"
                " implemented yet.")

    def readDataLines(self, inp, model):
        """no data line for *CONNECTOR ELASTICITY, ignore any data lines
        until next block starts
        """
        cnt = 0
        err = False
        for cnt, line in self.filteredLines(inp, 1):
            if line[0]=="*":
                break
            err = True
            line = None  # in case inp hits EOF: return empty line
        if err:
            msg("WARNING: *CONNECTOR ELASTICITY for behavior %s appears to"
                " have %d data line(s) on line %d. This is not yet"
                " implemented. The data will be ignored."
                % (self.superblock.name, cnt, model.linecnt))
        return line, cnt
BlockConnectorBehavior.validSubBlocks["CONNECTORELASTICITY"] = \
    BlockConnectorElasticity


class BlockConnectorPlasticity(Block):
    """class for reading an *CONNECTOR PLASTICITY block in the input file

    creates:
    model.connectorBehavior[...]['plasticity']
    """

    options = OptionList((
        ("component", "OneInt"),))

    def __init__(self, line, model, superblock):
        Block.__init__(self, "CONNECTORPLASTICITY", line, model)
        self.superblock = superblock

        try:
            plasticityList = self.superblock.behavior['plasticity']
        except KeyError:
            plasticityList = list()
            self.superblock.behavior['plasticity'] = plasticityList

        plastDispFrictForce = list()
        plasticityList.append((self.component, plastDispFrictForce))
        self.superblock.plastDispFrictForce = plastDispFrictForce

    def readDataLines(self, inp, model):
        """no data line for *CONNECTOR PLASTICITY, ignore any data lines
        until next block starts
        """
        cnt = 0
        err = False
        for cnt, line in self.filteredLines(inp, 1):
            if line[0]=="*":
                break
            err = True
            line = None  # in case inp hits EOF: return empty line
        if err:
            msg("WARNING: *CONNECTOR PLASTICITY for behavior %s appears to"
                " have %d data line(s) on line %d. The data will be ignored."
                % (self.superblock.name, cnt, model.linecnt))
        return line, cnt
BlockConnectorBehavior.validSubBlocks["CONNECTORPLASTICITY"] = \
    BlockConnectorPlasticity


class BlockConnectorHardening(Block):
    """class for reading an *CONNECTOR HARDENING block in the input file

    creates:

    items in the list "hardening table" in the following structure:
    model.connectorBehavior[...]['plasticity'] =
    [(component, hardening table), (component, hardening table), ...]
    """

    def __init__(self, line, model, superblock):
        Block.__init__(self, "CONNECTORHARDENING", line, model)
        self.superblock = superblock

        if not hasattr(self.superblock, "plastDispFrictForce"):
            msg("WARNING: *CONNECTOR HARDENING on line %d is lacking a"
                " preceeding *CONNECTOR PLASTICITY option. The data"
                " will be ignored."
                % model.linecnt)
            # just a dummy list that will be forgotten afterwards
            self.superblock.plastDispFrictForce = list()

    def readDataLines(self, inp, model):
        "read *CONNECTOR HARDENING data lines until next block starts"
        cnt = 0
        for cnt, line in self.filteredLines(inp, 1):
            if line[0]=="*":
                break
            data = [float(x) for x in line.split(",")]
            self.superblock.plastDispFrictForce.append(data)
            line = None  # in case inp hits EOF: return empty line
        return line, cnt
BlockConnectorBehavior.validSubBlocks["CONNECTORHARDENING"] = \
    BlockConnectorHardening


class BlockConnectorDamping(Block):
    """class for reading an *CONNECTOR DAMPING block in the input file

    creates:
    model.connectorBehavior[...]['damping']
    """

    options = OptionList((
        ("type", "OneName"),
        ("component", "OneInt"),))

    def __init__(self, line, model, superblock):
        Block.__init__(self, "CONNECTORDAMPING", line, model)
        self.superblock = superblock

        try:
            dampingList = self.superblock.behavior['damping']
        except KeyError:
            dampingList = list()
            self.superblock.behavior['damping'] = dampingList

        if self.type and self.type.upper()!="VISCOUS":
            raise ReadError(
                "*CONNECTOR DAMPING command with type=%s not"
                " implemented yet. TYPE must be VISCOUS."
                % self.type)

        self.dampingDataLines = list()
        dampingList.append((self.component, self.dampingDataLines))

    def readDataLines(self, inp, model):
        "read *CONNECTOR DAMPING data lines until next block starts"
        cnt = 0
        for cnt, line in self.filteredLines(inp, 1):
            if line[0]=="*":
                break
            data = [float(x) for x in line.split(",")]
            self.dampingDataLines.append(data)
            line = None  # in case inp hits EOF: return empty line
        return line, cnt
BlockConnectorBehavior.validSubBlocks["CONNECTORDAMPING"] = \
    BlockConnectorDamping


class BlockConnectorStop(Block):
    """class for reading an *CONNECTOR STOP block in the input file

    Stores a new (component, stop data lines) tuple in the list
    model.connectorBehavior[...]['stop']
    """

    options = OptionList((
        ("component", "OneInt"),))

    def __init__(self, line, model, superblock):
        Block.__init__(self, "CONNECTORSTOP", line, model)
        self.superblock = superblock

        try:
            stopList = self.superblock.behavior['stop']
        except KeyError:
            stopList = list()
            self.superblock.behavior['stop'] = stopList

        self.stopDataLines = list()
        stopList.append((self.component, self.stopDataLines))

    def readDataLines(self, inp, model):
        "read *CONNECTOR STOP data lines until next block starts"
        cnt = 0
        for cnt, line in self.filteredLines(inp, 1):
            if line[0]=="*":
                break
            data = [float(x) for x in line.split(",")]
            self.stopDataLines.append(data)
            line = None  # in case inp hits EOF: return empty line
        return line, cnt
BlockConnectorBehavior.validSubBlocks["CONNECTORSTOP"] = BlockConnectorStop


class BlockConnectorLock(Block):
    """class for reading an *CONNECTOR LOCK block in the input file

    Stores a new (component, lock argument, lock data lines) triple in the
    list model.connectorBehavior[...]['lock']
    """

    options = OptionList((
        ("component", "OneInt"),
        ("lock", "OneName"),))

    def __init__(self, line, model, superblock):
        Block.__init__(self, "CONNECTORLOCK", line, model)
        self.superblock = superblock

        try:
            lockList = self.superblock.behavior['lock']
        except KeyError:
            lockList = list()
            self.superblock.behavior['lock'] = lockList

        self.lockDataLines = list()
        lockList.append((self.component, self.lock, self.lockDataLines))

    def readDataLines(self, inp, model):
        "read *CONNECTOR LOCK data lines until next block starts"
        cnt = 0
        for cnt, line in self.filteredLines(inp, 1):
            if line[0]=="*":
                break
            data = [float(x) for x in line.split(",")]
            self.lockDataLines.append(data)
            line = None  # in case inp hits EOF: return empty line
        return line, cnt
BlockConnectorBehavior.validSubBlocks["CONNECTORLOCK"] = BlockConnectorLock


######################################################
###    Blocks: model definition - section definition
######################################################

class BlockAnySection(Block):
    """base class for reading of a *xxxx SECTION - block in the input file

    updates:
    model.properties, a PropertiesDict object: {elset: SectionProperty()}
    """
    def __init__(self, name, line, model):
        Block.__init__(self, name, line, model)

        # check: elset must be given
        if not self.elset:
            raise ReadError("Property definition without ELSET option.")

        # create properties-dict:
        #  - type=SOLIDSECTION, SHELLSECTION, ... (self.name)
        #  - append all other options
        prop = dict(type=self.name)
        for op in self.options.getNameList():
            val = getattr(self, op)
            if val is not None:
                prop[op] = val

        # store data in model
        self.prop = model.properties.updateSection(prop)

    def readDataLines(self, inp, model):
        "read *SECTION data lines until next block starts"
        cnt = 0
        line = None  # initialize for if filteredLines yields nothing
        for cnt, line in self.filteredLines(inp, 1):
            if line[0]=="*":
                break
            data = line.split(",")
            try:
                self.prop.addToSection(
                    [float(val) for val in data if val.strip()])
            except ValueError:
                raise ReadError(
                    "At least one parameter on the dataline of *%s %s"
                    " on line %d is not a number:\n%s"
                    % (self.name, self.elset, model.linecnt, line))
            line = None  # in case inp hits EOF: return empty line
        return line, cnt


class BlockSolidSection(BlockAnySection):
    """class for reading of a *SOLID SECTION - block in the input file
    """
    options = OptionList([("elset", "OneName", 1),
                          ("material", "OneName", 1),
                          ("controls", "OneName", 3),
                          ])

    def __init__(self, line, model):
        BlockAnySection.__init__(self, "SOLIDSECTION", line, model)
blockdict["SOLIDSECTION"] = BlockSolidSection


class BlockCohesiveSection(BlockAnySection):
    """class for reading of a *COHESIVE SECTION - block in the input file
    """
    options = OptionList([("elset", "OneName"),
                          ("material", "OneName", 1),
                          ("controls", "OneName", 1),
                          ("response", "OneName", 1),
                          ("thickness", "OneName", 1),
                          ])

    def __init__(self, line, model):
        BlockAnySection.__init__(self, "COHESIVESECTION", line, model)

blockdict["COHESIVESECTION"] = BlockCohesiveSection


class BlockShellSection(BlockAnySection):
    """class for reading of a *SHELL SECTION - block in the input file
    """
    options = OptionList([("elset", "OneName", 1),
                          ("material", "OneName", 1),
                          ("controls", "OneName", 3),
                          ("offset", "OneFloat"),
                          ])

    def __init__(self, line, model):
        BlockAnySection.__init__(self, "SHELLSECTION", line, model)
blockdict["SHELLSECTION"] = BlockShellSection


class BlockBeamSection(BlockAnySection):
    """class for reading of a *BEAM SECTION - block in the input file
    """
    options = OptionList([("elset", "OneName", 1),
                          ("material", "OneName", 1),
                          ("section", "OneName", 1),
                          ])

    def __init__(self, line, model):
        BlockAnySection.__init__(self, "BEAMSECTION", line, model)
blockdict["BEAMSECTION"] = BlockBeamSection


class BlockConnectorSection(BlockAnySection):
    """class for reading of a *CONNECTOR SECTION - block in the input file
    """
    options = OptionList([("elset", "OneName", 1),
                          ("behavior", "OneName", 1),
                          ])

    def __init__(self, line, model):
        BlockAnySection.__init__(self, "CONNECTORSECTION", line, model)

    def readDataLines(self, inp, model):
        "read *CONNECTOR SECTION data lines until next block starts"
        cnt = 0
        for cnt, line in self.filteredLines(inp, 1):
            if line[0]=="*":
                break

            # convert data line (i.e. comma separated string) to list of
            # strings with leading and trailing whitespaces stripped off.
            # Empty strings are replaced by Nones.
            # ... This is equivalent to:
            # data = map(str.strip, line.split(","))
            # for i, val in enumerate(data):
            #     if val=="":
            #         data[i]=None
            data = [
                {"": None}.get(d, d)
                for d in (dd.strip() for dd in line.split(","))]

            self.prop.addToSection(data)
            line = None  # in case inp hits EOF: return empty line
        return line, cnt
blockdict["CONNECTORSECTION"] = BlockConnectorSection


######################################################
###    Blocks: model definition - orientation
######################################################

class BlockOrientation(Block):
    """Class for reading of an *ORIENTATION block in the input file

    updates:
    model.orientation
    """

    options = OptionList((
        ("name", "OneName"),))

    def __init__(self, line, model):
        Block.__init__(self, "ORIENTATION", line, model)

        # check name
        if not self.name:
            raise ReadError("*ORIENTATION command without NAME-option.")

    def readDataLines(self, inp, model):
        "read *ORIENTATION data lines until next block starts"
        cnt = 0
        for cnt, line in self.filteredLines(inp, 1):
            if line[0]=="*":
                break

            data = [float(d) for d in line.split(",")]
            if cnt==1:
                numData = len(data)
                if numData < 6:
                    raise ReadError("*ORIENTATION data line with %d values. At"
                                    " least 6 are required." % numData)
                try:
                    thisOri = model.orientation[self.name]
                except KeyError:
                    thisOri = Orientation(name=self.name)
                    model.orientation[self.name] = thisOri
                thisOri.a = data[0:3]
                thisOri.b = data[3:6]

            line = None  # in case inp hits EOF: return empty line

        if cnt>1:
            msg("WARNING: *ORIENTATION, NAME=%s on line %d has more than one"
                " data lines. All but the first will be ignored."
                % (self.name, model.linecnt))
        return line, cnt
blockdict["ORIENTATION"] = BlockOrientation


######################################################
###    Blocks: model definition - amplitude
######################################################

class BlockAmplitude(Block):
    """Class for reading of an *AMPLITUDE block in the input file

    updates:
    model.amplitude
    """

    options = OptionList((
        ("name", "OneName"),
        ("time", "OneName"),
        ("smooth", "OneFloat"),
        ("definition", "OneName", 3, ),  # might be "smooth step", "tabular"
    ))

    def __init__(self, line, model):
        Block.__init__(self, "AMPLITUDE", line, model)

        # check name
        if not self.name:
            raise ReadError("*AMPLITUDE command without NAME-option.")

        # upper definition
        if self.definition is not None:
            self.definition = self.definition.upper().replace(" ", "")

        # check: smooth option only with definition=TABULAR
        if (self.smooth and
                self.definition is not None and self.definition!="TABULAR"):
            warnings.warn(
                "Amplitude definition %s has a smooth option and"
                " definition=%s. This is illegal. Ignoring smooth option."
                % (self.name, self.definition))

        # check: definition if present is TABULAR or SMOOTH STEP
        if (self.definition is not None
                and self.definition!="TABULAR"
                and self.definition!="SMOOTHSTEP"):
            warnings.warn(
                "Amplitude definition %s with definition=%s. This is not"
                " supported by the reader yet. Data might be meaningless."
                % (self.name, self.definition))

    def readDataLines(self, inp, model):
        "read data lines until next block starts"
        cnt = 0
        for cnt, line in self.filteredLines(inp, 1):
            if line[0]=="*":
                break
            data = [float(d) for d in line.split(",")]
            numData = len(data)
            if numData % 2 != 0:
                raise ReadError(
                    "*AMPLITUDE, NAME=%s on line %d has a data line with an"
                    " odd number (%d) of values."
                    % (self.name, model.linecnt, numData))
            if cnt==1:
                if self.name in model.amplitude:
                    msg("WARNING: *AMPLITUDE, NAME=%s on line %d already"
                        " exists in the model. Ignoring all but the last"
                        " definition."
                        % (self.name, model.linecnt))
                # create new amplitude object
                thisAmpl = Amplitude(
                    name=self.name, smooth=self.smooth,
                    time=self.time, definition=self.definition)
                model.amplitude[self.name] = thisAmpl

            # append data to amplitude object
            for i in range(0,numData,2):
                thisAmpl.append(data[i:i+2])

            line = None  # in case inp hits EOF: return empty line
        return line, cnt
blockdict["AMPLITUDE"] = BlockAmplitude


######################################################
###    Blocks: model definition - time points
######################################################

class BlockTimePoints(Block):
    """Class for reading of an *TIME POINTS block in the input file

    updates:
    model.timepoints
    """

    options = OptionList((
        ("name", "OneName"),))

    def __init__(self, line, model):
        Block.__init__(self, "TIMEPOINTS", line, model)

        # check name
        if not self.name:
            raise ReadError("*TIME POINTS command without NAME-option.")

    def readDataLines(self, inp, model):
        "read *TIME POINTS data lines until next block starts"
        cnt = 0
        for cnt, line in self.filteredLines(inp, 1):
            if line[0]=="*":
                break
            data = [float(d) for d in line.split(",")]

            if cnt==1:
                if self.name in model.timepoints:
                    msg("WARNING: *TIME POINTS, NAME=%s on line %d already"
                        " exists in the model. Ignoring all but the last"
                        " definition."
                        % (self.name, model.linecnt))
                # store time points in model
                timepoints = data
                model.timepoints[self.name] = timepoints
            else:
                timepoints.extend(data)

            line = None  # in case inp hits EOF: return empty line
        return line, cnt
blockdict["TIMEPOINTS"] = BlockTimePoints


######################################################
###    Blocks: model definition - initial condition
######################################################

class BlockInitCond(Block):
    """Class for reading of an *INITIAL CONDITIONS block in the input file

    updates:
    model.initCond
    """

    options = OptionList((
        ("type", "OneName"),
        ("geostatic", "NoArgs"),
        ("variable", "OneInt"),))

    def __init__(self, line, model):
        Block.__init__(self, "INITIALCONDITIONS", line, model)

        # check name
        if not self.type:
            raise ReadError(
                "*INITIAL CONDITIONS command without TYPE-option.")

        if self.type in initCondTypeMultLine:

            # current list of InitialCond objects, ignore all options
            self.dataList = model.initCond.addNewList(self.type)

            self.skip = False
        elif self.type in initCondTypeToClass:
            options = self.getOptionsDict()
            del options["type"]
            # current list of InitialCond objects
            self.dataList = model.initCond.addNewList(self.type, **options)

            self.skip = False
        else:
            msg("WARNING: *INITIAL CONDITIONS, TYPE=%s occured on line %d. This"
                " type has not been implemented so far."
                % (self.type, model.linecnt))
            self.skip = True

    def readDataLines(self, inp, model):
        "read *INITIAL CONDITIONS data lines until next block starts"

        if self.skip:
            cnt = -1
            line = None
            for cnt, line in enumerate(inp):
                if line[0]=="*" and line[1]!="*":
                    break
                line = None  # in case inp hits EOF: return empty line
            return line, cnt+1

        # initial condition for one elset or nset on multiple lines
        # e.g. *INITIAL CONDITIONS, TYPE=SOLUTION
        if self.type in initCondTypeMultLine:
            cnt = 0
            for cnt, line in self.filteredLines(inp, 1):
                if line[0]=="*":
                    break

                newData = line.split(",")
                try:
                    newData = [float(d) for d in newData]
                except ValueError:
                    # There must be a non float on the line, therefore we assume
                    # that this is a starting line with an elset or nset label
                    # as first item

                    # initialize new InitialCond object
                    # NOTE: this is completely different from list.append() !!!
                    # ... This append method returns a list (actually of type
                    #     InitialCondTypeSolution) as current list of values
                    currInitCond = self.dataList.append(
                        label=newData[0].strip())

                    # remove first item from line
                    try:
                        newData = [float(d) for d in newData[1:]]
                    except ValueError:
                        raise ReadError(
                            "Unexpected value could not be converted to float"
                            " on line %d:\n%s" % (model.linecnt, line))

                # in any case append values to current list
                currInitCond.extend(newData)

                line = None  # in case inp hits EOF: return empty line

        # only one initial condition per line
        # first item is element or node /-set label/number
        # subsequent items are floats
        else:
            cnt = 0
            for cnt, line in self.filteredLines(inp, 1):
                if line[0]=="*":
                    break

                data = line.split(",")
                d = data[0]
                try:
                    d = int(d)
                except ValueError:
                    d = d.strip().upper()
                self.dataList.append(label=d, listInit=map(float, data[1:]))

                line = None  # in case inp hits EOF: return empty line

        return line, cnt
blockdict["INITIALCONDITIONS"] = BlockInitCond


######################################################
###    Blocks: Loads, Steps, Reports
######################################################

# yet to be implemented....
