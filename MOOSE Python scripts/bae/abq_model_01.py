"""abq_model_01.py: DEPRECATED, use L{bae.abq_model_02}

This is for reading and writing abaqus input file data. Everything not
belonging to data I/O should go in other modules and classes.


Usage:
======


Initialization:
---------------

from bae.abq_model_01 import Model
m = Model()    # an optional parameter may state the already opened log file


Reading model data:
-------------------

m.read(inp_name)             # specify filename *with* extension .inp
print m.get_data_members()   # shows all newly loaded data attributes of m

-- or --
# read only *NODE and *ELEMENT - blocks
# no warnings about unknown commands in this case
m.read(inp_name, recognized_blocks=['NODE', 'ELEMENT'])

-- or --
# read all known block types except *BOUNDARY blocks
blocks = m.recognized_blocks
blocks.remove("BOUNDARY")
m.read(inp_name, blocks)

The recognized_blocks argument may consist of the following abaqus keywords
*without* the leading star.

ABAQUS keyword | created members of the Model object
---------------+------------------------------------
*NODE          | node_coord, possibly nset
*ELEMENT       | el_node, el_type, type_el, possibly elset
*NSET          | nset
*ELSET         | elset
*BOUNDARY      | node_bc
*SURFACE       | surface


Note:
- only those data members that are actually needed are created, for example
  there is no m.nset member at all (not even an empty dict or a None), if there
  are no nsets defined in the input file or you didn't chose them to be loaded.
- If there is at least one *ELEMENT block, m.el_node, m.el_type and m.type_el
  are created all three.
- Multiple *NSET commands (or *NODE, NSET=... commands) for the same *NSET
  will subsequently add nodes to the same nset. Same with elsets.
- No check for duplicity is performed. The same node may be in one nset several
  times. If the same node with its coordinates is defined multiple times, all
  but the last definition will be silently ignored. Same with elements.
- You can perform multiple read commands, new data will be added to the
  model subsequently.
- All arguments to specific command options are converted to uppercase, e.g.
  NSET-name "gaga" becomes "GAGA" in *NSET, NSET=gaga


Data structure / how to add and remove data manually:
-----------------------------------------------------

This data is created by the m.read() method. You can as well create the model
data (or part of it) from scratch and then export it (see section "Writing
model data" below).
See also the get_data_members() function in section "Other service functions"
below.

m = Model()

# Nodes
m.node_coord = NodesDict()  # basically a dict with a write method
m.node_coord[1031] = [12.7, 51.1, -4.2]  # {node number : coord-list}
...

# Elements
m.el_node = ElementNodesDict()
m.el_type = dict()
m.type_el = dict()              # this member is necessary for the write method

m.el_node[4623] = [10, 510, 11, 542] # {element number: list of node numbers}
m.el_type[4623] = "C3D4"             # {element number: element type string}
m.el_type["C3D4"].append(4623)       # {element type: list of element numbers}
...

# or simpler...  (replaces all six commands above!)
m.update_elem(4623, "C3D4", [10, 510, 11, 542])

# remove an element (also from all elsets but don't remove nodes)
m.remove_elem(4623)


# Nsets
m.nset = NsetsDict()
m.nset["ALL"].extend([10, 510, 11, 542]) # {nset name: list of node numbers}
# or ...
m.create_nset("ALL").extend([10, 510, 11, 542])
...

# Elsets
m.elset = ElsetsDict()
m.elset["ALL"].append(4623)           # {elset name: list of element numbers}
# or ...
m.create_elset("ALL").append(4623)
...

# Boundary conditions
m.node_bc[4623] = "PINNED"    # either {node number: boundary type string}
m.node_bc[4624] = [1,3,0.0]   # or {node num: [first dof, last dof, magnitude]}
...

# Boundary conditions per nset
# (this is not created automatically when reading from input file)
# useful only for writing
m.nset_bc = NsetBoundaryDict()
m.nset_bc["BOTTOM"] = "PINNED"     # either {node number: boundary type string}
m.nset_bc["RIGHT"] = [1, 1, 5.0]  # or {node: [first dof, last dof, magnitude]}
...


Writing model data to file:
---------------------------

m.write("mymodel.inp")
-- or --
m.write("mymodel.inp", ["NODE", "NSET"])
-- or --
fo = open("mynodelist.inp", "w")
fo.write("** Some extra things to write to this file.\n")
m.write(fo, ["NODE", "NSET"])
fo.close()

The second argument (named recognized_blocks) may consist of the following
abaqus keywords *without* the leading star. A single string or an iterable of
strings is accepted.

ABAQUS keyword | members of the Model object providing the data
---------------+-----------------------------------------------
*NODE          | node_coord
*ELEMENT       | el_node, type_el
*NSET          | nset
*ELSET         | elset
*BOUNDARY      | node_bc, nset_bc
*SURFACE       | surface


The following writing methods may also be used, for more flexibiliy.
All writing methods listed below take an open file object as first argument,
so it goes like this:

fo = open("mynodelist.inp", "w")
m.nset.write_one(fo, "MYNODESET")
fo.close()

... and that's what is in so far:

m.write_all_elems(fo)
m.node_coord.write(fo)    # write all node_coords

m.nset.write_one(fo, ...)      # write just one nset
m.nset.write_all(fo)
m.elset.write_one(fo, ...)     # write just one elset
m.elset.write_all(fo)
m.node_bc.write(fo)       # write all bc's

# m.nset_bc not created by m.read() !
m.nset_bc = NsetBoundaryDict()
m.nset_bc["BOTTOM"] = "PINNED"
m.nset_bc.write(fo)


Other service functions:
------------------------

m.get_data_members()
        Returns a set of strings naming the data attributes of the Model
        object. Useful to find out easily what's been loaded by the read()
        method or created manually so far.

m.get_nodes_from_elset(elset_name=None)
        Returns a set of nodes that are connected to the elset given as
        argument. If elset_name==None, return all nodes connected to any
        element in the model.


Adding new structures to read from the input file:
--------------------------------------------------

That's quite easy. For now see the other Block...classes and ask Gero.
... and please update the documentation.
"""
###############################################################################
###############################################################################
###############################################################################

_version_ = "1.8"

_version_history_ = """\
Versions:
=========

read_inp...
2.1 : use of module string removed
2.2 : switched to new style classes (without effect) and cut the overhead of
      block reading classes, complete translation
2.3 : added rudimentary *SURFACE Block reader, log-output may be suppressed
      by Model()-Argument logfile=None

abq_model_01
1.1 : GP: changed name to abq_model
      changed interface: Model() without filename, read() with filename-arg
1.2 : GP: read only selected blocks; added internal, instance-options,
      warn if they occur; get_data_objects() method
1.3 : GP: * some write methods: model.write(), model.write_all_elems(),
      model.node_coord.write(), model.nset.write_one() / ...write_all(),
      model.elset.write_one() / ...write_all() / model.node_bc.write(),
      * lots of more documentation
      * reading *BOUNDARY "direct" format (to specify magnitudes, still no step
        info though)
      * new methods to add data manually: create_elset, create_nset,
        create_surface
1.4 : GP surfaces -> surface (like all the other: elset, nset, el_node ...)
      corrected some bugs in write methods, improved error messages
      update_elem not tested yet
1.5 : GP fixed bugs in Model.update_elem(), ElementNodesDict.write_elem_lines()
1.6 : GP fixed a bug in Model().read()
      model.write() with only one ELSET
1.7 : GP write nsets and elsets ordered, accept string as recognized_blocks
      parameter in model.read()
      added: model.remove_elem(element_id)
      added: model.get_submodel_in_sphere()
      changed: now possible: m = Model().read('abaqusinputfile.inp')
      fixed: runs with python 2.4 (all builtin-function)
      no comment lines after each nset, elset
1.8 : GP added reading of some material options, section definition
      groupIter format argument may be function
      more tolerant to incorrect data
"""

_bugs_ = """\
known bugs/problems:
====================

- does not process ABAQUS parts/instances input files
  use command mdb.models[modelName].setValues(noPartsInputFile=ON)
  or in cae menu: Model -> Edit Attributes -> Model-1
  mark "do not use parts and assemblies in input file"
- only some blocks / ABAQUS-input-file-commands implemented so far;
  anything besides the pure mesh info is not implemented adequately.
- whitespace-handling not optimal, int(), float() conversion accept
  leading and trailing whitespaces! strip is necessary only for string
  arguments
- If the last command is *USER MATERIAL the CONSTANTS parameter is not
  being checked against the actual number of parameters

ideas:
======

- use csv-module for reading large datablocks? does it read numbers?
- in the Model-object at initialization, create one static Block-object
  for each Block-type. Blocks to be ignored are of type SkipBlock (or
  similar but specialized). Then those instances don't have to be created
  and destroyed while reading. Attributes like the options-warnflag may
  be changed after creating the Model-object before actually reading.
- speed optimizations: replace map(float..), map(int...) and re.split!
  for further ideas do the following:
  import cProfile
  from bae.abq_model_01 import Model
  m = Model()
  cProfile.run('m.read("somebigmodel.inp")')
"""


from sys import stdout, stderr
from os import access, R_OK
from re import compile, IGNORECASE

import types

from bae.vecmath_01 import *
from bae.future_01 import *

######################################################
######################################################
###   Model - Service classes for reading the input file
###   those classes don't actually hold the bulk data
###   they are just for interpreting the input file
######################################################
######################################################

class ReadError(Exception): pass
class UnknownOption(ReadError): pass

class CmdOption(object):
    """Descriptor class that finds one command option in the ABAQUS input file,
    something like NSET=myset or GENERATE, and creates an appropriate
    attribute/member in the owner class instance.

    Uses a regular expression to find a certain kind of command option
    in the string argument passed to the member function find_in.

    The first argument to the constructor "keyword" states the name of the
    option in the ABAQUS input file (e.g. GENEREATE) as well as the name of the
    attribute in the owner class instance. While the option name is case
    insensitive, the attribute name of the owner class is not, therefore it's
    suggested to use lowercase letters.

    There are different types of command options, stated by the second
    argument to __init__:
    typ == "NoArgs": option without argument, e.g. GENERATE
          its value is True (found) or False
    typ == "OneName": option with one argument, e.g. ELSET=myset
          its value is None (keyword not found) or the options argument string
          (e.g. "myset" in case of ELSET=myset)
    typ == "OneInt": option with one argument, e.g. DEPVAR=20
          its value is None (keyword not found) or integer specified
    typ == "OneFloat": option with one argument, e.g. ALPHA=0.5
          its value is None (keyword not found) or float specified

    If the member warnflag is set to True, each occurence of this option in the
    input file will trigger a warning in the log.
    """
    def __init__(self, keyword, typ="NoArgs"):
        self.typ = typ
        if self.typ == "NoArgs":
            self.regex = compile(keyword, IGNORECASE)
            self.default = False
        elif self.typ == "OneName":
            self.regex = compile(keyword + r"\s*=\s*(.*)", IGNORECASE)
            self.default = None
        elif self.typ == "OneInt":
            self.regex = compile(keyword + r"\s*=\s*(\d*)", IGNORECASE)
            self.default = None
        elif self.typ == "OneFloat":
            self.regex = compile(keyword + r"\s*=\s*(.*)", IGNORECASE)
            self.default = None
        else:
            raise AssertionError, \
                "ERROR: Unknown CmdOption type %s." % typ
        self.attr_name = keyword
        self.warnflag = False

    def init_member_attr(self, instance):
        setattr(instance, self.attr_name, self.default)
        
    def find_in(self, instance, arg, model):
        "set option if arg contains it. string values are converted uppercase"
        set_fnd = self.regex.match(arg)
        if set_fnd:
            if self.typ == "NoArgs":
                setattr(instance, self.attr_name, True)
            elif self.typ == "OneName":
                setattr(instance, self.attr_name, set_fnd.group(1).upper())
            elif self.typ == "OneInt":
                res = set_fnd.group(1)
                try:
                    resi = int(res)
                except ValueError:
                    raise ReadError(
                        "The value %s for option %s for command *%s on line %d"
                        " can not be converted to an integer."
                        % (res, self.attr_name.upper(), instance.name, model.linecnt))
                setattr(instance, self.attr_name, resi)
            elif self.typ == "OneFloat":
                res = set_fnd.group(1)
                try:
                    resf = float(res)
                except ValueError:
                    raise ReadError(
                        "The value %s for option %s for command *%s on line %d"
                        " can not be converted to a number."
                        % (res, self.attr_name.upper(), instance.name, model.linecnt))
                setattr(instance, self.attr_name, resf)

            # issue a warning, if warnflag is True
            if self.warnflag:
                model.logfile.write(
                    "WARNING: Option %s for command *%s occured on line %d.\n"
                    % (self.attr_name.upper(), instance.name, model.linecnt))
        return bool(set_fnd)

class OptionList(list):
    """List of CmdOption - objects

    just pass a list of (keyword, type) - tuples

    A third optional element of this tuple can be the string "warning". If so
    the attribute warnflag of the correspondig CmdOption-object will be set to
    True and a warning will be issued, when this option occurs in the input
    file.
    """
    def __init__(self, key_type_list):
        for opt_dat in key_type_list:
            key, typ = opt_dat[0:2]
            self.append(CmdOption(key, typ))
            if len(opt_dat)>2:
                self[-1].warnflag = (opt_dat[2][:4].lower()=='warn')
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

    Any child class that is to be used for reading must implement a readline
    method that recieves the current line from the input file for evaluation.
    It must return True if the block has not found its end. Or False if the end
    has been reached; the same line will then possibly again be evaluated for
    the next block.

    If the child class is a subblock of a superblock, its readline method must
    accept an additional argument superblock stating the superblock-object.

    Any block-class that needs special treatment for skipping this block (e.g.
    a *STEP-block does not end until its *END STEP counterpart), has to
    redefine the "SkipThisBlock" class member as a reference to a previously
    defined specialized class (possibly derived from SkipBlock).
    (example: class SkipMyBlock: ....
     then: class MyBlock(Block): SkipThisBlock = SkipMyBlock;def __init__ ....)
    The Block.SkipThisBlock class attribute is assigned after the definition of
    the class "SkipBlock".
    """

    re_split = compile(r"\s*,\s*")
    re_key = compile(r"\*(.*?)\s*(,.*)*$")
    options = list()

    def __init__(self, name, line, model, options = True):
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
        """extract the abaqus command key word from line

        line has to be stripped of spaces from the right end already:
           ... line=line.rstrip() already performed
        returns None or False if line is empty
        """
        if line == "":
            return None
        keymatch = self.re_key.match(line.replace(" ", ""))
        if not(keymatch):
            key = "???"
            raise ReadError(
                "Invalid beginning of block:\n%s" % line)
        key = keymatch.group(1).upper()
        return key

    def read_options(self, line, model):
        """extracts command options from the first line of a ABAQUS command
        block and assigns to appropriate member attributes"""
        for arg in self.re_split.split(line.replace(" ", ""))[1:]:
            if not arg: continue # ignore empty options (double comma or so)
            found = False
            for opt in self.options:
                found = opt.find_in(self, arg, model)
                if found: break
            # not a known option
            if not found:
                raise UnknownOption("Unknown option of the *%s-command:\n%s"
                                    % (self.name, arg))

    def datlist(self,line):
        return self.re_split.split(line.strip())

    def datlistfloat(self,line):
        data = self.datlist(line)
        try:
            data = map(float,data)
        except ValueError:
            raise ReadError("could not convert all items to floats for the"
                            " data line\n%s" % line)
        return data

class SkipBlock(Block):
    """class for skipping a block that is to be skipped, because
    it is not in the list recognized_blocks or because it is a block of
    unknown type.

    If a block contains subblocks (e.g. *STEP lasts until *END STEP)
    this blockclass has to provide its own specialized SkipBlock class
    with a special readline function. The skipblock-method of this block
    must then return a new object of this specialized SkipBlock class.
    """
    def __init__(self, line, model):
        Block.__init__(self, line.split(",")[0].strip().upper(),
                       line, model, options = False)

    def readline(self, line, model):
        "skip lines until next block starts"
        if line == "":
            return True
        return line[0] != "*"

# initialize default SkipThisBlock - class member as class SkipBlock
Block.SkipThisBlock = SkipBlock


class SuperBlock(Block):
    """base class for super blocks like *MATERIAL that have dependent subblocks
    following like *ELASTIC

    Any child class should call Block.__init__

    see description of class Block
    """
    validSubBlocks = dict()
    pass


class SkipSuperBlock(Block):
    """class for skipping a super block and all its sub blocks
    """
    def __init__(self, line, model):
        Block.__init__(self, line.split(",")[0].strip().upper(),
                       line, model, options = False)

    def readline(self, line, model):
        "skip lines until next block starts"
        if line == "":
            return True
        if line[0] != "*":
            return True
        else:
            key = self.read_key(line)
            return key in self.validSubBlocks

# initialize default SkipThisBlock - class member as class SkipBlock
SuperBlock.SkipThisBlock = SkipSuperBlock
# add a line like this for each SuperBlock-child class (here BlockXXXX)
# class SkipBlockXXXX(SkipSuperBlock):
#     validSubBlocks = BlockXXXX.validSubBlocks


######################################################
###    Blocks: general
######################################################


class BlockHeading(Block):
    "Class for reading the *HEADING - block"

    def __init__(self, line, model):
        Block.__init__(self, "HEADING", line, model, options = False)
    def readline(self, line, model):
        "read until next block starts"
        if line == "":
            return True
        return line[0] != "*"
blockdict["HEADING"] = BlockHeading


######################################################
###    Blocks: model definition - nodes
######################################################

class BlockNode(Block):
    """class for reading of a *NODE - block in the input file

    creates:
    model.node_coord (NodesDict nodenumber -> list of coords)
    possibly model.nset (NsetsDict nset-name -> list of nodenumbers)
    """

    options = OptionList([("nset", "OneName"),])

    def __init__(self, line, model):
        Block.__init__(self, "NODE", line, model)
        
        # make sure node_coord-attribute exists
        if not(hasattr(model, "node_coord")):
            model.node_coord = NodesDict()

        # find or create nset
        if self.nset:
            model.create_nset(self.nset, raise_exception = False)

    def readline(self, line, model):
        "read nodes line, return if line has been read"
        if line == "":
            return True
        if line[0] == "*":
            return False
        data = self.datlist(line)
        nodenum = int(data[0])
        model.node_coord[nodenum] = map(float,data[1:7])
        if self.nset:
            model.nset[self.nset].append(nodenum)
        if len(data) > 7:
            errstr = "more than 6 node coordinates. No 7 ff read:\n" \
                + ", ".join(data[7:])
            raise ReadError, errstr
        return True
blockdict["NODE"] = BlockNode    


class BlockNset(Block):
    """class for reading a *NSET - block in the input file

    creates:
    model.nset (NsetsDict nset-name -> list of nodenumbers)
    """

    options = OptionList((
            ("nset", "OneName"),
            ("generate", "NoArgs"),
            ("instance", "OneName", "warning"),
            ("internal", "NoArgs", "warning"),
            ("elset", "OneName", "warning"),
            ))

    def __init__(self, line, model):
        Block.__init__(self, "NSET", line, model)

        # check nset-name
        if not self.nset:
            raise ReadError, \
                "*NSET command without NSET-option."
        model.create_nset(self.nset, raise_exception = False)
            
    def readline(self, line, model):
        "read *NSET data line, return if line has been read"
        if line == "":
            return True
        if line[0] == "*":
            return False
        data = self.datlist(line)
        if self.generate:
            if len(data) < 2:
                raise ReadError, \
                      "*NSET,GENERATE with less than 2 data values."
            start = int(data[0])
            stop = int(data[1])
            if len(data) > 2:
                inc = int(data[2])
            else:
                inc = 1
            model.nset[self.nset].extend(range(start, stop+1, inc))
        else:
            for d in data:
                if d == "":
                    continue
                elif d.isdigit():
                    model.nset[self.nset].append(int(d))
                else:
                    d = d.upper()
                    if d in model.nset:
                        model.nset[self.nset].extend(model.nset[d])
                    else:
                        model.logfile.write(
                            ("WARNING: NSET %s specified as part of NSET %s, "
                             +"but it is not in the model.\n") % (d, self.nset))
        return True
blockdict["NSET"] = BlockNset


######################################################
###    blocks: model definition - elements
######################################################

class BlockElement(Block):
    """class for reading of an *ELEMENT block in the input file

    creates:
    model.el_node (ElementNodesDict: element number -> list of node numbers)
    model.el_type (dict: element number -> element type)
    model.type_el (dict: element type -> list of element numbers)
    possibly model.elset (ElsetsDict: elset name -> list of element numbers)

    If you want to change this BlockElement-implementation: It is probably
    tempting to use model.update_elem()-method. However this would probably be
    slower than this.
    """

    options = OptionList((
            ("type", "OneName"),
            ("elset", "OneName")))

    def __init__(self, line, model):
        Block.__init__(self, "ELEMENT", line, model)

        # make sure el_node-, el_type- and type_el-attributes exist
        if not(hasattr(model, "el_node")):
            model.el_node = ElementNodesDict()
        if not(hasattr(model, "el_type")):
            model.el_type = dict()
        if not(hasattr(model, "type_el")):
            model.type_el = dict()

        # find or create elset
        if self.elset:
            model.create_elset(self.elset, raise_exception = False)

        self.data = list() # temporary list of element numbers
            
    def readline(self, line, model):
        "read *ELEMENT data line, return if actually read (True) or not"
        if line == "":
            return True
        if line[0] == "*":
            return False
        if line[-1] == ",":  # nodenumbers to be continued
            self.data.extend(self.datlist(line)[:-1])
        else: # one element finished
            self.data.extend(self.datlist(line))
            elnum = int(self.data[0])

            model.el_node[elnum] = map(int, self.data[1:])
            model.el_type[elnum] = self.type

            if not model.type_el.has_key(self.type):
                model.type_el[self.type] = list()
            model.type_el[self.type].append(elnum)

            if self.elset:
                model.elset[self.elset].append(elnum)
            self.data = list()
        return True
blockdict["ELEMENT"] = BlockElement    


class BlockElset(Block):
    """class for reading an *ELSET - block

    creates:
    model.elset (ElsetsDict: elset name -> list of element numbers)
    """

    options = OptionList((
            ("elset", "OneName"),
            ("generate", "NoArgs"),
            ("instance", "OneName", "warning"),
            ("intarnal", "NoArgs", "warning")))

    def __init__(self, line, model):
        Block.__init__(self, "ELSET", line, model)

        # check: option ELSET mandatory
        if not self.elset:
            raise ReadError, \
                "*ELSET command without ELSET-option."

        model.create_elset(self.elset, raise_exception = False)
            
    def readline(self, line, model):
        "read *ELSET data line, return if actually read (True) or not"
        if line == "":
            return True
        if line[0] == "*":
            return False
        data = self.datlist(line)
        if self.generate:
            if len(data) < 2:
                raise ReadError, \
                      "*ELSET,GENERATE has less than 2 values in data line."
            try:
                start = int(data[0])
                stop = int(data[1])
            except ValueError:
                raise ReadError(
                    "*ELSET,GENERATE has non-integers as start ('%s') or stop"
                    " ('%s') values in the data line." % (data[0], data[1]))
            if len(data) > 2:
                try:
                    inc = int(data[2])
                except ValueError:
                    raise ReadError(
                        "*ELSET,GENERATE has the non-integer '%s' as increment"
                        " value in the data line." % data[2])
            else:
                inc = 1
            model.elset[self.elset].extend(range(start, stop+1, inc))
        else:
            for d in data:
                if d == "":
                    continue
                elif d.isdigit():
                    model.elset[self.elset].append(int(d))
                else:
                    d = d.upper()
                    if d in model.elset:
                        model.elset[self.elset].extend(model.elset[d])
                    else:
                        model.logfile.write(
                            ("WARNING: ELSET %s specified as part of ELSET %s, "
                             +"but it is not in the model.\n") % (d, self.elset))
        return True
blockdict["ELSET"] = BlockElset    


######################################################
###    Blocks: model definition - boundary conditions
######################################################

class BlockBoundary(Block):
    """Class for reading of a *BOUNDARY block in the input file

    This is not very useful yet, since it's not recognized whether the BC's
    are set within a certain step. The TYPE-option is also neglected.

    creates:
    model.boundary
    model.node_bc
    """

    options = OptionList((
            ("type", "OneName"),
            ("amplitude", "OneName"),
            ))

    def __init__(self, line, model):
        Block.__init__(self, "BOUNDARY", line, model)

        if not(hasattr(model, "node_bc")):
            model.node_bc = NodeBoundaryDict(model)
        if not(hasattr(model, "boundary")):
            model.boundary = BoundaryDict()
        if not(hasattr(self, "amplitude")):
            self.amplitude = ""
        if not(hasattr(self, "smooth")):
            self.smooth = None
            
    def readline(self, line, model):
        "read *BOUNDARY data line, return if actually read (True) or not"
        if line == "":
            return True
        if line[0] == "*":
            return False
        data = self.datlist(line)

        # which nodes / nsets
        if data[0].isdigit():
            key = int(data[0])
        else:
            key = data[0]

        # which flavour of BC definition?
        try:
            dof_1 = int(data[1])
        except ValueError:
            # zero-valued bc's using the "type" format (model data only)
            bc = data[1].upper()
        else:
            # using the "direct" format
            bc = [dof_1, dof_1, 0.0]
            if len(data)>2:
                bc[1] = int(data[2])
            if len(data)>3:
                bc[2] = float(data[3])

        # store data in dictionary
        model.boundary[key].append(BoundaryCondition(
                bc, amplitude=self.amplitude, smooth=self.smooth))
        return True
blockdict["BOUNDARY"] = BlockBoundary    


######################################################
###    blocks: model definition - surface
######################################################

class BlockSurface(Block):
    """class for reading of a *SURFACE block in the input file

    creates:
    model.surface (dict: surface name -> list of datalines of the *SURFACE command)

    this is not a very useful implementation yet. It is only useful for reproducing the
    *SURFACE command, not for accessing the information in detail.
    """

    options = OptionList((
            ("name", "OneName"),
            ("type", "OneName")))

    def __init__(self, line, model):
        Block.__init__(self, "SURFACE", line, model)

        # default type
        if not(self.type):
            self.type = "ELEMENT"

        # make sure surface[name] exists
        model.create_surface(self.name, raise_exception=False)

    def readline(self, line, model):
        "read *SURFACE data line, return if actually read (True) or not"
        if line == "":
            return True
        if line[0] == "*":
            return False
        model.surface[self.name].append(line)
        return 1
blockdict["SURFACE"] = BlockSurface


######################################################
###    blocks: model definition - materials
######################################################

class BlockMaterial(SuperBlock):
    """class for reading a *MATERIAL block in the input file

    creates:
    model.material (dict: material name ->
                       (dict: parameter name -> value/list of values))
    """

    options = OptionList((
            ("name", "OneName"),))

    def __init__(self, line, model):
        Block.__init__(self, "MATERIAL", line, model)

        # check: option NAME mandatory
        if not self.name:
            raise ReadError, \
                "*MATERIAL command without NAME-option."

        # make sure material-attribute exists
        if not(hasattr(model, "material")):
            model.material = dict()
        if self.name in model.material:
            model.logfile.write(
                "WARNING: *MATERIAL, NAME=%s is being redefined on"
                " line %d. The first definition will be ignored.\n"
                % (self.name, model.linecnt))
        self.material = dict()
        model.material[self.name] = self.material

    def readline(self, line, model):
        "no data line for *MATERIAL, ignore any data lines"
        if line == "":
            return True
        if line[0] == "*":
            return False
        else:
            model.logfile.write(
                "WARNING: *MATERIAL, NAME=%s appears to have a data line on"
                " line %d. It will be ignored.\n"
                % (self.name, model.linecnt))
            return True

blockdict["MATERIAL"] = BlockMaterial
class SkipBlockMaterial(SkipSuperBlock):
    validSubBlocks = BlockMaterial.validSubBlocks


class BlockElastic(Block):
    """class for reading an *ELASTIC block in the input file

    creates:
    model.material[...]['EMod']
    model.material[...]['Nue']
    """

    def __init__(self, line, model, superblock):
        Block.__init__(self, "ELASTIC", line, model)
        self.superblock = superblock

    def readline(self, line, model):
        if line == "":
            return True
        if line[0] == "*":
            return False
        # evaluate data line
        data = self.datlist(line)
        if ('EMod' in self.superblock.material
            or 'Nue' in self.superblock.material):
            model.logfile.write(
                "WARNING: Temperatur dependent data for *ELASTIC will be"
                "ignored for material %s on line %d.\n"
                % (self.superblock.name, model.linecnt))
        self.superblock.material['EMod'] = float(data[0])
        self.superblock.material['Nue'] = float(data[1])
        return True

BlockMaterial.validSubBlocks["ELASTIC"] = BlockElastic


class BlockDensity(Block):
    """class for reading an *DENSITY block in the input file

    creates:
    model.material[...]['Density']
    """

    def __init__(self, line, model, superblock):
        Block.__init__(self, "DENSITY", line, model)
        self.superblock = superblock

    def readline(self, line, model):
        if line == "":
            return True
        if line[0] == "*":
            return False
        # evaluate data line
        data = self.datlist(line)
        if 'Density' in self.superblock.material:
            model.logfile.write(
                "WARNING: Temperatur dependent data for *DENSITY will be"
                "ignored for material %s on line %d.\n"
                % (self.superblock.name, model.linecnt))
        self.superblock.material['Density'] = float(data[0])
        return True

BlockMaterial.validSubBlocks["DENSITY"] = BlockDensity


class BlockDepvar(Block):
    """class for reading an *DEPVAR block in the input file

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
        if hasattr(self, "delete") and self.delete!=None:
            self.superblock.material['UmatDelete'] = self.delete

    def readline(self, line, model):
        if line == "":
            return True
        if line[0] == "*":
            return False
        # evaluate data line
        data = self.datlist(line)
        self.superblock.material['UmatDepvarNum'] = float(data[0])
        return True

BlockMaterial.validSubBlocks["DEPVAR"] = BlockDepvar


class BlockDamping(Block):
    """class for reading an *DAMPING block in the input file

    *DAMPING, ALPHA=0.5

    creates:
    model.material[...]['Alpha']
    """

    options = OptionList((
            ("alpha", "OneFloat"),))

    def __init__(self, line, model, superblock):
        Block.__init__(self, "DAMPING", line, model)
        self.superblock = superblock
        if hasattr(self, "alpha") and self.alpha!=None:
            self.superblock.material['DampingAlpha'] = self.alpha
        else:
            model.logfile.write(
                "WARNING: No Parameter ALPHA for *DAMPING, material %s,"
                " line %d. *DAMPING command will be ignored.\n"
                % (self.superblock.name, model.linecnt))

    def readline(self, line, model):
        if line == "":
            return True
        if line[0] == "*":
            return False
        else:
            model.logfile.write(
                "WARNING: Datalines for *DAMPING, material %s, line %d,"
                " will be ignored.\n" % (self.superblock.name, model.linecnt))
            return True

BlockMaterial.validSubBlocks["DAMPING"] = BlockDamping


class BlockUserMaterial(Block):
    """class for reading an *User Material block in the input file

    *USER MATERIAL, CONSTANTS=29

    creates:
    model.material[...]['Umat'] = list of umat parameter
    """

    options = OptionList((
            ("constants", "OneInt"),))

    def __init__(self, line, model, superblock):
        Block.__init__(self, "USERMATERIAL", line, model)
        self.superblock = superblock
        if not(hasattr(self, "constants") and self.constants!=None):
            model.logfile.write(
                "WARNING: No Parameter CONSTANTS for *USER MATERIAL,"
                " material %s, line %d.\n"
                % (self.superblock.name, model.linecnt))
        self.para = list()
        self.superblock.material['Umat'] = self.para

    def readline(self, line, model):
        "read *USER MATERIAL data line, return if actually read (True) or not"
        if line == "":
            return True
        if line[0] == "*":
            if len(self.para) != self.constants:
                model.logfile.write(
                    "WARNING: Number of Umat parameters (=%d) differs from the"
                    " value for the CONSTANTS option (=%d), material %s, line"
                    " %d\n" % (len(self.para), self.constants,
                               self.superblock.name, model.linecnt))
            return False

        try:
            data = [float(d) for d in self.datlist(line)]
        except ValueError:
            raise ReadError(
                "At least one parameter of User Material %s is not a number."
                % (self.superblock.name))
        self.para.extend(data)
        return True

BlockMaterial.validSubBlocks["USERMATERIAL"] = BlockUserMaterial


######################################################
###    blocks: model definition - section definition
######################################################


class BlockAnySection(Block):
    """base class for reading of a *xxxx SECTION - block in the input file

    creates:
    model.properties (PropertiesDict object)
    """

    options = OptionList([("elset", "OneName"),
                          ("material", "OneName"),
                          ("controls", "OneName"),])

    def __init__(self, name, line, model):
        Block.__init__(self, name, line, model)
        
        # make sure node_coord-attribute exists
        if not(hasattr(model, "properties")):
            model.properties = PropertiesDict()

        self.prop = model.properties.addSection(self)

    def readline(self, line, model):
        "read line, return if line has been read"
        if line == "":
            return True
        if line[0] == "*":
            return False
        data = self.datlist(line)
        try:
            self.prop.addToSection(
                [float(val) for val in data
                 if val!=""])
        except ValueError:
            raise ReadError(
                "At least one parameter on the dataline of *%s %s"
                " is not a number." % (self.name, self.elset))
        return True

class BlockSolidSection(BlockAnySection):
    """class for reading of a *SOLID SECTION - block in the input file
    """
    def __init__(self, line, model):
        BlockAnySection.__init__(self, "SOLIDSECTION", line, model)
blockdict["SOLIDSECTION"] = BlockSolidSection

class BlockCohesiveSection(BlockAnySection):
    """class for reading of a *COHESIVE SECTION - block in the input file
    """
    options = OptionList([("elset", "OneName"),
                          ("material", "OneName"),
                          ("controls", "OneName"),
                          ("response", "OneName"),
                          ("thickness", "OneName"),
                          ])
    def __init__(self, line, model):
        BlockAnySection.__init__(self, "COHESIVESECTION", line, model)

blockdict["COHESIVESECTION"] = BlockCohesiveSection

class BlockShellSection(BlockAnySection):
    """class for reading of a *SHELL SECTION - block in the input file
    """
    def __init__(self, line, model):
        BlockAnySection.__init__(self, "SHELLSECTION", line, model)
blockdict["SHELLSECTION"] = BlockShellSection


class BlockBeamSection(BlockAnySection):
    """class for reading of a *BEAM SECTION - block in the input file
    """
    options = OptionList([("elset", "OneName"),
                          ("material", "OneName"),
                          ("controls", "OneName"),
                          ("section", "OneName"),
                          ])
    def __init__(self, line, model):
        BlockAnySection.__init__(self, "BEAMSECTION", line, model)
blockdict["BEAMSECTION"] = BlockBeamSection


######################################################
###    Blocks: model definition - amplitude
######################################################

class BlockAmplitude(Block):
    """Class for reading of an *AMPLITUDE block in the input file

    creates:
    model.amplitude
    """

    options = OptionList((
            ("name", "OneName"),
            ("time", "OneName"),
            ("smooth", "OneFloat"),
            ))

    def __init__(self, line, model):
        Block.__init__(self, "AMPLITUDE", line, model)

        if not(hasattr(model, "amplitude")):
            model.amplitude = defaultdict(list)

        # check name
        if not self.name:
            raise ReadError, \
                "*AMPLITUDE command without NAME-option."
            
    def readline(self, line, model):
        "read *AMPLITUDE data line, return if actually read (True) or not"
        if line == "":
            return True
        if line[0] == "*":
            return False
        data = self.datlistfloat(line)
        numData = len(data)
        if numData%2 != 0:
            raise ReadError("*AMPLITUDE data line with odd number (%d) of"
                            " values." % numData)
        for i in range(0,numData,2):
            model.amplitude[self.name].append(data[i:i+2])
        return True
blockdict["AMPLITUDE"] = BlockAmplitude


######################################################
###    Blocks: model definition - time points
######################################################

class BlockTimePoints(Block):
    """Class for reading of an *TIME POINTS block in the input file

    creates:
    model.timepoints
    """

    options = OptionList((
            ("name", "OneName"),
            ))

    def __init__(self, line, model):
        Block.__init__(self, "TIMEPOINTS", line, model)

        if not(hasattr(model, "timepoints")):
            model.timepoints = defaultdict(list)

        # check name
        if not self.name:
            raise ReadError, \
                "*TIME POINTS command without NAME-option."
            
    def readline(self, line, model):
        "read *TIME POINTS data line, return if actually read (True) or not"
        if line == "":
            return True
        if line[0] == "*":
            return False
        data = self.datlistfloat(line)
        model.timepoints[self.name].extend(data)
        return True
blockdict["TIMEPOINTS"] = BlockTimePoints


######################################################
###    Blocks: Loads, Steps, Reports
######################################################

# yet to be implemented....


######################################################
######################################################
###    Classes, functions for actually holding the
###    data and for writing them to a file
######################################################
######################################################

class ItemAlreadyThere(ValueError): pass


def groupIter(input, maxcounts=8, format=None):
    """yields lists of max maxcounts elements from the "input"-iterable object
    format may be a formating string in that case a list of strings is returned
    format may also be a function accepting a single item of input an returning
       the formated string for this item

    It's a generator function. You may use in a for clause:
    for group in groupIter(mylist, 3):
    """
    input = iter(input)
    if (format!=None and type(format) not in (str, types.FunctionType)):
        raise TypeError(
            "groupIter does not recognize the format argument %s of type %s."
            " Only strings and functions are valid."
            % (str(format), str(type(val))))
    while 1:
        output = list()
        for i in range(maxcounts):
            try:
                val = input.next()
            except StopIteration:
                break
            if type(format)==str:
                try:
                    val = (format % val)
                except TypeError:
                    raise TypeError(
                        "groupIter could not convert %s of type %s with"
                        " the format string '%s'"
                        % (str(val), str(type(val)), format))
            elif type(format) == types.FunctionType:
                try:
                    val = format(val)
                except TypeError:
                    raise TypeError(
                        "groupIter could not convert %s of type %s with"
                        " the given format function."
                        % (str(val), str(type(val))))
            # output val
            output.append(val)

        if output: # no empty lists at the end
            yield output
        else:
            break


class NodesDict(dict):
    """class for holding the node coordinates

    this is basically a dict, implements write method
    """
    __slots__ = ['fmt_nodenum','fmt_coord','write']

    fmt_nodenum = "%7d"
    fmt_coord = "%12f"

    @classmethod
    def fmt_coord_func(self, coord):
        return self.fmt_coord % coord

    def write(self, output_file, nset_name=None):
        """write the *node command with data in this object
        output_file is a file object of the already opened abaqus input file

        if nset_name is given, an ,NSET-option is given, that means all nodes
        belong to that nset, e.g. "ALL"
        """
        if nset_name:
            output_file.write("*NODE, NSET=%s\n" % nset_name)
        else:
            output_file.write("*NODE\n")
        ids = self.keys()
        ids.sort()
        for id in ids:
            output_file.write(
                (self.fmt_nodenum % id) + ", "
                + ", ".join(map(self.fmt_coord_func, self[id])) + "\n")
        return


class NsetsDict(dict):
    """class for holding the nsets

    implements write methods
    """
    __slots__ = ['fmt_nodenum','write', 'write_all']

    fmt_nodenum = "%7d"

    def write_one(self, output_file, nset_name):
        """write the nset with the name nset_name
        output_file is a file object of the already opened abaqus input file

        silently ignore empty nsets
        """
        if nset_name not in self:
            raise KeyError, (
                "No NSET %s, when asked to write it." % nset_name)
        if len(self[nset_name])==0:
            output_file.write("**  NSET %s contains no nodes.\n" % nset_name)
            return
        output_file.write("*NSET, NSET=%s\n" % nset_name)
        nodes = self[nset_name]
        nodes.sort()
        for line in groupIter(nodes,
                              maxcounts=8, format=self.fmt_nodenum):
            output_file.write(", ".join(line)+"\n")
        return

    def write_all(self, output_file):
        """write all nsets using the write() method in alphabetical order
        output_file is a file object of the already opened abaqus input file
        silently ignore empty nsets
        """
        nset_list = self.keys()
        nset_list.sort()
        for nset_name in nset_list:
            self.write_one(output_file, nset_name)


class ElementNodesDict(dict):
    """class for holding the elements nodes
    {element number: list of node numbers}

    implements write method
    """
    __slots__ = ['fmt_elnum','fmt_nodenum','write']

    fmt_elnum = "%7d"
    fmt_nodenum = "%7d"

    def write_elem_lines(self, output_file, elset):
        """write the element definitions of all elements listed in
        elset as datalines in an *ELEMENT-block.

        Important Note: the *ELEMENT-command line is not written!

        output_file is a file object of the already opened abaqus input file
        elset is a list of elementnumbers

        silently ignore invalid element numbers
        """
        for element in elset:

            # create linelist to look like that: ["4132", ", 10", ", 101", ", 11", ...]
            linelist = [(self.fmt_elnum % element)+", ",]
            try:
                nodes = self[element]
            except KeyError:
                continue
            for node in nodes[:-1]:
                linelist.append((self.fmt_nodenum % node)+", ")
            linelist.append(self.fmt_nodenum % nodes[-1])
            for line in groupIter(linelist, maxcounts=8):
                output_file.write("".join(line)+"\n")
                
        return


class ElsetsDict(dict):
    """class for holding the elsets
    {elset name: lsit of element numbers}

    implements write method
    """
    __slots__ = ['fmt_elnum', 'write', 'write_all']

    fmt_elnum = "%7d"

    def write_one(self, output_file, elset_name):
        """write the elset with the name elset_name
        output_file is a file object of the already opened abaqus input file

        silently ignore empty elsets
        """
        if elset_name not in self:
            raise KeyError, (
                "No ELSET %s, when asked to write it." % elset_name)
        if len(self[elset_name])==0:
            output_file.write("**  ELSET %s contains no elements.\n" % elset_name)
            return
        output_file.write("*ELSET, ELSET=%s\n" % elset_name)
        elements = self[elset_name]
        elements.sort()
        for line in groupIter(elements,
                              maxcounts=8, format=self.fmt_elnum):
            output_file.write(", ".join(line)+"\n")
        return

    def write_all(self, output_file):
        """write all elsets using the write() method
        output_file is an file object of the already opened abaqus input file
        silently ignore empty elsets
        """
        elset_list = self.keys()
        elset_list.sort()
        for elset_name in elset_list:
            self.write_one(output_file, elset_name)


class BoundaryCondition(list):
    """class for holding one boundary condition of the type [first dof,
    last dof, magnitude], like at the *BOUNDARY data line.

    additionally it can store the amplitude and smooth option
    """
    __slots__ = ['amplitude', 'smooth']
    def __init__(self, vals, amplitude='', smooth=None):
        list.__init__(self, vals)
        self.amplitude = amplitude
        self.smooth = smooth


class BoundaryDict(defaultdict):
    """class for holding the boundary conditions
    It's basically a dict {node number / nset label : bcList}, where bcList is
    a list of ABAQUS boundary type strings and/or BoundaryCondition objects.
    """
    __slots__ = ['fmt_nodenum', 'fmt_bc', 'write']

    fmt_nodenum = "%7d"
    fmt_bc = "%2d, %2d, %12G"

    def __init__(self):
        defaultdict.__init__(self, list)

    def write(self, output_file):
        """write the *boundary command with data in this object
        output_file is a file object of the already opened abaqus input file

        bugs: doesn't write amplitude and smooth options!
        """
        output_file.write("*BOUNDARY\n")
        ids = self.keys()
        ids.sort()
        for id in ids:
            bcList = self[id]
            for bc in bcList:
                if type(bc)==str:
                    bc = self.fmt_bc % tuple(bc)
                output_file.write((self.fmt_nodenum + ", %s\n") % (id, bc))
        return


class NodeBoundaryDict(object):
    """class for finding boundary conditions for a specific node
    It works like a dict {node number : bcList}, where bcList is a list of
    ABAQUS boundary type strings and/or BoundaryCondition objects.
    """
    def __init__(self, model):
        self.model = model

    def __getitem__(self, key):
        bcList = list()

        # BCs for single node
        try:
            bcList.extend(self.model.boundary[key])
        except KeyError:
            pass
        
        # BCs for nsets
        if hasattr(self.model, 'nset'):
            for nsetname, nsetbcList in self.model.boundary.iteritems():
                if type(nsetname)!=str:
                    continue
                try:
                    nset = self.model.nset[nsetname]
                    if key in nset:
                        bcList.extend(nsetbcList)
                except KeyError:
                    pass

        # if key not found -> raise KeyError
        if bcList:
            return bcList
        else:
            raise KeyError("No bC found for %s." % str(key))


class NsetBoundaryDict(BoundaryDict):
    """class for holding the boundary conditions for nsets
    It's basically a dict {nset name : bc}, where bc is either a ABAQUS
    boundary type string or a list [first dof, last dof, magnitude], like at
    the *BOUNDARY data line.

    At the moment this class has been implemented for writing only. Now it is 
    a reference to the new BoundaryDict class that does the same. This class
    will be removed in the next version of the abq_model module.
    """
    pass



class SurfacesDict(dict):
    """class for holding surfaces
    It's basically a dict {surface name : list of *SURFACE datalines}.

    At the moment this class stores the surface data only in the form of
    complete data line strings as they appear in the ABAQUS input file.
    """
    __slots__ = ['write_one', 'write_all']

    def write_one(self, output_file, surf_name, surf_type):
        """write the *surface command with data in this object
        output_file is a file object of the already opened abaqus input file
        surf_name specifies the name and surf_type is the type to appear on
        the command line in the input file.
        """
        output_file.write("*SURFACE, NAME=%s, TYPE=%s\n" % (surf_name, surf_type))
        for line in self[surf_name]:
            output_file.write(line+"\n")
        return

    def write_all(self, output_file, surf_type = "ELEMENT"):
        """write the *surface commands for all surfaces in this object
        assuming they all have the same type. Order them alphabetically.

        output_file is a file object of the already opened abaqus input file
        surf_type is the type to appear on the command line in the input file.
        """
        ids = self.keys()
        ids.sort()
        for surf_name in ids:
            self.write_one(output_file, surf_name, surf_type)
        return


class SectionProperty(dict):
    """dictionary for all properties of one section definition

    Data lines of the abaqus *... SECTION command are added to the 'data'-item.
    Parameters given as options to the abaqus *... SECTION command are added
    to the dictionary by the self.addSome() function. This creates an item of
    the same name as the abaqus command option (as stated in the attrnames
    argument to the addSome function).
    """
    __slots__ = ['addToSection', 'addSome']

    def addToSection(self, data):
        "add the list of parameters data to the 'data'-item of the prop."
        if 'data' in self:
            self['data'].extend(data)
        else:
            self['data'] = list(data)

    def addSome(self, source, attrnames):
        "add data members of source listed in attrnames, but only if not None"
        for atna in attrnames:
            try:
                val = getattr(source, atna)
            except AttributeError:
                continue
            if val != None:
                self[atna] = val

class PropertiesDict(dict):
    """dictionary with all property definitions in the model

    for section properties (*BEAM SECTION, *SHELL SECTION, *SOLID SECTION) the
    key value is the elset name, the value is an object of type
    SectionProperty.

    self.matElset ... dictionary {material name: set of elset names}
    """
    __slots__ = [
        # member attributes
        'matElset',
        # method attributes
        'addSection',
        ]

    def __init__(self):
        dict.__init__(self)
        self.matElset = defaultdict(set)

    def addSection(self, obj):
        prop = SectionProperty()
        prop['type'] = obj.__class__.__name__
        addPropList = obj.options.getNameList()
        try:
            addPropList.remove('elset')
        except ValueError:
            raise ReadError("Property def. w/o elset option. Program error.")
        prop.addSome(obj, addPropList)

        self[obj.elset] = prop
        if obj.elset and obj.material:
            self.matElset[obj.material].add(obj.elset)
        return prop


######################################################
######################################################
###    Model class, holding it all together
######################################################
######################################################

class QuietLog(object):
    def write(self, s):
        pass

class Model(object):
    """Reads an ABAQUS input file and provides the data as member attributes

    You may supply an open logfile for diagnostic output or None to suppress
    it to the constructor.
    """

    def __init__(self, logfile = stdout):
        "if logfile == None, no diagnostic output"
        if logfile is None:
            self.logfile = QuietLog()
        else:
            self.logfile = logfile
        self.linecnt = None
        self.recognized_blocks = blockdict.keys()

        # attribute of the obj models so far,
        # what comes later is considered a "data attribute", see data_objects()
        self.initial_attribs = None
        self.initial_attribs = set(dir(self))

    def get_data_members(self):
        """Returns a set of strings naming the data attributes of the Model
        object. Useful to find out easily what's been loaded by the read()
        method."""
        return set(dir(self)) - self.initial_attribs

    def get_nodes_from_elset(self, elset_name=None):
        """Returns a set of nodes that are connected to the elset given as
        argument. If elset_name==None, return all nodes connected to any
        element in the model."""

        assert(hasattr(self, "el_node"))

        if elset_name is None:
            elset = self.el_node.iterkeys()
        else:
            assert(hasattr(self, "elset"))
            elset = self.elset[elset_name]

        nodes = set()
        for element in elset:
            nodes.update(self.el_node[element])
        return nodes


    ######################################################
    ## read method

    def read(self, inp_name, recognized_blocks = None):
        """read an input file, adding the data to this Model object
        for further docs see the modules doc string

        return self, so the following is possible:
        m = Model().read('abaqusinputfile.inp')
        """

        if type(recognized_blocks)==str:
            if recognized_blocks.upper() not in self.recognized_blocks:
                raise ReadError(
                    "ERROR: The following input file command can not be"
                    " read yet: %s" % recognized_blocks)
            else:
                recognized_blocks = ( recognized_blocks.upper(), )
        elif recognized_blocks:
            recognized_blocks = set(
                [b.upper() for b in recognized_blocks])
            if not recognized_blocks.issubset(self.recognized_blocks):
                raise ReadError(
                    "ERROR: The following input file commands are not"
                    " implemented yet: %s"
                    % ", ".join(recognized_blocks-set(self.recognized_blocks)))

        self.logfile.write('Start reading input file %s ...\n' % inp_name)
        if not access(inp_name, R_OK):
            raise IOError("File %s can not be opened or read." % inp_name)
        self.inp_name = inp_name

        #--- read the File
        f = file(self.inp_name, "r")
        self.linecnt = 0

        # inblock is the object of the readblock that is currently being read
        # superblock is a block object the current position belongs to.
        # e.g. while reading an *ELASTIC block, superblock should be a
        # *MATERIAL block.
        inblock = None
        superblock = None
        for line in f:
            self.linecnt += 1

            # skip comment line
            if line[:2] == "**":
                continue

            try:
                line = line.rstrip()

                # already inside a block
                if inblock:
                    still_in_block = inblock.readline(line, self)
                    if not(still_in_block):
                        inblock = None
                        # note: superblock remains active

                # not yet or no longer in a block
                if not(inblock):
                    key = Block.read_key(line)
                    if not(key):
                        # empty line, loop: next line
                        continue

                    if superblock:
                        try:
                            inblock = superblock.validSubBlocks[key](
                                line, self, superblock)
                        except KeyError:
                            inblock = None
                            superblock = None
                        if inblock:
                            # loop: next line
                            continue

                    if ((recognized_blocks and (key in recognized_blocks))
                        or ((not recognized_blocks) and (key in blockdict))):
                        # load this block
                        inblock = blockdict[key](line, self)
                        if isinstance(inblock, SuperBlock):
                            superblock = inblock
                    elif (key in blockdict):
                        # skip this known block
                        inblock = blockdict[key].SkipThisBlock(line, self)
                    else:
                        # skip this unknown block
                        inblock = SkipBlock(line, self)

                    if ((not recognized_blocks) and (key not in blockdict)):
                        # issue a warning if block not known
                        self.logfile.write(
                            "WARNING: Unknown command:\n%s\n" % line)

            except ReadError, err:
                self.logfile.write("ERROR in the %s-block, in %s, line %d:\n"
                                   % (key, inp_name, self.linecnt))
                self.logfile.write(err.args[0]+'\n')
                raise
        self.logfile.write('Finished reading input file %s.\n' % inp_name)
        return self


    ######################################################
    ## methods for writing


    def write(self, output_file, recognized_blocks = None):
        """write the model data to an input file
        for further docs see the modules doc string

        output_file may be an open file object (which is left open by this
        method) or a filename (a string) in which case a file will be
        created and closed afterwards.

        recognized_blocks may be a string like "NODE", "ELSET" or a iterable
        of such strings.
        """

        if type(recognized_blocks)==str:
            if recognized_blocks.upper() not in self.recognized_blocks:
                raise ReadError, (
                    "ERROR: The following input file command can not be"
                    " written yet: %s" % recognized_blocks)
            else:
               recognized_blocks = ( recognized_blocks.upper(), )
        elif recognized_blocks:
            recognized_blocks = set(
                [b.upper() for b in recognized_blocks])
            if not recognized_blocks.issubset(self.recognized_blocks):
                raise ReadError, (
                    "ERROR: The following input file commands can not be"
                    " written yet:" +
                    ", ".join(recognized_blocks - set(self.recognized_blocks)))
        else:
            recognized_blocks = self.recognized_blocks

        if type(output_file) == str:
            output_file = open(output_file, "w")
            close_file = True
        else:
            close_file = False

        if (('NODE' in recognized_blocks)
            and hasattr(self, "node_coord")):
            assert type(self.node_coord) == NodesDict
            self.node_coord.write(output_file)

        if ('ELEMENT' in recognized_blocks
            and hasattr(self, "el_node")):
            self.write_all_elems(output_file)

        if (('NSET' in recognized_blocks)
            and hasattr(self, "nset")):
            assert type(self.nset) == NsetsDict
            self.nset.write_all(output_file)

        if (('ELSET' in recognized_blocks)
            and hasattr(self, "elset")):
            assert type(self.elset) == ElsetsDict
            self.elset.write_all(output_file)

        if ('BOUNDARY' in recognized_blocks):
            if hasattr(self, "nset_bc"):
                assert type(self.nset_bc) == NsetBoundaryDict
                self.nset_bc.write(output_file)
            if hasattr(self, "node_bc"):
                assert type(self.node_bc) == NodeBoundaryDict
                self.node_bc.write(output_file)

        if (('SURFACE' in recognized_blocks)
            and hasattr(self, "surface")):
            assert type(self.surface) == SurfacesDict
            self.surface.write_all(output_file)

        if close_file:
            output_file.close()

        return

    
    def write_all_elems(self, output_file, elset_name=None):
        """write all element definitions in the model

        output_file is a file object of the already opened abaqus input file

        if elset_name is given, all elements are assigned to that elset by
        means of an ELSET-option, useful for elset "ALL"
        """
        assert type(self.el_node) == ElementNodesDict

        if elset_name:
            elset_string = ", ELSET=" + elset_name
        else:
            elset_string = ""

        for eltype, elems in self.type_el.items():
            if len(elems)==0:
                continue
            output_file.write("*ELEMENT, TYPE=%s%s\n" % (eltype, elset_string))
            self.el_node.write_elem_lines(output_file, elems)

        return


    ######################################################
    ## methods for adding data manually

    def create_nset(self, nset_name, raise_exception = True):
        """Create a new nset in the model, create nset member if needed
        Return the nset whether just created or not
        """
        return self.create_data_entry(
            "nset", NsetsDict, nset_name, list, raise_exception)

    def create_elset(self, elset_name, raise_exception = True):
        """Create a new elset in the model, create elset member if needed
        Return the elset whether just created or not
        """
        return self.create_data_entry(
            "elset", ElsetsDict, elset_name, list, raise_exception)

    def create_surface(self, surf_name, raise_exception = True):
        """Create a new surface in the model.
        see create_data_entry method for details
        """
        return self.create_data_entry(
            "surface", SurfacesDict, surf_name, list, raise_exception)


    def create_data_entry(self,
                          attr_name, attr_type,
                          inst_name, inst_type,
                          raise_exception):
        """For internal use only. Do not call directly.

        Create a new entry in a certain data dict in the model object
        (i.e. self).
        The data dict is given by its attr_name and attr_type.
        inst_name is the key value in the data dict and inst_type is the type,
        in case it has yet to be created.
        
        See examples in create_nset or create_elset for how to use it.

        Return the member of the data dict (its value), whether it's just been
        created or not.

        If there is already a entry with this name in the dict, leave it as it
        is. Raise an ItemAlreadyThere exception if the flag raise_exception
        is True.
        """
        try:
            this_dict = getattr(self, attr_name)
        except AttributeError:
            this_dict = attr_type()
            setattr(self, attr_name, this_dict)

        if inst_name not in this_dict:
            this_dict[inst_name] = inst_type()
        elif raise_exception:
            raise ItemAlreadyThere(
                "When trying to create %s['%s'] found it in "
                "the model already." % (attr_name, inst_name))
        
        return this_dict[inst_name]


    def update_elem(self, el_num, el_type, nodes):
        """inserts a new element into the model
        updates el_node, el_type, type_el dicts but no nset

        Silently replaces an existing element with the same number.
        """
        # node connectivity
        self.create_data_entry("el_node", ElementNodesDict, el_num, list,
                               raise_exception = False)
        self.el_node[el_num] = nodes

        # type of this element
        if not(hasattr(self, "el_type")):
            self.el_type = dict()
        self.el_type[el_num] = el_type

        # elements of a certain type
        self.create_data_entry("type_el", dict, el_type, list,
                               raise_exception = False).append(el_num)


    def remove_elem(self, el_num, remove_nodes=False):
        """removes an element from the model
        updates el_node, el_type, type_el dicts and all elsets

        remove_nodes=True ... remove unused nodes. This is not very efficient
           if many elements are beeing deleted, in this case rather clean up
           after all elements have been removed using get_nodes_from_elset()

        Silently ignores missing element. Tries to remove as much from it as
        possible but doesn't search all type_el[...]-lists. Relies on the
        corresponding el_type[...] item.
        """

        # change elsets
        if hasattr(self, 'elset'):
            for els in self.elset.itervalues():
                try:
                    els.remove(el_num)
                except ValueError:
                    pass

        # node list of element nodes
        if remove_nodes:
            try:
                remove_nodes = self.el_node[el_num]
            except (KeyError, AttributeError):
                # el_num not in self.el_node or no attribute self.el_node
                remove_nodes = False

        # remove element
        try:
            thistype = self.el_type[el_num]
            del self.el_type[el_num]
        except (KeyError, AttributeError):
            # el_num not in self.el_type or no attribute self.el_type
            thistype = None
        try:
            self.type_el[thistype].remove(el_num)
            if len(self.type_el[thistype])==0:
                del self.type_el[thistype]
        except (KeyError, AttributeError, ValueError):
            # thistype==None, thistype not in self.type_el, no attribute
            # self.type_el or el_num not in self.type_el[thistype]
            pass
        try:
            del self.el_node[el_num]
        except (KeyError, AttributeError):
            # el_num not in self.el_node or no attribute self.el_node
            pass

        # remove nodes
        if remove_nodes and hasattr(self, 'node_coord'):
            remove_nodes = set(remove_nodes)
            for nodes in self.el_node.itervalues():
                remove_nodes.difference_update(nodes)
                if len(remove_nodes)==0:
                    break
            for node in remove_nodes:
                del self.node_coord[node]
            

        
    ########################################################################
    ### other element specific methods

    eltype_set_tet = set(('C3D4','C3D10M',))
    eltype_set_wedge = set(('COH3D6',))
    eltype_set_tri = set(('S3',))
    eltype_set_box = set(('C3D8R',))

    def calc_centroid(self, element_id):
        """calculate the elements centroid

        bugs: at the moment just take the averange of the corner nodes
        """

        ty = self.el_type[element_id]
        avg_n = 0
        if ty in self.eltype_set_tet:
            avg_n = 4
        elif ty in self.eltype_set_wedge:
            avg_n = 6
        elif ty in self.eltype_set_tri:
            avg_n = 3
        elif ty in self.eltype_set_box:
            avg_n = 8
        else:
            raise Exception('Element type %s not implemented in abq_model_01.'
                            ' Cannot calculate centroid.')

        if avg_n:
            for node in self.el_node[element_id][:avg_n]:
                coords=self.node_coord[node]
                x=x+coords[0]/avg_n
                y=y+coords[1]/avg_n
                z=z+coords[2]/avg_n
            return [x,y,z]


    def get_equiv_tets(self, element_id):
        """
        if the given element is not already a tet (return empty list in this
        case) return tet elements that fill the volume of this element.

        return an empty list also if the given element is not a 3D element.

        The tuple returned is suitable to insert the elements into the model
        quite easily:
        
        result = model.get_equiv_tets(element_id)
        if result:
           for el_tuple in result:
              model.update_elem(*el_tuple)
           model.remove_elem(element_id)

        bugs/limitaions:
        so far works for box elements only
        """
        this_type = self.el_type[element_id]

        if ((this_type in self.eltype_set_tet) or
            (this_type in self.eltype_set_tri)):
            return list()

        if this_type in self.eltype_set_box:
            el_node_idx = [
                (1, 3, 8, 6),
                (1, 3, 6, 2),
                (1, 8, 3, 4),
                (1, 6, 8, 5),
                (3, 8, 6, 7)]
        # eltype_set_wedge
        else:
            raise Exception('Model.get_equiv_tets: Element type %s not'
                            ' implemented yet.')

        result = list()
        last_el_id = max(self.el_node.iterkeys())
        oldnodes = self.el_node[element_id]
        for newtet in el_node_idx:
            nodelist = [oldnodes[n-1] for n in newtet]
            last_el_id += 1
            result.append((last_el_id, 'C3D4', nodelist))

        return result


        
    ########################################################################
    ### other methods extracting parts of the model

    def get_submodel_in_sphere(self, centre, radius):
        """return another model with just all elements that are completely
        within the given sphere.

        elsets and nsets are copied as well
        """
        assert(hasattr(self, 'node_coord')
               and hasattr(self, 'el_node')
               and hasattr(self, 'el_type'))
        
        mo = Model()
        mo.node_coord = NodesDict()
        for node, coord in self.node_coord.iteritems():
            if length(vector(coord,centre))>radius: continue
            mo.node_coord[node] = coord
        
        for elem, nodes in self.el_node.iteritems():
            if all(((nn in mo.node_coord) for nn in nodes)):
                mo.update_elem(elem, self.el_type[elem], nodes)

        connected = self.get_nodes_from_elset()
        
        for node in mo.node_coord.keys():
            if node not in connected:
                del mo.node_coord[node]

        # copy elsets
        if hasattr(self, 'elset'):
            allelems = set(mo.el_node)
            for elsname, elset in self.elset.iteritems():
                inset = allelems.intersection(elset)
                if inset:
                    mo.create_elset(elsname).extend(inset)

        # copy nsets
        if hasattr(self, 'nset'):
            allnodes = set(mo.node_coord)
            for setname, nset in self.nset.iteritems():
                inset = allnodes.intersection(nset)
                if inset:
                    mo.create_nset(setname).extend(inset)

        return mo
