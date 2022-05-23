"""Classes for holding the actual data of an abq_model_02.Model object.
They also do the bulk work when writing an Abaqus input file.

See docs and notes in bae/abq_model_02/__init__.py
"""

# no __version__ info here, it's in __init__.py

import warnings
from collections import defaultdict
from bae.misc_01 import groupIter


class ModelInconsistent(Exception):
    pass


class NodesDict(dict):
    """class for holding the node coordinates

    this is basically a dict, implements write method
    """

    fmtNodenum = "%7d"
    fmtCoord = "%12f"

    @classmethod
    def fmtCoordFunc(self, coord):
        return self.fmtCoord % coord

    def write(self, outputFile, nsetName=None):
        """write the *node command with data in this object
        outputFile is a file object of the already opened abaqus input file

        if nsetName is given, an ,NSET-option is given, that means all nodes
        belong to that nset, e.g. "ALL"
        """
        if len(self)==0:
            return
        if nsetName:
            outputFile.write("*NODE, NSET=%s\n" % nsetName)
        else:
            outputFile.write("*NODE\n")
        ids = self.keys()
        ids.sort()
        for id in ids:
            outputFile.write(
                (self.fmtNodenum % id) + ", "
                + ", ".join(map(self.fmtCoordFunc, self[id])) + "\n")
        return


class ElementNodesDict(dict):
    """class for holding the elements nodes
    {element number: list of node numbers}

    implements write method
    """

    fmtElnum = "%7d"
    fmtNodenum = "%7d"

    def writeElemLines(self, outputFile, elset):
        """
        Write the element definitions of all elements listed in
        elset as datalines in an *ELEMENT-block.

        Important Note: the *ELEMENT-command line is not written! The element
        order in elset is not changed.

        outputFile is a file object of the already opened abaqus input file
        elset is a list or other iterable of elementnumbers (the order of its
        items is not changed)

        silently ignore invalid element numbers
        """
        for element in elset:

            # create linelist to look like that:
            # ...["4132", ", 10", ", 101", ", 11", ...]
            linelist = [(self.fmtElnum % element)+", ",]
            try:
                nodes = self[element]
            except KeyError:
                continue
            for node in nodes[:-1]:
                linelist.append((self.fmtNodenum % node)+", ")
            linelist.append(self.fmtNodenum % nodes[-1])
            for line in groupIter(linelist, maxcounts=8):
                outputFile.write("".join(line)+"\n")

        return

    def update(self, *args):
        """Don't use self.update(), use Model.updateElem() instead.
        If you definitely know what you are doing use dict.update(self, ...)
        """
        raise ModelInconsistent(
            "Don't use ElementNodesDict.update(), this easily leads to"
            " inconsistent data. Instead use Model.updateElem().")


class AbqNameDict(dict):
    """Use instead of dict, when the key is a string in the ABAQUS input file
    converts all keys to uppercase
    """
    def __init__(self, *args, **kwargs):
        dict.__init__(self)
        if len(args)==1:
            if isinstance(args[0], dict):
                for key, value in args[0].iteritems():
                    self[key.upper()] = value
            else:
                for key, value in args[0]:
                    self[key.upper()] = value
        elif len(kwargs)>0:
            for key, value in kwargs.iteritems():
                self[key.upper()] = value

    def __contains__(self, key):
        "convert name to uppercase"
        return dict.__contains__(self, key.upper())

    def __setitem__(self, key, value):
        "convert name to uppercase"
        return dict.__setitem__(self, key.upper(), value)

    def __getitem__(self, key):
        "convert name to uppercase"
        return dict.__getitem__(self, key.upper())

    def __delitem__(self, key):
        "convert name to uppercase"
        return dict.__delitem__(self, key.upper())

    def update(self, iterable):
        if isinstance(iterable, dict):
            iterable = iterable.iteritems()
        for key, val in iterable:
            dict.__setitem__(self, key.upper(), val)

    def updateItem(self, *args, **kwargs):
        """
        Create a new item in the dictionary self if it did not exist so far
        otherwise silently replace the existing. The value of the new item will
        be of type self.typeOfItems. All arguments will be passed to the
        corresponding constructor. The first argument will be used as name/key
        of the new item.

        @param args: the first item of args must be the name (key) of the new
           item
        """
        newItem = self.typeOfItems(*args, **kwargs)
        try:
            name = kwargs["name"]
        except KeyError:
            try:
                name = args[0]
            except IndexError:
                raise ValueError(
                    "Did not find the name argument for the new item %s"
                    % str(newItem))
        self[name] = newItem

    def renameItem(self, oldName, newName):
        """Rename an item.

        This only works for subclasses that store objects with a name attribute
        like Amplitude, Orientation. Otherwise don't use this.

        @Note: If newName already exists it will be silently overwritten.
        """
        try:
            oldItem = self[oldName]
        except KeyError:
            raise ValueError("Could not rename %s-item in the %s-dictionary."
                             "Item not found." % (oldName, self.__class__))
        oldItem.name = newName
        dict.__setitem__(self, newName.upper(), oldItem)
        del self[oldName]

    def writeAll(self, outputFile):
        """
        Write the Abaqus input file commands for all items in this object.
        Order them alphabetically.

        outputFile is a file object of the already opened abaqus input file

        This method assumes that the individual items have a write method.
        Otherwise don't use this method or define a replacement for the derived
        child class.
        """
        ids = self.keys()
        ids.sort()
        for name in ids:
            self[name].write(outputFile)
        return


class BaseSetsDict(AbqNameDict):
    """base class for holding the elsets or nsets
    {set name: set of item labels (element or node numbers)}

    not used directly, the derived classes must have the following class
    members to differentiate between elset and nset:

    setType: "ELSET" or "NSET"
    itemType: "element" or "node"

    implements write methods

    arguments to the constructor:
     - a single argument to initialize the dictionary
     - writeCompressed = True (default): try to create blocks of nodes/elements
       which are then written with the GENERATE flavour of the *ELSET or *NSET
       command.
       self.writeCompressed attribute contains this flag and may be modified.
     - checkItems = True (default): before writing (other operations may be
       added here in the future) check that all items of a set are ints.
       self.checkItems attribute contains this flag and may be modified.

    You may assign other types than sets as values to the dictionary but they
    will possibly converted to a set. A warning is issued at least:
    >>> model.elset["ALL"] = 1
    Single integer given as element set ALL. Will create a set with that number
    as single item.
    >>> model.nset["ALL"] = range(1,11)
    The new node set ALL is not of type set. Converted it to a set.
    """

    fmtItemnum = "%7d"
    minBlockSizeToCompress = 4

    def __init__(self, arg=None, writeCompressed=True, checkItems=True, writeCompressedWithCommaOne=False):
        if arg is None:
            AbqNameDict.__init__(self)
        else:
            AbqNameDict.__init__(self, arg)
        self.writeCompressed = writeCompressed
        self.writeCompressedWithCommaOne = writeCompressedWithCommaOne
        self.checkItems = checkItems

    def __setitem__(self, key, val):
        if type(val)==int:
            warnings.warn("Single integer given as %s set %s. Will create a"
                          " set with that number as single item."
                          % (self.itemType, key.upper()))
            val=set((val,))
        elif type(val)!=set:
            try:
                val=set(val)
                warnings.warn("The new %s set %s is not of type set. Converted"
                              " it to a set."
                              % (self.itemType, key.upper()))
            except TypeError:
                raise TypeError("The new %s set %s is not of type set. Could"
                                " not convert it to a set."
                                % (self.itemType, key.upper()))

        AbqNameDict.__setitem__(self, key, val)

    def writeOne(self, outputFile, setName):
        """write the set with the name setName
        outputFile is a file object of the already opened abaqus input file

        silently ignore empty sets
        """
        if setName not in self:
            raise KeyError(
                "No %s %s, when asked to write it." % (self.setType, setName))
        if len(self[setName])==0:
            outputFile.write("**  %s %s contains no %ss.\n"
                             % (self.setType, setName, self.itemType))
            return

        items = list(self[setName])
        if self.checkItems and not(
                all(isinstance(i, (int, long)) for i in items)):
            raise ModelInconsistent(
                "%s %s contains non integer items."
                % (self.setType, setName))
        items.sort()

        if self.writeCompressed:
            if self.writeCompressedWithCommaOne:
                strGenerateCommaOne = ",1"
            else:
                strGenerateCommaOne = ""
            itemSingles = list()
            itemBlocks = list()
            start = items[0]
            last = start
            offset = start
            for i, eln in enumerate(items):
                supposedEln = i+offset
                if supposedEln != eln:
                    if supposedEln-start<self.minBlockSizeToCompress:
                        itemSingles.extend(range(start, supposedEln))
                    else:
                        itemBlocks.append((start, supposedEln-1))
                    offset = eln-i
                    start = eln
            # write also the last block
            if eln+1-start<self.minBlockSizeToCompress:
                itemSingles.extend(range(start, eln+1))
            else:
                itemBlocks.append((start, eln))

            # now write the blocks
            if len(itemBlocks)>0:
                outputFile.write("*%s, %s=%s, GENERATE\n"
                                 % (self.setType, self.setType, setName))
                for eb in itemBlocks:
                    outputFile.write(
                        (self.fmtItemnum+","+self.fmtItemnum+strGenerateCommaOne+"\n") % eb)
            # prepare to write also the single items
            items = itemSingles

        # write single items
        if len(items)>0:
            outputFile.write("*%s, %s=%s\n"
                             % (self.setType, self.setType, setName))
            for line in groupIter(items,
                                  maxcounts=8, format=self.fmtItemnum):
                outputFile.write(", ".join(line)+"\n")
        return

    def writeAll(self, outputFile):
        """write all sets using the writeOne() method
        outputFile is an file object of the already opened abaqus input file
        silently ignore empty sets
        """
        set_list = self.keys()
        set_list.sort()
        for setName in set_list:
            self.writeOne(outputFile, setName)


class ElsetsDict(BaseSetsDict):
    """
    class for holding the elsets {elset name: set of element numbers}

    implements write methods

    arguments to the constructor:
     - a single argument to initialize the dictionary
     - writeCompressed = True (default): try to create blocks of elements which
       are then written with the *ELSET, GENERATE flavour of the *ELSET command.

    self.writeCompressed attribute contains this flag and may be modified.

    You may assign other types than sets as values to the dictionary but they
    will possibly converted to a set. A warning is issued at least:

    >>> model.elset["ALL"] = 1

    Single integer given as element set ALL. Will create a set with that number
    as single item.

    >>> model.elset["ALL"] = range(1,11)

    The new element set ALL is not of type set. Converted it to a set.
    """

    setType = "ELSET"
    itemType = "element"


class NsetsDict(BaseSetsDict):
    """
    class for holding the nsets {nset name: set of node numbers}

    implements write methods

    arguments to the constructor:
     - a single argument to initialize the dictionary
     - writeCompressed = True (default): try to create blocks of nodes which
       are then written with the *NSET, GENERATE flavour of the *NSET command.

    self.writeCompressed attribute contains this flag and may be modified.

    You may assign other types than sets as values to the dictionary but they
    will possibly converted to a set. A warning is issued at least:

    >>> model.nset["ALL"] = 1

    Single integer given as node set ALL. Will create a set with that number
    as single item.

    >>> model.nset["ALL"] = range(1,11)

    The new node set ALL is not of type set. Converted it to a set.
    """

    setType = "NSET"
    itemType = "node"


class BoundaryCondition(list):
    """
    Class for holding one boundary condition of the type [first dof,
    last dof, magnitude], like at the *BOUNDARY data line.

    Additionally it can store the amplitude and type options.

    Objects of this type are stored in model.boundary of type L{BoundaryDict}

    @ivar type: the string-value of the TYPE parameter of the *BOUNDARY card.
       Typically "DISPLACEMENT", "VELOCITY", "ACCELERATION"
    @ivar amplitude: the string-value of the AMPLITUDE parameter of the
       *BOUNDARY card. Should correspond to a key in model.amplitude
       (L{AmplitudesDict}).

    @Note: Boundary conditions with an ABAQUS boundary type string (like
    "PINNED", "ENCASTRE") are stored as plain strings, not as objects of
    this type. They can not have the additional attributes mentioned above.
    """
    __slots__ = ['type', 'amplitude']

    def __init__(self, vals, type="", amplitude=""):
        if isinstance(vals, basestring):
            raise ValueError(
                "Improper data for BoundaryCondition object. ABAQUS boundary"
                " type strings should be stored as plain strings.")
        list.__init__(self, vals)
        if isinstance(type, basestring):
            self.type = type.upper()
        if isinstance(amplitude, basestring):
            self.amplitude = amplitude.upper()


class BoundaryDict(dict):
    """
    Class for holding the boundary conditions

    It's basically a dict {node number / nset label : bcList}, where bcList is
    a list of L{BoundaryCondition} objects and/or ABAQUS boundary type strings.
    The BoundaryCondition objects are basically lists with additional
    attributes (amplitude and type). The list consists of three values
    like those at the data line of the *BOUNDARY command.

    Example:
     >>> model.boundary.insert("NSET_WALL", "PINNED")
     >>> print model.boundary["NSET_WALL"]
     ['PINNED']
     >>> model.boundary.insert("NSET_NNE", [1, 1, 2.0], type="velocity")
     >>> model.boundary.insert("NSET_NNE", [2, 2, 5.0], amplitude="testamp")
     >>> print model.boundary["NSET_NNE"]
     [[1, 1, 2.0], [2, 2, 5.0]]
     >>> print model.boundary["NSET_NNE"][0]
     [1, 1, 2.0]
     >>> print model.boundary["NSET_NNE"][0].type
     VELOCITY
     >>> print model.boundary["NSET_NNE"][1].amplitude
     TESTAMP

    @note: The keys of this dictionary are converted to uppercase similar to
    what the class AbqNameDict does. But in contrast to this we might also have
    numerical keys here (node numbers).
    """

    fmtNodenum = "%7d"
    fmtBC = "%d, %d, %G"

    def __contains__(self, key):
        """
        convert name to uppercase

        """
        if isinstance(key, basestring):
            key = key.upper()
        return dict.__contains__(self, key)

    def __setitem__(self, key, value):
        "convert name to uppercase"
        if isinstance(key, basestring):
            key = key.upper()
        return dict.__setitem__(self, key, value)

    def __getitem__(self, key):
        "convert name to uppercase"
        if isinstance(key, basestring):
            key = key.upper()
        return dict.__getitem__(self, key)

    def __delitem__(self, key):
        "convert name to uppercase"
        if isinstance(key, basestring):
            key = key.upper()
        return dict.__delitem__(self, key)

    def update(self, iterable):
        if isinstance(iterable, dict):
            iterable = iterable.iteritems()
        for key, val in iterable:
            if isinstance(key, basestring):
                key = key.upper()
            dict.__setitem__(self, key, val)

    def write(self, outputFile):
        """
        write *boundary commands with data in this object

        outputFile is a file object of the already opened abaqus input file

        @note: The order of amplitudes might be rearranged.
        """
        if len(self)>0:
            # optionsIdDict is a dict { (type, amplitude) : (id, index) },
            #   where id is the nset label or node number and index is the
            #   index we need for self[id][index] (index=0 refers to the first
            #   BC stated for this id (i.e. nset or node), index=1 refers to
            #   the second BC for the same id ...).
            optionsIdDict = defaultdict(list)
            for id, values in self.iteritems():
                for indx, value in enumerate(values):
                    try:
                        optionsIdDict[(value.type, value.amplitude)].append(
                            (id,indx))
                    except AttributeError:
                        optionsIdDict[("", "")].append((id,indx))

            optList = optionsIdDict.keys()
            optList.sort()
            for typ, amplitude in optList:
                optstring = ""
                if typ:
                    optstring += ", TYPE=%s" % typ
                if amplitude:
                    optstring += ", AMPLITUDE=%s" % amplitude

                outputFile.write("*BOUNDARY%s\n" % optstring)

                idIndexList = optionsIdDict[(typ, amplitude)]
                idIndexList.sort()
                for id, indx in idIndexList:
                    if isinstance(id, basestring):
                        idString = id
                    else:
                        idString = self.fmtNodenum % id

                    bc = self[id][indx]
                    if isinstance(bc, BoundaryCondition):
                        bc = self.fmtBC % tuple(bc)
                    elif isinstance(bc, basestring):
                        pass
                    else:
                        raise ValueError(
                            "Wrong type <%s> of BC data for node/nset %s in"
                            " the model. Can't write that to input file.\n"
                            % (type(bc), idString))
                    outputFile.write("%s, %s\n" % (idString, bc))
        return

    def insert(self, key, value, type='', amplitude=''):
        """
        insert an existing BC

        @param key: an nset label or a node number

        @param value: Either a string which will be converted to uppercase or a
        list of three values like those at the data line of the *BOUNDARY
        command.

        @note: It is possible to insert contradicting BCs. No check is made if
        there is already a BC for this DOF in the model.
        """
        if isinstance(value, basestring):
            bc=value
        else:
            bc = BoundaryCondition(
                value, type=type, amplitude=amplitude)
        try:
            self[key].append(bc)
        except KeyError:
            self[key] = [bc]


class NodeBoundaryDict(object):
    """class for finding boundary conditions for a specific node
    It works like a dict {node number : bcList}, where bcList is a list of
    ABAQUS boundary type strings and/or BoundaryCondition objects.

    model.nodeBC is an object of this type.

    This class does not store boundary conditions or anything else. It only
    provides the means to retrieve the boundary conditions stored in the
    model for a particular node. It behaves like a dictionary
    {node number: list of boundary conditions}. model.boundary of type
    L{BoundaryDict} on the contrary can have node numbers *or* nset names
    as keys, like the *BOUNDARY card of Abaqus.

    To identify the boundary conditions for a particular node do:
     >>> mynode = 123456
     >>> bcList = model.nodeBC[mynode]

    bcList will now be a list of L{BoundaryCondition} objects and/or ABAQUS
    boundary type strings as stored in model.boundary (L{BoundaryDict}).
    It does not matter if mynode is a key in model.boundary or is part of one
    or more node sets whose names are keys in model.boundary. The result will
    be the combined list of the boundary conditions stored for the particular
    node number (if applicable) and all boundary conditions for all node sets
    mynode is part of.
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
            raise KeyError("No BC found for node %s." % str(key))


class MpcsList(list):
    """class for holding MPC data

    A list of tuples, each tuple consisting of:
      1. string defining the type (first item of the data line)
      2. list of node numbers
    """

    __slots__ = ("mode_option")

    def __init__(self, *args):
        list.__init__(self, *args)
        self.mode_option = "DOF"

    @staticmethod
    def writeOne(outputFile, typ, points):
        firstline = True
        for line in groupIter(points, maxcounts=15, format="%s"):
            if firstline:
                line[0:0] = ["%s" % str(typ),]
                firstline = False
            else:
                line[0:0] = [" 0",]
            outputFile.write(", ".join(line) + "\n")
        return

    def writeAll(self, outputFile):
        userMPC = [dat for dat in self
                   if type(dat[0])==int]
        nonUserMPC = [dat for dat in self
                      if type(dat[0])==str]

        for cmdOption, mpcList in (
                (", MODE=%s, USER" % self.mode_option, userMPC),
                ("", nonUserMPC)):
            if len(mpcList)==0:
                continue
            outputFile.write("*MPC%s\n" % cmdOption)
            for typ, points in mpcList:
                self.writeOne(outputFile, typ, points)
        return


class SurfacesDict(AbqNameDict):
    """Class for holding surfaces.
    It's basically a dict {surface name : surface data}.

    surface data may be a dict {face type : elsetName} in case of surface type
    == "ELEMENT".
    Or it's just a list of data lines, each line represented by a list (or
    tuple) of individual string values (which will be separated by comma in the
    input file).

    Example:
     >>> from bae.abq_model_02 import Model
     >>> model = Model().read("mesh.inp", "NODE,ELEMENT")
     >>> model.surface["Top"] = dict(
     >>>     "S1": set((47, 11)), "S2": set((8, 15)), ... )
     >>> model.write()

    ...Will generate a model with the mesh from mesh.inp plus elsets
    __TOP_S1, __TOP_S2, ... and the surface definition
     >>> *surface, name=Top, type=element
     >>> __TOP_S1, S1
     >>> __TOP_S2, S2
     >>> ...

    Bugs: This implementation is not very good: The type is not stored with the
    surface. The type of self[surfName] is not checked vs. surfType.

    @Note: The interface is very likely to change eventually: surface data
    objects will be of a specialized type incorporating a type attribute.
    For type=ELEMENT this will be derived from dict with very much the same
    data content. Other types need to be implemented separately.

    Surfaces of abitrary type with data lines as sarface data might be dropped
    completely. Because they don't store their type and require the type to
    be stated when exporting.
    """

    def writeOne(self, outputFile, surfName, surfType):
        """Write the *surface command with data in this object.
        @param outputFile: a file object of the already opened abaqus input file
        @param surfName: specifies the name
        @param surfType: the type to appear on the command line in the input
        file.
        """
        surfType = surfType.upper()
        thisSurf = self[surfName]
        if surfType=="ELEMENT":
            thisSurfMod=dict()
            tempElsetsDict=ElsetsDict()
            for faceId, elset in thisSurf.iteritems():
                if type(elset)==set:
                    if len(elset):
                        tempElsetName='__'+surfName+'_'+faceId
                        tempElsetsDict[tempElsetName]=elset
                        thisSurfMod[faceId]=tempElsetName
                else:
                    thisSurfMod[faceId]=elset
            thisSurf = thisSurfMod
            tempElsetsDict.writeAll(outputFile)
        outputFile.write("*SURFACE, NAME=%s, TYPE=%s\n" % (surfName, surfType))
        if surfType=="ELEMENT":
            for faceId, elset in thisSurf.iteritems():
                if faceId=="all exterior faces":
                    outputFile.write("%s\n" % elset)
                else:
                    outputFile.write("%s, %s\n" % (elset, faceId))
        else:
            for line in thisSurf:
                outputFile.write(", ".join(line)+"\n")
        return

    def writeAll(self, outputFile, surfType="ELEMENT"):
        """write the *surface commands for all surfaces in this object
        assuming they all have the same type. Order them alphabetically.

        @param outputFile: a file object of the already opened abaqus input file
        @param surfType: the type to appear on the command line in the input
        file.
        """
        ids = self.keys()
        ids.sort()
        for surfName in ids:
            self.writeOne(outputFile, surfName, surfType)
        return


class SectionProperty(AbqNameDict):
    """Dictionary for all properties of one section definition

    Data lines of the abaqus *... SECTION command are added to the 'data'-item.
    Parameters given as options to the abaqus *... SECTION command are added
    to the dictionary by the self.addSome() function. This creates an item of
    the same name as the abaqus command option (as stated in the attrnames
    argument to the addSome function).

    Items of the dictionary / keyword arguments:

     - TYPE: SOLIDSECTION, BEAMSECTION, COHESIVESECTION, CONNECTORSECTION
     - MATERIAL (required for SOLID-, BEAM-, COHESIVE SECTION): the name of
       the material stored in Model.material.
     - BEHAVIOR (optional for CONNECTOR SECTION): the name of the connector
       behavior stored in Model.connectorBehavior.
     - DATA (optional): a list of data lists, the latter corresponding to
       the data lines in the abaqus input file
     - miscellaneous other options to the *...SECTION command, e.g. CONTROLS,
       SECTION (for beams), ...

    Example: ...yet to be completed...
    """
    __slots__ = []

    def addToSection(self, data):
        """Add the list of parameters data to the 'data'-item of the prop.
        Mainly useful for internal purposes.

        @param data: is a list of values corresponding to one data line of
        the corresponding *...SECTION command
        """
        try:
            datalist = self['data']
        except KeyError:
            datalist = list()
            self['data'] = datalist
        # append the data (one data line of the *... SECTION command)
        self['data'].append(data)

    def write(self, outputFile, elsetName):
        """write this section definition to the abaqus input file

        @param outputFile: write to this file like object
        @param elsetName: this name will appear in the ELSET option
        """
        cmd = "*" + self['TYPE'].replace("SECTION", " SECTION")
        outputFile.write("%s, ELSET=%s" % (cmd, elsetName))
        options = [key for key in self
                   if key not in ('TYPE', 'DATA')]
        for key in options:
            outputFile.write(", %s=%s" % (key, self[key]))
        outputFile.write("\n")
        try:
            dataitems = self['DATA']
        except KeyError:
            return

        try:
            iter(dataitems)
        except TypeError:
            raise TypeError(
                "DATA option of the section property must be iterable."
                " It's not for type %s, elset %s"
                % (self['TYPE'], elsetName))
        for data in dataitems:
            try:
                iter(data)
            except TypeError:
                raise TypeError(
                    "DATA option of the section property must be list of"
                    " lists. For each data line a list of values."
                    " It's not for type %s, elset %s"
                    % (self['TYPE'], elsetName))
            outputFile.write(", ".join(map(str, data)) + "\n")


class PropertiesDict(AbqNameDict):
    """A dictionary {elset: SectionProperty()} with all property definitions
    in the model. Contains L{SectionProperty}-objects.

    Used as L{Model<bae.abq_model_02.Model>}.L{properties<bae.abq_model_02.Model.properties>}

    Usage:
     >>> m = Model()
     >>> m.properties.updateSection(type="SOLIDSECTION", elset="MySet",
     ...                            material="STEEL", data=[[1.0],])
     >>> m.properties.updateSection(
     ...     type="BEAMSECTION", elset="MyBeam", section="CIRC",
     ...     material="STEEL", data=[[radius,], sectionOrientation])
     >>> m.properties.updateSection(
     ...     type="SHELLSECTION", elset="MyShell",
     ...     material="STEEL", data=[[thickness, nbSectionPts],])

    For section properties (*BEAM SECTION, *SHELL SECTION, *SOLID SECTION,
    *CONNECTOR SECTION) the key value is the elset name, the value is an object
    of type L{SectionProperty<bae.abq_model_02.container.SectionProperty>}.

    Each material name should refer to an entry in
    L{Model<bae.abq_model_02.Model>}.L{material<bae.abq_model_02.Model.material>}
    of type L{Material<bae.abq_model_02.container.Material>}.

    Each behavior name should refer to an entry in
    L{Model<bae.abq_model_02.Model>}.L{connectorBehavior<bae.abq_model_02.Model.connectorBehavior>}
    of type L{ConnectorBehaviorDict<bae.abq_model_02.container.ConnectorBehaviorDict>}.

    @ivar matElset: ... dictionary {material name: set of elset names}
    @ivar behaviorElset: ... dictionary {behavior name: set of elset names}

    @Note: Use the L{updateSection}() and L{removeSection}() methods to modify
    this dictionary in order to preserve the structure (i.e. matElset).
    """

    def __init__(self):
        dict.__init__(self)
        self.matElset = defaultdict(set)
        self.behaviorElset = defaultdict(set)

    def updateSection(self, *args, **kwargs):
        """update (add or modify an existing) section definition

        Arguments may be either a dict or you may specify the contents of the
        dict as keyword arguments.

        Items of the dictionary / keyword arguments:

         - type: SOLIDSECTION, BEAMSECTION, COHESIVESECTION, CONNECTORSECTION
         - elset: the elset name
         - material: the name of the material
         - data (optional): a list of data lists, the latter corresponding to
           the data lines in the abaqus input file
         - miscellaneous other options to the *...SECTION command, e.g.
           CONTROLS, BEHAVIOR (connector section), SECTION (beam section)

        Usage:
         >>> m = Model()
         >>> # first possibility
         >>> m.properties.updateSection(type="SOLIDSECTION", elset="MySet",
         ...                            material="STEEL", data=[[1.0],])
         >>> # second possibility (same result)
         >>> m.properties.updateSection({
         ...      "type":"SOLIDSECTION",
         ...      "elset":"MySet",
         ...      "material":"STEEL",
         ...      "data":[[1.0],]})
         >>> # third possibility (same result)
         >>> m.properties.updateSection({
         ...      "type":"SOLIDSECTION",
         ...      "elset":"MySet",
         ...      "material":"STEEL"})
         >>> m.properties["MySet"].addToSection([1.0])
        """
        if len(args)==1:
            obj = args[0]
        elif len(kwargs)>0:
            obj = kwargs
        else:
            raise ValueError(
                "Improper arguments for PropertiesDict.updateSection()."
                "Must be either one positional argument or keyword arguments.")

        # default None. Either behaviorName or materialName will be updated
        behaviorName = None
        materialName = None

        if isinstance(obj, dict):
            # lowercase 'type': argument label; change type to uppercase
            sectionType = obj['type'].upper()
            obj['type'] = sectionType

            prop = SectionProperty(obj)
            elsetName = prop['ELSET']
            if (sectionType == "CONNECTORSECTION"):
                try:
                    behaviorName = prop['BEHAVIOR']
                except KeyError:
                    pass
            else:
                materialName = prop['MATERIAL']

            # remove the item ELSET
            del prop['ELSET']

            # upper all strings
            for key, value in prop.items():
                if isinstance(value, basestring):
                    prop[key] = value.upper()
        else:
            raise ValueError(
                "Improper argument type for PropertiesDict.updateSection(): %s."
                "Must be a dict")

        assert ('TYPE' in prop)

        self[elsetName] = prop
        if elsetName and materialName:
            self.matElset[materialName].add(elsetName)
        if elsetName and behaviorName:
            self.behaviorElset[behaviorName].add(elsetName)
        return prop

    def removeSection(self, elsetNames):
        """Remove a section definition.

        Usage:
         >>> m = Model().read("mymodel.inp")
         >>> m.properties.removeSection("VOLC_FAULTS")
         >>> m.properties.updateSection(
         ...     type="COHESIVESECTION", elset="VOLC_FAULT_CLAUDIA", ...)
         >>> m.properties.updateSection(
         ...     type="COHESIVESECTION", elset="VOLC_FAULT_GREEN", ...)

        Or, to remove all section definitions for elsets starting with "VOLC_":
         >>> m.properties.removeSection(m.regexpNames("elset", "VOLC_"))

        @param elsetNames: Either the elset name as a string or an iterable of
        elset names of the section definition(s) to be removed. Elsets given
        through this parameter that lack a section definition are being
        silently ignored.
        @note: You don't have to remove a section definition to change it as
        long as it is still meant for the same elset. Use updateSection() in
        this case.
        """
        if isinstance(basestring, elsetNames):
            elsetNames = [elsetNames]

        for elsetName in elsetNames:
            try:
                materialName = self[elsetName]['MATERIAL']
            except KeyError:
                pass
            else:
                self.matElset[materialName].remove(elsetName)

            try:
                behaviorName = self[elsetName]['BEHAVIOR']
            except KeyError:
                pass
            else:
                self.matElset[behaviorName].remove(elsetName)

            try:
                del self[elsetName]
            except KeyError:
                pass

    def writeAll(self, outputFile):
        elsetTypeList = sorted(
            (value["TYPE"], key)
            for key, value in self.iteritems())
        for typ, els in elsetTypeList:
            self[els].write(outputFile, els)


class Material(dict):
    r"""Object to store material properties for the *MATERIAL card and
    associates of an Abaqus input file.

    It's a dict with the following optional items:
     - 'DampingAlpha': argument of the *DAMPING, ALPHA=... option
     - 'Density': value for the *DENSITY option
     - 'EMod', 'Nue': first two values of the data line of the *ELASTIC option
     - 'Plastic': list of (stress, strain, temp, fld1, fld2, ...) tuples as
       in the *PLASTIC, HARDENING=ISOTROPIC option
     - 'UmatDelete': argument of the *DEPVAR,DELETE=... option
     - 'UmatDepvarNum': value of the data line of the *DEPVAR card
     - 'Umat': values for the data lines of the *USER MATERIAL card
     - 'UmatType': TYPE-option for the *USER MATERIAL card; can be "THERMAL",
       default (==not given) is "DISPLACEMENT"
     - 'UserDefinedField': exists and == True if the *USER DEFINED FIELD card
       is part of the material definition

    Additionally the Material object stores it's name (the key in the
    L{MaterialDict} repository) in self.name.

    Objects of this type are stored as values in the L{MaterialDict} repository
    L{bae.abq_model_02.Model.material}.

    Adding a new material to a model object:
     >>> model.material.updateItem("MYROCK", {
     >>>     "DampingAlpha": 0.5,
     >>>     "Density": 1000.0,
     >>>     "EMod": 1E6,
     >>>     "Nue": 0.3,
     >>>     "Plastic": [
     >>>         # stress, strain, temp, field1, field2, ...
     >>>         [ 20.0, 0, 0.4],   # not installed (temp=0 < 0.4)
     >>>         [ 540E6, 0, 0.6],  # installed (temp=1 > 0.6), yielding onset
     >>>         [ 635E6, 0.05, 0.6], # inst'ed (temp=1 > 0.6), ultimate force
     >>>         [ 400E6, 0.1, 0.6],  # inst'ed (temp=1 > 0.6), slowly degrading
     >>>         [ 2E6, 0.2, 0.6],    # inst'ed (temp=1 > 0.6), failed
     >>>         ],
     >>>     })
    """
    __slots__ = ['name']

    def __init__(self, name, data=(), **kwargs):
        if data:
            dict.__init__(self, data)
        elif kwargs:
            dict.__init__(self, kwargs)
        else:
            dict.__init__(self)
        self.name = name

    def write(self, outputFile):
        """
        Write the *material command with data in this object.

        @param outputFile: A file object of the already opened abaqus input
          file
        """

        outputFile.write("**\n")
        outputFile.write("**\n")
        outputFile.write("*MATERIAL, NAME=%s\n" % self.name)

        try:
            data = self['DampingAlpha']
        except KeyError:
            pass
        else:
            outputFile.write("*DAMPING, ALPHA=%g\n" % data)

        try:
            data = self['Density']
        except KeyError:
            pass
        else:
            outputFile.write("*DENSITY\n %g\n" % data)

        try:
            data = (self['EMod'], self['Nue'])
        except KeyError:
            pass
        else:
            outputFile.write("*ELASTIC\n %g, %g\n" % data)

        try:
            data = self['Plastic']
        except KeyError:
            pass
        else:
            outputFile.write("*PLASTIC\n")
            # outputFile.write("*PLASTIC, HARDENING=ISOTROPIC\n")
            # write data lines: stress, strain, temp, fld1, fld2, ...
            for line in data:
                if len(line)>8:
                    raise ValueError(
                        "*PLASTIC isotropic hardening: Dependence on more than"
                        " five field variables not implemented so far.")
                outputFile.write(", ".join("%g" % x for x in line)+"\n")

        try:
            data = self['UmatDepvarNum']
        except KeyError:
            pass
        else:
            # check if UmatDelete exists
            try:
                delStr = ", DELETE=%d" % self['UmatDelete']
            except KeyError:
                delStr = ""

            # write *DEPVAR card
            outputFile.write("*DEPVAR%s\n %g\n" % (delStr, data))

        if 'UserDefinedField' in self and self['UserDefinedField']:
            # if the option exists and is True...
            outputFile.write("*USER DEFINED FIELD\n")
            
        try:
            data = self['Umat']
        except KeyError:
            pass
        else:
            # check if type option exists
            try:
                typeStr = ", TYPE=%s" % self['UmatType']
            except KeyError:
                typeStr = ""

            # write umat-card
            outputFile.write(
                "*USER MATERIAL, CONSTANTS=%d%s\n" % (len(data), typeStr))
            try:
                for line in groupIter(data,
                                      maxcounts=8, format="%12.4E"):
                    outputFile.write(", ".join(line)+"\n")
            except TypeError, exc:
                text = (
                    "%s\nFor material %s the list of Umat-values seems to"
                    " not only contain floats:\n%s"
                    % (exc.args[0], self.name, data))
                raise TypeError(text)

        return


class MaterialDict(AbqNameDict):
    """Class for holding material definitions. It's a dict { material name:
    L{Material} object }

    Example:
     >>> model = Model()
     >>> model.material.updateItem("MYROCK", {
     >>>          "DampingAlpha": 0.5,
     >>>          "Density": 1000.0,
     >>>          "EMod": 1E6,
     >>>          "Nue": 0.3,
     >>>          })
     >>> model.write("gaga.inp")

    ...yields gaga.inp containing:
     >>> *MATERIAL, NAME=MYROCK
     >>> *DENSITY
     >>> 1000.0
     >>> *DAMPING, ALPHA=0.5
     >>> *ELASTIC
     >>> 1E6, 0.3
    """
    typeOfItems = Material


class ConnectorBehavior(dict):
    r"""Object to store connector behavior data. It's a dict with the following
    optional items:
     - 'elasticity': a list of (component nb, "rigid") tuples ("rigid" being
       the only option implemented so far)
     - 'plasticity': a list of (component nb, hardening table) tuples. The
       hardening table is a list in itself containing the data lines of the
       *CONNECTOR HARDENING option. Each item of the hardening table is a list
       of floats with all values of one data line. Empty items on the data line
       are denoted by a None value.
     - 'damping' : a list of ( component nb, dataline ) tuples, the
       dataline being a list of floats corresponding to the *CONNECTOR
       DAMPING option: [ damping coefficient (force or moment per relative
       velocity), None, temperature ... see Abaqus manual ]
     - 'stop': a list of (component nb, data lines) tuples. data lines is a
       list in itself containing the data lines of the *CONNECTOR STOP option.
     - 'lock': a list of (component nb, lock, data lines) triples. lock is the
       value of the "lock=" argument of the *CONNECTOR LOCK option. data lines
       is a list in itself containing the data lines of the *CONNECTOR LOCK
       option.

    Additionally the ConnectorBehavior object stores it's name (the key in the
    ConnectorBehaviorDict repository) in self.name.
    """
    __slots__ = ['name']

    def __init__(self, name, data=(), **kwargs):
        if data:
            dict.__init__(self, data)
        elif kwargs:
            dict.__init__(self, kwargs)
        else:
            dict.__init__(self)
        self.name = name

    def _formatDataLine(self, data, fmt="%g"):
        r"""Format a list of floats and/or None's to a string which can be used
        as Abaqus input file data line. Each float is converted with the string
        template given by the "fmt" argument, each None yields an empty value.
        """
        newdata = []
        for x in data:
            if x is None:
                newdata.append("")
            else:
                newdata.append(fmt % x)
        return ", ".join(newdata) + "\n"

    def write(self, outputFile):
        """
        Write the *connector behavior command with data in this object.

        @param outputFile: A file object of the already opened abaqus input
          file
        """

        outputFile.write("*CONNECTOR BEHAVIOR, NAME=%s\n" % self.name)

        try:
            data = self['elasticity']
        except KeyError:
            pass
        else:
            for component, option in data:
                if option == 'rigid':
                    outputFile.write("*CONNECTOR ELASTICITY, COMPONENT=%d,"
                                     " RIGID\n" % component)
                else:
                    warnings.warn(
                        "Connector elasticity other than rigid not implemented"
                        " so far. Ignoring option for connector behavior %s"
                        " component %d. (Found option <%s>)"
                        % (self.name, component, option))

        try:
            data = self['plasticity']
        except KeyError:
            pass
        else:
            for component, hardeningTab in data:
                outputFile.write(
                    "*CONNECTOR PLASTICITY, COMPONENT=%d\n"
                    "*CONNECTOR HARDENING\n"
                    % component)
                for line in hardeningTab:
                    outputFile.write(self._formatDataLine(line))

        try:
            data = self['damping']
        except KeyError:
            pass
        else:
            for component, dampingTab in data:
                outputFile.write(
                    "*CONNECTOR DAMPING, COMPONENT=%d, TYPE=VISCOUS\n"
                    % component)
                for line in dampingTab:
                    outputFile.write(self._formatDataLine(line))

        try:
            data = self['stop']
        except KeyError:
            pass
        else:
            for component, datalines in data:
                outputFile.write(
                    "*CONNECTOR STOP, COMPONENT=%d\n"
                    % component)
                for line in datalines:
                    outputFile.write(self._formatDataLine(line))

        try:
            data = self['lock']
        except KeyError:
            pass
        else:
            for component, lock, datalines in data:
                outputFile.write(
                    "*CONNECTOR LOCK, COMPONENT=%d, LOCK=%s\n"
                    % (component, lock))
                for line in datalines:
                    outputFile.write(self._formatDataLine(line))

        return


class ConnectorBehaviorDict(AbqNameDict):
    r"""Class for holding connector behavior definitions

    It's basically a dict {connector behavior name : L{ConnectorBehavior}
    object}.

    Example:
     >>> model = Model()
     >>> model.connectorBehavior.updateItem(
     >>>    "MYBEHAVIOR",
     >>>    elasticity=[(1, "rigid"),],
     >>>    plasticity=[(1, [
     >>>        [1E8, 0.0],
     >>>        [1E7, 0.01],
     >>>        ]),],
     >>>    stop=[(1, [...missing...]),],
     >>>    lock=[(1, "ALL", [...missing...]),],
     >>>    )
     >>> model.connectorBehavior.updateItem(
     >>>    "DAMPINGBEHAVIOR",
     >>>    damping=[
     >>>        # component, [[damp coeff, None, temp], [c2, None, T2], ...]
     >>>        (i, [[0.5]])
     >>>        for i in 1, 2, 3])
     >>> model.write("gaga.inp")

    ...yields gaga.inp containing:
     >>> *CONNECTOR BEHAVIOR, NAME=MYBEHAVIOR
     >>> *CONNECTOR ELASTICITY, RIGID, COMPONENT=1
     >>> *CONNECTOR PLASTICITY, COMPONENT=1
     >>> *CONNECTOR HARDENING
     >>>  .....
     >>> *CONNECTOR STOP, COMPONENT=1
     >>> , 0.0
     >>> *CONNECTOR LOCK, COMPONENT=1, LOCK=ALL
     >>>  .....
     >>> *CONNECTOR BEHAVIOR, NAME=DAMPINGBEHAVIOR
     >>> *CONNECTOR DAMPING, COMPONENT=1, TYPE=VISCOUS
     >>>  0.5
     >>> *CONNECTOR DAMPING, COMPONENT=2, TYPE=VISCOUS
     >>>  0.5
     >>> *CONNECTOR DAMPING, COMPONENT=3, TYPE=VISCOUS
     >>>  0.5
    """
    typeOfItems = ConnectorBehavior


class Orientation(object):
    r"""
    Class for holding one orientation. Objects of this type are being stored
    in Model.orientation which is of type L{OrientationsDict}

    @ivar a: point a
    @ivar b: point b
    """
    __slots__ = ['name', 'a', 'b']

    def __init__(self, name, a=[], b=[]):
        """
        @param name: name of the orientation in the input file
        @param a: point a (see Abaqus manual)
        @param b: point b
        """
        self.name = name
        self.a = a
        self.b = b

    def write(self, outputFile):
        """
        Write the *orientation command with data in this object.

        @param outputFile: A file object of the already opened abaqus input
          file
        @note: Only DEFINITION=COORDINATES type (default) recognized so far.
        """

        outputFile.write("*ORIENTATION, NAME=%s\n" % self.name)
        outputFile.write(", ".join("%g" % x for x in self.a+self.b)
                         + "\n")
        return


class OrientationsDict(AbqNameDict):
    r"""
    Class for holding orientations.
    It's basically a dict {orientation name : L{Orientation} object}.

    Used as Model.orientation
    """
    typeOfItems = Orientation


class Amplitude(list):
    r"""
    Class for holding one amplitude

    It is a list of (time, amplitude value) tuples with extra attributes
    corresponding to the options of the ABAQUS *AMPLITUDE command and a write
    method.

    @ivar definition: An uppercase string without whitespace. "TABULAR" or
    "SMOOTHSTEP".

    @note: Only DEFINITION=TABULAR or SMOOTH STEP type recognized so far.
    Otherwise the data might be corrupt!
    """
    __slots__ = ['name', 'time', 'smooth', 'definition']

    def __init__(self, name, listInit=[], smooth=None, time=None,
                 definition=None):
        r"""
        Constructor:
        @param name: name of the amplitude in the input file
        @param definition: Defaults to None which is equivalent to TABULAR.
          Should be None, "TABULAR" or "SMOOTH STEP", other values are not
          supported but might not cause error messages! Converted to upper
          case and white space being removed.
        @param listInit: a list of (time, amplitude value) tuples
        @param smooth: A float stating the value of the SMOOTH option of the
          *AMPLITUDE command.
        @param time: Either "STEP TIME" or "TOTAL TIME". See the ABAQUS manual.
          If not specified Abaqus assumes STEP TIME.
        """
        list.__init__(self, listInit)
        self.name = name
        self.smooth = smooth
        self.time = time
        if definition is None:
            self.definition = definition
        else:
            self.definition = definition.upper().replace(" ", "")

    def write(self, outputFile):
        """
        Write the *amplitude command with data in this object.

        @param outputFile: A file object of the already opened abaqus input
          file
        @note: Only DEFINITION=TABULAR or SMOOTH STEP type recognized so far.
        Otherwise the datalines might be corrupt.
        """

        if self.definition == "SMOOTHSTEP":
            defi = ", DEFINITION=SMOOTH STEP"
        elif self.definition:
            defi = ", DEFINITION=%s" % self.definition
        else:
            defi = ""

        try:
            smooth = ", SMOOTH=%g" % float(self.smooth)
        except (ValueError, TypeError):
            smooth = ""

        if self.time:
            ti = ", TIME=%s" % self.time
        else:
            ti = ""

        outputFile.write("*AMPLITUDE, NAME=%s%s%s%s\n"
                         % (self.name, defi, ti, smooth))
        try:
            for quadruple in groupIter(self, 4):
                outputFile.write(", ".join([
                    "%g, %g" % tuple(tiVal)
                    for tiVal in quadruple]) + "\n")
        except TypeError, exc:
            raise TypeError(
                "%s\nThis error occurs if the amplitude list is a flat list"
                " of time and amplitude values instead of a list of (time,"
                " amplitude) tuples." % exc.args[0])
        return


class AmplitudesDict(AbqNameDict):
    r"""
    Class for holding amplitudes.
    It's basically a dict {amplitude name : list of amplitude tuples} with case
    insensitive keys. The amplitude value is of type L{Amplitude} and is
    basically a list of (time, amplitude value) tuples.

    The Model.amplitude attribute is of this type.
    """

    typeOfItems = Amplitude

    # deprecated alias
    def updateAmp(self, *args, **kwargs):
        "Deprecated alias for L{AbqNameDict.updateItem}"
        return AbqNameDict.updateItem(self, *args, **kwargs)


class TimePointsDict(AbqNameDict):
    r"""
    Class for holding time points.
    It's basically a dict {time points name : list of time points}.

    This is the class of the Model.timepoints attribute of the Model class.

    arguments to the constructor:
     - a single argument to initialize the dictionary
     - numberFormat (default='%6.3f'): format of the numbers of the time points
       when being written.

    @ivar numberFormat: Format of the numbers of the time points when being
    written. Initialized by the corresponding argument to the constructor.
    """

    def __init__(self, arg=None, numberFormat='%6.3f'):
        if arg is None:
            AbqNameDict.__init__(self)
        else:
            AbqNameDict.__init__(self, arg)
        self.numberFormat = numberFormat

    def writeOne(self, outputFile, tpName):
        """
        write one time points list to the given output file
        """
        tpList = self[tpName]
        outputFile.write("*TIME POINTS, NAME=%s\n" % tpName)
        for line in groupIter(tpList, maxcounts=8, format=self.numberFormat):
            outputFile.write(", ".join(line)+"\n")

    def writeAll(self, outputFile):
        """
        write the *amplitude commands for all amplitudes in this object.
        Order them alphabetically.

        outputFile is a file object of the already opened abaqus input file
        """
        ids = self.keys()
        ids.sort()
        for name in ids:
            self.writeOne(outputFile, name)
        return

#------------------------------------------------------------------------------
#-- store initial conditions
#------------------------------------------------------------------------------

initCondTypeToClass = dict()
initCondTypeMultLine = set()


class InitialCond(object):
    "Base class for all initial condition classes"

    def __deepcopy__(self, memoDict):
        "For copy.deepcopy"
        return InitialCondTypeTemperature(self.label, self)


class InitialCondTypeSolution(list, InitialCond):
    """
    Class for holding one set of initial conditions of type solution for one
    element or one elset.

    See also the description for class L{InitialCondList}.

    @ivar label:  An elset label (the first value on the first data line)
    """
    __slots__ = ['label']
    type = "SOLUTION"

    def __init__(self, element, listInit=[]):
        r"""
        Constructor:
        @param element: An element number or elset label (the first value on the
        first data line)
        """
        list.__init__(self, listInit)
        self.label = element

    @property
    def element(self):
        """For compatibility reasons. Attribute deprecated."""
        return self.label

    def writeDataLines(self, outputFile):
        """function to write the data line(s) of this initial condition to the
        specified Abaqus input file
        """
        for line in groupIter([self.label,]+self, maxcounts=8, format="%s"):
            outputFile.write(", ".join(line)+"\n")

initCondTypeToClass["SOLUTION"] = InitialCondTypeSolution
initCondTypeMultLine.add("SOLUTION")


class InitialCondTypeStress(list, InitialCond):
    """
    Class for holding one set of initial conditions of type solution for one
    element or one elset.

    See also the description for class L{InitialCondList}.

    @ivar label:  An element number or elset label (the first value on the
       first data line)
    """
    __slots__ = ['label']
    type = "STRESS"

    def __init__(self, element, listInit=[]):
        r"""
        Constructor:
        @param element: An element number or elset label (the first value on the
        first data line)
        """
        list.__init__(self, listInit)
        self.label = element

    @property
    def element(self):
        """For compatibility reasons. Attribute deprecated."""
        return self.label

    def writeDataLines(self, outputFile):
        """function to write the data line(s) of this initial condition to the
        specified Abaqus input file
        """
        outputFile.write("%s, %s\n" % (
            self.label,
            ", ".join("%s" % x for x in self)))

initCondTypeToClass["STRESS"] = InitialCondTypeStress
# initCondTypeMultLine.add()  ... no it's a single data line version


class InitialCondTypeField(float, InitialCond):
    """
    Class for holding one set of initial conditions of type field for one
    node or one nset.

    See also the description for class L{InitialCondList}.

    @ivar label:  An node number or nset label (the first value on the
       first data line)
    """
    __slots__ = ['label']
    type = "FIELD"

    def __new__(cls, node, listInit=[]):
        r"""
        Constructor:
        @param node: A node number or nset label (the first value on the
        first data line
        """
        try:
            # if called with a list argument representing the data line
            # (model.read() does this)
            self = float.__new__(cls, listInit[0])
        except TypeError:
            # if called with any other type of argument (string, float)
            self = float.__new__(cls, listInit)
        self.label = node
        return self

    @property
    def node(self):
        """For compatibility reasons. Attribute deprecated."""
        return self.label

    def writeDataLines(self, outputFile):
        """function to write the data line(s) of this initial condition to the
        specified Abaqus input file
        """
        outputFile.write("%s, %s\n" % (self.label, self))

initCondTypeToClass["FIELD"] = InitialCondTypeField
# initCondTypeMultLine.add()  ... no it's a single data line version


class InitialCondTypeTemperature(float, InitialCond):
    """
    Class for holding one set of initial conditions of type temperature for one
    node or one nset.

    See also the description for class L{InitialCondList}.

    @ivar label:  An node number or nset label (the first value on the
       first data line)
    """
    __slots__ = ['label']
    type = "TEMPERATURE"

    def __new__(cls, node, listInit):
        r"""
        Constructor:
        @param node: A node number or nset label (the first value on the
        first data line
        """
        try:
            # if called with a list argument representing the data line
            # (model.read() does this)
            self = float.__new__(cls, listInit[0])
        except TypeError:
            # if called with any other type of argument (string, float)
            self = float.__new__(cls, listInit)
        self.label = node
        return self

    @property
    def node(self):
        """For compatibility reasons. Attribute deprecated."""
        return self.label

    def writeDataLines(self, outputFile):
        """function to write the data line(s) of this initial condition to the
        specified Abaqus input file
        """
        outputFile.write("%s, %s\n" % (self.label, self))

initCondTypeToClass["TEMPERATURE"] = InitialCondTypeTemperature
# initCondTypeMultLine.add()  ... no it's a single data line version


class InitialCondList(list):
    """
    This list holds the InitialCond items for one abaqus input file command.

    This object also holds all options data and stores the type of the initial
    condition in self.type. All items in the list are objects of the
    corresponding subclass of L{InitialCond}. See the global dict
    L{initCondTypeToClass} that contains the Abaqus type keyword as key and
    the corresponding subclass of L{InitialCond} as value, for example
    L{InitialCondTypeSolution}, L{InitialCondTypeStress},
    L{InitialCondTypeField} or L{InitialCondTypeTemperature}.

    L{InitialCond}-objects typically have a label attribute that identifies the
    elset name or node label or similar. The InitialCond-object itself is also
    derived from an appropriate type to hold the actual data. For example:
    L{InitialCondTypeSolution} is a list of SDV values (floats),
    L{InitialCondTypeStress} is a list of stress components (floats),
    L{InitialCondTypeField} and L{InitialCondTypeTemperature} are floats.

    Example:
     >>> m = Model().read("mymodel.inp")
     >>> for initCondContainer in m.initCond:
     >>>     if initCondContainer.type!="SOLUTION":
     >>>         print "skipping init cond type %s" % initCondContainer.type
     >>>         continue
     >>>     for initCond in initCondContainer:
     >>>         elsetName = initCond.label
     >>>         timing = initCond[16:20]
     >>>         dosomething(elsetName, timing)

    @ivar type: Uppercase string of the Abaqus type keyword. This is the key in
      the global L{initCondTypeToClass} dictionary.
    @ivar options: A dictionary of other keyword-value pairs of the Abaqus
      command. Options without a value (like the GEOSTATIC option of type
      STRESS) will have True as value.
      Note that these items are also duplicated as attributes of the
      L{InitialCondList} object itself.

      Example: Suppose you supply an option 'GEOSTATIC' in the (one and
      only) *INITIAL CONDITION command of mymodel.inp:
       >>> m = Model().read("mymodel.inp")
       >>> icList = m.initCond[0]
       >>> print icList.options
       { "GEOSTATIC":True }
       >>> print icList.GEOSTATIC
       True
    """
    def __init__(self, type, **optionsDict):
        """
        @param type: must be a key in the initCondTypeToClass dictionary. It is
          converted to uppercase and stored in self.type. The list items are
          objects of the corresponding subclass of InitialCond.
        @param optionsDict: All other arguments must be keyword arguments and
          are stored as additional attributes of this object (self) and a
          reference (not a copy!) is stored in self.options.
        """
        self.type = type.upper()
        try:
            self.initCondClass = initCondTypeToClass[self.type]
        except KeyError:
            raise ValueError(
                "*INITIAL CONDITION, TYPE=%s has not been implemented so far.\n"
                % type)
        self.options = optionsDict
        for key, val in optionsDict.iteritems():
            setattr(self, key, val)

    def __deepcopy__(self, memoDict):
        "For copy.deepcopy"
        newList = InitialCondList(self.type, **self.options)
        for initCond in self:
            newList.append(initCond.label, initCond)
        return newList

    def append(self, label, listInit=[]):
        """
        Create a new object of the appropriate subclass of InitialCond and
        append it to the list (to self).

        @param label: element number, elset label, node number or nset label.
        @param listInit: initial condition values. Can be a list in any case
           representing the Abaqus input file data line. For certain types
           (e.g. TEMPERATURE, FIELD) does not need to be a list but can be
           the actual value. This parameter will be passed over to the
           constructor of the corresponding InitialCond class.
        @returns: the new InitialCond object.
        """
        newInitCond = self.initCondClass(label, listInit)
        list.append(self, newInitCond)
        return newInitCond

    def labelValueIterator(self):
        """yield a (label, value) tuple for each item in self.

        label is either the node or element attribute of each item,
        value is the item itself (an InitialCond-object)
        """

        # attr is the name of the label attribute of the corresponding
        # InitCond objects. It is either "node" or "element".
        attr = self.initCondClass.__slots__[0]

        for initCond in self:
            yield (getattr(initCond, attr), initCond)

    def write(self, outputFile):
        """
        Write the *initial condition command with data in this object.

        @param outputFile: A file object of the already opened abaqus input
          file
        """
        optionsString = ", TYPE=%s" % self.type
        for key, val in self.options.iteritems():
            if val is True:
                optionsString += ", %s" % key
            elif val is False:
                pass
            else:
                optionsString += ", %s=%s" % (key, val)

        outputFile.write("*INITIAL CONDITIONS%s\n" % optionsString)
        for initCond in self:
            initCond.writeDataLines(outputFile)

        return


class InitialCondContainer(list):
    """
    This is basically a list of InitialCondList items.
    abq_model_02.Model.initCond is of this type.

    Example:
     >>> model = Model()
     >>> newInitCondList = model.initCond.addNewList("TEMPERATURE")
     >>> newInitCondList.append("INSTALLER", 1.0)
     >>> model.write("gaga.inp")

    ...yields gaga.inp containing:
     >>> *INITIAL CONDITION, TYPE=TEMPERATURE
     >>> INSTALLER, 1.0

    @Note: writing only works for a limited set of varieties. Always check a
    specific type of output (and add a corresponding test function to the test
    suite).
    """

    def addNewList(self, type, **optionsDict):
        """return a new (empty) InitialCond object that is stored in this
        container together with the given options. The options argument must
        contain a "type" item.

        @param type: The type of the initial condition as in the
          initCondTypeToClass dictionary.

        @param optionsDict: All other arguments must be keyword arguments and
          are stored as additional atttributes of this object (self).
        """
        newCondList = InitialCondList(type, **optionsDict)
        self.append(newCondList)
        return newCondList

    def ofType(self, type):
        """Return an iterator, yielding (label, value, initCondList)-tuples of
        all items of self with the specified type.

         - "label" is the element number, elset label, node number or nset label
         - "value" is the InitialCond-object. It might be a list of values as
           specified on the data line of the *INITIAL CONDITION command. Note
           that value is an InitialCond-object that actually contains label as
           either its node or element attribute.
         - "initCondList" is the InitialCondList item this very initial
           condition is stored in. initCondList.options is a dictionary of
           options of the *INITIAL CONDITION command. E.g. might be
           { "GEOSTATIC": True }

        Example:
         >>> initCondTemp = dict(
         >>>     (label, float(val))
         >>>     for label, val, icl in model.initCond.ofType("TEMPERATURE")
         >>>     )
         >>> for node in mynodelist:
         >>>     print "initial temp of node %d:%g"%(node, initCondTemp[node])
        """
        for icl in self:
            if icl.type != type:
                continue
            # attr is the name of the label attribute of the corresponding
            # InitCond objects. It is either "node" or "element".
            attr = icl.initCondClass.__slots__[0]

            for initCond in icl:
                yield (getattr(initCond, attr), initCond, icl)

    def writeAll(self, outputFile):
        "write all initial conditions of the model to the given file"
        for initCondList in self:
            initCondList.write(outputFile)


######################################################
######################################################
###    Classes for making Model a subclass of Mesh
######################################################
######################################################

class FakeElShape(object):
    """Objects of this class look like a dict but they look up
    a different dict.

    It is used in the Model.__init__() function.
    """

    def __init__(self, elType, elTypeToShape):
        self.elType = elType
        self.elTypeToShape = elTypeToShape

    def __getitem__(self, el):
        return self.elTypeToShape[self.elType[el]]

    def __iter__(self):
        return self.elType.__iter__()

    def iterkeys(self):
        return self.__iter__()

    def iteritems(self):
        return ((el, self.elTypeToShape[type_])
                for el, type_ in self.elType.iteritems())


class FakeShapeEl(object):
    """Objects of this class look like a dict but they look up
    a different dict.

    It is used in the L{bae.abq_model_02.Model}.__init__() function.
    L{bae.abq_model_02.Model.shapeEl} is an object of this type.

    Note that L{FakeShapeEl} doesn't have the keys() attribute that ordinary
    dictionaries have. To list all element shapes in the model do the
    following:
     >>> print sorted(model.shapeEl)
    """

    def __init__(self, typeEl, elShapeToTypes, elTypeToShape):
        self.typeEl = typeEl
        self.elShapeToTypes = elShapeToTypes
        self.elTypeToShape = elTypeToShape

    def getTypesFromShape(self, shape):
        """Get element types existing in self corresponding to the given
        element shape.

        This is basically a service function for self.__getitem__() and
        self.get().

        @Returns: A set of Abaqus element types (strings)
        """
        if "_" in shape:
            # shape is a specific linear or quad version
            try:
                types = self.elShapeToTypes[shape]
            except KeyError:
                raise KeyError(
                    "Element shape %s not defined in"
                    " bae.abq_model_02.Model.elShapeToTypes ."
                    " Please update this dictionary." % shape)
        else:
            # shape is a general version
            types = self.elShapeToTypes.get("%s_L" % shape, set()).union(
                self.elShapeToTypes.get("%s_Q" % shape, set()))
            if not types:
                raise KeyError(
                    "Neither element shape %s_L nor %s_Q is defined in"
                    " bae.abq_model_02.Model.elShapeToTypes ."
                    " Please update this dictionary." % (shape, shape))
        return types.intersection(self.typeEl)

    def __getitem__(self, shape):
        """get set of element numbers of the specified shape.
        Both specifically linear or quad versions (like "TET_L") as well as
        the general version (like "TET") might be queried.

        @raises KeyError: if this shape is not present in the model.
        """
        types = self.getTypesFromShape(shape)
        if len(types)==0:
            raise KeyError("No Element of shape %s" % shape)
        elif len(types)==1:
            res = self.typeEl[iter(types).next()]
        else:
            res = set()
            for tt in types:
                res.update(self.typeEl[tt])
        return res

    def get(self, shape, default=None):
        """get set of element numbers of the specified shape.
        Both specifically linear or quad versions (like "TET_L") as well as
        the general version (like "TET") might be queried.

        @returns: a set of elements or the value of the default argument if this
        shape is not present in the model.
        """

        types = self.getTypesFromShape(shape)
        if len(types)==0:
            return default
        elif len(types)==1:
            res = self.typeEl[iter(types).next()]
        else:
            res = set()
            for tt in types:
                res.update(self.typeEl[tt])
        return res

    def __iter__(self):
        shapes = set(self.elTypeToShape[typ] for typ in self.typeEl)
        return iter(shapes)

    def iteritems(self):
        shapeTypes = defaultdict(list)
        for typ in self.typeEl:
            shapeTypes[self.elTypeToShape[typ]].append(typ)
        for shape, types in shapeTypes.iteritems():
            if len(types)==1:
                res = self.typeEl[iter(types).next()]
            else:
                res = set()
                for tt in types:
                    res.update(self.typeEl[tt])
            yield (shape, res)
