"""abq_model_02/model.py

This module is for reading, writing and manipulating abaqus input file data.

Usage
=====
  >>> from bae.abq_model_02 import Model
  >>> m = Model() # an optional parameter may state the already opened log file
  >>> help(m) for more information
"""

###############################################################################
###############################################################################
###############################################################################

# no __version__ string here, version info is kept in __init__.py


from sys import stdout
import re, copy, csv
from itertools import izip, count
from cStringIO import StringIO

from bae.vecmath_01 import length, vector, norm, mat_multvec, \
    vector_rot_x, vector_rot_y, vector_rot_z, vector_minus, \
    vector_modif_add, vector_scale, vector_modif_scale
from bae.future_01 import defaultdict

from bae.mesh_01 import Mesh
from . import container
from .container import ModelInconsistent
from . import reading
from .reading import \
    ReadError, UnionSetInvalidItem, UnknownElType

from bae.misc_01 import quietlog, findDuplicatePoints
from bae.log_01 import msg, MsgTicker

try:
    import numpy as np
except ImportError:
    np = None

######################################################
######################################################
###    Service funtion, check Version of model object
######################################################
######################################################

def checkModelModuleVersion(model, callingfunction):
    """Check that model is a Model object of this module (and not abq_model_01)

    raises an exception if not.

    callingfunction is the name of the calling function (for diagnostic output)

    typical call sequence (from a classes __init__ function) is

    >>> from bae.abq_model_02 import checkModelModuleVersion
    >>> checkModelModuleVersion(model, "%s.%s.__init__" % (__name__, self.__class__.__name__))
    """
#     model_module = model.__class__.__module__.split('.')[-1]
#     this_module = __name__.split('.')[-1]
    model_module = model.__class__.__module__
    this_module = __name__

    if not isinstance(model, Model):
        raise ValueError(
            "%s: The supplied model is an <%s> object from the module %s."
            " It's got to be <Model> from %s."
            % (callingfunction, str(type(model)), model_module, this_module))

    if model_module!=this_module:
        raise ValueError(
            "%s: The supplied model is an object from the module %s."
            " It's got to be from %s."
            % (callingfunction, model_module, this_module))
    return True


######################################################

class Model(Mesh):
    """
    Reads and writes ABAQUS input files, provides the data as member
    attributes and has some functionality to manipulate this data.

    Usage:

      >>> from bae.abq_model_02 import Model
      >>> m = Model()
      >>> m.read(inp_name)  # specify filename *with* extension .inp

    @ivar inpName: File name of the file the data has been read from.
        Will be created by L{read}. Might not be there in all cases.
    @ivar elNodes: dict {element number: [node numbers]}, B{do not modify
        directly},
        use L{updateElems()<bae.abq_model_02.Model.updateElems>},
        L{insertElems()<bae.abq_model_02.Model.insertElems>}.
        L{removeElem()<bae.abq_model_02.Model.removeElem>}.
    @ivar elType: dict {element number: element type string}, B{do not modify
        directly},
        use L{updateElems()<bae.abq_model_02.Model.updateElems>},
        L{insertElems()<bae.abq_model_02.Model.insertElems>}.
        L{removeElem()<bae.abq_model_02.Model.removeElem>}.
    @ivar typeEl: defaultdict {element type string: set of element numbers},
        B{do not modify directly},
        use L{updateElems()<bae.abq_model_02.Model.updateElems>},
        L{insertElems()<bae.abq_model_02.Model.insertElems>}.
        L{removeElem()<bae.abq_model_02.Model.removeElem>}.
    @ivar elShape: pseudo-dict {element number: element shape string},
        read-only attribute,
        type: L{container.FakeElShape<bae.abq_model_02.container.FakeElShape>}
    @ivar shapeEl: pseudo-dict {element shape string: set of element numbers}
        read-only attribute, contains general keys (like "TET") as well as
        specialized (like "TET_Q"),
        type: L{container.FakeShapeEl<bae.abq_model_02.container.FakeShapeEl>}
    @ivar elset: dict {elset name: set of element numbers},
        type: L{container.ElsetsDict<bae.abq_model_02.container.ElsetsDict>}.
    @ivar nset: dict {nset name: set of node numbers},
        type: L{container.NsetsDict<bae.abq_model_02.container.NsetsDict>}.
    @ivar boundary: dict {node number or nset name: list of
        L{BoundaryCondition<bae.abq_model_02.container.BoundaryCondition>}
        objects}, type:
        L{container.BoundaryDict<bae.abq_model_02.container.BoundaryDict>}.
    @ivar nodeBC: dict {node number: },
        type: L{container.NodeBoundaryDict<bae.abq_model_02.container.NodeBoundaryDict>}.
    @ivar nodeCoords: dict {node number: [x,y,z]},
        type: L{container.NodesDict<bae.abq_model_02.container.NodesDict>}.
    @ivar surface: element based surfaces,
        type: L{container.SurfacesDict<bae.abq_model_02.container.SurfacesDict>}.
    @ivar mpc: multi point constraints
        type: L{container.MpcsList<bae.abq_model_02.container.MpcsList>}.
    @ivar properties: *SOLID SECTION, *BEAM SECTION and so on,
        type: L{container.PropertiesDict<bae.abq_model_02.container.PropertiesDict>}.
    @ivar orientation: dict { orientation name: L{container.Orientation<bae.abq_model_02.container.Orientation>} object }
    @ivar amplitude: dict { amplitude name: Amplitude object, a list of
        (time, amplitude)-tuples },
        type: L{container.AmplitudesDict<bae.abq_model_02.container.AmplitudesDict>}.
         >>> m.amplitude.updateAmp("GRAVAMP", [
         ...     (0.0, 0.0), (40.0, 1.0), (50.0, 1.0)],
         ...     smooth=None, time="STEP TIME")
         >>> print m.amplitude["GRAVAMP"]
         [(0.0, 0.0), (40.0, 1.0), (50.0, 1.0)]

    @ivar timepoints: dict { time points name: list of time points }
        type: L{container.TimePointsDict<bae.abq_model_02.container.TimePointsDict>}.
    @ivar initCond: [<type, options, [ <label, [value, value,...]>,...]>,...]
        types:
         - L{container.InitialCondContainer<bae.abq_model_02.container.InitialCondContainer>}:
             a list of ...
         - L{container.InitialCondList<bae.abq_model_02.container.InitialCondList>}:
             a list of childs of InitialCond, e.g. InitialCondTypeSolution
             objects.

    @ivar material: an L{container.MaterialDict} {name : material} created on
        demand (if found in the input file). name is the NAME parameter of the
        *MATERIAL command and material is an
        L{container.Material<bae.abq_model_02.container.Material>} object,
        basically a dictionary that may contain:
         - 'DampingAlpha' : *DAMPING, ALPHA=... (a float)
         - 'Density' : a float
         - 'EMod' : Youngs modulus
         - 'Nue' : Poissons ratio
         - 'Umat' : list of umat parameters

    @ivar connectorBehavior: an L{container.ConnectorBehaviorDict}
        {name : parameter dict}. name is the NAME parameter of the
        *CONNECTOR BEHAVIOR command and the parameter dictionary of type
        L{container.ConnectorBehavior} contains:
         - 'elasticity' : a list of ( component nb, 'rigid' ) tuples (only the
           rigid option has been implemented so far)
         - 'plasticity' : a list of ( component nb, datalines ) tuples, the
           datalines being those of the *CONNECTOR HARDENING option. Each
           dataline being a list of floats: [ friction force, sliding distance,
           ... see Abaqus manual ]
         - 'damping' : a list of ( component nb, dataline ) tuples, the
           dataline being a list of floats corresponding to the *CONNECTOR
           DAMPING option: [ damping coefficient (force or moment per relative
           velocity), None, temperature ... see Abaqus manual ]
         - 'stop' : a list of (component, stop data lines) tuples
         - 'lock' : a list of (component, lock argument, lock data lines)
           triples

    @note: docu on boundary, nodeBC, surface, mpc and properties need to be
        updated.
    """
    #: all known blocks
    recognizedBlocks = set(reading.blockdict.keys())
    # some subsets:
    sectionBlocks = set([key for key in recognizedBlocks
                         if key.endswith("SECTION")])
    del key

    def __init__(self, *args, **kwargs):
        """
        @param args: Specify one postional argument. It should be of type Model
        to initialize this object (use as copy constructor). It might also be
        the logfile argument for compatibility reasons, see below.

        @keyword logfile: deprecated argument, only for compatibility reasons.
        No effect.

        @Note: The interface is a bit weird for compatibility reasons.
        The following is the recommended way of using:
         - supply an (optional) logfile argument as keyword argument.
         - supply an (optional) positional argument of type Model to initialize
           the new model object. (Use as copy contructor)

        @Note: When using as a copy constructor note that the point search
        initialization is not being copied from the old model. This seems to
        not have been missed so far.
        """

        # default if no initializer is specified
        oldmodel = None

        # process arguments
        if len(args) > 1:
            raise ValueError("Model.__init__: may only take one positional"
                             " argument. %d specified." % len(args))
        if (len(args)==1
                and not(isinstance(args[0], Model))
                and hasattr(args[0], "write")):
            # logfile is the first and only postional argument
            self.logfile = args[0]
        elif (len(args)==1 and args[0] is None):
            self.logfile = quietlog
        else:
            try:
                logfile = kwargs["logfile"]
            except KeyError:
                self.logfile = stdout
            else:
                if hasattr(logfile, "write"):
                    self.logfile = logfile
                elif logfile is None:
                    self.logfile = quietlog
                del logfile

            if len(args)==1:
                # now this positional argument can't be a logfile
                oldmodel = args[0]

        # for reading
        self.linecnt = None

        if oldmodel:
            self.inpName = oldmodel.inpName

            # data members
            self.nodeCoords = copy.deepcopy(oldmodel.nodeCoords)
            self.elNodes = copy.deepcopy(oldmodel.elNodes)
            self.elType = copy.deepcopy(oldmodel.elType)
            self.typeEl = copy.deepcopy(oldmodel.typeEl)
            self.elShape = copy.deepcopy(oldmodel.elShape)
            self.shapeEl = copy.deepcopy(oldmodel.shapeEl)

            self.elset = copy.deepcopy(oldmodel.elset)
            self.nset = copy.deepcopy(oldmodel.nset)

            self.boundary = copy.deepcopy(oldmodel.boundary)
            self.nodeBC = copy.deepcopy(oldmodel.nodeBC)

            self.surface = copy.deepcopy(oldmodel.surface)
            self.mpc = copy.deepcopy(oldmodel.mpc)
            self.properties = copy.deepcopy(oldmodel.properties)
            self.orientation = copy.deepcopy(oldmodel.orientation)
            self.amplitude = copy.deepcopy(oldmodel.amplitude)

            self.connectorBehavior = copy.deepcopy(oldmodel.connectorBehavior)
            self.material = copy.deepcopy(oldmodel.material)

            self.timepoints = copy.deepcopy(oldmodel.timepoints)
            self.initCond = copy.deepcopy(oldmodel.initCond)

        else:
            self.inpName = ""

            # data members, they are always there
            self.nodeCoords = container.NodesDict()
            self.elNodes = container.ElementNodesDict()
            self.elType = dict()
            self.typeEl = defaultdict(set)
            self.elShape = container.FakeElShape(
                self.elType, self.elTypeToShapeSep)
            self.shapeEl = container.FakeShapeEl(
                self.typeEl, self.elShapeToTypes, self.elTypeToShapeSep)

            self.elset = container.ElsetsDict()
            self.nset = container.NsetsDict()

            self.boundary = container.BoundaryDict()
            self.nodeBC = container.NodeBoundaryDict(self)

            self.surface = container.SurfacesDict()
            self.mpc = container.MpcsList()
            self.properties = container.PropertiesDict()
            self.orientation = container.OrientationsDict()
            self.amplitude = container.AmplitudesDict()

            self.connectorBehavior = container.ConnectorBehaviorDict()
            self.material = container.MaterialDict()

            self.timepoints = container.TimePointsDict()
            self.initCond = container.InitialCondContainer()

        # additional stuff from the Mesh constructor
        self.point_search_initialized = False

    ######################################################
    #{ consistency checks
    def checkElType(self):
        """Check consistency of the Model object.

        In detail: check if set(elType) == set(elNodes) and check typeEl
        """
        if set(self.elType) != set(self.elNodes):
            raise ModelInconsistent(
                "Different elements in Model.elNodes (#=%d) and Model.elType"
                " (#=%d)." % (len(self.elNodes), len(self.elType)))

        in_typeEl = set()
        for typ, els in self.typeEl.iteritems():
            if len(in_typeEl.intersection(els))!=0:
                raise ModelInconsistent(
                    "%d element(s) found more than once in Model.typeEl"
                    % len(in_typeEl.intersection(els)))
            in_typeEl.update(els)
            if not all([typ == self.elType[el] for el in els]):
                wrongElTypeEls = [el for el in els
                                  if typ != self.elType[el]]
                nbWrong = len(wrongElTypeEls)
                wrongElTypeEls.sort()
                wrongElTypeEls = wrongElTypeEls[:4]
                raise ModelInconsistent(
                    "Inconsistency between typeEl and elType: %d elements in"
                    " typeEl['%s'] have other types: Some examples:\n%s"
                    % (nbWrong, typ,
                       "; ".join(["element %d: %s" % (el,self.elType[el])
                                  for el in wrongElTypeEls])))

        if (in_typeEl != set(self.elType)):
            raise ModelInconsistent(
                "all(typeEl.values()) != set(elType)")
        return
    #} end of consistency checks

    ######################################################
    #{ diagnostic output
    def __str__(self):
        """Some diagnostic output. What's in the model?
        """
        res = []

        nb = len(self.nodeCoords)
        if nb:
            res.append("%d nodes" % nb)

        nb = len(self.elNodes)
        if nb:
            res.append("%d elements" % nb)

        nb = len(self.typeEl)
        if nb>1:
            res.append("%d element types" % nb)

        nb = len(self.elset)
        if nb:
            res.append("%d element sets" % nb)

        nb = len(self.nset)
        if nb:
            res.append("%d node sets" % nb)

        nb = len(self.boundary)
        if nb:
            res.append("%d BCs" % nb)

        nb = len(self.surface)
        if nb:
            res.append("%d surfaces" % nb)

        nb = len(self.mpc)
        if nb:
            res.append("%d MPCs" % nb)

        nb = len(self.properties)
        if nb:
            res.append("%d property/section definitions" % nb)

        nb = len(self.orientation)
        if nb:
            res.append("%d orientations" % nb)

        nb = len(self.amplitude)
        if nb:
            res.append("%d amplitudes" % nb)

        nb = len(self.timepoints)
        if nb:
            res.append("%d sets of time points" % nb)

        nb = len(self.initCond)
        if nb:
            res.append("%d initital conditions" % nb)

        nb = len(self.material)
        if nb:
            res.append("%d material definitions" % nb)

        if res:
            return "Abaqus model data: %s" % (", ".join(res))
        else:
            return "Abaqus model data (empty)"
    #} end of diagnostic output

    ######################################################
    #{ file I/O

    def read(self, inputFile, recognizedBlocks=None):
        r"""Reads an input file, adding the data to this Model object.

        Examples, Usage
        ===============

        >>> from bae.abq_model_02 import Model
        >>> m = Model()
        >>> m.read(inp_name)  # specify filename *with* extension .inp

        to read only *NODE and *ELEMENT - blocks (no warnings about unknown
        commands in this case)
         >>> m.read(inp_name, recognizedBlocks=['NODE', 'ELEMENT'])

        to read all known block types except *BOUNDARY blocks do:
        (duplicate m.recognizedBlocks or this setting applies to all Model
        objects)
         >>> blocks = set(m.recognizedBlocks)
         >>> blocks.remove("BOUNDARY")
         >>> m.read(inp_name, blocks)

        returns self, so the following is possible:
         >>> m = Model().read('abaqusinputfile.inp')


        Recognized keywords, relation to Model object
        =============================================

        Keyword related members::

           ABAQUS keyword | created members of the Model object
           ---------------+------------------------------------
           *NODE          | nodeCoords, possibly nset
           *ELEMENT       | elNodes, elType, typeEl, possibly elset
           *NSET          | nset
           *ELSET         | elset
           *BOUNDARY      | boundary, nodeBC
           *SURFACE       | surface
           *TIME POINTS   | timepoints
           ... needs update


        Interface
        =========

        @param inputFile: May be a filename of
           1. an abaqus input file (including the extension .inp), or
           2. a gzipped abaqus input file (including the extension .gz), or
           3. a zipped abaqus input file (including the extension .zip).
        It might be an open file or any other object that has a next() method
        yielding input lines or is convertible into an appropriate iterator.
        If this object has a name attribute this is used instead of a filename
        for diagnostic output.

        If this open file object has a name attribute that ends in .gz or .zip
        it is treated as an already open compressed file and treated
        accordingly.

        If inputFile is a string not ending in any of the recognized extensions
        .inp, .gz, .zip but beginning with a star ('*'), the inputFile string
        itself is taken as the content to evaluate as input file.

        @param recognizedBlocks: String or list of strings. May contain abaqus
          keywords *without* the leading star, refer to the table above for
          possible ABAQUS keywords. "NODE, TIME POINTS", ["NODE",
          "TIME POINT"], "NODE" would be possible values. Case and all spaces
          are ignored: "Time Points", "time points", "TIMEPOINTS" are all the
          same.

          If NSET is specified, nset definitions from *NSET commands as well as
          from *NODE, NSET=... commands are being recognized. (And still it's
          faster than also reading the node coordinates!) The same applies for
          ELSET definitions in the form *ELEMENT, ELSET=..., accordingly.

        @note:
         - Multiple *NSET commands (or *NODE, NSET=... commands) for the
           same *NSET will subsequently add nodes to the same nset. Same with
           elsets.
         - No check for duplicity is performed. If the same node with its
           coordinates is defined multiple times, all but the last definition
           will be silently ignored. Same with elements.
         - You can perform multiple read commands, new data will be added
           to the model subsequently.
         - All arguments to specific command options are converted to
           uppercase, e.g. NSET-name "gaga" becomes "GAGA" in *NSET, NSET=gaga.
         - *INITIAL CONDITION, TYPE=SOLUTION are not recognized if the elements
           are addressed by element number, only if elset names are used.
         - If the last command is *USER MATERIAL the CONSTANTS parameter is not
           being checked against the actual number of parameters
         - does not process ABAQUS parts/instances input files
           use command mdb.models[modelName].setValues(noPartsInputFile=ON)
           or in cae menu: Model -> Edit Attributes -> Model-1
           mark "do not use parts and assemblies in input file"
         - Abbreviations of keywords not supported. However some command options
           might be abbreviated, e.g. "*BOUNDARY, AMP=...". (Implementation
           note: Look at the OptionList class attribute of the corresponding
           BlockXXX class, search for OptionList in abq_model_02.reading. The
           optional third item in each (keyword, type, ...) - tuple states the
           number of significant letters of the keyword.)
         - only some blocks / ABAQUS-input-file-commands implemented so far;
           anything besides the pure mesh info is not implemented adequately.
         - *SYSTEM is ignored possibly leading to wrong node coordinates.
        """

        # in case it's a zip archive: find the right member
        def openZipArchive(archive):
            archnamelist = [x.filename for x in archive.infolist()
                            if x.filename[-4:].lower() == ".inp"]
            if len(archnamelist)<1:
                raise ReadError(
                    "Zip file <%s> does not contain a .inp file."
                    % self.inpName)
            elif len(archnamelist)>1:
                raise ReadError(
                    "Zip file <%s> contains more than one .inp file."
                    " Reading multiple files in one go not implemented"
                    " (yet)." % self.inpName)
            return archive.open(archnamelist[0], "r")

        if isinstance(inputFile, basestring):
            self.inpName = inputFile
            if self.inpName.lower().endswith(".inp"):
                inp_file = open(self.inpName, "r")
            elif self.inpName.lower().endswith(".inp.gz"):
                import gzip
                inp_file = gzip.open(self.inpName, "rb")
            elif self.inpName.lower().endswith(".zip"):
                import zipfile, os
                if not os.access(self.inpName, os.R_OK):
                    raise ReadError("Input file <%s> is not readable."
                                    % self.inpName)
                if not(zipfile.is_zipfile(self.inpName)):
                    raise ReadError("Input file <%s> is not a valid zip file."
                                    % self.inpName)
                inp_file = openZipArchive(zipfile.ZipFile(self.inpName, "r"))
            elif self.inpName.startswith("*"):
                inp_file = StringIO(self.inpName)
                self.inpName = 'inputstream'
            else:
                raise ReadError("Unrecognized input file name extension <%s>"
                                % self.inpName)
        else:
            if hasattr(inputFile, 'name'):
                self.inpName = inputFile.name
                if self.inpName.lower().endswith(".gz"):
                    import gzip
                    inp_file = gzip.GzipFile(fileobj=inputFile)
                elif self.inpName.lower().endswith(".zip"):
                    import zipfile
                    inp_file = openZipArchive(zipfile.ZipFile(inputFile, "r"))
                else:
                    inp_file = inputFile
            else:
                self.inpName = 'inputstream'
                inp_file = inputFile

        recognizedBlocks = self._checkRecognizedBlocks(recognizedBlocks)

        msg('Start reading input file %s ...' % self.inpName)


        #--- read the File
        self.linecnt = 0
        ticker = MsgTicker("... %d lines read and processed.")
        try:

            # inblock is the object of the readblock that is currently being
            # read superblock is a block object the current position belongs to.
            # e.g. while reading an *ELASTIC block, superblock should be a
            # *MATERIAL block.
            inblock = None
            superblock = None

            # initially: no line from previous readDataLine-call
            line = None
            while 1:

                # we always start a new block here
                # (but we might stay in the current superblock, btw.)
                inblock = None

                # get next line
                while not line:
                    line = inp_file.next()
                    self.linecnt += 1

                    # skip comment line
                    if line[:2] == "**":
                        line = None
                        continue

                    # remove trailing whitespace (possibly leaving nothing)
                    line = line.rstrip()

                # ... found something
                ticker.msg(self.linecnt)
                # print "got line: <%s>" % line   ###### DEBUG

                # identify next block (or find end of file)
                key = reading.Block.read_key(line)

                if superblock:
                    try:
                        inblock = superblock.validSubBlocks[key](
                            line, self, superblock)
                    except KeyError:
                        inblock = None
                        superblock = None

                if not inblock:
                    if (key in recognizedBlocks):
                        if (key=="NODE" and ("NSET" not in recognizedBlocks)
                            or (key=="ELEMENT"
                                and ("ELSET" not in recognizedBlocks))):
                            # special procedure: load node/element definition
                            # ...but ignore sets
                            inblock = reading.blockdict[key](
                                line, self, noSets=True)
                        else:
                            # general procedure: load this block
                            inblock = reading.blockdict[key](line, self)
                        if isinstance(inblock, reading.SuperBlock):
                            superblock = inblock
                    elif ((key=="NODE" and ("NSET" in recognizedBlocks))
                          or (key=="ELEMENT"
                              and ("ELSET" in recognizedBlocks))):
                        # load only the node/element numbers in this block
                        inblock = reading.blockdict[key](
                            line, self, onlySet=True)
                        if key=="NODE":
                            setName = inblock.nset
                        else:
                            setName = inblock.elset
                        if not(setName):
                            inblock = inblock.SkipThisBlock(line, self)
                    elif (key in reading.blockdict):
                        # skip this known block
                        inblock = reading.blockdict[key].SkipThisBlock(
                            line, self)
                    else:
                        # skip this unknown block
                        inblock = reading.SkipBlock(line, self)

                    if (recognizedBlocks==reading.blockdict
                            and (key not in reading.blockdict)):
                        # issue a warning if block not known
                        # only if recognizedBlocks was not specified
                        msg("WARNING: Unknown command:")
                        msg(line)

                #--- read data lines of block
                line, cnt = inblock.readDataLines(inp_file, self)
                # print "read %d data lines of block %s: %s" % (cnt, type(inblock), line)   ###### DEBUG
                self.linecnt += cnt
                ticker.msg(self.linecnt)

        except StopIteration:
            # end of file reached, fini
            pass

        except ReadError, err:
            if inblock:
                msgtext = "ERROR in the %s-block" % inblock.name
            else:
                msgtext = "ERROR"
            msg("%s in %s, line %d:" % (msgtext,self.inpName,self.linecnt))
            msg(err.args[0])
            raise

        del ticker
        msg('Finished reading input file %s.' % self.inpName)
        return self

    def write(self, outputFile, recognizedBlocks=None,
              header="", withSummary=False):
        r"""Writes the model data to an input Abaqus file.

        Examples, Usage
        ===============

        >>> m = Model()
        >>> m.write("mymodel.inp")

        >>> m.write("mymodel.inp", ["NODE", "NSET"])

        >>> fo = open("mynodelist.inp", "w")
        >>> fo.write("** Some extra things to write to this file.")
        >>> m.write(fo, ["NODE", "NSET"])
        >>> fo.close()


        Alternative more flexible methods
        =================================

        There are more flexible methods at hand (defined in
        bae.abq_model_02.container.py) for writing parts of the model
        if you need it, for example:

        >>> m.writeAllElems(fo)
        >>> m.nodeCoords.write(fo)        # write all nodeCoordss
        >>> m.nset.writeOne(fo, ...)      # write just one nset
        >>> m.nset.writeAll(fo)
        >>> m.elset.writeOne(fo, ...)     # write just one elset
        >>> m.elset.writeAll(fo)
        >>> m.boundary.write(fo)          # write all bc's

        However I recommend to create a model with just the data you need and
        then write that whole model. Only the use of the writeAllElems() method
        might be useful if you want to use the ELSET option of the *ELEMENT
        input file command. Though in general because of the compressed
        writing of the elsets even this is of limited use, see below.


        Options for writing ELSET and NSET
        ==================================

        When writing elsets and nsets this function usually tries to create
        blocks of nodes/elements which are then written with the GENERATE
        flavour of the *ELSET or *NSET input file command. This default
        behaviour may be altered:
         >>> m.elset.writeCompressed = False # write each single element number
         >>> m.nset.writeCompressed = False  # write each single node number


        Options for writing TIME POINTS
        ===============================

        When writing time points you can specify a nonstandard format for the
        time point values (otherwise the default is '%6.3f'):
         >>> m.timepoints.numberFormat = '%7.4f'


        Interface
        =========

        @param outputFile: May be an open file object (which is left open by
        this method) or a filename (a string) in which case a file will be
        created and closed afterwards. It silently overwrites an already
        existing file with the same name.

        @param recognizedBlocks: May be an abaqus keyword *without* the leading
        star like "NODE", "ELSET" or a iterable of such strings. See the
        recognizedBlocks argument to the read() method for more info.

        @param header: Optional string to be written before the actual data
        from self and before the optional summary (see withSummary option).
        Can include newline characters for multiple lines. Leading and trailing
        white-space will be chopped off and a newline will be added
        automatically.

        @param withSummary: Add a summary of the contents will be written as a
        comment after the optional header and before the actual data:
        nb of nodes, elements, elsets

        @note: It is not possible to selectivly write e.g. *SOLID SECTION and
        not *BEAM SECTION. If any *...SECTION command is to be written, all of
        them are being written.
        """

        recognizedBlocks = self._checkRecognizedBlocks(recognizedBlocks)

        if type(outputFile) == str:
            outputFile = open(outputFile, "w")
            close_file = True
        else:
            assert hasattr(outputFile, "write")
            close_file = False

        if header:
            outputFile.write(header.strip())
            outputFile.write("\n")

        if withSummary:
            # list contents with only number of items
            contents = list()
            if 'NODE' in recognizedBlocks and self.nodeCoords:
                contents.append("%d nodes" % len(self.nodeCoords))
            if 'ELEMENT' in recognizedBlocks and self.elNodes:
                contents.append("%d elements" % len(self.elNodes))
            if contents:
                outputFile.write(
                    "** Contains %s.\n**\n" % " and ".join(contents))
            del contents

            # list nsets and elsets (with list of items)
            if 'NSET' in recognizedBlocks and self.nset:
                outputFile.write(
                    "** Contains these %d nsets:\n**\n" % len(self.nset))
                for name in sorted(self.nset):
                    nb = len(self.nset[name])
                    outputFile.write("**   %-28s with %7d items\n" % (name, nb))
                outputFile.write("**\n")
            if 'ELSET' in recognizedBlocks and self.elset:
                outputFile.write(
                    "** Contains these %d elsets:\n**\n" % len(self.elset))
                for name in sorted(self.elset):
                    nb = len(self.elset[name])
                    outputFile.write("**   %-28s with %7d items\n" % (name, nb))
                outputFile.write("**\n")

        if 'NODE' in recognizedBlocks:
            assert type(self.nodeCoords) == container.NodesDict
            self.nodeCoords.write(outputFile)

        if 'ELEMENT' in recognizedBlocks:
            self.writeAllElems(outputFile)

        if 'NSET' in recognizedBlocks:
            assert type(self.nset) == container.NsetsDict
            self.nset.writeAll(outputFile)

        if 'ELSET' in recognizedBlocks:
            assert type(self.elset) == container.ElsetsDict
            self.elset.writeAll(outputFile)

        if 'BOUNDARY' in recognizedBlocks:
            assert type(self.boundary) == container.BoundaryDict
            self.boundary.write(outputFile)

        if 'SURFACE' in recognizedBlocks:
            assert type(self.surface) == container.SurfacesDict
            self.surface.writeAll(outputFile)

        if 'MPC' in recognizedBlocks:
            assert type(self.mpc) == container.MpcsList
            self.mpc.writeAll(outputFile)

        # all ...SECTION commands
        if len(self.sectionBlocks.intersection(recognizedBlocks)):
            assert type(self.properties) == container.PropertiesDict
            self.properties.writeAll(outputFile)

        if 'MATERIAL' in recognizedBlocks:
            assert type(self.material) == container.MaterialDict
            self.material.writeAll(outputFile)

        if 'ORIENTATION' in recognizedBlocks:
            assert type(self.orientation) == container.OrientationsDict
            self.orientation.writeAll(outputFile)

        if 'CONNECTORBEHAVIOR' in recognizedBlocks:
            assert type(self.connectorBehavior)==container.ConnectorBehaviorDict
            self.connectorBehavior.writeAll(outputFile)

        if 'AMPLITUDE' in recognizedBlocks:
            assert type(self.amplitude) == container.AmplitudesDict
            self.amplitude.writeAll(outputFile)

        if 'TIMEPOINTS' in recognizedBlocks:
            assert type(self.timepoints) == container.TimePointsDict
            self.timepoints.writeAll(outputFile)

        if 'INITIALCONDITIONS' in recognizedBlocks:
            assert type(self.initCond) == container.InitialCondContainer
            self.initCond.writeAll(outputFile)

        if close_file:
            outputFile.close()

        return

    def writeAllElems(self, outputFile, elsetName=None):
        """
        writes all element definitions in the model to an Abaqus input file

        outputFile is a file object of the already opened Abaqus input file

        if elsetName is given, all elements are assigned to that elset by
        means of an ELSET-option, useful for elset "ALL"
        """
        assert type(self.elNodes) == container.ElementNodesDict

        if elsetName is not None:
            elset_string = ", ELSET=" + elsetName
        else:
            elset_string = ""

        for elType, elems in self.typeEl.iteritems():
            if len(elems)==0:
                continue
            outputFile.write("*ELEMENT, TYPE=%s%s\n" % (elType, elset_string))
            elems = list(elems)
            elems.sort()
            self.elNodes.writeElemLines(outputFile, elems)

        return

    def writeAllElemsOfCertainType(self, outputFile, elType, elsetName=None):
        """writes all element definitions in the model to an Abaqus input file

        @param outputFile: is a file object of the already opened Abaqus input
           file

        @param elType: the elementType, all elements of this type are written
           to outputFile

        @param elsetName: if given, all elements are assigned to that elset by
           means of an ELSET-option, useful for elset "ALL"
        """
        assert type(self.elNodes) == container.ElementNodesDict

        if elsetName is not None:
            elset_string = ", ELSET=" + elsetName
        else:
            elset_string = ""

        try:
            elems = self.typeEl[elType]
        except:
            raise UnknownElType('Element type %s not implemented in %s.'
                                ' Cannot write corresponding elements.'
                                % (elType, self.__module__))

        if len(elems)>0:
            outputFile.write("*ELEMENT, TYPE=%s%s\n" % (elType, elset_string))
            elems = list(elems)
            elems.sort()
            self.elNodes.writeElemLines(outputFile, elems)

        return

    def writeElsetNameCsv(self, outputFile=None, ignoreEmptySets=True):
        """Export all elset names alphabetically sorted to a csv file.

        @param outputFile: Where to write the elset list?

        May be a filename (a string) in which case a file will be
        created and closed afterwards. It silently overwrites an already
        existing file with the same name.

        If None then a generic file name will be deducted from self.inpName
        (the name of the input file read last) This file will be .

        May be an open file object, then this file is left open after writing
        by this method.

        @param ignoreEmptySets: If True (the default) then empty sets will be
        ignored.
        """

        if outputFile is None:
            outputFile = self.inpName.split(".inp")[0]+".csv"

        if isinstance(outputFile, basestring):
            outputFile = open(outputFile, "wb")
            close_file = True
        else:
            assert hasattr(outputFile, "write")
            close_file = False

        output = csv.writer(outputFile)

        if ignoreEmptySets:
            elsetList = [elsetName
                         for elsetName, elset in self.elset.iteritems()
                         if len(elset)>0]
        else:
            elsetList = list(self.elset)

        for elsetName in sorted(elsetList):
            output.writerow([elsetName,])

        if close_file:
            outputFile.close()

        return
    #} end of file I/O

    ######################################################
    #{ edit / modify
    def updateElems(self, elNodes, elType):
        """Inserts new elements into the mesh or replaces existing.

        Use this method if you want to specify the element numbers of the
        elements to be created.

        See also: L{insertElems()<bae.abq_model_02.Model.insertElems>}

        @param elNodes: dictionary {element number: node number
        list} of elements to be updated or inserted
        @param elType: Is the string of the abaqus element type.

        @note: The method updateElems of the base class Mesh takes a shape
        argument (something like "TET_L") instead of the abaqus element type.
        """

        # possibly remove elements from self.typeEl
        overwriteElems = set(elNodes).intersection(self.elNodes)
        for elem in overwriteElems:
            self.typeEl[self.elType[elem]].remove(elem)

        # update Model object attributes
        dict.update(self.elNodes, elNodes)
        self.elType.update(dict.fromkeys(elNodes, elType))
        if elNodes:
            # only touch defaultdict typeEl if there is anything to add
            self.typeEl[elType].update(elNodes)

    def insertNodes(self, coords, firstNode=1):
        """Inserts new nodes.

        This can as well be achieved easily by just assigning a new key, value
        pair to self.nodeCoords. The only benefit of this method is that it
        looks for free node numbers starting with firstNode. This avoids
        overwriting existing nodes.

        The method copies the coord triples into a new list (to make sure it is
        a list, inserting e.g. tuples into self.nodeCoords might create very
        hard to find errors). As a side effect, the coords list may be altered
        afterwards without effect on the model self.

        @param coords: List or other iterable of coordinate triples.
        @param firstNode: first possible node number.

        @returns: List of newly created node numbers.
        """
        insertedNodes = list()
        for pt in coords:
            while firstNode in self.nodeCoords:
                firstNode+=1
            self.nodeCoords[firstNode] = list(pt)
            insertedNodes.append(firstNode)
            firstNode+=1
        return insertedNodes

    def insertElems(self, elNodes, elType, firstElem=1):
        """Inserts new elements. This method never replaces existing elements.

        Use this method if you don't want to specify the element numbers of the
        elements to be created.

        See also: L{updateElems()<bae.abq_model_02.Model.updateElems>}

        @param elNodes: list or iterable of node number lists. One node number
        list per elements to be inserted into the model

        @param elType: Is the string of the abaqus element type.

        @param firstElem: Lowest possible element number for the elements to be
        created. If there is no element with the element number firstElem in
        the model yet, then the first new element will get this number.
        Otherwise the next higher number not being used already will be chosen.

        @Returns: a list of element numbers of the newly creted elements.

        @Raises ValueError: if the items of elNodes are not lists or if it's
        not possible to convert each item to a list.

        @note: There should be a corresponding method defined in the base
        class Mesh taking a shape argument (something like "TET_L") instead of
        the abaqus element type. It's not yet implemented.
        """

        elNodesDict = dict()
        inserted = list()
        for nodes in elNodes:
            while firstElem in self.elNodes:
                firstElem += 1
            inserted.append(firstElem)
            if not isinstance(nodes, list):
                try:
                    nodes = list(nodes)
                except TypeError:
                    raise ValueError(
                        "The elNodes argument of Model.insertElems() must be a"
                        " list of node lists.")
            elNodesDict[firstElem] = nodes
            firstElem += 1

        # update Model object attributes
        dict.update(self.elNodes, elNodesDict)
        self.elType.update(dict.fromkeys(inserted, elType))
        self.typeEl[elType].update(inserted)

        return inserted

    def insertModel(self, otherModel, mergeNodes={}, recognizedBlocks=None,
                    orientationTolerance=0.01):
        """Add the contents of another model (nodes, elements, elsets, ...)
        to self.

        This has been developed with ground support models in mind: Create a
        bolt model, duplicate and modify the bolt model then add it to the all
        encompassing ground support model that already includes the gs mesh.

        Node and element number of the otherModel will be conserved if there is
        no node resp. element in self falling in between the those new numbers.
        I.e. if there is node 1, 2, 10, 11 in self and nodes 3, 4, 5 in
        otherModel then they will retain their numbers upon import to self and
        will be added to self as nodes 3, 4, 5. If there is again node 1, 2,
        10, 11 in self and now nodes with numbers 7, 12, 17 in otherModel then
        they will get new numbers and become the new nodes 12, 13, 14 in self.
        (The same applies to element numbers accordingly.)

        Nodes that are not defined in otherModel.nodeCoords but used in the
        element definitions otherModel.elNodes will keep their old numbers.
        This is a feature to enable import of elements connected to already
        defined nodes with this function. You might want to check beforehand
        if all nodes are well defined if you don't need this feature.

        Node coordinates are deep copied, i.e. if you later change individual
        node coordinates in otherModel then the node coordinates of self are
        not affected. Same with the element connectivity elNodes.

        @param otherModel: Model-instance to be added to self.

        @param mergeNodes: A dictionary {nodeOther : nodeThis}. nodeOther
          is the node number in otherModel of the node to be replaced by
          the nodeThis-node in self.

          I.e. say mergeNodes={1000:10, ...}. Node 10 is a node in self. All
          occurences of node 1000 in otherModel are replaced by 10 during the
          process of adding otherModel to self. The node 1000 from otherModel
          won't be imported to self.

        @param recognizedBlocks: String or list of strings. May contain abaqus
          keywords *without* the leading star. "NODE, TIME POINTS", ["NODE",
          "TIME POINT"], "NODE" would be possible values. Case and all spaces
          are ignored: "Time Points", "time points", "TIMEPOINTS" are all the
          same.

        @param orientationTolerance: There are two axes a and b defining each
          orientation. When comparing two (possibly identical or similar)
          orientation the (vector-)difference between the axes a is being
          computed. It's length is approximately the angle between the two a's
          in radians as long as the angle is small. The corresponding length
          is being calculated for the b axes as well. If the sum of those two
          lengths (approximating the angles) is smaller than this parameter
          orientationTolerance then the two orientations are considered
          equal.

          A number of 0.01 approximately corresponds to the sum of the angles
          between the two corresponding axes being less than 0.5 degrees.

        @Return: A tuple (nodesOldToNew, elemsOldToNew):
          - nodesOldToNew: {node number in otherModel, new node num}
          - elemsOldToNew: {elNum in otherModel, new elNum}

        @Note: Currently only the following items are recognized (as may be
          stated in the recognizedBlocks argument):
          NODE, ELEMENT, NSET, ELSET, MATERIAL, ORIENTATION, CONNECTORBEHAVIOR,
          (SOLID, COHESIVE, SHELL, BEAM, CONNECTOR)-SECTION

        @Note: For material, orientation and connector behavior only new items
            will be inserted. Items with a name that already exists in self are
            not being imported.

        @Note: Of the material section definitions only those will be imported
            that refer to an elset that does not already have a section
            definition in self.
        """

        recognizedBlocks = self._checkRecognizedBlocks(recognizedBlocks)

        otherNodesToAdd = sorted( node for node in otherModel.nodeCoords
                                  if node not in mergeNodes )

        # determine new node numbers
        if len(otherNodesToAdd)==0:
            newNodes = []
        elif set(self.nodeCoords).intersection(
                range(otherNodesToAdd[0], otherNodesToAdd[-1]+1)):
            # conflicting node number ranges
            firstNode = max(self.nodeCoords)+1
            newNodes = range(firstNode, firstNode+len(otherNodesToAdd))
        else:
            newNodes = otherNodesToAdd

        # insert nodes
        if ("NODE" in recognizedBlocks):
            self.nodeCoords.update(
                (newN, otherModel.nodeCoords[oldN][:])
                for newN, oldN in izip(newNodes, otherNodesToAdd))

        # initialize nodesOldToNew: {node number in otherModel, new node num}
        nodesOldToNew = dict(mergeNodes)
        nodesOldToNew.update(izip(otherNodesToAdd, newNodes))
        del otherNodesToAdd, newNodes

        # update nsets
        if ("NSET" in recognizedBlocks):
            for key, items in otherModel.nset.iteritems():
                self.forceNset(key).update(
                    nodesOldToNew[node] for node in items)

        # determine new element numbers
        oldElems = sorted(otherModel.elNodes)
        if set(self.elNodes).intersection(
                range(oldElems[0], oldElems[-1]+1)):
            # conflicting element number ranges
            firstElem = max(self.elNodes)+1
            newElems = range(firstElem, firstElem+len(oldElems))
        else:
            newElems = oldElems

        # initialize elemsOldToNew: {elNum in otherModel, new elNum}
        elemsOldToNew = dict(izip(oldElems, newElems))
        del oldElems, newElems

        # insert elements
        if ("ELEMENT" in recognizedBlocks):
            for elType, elems in otherModel.typeEl.iteritems():
                elNodesDict = dict(
                    (elemsOldToNew[elem],
                     [nodesOldToNew.get(node, node)
                      for node in otherModel.elNodes[elem]])
                    for elem in elems)
                # update Model object attributes
                dict.update(self.elNodes, elNodesDict)
                self.elType.update(dict.fromkeys(elNodesDict, elType))
                self.typeEl[elType].update(elNodesDict)

        # update elsets
        if ("ELSET" in recognizedBlocks):
            for key, items in otherModel.elset.iteritems():
                self.forceElset(key).update(
                    elemsOldToNew[elem] for elem in items)

        # update material section definitions (all ...SECTION commands)
        if (len(self.sectionBlocks.intersection(recognizedBlocks))
                and len(otherModel.properties)):
            assert (type(self.properties) == type(otherModel.properties)
                    == container.PropertiesDict)
            newData = set(otherModel.properties).difference(self.properties)
            for name in newData:
                self.properties[name] = copy.deepcopy(
                    otherModel.properties[name])

        # update material
        if 'MATERIAL' in recognizedBlocks and len(otherModel.material):
            newData = set(otherModel.material).difference(self.material)
            for name in newData:
                self.material[name] = copy.deepcopy(
                    otherModel.material[name])

        # update orientations
        if 'ORIENTATION' in recognizedBlocks and len(otherModel.orientation):
            for name, ori in otherModel.orientation.iteritems():
                try:
                    oldOri = self.orientation[name]
                except KeyError:
                    # orientation with new name will be copied
                    self.orientation[name] = copy.deepcopy(ori)
                else:
                    # orientation with already existing name
                    if ((length(vector(norm(ori.a), norm(oldOri.a)))
                         + length(vector(norm(ori.b), norm(oldOri.b))))
                            > orientationTolerance):
                        # orientations are different
                        msg("WARNING: There is a different orientation of the"
                            " same name %s in the model to be inserted. It"
                            " won't be inserted, the old definition prevails."
                            % name)
                        # idea: store under a different name.
                        # problem: need to update section definition and then
                        # possibly rename that as well...
                        # newname = incrementName(name) # to be found in misc_01
                        # self.orientation[newname] = copy.deepcopy(ori)
                        # oldToNewOri[name] = newname

        # update connectorBehavior
        if ('CONNECTORBEHAVIOR' in recognizedBlocks
                and len(otherModel.connectorBehavior)):
            newData = set(otherModel.connectorBehavior).difference(
                self.connectorBehavior)
            for name in newData:
                self.connectorBehavior[name] = copy.deepcopy(
                    otherModel.connectorBehavior[name])

        # update initial conditions
        if 'INITIALCONDITIONS' in recognizedBlocks:
            for newInitCondList in otherModel.initCond:

                # which (node/element)-labels have already got an initial
                # condition in self
                alreadyDefinedLabels = dict(
                    (label, (value, icl))
                    for label, value, icl
                    in self.initCond.ofType(newInitCondList.type))
                for newICLabel, newICValue in newInitCondList.labelValueIterator():
                    if newICLabel in alreadyDefinedLabels:
                        if len(newInitCondList.options):
                            msg("WARNING: When inserting the new model the"
                                " initial condition of type %s is being"
                                " ignored for node/element %s because it's"
                                " got additional options. Feature not"
                                " implemented yet."
                                % (newInitCondList.type, newICLabel))
                            continue
                        elif len(alreadyDefinedLabels[newICLabel][1].options):
                            msg("WARNING: When inserting the new model the"
                                " initial condition of type %s is being"
                                " ignored for node/element %s because there"
                                " is already a corresponding initial"
                                " condition with options in self. Feature"
                                " not implemented yet."
                                % (newInitCondList.type, newICLabel))
                            continue
                        else:
                            # no warning if there is no option involved
                            # new initial condition is silently discarded
                            continue
                    else:  # newICLabel not in alreadyDefinedLabels
                        try:
                            oldInitCondList = alreadyDefinedLabels[newICLabel][1]
                        except KeyError:
                            oldInitCondList = self.initCond.addNewList(
                                newInitCondList.type, **newInitCondList.options)

                        # store additional initial condition in the
                        # corresponding list
                        oldInitCondList.append(
                            newICLabel, copy.deepcopy(newICValue))

        # return some of the transfer parameters
        return (nodesOldToNew, elemsOldToNew)

    def updateElem(self, elNum, elType, nodes):
        """Try to avoid this function, it may be removed, instead use
        L{updateElems()<bae.abq_model_02.Model.updateElems>}
        if you want to specify element numbers one by one or
        L{insertElems()<bae.abq_model_02.Model.insertElems>}
        if you don't want to specify each element number. The
        interface of updateElem is messed up somehow.

        Inserts a new element into the model or replaces an existing one.

        Updates elNodes, elType, typeEl dicts.

        See also: L{insertElems()<bae.abq_model_02.Model.insertElems>}

        @param nodes: It's either a list (or iterable) of node numbers or a
        list of node number lists or a dictionary {element number: node number
        list}

        If nodes is a list (or iterable) and its first item is not a list then
        nodes is treated as the node list for one single element.

        If nodes is a list (or iterable) and its first item is a list then
        nodes is treated as a list of node lists for as many elements of type
        elType as there are node lists in the nodes argument.

        If nodes is a dictionary it is interpreted as {element number : node
        list}. In this case elNum ist being ignored and element numbers and
        corresponding node lists are taken from the dictionary.

        @param elType: Is the string of the abaqus element type.

        @param elNum: Is the first element number used by replacing
        or adding one or more elements. If elNum == None, select next free
        element number. If many elements are to be inserted (i.e. if nodes
        is a list of node lists) elNum is the number of the first element,
        elNum+1 the next and so on.

        @note:
        Silently replaces an existing element with the same number. (This is
        also removing the element from the old self.typeEl set if the type
        changes, so it's save to use this method to modify an element.)

        @note:
        The actual list of nodes suplied in either way through the nodes
        argument are not copied but referenced. If you modify those lists it
        will affect self's connectivity as well.

        @note:
        If you want to update / insert many elements then it is way(*) faster
        to first create a list of node number lists or a dictionary {element
        number: list of node numbers} and supply this to the nodes argument,
        thus calling this method only once. (Or once for each element type.)

        (*) "way faster" means: It is not feasible to insert 1 Mio Elements by
        calling this method for each element, it would take hours. Creating a
        dict or list and then calling this method once is a matter of seconds.

        @return: First used element number (elNum if specified)
        """

        if isinstance(nodes, dict):
            newElNodes = nodes
        else:
            if elNum is None:
                try:
                    elNum = max(self.elNodes) + 1
                except ValueError:
                    # if len(self.elNodes)==0: "max() arg is an empty sequence"
                    elNum = 1

            elNodesIterator = iter(nodes)
            thisNodes = elNodesIterator.next()
            if isinstance(thisNodes, list):
                newElNodes = dict(izip(count(elNum+1), elNodesIterator))
                newElNodes[elNum] = thisNodes
            elif isinstance(thisNodes, (int, long)) and isinstance(nodes, list):
                newElNodes = {elNum: nodes}
            elif isinstance(thisNodes, (int, long)):
                newElNodes = {elNum: list(nodes)}
            else:
                raise ValueError(
                    "Invalid type %s for nodes argument to Model.updateElem()"
                    % type(nodes))

        self.updateElems(newElNodes, elType)
        return elNum

    def removeElems(self, elems, removeNodes=False, updateElsets=True):
        """Removes elements from the model.
        Updates elNodes, elType, typeEl dicts and all elsets.

        @param elems: A list or a set of element numbers (also accepts a dict,
        in that case takes its keys.)

        @param removeNodes: if True remove unused nodes. This is not very
        efficient if many elements are beeing deleted, in this case rather
        clean up after all elements have been removed using
        getNodesFromElset():
          >>> model.removeElems(elems=myBigElset)
          >>> allNodes = model.getConnectedNodes()
          >>> removeNodes = set(model.nodeCoords).difference(allNodes)
          >>> for node in removeNodes:
          ...    del model.nodeCoords[node]

        @param updateElsets: If True (default) also remove the element from
        all elsets

        @note: Silently ignores missing elements.

        @note: It is save to do the following (it is not only save but the
        recommended way to accomplish those tasks):
          >>> model.removeElems(model.elNodes)  # remove all elements
          >>> model.removeElems(model.typeEl['S3']) # remove one type

        @Note: Historically abq_model_02.Model.removeElem was first. Next was
        mesh_01.Mesh.removeElems (note the trailing "s"). The latter was
        simpler and presumably faster. There was a plcaeholder that prevented
        unintended use of Mesh.removeElems by a Model-instance.
        """

        # special treatment if you want to remove all, wouldn't work otherwise
        if elems is self.elNodes:
            self.elNodes.clear()
            self.typeEl.clear()
            self.elType.clear()
            if updateElsets:
                self.elset.clear()
            return

        # change elsets
        # remove elements from all elsets
        if updateElsets:
            for els in self.elset.itervalues():
                if els is elems:
                    # special case: elems is one of the elsets
                    # make a copy before clearing this set
                    elems = set(elems)
                    els.clear()
                els.difference_update(elems)

        # node list of element nodes
        if removeNodes:
            removeNodesSet = set()
            for elNum in elems:
                try:
                    removeNodesSet.update(self.elNodes[elNum])
                except KeyError:
                    # elNum not in self.elNodes
                    pass

        # remove elements from self.elNodes and self.elType
        for elNum in elems:
            try:
                del self.elType[elNum]
            except KeyError:
                pass

            try:
                del self.elNodes[elNum]
            except KeyError:
                pass

        # remove elements from self.typeEl
        # use a copy (self.typeEl.items()) in the loop
        for thistype, els in self.typeEl.items():
            els.difference_update(elems)
            if len(els)==0:
                del self.typeEl[thistype]

        # remove nodes
        if removeNodes:
            for nodes in self.elNodes.itervalues():
                removeNodesSet.difference_update(nodes)
                if len(removeNodesSet)==0:
                    break
            for node in removeNodesSet:
                del self.nodeCoords[node]

    def removeElem(self, elem, removeNodes=False, updateElsets=True):
        """DEPRECATED! Use L{removeElems} instead.
        Removes elements from the model.
        Updates elNodes, elType, typeEl dicts and all elsets.

        @param elem: Might be a single element number, a list or a set of
        element numbers (also accepts a dict, in that case takes its keys.)

        @param removeNodes: if True remove unused nodes. This is not very
        efficient if many elements are beeing deleted, in this case rather
        clean up after all elements have been removed using
        getNodesFromElset():
          >>> model.removeElem(elem=myBigElset)
          >>> allNodes = model.getNodesFromElset()
          >>> removeNodes = set(model.nodeCoords).difference(allNodes)
          >>> for node in removeNodes:
          ...    del model.nodeCoords[node]

        @param updateElsets: If True (default) also remove the element from
        all elsets

        @note: Silently ignores missing elements.
        @note: It is save to do the following (it is not only save but the
        recommended way to accomplish those tasks):
          >>> model.removeElem(elem=model.elNodes)  # remove all elements
          >>> model.removeElem(elem=model.typeEl['S3']) # remove one type
        """

        if type(elem) == int:
            elem = [elem,]

        Model.removeElems(self, elem, removeNodes, updateElsets)

    def mergeNodesOnElems(self, elset=None, tolerance=0.01):
        """Reconnects elements to one single node if nodes are very close
        together.

        To get rid of nodes no longer needed do the following. This assumes
        that none of the merged nodes is used in any other element.
        Use with Caution!
          >>> oldToNew = model.mergeNodesOnElems()
          >>> for node in oldToNew:
          ...    del model.nodeCoords[node]

        If you there might be other elements still be attached to the same
        nodes then use the following alternative to get rid of nodes that are
        no longer connected to any element:
          >>> model.mergeNodesOnElems()
          >>> allNodes = model.getNodesFromElset()
          >>> removeNodes = set(model.nodeCoords).difference(allNodes)
          >>> for node in removeNodes:
          ...    del model.nodeCoords[node]

        For diagnostic output:
          >>> oldToNew = model.mergeNodesOnElems()
          >>> nbNew = len(set(oldToNew.itervalues()))
          >>> nbOld = len(oldToNew)+nbNew
          >>> msg("Merged %d nodes into %d." % (nbOld, nbNew))

        @param elset: may be an element set name or a list of element set
           names (everything accepted by L{getUnionSet()}. If not given
           defaults to the whole model.
        @param tolerance: Approximate tolerance for two nodes considered
           identical. For a more precise explanation look at
           L{findDuplicatePoints<bae.misc_01.findDuplicatePoints>}.
        @returns: oldToNew a dict {old node number : new node number}
        """

        if elset is None:
            elset = set(self.elNodes)
        else:
            elset = self.getUnionSet("elset", elset)

        # find duplicate nodes
        allNodes = list(self.getNodesFromElset(elset))
        coords = [self.nodeCoords[n] for n in allNodes]
        samePts = findDuplicatePoints(coords, tolerance=tolerance)
        oldToNew = dict(
            (allNodes[other], allNodes[ids[0]])
            for ids in samePts
            for other in ids[1:])

        # update element connectivity
        for elem in elset:
            self.elNodes[elem] = [
                oldToNew.get(node, node) for node in self.elNodes[elem]]

        return oldToNew

    def updateSurf(self, surfacename, surfaceelementlists):
        """Create a new surface from scratch, silently overwriting an already
        existing one with the same name.

        not implemented yet
        """
        raise Exception("Model.updateSurf() not implemented yet.")
        return

    def forceElset(self, elsetName):
        """return self.elset[elsetName] or create a new elset with that name if
        it didn't exist so far.

        This is very handy for adding elements to the set without prior check
        if the set already exists, e.g.:
          >>> m = Model()
          >>> newSets = {"ELSET_1":(45, 711, 7612), "ELSET_2":(1, 204, 1877)}
          >>> for elset, element in newSets.iteritems():
          >>>     m.forceElset(elset).update(element)
        """
        try:
            return self.elset[elsetName]
        except KeyError:
            newset = set()
            self.elset[elsetName] = newset
            return newset

    def forceNset(self, nsetName):
        """return self.nset[nsetName] or create a new nset with that name if
        it didn't exist so far."""
        try:
            return self.nset[nsetName]
        except KeyError:
            newset = set()
            self.nset[nsetName] = newset
            return newset

    def renumber(self, nodeStart=1, elemStart=1):
        """renumber elements and nodes

        Examples:
          >>> # renumber node and element numbers consecutively starting at 1
          >>> model.renumber()

          >>> # renumber only the nodes to start at 500001
          >>> model.renumber(nodeStart=500001, elemStart=None)

        @param nodeStart: New start for node numbers. If None, nodes are not
        renumbered.
        @param elemStart: New start for element numbers. If None, elements are
        not renumbered.

        @Note: The relative order of nodes and elements is not touched.

        @Note: If the element numbers are touched (elemStart != None), the
        point search must eventually be initialized again. This could be changed
        (i.e. such that this method also updates the point search initialization
        but this is not done yet.)

        @Note: At the moment besides the actual node and element information
        only nsets and elsets are updated. I.e. neither boundary conditions,
        surface definitions nor all the rest.

        @Returns: Tuple of two dictionaries (nodesOldToNew, elemsOldToNew)
        maping old node/element numbers to the new ones.

        @Raises KeyError: If a node or an element is found in e.g. an elset,
        nset, typeEl, elType or other model attributes which is not defined in
        self.nodeCoords or self.elNodes respectively. Even if you catch that
        exception, the model is likely to be inconsistent. Better check
        beforehand.
        """

        (nodesOldToNew, elemsOldToNew) = Mesh.renumber(
            self, nodeStart, elemStart)

        # change self.nodeCoords and self.elNodes to correct types
        if nodeStart is not None:
            holdit = self.nodeCoords
            self.nodeCoords = container.NodesDict()
            dict.update(self.nodeCoords, holdit)

            # update self.nset
            try:
                newNset = [
                    (setName, set(nodesOldToNew[oldNode] for oldNode in nodes))
                    for setName, nodes in self.nset.iteritems()]
            except KeyError:
                notDefined = list(nodes.difference(nodesOldToNew))
                someNodesTxt = ",".join(map(str, notDefined[:10]))
                msgtext = (
                    "self.nset['%s'] contains %d nodes that are not defined in"
                    " the model, e.g. %s."
                    " This would leave the model in a corrupted state. Bailing"
                    " out." % (setName, len(notDefined), someNodesTxt))
                msg(msgtext)
                raise KeyError(msgtext)
            else:
                self.nset.clear()
                dict.update(self.nset, newNset)

        if elemStart is not None:
            holdit = self.elNodes
            self.elNodes = container.ElementNodesDict()
            dict.update(self.elNodes, holdit)

            # update elType
            try:
                newElType = dict(
                    (elemsOldToNew[oldElem], typ)
                    for oldElem, typ in self.elType.iteritems())
            except KeyError:
                notDefined = set(self.elType).difference(elemsOldToNew)
                someElemsTxt = ",".join(map(str, notDefined[:10]))
                msgtext = (
                    "self.elType contains %d elements that are not defined in"
                    " the model, e.g. %s."
                    " This would leave the model in a corrupted state. Bailing"
                    " out." % (len(notDefined), someElemsTxt))
                msg(msgtext)
                raise KeyError(msgtext)
            else:
                self.elType = newElType

            # update typeEl, this is treated differently...
            # ...because it's a defaultdict
            try:
                newTypeEl = [
                    (typ, set(elemsOldToNew[oldElem] for oldElem in elems))
                    for typ, elems in self.typeEl.iteritems()]
            except KeyError:
                notDefined = list(elems.difference(elemsOldToNew))
                someElemsTxt = ",".join(map(str, notDefined[:10]))
                msgtext = (
                    "self.typeEl['%s'] contains %d elements that are not"
                    " defined in the model, e.g. %s."
                    " This would leave the model in a corrupted state. Bailing"
                    " out." % (typ, len(notDefined), someElemsTxt))
                msg(msgtext)
                raise KeyError(msgtext)
            else:
                self.typeEl.clear()
                dict.update(self.typeEl, newTypeEl)

            # reinstall self.elShape and self.shapeEl
            self.elShape = container.FakeElShape(
                self.elType, self.elTypeToShapeSep)
            self.shapeEl = container.FakeShapeEl(
                self.typeEl, self.elShapeToTypes, self.elTypeToShapeSep)

            # update self.elset
            try:
                newElset = [
                    (setName, set(elemsOldToNew[oldElem] for oldElem in elems))
                    for setName, elems in self.elset.iteritems()]
            except KeyError:
                notDefined = list(elems.difference(elemsOldToNew))
                someElemsTxt = ",".join(map(str, notDefined[:10]))
                msgtext = (
                    "self.elset['%s'] contains %d elements that are not"
                    " defined in the model, e.g. %s."
                    " This would leave the model in a corrupted state. Bailing"
                    " out." % (setName, len(notDefined), someElemsTxt))
                msg(msgtext)
                raise KeyError(msgtext)
            else:
                self.elset.clear()
                dict.update(self.elset, newElset)

        # Not yet converted, not touched:
        # self.boundary = container.BoundaryDict()
        # self.nodeBC = container.NodeBoundaryDict(self)

        # self.surface = container.SurfacesDict()
        # self.mpc = container.MpcsList()

        # self.initCond = container.InitialCondContainer()

        return (nodesOldToNew, elemsOldToNew)

    #} end of edit / modify

    ########################################################################
    #{ transform geometry
    def _rotateXXX(self, refPoint, rotAmount, meshRotFunc, vecRotFunc):
        "used by rotate, rotateX, ..."
        # rotate nodes
        meshRotFunc(self, refPoint, rotAmount)
        # rotate orientations
        for ori in self.orientation.itervalues():
            ori.a = vecRotFunc(ori.a, rotAmount)
            ori.b = vecRotFunc(ori.b, rotAmount)
        # rotate beam cross section
        for prop in self.properties.itervalues():
            if prop['TYPE'] == "BEAMSECTION":
                sectiontype = prop.get('SECTION', "undefined")
                if sectiontype in ("BOX", "CIRC", "HEX", "I", "L", "PIPE",
                                   "RECT", "TRAPEZOID"):
                    beamAxis1 = prop['DATA'][1]
                elif sectiontype == "ARBITRARY":
                    beamAxis1 = prop['DATA'][-1]
                else:
                    msg("WARNING: BEAM SECTION, SECTION=%s is"
                        " not being rotated properly."
                        % (sectiontype))
                    continue
                beamAxis1[:] = vecRotFunc(beamAxis1, rotAmount)

    @staticmethod
    def _vector_rot_mat(vec, rotMat):
        return mat_multvec(rotMat,vec)

    def rotate(self, refPoint, rotMat):
        """Rotate the model according to the specified transformation matrix
        rotMat relative to the given point refPoint.

        @param refPoint: point around which to rotate the mesh
        @param rotMat: transformation matrix. No check is performed but the
           matrix should perform a pure rotation.

        @Note: Besides rotating the mesh (nodes) the following items will be
        treated:
         - orientations will be rotated (the two points a and b)
         - beam sections for section type BOX, CIRC, HEX, I, L, PIPE, RECT,
           TRAPEZOID, and ARBITRARY will be rotated. For other BEAM SECTION
           types a warning will be issued.

        @Note: Not being rotated:
         - initial conditions of type stress
        """
        self._rotateXXX(refPoint, rotMat, Mesh.rotate, Model._vector_rot_mat)

    def rotateX(self, refPoint, alpha):
        """Rotate the model around a line parallel to the x axis.

        @param refPoint: point around which to rotate the mesh
        @param alpha: angle in degree, pos == clockwise when looking in
        x-direction.

        @Note: See L{rotate} for further note on limitations and features.
        """
        self._rotateXXX(refPoint, alpha, Mesh.rotateX, vector_rot_x)

    def rotateY(self, refPoint, alpha):
        """Rotate the mesh around a line parallel to the y axis.

        @param refPoint: point around which to rotate the mesh
        @param alpha: angle in degree, pos == clockwise when looking in
        y-direction.

        @Note: See L{rotate} for further note on limitations and features.
        """
        self._rotateXXX(refPoint, alpha, Mesh.rotateY, vector_rot_y)

    def rotateZ(self, refPoint, alpha):
        """Rotate the mesh around a line parallel to the z axis.

        @param refPoint: point around which to rotate the mesh
        @param alpha: angle in degree, pos == anticlockwise when seen from
        above (i.e. looking in neg. z-direction).

        @Note: See L{rotate} for further note on limitations and features.
        """
        self._rotateXXX(refPoint, alpha, Mesh.rotateZ, vector_rot_z)
    #} end of transform geometry

    ########################################################################
    #{ other data extraction methods

    def getConnectedNodes(self, elems=None):
        """Returns a set of nodes that are connected to the elements given as
        argument.

        No check is performed if the nodes are actually defined in the model.
        If you want an error then check the elems-argument against self.elNodes
        beforehand like this:
         >>> elems = set(elementNumbers)
         >>> if elems.difference(model.elNodes):
         >>>     msg("Warning: Elements not defined in the model.")
         >>> else
         >>>     nodes = model.getConnectedNodes(elems)

        @param elems: A iterator of element numbers (also accepts a dict,
          in that case takes its keys.) Or anything that L{getUnionSet}
          accepts as input. Elements not defined are silently ignored.

          If elems==None, return all nodes connected to any element in the
          model.

        @Note: It's not checked if the nodes are actually defined in the model.
        """
        if not(elems is None or isinstance(elems, dict)):
            elems = self.getUnionSet("elset", elems)
        return Mesh.getConnectedNodes(self, elems)

    def getNodesFromElset(self, elements=None, raiseKeyError=True):
        """DEPRECATED: Use L{getConnectedNodes} instead. If you want an error
        then check the elements-argument against self.elNodes beforehand.

        Returns a set of nodes that are connected to the elset given as
        argument.

        elements may be an elset name or a list (or other iterable) of elset
        names and/or element numbers. (Even a mixed list of elset names and
        element numbers is accepted.)

        If elements==None, return all nodes connected to any
        element in the model.

        If raiseKeyError==True, raises a KeyError if an element in any of the
        given sets is not in self.elNodes. Otherwise just a warning is written
        to the logfile.
        """

        if elements is None:
            elset = self.elNodes.iterkeys()
        else:
            elset = self.getUnionSet("elset", elements)

        nodes = set()
        elemsNotFound = list()
        for element in elset:
            try:
                nodes.update(self.elNodes[element])
            except KeyError, err:
                if elements is None:
                    raise KeyError(
                        "Very strange error in Model.getNodesFromElset():"
                        " Element %d is not defined in the model (not in"
                        " self.elNodes). Try again or look into the code."
                        % err.args[0])
                else:
                    if raiseKeyError:
                        maxElemArgStrLen=60
                        elemArgStr = str(elements)
                        if len(elemArgStr)>maxElemArgStrLen:
                            elemArgStr = "%s..." % elemArgStr[:maxElemArgStrLen]
                        msgtext = (
                            "In Model.getNodesFromElset(): Element %d from the"
                            " given elements argument %s is not defined in the"
                            " model (not in self.elNodes)"
                            % (err.args[0], elemArgStr))
                        raise KeyError(msgtext)
                    else:
                        elemsNotFound.append(err.args[0])

        # if an error occurred:...
        if len(elemsNotFound) and not raiseKeyError:
            maxElemArgStrLen=60
            maxNotFoundLen=100

            elemArgStr = str(elements)
            if len(elemArgStr)>maxElemArgStrLen:
                elemArgStr = "%s..." % elemArgStr[:maxElemArgStrLen]

            elemNotFoundStr = str(elemsNotFound)
            if len(elemNotFoundStr)>maxNotFoundLen:
                elemNotFoundStr = "%s..." % elemNotFoundStr[:maxNotFoundLen]

            msg("WARNING: In Model.getNodesFromElset(): %d elements from the"
                " given elset %s are not defined in the model (not in"
                " self.elNodes):" % (len(elemsNotFound), elemArgStr))
            msg(elemNotFoundStr)

        return nodes

    def regexpNames(self, attribute, *args):
        r"""
        Returns a sorted list of all elset, nset or surface names that match
        one of the given the regular expressions.

        >>> m = Model().read("mymodel.inp")
        >>> allSeqElsetNames = m.regexpNames("elset", "DEV_", "STOPE_")
        >>> allBcNsetNames = m.regexpNames("nset", "NODES_..")

        You can use the *-expression if you already have a list of regular
        expressions:
         >>> setNames = ["DEV_", "STOPE_", "UCUT_"]
         >>> otherSetNames =  m.regexpNames("elset", *setNames)

        @param attribute: may be "elset", "nset", "surface"

        @param args: (list of all further arguments) state regular expression
          search string(s).

          dot (".") serves as wildcard for a single character, ".*" is a
          wildcard for any number of characters including zero, ".+" same
          but at least one character

          Some examples:
           - "DEV_" matches "DEV_STP", "dev_ug_1" but not "DEV01" and not
             "UDEV_1"
           - ".*M02" matches "DEV_M02", "M02", "M0240"
           - ".+M02" matches "DEV_M0240" but not "M0240"
           - ".+M02$" matches "DEV_M02" but not "DEV_M0240" and not "M02"
           - "..._UG" matches "DEV_UG1", "PRD_UG2" but not "PROD_UG3"
           - "" (empty string) matches all

          for further info refer to:
          http://docs.python.org/library/re.html#regular-expression-syntax

        @note: regexpNames() uses the re.match() method. It's case insensitive.
        """
        res = set()
        if not(hasattr(self, attribute)):
            raise AttributeError("Invalid attribute argument %s to regexpNames"
                                 % attribute)
        for regexp in args:
            rec = re.compile(regexp, re.IGNORECASE)
            for name in getattr(self, attribute):
                if rec.match(name):
                    res.add(name)
        res = sorted(res)
        return res

    def getUnionSet(self, settype, items, onlyValid=False):
        r"""
        Return a set that is the union of all items.

        @param settype: May be "nset" or "elset" (case matters)

        @param items: May be a set name or a list (or other iterable) of
          set names or element/node numbers. Set names and element/nodes
          numbers may be mixed.

          If items is not a string but iterable and any of its items
          is neither a string (i.e. elset name) nor an int (i.e. element
          number), an UnionSetInvalidItem exception is raised with the
          questionable item as argument.

        @param onlyValid: If True check the result against self.elNodes or
          self.nodeCoords and return only those items that are defined in the
          model. Note that in that case it returns an empty set if self only
          contains the sets at all.

        @return: A set that is the union of all items. The set returned is a
          new set (not a reference to any other elset) in any case.
          If no elements/nodes are found (unknown set name, empty set) an empty
          set is returned. No warning in this case!
        """

        try:
            attr = getattr(self, settype)
        except AttributeError:
            raise ValueError("Invalid settype argument %s to getUnionSet"
                             % settype)

        if isinstance(items, basestring):
            try:
                result = set(attr[items])
            except KeyError:
                result = set()
        else:
            result = set()
            for item in items:
                if isinstance(item, basestring):
                    try:
                        result.update(attr[item])
                    except KeyError:
                        pass
                elif isinstance(item, int):
                    result.add(item)
                elif isinstance(item, long):
                    result.add(int(item))
                else:
                    raise UnionSetInvalidItem(
                        "%s of type %s is an invalid item for the items"
                        " argument of %s.%s.getUnionSet()."
                        % (str(item), str(type(item)),
                           self.__class__.__module__, self.__class__.__name__))

        if onlyValid:
            if settype=="elset":
                result.intersection_update(self.elNodes)
            elif settype=="nset":
                result.intersection_update(self.nodeCoords)
        return result

    def getUnionElsetREs(self, *args):
        """returns a set of element numbers from the union of all elsets whose
        names match the given regular expressions

        You can use the *-expression if you already have a list of regular
        expressions:
         >>> setNames = ["DEV_", "STOPE_", "UCUT_"]
         >>> otherElset =  m.getUnionElsetREs(*setNames)

        shortcut for self.L{getUnionSet}("elset",self.L{regexpNames}("elset",
        *args))
        """
        assert all(isinstance(val, basestring) for val in args), \
            "All arguments to Model.getUnionElsetREs must be strings."
        return self.getUnionSet("elset",self.regexpNames("elset",*args))

    def getUnionElset(self, *args):
        """Deprecated alias for getUnionElsetREs() kept for compatibility
        reasons. Use getUnionElsetREs()!
        """
        return self.getUnionElsetREs(*args)

    def _checkRecognizedBlocks(self, recognizedBlocks):
        """for internal use only
        check the recognizedBlocks argument of the read and write methods
        or set all recognized blocks as default.
        """

        if isinstance(recognizedBlocks, basestring):
            recognizedBlocks = [b.strip() for b in recognizedBlocks.split(",")]
        elif (recognizedBlocks
              and not(isinstance(recognizedBlocks, (list, tuple, set))
                      or (hasattr(recognizedBlocks, "next")
                          and callable(recognizedBlocks.next)))):
            raise ValueError(
                "The recognizedBlocks argument to the abq_model_02.Model"
                ".checkRecognizedBlocks() method must be a list or another"
                " iterable object. It was found to be of type %s: %s"
                % (str(type(recognizedBlocks)), str(recognizedBlocks)))

        if recognizedBlocks:
            # convert to uppercase and remove all spaces
            recognizedBlocks = set([b.upper().replace(" ","")
                                    for b in recognizedBlocks])
            if not recognizedBlocks.issubset(self.recognizedBlocks):
                raise ReadError(
                    "ERROR: The following input file command(s) is/are not"
                    " implemented (yet) and can not be processed: %s"
                    % ", ".join(recognizedBlocks - set(self.recognizedBlocks)))
        else:
            recognizedBlocks = self.recognizedBlocks

        return recognizedBlocks

    def freeSurfElemsFromElset(self, elset):
        """DEPRECATED! Use L{bae.surface_02.ElemFaceSurface().addFreeSurfFromElset(...)<bae.surface_02.ElemFaceSurface.addFreeSurfFromElset>}

        Returns a dict of all elements on the perimeter surface / exterior
        surface of the given elset. Keys of dict are element faceIds.

        elset may be an element set name or a list of element set names.
        """

        elset = self.getUnionSet("elset", elset)

        faceToElem = defaultdict(list)
        for thisElem in elset:
            for faceNum, face in enumerate(self.elemFacesIter(thisElem)):
                faceToElem[face].append([thisElem,faceNum])

        freeSurfElems = defaultdict(set)
        for face, elems in faceToElem.iteritems():
            if len(elems)==1:
                freeSurfElems['S%d' %(elems[0][1]+1)].add(elems[0][0])
            else:
                pass

        return dict(freeSurfElems)

    def freeSurfNodesFromElset(self, elset=None, includeMidNodes=True):
        """Returns a set of all nodes on the perimeter surface of the given
        elset.

        If there are quadratic tets and linear cohesives in the model (and in
        elset) then also consider their non-conforming connection. I.e. don't
        make faults part of the exteriour surface. The parameter
        includeMidNodes is forcedly True in this case.

        @param elset: may be a set of element numbers or an element set name
           or a list of element set names (everything accepted by
           L{getUnionSet()}. If not given defaults to the whole model.
        @param includeMidNodes: If False then return corner nodes only,
           otherwise all nodes on the surface will be returned.
        """

        if elset is None:
            elset = self.elNodes
        elif not(isinstance(elset, set)):
            elset = self.getUnionSet("elset", elset)

        if ( ("COH3D6" in self.typeEl) or ("C3D10M" in self.typeEl)
             and self.typeEl["COH3D6"].intersection(elset)):

            #--- some preparation --------------------
            # dict edge (frozenset of corner nodes) -> corresponding mid node
            # nodes identified by index (0,1,2,....,9)
            # i.e.: { (0,1) : 5, (1,2) : 6, ... }
            edgeToMidNode = dict(
                (frozenset(edge), idx+4)
                for idx, edge in enumerate(Model.quadMidNodeEdges["TET"]))

            def getCornerFaceQ(p0, p1, p2):
                """Get part of face that is (potentially) connected to one
                cohesive face (ids of three nodes) on the tet face identified
                by p0,p1,p2. Get the face part that is connected to p0.
                """
                return [ p0,
                         edgeToMidNode[frozenset((p0, p1))],
                         edgeToMidNode[frozenset((p2, p0))] ]

            def getMiddleFaceQ(p0, p1, p2):
                """Get part of face that is (potentially) connected to one
                cohesive face (ids of three nodes) on the tet face identified
                by p0,p1,p2. Get the part that is connected to all three mid
                nodes.
                """
                return [ edgeToMidNode[frozenset((p0, p1))],
                         edgeToMidNode[frozenset((p1, p2))],
                         edgeToMidNode[frozenset((p2, p0))] ]

            # faceLfaceQ ... dict { tet face : face parts } stating all face
            # parts on each face of a tet. A face part is a list of node ids
            # that are (potentially) connected to one cohesive face (ids of
            # three nodes like [0,5,6], [5,6,7], ..)
            faceLfaceQ = dict()
            for id, faceL in enumerate(Model.elShapeFaceNodes["TET"]):
                faceLfaceQ["S%d" % (id+1)] = [
                    getCornerFaceQ(faceL[0], faceL[1], faceL[2]),
                    getCornerFaceQ(faceL[1], faceL[2], faceL[0]),
                    getCornerFaceQ(faceL[2], faceL[0], faceL[1]),
                    getMiddleFaceQ(faceL[0], faceL[1], faceL[2]) ]
            #--- end of preparation --------------------

            # all external surface of tets as dict face id -> tet elements
            msg("Finding exteriour faces for tets.")
            faceToElem = defaultdict(list)
            for thisElem in self.typeEl["C3D10M"].intersection(elset):
                for faceNum, face in enumerate(self.elemFacesIter(thisElem)):
                    faceToElem[face].append([thisElem,faceNum])

            faceIdElset = defaultdict(set)
            for face, elems in faceToElem.iteritems():
                if len(elems)==1:
                    faceIdElset['S%d' %(elems[0][1]+1)].add(elems[0][0])
                else:
                    pass
            # from now on behave like an ordinary dict
            faceIdElset.default_factory = None

            msg("Finding face parts for tet elements.")
            # all face parts (potentially connected to one cohesive face)
            # set of face parts (three node labels)
            tetFaceParts = set()
            ticker = MsgTicker("... %d elements done")
            cnt = 0
            for faceId in ["S1", "S2", "S3", "S4"]:
                fourFacePartIds = faceLfaceQ[faceId]
                for elem in faceIdElset[faceId]:
                    nodes = self.elNodes[elem]
                    for oneFacePartIds in fourFacePartIds:
                        facePart = frozenset(
                            nodes[idx] for idx in oneFacePartIds)
                        if facePart in tetFaceParts:
                            msg("ERROR: facePart %s on elem %d is double"
                                " connected but it's on an exteriour face."
                                % (facePart, elem))
                        tetFaceParts.add(facePart)
                    cnt += 1
                    ticker.msg(cnt)
            del ticker
            msg("Found %d exteriour face parts on tets." % len(tetFaceParts))

            msg("Collecting cohesive faces.")
            cohsFaces = set()
            ticker = MsgTicker("... %d elements done")
            cnt = 0
            for elem in self.typeEl["COH3D6"].intersection(elset):
                nodes = self.elNodes[elem]
                cohsFaces.add(frozenset(nodes[0:3]))
                cohsFaces.add(frozenset(nodes[3:6]))
                cnt += 1
                ticker.msg(cnt)
            del ticker
            msg("Found %d exteriour faces on cohs." % len(cohsFaces))

            msg("Identifying connected faces.")
            ticker = MsgTicker(
                "... tested %s/%d cohs faces" % ("%d", len(cohsFaces)))
            for face in list(cohsFaces):
                try:
                    tetFaceParts.remove(face)
                    cohsFaces.remove(face)
                except KeyError:
                    pass
                cnt += 1
                ticker.msg(cnt)
            del ticker
            msg("Not connected faces left: %d on cohesives and %d on tets"
                % (len(cohsFaces), len(tetFaceParts)))

            msg("Collecting nodes on not connected faces.")
            freeSurfNodes = set()
            for face in cohsFaces:
                freeSurfNodes.update(face)
            for face in tetFaceParts:
                freeSurfNodes.update(face)
            msg("Found %d exteriour nodes." % len(freeSurfNodes))

        else:  # if not(("COH3D6" in self.typeEl) or ("C3D10M" in self.typeEl))
            faceToElem = defaultdict(list)
            for thisElem in elset:
                for faceNum, face in enumerate(
                        self.elemFacesIter(thisElem, sepLQ=includeMidNodes)):
                    faceToElem[face].append([thisElem,faceNum])

            freeSurfNodes = set()
            for face, elems in faceToElem.iteritems():
                if len(elems)==1:
                    freeSurfNodes.update(face)

        return freeSurfNodes

    def getSurfElems(self, surf):
        """DEPRECATED! Use L{bae.surface_02.ElemFaceSurface().updateFromModel(...)<bae.surface_02.ElemFaceSurface.updateFromModel>}

        Returns a dict of elements that define the surface surf.

        @param surf: Surface name of type string or a dict of all elements
          defining the surface. Keys of dict are element faceIds.

        @note: This returns references of the element sets. If you modify this
          data you may also modify the elset the surface referes to.
        """

        if type(surf)==str:
            surfElemData=dict()
            try:
                surfElemSets=self.surface[surf]
                for faceId in surfElemSets:
                    surfElemData[faceId]=self.elset[surfElemSets[faceId]]
            except KeyError:
                raise KeyError(
                    "Surface %s is not defined in model. Input of a surface"
                    " string name requires the existance of a similarily named"
                    " surface in the model." % surf)
        elif type(surf)==dict:
            surfElemData=surf
        else:
            raise ValueError("Surface not of correct type.")

        return surfElemData

    def surfaceAdd(self, surf1, surf2):
        """DEPRECATED! Use L{bae.surface_02.ElemFaceSurface}

        Constructs a boolean add operator for Abaqus continuum element surfaces.

        @param surf1 surf2: Can be surface names of type string or a dict of
        all elements defining the surfaces. Keys of dict are element faceIds.

        @return: Dict of all elements defining the surface. Keys of dict are
        element faceIds.
        """

        surfs=[surf1, surf2]

        surfElemData = dict()
        for surfNum in range(len(surfs)):
            surfElemData[surfNum+1]=self.getSurfElems(surfs[surfNum])

        faceIds=set(surfElemData[1].keys()).union(set(surfElemData[2].keys()))

        booleanSurfElems=defaultdict(set)
        for faceId in faceIds:
            if faceId in surfElemData[1]:
                booleanSurfElems[faceId].update(surfElemData[1][faceId])
            if faceId in surfElemData[2]:
                booleanSurfElems[faceId].update(surfElemData[2][faceId])

        return booleanSurfElems

    def surfaceSubtract(self, surf1, surf2):
        """DEPRECATED! Use L{bae.surface_02.ElemFaceSurface}

        Constructs a boolean subtract operator for Abaqus continuum element
        surfaces. Returns a dict of all elements defining the boolean surface.
        Keys of dict are element faceIds.

        @param surf1 surf2: can be surface names of type string or a dict of
        all elements defining the surfaces. Keys of dict are element faceIds.
        surf2 is subtracted from surf1.

        @return: Dict of all elements defining the surface. Keys of dict are
        element faceIds.

        """

        surfs=[surf1, surf2]

        surfElemData=dict()
        for surfNum in range(len(surfs)):
            surfElemData[surfNum+1]=self.getSurfElems(surfs[surfNum])

        faceIds=set(surfElemData[1].keys()).union(set(surfElemData[2].keys()))

        booleanSurfElems=defaultdict(set)
        for faceId in faceIds:
            if faceId in surfElemData[1]:
                booleanSurfElems[faceId].update(surfElemData[1][faceId])
            if faceId in surfElemData[2]:
                booleanSurfElems[faceId].difference_update(
                    surfElemData[2][faceId])

        return booleanSurfElems

    def surfaceIntersect(self, surf1, surf2):
        """DEPRECATED! Use L{bae.surface_02.ElemFaceSurface}

        Constructs a boolean intersect operator for Abaqus continuum element
        surfaces. Returns a dict of all elements defining the boolean surface.
        Keys of dict are element faceIds.

        @param surf1 surf2: can be surface names of type string or a dict of
        all elements defining the surfaces. Keys of dict are element faceIds.

        @return: Dict of all elements defining the boolean surface. Keys of
        dict are element faceIds.
        """

        surfs=[surf1, surf2]

        surfElemData=dict()
        for surfNum in range(len(surfs)):
            surfElemData[surfNum+1]=self.getSurfElems(surfs[surfNum])

        faceIds=set(surfElemData[1].keys()).union(set(surfElemData[2].keys()))

        booleanSurfElems=defaultdict(set)
        for faceId in faceIds:
            if faceId in surfElemData[1]:
                booleanSurfElems[faceId].update(surfElemData[1][faceId])
            if faceId in surfElemData[2]:
                booleanSurfElems[faceId].intersection_update(
                    surfElemData[2][faceId])

        return booleanSurfElems
    #} end of other data extraction methods

# this is not working, because we cannot import the volume_02 module because
# this in turn is importing the abq_model_02 module.
#
#     def pointsToTetLinShapeFunc(self, pointCoords):
#         """The function returns a tuple (element number, [SF1, SF2, SF3, SF4])
#         for each point in pointList. SF1...SF4 are the values of the (linear
#         tetrahedron) shape function associated with the four nodes of the tet
#         at the given point. Each of the SF1...SF4 values lies between 0 and 1.
#         They can be used to linearily interpolate values at the corner nodes.

#          >>> points = [[0,0,0], [0,2,0], [2,0,2]]
#          >>> elemShapeFunc = list(model.pointsToTetLinShapeFunc(points))
#          >>> print elemShapeFunc
#          [(None, [0, 0, 0, 0]), (3, [1, 0, 0, 0]), (4, [0, 0.5, 0.5, 0])]

#         @Returns: An iterator that yields a tuple as described above. If a given
#         point is not within any of the tet elements of the model the element
#         number item of the returned tuple is None.
#         """

#         # find min, max for all point coordinates
#         box = BoundingBox()
#         box.update(pointCoords)
#         msg("Requested element numbers and linear shape function values for"
#             " %d points in the following box:\n%s" % (len(pointCoords), box))

#         # create a volume instance to perform the element lookup
#         allTets = set()
#         for typ in self.elShapeToTypes["TET"].intersection(self.typeEl):
#             allTets.update(self.typeEl[typ])
#         vol = VolumeFromMesh(self, elset=allTets)
#         vol.initializePointSearchInBox(box=box)

#         # findig elements and shape function values for each point
#         # elemShapeFunc = [vol.pointIsInside(point, returnLinShapeFunc=True)
#         #                  for point in pointCoords]
#         msgTemplate = ("Finding elements for points (%s/%d done)"
#                        % ("%d", len(pointCoords)))
#         ticker = MsgTicker(msgTemplate=msgTemplate)
#         for cnt, point in enumerate(pointCoords):
#             res = vol.pointIsInside(point, returnLinShapeFunc=True)
#             ticker.msg(cnt+1)
#             if res==False:
#                 res = (None, [0,0,0,0])
#             yield res
#         del ticker


#     def pointsToTetQuadShapeFunc(self, pointCoords):
#         """The function returns a tuple (element number, [SF1, SF2, ..., SF10])
#         for each point in pointList. SF1...SF10 are the values of the shape
#         function associated with the ten nodes of the tet at the given point.
#         Each of the SF1...SF10 values lies between 0 and 1 and their sum is 1.
#         They can be used to interpolate values at the nodes in a quadratic tet
#         element.

#          >>> points = [[0,0,0], [0,2,0], [2,0,2]]
#          >>> elemShapeFunc = list(model.pointsToTetQuadShapeFunc(points))

#         @Returns: An iterator that yields a tuple as described above. If a given
#         point is not within any of the tet elements of the model the element
#         number item of the returned tuple is None.
#         """

# # node  1.0, xi, eta, zeta, xi**2, eta**2, zeta**2, xi*eta, xi*zeta, eta*zeta
# # 1   [[ 1.  -3.  -3.  -3.    2.     2.       2.      4.      4.       4.]
# # 2    [ 0.  -1.   0.   0.    2.     0.       0.      0.      0.       0.]
# # 3    [ 0.   0.  -1.   0.    0.     2.       0.      0.      0.       0.]
# # 4    [ 0.   0.   0.  -1.    0.     0.       2.      0.      0.       0.]
# # 5    [ 0.   4.   0.   0.   -4.     0.       0.     -4.     -4.       0.]
# # 6    [ 0.   0.   0.   0.    0.     0.       0.      4.      0.       0.]
# # 7    [ 0.   0.   4.   0.    0.    -4.       0.     -4.      0.      -4.]
# # 8    [ 0.   0.   0.   4.    0.     0.      -4.      0.     -4.      -4.]
# # 9    [ 0.   0.   0.   0.    0.     0.       0.      0.      4.       0.]
# # 10   [ 0.   0.   0.   0.    0.     0.       0.      0.      0.       4.]]

#         # find min, max for all point coordinates
#         box = BoundingBox()
#         box.update(pointCoords)
#         msg("Requested element numbers and linear shape function values for"
#             " %d points in the following box:\n%s" % (len(pointCoords), box))

#         # create a volume instance to perform the element lookup
#         allTets = set()
#         for typ in self.elShapeToTypes["TET"].intersection(self.typeEl):
#             allTets.update(self.typeEl[typ])
#         vol = VolumeFromMesh(self, elset=allTets)
#         vol.initializePointSearchInBox(box=box)

#         # findig elements and shape function values for each point
#         # elemShapeFunc = [vol.pointIsInside(point, returnLinShapeFunc=True)
#         #                  for point in pointCoords]
#         msgTemplate = ("Finding elements for points (%s/%d done)"
#                        % ("%d", len(pointCoords)))
#         ticker = MsgTicker(msgTemplate=msgTemplate)
#         sf0coeff = [1., -3., -3., -3., 2., 2., 2., 4., 4., 4.]
#         for cnt, point in enumerate(pointCoords):
#             res = vol.pointIsInside(point, returnLinShapeFunc=True)
#             ticker.msg(cnt+1)
#             if res==False:
#                 res = (None, [0,]*10)
#             else:
#                 xxx = (1.0, res[1][1], res[1][2], res[1][3],
#                        res[1][1]**2, res[1][2]**2, res[1][3]**2,
#                        res[1][1]*res[1][2], res[1][1]*res[1][3],
#                        res[1][2]*res[1][3])
#                 res = (res[0], [
#                         sum((a*b for a,b in izip(sf0coeff, xxx))),
#                         2.0*xxx[4] - xxx[1],
#                         2.0*xxx[5] - xxx[2],
#                         2.0*xxx[6] - xxx[3],
#                         4.0*xxx[1] - 4.0*xxx[4] - 4.0*xxx[7] - 4.0*xxx[8],
#                         4.0*xxx[7],
#                         4.0*xxx[2] - 4.0*xxx[5] - 4.0*xxx[7] - 4.0*xxx[9],
#                         4.0*xxx[3] - 4.0*xxx[6] - 4.0*xxx[8] - 4.0*xxx[9],
#                         4.0*xxx[8],
#                         4.0*xxx[9],
#                         ])
#             yield res
#         del ticker


    ########################################################################
    #{ some service functions: iterator function

    def elemFacesIter(self, elNum, sepLQ=False):
        """Returns an iterator of the faces of the given element

        >>> myElem = 5731
        >>> for faces in elemFacesIter(myElem):
        >>>    ...

        faces will be:
          - frozenset((3, 43, 12)), then
          - frozenset((43, 3, 36)) etc.

        @param elNum: element label
        @param sepLQ: "use seperate Linear/Quadratic element shape?"
           This option only makes a difference for quadratic elements.
           If True then return frozensets of all nodes on the respective face.
           If False then return frozensets of only the linear, i.e. corner-
           nodes on each face.

        @returns: An iterator that yields all faces of the element. Each face
           is represented by a frozenset of the node labels on this face. For
           a linear tet element this would be all combinations of three out of
           the four corner nodes of this tet.
        """
        try:
            if sepLQ:
                elShape = self.elTypeToShapeSep[self.elType[elNum]]
            else:
                elShape = self.elTypeToShape[self.elType[elNum]]
        except KeyError:
            raise KeyError("elType %s of element %d not yet implemented in"
                           " elTypeToShape: %s."
                           % (self.elType[elNum], elNum,
                              self.elTypeToShape.keys()))
        try:
            thisNodeIdx = self.elShapeFaceNodes[elShape]
        except KeyError:
            raise KeyError("elShape %s not equivalent to one of the known:"
                           " %s." % (elShape, self.elShapeFaceNodes.keys()))

        nodesOfThisElem = self.elNodes[elNum]
        for i in thisNodeIdx:
            yield frozenset([nodesOfThisElem[j] for j in i])
    #} end of some service functions: iterator function


    ########################################################################
    #{ some structural data objects (static class members)

    #: element types available for each element shape
    elShapeToTypes = {
        'LINE': frozenset(('B31','T3D2', 'CONN3D2')),
        'TRI': frozenset(('CPE3','S3','S3R','DS3','M3D3',
                          'CPS6M','CPE6','STRI65','DS6')),
        'QUAD': frozenset(('CPS4R','CPE4','CPE4R','S4','S4R','S4RS',
                           'CPS8','CPS8R','CPE8','CPE8R','CPE8RP')),
        'TET': frozenset(('C3D4','C3D4H','C3D4P','C3D4T',
                          'C3D10','C3D10M','C3D10H','C3D10MH','C3D10MP')),
        'TETHT': frozenset(('DC3D4',)),
        'WEDGE': frozenset(('COH3D6','COH3D6P',
                            'C3D6','C3D6R','C3D6P','DC3D6','C3D15','SC6R')),
        'HEX': frozenset(('C3D8','C3D8P','C3D8R','C3D20','C3D20R','SC8R')),

        'LINE_L': frozenset(('B31','T3D2', 'CONN3D2')),
        'TRI_L': frozenset(('CPE3','S3','S3R','DS3','M3D3')),
        'QUAD_L': frozenset(('CPS4R','CPE4','CPE4R','S4','S4R','S4RS')),
        'TET_L': frozenset(('C3D4','C3D4H','C3D4P','C3D4T')),
        'TETHT_L': frozenset(('DC3D4',)),
        'WEDGE_L': frozenset(('COH3D6','COH3D6P',
                              'C3D6','C3D6R','C3D6P','DC3D6','SC6R')),
        'HEX_L': frozenset(('C3D8','C3D8P','C3D8R','SC8R')),

        'LINE_Q': set(),
        # For LINE_Q 'B32' would be a candidate but is problematic: Because
        # of the node ordering it does not fit into the getCentroid()-algorithm
        'TRI_Q': frozenset(('CPS6M','CPE6','STRI65','DS6')),
        'QUAD_Q': frozenset(('CPS8','CPS8R','CPE8','CPE8R', 'CPE8RP')),
        'TET_Q': frozenset(('C3D10','C3D10M','C3D10H','C3D10MH','C3D10MP')),
        'WEDGE_Q': frozenset(('C3D15',)),
        'HEX_Q': frozenset(('C3D20','C3D20P','C3D20R')),
        }

    #: dicts elType to shape like {'C3D4' : 'TET', ...}
    #: not distinguishing order of elements, yields 'TRI', 'TET', ...
    elTypeToShape = dict()
    #: seperate shapes for different order: 'TRI_L', 'TRI_Q', 'TET_L', ...
    elTypeToShapeSep = dict()
    for shape, elTypes in elShapeToTypes.iteritems():
        if shape.endswith('_L') or shape.endswith('_Q'):
            currDict = elTypeToShapeSep
        else:
            currDict = elTypeToShape
        for t in elTypes:
            currDict[t] = shape
    del shape, elTypes, t, currDict

    #: standard type conversion to quadratic, see L{Model.convertToQuadratic}()
    elTypeLinToQuad = {
        'C3D4': 'C3D10M',
        'C3D4P': 'C3D10MP',
        'S3': 'DS6',
        'COH3D6': 'COH3D6',
        'C3D6': 'C3D15',
        }

    # deprecated symbols (for backwards compatibility only):
    # eltypesTet, eltypesWedge, eltypesTri, eltypesBrick
    # don't use them any more, use elShapeToType['TRI'] ... instead
    eltypesTet = elShapeToTypes['TET']
    eltypesWedge = elShapeToTypes['WEDGE']
    eltypesTri = elShapeToTypes['TRI']
    eltypesBrick = elShapeToTypes['HEX']

    #: when converting arbitrary meshes to Abaqus models select an appropriate
    #: element type
    defaultMapShapeToType = {
        'LINE_L': 'B31',
        'TRI_L': 'S3',
        'QUAD_L': 'S4R',
        'TET_L': 'C3D4',
        'WEDGE_L': 'C3D6',
        'HEX_L': 'C3D8R',
        'LINE_Q': 'B32',
        'TRI_Q': 'DS6',
        'QUAD_Q': 'CPE8R',
        'TET_Q': 'C3D10M',
        'WEDGE_Q': 'C3D15',
        'HEX_Q': 'C3D20R',
        }
    #} end of some structural data objects (static class members)

    ########################################################################
    #{ other element specific methods

    def getCentroid(self, elementNumber):
        """calculate the elements centroid

        bugs: At the moment just takes the average of the position vectors of
        the corner nodes
        """

        try:
            ty = self.elType[elementNumber]
            nodes = self.elNodes[elementNumber]
        except KeyError:
            raise KeyError('Model.getCentroid(): element %d'
                           ' is not in the model.' % elementNumber)

        try:
            shape = self.elTypeToShape[ty]
            avg_n = self.elShapeToNbCornerNodes[shape]
        except KeyError:
            raise UnknownElType('Element type %s not implemented in %s.'
                                ' Cannot calculate centroid.'
                                % (ty, self.__module__))

        try:
            x=y=z=0.0
            for node in nodes[:avg_n]:
                coords=self.nodeCoords[node]
                x=x+coords[0]/avg_n
                y=y+coords[1]/avg_n
                z=z+coords[2]/avg_n
        except KeyError:
            raise KeyError('Model.getCentroid(elementNumber=%d): node %d'
                           ' is not in the model.' % (elementNumber, node))
        return [x,y,z]

    def _np_getElCentroidDict(self, elset=None):
        """Not used, use L{getElCentroidDict}() instead!

        This was to try if numpy makes it faster. It does not (yet). It's a
        drop in replacement for getElCentroidDict()

        The idea of running it elType by elType to avoid the type lookup for
        each element could also be introduced into the non numpy version.

        >>> from bae.abq_model_02 import Model
        >>> from time import clock
        >>>
        >>> m = Model().read("/work/RWAY_meshL5_onlyTets.inp")
        >>>
        >>> def getNew(els):
        >>>     cl = clock()
        >>>     centroids_new = m._np_getElCentroidDict(elset=els)
        >>>     print "new, centroids for %d elements: %g secs" % (len(els), clock()-cl)
        >>>     return centroids_new
        >>>
        >>> def getOld(els):
        >>>     cl = clock()
        >>>     centroids_old = m.getElCentroidDict(elset=els)
        >>>     print "old, centroids for %d elements: %g secs" % (len(els), clock()-cl)
        >>>     return centroids_old
        >>>
        >>> for N in (100, 1000, 10000, 100000, None):
        >>>     if N is None:
        >>>         someelems = m.elNodes.keys()
        >>>     else:
        >>>         someelems = list(m.elNodes)[:N]
        >>>
        >>>     getNew(someelems)
        >>>     getOld(someelems)
        """
        from time import clock
        cl = clock()
        if elset is not None:
            elset = self.getUnionSet("elset", elset, onlyValid=True)
        print "--1--: %g sec" % (clock()-cl)

        el_pts = dict()
        for ty, elems in self.typeEl.iteritems():
            if elset is not None:
                elems = elems.intersection(elset)
            if len(elems)==0:
                continue
            try:
                shape = self.elTypeToShape[ty]
                avg_n = self.elShapeToNbCornerNodes[shape]
            except KeyError:
                raise UnknownElType('Element type %s not implemented in %s.'
                                    ' Cannot calculate centroid.'
                                    % (ty, self.__module__))
            print "--2-- (type %s): %g sec" % (ty, clock()-cl)

            elems = list(elems)  # fix the element order
            print "--3-- (type %s): %g sec" % (ty, clock()-cl)
            # coords shape: (nb elements, avg_n, 3)
            coords = np.array(
                [[self.nodeCoords[node]
                  for node in self.elNodes[el][:avg_n]]
                 for el in elems],
                dtype=float)
            print "--3b-- (type %s): %g sec" % (ty, clock()-cl)
            # coords shape: (nb elements, avg_n, 3)
            coords = np.empty((len(elems), avg_n, 3), dtype=float)
            for i, el in enumerate(elems):
                for j, node in enumerate(self.elNodes[el][:avg_n]):
                    coords[i,j] = self.nodeCoords[node]
            print "--4-- (type %s): %g sec" % (ty, clock()-cl)
            centroids = coords.sum(axis=1) / avg_n
            print "--5-- (type %s): %g sec" % (ty, clock()-cl)
            # store results in el_pts
            el_pts.update(izip(elems, centroids))
            print "--6-- (type %s): %g sec" % (ty, clock()-cl)
        return el_pts

    def getElCentroidDict(self, elset=None):
        """Deprecated, use self.L{getElCentroids}() instead! This is just there
        for compatibility reasons.

        Calculate element centroids for all specified elements and return a
        dict {element: centroid}

        @param elset: Anything that self.L{getUnionSet}() accepts as input or
        None for all elements.
        """
        if elset is None:
            # caution, this is fast but dangerous if you modify the code
            # in this case here, elset is a dictionary, not a set
            elset = self.elNodes
        else:
            elset = self.getUnionSet("elset", elset, onlyValid=True)

        # create list with tuples (el numer, [points to check])
        ticker = MsgTicker("... calculated %s/%d element"
                           " centroids so far." % ("%d", len(elset)))
        el_pts = dict()
        for cnt, el in enumerate(elset):
            el_pts[el] = self.getCentroid(el)
            ticker.msg(cnt+1)
        del ticker
        return el_pts

    def getElCentroids(self, elems=None):
        """
        Calculate element centroids for all specified elements and return a
        dict {element: centroid}

        @param elems: Elements to be found the centroid for. Anything that
        self.L{getUnionSet}() accepts as input or None for all elements.
        Elements not in the model are silently ignored.
        """
        if elems is not None:
            elems = self.getUnionSet("elset", elems, onlyValid=True)
        return Mesh.getElCentroids(self, elems)


    def getQuadGaussPointsDict(self, elset=None):
        """
        Calculate the APPROXIMATE gauss points for the quadratic version of all
        specified elements and return a dict {element: [list of gauss points]}
        (That means: tets will always have 4 gausspoints. Those of the linear
        type C3D4 as well.)

        See also L{bae.mesh_01.Mesh.getElGaussPoints}.

        @param elset: Anything that self.L{getUnionSet}() accepts as input or
        None for all elements.

        @note: Only tested for tet elements so far.

        @Note: Just takes the average of the position vectors of the corner
        nodes as centroid and assumes half the way to each corner node as gauss
        point. THIS IS ONLY APPROXIMATE!

        For quadratic tet elements instead of 0.5*(centre+corner) taken here
        you would get the right value by (1-a)*centre+a*corner; a=sqrt(5)/5.
        But this would probably not be correct for bricks (and others)...
        """
        if elset is None:
            # caution, this is fast but dangerous if you modify the code
            # in this case here, elset is a dictionary, not a set
            elset = set(self.elNodes)
        else:
            elset = self.getUnionSet("elset", elset, onlyValid=True)

        # create list with tuples (el numer, [points to check])
        el_pts = dict()
        ticker = MsgTicker("... calculated %s/%d element"
                           " gauss points so far." % ("%d", len(elset)))
        cnt = 0
        elNodes = self.elNodes
        nodeCoords = self.nodeCoords
        for ty, typElems in self.typeEl.iteritems():
            if len(elset)==0:
                break

            elems = elset.intersection(typElems)
            if len(elems)==0:
                continue

            try:
                shape = self.elTypeToShape[ty]
                avg_n = self.elShapeToNbCornerNodes[shape]
            except KeyError:
                raise UnknownElType('Element type %s not implemented in %s.'
                                    ' Cannot calculate gauss points.'
                                    % (ty, self.__module__))
            for el in elems:
                cnt += 1
                nodes = elNodes[el]
                try:
                    points = [nodeCoords[node] for node in nodes[:avg_n]]
                except KeyError:
                    raise KeyError('Model.getCentroid(elementNumber=%d): node'
                                   ' %d is not in the model.' % (el, node))
                centroid = [float(sum(coords))/avg_n
                            for coords in izip(*points)]
                points = [[0.5*(xp+xc) for xp, xc in izip(pt, centroid)]
                          for pt in points]
                el_pts[el] = points
                ticker.msg(cnt)
            elset.difference_update(elems)

        del ticker
        return el_pts


    def getCentroidVolume(self, elementNumber):
        """calculate the elements centroid and volume

        >>> elem = 1623534
        >>> centroid, volume = model.getCentroidVolume(elem)

        bugs: At the moment works only for tets

        @Note: see also L{bae.mesh_01.Mesh.getElemVolField}
        """
        ty = self.elType[elementNumber]
        try:
            shape = self.elTypeToShape[ty]
            avg_n = 4
        except KeyError:
            raise UnknownElType('Element type %s not implemented in %s.'
                                ' Cannot calculate centroid.'
                                % (ty, self.__module__))
        if shape!='TET':
            raise UnknownElType(
                'Cannot calculate volume of element shape %s yet.' % shape)

        points = [self.nodeCoords[node]
                  for node in self.elNodes[elementNumber][:avg_n]]

        # centroid
        x=y=z=0.0
        for coords in points:
            x=x+coords[0]/avg_n
            y=y+coords[1]/avg_n
            z=z+coords[2]/avg_n

        # volume
        vec = [vector_minus(pt, points[0])
               for pt in points[1:]]
        vol = ((
              vec[0][0]*vec[1][1]*vec[2][2]
            + vec[0][1]*vec[1][2]*vec[2][0]
            + vec[0][2]*vec[1][0]*vec[2][1]
            - vec[0][2]*vec[1][1]*vec[2][0]
            - vec[0][1]*vec[1][0]*vec[2][2]
            - vec[0][0]*vec[1][2]*vec[2][1])
               / 6.0)

        return [x,y,z], vol

    def getCentroidTetElset(self, elset=None):
        """
        calculate the centre of the volume of the given elset

        @param elset: Must be a valid input to self.getUnionSet(),
          i.e. a list of element numbers, elset names,
          a single elset name,...
          If elset is None (the default) return the centre of the
          whole model.

        @returns: the centroid as a list of three coordinates.

        @note: Ignores all but tet elements
        """

        # make elset a list of all elements to be converged
        if elset is None:
            elset = self.elNodes
        else:
            elset = self.getUnionSet("elset", elset)

        assert len(elset)>0, "Can't calculate the centroid of an empty elset."

        allVol = 0.0
        allCentre = [0.0, 0.0, 0.0]
        for el in elset:
            try:
                ty = self.elType[el]
                shape = self.elTypeToShape[ty]
                if shape!='TET':
                    continue
            except KeyError:
                continue

            centre, vol = self.getCentroidVolume(el)
            allVol += vol
            vector_modif_add(allCentre,
                             vector_scale(centre, vol))

        if allVol==0.0:
            raise Exception(
                "Model.getCentroidTetElset(): given elset has no volume.\n"
                "Either there are no tet elements in this elset or none of"
                " of the elements is defined in the mesh.")
        allCentre = vector_scale(allCentre, 1.0/allVol)
        return allCentre


    def getEquivTets(self, elementNumber):
        """Return tet elements that fill the volume of this element.

        If the given element is already a tet or not a 3D element, return None

        Returns a list of nodelists for C3D4 linear tet elements:
        [[nodeNums 1. elem], [nodeNums 2. elem], ...]

        The tuple returned is suitable to insert the elements into the model
        quite easily:

        >>> result = model.getEquivTets(elementNumber)
        >>> nextElemNb = max(model.elNodes.iterkeys())+1
        >>> if result:
        >>>    model.insertElems(result, "C3D4", nextElemNb)
        >>>    model.removeElem(elementNumber)

        @bug: so far works for linear brick and wedge elements only.

        @todo: generalize with the class members elTypeToShape and so on. Needs
          a new class member dict for the translation of other elements instead
          of the if - elif construct.
        """
        this_type = self.elType[elementNumber]

        if ((this_type in self.eltypesTet)
                or (this_type in self.eltypesTri)):
            return None

        if this_type in self.eltypesBrick:
            el_nodes_idx = [
                (1, 3, 8, 6),
                (1, 3, 6, 2),
                (1, 8, 3, 4),
                (1, 6, 8, 5),
                (3, 8, 6, 7)]
        elif this_type in self.eltypesWedge:
            el_nodes_idx = [
                (1,2,3,4),
                (2,3,4,5),
                (3,4,5,6)]
        else:
            raise UnknownElType('Model.getEquivTets: Element type %s not'
                                ' implemented yet.' % this_type)

        oldnodes = self.elNodes[elementNumber]
        result = list()
        for newtet in el_nodes_idx:
            try:
                nodelist = [oldnodes[n-1] for n in newtet]
            except IndexError:
                raise IndexError(
                    "In abq_model_02.Model.getEquivTets(): Element %d of type"
                    " %s has only %d nodes."
                    % (elementNumber, this_type, len(oldnodes)))
            result.append(nodelist)

        return result

    def convertToTets(self, elset=None,
                      firstNewElemNb=None, updateElsets=True,
                      tetType="C3D4"):
        """Convert all non-tet volume elements to tets.

        All non-tet elements are filled with tets independently of each other.
        New diagonal edges cutting four-node-faces do not necessarily conform
        to those of the neighbours.

        Not suitable if you want to use the mesh for calculations or conversion
        to quadratic. Use self.convertToConformingTets() in that case.

        @param elset: If argument elset is specified must be a valid input to
           self.getUnionSet() i.e. a list of element numbers, elset names,
           a single elset name,...

           If None (the default) is specified: Convert all elements.

        @param firstNewElemNb: If not None: States the number of the first
           element that is being inserted. The others get successive numbers.
           If None just takes the last number + 1.

        @param updateElsets: If True (the default): Add the newly created
           elements to all sets the old (non tet) elements belonged to and
           remove the old from there.

        @param tetType: The type of the new elements to insert, default is
           C3D4.

        @returns: dictionary {old elNum : list of new elNums}
        """

        # nothing to convert if elNodes is empty
        # (we later expect it to contain something!)
        if len(self.elNodes)==0:
            return dict()

        # make elset a list of all elements to be converged
        if elset is not None:
            elset = self.getUnionSet("elset", elset)

        if firstNewElemNb is None:
            nextElemId = max(self.elNodes) + 1
        else:
            nextElemId = firstNewElemNb

        typesToConvert = [
            typ for typ in self.typeEl
            if ((typ.startswith('C3D') or typ.startswith('COH3D'))
                and typ not in self.elShapeToTypes['TET'])]

        oldToNew = dict()  # {old elNum : list of new elNums}
        newElNodes = dict()  # {elNum: list of nodeNums}
        for typ in typesToConvert:
            nb_notConv = 0
            oldElems = self.typeEl[typ]  # only find elements if necessary
            if elset is not None:
                # make a copy, don't modify self.typeEl[typ]!
                oldElems = oldElems.intersection(elset)
            for oldEl in oldElems:
                thisNodelist = self.getEquivTets(oldEl)
                if len(thisNodelist):
                    newElNums = range(
                        nextElemId, nextElemId+len(thisNodelist))
                    nextElemId += len(thisNodelist)
                    oldToNew[oldEl] = newElNums
                    newElNodes.update(izip(newElNums, thisNodelist))
                else:
                    nb_notConv += 1

            if nb_notConv>0:
                msg("WARNING: %d elements of type %s not converted to"
                    " tets abq_model_02.Model.convertToTets()."
                    % (nb_notConv, typ))

        # insert new elements
        if len(newElNodes):
            self.updateElem(None, tetType, newElNodes)

        # update elsets: add new and remove old elements
        if updateElsets:
            for els in self.elset.itervalues():
                oldInEls = els.intersection(oldToNew)
                for oldEl in oldInEls:
                    els.update(oldToNew[oldEl])
                els.difference_update(oldInEls)

        # remove old elements
        self.removeElem(oldToNew, updateElsets=False)

        # fini
        return oldToNew


    class _VertEdgeData(object):
        __slots__ = "nodes", "connWedges", "topIndex"
        def __init__(self):
            self.connWedges = set()
            self.topIndex = 0
    class _OrientedEdgesList(set):
        def extend(self, newEdges, orientation):
            for nextNodes in newEdges:
                orderedEdge = (nextNodes[orientation], nextNodes[1-orientation])
                if orderedEdge[::-1] in self:
                    raise Exception("Ambigous edge orientation")
                self.add(orderedEdge)

    def convertToConformingTets(self, elset=None,
                                firstNewElemNb=None, updateElsets=True,
                                tetType="C3D4", extraConnections=[]):
        """Convert all non-tet volume elements to tets.

        New diagonal edges cutting four-node-faces conform to those of the
        neighbours (other that what self.convertToTets() does).

        If you don't want to use the mesh for calculations or conversion to
        quadratic self.convertToTets() may be faster.

        @note: Converts only six node wedge elements so far.

        @param elset: If argument elset is specified must be a valid input to
               self.getUnionSet() i.e. a list of element numbers, elset names,
               a single elset name,...

               If None (the default) is specified: Convert all elements.

        @param firstNewElemNb: If not None: States the number of the first
        element that is being inserted. The others get successive numbers. If
        None just takes the last number + 1.

        @param updateElsets: If True (the default): Add the newly created
        elements to all sets the old (non tet) elements belonged to and remove
        the old from there.

        @param tetType: The type of the new elements to insert, default is C3D4.
        @param extraConnections: A list of tuples of node numbers. The nodes in
        each of the tuples are considered to be connected as if they were one,
        e.g. by an MPC.
        """

        # make elset a list of all elements to be converged
        if elset is not None:
            elset = self.getUnionSet("elset", elset)

        if firstNewElemNb is None:
            nextElemId = max(self.elNodes.iterkeys()) + 1
        else:
            nextElemId = firstNewElemNb

        typesToConvert = set([
            typ for typ in self.typeEl
            if ((typ.startswith('C3D') or typ.startswith('COH3D'))
                and typ not in self.elShapeToTypes['TET'])])

        ##--- deal with all wedge shaped elements
        linWedgeTypesInModel = typesToConvert.intersection(
            self.elShapeToTypes['WEDGE_L'])
        typesToConvert -= linWedgeTypesInModel

        allWedgeElems = set()
        for typ in linWedgeTypesInModel:
            if elset is None:
                allWedgeElems.update(self.typeEl[typ])
            else:
                allWedgeElems.update(elset.intersection(
                    self.typeEl[typ]))

        # vertEdgeToData :
        # {vertical wedge edge: _VertEdgeData struct}
        # note: vertical edges are those edges between two tri faces, they are
        # connecting node 1-4, 2-5 and 3-6 in an abaqus wedge element
        vertEdgeToData = defaultdict(Model._VertEdgeData)
        for wedgeElem in allWedgeElems:
            elNodes = self.elNodes[wedgeElem]
            for i1 in range(3):
                i2 = i1+3
                edge = frozenset((elNodes[i1], elNodes[i2]))
                vertEdgeToData[edge].connWedges.add(wedgeElem)

        # nodesToVertEdges: {nodes in vertical edges : set of vertical edges}
        # each value of the dict (set of edges) contains at least one edge
        nodesToVertEdges = defaultdict(set)
        for edge in vertEdgeToData:
            for node in edge:
                nodesToVertEdges[node].add(edge)

        # xConNodesToOtherNodeOnVertEdges: { node : list of nodes }
        xConNodesToOtherNodeOnVertEdges = defaultdict(list)
        for xConNodes in extraConnections:
            for node in xConNodes:
                for edge in nodesToVertEdges[node]:
                    othernode = iter(edge.difference((node,))).next()
                    xConNodesToOtherNodeOnVertEdges[node].append(othernode)

        # xConEdges:
        # { edge : set of ordered edges connected via extra connections }
        extraConnections = set([frozenset(xConNodes)
                                for xConNodes in extraConnections])
        xConEdges = dict()
        for xConNodes in extraConnections:
            xConNodesOrdered = tuple(xConNodes)
            allOtherNodes = [xConNodesToOtherNodeOnVertEdges[node]
                             for node in xConNodesOrdered]
            otherNodes = [0,0]
            for otherNodes[0] in allOtherNodes[0]:
                for otherNodes[1] in allOtherNodes[1]:
                    if (otherNodes[0]==otherNodes[1]
                            or frozenset(otherNodes) in extraConnections):
                        xConEdgesPair = [((xConNodesOrdered[0],otherNodes[0])),
                                         ((xConNodesOrdered[1],otherNodes[1]))]
                        for edge_i in range(2):
                            found = False
                            for flip_i in range(2):
                                orderedEdge = (xConEdgesPair[edge_i][flip_i],
                                               xConEdgesPair[edge_i][1-flip_i])
                                found = orderedEdge in xConEdges
                                if found:
                                    break
                            if found:
                                xConEdges[orderedEdge].add((
                                    xConEdgesPair[1-edge_i][flip_i],
                                    xConEdgesPair[1-edge_i][1-flip_i]))
                            else:
                                xConEdges[orderedEdge[::-1]]=set(xConEdgesPair)

        # xConEdges: transform key to frozenset
        xConEdges = dict([ (frozenset(key), value)
                           for key, value in xConEdges.iteritems()])

        # wedgeOrient: {wedge elem nb : 0/1 flag marks wedge to be flipped}
        wedgeOrient = dict()
        wedgesWOOrient = set(allWedgeElems)  # wedges w/o orientation
        # make a copy. Items will be consumed upon use
        notUsedXConEdges = dict(xConEdges)
        while 1:
            # initialize first wedge (no other oriented wedge attached so far)
            try:
                wedgeElem = wedgesWOOrient.pop()
                elNodes = self.elNodes[wedgeElem]
            except KeyError:
                break

            wedgeOrient[wedgeElem] = 0

            # edges with known orientation not yet processed in next loop
            orientedEdges = Model._OrientedEdgesList([
                (elNodes[i1], elNodes[i1+3]) for i1 in range(3)])
            while len(orientedEdges)>0:
                nodes = orientedEdges.pop()
                edge = frozenset(nodes)

                # append edges in notUsedXConEdges to orientedEdges
                try:
                    xConOrientedEdges = notUsedXConEdges[edge]
                except KeyError:
                    pass
                else:
                    if nodes in xConOrientedEdges:
                        newEdges = xConOrientedEdges.difference((nodes,))
                        flip = 0
                    elif nodes[::-1] in xConOrientedEdges:
                        newEdges = xConOrientedEdges.difference((nodes[::-1],))
                        flip = 1
                    else:
                        raise Exception(
                            "Program error: nodes %s not in xConOrientedEdges"
                            " %s." % (str(nodes), str(xConOrientedEdges)))
                    orientedEdges.extend(newEdges=newEdges, orientation=flip)
                    for edgeToRemove in xConOrientedEdges:
                        del notUsedXConEdges[frozenset(edgeToRemove)]

                # for all wedges on this edge:
                # determine wedge orientation and add new edges to orientedEdges
                edgeData = vertEdgeToData[edge]
                edgeData.nodes = (nodes)
                edgeNeighbourWedges = set((edgeData.connWedges))
                edgeNeighbourWedges.intersection_update(wedgesWOOrient)
                for nextWedge in edgeNeighbourWedges:
                    nextElNodes = self.elNodes[nextWedge]
                    nextOrient = None
                    newEdges = list()
                    for ii in range(3):
                        nextNodes = (nextElNodes[ii], nextElNodes[ii+3])
                        nextEdge = frozenset(nextNodes)
                        if nextEdge==edge:
                            if nextNodes == nodes:
                                nextOrient = 0
                            else:
                                nextOrient = 1
                        else:
                            newEdges.append(nextNodes)
                    assert nextOrient is not None
                    wedgeOrient[nextWedge] = nextOrient
                    wedgesWOOrient.remove(nextWedge)
                    orientedEdges.extend(newEdges=newEdges,
                                         orientation=nextOrient)

                ## end while len(orientedEdges)>0 ...
                del edge
        ## end while 1, wedgesWOOrient now empty

        # create node connectivity for new tets
        oldToNew = defaultdict(list)  # {old elNum : list of new elNums}
        newElNodes = dict()           # {elNum: list of nodeNums}
        vertEdgesToProcess = set(vertEdgeToData)
        while len(vertEdgesToProcess)>0:
            # take any edge not processed so far
            startEdge = vertEdgesToProcess.pop()
            centreEdgeAndXConns = [startEdge]
            # add edges connected via extraConnections
            try:
                xConOrientedEdges = xConEdges[startEdge]
            except KeyError:
                pass
            else:
                currentXConEdges = [frozenset(edge)
                                    for edge in xConOrientedEdges
                                    if frozenset(edge) != startEdge]
                # so this is the replacement that works (for the moment?)
                currentXConEdges = [
                    edge for edge in map(frozenset, xConOrientedEdges)
                    if edge != startEdge]
                centreEdgeAndXConns.extend(currentXConEdges)
                vertEdgesToProcess.difference_update(currentXConEdges)
                del xConOrientedEdges, currentXConEdges, edge  # local variables

            # now process those edges
            for centreEdge in centreEdgeAndXConns:
                edgeData = vertEdgeToData[centreEdge]
                for wedgeElem in edgeData.connWedges:
                    wedgeNodes = self.elNodes[wedgeElem]
                    tetnodes = list()
                    # base triangle for the new tet
                    for ii in range(3):
                        otherEdge = frozenset((wedgeNodes[ii],wedgeNodes[ii+3]))
                        otherData = vertEdgeToData[otherEdge]
                        tetnodes.append(otherData.nodes[otherData.topIndex])
                    if wedgeOrient[wedgeElem] == 1:
                        tetnodes = tetnodes[::-1]
                    # top corner of the new tet
                    tetnodes.append(edgeData.nodes[1])
                    newElNodes[nextElemId] = tetnodes
                    oldToNew[wedgeElem].append(nextElemId)
                    nextElemId += 1
                edgeData.topIndex = 1  # from now take the top node of this edge

        ##--- deal with all other elements (not treated)
        for typ in typesToConvert:
            nb_notConv = len(self.typeEl[typ])
            if nb_notConv>0:
                msg("WARNING: %d elements of type %s not converted to"
                    " tets abq_model_02.Model.convertToTets()."
                    % (nb_notConv, typ))

        # insert new elements
        if len(newElNodes):
            self.updateElem(None, tetType, newElNodes)

        # update elsets: add new and remove old elements
        if updateElsets:
            for els in self.elset.itervalues():
                oldInEls = els.intersection(oldToNew)
                for oldEl in oldInEls:
                    els.update(oldToNew[oldEl])
                els.difference_update(oldInEls)

        # remove old elements
        self.removeElem(oldToNew, updateElsets=False)

        # fini
        return oldToNew


    def convertToQuadratic(self,
                           firstNewNodeNb=None, firstNewElemNb=None,
                           updateNsets=False, updateElsets=True):
        """
        convert all elements to quadratic elements, insert all needed
        nodes

        converts only element types listed in self.elTypeLinToQuad

        firstNewNodeNb, firstNewElemNb (if not None) states the number of the
        first node/element that is being inserted to convert the elements to
        quadratic. The others get successive numbers. If None just takes the
        last number + 1.

        updateNsets ... add new nodes to the same nsets both other nodes on the
        edge are in. Not implemented yet.
        updateElsets ... (meaningful only for specially treated elements like
        cohesives) add new elements all elsets the original elements belonged
        to.

        In order to achieve different conversions you may overwrite the
        conversion dictionaries. They are initialized as class instances, i.e.
        Model.elTypeLinToQuad but referenced as self.elTypeLinToQuad, e.g.:

        >>> m = Model().read(...)
        >>> m.elTypeLinToQuad = dict(Model.elTypeLinToQuad)  # make a copy!
        >>> del m.elTypeLinToQuad['C3D6']  # don't convert C3D6 elements
        >>> m.convertToQuadratic()

        conversion dictionaries used:
        (Those are class members, better make a copy before modifying them)

         - elTypeLinToQuad ... a dict {old (linear) type: new (quadratic) type}

         - elTypeToShapeSep ... a dict elType to shape distinguishing the order
           of elements: {'C3D4':'TET_L', 'C3D10M':'TET_Q', 'S3','TRI_L', ...}
           order indicator: '_L' stands for linear, '_Q' for quadratic

         - quadMidNodeEdges ... mid nodes for quadratic elems:
            {shape:edges list} shape is one of 'TRI', 'TET', ... (omitting the
            order indicator) The first item of the edges list is a list of node
            indices identifying the corner nodes between which the first
            quadratic mid node lies, the second item represents the edge for
            the second mid node and so on.

        @warning: updateNsets option does not work

        @returns: edgeToNewMidNode, oldElemToReplacements
            edgeToNewMidNode {edge:midNode}, edge is a frozenset of two nodes
            labels, midNode is the node label of the newly created mid node
            oldElemToReplacements is a dict containing only cohesive elements
            {old lin elem label: list of new elem labels}
        """

        if firstNewNodeNb is None:
            nextNodeId = max(self.nodeCoords) + 1
        else:
            nextNodeId = firstNewNodeNb

        if firstNewElemNb is None:
            nextElemId = max(self.elNodes) + 1
        else:
            nextElemId = firstNewElemNb

        # for cohesives: {old elNum : new elNums replacing old}
        oldElemToReplacements = defaultdict(list)
        edgeToNewMidNode = dict()  # {edge:midNodeID}
        for thisLinType in list(self.typeEl):  # make a copy, will be modified
            linTypesElements = self.typeEl[thisLinType]
            lenLinTypesElements = len(linTypesElements)

            try:
                thisQuadType = self.elTypeLinToQuad[thisLinType]
            except KeyError:
                msg('No matching QuadType for %d elements of LinType %s.'
                    % (lenLinTypesElements, thisLinType))
                continue
            quadTypesElements = self.typeEl[thisQuadType]

            elShape = self.elTypeToShapeSep[thisLinType]
            if not(elShape.endswith('_L')):
                msg('%d elements of type %s are not considered linear and are'
                    ' not converted therefore.'
                    % (lenLinTypesElements, thisLinType))
                continue
            try:
                # list of edges ordered by the corresponding mid node
                # represented by node index
                edgeIndexList = self.quadMidNodeEdges[elShape[:-2]]
                nbCornerNodes = self.elShapeToNbCornerNodes[elShape[:-2]]
            except KeyError:
                msg('Conversion data missing for %d elements of LinType %s.'
                    ' Those elements are not being converted.'
                    % (lenLinTypesElements, thisLinType))
                continue

            if thisQuadType=='COH3D6':
                # special treatment for cohesive elements
                # insert new nodes only along the edges of the top and bottom
                edgeIndexList = edgeIndexList[:6]
                # nodes indices for the four new wedge element replacing an old
                nodesIdx4CohWedges = [[0,6,8,3,9,11], [6,1,7,9,4,10],
                                      [7,2,8,10,5,11], [6,7,8,9,10,11]]
                newElNodes = dict()

            msg('%d elements of LinType %s will be converted to QuadType %s.'
                ' Calculating new connectivity now...'
                % (lenLinTypesElements, thisLinType, thisQuadType))
            msgTemplate = ("... processed %d of "+str(lenLinTypesElements)
                           +" elements.")
            ticker = MsgTicker(msgTemplate)
            # copy the set linTypesElements, because it will be modified
            for iCount, element in enumerate(list(linTypesElements)):
                thisElNodes = self.elNodes[element]
                newMidNodes = list()
                for thisEdgeIdx in edgeIndexList:
                    thisEdge = frozenset(thisElNodes[i] for i in thisEdgeIdx)

                    # find or create new mid node
                    try:
                        centreNode = edgeToNewMidNode[thisEdge]
                    except KeyError:
                        centreCoords = [0.0, 0.0, 0.0]
                        for node in thisEdge:
                            vector_modif_add(centreCoords,self.nodeCoords[node])
                        vector_modif_scale(centreCoords, 0.5)

                        centreNode = nextNodeId
                        nextNodeId += 1

                        self.nodeCoords[centreNode] = centreCoords
                        edgeToNewMidNode[thisEdge] = centreNode
                    # store new mid node
                    newMidNodes.append(centreNode)

                # convert the actual element to quadratic: insert new nodes
                # even if this type will be treated specially: store the nodes
                thisElNodes[nbCornerNodes:] = newMidNodes

                if thisQuadType=='COH3D6':
                    # special treatment for cohesive elements
                    # insert 4 wedges in place of the one old
                    newElemsList = list()
                    for wedgeNodesIdx in nodesIdx4CohWedges:
                        wedgeNodes = [thisElNodes[i] for i in wedgeNodesIdx]
                        newElNodes[nextElemId] = wedgeNodes
                        newElemsList.append(nextElemId)
                        nextElemId += 1
                    oldElemToReplacements[element] = newElemsList

                elif thisLinType!=thisQuadType:
                    # after standard conversion: update element type information
                    linTypesElements.remove(element)
                    quadTypesElements.add(element)
                    self.elType[element] = thisQuadType
                else:
                    raise ValueError(
                        "ERROR in Model.convertToQuadratic(): Special"
                        " conversion of elements of type %s not yet"
                        " implemented." % thisLinType)
                ticker.msg(iCount+1)
            del ticker
            if thisQuadType=='COH3D6':
                msg('Inserting %d new elements of QuadType %s into the model'
                    ' now...' % (len(newElNodes), thisQuadType))
                self.updateElems(newElNodes, thisQuadType)

            # remove this linear type from self.typeEl (if applicable)
            if len(linTypesElements) == 0:
                del self.typeEl[thisLinType]

        # from now on behave like an ordinary dict
        oldElemToReplacements.default_factory = None

        # if new elements have been inserted (only special cases)
        if len(oldElemToReplacements):
            msg('Removing %d elements of LinType %s that have been replaced'
                ' now...' % (lenLinTypesElements, thisLinType))
            # update elsets: add new and remove old elements
            if updateElsets:
                for els in self.elset.itervalues():
                    oldInEls = els.intersection(oldElemToReplacements)
                    for oldEl in oldInEls:
                        els.update(oldElemToReplacements[oldEl])
                    els.difference_update(oldInEls)

            # remove old elements
            self.removeElem(oldElemToReplacements, updateElsets=False)

        return (edgeToNewMidNode, oldElemToReplacements)
    #} end of other element specific methods



    ########################################################################
    #{ other methods extracting parts of the model

    def copy(self):
        """Return an exact copy of self. This is calling the constructor of
        self's type with self as initializer-argument. Which might be a
        subclass of Model.
        """
        return type(self)(self)


    def copySubmodelInSphere(self, centre, radius):
        """return another model with just all elements that are completely
        within the given sphere.

        elements and nodes are copied (a deep copy, not a reference)

        elsets and nsets are copied as well, but nothing else
        """
        mo = Model()
        for node, coord in self.nodeCoords.iteritems():
            if length(vector(coord,centre))>radius:
                continue
            mo.nodeCoords[node] = coord

        for elem, nodes in self.elNodes.iteritems():
            if all([(nn in mo.nodeCoords) for nn in nodes]):
                mo.updateElem(elem, self.elType[elem], nodes)

        connected = self.getNodesFromElset()

        for node in mo.nodeCoords.keys():
            if node not in connected:
                del mo.nodeCoords[node]

        # copy elsets
        allelems = set(mo.elNodes)
        for setname, elset in self.elset.iteritems():
            inset = allelems.intersection(elset)
            if inset:
                try:
                    mo.elset[setname].update(inset)
                except KeyError:
                    mo.elset[setname] = inset

        # copy nsets
        allnodes = set(mo.nodeCoords)
        for setname, nset in self.nset.iteritems():
            inset = allnodes.intersection(nset)
            if inset:
                try:
                    mo.nset[setname].update(inset)
                except KeyError:
                    mo.nset[setname] = inset

        return mo

    def copySubmodelInVolume(self, vol):
        """Return another model with just all elements that are completely
        within the given volume. More precisely: All nodes are being tested.
        For an element to be considered in the volume all of it's nodes must
        be in that volume.

        Elements and nodes are copied (a deep copy, not a reference).

        Elsets and nsets are copied as well, but nothing else.

        @param vol: An object of a subclass of L{Volume<bae.volume_02.Volume>}.
        """
        mo = Model()

        nodesInVol = vol.intersectionIdPointTuples(self.nodeCoords)
        elemsInVol = set(
            elem
            for elem, nodes in self.elNodes.iteritems()
            if all((nn in nodesInVol) for nn in nodes))
        for type_, elems in self.typeEl.iteritems():
            elNodes = dict( (el, self.elNodes[el])
                            for el in elems.intersection(elemsInVol) )
            if elNodes:
                mo.updateElems(elNodes, type_)

        connected = mo.getConnectedNodes()
        mo.nodeCoords.update(
            (node, self.nodeCoords[node])
            for node in connected)

        # copy elsets
        allelems = set(mo.elNodes)
        for setname, elset in self.elset.iteritems():
            inset = allelems.intersection(elset)
            if inset:
                mo.elset[setname] = inset

        # copy nsets
        allnodes = set(mo.nodeCoords)
        for setname, nset in self.nset.iteritems():
            inset = allnodes.intersection(nset)
            if inset:
                mo.nset[setname] = inset

        return mo

    def copySubmodelFromElsets(self, elements, recognizedBlocks=None):
        """Return another model with just all elements in one of the given
        elsets.

        All items are copied (a deep copy, not a reference).

        If no elements are found (unknown elset name, empty elset) an empty
        model is returned. No warning!
        If there are no element or no node definitions to the elements in the
        elsets this is silently ignored.

        @param elements: may be an elset name or a list (or other iterable) of
        elset names and/or element numbers. (Even a mixed list of elset names
        and element numbers is accepted.)

        @param recognizedBlocks: If something like NODE or ELEMENT or a list of
        such (see the correspondig argument of self.read()) only the specified
        items are copied to the new model. (And of course only those items that
        are connected to the given elements.)

        @note: There is a similar method for the more general Mesh class that
        does only a shallow copy for most of the items:
        L{Mesh.partMeshFromElements}().

        @note: Only those materials and connector behaviors are being copied
        that have a correspondint section property definition. Section property
        definitions are only copied if the corresponding elset exists.

        @note: Bugs/Limitations:
         - Only the following items are copied so far:
            elements, nodes, elset, nsets, properties, materials, connector
            behaviors, orientations, initial conditions
         - Initial conditions are copied completely regardless whether the
           addressed elements or nodes actuallay exist in the submodel.
           This is a real flaw and should be corrected, asap.
         - very slow for big models or big elsets (still the case?)
        """
        mo = Model()

        try:
            elemsToCopy = self.getUnionSet("elset", elements)
        except UnionSetInvalidItem, err:
            raise UnionSetInvalidItem(
                "copySubmodelFromElsets() invalid item in the elements"
                " list argument: %s" % err.args[0])

        recognizedBlocks = self._checkRecognizedBlocks(recognizedBlocks)

        if (('NODE' in recognizedBlocks)
                or ('NSET' in recognizedBlocks)):
            # initialize nodestocopy
            ticker = MsgTicker(
                "Collecting nodes to copy. Nb of elements already"
                " processed: %s/%d." % ("%d", len(elemsToCopy)))
            nodestocopy = set()
            for cnt, elem in enumerate(elemsToCopy):
                try:
                    nodes = self.elNodes[elem]
                    nodestocopy.update(nodes)
                except KeyError:
                    continue
                ticker.msg(cnt+1)
            del ticker

        if 'ELEMENT' in recognizedBlocks:
            msg("Copying %d elements to new submodel." % len(elemsToCopy))
            ticker = MsgTicker("... finished %s/%d." % ("%d",len(elemsToCopy)))
            cnt = 0
            for thisType, thisElems in self.typeEl.iteritems():
                elNodes = dict()
                for elem in thisElems.intersection(elemsToCopy):
                    cnt += 1
                    elNodes[elem] = list(self.elNodes[elem])
                    ticker.msg(cnt)
                if elNodes:
                    mo.updateElems(elNodes, thisType)
            del ticker

        if 'NODE' in recognizedBlocks:
            msg("Copying %d nodes to new submodel." % len(nodestocopy))
            ticker = MsgTicker("... finished %s/%d." % ("%d",len(nodestocopy)))
            for cnt, node in enumerate(nodestocopy):
                ticker.msg(cnt)
                try:
                    coords = self.nodeCoords[node]
                    mo.nodeCoords[node] = list(coords)
                except KeyError:
                    pass

        # copy elsets
        if 'ELSET' in recognizedBlocks and len(self.elset):
            msg("Copying elsets to new submodel.")
            for setname, elset in self.elset.iteritems():
                inset = elemsToCopy.intersection(elset)
                if inset:
                    mo.elset[setname] = inset

        # copy nsets
        if 'NSET' in recognizedBlocks and len(self.nset):
            msg("Copying nsets to new submodel.")
            for setname, nset in self.nset.iteritems():
                inset = nodestocopy.intersection(nset)
                if inset:
                    mo.nset[setname] = inset

        # copy properties
        if (len(self.sectionBlocks.intersection(recognizedBlocks))
                and len(self.properties)):
            msg("Copying section properties to new submodel.")
            for setname, property in self.properties.iteritems():
                if setname in mo.elset:
                    mo.properties[setname] = copy.deepcopy(property)
                    try:
                        matname = property['MATERIAL']
                    except KeyError:
                        behaviorname = property['BEHAVIOR']
                        mo.properties.behaviorElset[behaviorname] = setname
                    else:
                        mo.properties.matElset[matname] = setname

        # copy material properties
        if ('MATERIAL' in recognizedBlocks and len(self.material)):
            msg("Copying material properties to new submodel.")
            for name, item in self.material.iteritems():
                if name in mo.properties.matElset:
                    mo.material[name] = copy.deepcopy(item)

        # copy orientations
        if 'ORIENTATION' in recognizedBlocks and len(self.orientation):
            msg("Copying orientations to new submodel.")
            for name, item in self.orientation.iteritems():
                mo.orientation[name] = copy.deepcopy(item)

        # update connectorBehavior
        if ('CONNECTORBEHAVIOR' in recognizedBlocks
                and len(self.connectorBehavior)):
            ("Copying connectorBehaviors to new submodel.")
            for name, item in self.connectorBehavior.iteritems():
                if name in mo.properties.behaviorElset:
                    mo.connectorBehavior[name] = copy.deepcopy(item)

        # update initial conditions
        if ('INITIALCONDITIONS' in recognizedBlocks
                and len(self.initCond)):
            msg("Copying initial conditions to new submodel.")
            mo.initCond = copy.deepcopy(self.initCond)

        return mo
    #} end of other methods extracting parts of the model


def mesh_asAbqModel(self, mapShapeToType=Model.defaultMapShapeToType):
    """Returns self as an instance of abq_model_02.Model, so it can perform
    any actions specific to this class like reading and writing.

    @Note: This method is only available if you also import the module
    bae.abq_model_02
    """

    model = Model()

    # data members, they are always there
    model.nodeCoords = container.NodesDict(self.nodeCoords)
    model.elNodes = container.ElementNodesDict(self.elNodes)
    model.elType = dict(
        (el, mapShapeToType[shape])
        for el, shape in self.elShape.iteritems())
    for shape, elems in self.shapeEl.iteritems():
        model.typeEl[mapShapeToType[shape]] = elems
    model.elShape = container.FakeElShape(model.elType, model.elTypeToShapeSep)
    model.shapeEl = container.FakeShapeEl(
        model.typeEl, model.elShapeToTypes, model.elTypeToShapeSep)
    return model

Mesh.asAbqModel = mesh_asAbqModel
