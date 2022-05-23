"""For reading, writing and manipulating abaqus input file data.

Intended to be imported like a module.

Main content
============
B{L{Model} class}


Usage
=====
  >>> from bae.abq_model_02 import Model
  >>> m = Model() # an optional parameter may state the already opened log file
  >>> help(m) for more information

Remarks
=======
  New structures to be read from an input file can be defined in
  module abq_model_02.reading.

Known bugs/problems
===================
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
  - Several *MPC,USER commands will be unified to one with only the last
    specified MODE option being taken for all those MPCs. On output all
    MPCs with a numeric type will be grouped under a *MPC,USER command
"""


__version__ = "2.31"

_version_history_ = """
Versions:
=========

abq_model_02
2.0 : GP: derived from Version 1.8
      changed interface: other attribute names, standard to be for all modules;
      elset, nset -> python sets; elset, nset names case insensitive (class
      AbqNameDict); attributes of Model-object always there; Model.read() from
      already open file or anything iterable; removed nset_BC attribute;
      remove_elem() for many elements; update_elem() generates element number;
      documentation moved to class Model - docstring
      added: copy_submodel_from_elsets();
      renamed: many functions, arguments, variables from this_naming_form to
      thatNamingForm
      added: write elsets in blocks with the *ELSET, GENERATE command if
      possible

2.1 : GP split internals into abq_model_internals_02.py for clearity
      Model.read("...", "ELSET") now also reads elsets from *ELEMENT, ELSET=...
      and similar for NSETs
2.2 GP added: Model.convertToTets(), checkModelModuleVersion(),
           Model.checkElType(), Model.forceElset(), Model.forceNset()
       changed: copySubmodelFromElsets() accepts recognizeBlocks argument,
           getNodesFromElset() now accepts everything suitable for
           getUnionSet() as elements argument. It does no longer raise an
           exception but writes a warning to the log file if an element in
           the sets is not defined in the model (optional).
2.3 : BF added functionality to define a surface that comprises exterior surface
         of an single or group of elsets.
      GP fixed: sensible error message on reading rubbish from input file
      BF added functionality to boolean two surfaces.
      BF added functionailty to surface that comprises exterior surface nodes
         of an single or group of elsets.
      BF redefined freeSurfFromElset(), freeSurfNodesFromElset() and
         booleanSurface().
         freeSurfFromElset() replaced with freeSurfElemsFromElset().
         freeSurfNodesFromElset() now returns a set of nodes labels
         booleanSurface() decomposed into surfaceAdd(), surfaceSubtract() and
         surfaceIntersect()
         automatic creation of new surfaces removed from
         freeSurfElemsFromElset(), surfaceAdd(), surfaceSubtract() and
         surfaceIntersect(). Output from these functions can be directly used
         with
         Model().read('myModel.inp').surface['surfName']=freeSurfElemsFromElset()
         Model().read('myModel.inp').surface['surfName']=surfaceAdd()
         ...
      GP added: convertToQuadratic(), Model.mpc
      GP fixed: copySubmodelFromElset faster, with diagnostic output
2.4 : GP added getElCentroidDict(), getElGaussPointsDict()
2.5 : GP changed: make it a subclass of Mesh
2.6 : GP added: renumber
2.7 : Vlad and Gero added: writeElsetNameCsv
2.8 : GP added: insertElems
2.9 : GP added: insertNodes, __str__
2.10 : GP added: insertModel, Orientation, Connector Behavior and suboptions,
                 material is now a permanent option
          added: rotate, rotateX, ..., Model() as copy constructor
          changed: InitialCond-classes now have a label attribute instead of
                 node or element and a read-only property to get the label as
                 self.node or self.element.
2.11 : GP added: orientationTolerance argument to Model.insertModel()
2.12 : AF added: writeAllElemsOfCertainType to Model
2.13 : GP changed to bae.log_01
2.14 : GP added (subclass-aware) copy method.
2.15 : GP added: Model.read() can read zip files
2.16 : GP added: includeMidNodes option for Model.freeSurfNodesFromElset()
2.17 : GP fixed: make elShapeToTypes values immutable
                 (i.e. "TET_L" : frozenset(('C3D4', )) ... )
2.18 : GP added: mergeNodesOnElems
2.19 : GP added header argument to Model.write()
2.20 GP fixed: recognizedBlocks argument of Model.read()
2.21 GP replaced getNodesFromElset() by getConnectedNodes()
2.22 GP added copySubmodelInVolume(), improved freeSurfNodesFromElset()
2.23 GP added *PLASTIC option (read and write)
2.24 GP added withSummary argument to Model.write()
2.25 GP fixed Model.getQuadGaussPointsDict() for integer node coords
2.26 GP added Model.__str__ now includes nb of material definitions
2.27 GP speed up reading, last version with abq_model_02 being a module

versions of abq_model_internal_02:
2.0 : GP: derived from Version 1.8 ...
      added: write elsets in blocks with the *ELSET,GENERATE command if possible
      added: reading and writing MPCs
      added: writing of time points, amplitudes, section properties
2.1 : GP: added writing of BCs
2.2 : GP: added removeSection
2.3 : GP: added InitialCondTypeField
2.4 : GP: added Orientation, ConnectorBehavior
          changed Model.material from being a dict to being an AbqNameDict.
2.5 : GP: added Amplitude, definition=smooth step
2.6 : GP changed to bae.log_01
2.7 : GP fixed writing of initial conditions with GEOSTATIC=False
2.8 : GP fixed: FakeShapeEl raises KeyError if no element of given shape,
         added FakeShapeEl.get()
2.9 : GP fixed: don't modify FakeShapeEl.elShapeToTypes values on querying
         FakeShapeEl anymore.
         changed: _version_ to __version__
2.10: GP added: both element shape versions, specific like "TET_L"
         as well as general like "TET" can be queried with FakeShapeEl
2.11 GP added FakeElShape-methods __iter__(), iterkeys(), iteritems()
        ... needed for elShape.update(other.elShape.iteritems()) in
        DynaModel.__init__
2.12 GP speed up reading: replace dataline-methods by readDataLines
        changed (probably: fixed): connector hardening, connector damping,
        connector stop, connector lock previously had the data lines not
        converted to float.

common version recorded in __init__.py
2.28 GP Conversion to package with submodules.
2.29 TR Speed up reading by:
     - passing some string processing to StringIO before parsing
     - using pandas or numpy (if available) for stringToNum conversion
2.30 GP added TYPE-option for USERMATERIAL and *USER DEFINE FIELD card
2.31 GP added new elshape TETHT_L for DC3D4 linear heat tranfer tet element
    with four Gauss points
"""

_todo_ = """
known bugs/problems:
====================

- none, of course.

ideas:
======

- For reading speed optimizations: see module reading.py
- Interface improvements (make it cleaner, easier to understand and use):
  . The constructor Model() should take one optional positional argument to
    initialize (i.e. another Model object or an input file or input string,
    maybe a file name as well if it ends in .inp). The logfile argument should
    be an optional keyword argument. I.e.
     >>> def __init__(init=None, logfile=None):
  . simplify the storage of initial conditions:
    { type : { label : value(s), label2 : values2, ... }, type2: ... }
    I.e.: { "TEMPERATURE" : { "MyNSET":1.0, 500142:0.0 },
    "STRESS" : { "MyELSET":[1E6, 2E7, ...], 20053, [...] } }
    Problem: options. Store together with values? Or make different options
    induce different types?
    Problem: order may matter, i.e. if an elset gets stresses then individual
    elements out of this elset get other values. Or have intersecting elsets.
    Why? To simplify handling and merging (insertModel).
    And: Don't subclass floats.
"""

from .model import Model, checkModelModuleVersion,\
    ReadError, UnionSetInvalidItem, ModelInconsistent, UnknownElType
