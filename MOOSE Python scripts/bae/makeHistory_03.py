r"""
Processing the model history input and generating an input file as output.

For plans, design flaws, bugs ... see the __todo__ -string.
"""

## ---------------------------------------------------------------------------
## version control
__version__ = '3.30'

_version_history_ = """
3.01 GP new - created new from 2.??
3.02 GP added getLoadHistRow(rowIdx)
        fixed the double elset in sequence csv issue.
3.03 GP incompatible change: removed logfile argument, changed to bae.log_01
        functionality (need to remove Model()-logfile argument once
        bae.abq_model_02 has been updated as well, search for logfile)
3.04 GP added colToAmplitude function
        incompatible change: getStrTimePointsDef returns timepoints for all odb
        steps if odbStep=None (or not specified) (This feature has most
        probably not yet been used.)
3.05 GP fixed: elsets in last frame for seqElsets now enclosed in quotes
3.06 GP added: sequenceElsetsIter
3.07 GP fixed: ignore case of elset names
3.08 GP added: possibility to leave away disjointSetsFileName argument of
        LoadHistory.updateDisjointElsets() for no elsets output.
3.09 GP fixed: recognize the optional elsetColumns argument to
        LoadHistory.disjointSeqElsets() in the call to removeElsetFromStep() as
        well. Unfortunately for speed considerations this introduced an
        incompatible interface change to removeElsetFromStep(). If you encounter
        this as a problem give self.elsetColumns as additional argument to
        removeElsetFromStep().
        added: getElsetsAmp() for temperature amplitudes from nset columns
3.10 GP: fixed ignore double elsets in same step
3.11 GP: added exportSeqSurfOnePerStep from version 1.3 of
        py_abq_post/exportSequence/exportSeqSurfFromAbqModel.py
3.12 GP: changed the frame identifier in the output file name and object name
        from something like S2F005 to something like F2005.
3.13 GP: added generateAmplitude()
3.14 GP changed: exportSeqSurfOnePerStep now uses the objectNameTemplate
        argument for the filename ending as well.
3.15 GP added dynaAddSequence
3.16 GP incompatible change to exportSeqSurfOnePerStep: renamed
        postDataFileName to postData
        changed exportSeqSurfOnePerStep to use sequenceElsetsIter and added
        ignoreElsets, ignoreElems, addExtraSurfaces arguments
3.17 GP added colsToInitStressSolid argument to LoadHistory.dynaAddSequence()
3.18 GP changed: elset names and elset columns are all converted to upper case
        and are therefore not case sensitive anymore.
3.19 GP getStrSdvIni can now take two columnNames in the colToSDVini argument
        to create split elsets.
3.20 GP added loadHistory.updateElsetNames()
3.21 GP added sequenceElsetsIterFromCsv()
3.22 GP added restartFrom-argument to writePostDataFile()
3.23 GP added loadHistory.getFrameTimes()
3.24 GP added loadHistory.writeFrameTimes()
3.25 GP fixed loadHistory.intersectColumns: can now be called multiple times
3.26 GP changed sequenceElsetsIter to also consider elsets assigned to frame 0
3.27 GP added LoadHistory.odbStepLastFrameDict and ....getStepFrameToATN()
              write frameLists and stepFrameToCsStep to postData
3.28 GP added LoadHistory.frameRowIter, LoadHistory.writeCaveCouplingTiming
              and writeCaveCouplingRelationTable
3.29 GP added frameListHalfYearly to postData
3.30 GP added added useStepName (new "Step-2"-name convention) and
        frameListQuarterly to writePostDataFile
"""

todo = """
ERROR: (bugs, design flaws...)
==============================

 - frame names in the loadHistory-csv file that can be converted to floats
   may lead to nonsense results.

 - LoadHistory.elsetColumns identifies columns
   (a) with elsets to be considered for the LoadHistory.disjointSeqElsets()
       functionality,
   (b) with names that are to be treated case insensitive (might also be nsets)
   (c) that are to be treated as "listColumns": each item is a list even if
       there is only one column of that name.
   Not all of the features are wanted at the same time in every case! For
   example a column might contain elsets or nsets that will never be made
   disjoint, for example ground support installer nsets.
 - By placing an elset name in a particular cell in the loadhistory table it is
   assigned to a "column" in the sense of a type. A particlar entry in the
   colToSDVini and colToGravAmp arguments is chosen according to the column.
   At the same time a particular mining step is identified by the row.

   But: An elset might be of a particular type and belong to a particular
   column in that sense without being assigned to a particular mining step.
   E.g. NOTMINED.

   And: An elset (or nset) might ba associated with a particular mining step
   without necessarily belonging to a particular type. E.g. for amplitudes for
   particular elsets. ???

 - Handling of the elset names and elset column names is a mess. There should
   be a distinct, common, reasonable and intuitive point of conversion to
   uppercase (I think). Column names also serve as variables. For that we may
   want to conserve case. So we might have to have both: the original name with
   case preserved and a upper case converted copy.


Plans:
======

from discussions mid June 2020:
 - streamline column names, define standard columns (?)
   . outputNb (sequential number in e.g. vtk file names "_F001.vtk",
   . outputName (most of the times a date, to be used in seq-iv-file names),
   . odbStep, odbFrame
   . stepLength == miningStepLength == varDT_frame (uuuh)
 - isOdbFrame argument shall work automatically ?


 - We might want to move the makeHistory modules to bae.utils.
   For compatibility we'd still want it to be available as bae.makeHistory_03
   this seems to be achievable by adding the following line to bae.__init__.py:
   __path__.append("bae/utils")
   Then *all* modules in bae/utils would be available as bae/....

   But: readPostData and sequenceElsetsIter don't fall into that category,
   bae.utils.... doesn't suit them well.

   Also note the exportSeqSurfOnePerStep function.

   We might want to step further to makeHistory_04 when we split off the
   multicolumn colToSDVini-argument-functionality as suggested below.
   We might want to do some cleanup with amplitudes-functions, anyway?
   We might want to change the name alltogether....?
   ... something like bae.utils.sequence?

 - Or instead move (at least some core functionality) to a new sub-module
   abq_model_02.sequence? Maybe have a Sequence (or LoadSequence or other)
   class there with methods: sequenceElsetsIter(), getStrSdvIni(),
   getGravAmpDload() and other. Let this class get the seq. info from different
   sources, e.g. loadHistory-csv, just sequence elsets, Abaqus input files.
   Maybe also include the drawschedule-stuff?

 - getVelPatternBCFromGivenDispl()...
   . velocity_pattern = list of (time, amplitude)-tuples each item can be
     evaluated as always with makehistory...
   . amplitude values are scaled such that the resulting displacement is as
     prescribed (including the scale factor)
   . remove the middle of three identical sucessive values
   . note that the resulting times is a successive list
   . see Prominent Hill GS 2013/ runs earlier than R13, there is a
     makeHistory_BCAmpl_00.py module with some starters

 - split elsets functionality:
   . it's done for getStrSdvIni(), getGravAmpDload() and writePostDataFile()
     with a multicolumn colToSDVini and colToGravAmp-argument.
   . not yet implemented for updateDisjointElsets (do we need it? It might not
     be appropriate, not flexible enough...)
   . Currently disjointSeqElsets has to be called with warnNotUsedElsets=False.
     Can we get rid of that?
     I.e. add warning about element sets not in the load history from
     self.disjointSeqElsets() or updateDisjointElsets().

 - Don't use DictTableFromCsvFile anymore.
   Frame names are supposed to be strings in the csv file. DictTableFromCsvFile
   converts all columns to float if possible. It'd be better to only convert
   specified ones.
   Elset columns are list-columns. That's already implemented in DictTable...
   Make it less flexible and simpler.
   But the evaluateString - function is a nice feature....

 - Make the default usage such that it works without elset-list in __init__.
   Then optionally read a loadhistory table with optional elset columns.
   Elsets could also be taken from the seqElsets input file. And their names be
   interpreted. That should be the default behaviour/use case.

 - have all column labels for elsets start with a certain prefix like "elsets_"
   or "e_" in order to automatically identify them as such. Only those columns
   shall be treated as multi-column-columns. I.e. save the elsetColumns
   argument
"""

## ---------------------------------------------------------------------------

import os
import re
from cStringIO import StringIO
from itertools import izip
from bae.future_01 import defaultdict, OrderedDict
from bae.misc_01 import DictTableFromCsvFile, RowIterFromCsvFile,\
    groupIter, Container, filterDateListYearly, filterDateListInterval
from bae.abq_model_02 import Model
from bae.abq_model_02.container import Amplitude, TimePointsDict, ElsetsDict
from bae.volume_02 import VolumeFromMesh  # for exportSeqSurf()

from bae.log_01 import msg, log


## ---------------------------------------------------------------------------
class LoadHistory(object):
    r"""An object that stores the load history for the analysis. Meant to
    generate initial conditions for SDVs (SDVini), amplitudes, dloads and so
    on.


    Value evaluation
    ================

    Some strings (e.g. the miningStepLength argument if it's a string) are
    evaluated as python expression. In that expression all the values of
    the current row in the load history table (aka the load history csv file)
    are accessible as variables named like the column.
    I.e. if you have a column labeled 'stopetime' with a value 5 in the
    current row then the variable stopetime evaluates to 5 and e.g. an
    expression '0.5*stopetime' would be a valid expression for the
    miningStepLength argument yielding 2.5 in this case.

    Some extra columns are being added by __init__ automatically. See the
    description of the instance variable L{loadHistTable} for details.

    Additionally you can access other variables contained in the
    extraLocals dict. Specify extraLocals=locals() in the argument list
    to make available all current local variables.

    Probably the most frequent use will be just to specify a string
    identifying the column in the csv file from which to read the value.

    Note that spaces in column names will be replaced by an equal number
    of underscores ("_") to make them valid python identifiers.

    Hint: To convert a possible zero value into some sensible minimum you can
    use the expression int(not(...)). For example suppose the miningStepLength
    shall be 0.5*stopetime but should be set to 5 in case stopetime is zero
    then set miningStepLength to "0.5*stopetime + 5.0*int(not(stopetime))".

    This way you can also introduce more complicated conditional expressions
    but remember that each term contributes to the sum in any case, so you
    might need a conditional factor for each single term. For example:
    "0.5*stopetime*int(stopetime>5.0) + 5.0*int(stopetime<=5.0)" will yield
    the maximum of 0.5*stopetime and 5.0. Which you could also achieve with
    a max(...) expression, btw.


    Subclassing
    ===========

    In case there are special variations on how the load history should be
    processed subclassing this is the recommended way of doing it.

    Whenever self.loadHistTable attribute is being accessed it should be done
    using its rowIterColToList() or getRowColToList() methods with
    self.elsetColumns as listColumns argument and removeNones=True. I.e.:
     >>> for row in self.loadHistTable.rowIterColToList(
     ...     listColumns=self.elsetColumns, removeNones=True):
     >>>
     >>> specialRow = self.loadHistTable.getRowColToList(
     ...     rowIdx=5, listColumns=self.elsetColumns, removeNones=True)

    The methods self.getLoadHistTableIter() and self.getLoadHistRow(rowIdx)
    have been supplied for this purpose for convenience. So instead of the
    former you would now do:
     >>> for row in self.getLoadHistTableIter():
     ...     doSomething(row)
     >>>
     >>> specialRow = self.getLoadHistRow(5)


    Split elsets feature
    ====================

    This is for the case that excavation (T1 and T2) are sequenced independent
    from fill (T3) or deletion time (T4). The LoadHistory-object has methods to
    split/intersect elsets in order to assign the correct timings to each part.

    Example: It's been used for the Prominent Hill project PH2017099 for the
    dumps in the pit. We have pit elsets to be switched to void early in the
    sequence. And dumps elsets containing only elements that also belong to
    some pit elsets. Those dumps elsets represent the part of the pit that is
    being back filled at a much later time point. Of course there is also a
    (large) part of the pit which is not refilled. The split elset feature
    creates intersection elsets of the pit and dump elsets and assigns them
    their individual gravity amplitude and T1-T4 timings via SDV-initial
    values. Remaining pit-only elsets are also treated separately.

    Second example: Kidd Creek project KIDD2019022 has stopes that are refilled
    after having stood open for a certain period of time variing from stope to
    stope. So for each stope elset we have a mining step for the excavation and
    independent from that another ming step in which the refill takes place.
    Now we have stope-excavation elsets and stope-refill elsets that all
    together cover the same region. The split elset feature again creates
    intersection elsets and proper timing data.

    Procedure: There is supposed to be a collection of elsets that together
    contain all elements that are being mined in any way. In the PH example
    above that would be the pit elsets and all drives and stopes but not the
    dumps in the pit. We will be making those elsets disjoint and derive the
    remaining "NOTMINED" elset from that. In a second step we'll look at the
    dump elsets alone and make them disjoint separately from the rest. This is
    achieved by two separate calls to the method L{disjointSeqElsets}.

    For the first call the parameter checkOnlyColumns contains the list of
    column names of pit, stopes, drives and so on, all except the dumps.
    Additionally we have the elsets "NOTMINED" passed through the
    additionalElsets parameter. For the second call we have only the dumps
    given in the parameter checkOnlyColumns and nothing for the parameter
    additionalElsets.

    After making the elsets disjoint we're creating the required intersection
    elsets by means of the method L{intersectColumns}. We want to intersect
    the pit and dumps elsets and therefore pass the corresponding column names
    to the parameter elsetColumns. Those intersection elsets are also being
    recorded in a proper way to later appear in the seqElsets of the postData
    output and in the gravity amplitudes and sdvini creation.

    From now we follow the standard procedure. The methods L{getGravAmpDload},
    L{getStrSdvIni}, and L{writePostDataFile} can be used as usual and their
    output takes into account the split elsets.

    The colToGravAmp and colToSDVini arguments to L{getGravAmpDload} and
    L{getStrSdvIni} methods can use the "value evaluation" feature described
    above: values are derived from expressions which might contain variables
    named like a column of the history table (a.k.a. loadhistory csv).
    Thereby values may differ between mining steps. With the split elsets
    feature we have intersection elsets automatically generated by
    L{intersectColumns}. For these intersection elsets there are two different
    mining steps from which values can be taken. The corresponding variable
    names are "elset column label"."value column label". See the example below
    and note the colToGravAmp-item for the combined column labels, e.g.
    ["excav_XC", "fill_concrete"]. Note that the variable names are case
    sensitive: The prefix in the expression must match exactly the spelling
    in the column-label-pair of the colToGravAmp resp. colToSDVini arguments.
    Whereas the the column lookup in the loadhistory-table is done disregarding
    upper and lower case.

    Note: The elsetColumns argument to L{intersectColumns} must meet an
    identical double-column-item in the colToGravAmp argument to
    L{getGravAmpDload}() and an identical double-column-item in the colToSDVini
    argument to L{getStrSdvIni}(). I.e. if you have ["exc_XC", "fill"] as
    elsetColumns for L{intersectColumns} then you also need a
    ["exc_XC", "fill"] entry in the colToSDVini parameter for getStrSdvIni()
    and a ["exc_XC", "fill"] entry in the colToGravAmp argument of
    L{getGravAmpDload}(). Otherwise the grav amps and sdvini lines for the
    split elsets will be missing.

    Note: With the current implementation and the split elsets feature you
    should suppress warnings about element sets not in the load history from
    L{disjointSeqElsets}(). Use warnNotUsedElsets=False argument. The reason is
    that you are calling L{disjointSeqElsets}() only with a subset of all
    elsets at a time.

    Note: With the split elsets feature each mining step with elsets assigned
    needs a name even if there is no odb-frame. This is because combined elsets
    get names derived from the mining step name. (And of course as usual: The
    mining step names must be unique.)

    Note: With the split elsets feature L{updateDisjointElsets}() does not
    work. L{updateDisjointElsets}() is a convenience wrapper for
    L{disjointSeqElsets}() but here we need to call it twice...

    Note: If there is no collection of elsets that together contain all elements
    to be mined then the procedure of making elsets disjoint is flawed. For
    example if we have stopes that are eventually being eaten up by a pit
    then there are intersections between pit and stopes, i.e. elements that are
    stopes and later pit. And we have stopes that are never becoming pit and we
    have pit that has never been a stope. In this case the described procedure
    would not be able to make drive elsets disjoint with pit and with stopes as
    well. This certainly can be resolved but it's not covered yet...

    Note: The intersection elsets created by L{intersectColumns} can't be
    intersected again with another set of elsets because those new elsets don't
    appear in the elset columns of the loadhistory table. Intersections can
    only be made between elsets that were in the loadhistory table from the
    very beginning.

    Example:
     >>> colToGravAmp = [
     >>>     (["excav_XC", "fill_concrete"],
     >>>      [("excav_XC.totaltime", 1.0),
     >>>       ("excav_XC.totaltime+excavDuration", 0.0),
     >>>       ("fill_concrete.totaltime+0.8", 0.0),
     >>>       ("fill_concrete.totaltime+0.9", 0.85),
     >>>      ]),
     >>>     (["excav_XC", "fill_cave"],
     >>>      [ ...
     >>>      ]),
     >>>     ("excav_XC",
     >>>      [("totaltime", 1.0),
     >>>       ("totaltime+excavDuration", 0.0),
     >>>      ]),
     >>>     ...
     >>>     ]
     >>> colToSDVini = [
     >>>     (["excav_XC", "fill_concrete"], [
     >>>         0,0,0,0,0,0,0,1,1,
     >>>         "excav_XC.totaltime",
     >>>         "excav_XC.totaltime+1.6",
     >>>         "fill_concrete.totaltime+0.8",
     >>>         9.9E6]),
     >>>     (["excav_XC", "fill_cave"], [
     >>>         0,0,0,0,0,0,0,1,1,
     >>>         "excav_XC.totaltime",
     >>>         "excav_XC.totaltime+1.6",
     >>>         "fill_cave.totaltime+0.8",
     >>>         9.9E6]),
     >>>     ("excav_XC",[
     >>>         0,0,0,0,0,0,0,1,1,
     >>>         "totaltime",
     >>>         "totaltime+1.6",
     >>>         9.9E6,
     >>>         "totaltime+1.6"]),  ]
     >>> extraElsetSDVList=[
     >>>     ("NOTMINED", [
     >>>         0,0,0,0,0,0,0,1,1,
     >>>         9E9,9E9,9E9,9E9]),
     >>>     ]
     >>> ...
     >>> loadhistory = LoadHistory(
     >>>     loadhistoryCsv=loadhistoryCsvName,
     >>>     elsetColumns=["pit", "drive", "excav_XC", "fill_concrete",
     >>>                   "fill_cave", ...],
     >>>     ..., )
     >>> ...
     >>> modelSeqElsets = Model()
     >>> modelSeqElsets.read(seqElsetsBasePath)
     >>> seqElsets = modelSeqElsets.elset
     >>> msg("Got %d sequence elsets." % (len(seqElsets))
     >>>
     >>> # create disjoint elsets, remove empty elsets from the load history
     >>> disjointSets = dict()
     >>> disjointSets.update(loadhistory.disjointSeqElsets(
     >>>     additionalElsets=["NOTMINED"],
     >>>     checkOnlyColumns=["pit", drive", "excav_XC", ...],
     >>>     elsetsDict=seqElsets, warnNotUsedElsets=False, ...))
     >>> disjointSets.update(loadhistory.disjointSeqElsets(
     >>>     checkOnlyColumns=["fill_concrete", "fill_cave",],
     >>>     elsetsDict=seqElsets, warnNotUsedElsets=False, ...))
     >>>
     >>> # create not mined elset(s)
     >>> allMined = set().union(* disjointSets.itervalues())
     >>> modelSeqElsets = Model()
     >>> modelSeqElsets.read(notMinedElsetPath)
     >>> notMinedElsets = modelSeqElsets.elset
     >>> for elset in notMinedElsets.itervalues():
     >>>     elset.difference_update(allMined)
     >>> disjointSets.update(notMinedElsets)
     >>>
     >>> # split elsets
     >>> disjointSets = loadhistory.intersectColumns(
     >>>     disjointSets, elsetColumns=("excav_XC", "fill_concrete"))
     >>> disjointSets = loadhistory.intersectColumns(
     >>>     disjointSets, elsetColumns=("excav_XC", "fill_cave"))
     >>>
     >>> # export new elsets to Abaqus input file
     >>> mout = Model()
     >>> mout.elset.update(disjointSets)
     >>> mout.write(seqElsetsSplitPath, withSummary=True)
     >>> msg("Wrote %s" % seqElsetsSplitPath)
     >>>
     >>> # grav amplitudes and dload, include new split elsets
     >>> gravAmp, gravDload = loadhistory.getGravAmpDload(
     >>>     colToGravAmp, extraLocals=locals(), time="'TOTAL'")
     >>> sdvIni = loadhistory.getStrSdvIni(
     >>>     colToSDVini, extraLocals=locals(),
     >>>     extraElsetSDVList=extraElsetSDVList)
     >>>
     >>> # create post data, uses newSplitElsets
     >>> loadhistory.writePostDataFile(postDataFilePath, postTypeCols)
     >>> msg("Wrote %s" % postDataFilePath)

    @ivar loadHistTable: A L{bae.misc_01.DictTableFromCsvFile}-object containing
       all data from the loadhistory csv file.

       Additionally to all columns in the csv file there are four extra columns
       in this table: odbStep, odbFrame, steptime, totaltime, miningStepLength,
       miningStepName.

       odbStep and odbFrame hold the odb step and frame number of the mining
       step. odbStep=2 corresponds to "Step-2" in the odb. All mining steps
       without a dedicated odb frame (i.e. isOdbFrame==False) will be assigned
       to the following odb frame number. (Because the effects will be visible
       in the following odb frame.)

       steptime and totaltime hold the time at the beginning of the
       corresponding mining step (=odb frame) relative to the beginning of
       the odb step and relative to the beginning of the analysis respectively.

       miningStepLength contains the actual length of the mining step. This is
       being computed from the corresponding argument to L{__init__} and may
       be a duplicate of another column already contained in the csv file. But
       it may as well contain computed values. See self. L{__init__} for
       further explanations.

       miningStepName contains the names of each mining step. This is
       being computed from the corresponding argument to L{__init__}. See
       self. L{__init__} for further explanations.

       If there were already columns with those names in the original csv file
       those columns will be renamed in the resulting loadHistTable: The names
       will be preceded by "orig_". I.e. suppose there is a column labelled
       "miningStepLength" in the loadhistory csv file. Then there will be a
       column "orig_miningStepLength" containing the original content of this
       column and a column "miningStepLength" as described above.

       Note: Elset columns whose names only differ in case will be merged.
       I.e. columns "Pit" and "PIT" will be merged under the name "PIT".

    @ivar elsetColumns: A set of column names containing the elsets for the
       loadhistory. The column names are stored as uppercase, regardless of
       the case of the items in the elsetColumns argument to the constructor.

    @ivar odbStepLengthDict: A dict {odb step number: step length},
       for each odb step found in the loadhistory csv file.

       In case you need a list of all odb steps in the load history
       sorted(odbStepLengthDict) would do.

    @ivar odbStepLastFrameDict: A dict {odb step number: last frame number},
       for each odb step found in the loadhistory csv file.

    @ivar totalAnalysisTime: A float. The sum of all mining step lengths.

    @ivar stepColToSplitElsets: dict {(stepName, column label): [list of names
       of split elsets]} for the split elsets feature, only available after
       call to self.{intersectColumns}.
       L{writePostDataFile} reads this data to add the split elsets to the
       sequence elsets at the right mining step.

    @ivar colsStep1ToStep2Elsets: dict { (elset-columns-tuple, first rowIdx):
       list of (second rowIdx, combined elset name)-tuples } for the split
       elsets feature, only available after call to self.{intersectColumns}.
       L{getGravAmpDload} and L{getStrSdvIni} use this data to find combined
       elsets (from the intersection of one from each column in the
       elset-columns-tuple).

    @Note: B{IMPORTANT} Frame names a.k.a. mining step names in the loadhistory
       csv file must not be convertible to floats. Otherwise there is trouble
       ahead. This is a bug and not a feature. But it's rather difficult to
       correct...
    """

    #------------------------------------------------------------------
    #{ basic functions

    def __init__(
        self, loadhistoryCsv,
        elsetColumns=["Drives", "Cave", "Stopes", "Pit"],
        odbStepNb="StepNb",miningStepLength="Steplength",miningStepName="Name",
        isOdbFrame=True,
        extraLocals={},
        startTotalTime=0.0
        ):
        r"""
        Constructor. Reads the data for this loadhistory from the given csv
        file.

        @param loadhistoryCsv: file name or open file containing the data.
          Note: Don't have column names only differing in case! I.e. having
          a column "Pit" *or* "PIT" is ok but never have both at the same time!

        @param elsetColumns: A list of columns labels that identify all columns
          in the loadhistoryCsv-file that contain sequence elset names. Those
          names are treated in a case insensitive manner.

          It is essential that all columns that contain elset names are listed
          in this parameter. Elset names are case insensitive.

          B{Missing elset column labels may result in columns being ignored or
          elsets not being recognized.}

          Note: A warning will be issued if an elset column mentioned in this
          parameter elsetColumns actually does not exist in the loadhistory-csv
          file.

        @param odbStepNb: step nb in the odb

        @param miningStepLength: Either the step length (a float), in this case
          this is the length of all steps.  Or it may be a string evaluating to
          the step length. See the note on L{value evaluation<LoadHistory>}.

        @param miningStepName: A string evaluating to the step name. See the
          note on L{value evaluation<LoadHistory>}. If you want to take it from
          a column then supply that columns name. If you want the same name for
          all steps supply the name enclosed in quotation marks, i.e.:
           >>> ... miningStepName = "'OZM'"

          Here is a real fancy example for this parameter. Given elset names
          in "Pit" and "Fill" columns that contain a date after some prefix
          ending in an underscore like "PIT_YR2010_M09" or "DRV_YR2008_M12".
          Then the following string will extract the date from the first elset
          name in the "Pit" and "Fill" columns. If there is no elset in those
          columns for a particular step it defaults to "XXXX". Try to read
          this:
           >>> ... miningStepName="(Pit+Fill+['y_XXXX'])[0].split('_',1)[1]"

        @param isOdbFrame: This may be or evaluate to (see the note about
          L{value evaluation<LoadHistory>}) True, a non zero number, "yes" or
          "on", in that case the corresponding load history step will be added
          to a time points list named LoadHistory.defaultTPointsName. If it
          evaluates to False, zero, "no" or "off" this step is not
          included in the list.

        @param extraLocals: A namespace dictionary that contains additional
          variables that should be accessible during value evaluation, see
          L{value evaluation<LoadHistory>}. Specify extraLocals=locals()
          in the argument list to make available all current variables.

        @param startTotalTime: optional offset for totaltime

        @note: The order of columns with different names is not significant.
        """

        # initialize the data table form the csv file
        # try to determine the file name for diagnostic output
        self.loadHistTable = DictTableFromCsvFile(
            loadhistoryCsv, replaceSpaceBy="_", warnIgnoreCols="message")
        if isinstance(loadhistoryCsv, basestring):
            loadhistoryName = "the loadhistory table %s" % loadhistoryCsv
        else:
            try:
                loadhistoryName = (
                    "the loadhistory table %s" % (loadhistoryCsv.name))
            except AttributeError:
                loadhistoryName = "the loadhistory table"

        # store elsetColumns permanently (disregarding case)
        self.elsetColumns = set(c.upper() for c in elsetColumns)
        # filter elsetColumns: only take what is there in the csv file
        missing = self.elsetColumns.difference(
            c.upper() for c in self.loadHistTable.iterkeys())
        if missing:
            self.elsetColumns.difference_update(missing)
            msg("WARNING: Elset column(s) %s not found in %s. Please check"
                " spelling and possibly adapt the elsetColumns parameter to"
                " the LoadHistory-constructor."
                % (",".join(sorted(missing)), loadhistoryName))
        # rename elsetColumns to upper case
        for colKey in dict.keys(self.loadHistTable):
            if colKey is None:
                continue  # ignore columns without title
            colKeyUp = colKey.upper()
            if colKeyUp not in self.elsetColumns:
                continue  # ignore all columns but elsetColumns

            if colKey != colKeyUp:
                if colKeyUp in self.loadHistTable:
                    # merge double columns like "PIT" and "pit" in the csv file
                    merged = self.loadHistTable.mergeCols((colKey, colKeyUp))
                    del self.loadHistTable[colKey]
                    dict.__setitem__(self.loadHistTable, colKeyUp, merged)
                else:
                    # just rename otherwise
                    self.loadHistTable.renameCol(colKey, colKeyUp)

        # rename elsets in elsetColumns to upper case
        def upper(s):
            if s is None:
                return None
            else:
                return s.upper()
        for colKey, column in self.loadHistTable.iteritems():
            if colKey not in self.elsetColumns:
                continue  # ignore all columns but elsetColumns
            if isinstance(column[0], list):
                column[:] = (map(upper, item) for item in column)
            else:
                column[:] = (upper(item) for item in column)

        #-- sum up timing, evaluate mining step names
        totaltime = startTotalTime
        steptime = 0.0
        oldOdbStep = 0
        odbFrame = 1

        steptimeList = list()
        totaltimeList = list()
        stepLengthList = list()
        odbStepList = list()
        odbFrameList = list()
        stepNameList = list()
        self.odbStepLengthDict = dict()
        self.odbStepLastFrameDict = dict()
        for rowCnt, row in enumerate(self.getLoadHistTableIter()):

            errorLocStr = "the %dth data row of %s" % (rowCnt+1,loadhistoryName)

            # the locals argument to the eval() call to evaluate non numeric
            # values in several arguments. E.g. stepLength
            evalNamespaceDict = dict(extraLocals)
            evalNamespaceDict.update(row)

            # evaluate stepLength
            stepLengthVal = self.evaluateString(
                miningStepLength, evalNamespaceDict,
                "the step length in %s" % errorLocStr)
            try:
                stepLengthVal = float(stepLengthVal)
            except (ValueError, TypeError):
                raise ValueError(
                    "The step length in %s does not evaluate to a numeric"
                    " value (float): string <%s> evaluates to <%s>"
                    % (errorLocStr, miningStepLength, stepLengthVal))

            # evaluate odbStepNb
            odbStepVal = self.evaluateString(
                odbStepNb, evalNamespaceDict, "the odb step in %s" %errorLocStr)

            # evaluate isOdbFrame
            isOdbFrameFlag = self.evaluateStringAsFlag(
                isOdbFrame, evalNamespaceDict,
                "the odb frame nb in %s" % errorLocStr)

            # evaluate stepName
            stepNameVal = self.evaluateString(
                miningStepName, evalNamespaceDict,
                "the step name in %s" % errorLocStr)
            if stepNameVal is None:
                stepNameVal = ""

            # determine currentOdbStep from odbStepVal
            try:
                currentOdbStep = int(odbStepVal)
            except (ValueError, TypeError):
                if oldOdbStep==0:
                    if rowCnt>0:
                        raise ValueError(
                            "Improper value for currentOdbStep %d when"
                            " processing %s." % (currentOdbStep, errorLocStr))
                    if totaltime>0.0:
                        currentOdbStep = 2
                    else:
                        currentOdbStep = 1
                    msg("WARNING: No proper odb step number given on the first"
                        " data line of %s. Assuming %d as first odb step nb."
                        % (loadhistoryName, currentOdbStep))
                elif (odbStepVal is None) or (odbStepVal==""):
                    # if not specified then just leave as before
                    pass
                else:
                    # if it cannot be interpreted then leave as before
                    msg("WARNING: odb step number could not be properly"
                        " evaluated in %s. Continuing the step %d."
                        % (errorLocStr, currentOdbStep))

            # check odb step nb did not decrease
            if currentOdbStep<oldOdbStep:
                msg("WARNING: odb step number decreased from %d to %d in %s."
                    " This will be ignored and the old step %d will be cont'd."
                    % (oldOdbStep, currentOdbStep, errorLocStr, currentOdbStep))
                currentOdbStep = oldOdbStep

            # store totaltime and steptime
            totaltimeList.append(totaltime)
            if currentOdbStep>oldOdbStep:
                steptimeList.append(0.0)            # odb step change
            else:
                steptimeList.append(steptime)

            # if odb step nb changed:
            #  ... store current step in self.odbStepLengthDict
            #  ... reset steptime and odbFrame
            if currentOdbStep>oldOdbStep and oldOdbStep>0:
                self.odbStepLengthDict[oldOdbStep] = steptime
                self.odbStepLastFrameDict[oldOdbStep] = odbFrameList[-1]
                steptime = 0.0
                odbFrame = 1

            # store odb step and frame nb and steplength and stepname
            odbStepList.append(currentOdbStep)
            odbFrameList.append(odbFrame)
            stepLengthList.append(stepLengthVal)
            stepNameList.append(stepNameVal)

            # update steptime and totaltime for next iteration
            steptime += stepLengthVal
            totaltime += stepLengthVal

            # increase odbFrame counter for next iteration
            if isOdbFrameFlag:
                odbFrame += 1

            # store current as old for next iteration
            oldOdbStep = currentOdbStep


        # store last step processed in self.odbStepLengthDict
        self.odbStepLengthDict[oldOdbStep] = steptime
        self.odbStepLastFrameDict[oldOdbStep] = odbFrameList[-1]

        #-- store odbStep, odbFrame, steptime, totaltime, stepname, steplength
        #   in self.loadHistTable
        for colKey, data in [
            ("odbStep", odbStepList),
            ("odbFrame", odbFrameList),
            ("steptime", steptimeList),
            ("totaltime", totaltimeList),
            ("miningStepLength", stepLengthList),
            ("miningStepName", stepNameList)
        ]:

            # rename conflicting columns
            bakColKey = colKey
            while bakColKey in self.loadHistTable:
                bakColKey = "orig_%s" % bakColKey
            if bakColKey != colKey:
                self.loadHistTable.renameCol(colKey, bakColKey)

            # insert new column
            self.loadHistTable.insertCol(colKey, data)

        # store some data permanently
        self.totalAnalysisTime = totaltime

    @staticmethod
    def filterDateAndColumn_Type_Date(elsetName):
        """Example function for the filterDateAndColumn argument of
        L{LoadHistory.updateElsetNames}. Takes a single elset name as argument
        and returns a (type, date) tuple. For elset names not to be considered
        returns None.

        This is a rather simple example and you might want to supply a more
        sophisticated version for your actual problem.

        Tries to locate the string "_Y" (assuming the "Y" meaning "year") and
        then returns everything before as type and and everything after and
        including the "Y" as date.

        @param elsetName: ... guess
        @returns: (type, date) tuple. The type part will serve as column name
           and the date part will be looked up in the miningStepName column to
           identify the correct mining step.
           For elset names not to be considered returns None.
        """
        type_, date = elsetName.rsplit("_Y", 1)
        if not date:
            return None
        return (type_, "Y"+date)

    def updateElsetNames(self, elsetNames, filterDateAndColumn=None):
        """Add given elsets to the load history table. ElsetNames are expected
        to be self-explanatory in the sense that they carry the type /
        elset-column and date informtation coded into themselves. A function is
        provided by the argument filterDateAndColumn that extracts type and
        date from the elset name.

        This can be used with an almost empty loadHistory csv file.
        (For efficiency this example does not use
        L{LoadHistory.updateDisjointElsets} but reads the sequence elsets
        "manually" only once for updateElsetNames and disjointSeqElsets.)
         >>> loadhistory = LoadHistory(
         >>>     loadhistoryCsv=loadhistoryCsvName,
         >>>     elsetColumns=set(x[0] for x in colToGravAmp+colToSDVini),
         >>>     odbStepNb=2,
         >>>     miningStepLength="stepLength", miningStepName='Step_Name',
         >>>     isOdbFrame="Odb_Frame", startTotalTime=stepTimeEq,
         >>>     extraLocals=locals())
         >>>
         >>> # read the sequence elset definition files
         >>> setsSequenceModel = Model()
         >>> for setSequenceFile in setsSequenceFnames:
         >>>     setsSequenceModel.read(
         >>>         setSequenceFile, recognizedBlocks="ELSET")
         >>>
         >>> # add sequence elsets to loadhistory
         >>> loadhistory.updateElsetNames(sorted(setsSequenceModel.elset))
         >>>
         >>> # create disjoint elsets (boolElsets) and
         >>> # ...  remove empty elsets from the load history
         >>> disjointSets = loadhistory.disjointSeqElsets(
         >>>     elsetsDict=setsSequenceModel.elset,
         >>>     updateLoadHistory=True,
         >>>     additionalElsets=additionalElsetNames,
         >>>     reportWrongElsetNb=200)
         >>>
         >>> # store disjoint elsets in Abaqus input file
         >>> disjointSetsModel = Model()
         >>> disjointSetsModel.elset.update(disjointSets)
         >>> disjointSetsModel.write(setsDisjointSeqPath)

        @param elsetNames: An iterable of elsetNames.
        @param filterDateAndColumn: A function taking a single elset name as
           argument and returning a (type, date) tuple. The type part will
           serve as column name and the date part will be looked up in the
           miningStepName column to identify the correct mining step.
           For elset names to be discarded this function returns None.

           If this function is not provided it defaults to
           L{filterDateAndColumn_Type_Date} which assumes the common pattern
           like "PIT_Y2017_M01".
        """
        # default argument
        if filterDateAndColumn is None:
            filterDateAndColumn=self.filterDateAndColumn_Type_Date

        # colElsetNames: {column : {date : [list of elsets] } }
        class DefaultDictListDict(defaultdict):
            @staticmethod
            def _newListDict():
                return defaultdict(list)
            def __init__(self):
                defaultdict.__init__(self)
                self.default_factory = self._newListDict
        colElsetNames = DefaultDictListDict()

        # initialize colElsetNames with the elsetNames argument
        newSetsNumber = 0
        ignoredElsetNames = list()
        for elsetName in elsetNames:
            typeDate = filterDateAndColumn(elsetName)
            if not typeDate:
                ignoredElsetNames.append(elsetName)
                continue
            colElsetNames[typeDate[0].upper()][typeDate[1].upper()].append(
                elsetName.upper())
            newSetsNumber += 1

        msg("Updating loadhistory with %d new elsets of %d types/columns..."
            % (newSetsNumber, len(colElsetNames)))
        if ignoredElsetNames:
            # print a warning
            # some heuristics to group elsets in ignoredElsetNames
            typedElsets = defaultdict(list)
            for e in ignoredElsetNames:
                typedElsets[e.split("_", 1)[0]].append(e)
            cnt = int(100 / len(typedElsets))+1
            typedElsets = sorted(typedElsets.iteritems())
            elsetNames = [
                x for t, xx in typedElsets for x in xx[:cnt]]
            msg("LoadHistory.updateElsetNames(): %d elsets from the list of"
                " elsets to be added didn't match the elset name pattern"
                " and have been ignored. Examples: %s"
                % (len(ignoredElsetNames), elsetNames[:100]))
            del cnt, typedElsets, elsetNames, ignoredElsetNames
        # newData = DictTableFromCsvFile()  ... would be nice:
        # ... to store added columns in a DictTableFromCsvFile to (optionally
        #     write this to a csv file.
        #     currently DictTableFromCsvFile needs some data in the csv on
        #     initialization.
        for colKey, dateElsets in colElsetNames.iteritems():
            if colKey not in self.loadHistTable:
                # a new column must also be added to the set of elsetColumns
                self.elsetColumns.add(colKey)
            self.loadHistTable.insertCol(
                colKey, [
                    dateElsets.pop(date, [])
                    for date in self.loadHistTable["miningStepName"]])
            if dateElsets:
                # print a warning
                # some heuristics to group elsets in dateElsets
                cnt = int(100 / len(dateElsets))+1
                elsetNames = [x for t, xx in dateElsets.iteritems()
                              for x in xx[:cnt]]
                cnt = sum(len(x) for x in dateElsets.itervalues())
                msg("WARNING: LoadHistory.updateElsetNames() did not find a"
                    " matching date in the loadhistory table for %d elsets."
                    " Examples: %s" % (cnt, elsetNames[:100]))
                del cnt, elsetNames


    def getLoadHistTableIter(self):
        r"""Convenience function, equivalent to
        self.loadHistTable.rowIterColToList(listColumns=self.elsetColumns,
        removeNones=True)

        @Returns: an iterator over all rows in the load history table.
        All elset columns (self.elsetColumns) are converted to lists with just
        the applicable elset names (no Nones for empty columns).
        """
        return self.loadHistTable.rowIterColToList(
            listColumns=self.elsetColumns, removeNones=True)

    def getLoadHistRow(self, rowIdx):
        r"""Convenience function, equivalent to
        self.loadHistTable.getRowColToList(rowIdx=rowIdx,
        listColumns=self.elsetColumns, removeNones=True)

        @Returns: the specified row of the loadHistoryTable.
        All elset columns (self.elsetColumns) are converted to lists with just
        the applicable elset names (no Nones for empty columns).

        @Note: rowIdx is a zero based index! No matter how your steps are
        counted.
        """
        return self.loadHistTable.getRowColToList(
            rowIdx=rowIdx, listColumns=self.elsetColumns, removeNones=True)
    #} end of basic functions


    #------------------------------------------------------------------
    #{ dynamic value evaluation

    @staticmethod
    def evaluateString(value, evalNamespaceDict, description):
        r"""
        If value is a string it is tried to evaluate as a python
        expression with the evalNamespaceDict as local name space.

        >>> print self.evaluateString("time*0.5", {"time":10}, "a test")
        5.0

        This method is intended to be used to evaluate an expression in the
        context of the current row of the load history table. I.e. like that:
        (extraLocals is a dict like locals(), value is the expression to be
        evaluated, self is the LoadHistory object)
         >>> for row in self.getLoadHistTableIter():
         >>>     evalNamespaceDict = extraLocals.copy()
         >>>     evalNamespaceDict.update(row)
         >>>     val = self.evaluateString(
         >>>         value, evalNamespaceDict,
         >>>         "some error string for diagnostic purposes")

        See the note on L{value evaluation<LoadHistory>}.

        @param description: description of the value for error messages.
        """
        if isinstance(value, basestring):
            try:
                value = eval(value,globals(), evalNamespaceDict)
            except Exception as exc:
                raise ValueError(
                    "Could not evaluate the expression '%s' for %s:\n%s"
                    % (value, description, exc.args[0]))
        return value

    @staticmethod
    def evaluateStringAsFlag(value, evalNamespaceDict, description):
        r"""
        If value is a string it is tried to evaluate as a python
        expression with the evalNamespaceDict as local name space.

        >>> print self.evaluateStringAsFlag("flag", "flag":True}, "a test")
        True
        >>> print self.evaluateStringAsFlag(False, "flag":True}, "a test")
        False
        >>> print self.evaluateStringAsFlag("flag", "flag":"On"}, "a test")
        True

        See the note on L{value evaluation<LoadHistory>}.

        @param description: description of the value for error messages.
        """
        if isinstance(value, basestring):
            try:
                value = eval(value, globals(), evalNamespaceDict)
            except Exception as exc:
                raise ValueError(
                    "Could not evaluate the expression '%s' for %s:\n%s"
                    % (value, description, exc.args[0]))

        value = (value is True
                 or (isinstance(value,int) and value!=0)
                 or (isinstance(value, float)
                     and not (-1E-6 < value < 1E-6))
                 or (isinstance(value, basestring)
                     and value.lower() in ("yes", "on", str(1)))
                 )

        return value
    #} end of dynamic value evaluation


    #------------------------------------------------------------------
    #{ disjoint elsets

    def removeElsetFromSteps(self, elsetName, stepIds=None, elsetColumns=None):
        """remove the given elset from the steps listed (by index) in
        the iterable stepIds.

        This is basically a service function for internal use.

        @param stepIds: iterable of step indexes. If None then remove elset
          from all steps. (This may be slow.)

        @param elsetColumns: list or set of columns labels specifying all
          columns in the loadhistoryCsv-file containing all elsets that should
          be considered. If not specified then defaults to self.elsetColumns.
          Case insensitive.
        """
        if elsetColumns is None:
            elsetColumns = self.elsetColumns
        else:
            # if argument passed in then convert to upper case
            elsetColumns = [x.upper() for x in elsetColumns]

        if stepIds is None:
            stepIds = range(self.loadHistTable.getColumnLength())

        # find this empty/not existing elset in self.loadHistTable
        indicesToRemove = list()
        for rowCnt in stepIds:
            indicesToRemove.extend(
                (rowCnt, colKey, vIdx)
                for colKey, vIdx in self.loadHistTable.findInRow(
                    elsetName, rowCnt, inColumns=elsetColumns,
                    singleColValIndex=None))
        # remove this empty/not existing elset from the steps
        for rowCnt, colKey, vIdx in indicesToRemove:
            if vIdx is None:
                self.loadHistTable[colKey][rowCnt] = None
            else:
                self.loadHistTable[colKey][rowCnt][vIdx] = None

    def removeElsetFromStep(self, elsetName, stepId, elsetColumns):
        """remove the given elset from the step given by its index.

        This is basically a service function for internal use.

        @param elsetColumns: list or set of columns labels specifying all
          columns in the loadhistoryCsv-file containing all elsets that should
          be considered. Case insensitive.
        """
        # find this empty/not existing elset in self.loadHistTable
        elsetColumns = [x.upper() for x in elsetColumns]  # convert to up case
        indicesToRemove = self.loadHistTable.findInRow(
            elsetName, stepId, inColumns=elsetColumns,
            singleColValIndex=None)
        # remove this empty/not existing elset from the step
        for colKey, vIdx in indicesToRemove:
            if vIdx is None:
                self.loadHistTable[colKey][stepId] = None
            else:
                self.loadHistTable[colKey][stepId][vIdx] = None

    def disjointSeqElsets(self, elsetsDict, updateLoadHistory=True,
                          additionalElsets=[], elsetColumns=None,
                          checkOnlyColumns=None,
                          reportWrongElsetNb=20, warnNotUsedElsets=True):
        r"""
        Does what the boolElsets script did: Modify the given elsets in such a
        way that they are disjoint according to the order those elsets appear
        in self.loadHistTable and after that in additionalElsets. Elsets that
        don't appear in the load history (or in the additionalElsets list) or
        that are empty after this procedure are not considered in the result.
        At the same time it removes empty elsets from the load history steps
        if the argument updateLoadHistory is set to True.

        @param elsetsDict: A dictionary like abq_model_02.Model.elset. This is
          not touched. Instead a modified copy is returned

        @param updateLoadHistory: If True remove all elsets not found in the
          resulting elsets dictionary from the self.loadHistTable items. I.e.
          elsets that are empty after this process or have not been there from
          the beginning won't be assigned a gravity amplitude anymore nor will
          they appear in the initial conditions lists.

        @param additionalElsets: List of elset names that appear in elsetsDict
          and are supposed to appear in the output (if not empty). They are
          being made disjoint from all others as well. This list acts as if
          there where additional steps at the end of the loadhistory, one for
          each elset in the order of this list. (This should be used for the
          NOTMINED / GRAV-01 elsets.)

        @param elsetColumns: DEPRECATED argument. If present it's being
          compared to self.L{elsetColumns} and raises a ValueError if the
          argument contains an item that is not present in self.L{elsetColumns}.

          Was intended to specify all columns in the loadhistoryCsv-file
          containing all elsets that should be considered. This is now the
          purpose of self.L{elsetColumns} instance variable which is set by
          the elsetColumns argument to the constructor.

        @param checkOnlyColumns: optional list of elset column names to be
          considered for this operation. If not specified then take
          self.elsetColumns. Needed for the split elset feature.

        @param reportWrongElsetNb: how many elsets should be displayed that
          don't appear in the loadhistory table but are defined in the input
          file.

        @param warnNotUsedElsets: If True (the default) warn about elsets that
          don't appear in the load history table.
          This is for internal use only. And currently it's recommended to be
          set to False when using the "split elsets feature" because otherwise
          you'd get meaningless warnings. May be removed at a later stage.

        @returns: A modified copy of the input elsetsDict with all the elsets
          made disjoint according to the order those elsets appear in
          self.loadHistTable and after that in additionalElsets.

        @Note: The order of elsets that occur in the same step for the first
          time is random. (Maybe this is not true. I now think the order is
          determined first by elsetColumns (or self.elsetColumns) and second
          --for elsets in columns of the same header-field/name-- by the order
          they appear in the csv file from left to right.)

        @Note: Someone has to explain why this function is so complicated where
          what it does is so simple. Thinks Gero.
        """

        # deprecated functionality!
        if elsetColumns is not None:
            missing = set(c.upper() for c in elsetColumns) \
                .difference(self.elsetColumns)
            if missing:
                raise ValueError(
                    "LoadHistory.disjointSeqElsets: Elset column(s) %s not"
                    " recognized as elsetColums (as given by the corresponding"
                    " argument to LoadHistory.__init__). Please update the"
                    " elsetColumns argument."
                    % ",".join(sorted(missing)))

        # determine elset columns to be considered
        if checkOnlyColumns is None:
            checkOnlyColumns = self.elsetColumns
        else:
            checkOnlyColumns = [x.upper() for x in checkOnlyColumns]

        # convert elsetsDict to a ElsetsDict in order to ignore upper/lower case
        elsetsDict = ElsetsDict(elsetsDict)

        # elsetNamesToSteps is a dict {elset name: set of step indexes}. A step
        # index is the index in self.loadHistTable of the step this elset occurs
        # in any of the columns listed in elsetColumns
        elsetNamesToSteps = defaultdict(list)
        # elsetsOrder: elset names in the order they appear in
        # self.loadHistTable
        elsetsOrder = list()

        allElsetNames = set()  # all elset names in a set

        rowIterator = self.loadHistTable.rowIterColToList(
            listColumns=checkOnlyColumns, removeNones=True)
        for rowCnt, row in enumerate(rowIterator):
            # sum up all lists of lists of elsets and convert to uppercase
            # find and remove duplicates
            elsetsOfThisStep = dict()
            for col in checkOnlyColumns:
                for idx, elsetName in enumerate(row[col]):
                    elsetName = elsetName.upper()
                    if elsetName in elsetsOfThisStep:
                        msg("WARNING: On line %d of the loadhistory table the"
                            " elset %s occurs twice. All but the first"
                            " occurence will be ignored."
                            % ((rowCnt+1), elsetName))
                        if updateLoadHistory:
                            self.loadHistTable[col][rowCnt][idx] = None
                    else:
                        elsetsOfThisStep[elsetName] = len(elsetsOfThisStep)

            # sorted list (sorted because order matters if elsets intersect)
            elsetsOfThisStep = [x[1] for x in sorted(
                (i,y) for y, i in elsetsOfThisStep.iteritems())]

            # add all elsets to elsetNamesToSteps {elset: list of steps}
            for elsetName in elsetsOfThisStep:
                elsetNamesToSteps[elsetName].append(rowCnt)

            # check for elsets already mined earlier
            doubleElsets = allElsetNames.intersection(elsetsOfThisStep)
            if len(doubleElsets):
                msg("WARNING: On line %d of the loadhistory table the"
                    " following elset(s) occur for the second time in this"
                    " sequence: %s. All but the first occurence will be"
                    " ignored."
                    % ((rowCnt+1), sorted(doubleElsets)[:50]))
                elsetsOfThisStep = [elset for elset in elsetsOfThisStep
                                    if elset not in doubleElsets]
                if updateLoadHistory:
                    for elsetName in doubleElsets:
                        self.removeElsetFromStep(
                            elsetName, rowCnt, checkOnlyColumns)

            # add to the ordered list for making them disjoint (only once)
            elsetsOrder.extend(elsetsOfThisStep)
            allElsetNames.update(elsetsOfThisStep)

        elsetsOrder.extend(additionalElsets)  # order of addi. preserved so far
        additionalElsets = set(additionalElsets)  # convert for speed

        # create:
        #  - disjointElsets: dictionary of disjoint elsets
        #  - allElements: set of all elements in any of the elsets
        disjointElsets = dict()
        allElements = set()
        emptyElsets = list()
        for elsetName in elsetsOrder:
            removeThisElset = False
            newElset = []
            try:
                newElset = elsetsDict[elsetName].difference(allElements)
            except KeyError:
                if elsetName in additionalElsets:
                    # elsetName from additionalElsets, issue warning then
                    # ignore this and go on with the next.
                    msg("WARNING: 'Additional' elset %s not found. Will be"
                        " ignored." % elsetName)
                    continue
                # elsetName from the loadhistory, issue warning and remove from
                # load history
                msg("WARNING: Elset %s listed in the loadhistory table"
                    " not found. Will be ignored." % elsetName)

                removeThisElset = True
            else:
                # elset found in elsetsDict, already disjoint to all processed
                # so far
                if len(newElset):
                    allElements.update(newElset)
                    disjointElsets[elsetName] = newElset
                else:
                    emptyElsets.append(elsetName)
                    removeThisElset = True

            # if elset became empty or has not been in elsetsDict initially
            # then remove from the load history
            if (removeThisElset and updateLoadHistory
                and elsetName in elsetNamesToSteps):
                for stepId in elsetNamesToSteps[elsetName]:
                    self.removeElsetFromStep(
                        elsetName, stepId, checkOnlyColumns)

        if len(emptyElsets):
            if len(emptyElsets)>20:
                elsetsText = "%s..." % (", ".join(emptyElsets[:20]))
            else:
                elsetsText = ", ".join(emptyElsets)
            msg("WARNING: %d elsets have been removed because"
                " they are empty: %s"
                % (len(emptyElsets), elsetsText))

        # element sets that don't appear in the load history
        notUsedElsets = set(elsetsDict).difference(emptyElsets) \
            .difference(disjointElsets)
        if len(notUsedElsets) and warnNotUsedElsets:
            if len(notUsedElsets)>reportWrongElsetNb:
                elsetsText = ("%s... and %d more ..." % (
                    ", ".join(sorted(notUsedElsets)[:reportWrongElsetNb]),
                    len(notUsedElsets)-reportWrongElsetNb ))
            else:
                elsetsText = ", ".join(sorted(notUsedElsets))
            msg("WARNING: %d elsets do not appear in the load history. They"
                " won't appear as one of the resulting disjoint elsets either:"
                " %s" % (len(notUsedElsets), elsetsText))

        return disjointElsets

    def intersectColumns(self, elsetsDict, elsetColumns,
                         updateLoadHistory=True):
        r"""
        Intersect two columns of the loadhistory table for the split elsets
        feature.

        Appart from the split elsets as directly returned result (see below)
        this function also creates two dictionaries with data for later
        processing:
         - self.stepColToSplitElsets for L{writePostDataFile} is a dict
           {(row index, elset column name): list of names of intersection
           elsets}. E.g. { (7, "stopes"): ["STOPES_Y2012_FILL_Y2014",
           "STOPES_Y2012_FILL_Y2015"], (9, "fill"):
           ["STOPES_Y2012_FILL_Y2014",], (10, "fill"):
           ["STOPES_Y2012_FILL_Y2015",], ... }
         - self.colsStep1ToStep2Elsets for L{getGravAmpDload} and
           L{getStrSdvIni} is a dict {((pair of elset columns), row index for
           first elset column): list of (row index for second elset column,
           combined name of the intersection-elset)}
           E.g. {(("stope","fill"), 7): [(9, "STOPES_Y2012_FILL_Y2014"),
           (10, "STOPES_Y2012_FILL_Y2015")], ...}

        @param elsetsDict: A dictionary like abq_model_02.Model.elset. This is
          not touched. Instead a modified copy is returned

        @param elsetColumns: tuple of two elset column names to be
          considered for this operation. This tuple (including the order of
          the items) must be identical to the corresponding items in the
          colToGravAmp-argument of L{getGravAmpDload} and in the colToSDVini-
          argument of L{getStrSdvIni}. Those elset column names are treated in
          a case insensitive manner.

        @param updateLoadHistory: If True remove all elsets that became empty
          from self.loadHistTable. I.e. elsets that are empty after this
          process or have not been there from the beginning won't be assigned
          a gravity amplitude anymore nor will they appear in the initial
          conditions lists.

        @returns: A modified copy of the input elsetsDict with new intersection
          elsets according to self.loadHistTable.
        """

        elsetColumns = tuple(x.upper() for x in elsetColumns)

        if not hasattr(self, "stepColToSplitElsets"):
            self.stepColToSplitElsets = defaultdict(list)
        else:
            self.stepColToSplitElsets.default_factory = list
        if not hasattr(self, "colsStep1ToStep2Elsets"):
            self.colsStep1ToStep2Elsets = defaultdict(list)
        else:
            self.colsStep1ToStep2Elsets.default_factory = list
        splitElsetDict = dict()
        elsetsSplit = set()

        # iteration over mining steps
        for rowCnt, row in enumerate(self.getLoadHistTableIter()):

            # get elset of first column (all elsets combined)
            elsetListFirstCol = self._getElsetListFromColumn(
                elsetColumns[0], row, "intersectColumns")
            if not(elsetListFirstCol):
                continue
            firstElset = set().union(*(elsetsDict[e]
                                       for e in elsetListFirstCol))

            # iterate over all steps for the second column
            for secRowCnt, secRow in enumerate(
                    self.getLoadHistTableIter()):
                # determine intersection of elset of first column with elsets
                # of second column
                elsetListSecondCol = self._getElsetListFromColumn(
                    elsetColumns[1], secRow, "intersectColumns")
                if not(elsetListSecondCol):
                    continue
                secElset = set().union(*(elsetsDict[e]
                                         for e in elsetListSecondCol))
                combinedElset = firstElset.intersection(secElset)
                if not combinedElset:
                    continue

                # store new elset
                combinedName = "_".join((
                    elsetColumns[0].upper(), row["miningStepName"],
                    elsetColumns[1].upper(), secRow["miningStepName"]))
                splitElsetDict[combinedName] = combinedElset

                # store split elset names for later use (a,b,c)
                # (a) for combining old and new elsets
                elsetsSplit.update(elsetListFirstCol)
                elsetsSplit.update(elsetListSecondCol)

                # (b) for postData
                self.stepColToSplitElsets[
                    (rowCnt, elsetColumns[0])].append(combinedName)
                self.stepColToSplitElsets[
                    (secRowCnt, elsetColumns[1])].append(combinedName)

                # (c) for getGravAmpDload and getStrSdvIni
                self.colsStep1ToStep2Elsets[
                    (tuple(elsetColumns), rowCnt)].append(
                    (secRowCnt, combinedName))

        # convert from defaultdict(list) to ordinary dict
        self.stepColToSplitElsets.default_factory = None
        self.colsStep1ToStep2Elsets.default_factory = None

        # combine old and new elsets in result
        allNewElems = set.union(*(splitElsetDict.itervalues()))
        msg("Created %d new split elsets with %d elements."
            % (len(splitElsetDict), len(allNewElems)))

        # print "elsetsSplit:", elsetsSplit  #######
        # print "splitElsetDict:", splitElsetDict  #######
        # print "allNewElems:", allNewElems  #######
        # print ":",   #######
        emptysets = list()
        for eln, els in elsetsDict.iteritems():
            if eln in elsetsSplit:
                croppedEls = els.difference(allNewElems)
                if croppedEls:
                    # store remains in splitElsetDict
                    splitElsetDict[eln] = croppedEls
                else:
                    emptysets.append(eln)
            else:
                # just copy elset not involved in splitting
                splitElsetDict[eln] = els

        if updateLoadHistory:
            for eln in emptysets:
                self.removeElsetFromSteps(eln)
        return splitElsetDict

    def updateDisjointElsets(self, setsSequenceFileName,
                             disjointSetsFileName=None, additionalElsets=[],
                             checkOnlyColumnsIter=None,
                             reportWrongElsetNb=20):
        """
        Removes empty elsets from the load history and (re)creates the Abaqus
        input file that contains the disjoint elsets if it appears to be
        necessary.

        This is basically a wrapper around self.disjointSeqElsets() that reads
        and writes the necessary files.

        Note: Generally don't use this method for the split elset functionality
        because after making the elsets disjoint you still have to split them
        which is done by self.L{getStrSdvIni}().

        @param setsSequenceFileName: file name or list of file names of the
          Abaqus input file(s) containing the elsets before they have been made
          disjoint. Instead of file names you can also provide open file
          objects.

        @param disjointSetsFileName: file name of the Abaqus input file
          to be written on output containing the modified elsets.
          If None (default) then elsets are not being written. (Just removes
          empty elsets from the load history in this case.)

        @param additionalElsets: May be a list or a dictionary of additional
          elsets not in the load history steps that shall be made disjoint from
          the rest as if they appeared at the end of the load history. (This
          should be used for the NOTMINED --formerly called GRAV-01-- elsets.)

          If it is a list of elset names then those elsets must appear in the
          Abaqus input file mentioned as setsSequenceFileName.

          If the elsets don't appear in this Abaqus input file, you have to
          supply a dictionary {elset name: set of elements} e.g. like
          abq_model_02.Model.elset.

        @param checkOnlyColumnsIter: optional. An iterator of lists
          of elset column names to be considered for subsequent calls to
          self.disjointSeqElsets(). Needed for the split elset feature.

          If not specified then call self.disjointSeqElsets() only once with
          checkOnlyColumns=self.elsetColumns.

          Note that the parameter additionalElsets is only considered for the
          first item of checkOnlyColumnsIter. I.e. The NOTMINED elset will
          be shaped by the elsets in columns listed in the first item of this
          iterator. Following items only make sure that the secondary elsets
          used for elset splitting are free of overlaps in themselves.

        @param reportWrongElsetNb: how many elsets should be displayed that
          don't appear in the loadhistory table but are defined in the input
          file.
        """

        if type(setsSequenceFileName)==list:
            setsSequenceModel = Model()
            for setSequenceFile in setsSequenceFileName:
                setsSequenceModel.read(setSequenceFile)
        else:
            setsSequenceModel = Model().read(setsSequenceFileName)

        if isinstance(additionalElsets, dict):
            additionalElsetNames = additionalElsets.keys()
            additionalElsetNames.sort()
        elif isinstance(additionalElsets, list):
            additionalElsetNames = additionalElsets
            additionalElsets = dict()
        else:
            raise ValueError(
                "The additionalElsets argument must be either a dictionary or"
                " a list. (It's a %s.)" % type(additionalElsets))

        if checkOnlyColumnsIter is None:
            oldSeqElsets = dict(setsSequenceModel.elset)
            oldSeqElsets.update(additionalElsets)
            disjointSets = self.disjointSeqElsets(
                elsetsDict=oldSeqElsets,
                updateLoadHistory=True,
                additionalElsets=additionalElsetNames,
                reportWrongElsetNb=reportWrongElsetNb)
        else:
            disjointSets = dict()
            firstIteration = True
            for checkOnlyColumns in checkOnlyColumnsIter:
                if firstIteration:
                    oldSeqElsets = dict(setsSequenceModel.elset)
                    oldSeqElsets.update(additionalElsets)

                disjointSets.update(self.disjointSeqElsets(
                    elsetsDict=oldSeqElsets,
                    updateLoadHistory=True,
                    additionalElsets=additionalElsetNames,
                    checkOnlyColumns=checkOnlyColumns,
                    reportWrongElsetNb=reportWrongElsetNb,
                    warnNotUsedElsets=False))

                if firstIteration:
                    # for next iteration only take elsets from the model
                    oldSeqElsets = dict(setsSequenceModel.elset)
                    additionalElsetNames=[]
                    firstIteration = False

        if disjointSetsFileName is not None:
            disjointSetsModel = Model()
            disjointSetsModel.elset.update(disjointSets)
            disjointSetsModel.write(disjointSetsFileName)
    #} end of disjoint elsets


    #------------------------------------------------------------------
    #{ Abaqus input file commands

    def _getElsetListFromColumn(self, elsetColumn, row, caller):
        """Service function: get list of elset names from the given elsetColumn
        in the current row of the loadhistory.
        Caller is the method name of the calling function for dignostic output.

        To be called from L{getGravAmpDload} and L{getStrSdvIni}

        @param elsetColumn: must be uppercase!
        """
        if elsetColumn not in self.elsetColumns:
            raise ValueError(
                "LoadHistory.%s: Elset column %s not"
                " recognized as elsetColum (as given by the"
                " corresponding argument to LoadHistory.__init__)."
                " Please update the elsetColumns argument."
                % (caller, elsetColumn))

        try:
            elsetList = row[elsetColumn]
        except KeyError:
            raise KeyError(
                "LoadHistory.%s: Elset column %s not"
                " found in the loadhistory table. Please check."
                % (caller, elsetColumn))

        return elsetList

    def getGravAmpDload(self, colToGravAmp, extraLocals={},
                        odbStep=None, time="'STEP'"):
        r"""
        Return an (amplitude string, dload string) tuple.

        The amplitude string contains the *amplitude commands and the dload
        string contains the *dload command needed to apply the gravity loads
        in an abaqus input file.

        The amplitudes are called "GRAVAMP_%s" % elsetName

        Examples:
        =========

        Same timing for all elsets in the "Pit" column and another for
        stopes. Note: the column "steptime" does not have to be defined in the
        load history csv file, it is generated automatically, see L{__init__}:
         >>> self.getGravAmpDload([
         ...     ("Pit",    [("steptime",     1.0),
         ...                 ("steptime+0.5", 0.0)]),
         ...     ("Stopes", [("steptime",     1.0),
         ...                 ("steptime+0.5", 0.0),
         ...                 ("steptime+1.5", 0.0),
         ...                 ("steptime+1.6", 0.8)]),
         ...     ])

        Variable timing for Pit, being defined in the column T_EndExc. And
        timing relative to the mining step length for Stopes. Note: the columns
        "steptime" and "miningStepLength" do not have to be defined in the
        load history csv file, they are generated automatically,
        see L{__init__}:
         >>> self.getGravAmpDload([
         ...     ("Pit",    [("steptime", 1.0),
         ...                 ("steptime+T_EndExc", 0.0)]),
         ...     ("Stopes", [("steptime", 1.0),
         ...                 ("steptime+0.25*miningStepLength", 0.0),
         ...                 ("steptime+1.5", 0.0),
         ...                 ("steptime+1.6", 0.8)]),
         ...     ])

        Variable amplitude value according to the column EndGrav in the load
        history csv file and a variable drivesFact defined locally:
         >>> drivesFact = 0.1
         >>> self.getGravAmpDload(
         ...     [ ("Stopes", [("steptime", 1.0),
         ...                   ("steptime+0.5", 0.0),
         ...                   ("steptime+1.5", 0.0),
         ...                   ("steptime+1.6", "EndGrav")]),
         ...       ("Pit",    [("steptime", 1.0),
         ...                   ("steptime+T_EndExc", "drivesFact*endGrav")]),
         ...     ],
         ...     extraLocals=locals()
         ...     )

        If you split the mining history over many odb steps and want to define
        the gravity dload with the amplitudes just in the first of those
        (probably Step-2) you'd have to use TOTAL TIME - amplitudes.
        Don't forget to switch from steptime to totaltime in the times of the
        colToGravAmp argument!
         >>> self.getGravAmpDload(
         ...     colToGravAmp=[
         ...         ("Stopes", [("totaltime", 1.0),
         ...                     ("totaltime+0.5", 0.0),
         ...                     ("totaltime+1.5", 0.0),
         ...                     ("totaltime+1.6", "EndGrav")]),
         ...         ],
         ...     time="'TOTAL'",
         ...     )

        @param colToGravAmp: A list of (column label, gravity amplitude list)
          tuples.

          Each column label identifies the column(s) under which the elsets are
          listed to which this timing is supposed to apply to. E.g. 'Pit',
          'Cave', 'Fill'. Case insensitive.

          For the "split elsets feature" a tuple of two elset names is given
          as column label to identify the corresponding intersection elset.
          This tuple (including the order of the items) must be identical to
          that used in the elsetColumns argument of the preceeding call to
          L{intersectColumns}. Otherwise this item will not be recognized.

          The gravity amplitude list is a list of (time, amplitude) tuples
          defining the gravity amplitude.

          Each of the time and amplitude values may be floats, in that case
          this value applies to all loadhistory steps. Or it may be a string
          evaluating to the corresponding value.
          See the note on L{value evaluation<LoadHistory>}.

          The gravity amplitude list might also be None. In this case the
          corresponding elset column is to be ignored for this purpose.
          A warning message about missing column labels is being suppressed in
          this case.

          Note: All column labels are converted to upper case. It is essential
          to have had all column labels in the elsetColumns argument of
          self.__init__. Otherwise this will raise a ValueError.

        @param extraLocals: A namespace dictionary that contains additional
          variables that should be accessible during value evaluation, see
          L{value evaluation<LoadHistory>}. Specify extraLocals=locals()
          in the argument list to make available all current variables.

        @param odbStep: return amplitudes and dloads for the given odb step
          number only. If not specified take all steps in the loadhistory.

        @param time: Defines the time parameter of the *AMPLITUDE command.
          Might be "'STEP'" or "'TOTAL'" or anything that evaluates to one of
          those. See the note on L{value evaluation<LoadHistory>}.
          And remember to put extra quotes on literal strings, otherwise the
          scripts looks for a TOTAL column in the csv file for instance.

          WARNING: This does not affect the time values in the amplitude
          definition! Most certainly you'll have to adapt the colToGravAmp
          parameter as well if you change from STEP TIME to TOTAL TIME!
        """
        # check parameter
        if not isinstance(colToGravAmp, list):
            raise ValueError(
                "LoadHistory.getGravAmpDload(): expect list as colToGravAmp"
                " argument, got %s." % type(colToGravAmp))

        # further check: all elsetColumns covered in colToGravAmp?
        elsetColumnsInPara = set()
        for elsetColumn, amplist in colToGravAmp:
            if isinstance(elsetColumn, basestring):
                elsetColumnsInPara.add(elsetColumn.upper())
            elif isinstance(elsetColumn, (list, tuple)):
                elsetColumnsInPara.update(e.upper() for e in elsetColumn)
        missing = self.elsetColumns.difference(elsetColumnsInPara)
        if missing:
            msg("WARNING: The following elset columns are not mentioned in the"
                " colToGravAmp parameter to LoadHistory.getGravAmpDload(): %s"
                "\nIf this is intended and you want to avoid this warning then"
                " specify those columns with None instead of the gravity"
                " amplitude list."
                % ",".join(sorted(missing)))
            # Note: The reverse condition (colToGravAmp mentions column that
            # is not in self.elsetColumns) leads to an exception further down.
            # Elset columns must be recognized as such for the class to work
            # correctly.

        ampStr = StringIO()
        dloadStr = StringIO()
        checkCombinedElsets = set()

        # iteration over mining steps
        for rowCnt, row in enumerate(self.getLoadHistTableIter()):

            if odbStep is not None and row["odbStep"] < odbStep:
                # right step not found yet
                continue
            elif odbStep is not None and row["odbStep"] > odbStep:
                # reached end of this step
                break

            ampStr.write("**\n** mining step %03d %s (odb step %d, frame %d)\n"
                         % (rowCnt+1,row["miningStepName"],
                            row["odbStep"], row["odbFrame"]))

            evalNamespaceDict = dict(extraLocals)
            evalNamespaceDict.update(row)
            errorLocTmp =("%%s value in the %dth data row of the load history"
                          % (rowCnt+1))

            # iteration over elsetColumns (like "Cave", "Stopes", "Drives")
            for elsetColumn, amplist in colToGravAmp:

                if amplist is None:
                    # elset column specified with None as amplist part are
                    # interpreted as columns to be ignored here
                    # (they may serve other purposes like generating
                    # non-gravity amplitudes)
                    continue

                if isinstance(elsetColumn, basestring):
                    # ordinary elset column

                    elsetList = self._getElsetListFromColumn(
                        elsetColumn.upper(), row, "getGravAmpDload")
                    del elsetColumn  # don't use anymore may be lower case
                    if not(elsetList):
                        continue

                    # evaluation of amplitude values
                    amplist = [
                        (self.evaluateString(
                            ti, evalNamespaceDict, errorLocTmp %"time"),
                         self.evaluateString(
                             amp, evalNamespaceDict, errorLocTmp %"amplitude"))
                        for ti, amp in amplist]

                    # check amplitude times increasing monotonically
                    if not all(amplist[i+1][0] - amplist[i][0]
                               for i in range(len(amplist)-1)):
                        ampName = ", ".join("GRAVAMP_%s" % elsetName
                                            for elsetName in elsetList)
                        msg("ERROR: The tabular amplitude(s) %s is/are defined"
                            " with time not monotonically increasing. %dth data"
                            " row of the load history."
                            "\nProceeding now, but the Abaqus preprocessor"
                            " will not accept this data." % (ampName, rowCnt+1))
                        msg("WARNING: Also check sdv-initialization. There is"
                            " no warning being generated for sdv initial"
                            " values!")
                    # evaluation of time parameter
                    timepara = ("%s TIME" % self.evaluateString(
                        time, evalNamespaceDict,
                        errorLocTmp %"STEP/TOTAL TIME"))

                    # iteration over elsets
                    for elsetName in elsetList:
                        ampName = "GRAVAMP_%s" % elsetName
                        Amplitude(ampName, amplist,
                                  smooth=None, time=timepara).write(ampStr)
                        dloadStr.write(
                            "*DLOAD, OP=NEW, AMPLITUDE=%s\n"
                            "%s, GRAV, 9.81, 0., 0., -1.\n"
                            % (ampName, elsetName))

                else:
                    # two column entry

                    # column labels are in general case insensitive
                    # but not for value eval
                    elsetColumnsUpper = tuple(x.upper() for x in elsetColumn)

                    # check if there is a split elset for current two-columns
                    # in current step
                    try:
                        step2ElsetsList = self.colsStep1ToStep2Elsets[
                            (elsetColumnsUpper, rowCnt)]
                    except KeyError:
                        continue

                    # make available values from extraLocals and then both
                    # elsetColumns
                    combinedEvalDict = dict(extraLocals)

                    # make available values from the first elsetColumn
                    # This is done in a case ***sensitive*** manner.
                    combinedEvalDict[elsetColumn[0]] = Container(**row)

                    for secRowCnt, combinedElset in step2ElsetsList:

                        # check for duplicates
                        if combinedElset in checkCombinedElsets:
                            msg("WARNING: The split-off elset <%s> appears"
                                " twice when generating amplitude curves"
                                " for it. This may cause invalid output!"
                                "\nPlease check the documentation that all"
                                " parameters are well chosen."
                                "\nIn particular check that the same pairs of"
                                " coulumns are listed in the colToSDVini"
                                " argument to getStrSdvIni and in the"
                                " colToGravAmp argument to getGravAmpDload."
                                " And check that no elset is listed more than"
                                " once." % combinedElset)
                        else:
                            checkCombinedElsets.add(combinedElset)

                        # make available values from the second elsetColumn
                        secRow = self.getLoadHistRow(secRowCnt)
                        combinedEvalDict[elsetColumn[1]] = Container(**secRow)
                        ampName = "GRAVAMP_%s" % combinedElset

                        # evaluation of amplitude values
                        secErrorLocTmp = (
                            "%%s value in the %dth or %dth data row of the"
                            " load history" % (rowCnt+1, secRowCnt+1))
                        amplistEvaled = [
                            (self.evaluateString(
                                ti, combinedEvalDict,
                                secErrorLocTmp % "time"),
                             self.evaluateString(
                                 amp, combinedEvalDict,
                                 secErrorLocTmp % "amplitude"))
                            for ti, amp in amplist]

                        # check amplitude times increasing monotonically
                        if not all(amplistEvaled[i+1][0] - amplistEvaled[i][0]
                                   for i in range(len(amplistEvaled)-1)):
                            msg("ERROR: The tabular amplitude %s is defined"
                                " with time not monotonically increasing. %dth"
                                " and %dth data row of the load history."
                                "\nProceeding now, but the Abaqus preprocessor"
                                " will not accept this data."
                                % (ampName, rowCnt+1, secRowCnt+1))
                            msg("WARNING: Also check sdv-initialization. There"
                                " is no warning being generated for sdv"
                                " initial values!")
                        # evaluation of time parameter
                        timepara = ("%s TIME" % self.evaluateString(
                            time, combinedEvalDict,
                            errorLocTmp %"STEP/TOTAL TIME"))
                        Amplitude(ampName, amplistEvaled,
                                  smooth=None, time=timepara).write(ampStr)
                        dloadStr.write(
                            "*DLOAD, OP=NEW, AMPLITUDE=%s\n"
                            "%s, GRAV, 9.81, 0., 0., -1.\n"
                            % (ampName, combinedElset))

        # return results
        return (ampStr.getvalue(), dloadStr.getvalue())

    def getStrSdvIni(self, colToSDVini, extraLocals={}, extraElsetSDVList=[],
                     elsetsToSplit=None):
        r"""
        Return a string containing the *INITIAL CONDITIONS,TYPE=SOLUTION
        command.

        @param colToSDVini: A list of (column label, list of SDV values)
          tuples.

          Each column label identifies the column(s) under which the elsets are
          listed to which those SDVs are applied. E.g. 'Pit', 'Cave', 'Fill'.
          Column labels are considered case insensitive.

          "column label" may also be a tuple of two column labels identifying
          the corresponding intersection elset. This is called the
          "split elsets feature". Elsets should have been split before by a
          call to L{intersectColumns}.

          The "column label" tuple must be identical to that used in the
          elsetColumns argument of the preceeding call to L{intersectColumns},
          including the order of the items. Otherwise this item will not be
          recognized.

          Contents of the column identified by the miningStepName argument
          to self.__init__() must be unique (and not empty) for this to work
          because this name is used for the new elset names.

          String evaluation is slightly modified with this feature: for example
          instead of variable "totaltime" you'll have to write "Stope.totaltime"
          to get the corresponding value from the column "Stope".

          See also the example for the split elsets feature in the class
          description.

          The corresponding list of SDV values contains floats or string
          values. Strings are evaluated to yield the corresponding SDV value.
          See the note on L{value evaluation<LoadHistory>}.

          There must be as many items in the list as SDVs in the definition
          of the user material (and the corresponding vumat). No check is being
          performed at this stage.

          The list of SDV values might also be None. In this case the
          corresponding elset column is to be ignored for this purpose.
          A warning message about missing column labels is being suppressed in
          this case.

          Note: All column labels are converted to upper case. It is essential
          to have had all column labels in the elsetColumns argument of
          self.__init__. Otherwise this will raise a ValueError.

        @param extraLocals: A namespace dictionary that contains additional
          variables that should be accessible during value evaluation, see
          L{value evaluation<LoadHistory>}. Specify extraLocals=locals()
          in the argument list to make available all current variables.

        @param extraElsetSDVList: A list of (elset name, list of SDV values)
          tuples. Each elset name identifies an elset not contained in the
          load history. The list of SDV values must contain all initial values
          of the SDVs for those element sets.

          This list if specified results in extra line(s) for elsets not
          listed in the load history, typically this would be the elset
          'NOTMINED'.

          Example:
           >>> extraElsetSDVList=[
           >>>     ("NOTMINED",[0,0,0,0,0,0,0,1,1,9.9E6,9.9E6,9.9E6,9.9E6]),]

        @param elsetsToSplit: DEPRECATED argument, will be ignored.
        """
        # check parameter
        if not isinstance(colToSDVini, list):
            raise ValueError(
                "LoadHistory.getStrSdvIni(): expect list as colToSDVini"
                " argument, got %s." % type(colToSDVini))
        if not isinstance(extraElsetSDVList, list):
            raise ValueError(
                "LoadHistory.getStrSdvIni(): expect list as extraElsetSDVList"
                " argument, got %s." % type(extraElsetSDVList))

        # further check: all elsetColumns covered in colToSDVini?
        allElsetColumns = list()
        for elsetColumn, sdvList in colToSDVini:
            if isinstance(elsetColumn, basestring):
                allElsetColumns.append(elsetColumn)
            else:
                allElsetColumns.extend(elsetColumn)
        allElsetColumns = map(str.upper, allElsetColumns)
        missing = self.elsetColumns.difference(allElsetColumns)
        if missing:
            msg("WARNING: The following elset columns are not mentioned in the"
                " colToSDVini parameter to LoadHistory.getStrSdvIni(): %s"
                "\nIf this is intended and you want to avoid this warning then"
                " specify those columns with None instead of the list of SDV"
                " values."
                % ",".join(sorted(missing)))
            # Note: The reverse condition (colToSDVini mentions column that
            # is not in self.elsetColumns) leads to an exception further down.
            # Elset columns must be recognized as such for the class to work
            # correctly.

        output = StringIO()
        checkCombinedElsets = set()

        # write command for extra elset (from extraElsetSDVList)
        output.write("*INITIAL CONDITIONS,TYPE=SOLUTION\n")
        for elsetName, sdvList in extraElsetSDVList:
            for line in groupIter(
                    [elsetName] + ["%g" %val for val in sdvList], 8):
                output.write(", ".join(line) + "\n")

        # iteration over mining steps
        for rowCnt, row in enumerate(self.getLoadHistTableIter()):

            output.write("**\n** mining step %03d %s (odb step %d, frame %d)\n"
                         % (rowCnt+1,row["miningStepName"],
                            row["odbStep"], row["odbFrame"]))

            evalNamespaceDict = dict(extraLocals)
            evalNamespaceDict.update(row)
            errorLocStr = ("SDVini value in the %dth data row of the load"
                           " history" % (rowCnt+1))

            # iteration over elsetColumns (like "Cave", "Stopes", "Drives")
            for elsetColumn, sdvList in colToSDVini:

                if sdvList is None:
                    # elset column specified with None as sdvList part are
                    # interpreted as columns to be ignored here
                    # (they may serve other purposes like generating amplitudes)
                    continue

                if isinstance(elsetColumn, basestring):
                    # ordinary elset column

                    elsetList = self._getElsetListFromColumn(
                        elsetColumn.upper(), row, "getStrSdvIni")
                    del elsetColumn  # don't use anymore may be lower case
                    if not(elsetList):
                        continue  # ignore if no elsets

                    # evaluation of sdvinit values
                    sdvList = [
                        "%g" % self.evaluateString(
                            val,evalNamespaceDict,errorLocStr)
                        for val in sdvList]

                    # iteration over elsets
                    for elsetName in elsetList:
                        for line in groupIter([elsetName] + sdvList, 8):
                            output.write(", ".join(line) + "\n")

                else:
                    # elsetColumn is a pair of two elset column labels
                    # case insensitive in general (but not for value eval)
                    elsetColumnsUpper = tuple(x.upper() for x in elsetColumn)

                    # check if there is a split elset for current two-columns
                    # in current step
                    try:
                        step2ElsetsList = self.colsStep1ToStep2Elsets[
                            (elsetColumnsUpper, rowCnt)]
                    except KeyError:
                        continue

                    # make available values from extraLocals and then both
                    # elsetColumns
                    combinedEvalDict = dict(extraLocals)

                    # make available values from the first elsetColumn
                    # This is done in a case ***sensitive*** manner.
                    combinedEvalDict[elsetColumn[0]] = Container(**row)

                    for secRowCnt, combinedElset in step2ElsetsList:

                        # check for duplicates
                        if combinedElset in checkCombinedElsets:
                            msg("WARNING: The split-off elset <%s> appears"
                                " twice when generating sdvinit lines"
                                " for it. This may cause invalid output!"
                                "\nPlease check the documentation that all"
                                " parameters are well chosen."
                                "\nIn particular check that the same pairs of"
                                " coulumns are listed in the colToSDVini"
                                " argument to getStrSdvIni and in the"
                                " colToGravAmp argument to getGravAmpDload."
                                " And check that no elset is listed more than"
                                " once." % combinedElset)
                        else:
                            checkCombinedElsets.add(combinedElset)

                        # make available values from the second elsetColumn
                        secRow = self.getLoadHistRow(secRowCnt)
                        combinedEvalDict[elsetColumn[1]] = Container(**secRow)
                        secErrorLocStr = (
                            "SDVini value in the %dth or %dth data row of the"
                            " load history" % (rowCnt+1, secRowCnt+1))

                        sdvListEvaled = [
                            "%g" % self.evaluateString(
                                val, combinedEvalDict, secErrorLocStr)
                            for val in sdvList]

                        # write sdvini lines, iteration over elsets
                        for line in groupIter([combinedElset]+sdvListEvaled, 8):
                            output.write(", ".join(line) + "\n")

        return output.getvalue()


    def frameRowIter(self):
        """Iterate over odb frames and yield the last row corresponding to each
        frame. Ignore the equilibrium step and frame 0 of each step.
        """
        oldStep = None
        lastRow = None

        # iteration over mining steps
        for row in self.getLoadHistTableIter():
            currentStep = row["odbStep"]

            # ignore equilibrium step
            if currentStep<2:
                continue

            if oldStep!=currentStep:

                # add end of last step
                if oldStep is not None:
                    yield lastRow

                oldStep = currentStep
                oldframe = 1  # don't consider frame 0

            curframe = row["odbFrame"]
            if curframe > oldframe:
                # note: curframe can't be larger than 1 for the first row
                # in each odb step
                yield lastRow

            oldframe = curframe
            lastRow = row

        # add end of step
        if lastRow is not None:
            yield row
        return


    def getFrameTimes(self):
        r"""Retrieve time points (total time) of all time frames / end points
        of sequence periods.

        See also L{writeFrameTimes}.

        @returns: a list of frame times. A.k.a. time points at the end of
          sequence periods.

        @Note: The end of the step is always included, regardless whether it
          has actually been requested.
        """

        timepoints = [
            row["totaltime"] + row["miningStepLength"]
            for row in self.frameRowIter()]
        return timepoints


    def writeFrameTimes(self, outputFile="vusubs_param_timepoints.py"):
        """Write frame times into python script file to be imported into the
        vusubs_param.py parameter file for the Abaqus user subroutines.

        @param outputFile: can be a file name (or path) or an open file.
           Defaults to vusubs_param_frameTimes.py which is the default name
           included by vusubs_param.py if you activate
            >>> useModule_seqPeriodCounter = "varLength"
        """
        if not hasattr(outputFile, "write"):
            outputFile = open(outputFile, "w")

        outputFile.write(
            "\n"
            "# Time points for time frames for varLength-seqPeriodCounter\n"
            "#  . Will be included by vusubs_param.py if\n"
            "#    useModule_seqPeriodCounter = 'varLength' has been selected.\n"
            "#  . To create this file automatically from makeHistory add:\n"
            "#    loadhistory.writeFrameTimes()\n"
            "\n"
            "timepoints = [\n")
        timepoints = self.getFrameTimes()
        for line in groupIter(timepoints, 10, format="%g"):
            outputFile.write("    " + ", ".join(line) + ",\n")
        outputFile.write("    ]\n")
        return


    def writeCaveCouplingTiming(
            self, outputFile,
            colNameSeqPeriod="SeqPeriod", colNameATN="ATN",
            dUEvalStartMarker="dUEvalStart"):
        """Write the sequence-period to Abaqus Transfer Number relation for the
        Abaqus user subroutines. I.e. write the variable seqPeriodToATN in a
        form suitable for vusubs_param_timepoints.py to be imported into
        vusubs_param.py.

        @param outputFile: file name or file-like object, e.g.
           "vusubs_param_timepoints.py".
        @param colNameSeqPeriod: The column with that name contains
           the sequence period number. It's being checked if the numbers
           are 1,2,3,4,... (consecutive) through all the frames.
           You can specify colNameSeqPeriod=None to avoid this check and a
           possible warning that the sequence period column could not be found.
        @param colNameATN: column name for the Abaqus Transfer Numbers
        @param dUEvalStartMarker: marker-text in the column identified by
           colNameATN. This text marks the sequence period from which
           incremental dsiplacements dU3 are supposed to be recorded. If this
           marker is not found then the first occurence of a number (usually 1)
           serves as this marker.

           The variable dUEvalStartSeqPeriodIdx is set to the corresponding
           sequence period number.

           In the common case that we want to start recording dU3 for the very
           same sequence period at whose end we want the first coupling data
           exchange to happen there is no special marker. In this case the first
           entry "1" in this column will be treated as "dUEvalStartMarker" and
           dUEvalStartSeqPeriodIdx will be set to cavingStartSeqPeriodIdx-1.
        """

        cavingStartSeqPeriodIdx = None
        dUEvalStartSeqPeriodIdx = None
        seqPeriodToATN = dict()
        if colNameSeqPeriod is None:
            checkSeqPeriod = False
        else:
            checkSeqPeriod = colNameSeqPeriod in self.loadHistTable
            if not checkSeqPeriod:
                msg("WARNING: No column <%s> found for the sequence period."
                    "\nIf there is no column with a sequence period then"
                    " pass colNameSeqPeriod=None to suppress this warning."
                    "\nOtherwise rename the column (recommended) or pass the"
                    " correct column name as colNameSeqPeriod argument."
                    % colNameSeqPeriod)
        for seqPeriod, row in enumerate(self.frameRowIter(), start=1):
            if checkSeqPeriod:
                try:
                    seqPeriodCheck = int(round(row[colNameSeqPeriod]))
                except (TypeError, ValueError):
                    msg("ERROR: Illeagal value <%s> in column %s in the"
                        " loadHistory table when I expected %d."
                        % (row[colNameSeqPeriod], colNameSeqPeriod, seqPeriod))
                if (seqPeriod != seqPeriodCheck):
                    raise ValueError(
                        "ERROR: The %d-th sequence period (output frame) has a"
                        " different number %d assigned in the column <%s>."
                        " If you don't care then just remove this column."
                        % (seqPeriod, seqPeriodCheck, colNameSeqPeriod))
            atn = row[colNameATN]
            if atn is None:
                atn = -1
            else:
                try:
                    atn = int(round(float(atn)))
                except ValueError:
                    if atn != dUEvalStartMarker:
                        msg("WARNING: Found unexpected dU3EvalStart - marker:"
                            " <%s> in sequence period %d. Expected <%s> as"
                            " marker."
                            % (atn, seqPeriod, dUEvalStartMarker))
                    if dUEvalStartSeqPeriodIdx is None:
                        dUEvalStartSeqPeriodIdx = seqPeriod
                    else:
                        msg("WARNING: Found multiple dU3EvalStart - marker, or"
                            " unrecognizable Abaqus Transfer Number <%s> in"
                            " sequence period %d." % (atn, seqPeriod))
                    atn = -1
            if cavingStartSeqPeriodIdx is None and atn>0:
                cavingStartSeqPeriodIdx = seqPeriod+1
                if atn!=1:
                    msg("WARNING: First Abaqus Transfer Number in column %s"
                        "is %d instead of 1, sequence period %d."
                        % (colNameATN, atn, seqPeriod))
                if dUEvalStartSeqPeriodIdx is None:
                    dUEvalStartSeqPeriodIdx = seqPeriod

            seqPeriodToATN[seqPeriod] = atn

        seqPeriodToATN = [seqPeriodToATN.get(cnt+1, -1)
                          for cnt in range(max(seqPeriodToATN))]
        if not hasattr(outputFile, "write"):
            outputFile = open(outputFile, "w")
        outputFile.write(
            "\n"
            "# sequence period to Abaqus Transfer Number translation table:\n"
            "# Python: ATN = seqPeriodToATN[seqPeriod-1]\n"
            "# Fortran: ATN = seqPeriodToATN(seqPeriod)\n"
            "\n"
            "seqPeriodToATNLen = %d\n"
            "seqPeriodToATN = [\n" % len(seqPeriodToATN))
        for line in groupIter(seqPeriodToATN, 10, format="%d"):
            outputFile.write("    " + ", ".join(line) + ",\n")
        outputFile.write("    ]\n")

        outputFile.write("""
# The cave criterion
# will be evaluated for the first time at the beginning of the sequence
# period with seqPeriodIdx=cavingStartSeqPeriodIdx.
# cavingStartSeqPeriodIdx must be >= seqPeriodIdxStart+1 because dU3 can
# only be calculated at the end of seqPeriodIdx==seqPeriodIdxStart and this
# will then be effective for seqPeriodIdx==seqPeriodIdxStart+1 earliest.
# Example with eqStepTime=40, seqPeriodLength=2.5, monthly step starting
# Jan 2006:
# First rings are being fired Apr 2006 which is steptime 10 = totaltime 50,
# odb frame 4. This is the first coupling frame, i.e. frame 4 is coupling
# step 1 aka Abaqus Transfer Number ATN=1. At the beginning of the
# following sequence period with seqPeriodIdx=5 rock can be switched to
# cave for the first time: set cavingStartSeqPeriodIdx=5.

# DU3 will be calculated for the first time at the end of this sequence
# period (i.e. the du in this sequence period is recorded).
dUEvalStartSeqPeriodIdx = %d

# Sequence period during which rock can turn into cave for the first time.
# Usually one sequence period after dUEvalStartSeqPeriodIdx
# - cavingStartSeqPeriodIdx>seqPeriodIdxStart must always apply
# - cavingStartSeqPeriodIdx>dUEvalStartSeqPeriodIdx must always apply
# - cavingStartSeqPeriodIdx must correspond to frame1stStop+1 in
#   jobInfoCoupling.py:
#   . For non-restart runs: cavingStartSeqPeriodIdx=frame1stStop+1
#   . For restart runs: cavingStartSeqPeriodIdx=frame1stStop+seqPeriodIdxStart.
# cavingStartSeqPeriodIdx = dUEvalStartSeqPeriodIdx + 1
cavingStartSeqPeriodIdx = %d
""" % (dUEvalStartSeqPeriodIdx, cavingStartSeqPeriodIdx))
        return


    def writeCaveCouplingRelationTable(
            self, outputFile, colNameATN="ATN",
            odbStepConversion=lambda x: x-2 if x>2 else x):
        """Write the Abaqus frame to FS4 mining step relation into the variable
        fileNumFrameList to be imported into jobInfoCoupling.py.

        @param outputFile: file name or file-like object, e.g.
           "jobInfoCoupling.py".
        @param colNameATN: column name for the Abaqus Transfer Numbers

        @Note: The variable fileNumFrameList distinguishes between FS4 file
        number and Abaqus Transfer Number. Typically these numbers are the same
        in our coupling framework. However in general it can be used (and has
        been used) to let the Abaqus and FS4/Cavesim schedules deviate.

        This function always assigns identical FS4 file numbers and ATNs.
        """

        # initialize output
        if not hasattr(outputFile, "write"):
            outputFile = open(outputFile, "w")
        outputFile.write("""
# list of tuples for each coupling sequence period:
# (FS4 file number (handshake number), Abaqus transfer number (ATN),
#  odb step number, frame number)
fileNumFrameList = [
    # [FS4 nb, ATN, odb step, odb frame]
""")

        for cnt, row in enumerate(self.frameRowIter()):
            atn = row[colNameATN]
            try:
                atn = int(atn)
            except (ValueError, TypeError):
                continue

            odbStep = row["odbStep"]
            if callable(odbStepConversion):
                odbStep = odbStepConversion(odbStep)

            outputFile.write(
                "    (%3d, %3d, %d, %3d),\n"
                % (atn, atn, odbStep, row["odbFrame"]))

        # finalize output
        outputFile.write("    ]\n")


    def getStrTimePointsDef(self, timePointsName=None, odbStep=None):
        r"""
        Return a string containing the *TIME POINTS commands for the field
        output requests in an Abaqus input file.

        @param timePointsName: name argument for the *TIME POINTS command.
          Defaults to "TPOINT_X" where X is the odb step number. Only specify
          this if there is just one step (or odbStep is specified).

        @param odbStep: return time points for this odb step. If not specified
          take the first step in the loadhistory.

        @Note: The current implementation always includes the end of the step
          in the timepoints, regardless whether it has actually been requested.
          This does not matter because Abaqus writes results at the end of the
          step no matter if a corresponding time point has been specified.

        @Raises ValueError: If there are more than one steps and
          timePointsName has been specified.
        """

        tpDict = TimePointsDict()
        oldStep = None
        timeAtEnd = None

        # iteration over mining steps
        for rowCnt, row in enumerate(self.getLoadHistTableIter()):
            currentStep = row["odbStep"]

            if (odbStep is not None) and (currentStep < odbStep):
                # right step not found yet
                continue
            elif (odbStep is not None) and (currentStep > odbStep):
                # reached end of this step
                break

            if oldStep != currentStep:

                # add end of last step
                if oldStep is not None:
                    timepoints.append(timeAtEnd)

                # Abaqus name for the time points list
                if timePointsName is None:
                    curTPName = "TPOINTS_%d"%currentStep
                elif (oldStep is not None):
                    # this is not the first step, then one timePointsName is
                    # not suitable
                    raise ValueError(
                        "LoadHistory.getStrTimePointsDef() only accepts"
                        " timePointsName if there is no more than one odbStep.")
                else:
                    curTPName = timePointsName

                # initialize new time points list
                timepoints = list()
                tpDict[curTPName] = timepoints
                oldframe = 1  # don't write time point 0.0
                oldStep = currentStep

            curframe = row["odbFrame"]
            if curframe > oldframe:
                timepoints.append(timeAtEnd)

            oldframe = curframe
            timeAtEnd = row["steptime"] + row["miningStepLength"]

        # add end of step
        if timeAtEnd is not None:
            timepoints.append(timeAtEnd)

        output = StringIO()
        tpDict.writeAll(output)

        return output.getvalue()


    def generateAmplitude(
            self, name, rowToTimeValTuples,
            outputAsString=True,
            **kwargs):
        r"""Generates an Abaqus amplitude from the loadHistory data by means of
        a user supplied function.

        This method can generate an arbitrary number of (time, value)-pairs
        for each mining step. It's a more flexible successor to the older
        method L{colToAmplitude}.

        The Amplitude object is basically a list of (time, amplitude value)
        tuples with extra attributes corresponding to the options of the ABAQUS
        *AMPLITUDE command. It's got a write method with the open output file
        object as argument.

        Example of a sawtooth shaped velocity amplitude. Suppose there is a
        column "BC_DISP" in the sequence-csv file that specifies the relative
        displacement in each step:
         >>> def ampVals(row, amp):
         >>>     if row["BC_DISP"] is None:
         >>>         return []
         >>>     steptime = row["steptime"]
         >>>     stepLength = 0.6*row["miningStepLength"]
         >>>     maxVel = row["BC_DISP"] *2.0/stepLength
         >>>     return [ [steptime               , 0.0],
         >>>              [steptime+0.5*stepLength, maxVel],
         >>>              [steptime+stepLength    , 0.0] ]
         >>>
         >>> loadcase = LoadHistory("sequence_Q01.csv", ...)
         >>> displacementAmp = loadcase.generateAmplitude(
         >>>     "BC_VEL", ampVals, outputAsString=False)
         >>> ...
         >>> output = open(abqMainName, "w")
         >>> ...
         >>> displacementAmp.write(output)

        Second example, using the second argument "amp" to the generator
        function. Here the column "ampModif" contains the offsets by which
        the amplitude is to be ramped up during the first half of the mining
        step. Shown here is just the generator function to be supplied to the
        rowToTimeValTuples argument:
         >>> def ampVals(row, amp):
         >>>     if row["ampModif"] is None:
         >>>         return []
         >>>
         >>>     if len(amp)==0:
         >>>         prevVal = 0.0
         >>>     else:
         >>>         prevVal = amp[-1][1]
         >>>
         >>>     return [ [row["steptime"], prevVal],
         >>>              [row["steptime"]+0.5*row["miningStepLength"],
         >>>               prevVal+row["ampModif"]] ]

        @param name: Name for the Amplitude object to be returned.

        @param rowToTimeValTuples: A generator function that is called for
          each row in the loadhistory.

          It's supplied with two arguments:
          Firstly a dictionary that contains
          the column names as keys and the corresponding values from the
          current row. This row-dictionary includes the extra columns
          automatically generated by self.__init__(). See the description of
          the instance variable L{loadHistTable} for details.

          The second argument is the amplitude itself as far as it's been
          generated so far. This might be used to make use of previously
          generated values. Suppose the second argument is called "amp"
          then amp[-1] will be the preceding [time, value]-tuple (actually
          a list with two items). This could also be used to retrospectively
          modify earlier values.

          The function is expected to return a sequence of [time, value]-
          lists to be appended to the amplitude as it's being generated. The
          sequence might of course be empty as well. Note that the time-values
          must be steptime or totaltime corresponding to the optional time
          argument to the amplitude, defaulting to steptime.

          Note: The global variable name space for this function is the globals
          namespace where the function rowToTimeValTuples has been defined.
          Those global variables of course have the current values at the time
          the function is actually called, not at the time the function has
          been defined. So it's straightforward to define some extra parameter
          for the function rowToTimeValTuples as global variables.

        @param outputAsString: if True return the Abaqus input file string
          corresponding to the amplitude. This includes the trailing newline.

        @param kwargs: arguments passed to the
          L{Amplitude<bae.abq_model_02.container.Amplitude>} object that is to
          be returned.

        @kwarg time: Either "STEP TIME" or "TOTAL TIME". See the ABAQUS manual.
          If not specified Abaqus assumes STEP TIME.

        @returns: an L{Amplitude<bae.abq_model_02.container.Amplitude>}-object
          with the values from the given column in the loadhistory table.

        @Note: Specify all but the first two arguments as keyword argument to
          allow for easier upgrade of the function interface!
        """

        amp = Amplitude(name=name, **kwargs)

        for row in self.getLoadHistTableIter():

            amp.extend(rowToTimeValTuples(row, amp))

        # return amplitude
        if outputAsString:
            output = StringIO()
            amp.write(output)
            return output.getvalue()
        else:
            return amp

    def colToAmplitude(self, value, name=None, filterVal=None,
                       timeOffsetFraction=0.0, extraLocals={},
                       outputAsString=True,
                       **kwargs):
        r"""Returns an Abaqus amplitude corresponding to value evaluated for
        each row in the loadhistory table.

        This is the simpler (and older) method to generate an amplitude card
        from the loadhistory table. It creates (at most) one
        (time, value)-point for each mining step (row in the loadhistory
        table).

        You may want to consider the method L{generateAmplitude}() if you need
        more flexibility.

        The Amplitude object is basically a list of (time, amplitude value)
        tuples with extra attributes corresponding to the options of the ABAQUS
        *AMPLITUDE command. It's got a write method with the open output file
        object as argument.

        A simple example:
         >>> loadcase = LoadHistory("sequence_Q01.csv", ...)
         >>> displacementAmp = loadcase.colToAmplitude("BC_DISP")
         >>> ...
         >>> output = open(abqMainName, "w")
         >>> ...
         >>> displacementAmp.write(output)

        A more elaborate example:
         >>> displacementAmp = loadcase.colToAmplitude("BC_DISP")
         >>> # prepare other values from a "factor"-column and the
         >>> # automacally generated columns steptime and miningStepLength
         >>> secondValDict = dict(loadcase.colToAmplitude(
         >>>         "factor * (steptime+miningStepLength)"))
         >>> modifiedAmp = Amplitude(name="modamp", listInit=[
         >>>         (time, val * secondValDict[time])
         >>>         for time, val in displacementAmp])

        Another example (maybe more sensible): Computes a piecewise linear
        velocity amplitude that reaches the displacements specified in the
        BC_DISP column of the loadhistory.

        This example creates an amplitude with more than one point per mining
        step. This is done by separately processing the initial amplitude
        extracted by colToAmplitude(). Gero thinks that it makes more sense
        to use the method L{generateAmplitude}() in this use case.
         >>> displacementAmp = loadcase.colToAmplitude(
         >>>     value="AMP_BC_DISP",
         >>>     filterVal="odbStep==2",
         >>>     timeOffsetFraction=1.0)
         >>> peakTimeFraction = 0.5 # must be greater zero and less than one!
         >>> lastStepTime = 0.0
         >>> lastVal = 0.0
         >>> velocityAmp = Amplitude(name="AMP_BC_VELOCITY")
         >>> velocityAmp.append( [0.0, 0.0] )
         >>> for i, (time, val) in enumerate(displacementAmp):
         >>>     dt = time - lastStepTime
         >>>     dx = val - lastVal
         >>>     velocityAmp.append( [lastStepTime+peakTimeFraction*dt,
         >>>                          2.0*dx/dt] )
         >>>     velocityAmp.append( [time, 0.0] )
         >>>     lastStepTime = time
         >>>     lastVal = val
         >>> ...
         >>> velocityAmp.write(output)

        @param value: column name / expression serving as value for the
          amplitude. Might be a column name in the loadhistory csv file or
          might be an expression containing such column names and other
          variables to be taken from extraLocals.

          If the expression evaluates to None the corresponding row is being
          ignored and simply does not lead to any new step in the resulting
          amplitude.

          Note that you need to make sure that the expression is valid for each
          single row. For instance a value "ampVal" is acceptable even if the
          column ampVal only contains values in some of its rows. (Only those
          rows will generate time steps in the amplitude.) But a value
          "5*ampVal" will generate an error saying that you can't multiply an
          integer (5) to None which is the result of an empty ampVal-row. You
          would need to make use of the filterVal argument in this case:
          filterVal="ampVal is not None". The simpler filterVal="ampVal" might
          do as well but would also filter out every zero value in ampVal.

        @param filterVal: column name / expression serving as filter value.
          There will only be a value in the resulting amplitude object for
          mining steps (i.e. rows in the loadhistory table) at which the
          filter-expression evaluates to True (after conversion to bool).

        @param name: Name for the Amplitude object to be returned. If not
          specified defaults to the value argument (This might not be valid
          Abaqus identifier!).

        @param timeOffsetFraction: controls the time values in the amplitude
          definition. If 0.0 then the times in the amplitude will be the start
          time of each mining step. If 1.0 then it will be the end times of each
          mining step. 0.5 is in the middle of the mining step. Defaults to 0.0
          i.e. start of mining step.

        @param extraLocals: A namespace dictionary that contains additional
          variables that should be accessible during value evaluation, see
          L{value evaluation<LoadHistory>}. Specify extraLocals=locals()
          in the argument list to make available all current variables.

        @param outputAsString: if True return the Abaqus input file string
          corresponding to the amplitude

        @param kwargs: arguments passed to the
          L{Amplitude<bae.abq_model_02.container.Amplitude>} object that is to
          be returned.

        @kwarg time: Either "STEP TIME" or "TOTAL TIME". See the ABAQUS manual.
          If not specified Abaqus assumes STEP TIME.

        @returns: an L{Amplitude<bae.abq_model_02.container.Amplitude>}-object
          with the values from the given column in the loadhistory table.

        @Note: Specify all but the value argument as keyword argument to
          allow for easier upgrade of the function interface!
        """

        if name is None:
            name = value

        timeCol = kwargs.get("time", "steptime").replace(" ", "").lower()

        amp = Amplitude(name=name, **kwargs)

        for row in self.getLoadHistTableIter():

            # prepare evaluating values
            evalNamespaceDict = extraLocals.copy()
            evalNamespaceDict.update(row)

            # check filterVal
            if filterVal is not None:
                if not(self.evaluateString(
                        filterVal, evalNamespaceDict,
                        "the filter value for amplitude %s at totaltime %g"
                        % (name, row["totaltime"]))):
                    continue

            # evaluate value
            val = self.evaluateString(
                value, evalNamespaceDict,
                "the value for amplitude %s at totaltime %g"
                % (name, row["totaltime"]))

            # assign to amplitude
            if val is not None:
                time = row[timeCol]+timeOffsetFraction*row["miningStepLength"]
                amp.append([time, val])

        # return amplitude
        if outputAsString:
            output = StringIO()
            amp.write(output)
            return output.getvalue()
        else:
            return amp

    def getElsetsAmp(self, colToAmp, ampNameTemplate="AMP_%s",
                     extraLocals={}, outputAsString=True, **kwargs):
        r"""Generate amplitudes from column(s) of set labels in the loadhistory
        csv file. This can be used for example to generate *temperature
        amplitudes for nsets listed in the loadhistory.

        Example:
         >>> ampNameTemplate="TEMPAMP_%s"
         >>> nsetNames, ampString = loadhistory.getElsetsAmp(
         >>>     { "GS_Installer" : [
         >>>         ("steptime+0.5*miningStepLength", 0), # start time, value
         >>>         ("steptime+0.6*miningStepLength", 1), # end time, value
         >>>         ], },
         >>>     ampNameTemplate )
         >>> output.write(ampString)
         >>> for nset in nsetNames:
         >>>     ampName = ampNameTemplate % nset
         >>>     output.write("*TEMPERATURE, AMPLITUDE=%s\n" % ampName)
         >>>     output.write("%s, 1.\n" % nset)

        Another one:
         >>> ampNameTemplate="TEMPAMP_%s"
         >>> nsetNames, tempAmpString = loadhistory.getElsetsAmp(
         >>>     { "BOLTS" : [
         >>>         ("totaltime+0.8*miningStepLength", 0), # start time, value
         >>>         ("totaltime+0.9*miningStepLength", 1), # end time, value
         >>>         ], },
         >>>     ampNameTemplate, extraLocals=locals(), time="TOTAL TIME")
         >>> tempCardsStr = "".join(
         >>>     "*TEMPERATURE, AMPLITUDE=%s\n%s, 1.\n"
         >>>     % (ampNameTemplate%nset, nset)
         >>>     for nset in nsetNames)
         >>> del ampNameTemplate, nsetNames

        @param colToAmp: a dict {column label: list of (time, value)-tuples};
          times and values are subject to L{value evaluation<LoadHistory>}.

          Note: All column labels are converted to upper case. It is essential
          to have had all column labels in the elsetColumns argument of
          self.__init__. Otherwise this will raise a ValueError.

        @param ampNameTemplate: template for the amplitude name corresponding
          to the label found in the loadhistory csv file.

        @param outputAsString: if True return the Abaqus input file string
          corresponding to the amplitudes, otherwise return list of
          L{Amplitude<bae.abq_model_02.container.Amplitude>}-objects

        @param extraLocals: A namespace dictionary that contains additional
          variables that should be accessible during value evaluation, see
          L{value evaluation<LoadHistory>}. Specify extraLocals=locals()
          in the argument list to make available all current variables.

        @param kwargs: arguments passed to the Amplitude object that is to
          be returned.

        @kwarg time: Either "STEP TIME" or "TOTAL TIME". See the ABAQUS manual.
          If not specified Abaqus assumes STEP TIME.

        @returns: a (labelNames, amplitudes) tuple. labelNames is a list of
          (possibly elset- or nset-) names as listed in the corresponding
          columns in the loadhistory csv file. amplitudes is a list of
          L{Amplitude<bae.abq_model_02.container.Amplitude>}-objects or an
          Abaqus input file string corresponding to the amplitudes, see the
          outputAsString argument. labelNames and amplitudes are ordered
          chronologically.
        """

        ampDict = OrderedDict()
        labelNames = OrderedDict()  # actually needed as ordered set only

        # convert all column labels to upper case
        colToAmp = dict(
            (k.upper(), v) for k,v in colToAmp.iteritems())

        # columns to be treated as multi value columns:
        if set(colToAmp).difference(self.elsetColumns):
            raise ValueError(
                "Elset column(s) %s not recognized as elsetColums (as given"
                " by the corresponding argument to LoadHistory.__init__)."
                " Please update the elsetColumns argument."
                % ",".join(sorted(set(colToAmp).difference(self.elsetColumns))))

        # iterate over steps in the load history
        for rowCnt, row in enumerate(
            # self.getLoadHistTableIter() not feasible because of other columns
            self.loadHistTable.rowIterColToList(
                listColumns=self.elsetColumns, removeNones=True)):

            evalNamespaceDict = dict(extraLocals)
            evalNamespaceDict.update(row)
            errTmpTime =(
                "%%d-th time value in the %dth data row of the load history"
                % (rowCnt+1))
            errTmpVal = (
                "%%d-th value in the %dth data row of the load history"
                % (rowCnt+1))

            # collect labels, iteration over labelColumns (like "GS")
            for labelColumn, ampTemplate in colToAmp.iteritems():
                try:
                    labelList = row[labelColumn]
                except KeyError:
                    continue
                if len(labelList)==0:
                    continue   # only proceed if labels exist

                # new (time, value)-tuples
                newValues = [
                    (self.evaluateString(
                        time, evalNamespaceDict, errTmpTime % ampCnt),
                     self.evaluateString(
                         val, evalNamespaceDict, errTmpVal % ampCnt))
                    for ampCnt, (time, val) in enumerate(ampTemplate)]

                # add those to the correspoding amplitudes
                for label in labelList:
                    # store label
                    labelNames[label] = None

                    # store new values to amplitude
                    ampName = ampNameTemplate % label
                    try:
                        amp = ampDict[ampName]
                    except KeyError:
                        amp = list()
                        ampDict[ampName] = amp
                    amp.extend(newValues)

        # prepare final labelNames
        labelNames = labelNames.keys()

        # sort amplitude tuples and prepare final amplitude list
        amplitudes = [
            Amplitude(name=name, listInit=sorted(amp), **kwargs)
            for name, amp in ampDict.iteritems()]

        if outputAsString:
            output = StringIO()
            for amp in amplitudes:
                amp.write(output)
            amplitudes = output.getvalue()

        return (labelNames, amplitudes)

    def getStepFrameToATN(self, colNameATN="ATN"):
        """Creates the odb-frame to cavesim-step association required for
        combined vtk files.
        See L{utils.odbToPointsData_02.Converter_CombinedVtk_FS4} and friends.

        Creates the equivalent to this old formulation in the config file:
         >>> stepFrameToCsStep = {
         >>>     2: [-1,]*(frame1stStop+1) + range(1,nbCouplingSteps+1)
         >>>        # if there are further steps without coupling:
         >>>        + [nbCouplingSteps]*999 }

        @param colNameATN: column label of the Abaqus Transfer Number (ATN)
           values for the coresponding odb-frame.
        """

        stepFrameATN = defaultdict(dict)
        for row in self.getLoadHistTableIter():
            stepFrameATN[row["odbStep"]][row["odbFrame"]] = row[colNameATN]
        stepFrameATN.default_factory = None  # convert into ordinary dict

        stepFrameToCsStep = dict()
        lastCsStep = -1
        for odbStep in sorted(stepFrameATN):
            lastFrame = self.odbStepLastFrameDict[odbStep]
            # initialize with defaults
            csStepList = [lastCsStep,]*(lastFrame + 1)
            for odbFrame in range(1, lastFrame):
                csStep = stepFrameATN[odbStep][odbFrame]
                try:
                    csStep = int(csStep)
                    if csStep<1:
                        raise ValueError("Just ignore CS steps<1")
                except (ValueError, TypeError):
                    csStepList[odbFrame+1] = lastCsStep
                else:
                    csStepList[odbFrame+1] = csStep
                    lastCsStep = csStep
            stepFrameToCsStep[odbStep] = csStepList

        return stepFrameToCsStep

    #} end of Abaqus input file commands

    #------------------------------------------------------------------
    #{ LS-Dyna keyword file commands

    def dynaAddSequence(
            self, model, elsetsDict,
            colsToLoadStiffen=[], colsToLoadRemove=[],
            colsToInitStressSolid=[],
            newNameTemplate="{oldName}_{elsetColumns}_{miningStepName}",
            extraLocals={}):
        r"""
        Add LOAD_STIFFEN_PART commands and so on to the given DynaModel.

        Examples:
        =========

        Variable timing for Pit, being defined in the column T_EndExc. And
        timing relative to the mining step length for Stopes. Note: the columns
        "totaltime" and "miningStepLength" do not have to be defined in the
        load history csv file, they are generated automatically,
        see L{__init__}:
         >>> self.dynaAddSequence([
         ...     ("Pit",    [(0.0, 1.0),
         ...                 ("totaltime", 1.0),
         ...                 ("totaltime+T_EndExc", 1e-5),
         ...                 (999.9, 1e-5),
         ...                ]),
         ...     ("Stopes", [(0.0, 1.0),
         ...                 ("totaltime", 1.0),
         ...                 ("totaltime+0.25*miningStepLength", 1e-5),
         ...                 ("totaltime+1.5", 1e-5),
         ...                 ("totaltime+1.6", 0.8),
         ...                 (999.9, 0.8)]),
         ...     ])

        Variable amplitude value according to the column EndGrav in the load
        history csv file and a variable drivesFact defined locally:
         >>> drivesFact = 0.1
         >>> self.dynaAddSequence(
         ...     [ ...
         ...       ("Pit", [...,
         ...                ("totaltime+T_EndExc", "drivesFact*endGrav"),
         ...                ...]),
         ...     ],
         ...     extraLocals=locals()
         ...     )

        @param model: a L{bae.dyna_model_01.DynaModel} object.

        @param elsetsDict: A dictionary like L{bae.abq_model_02.Model}.elset.

        @param colsToLoadStiffen: A list of ([list of column labels],
          [list of amplitude tuples]) tuples.

          Each column label identifies the column(s) under which the elsets are
          listed to which this timing is supposed to apply to. E.g. 'Pit',
          'Cave', 'Fill'. Instead of a list of column labels a single label is
          accepted as well.

          The list of amplitude tuples is a list of (x,y) or (time, value)
          tuples defining the load stiffen curve.

          Each of the x- and y-values may be floats, in that case they apply to
          all loadhistory steps. Or they may be a strings evaluating to the
          corresponding value. The time value will very rarely be a float but
          rather a string like "totaltime+xyz".
          See the note on L{value evaluation<LoadHistory>}.

          Note: All column labels are converted to upper case.

        @param colsToLoadRemove: A list of ([list of column labels], start
          time, end time ) tuples.

          Each column label identifies the column(s) under which the elsets are
          listed to which this timing is supposed to apply to. E.g. 'Pit',
          'Cave', 'Fill'. Instead of a list of column labels a single label is
          accepted as well.

          The start and end time define the start and end time of the removal.
          They will usually be strings evaluating to the corresponding time
          value. See the note on L{value evaluation<LoadHistory>}.

          Note: All column labels are converted to upper case.

        @param colsToInitStressSolid: A list of ([list of column labels], nint,
          nhisv, nthint, nthhsv, initStress, initEps, initHisv) tuples.

          nint is the number of integration points, nhisv the number of history
          varibles per integration point, nthint the number of thermal
          integration points per element, nthhsv the number of thermal history
          variables per integration point. The last two must be zero currently.

          initStress is a list of nint 6-tuples of stress values or an
          object of a subclass of L{StressDistribution}, e.g.
          L{StressGeostatic}.

          initEps is a list of nint floats, the effective plastic strains
          initHisv is a list of nint tuple of nhisv floats to initialize the
          history variables

          Note: All column labels are converted to upper case.

        @param newNameTemplate: Template for the name / header of new parts.
          May contain the template strings {oldName}, {elsetColumns} and
          {miningStepName} for substitution with str.format().
          defaults to "{oldName}_{elsetColumns}_{miningStepName}"

        @param extraLocals: A namespace dictionary that contains additional
          variables that should be accessible during value evaluation, see
          L{value evaluation<LoadHistory>}. Specify extraLocals=locals()
          in the argument list to make available all current variables.
        """
        # check parameter, convert single string elsetColumns to list
        for colsToXXX in (colsToLoadStiffen, colsToLoadRemove,
                          colsToInitStressSolid):
            if not isinstance(colsToXXX, list):
                raise ValueError(
                    "LoadHistory.dynaAddSequence(): expect list as"
                    " each of the colsTo... arguments, got %s."
                    % type(colsToXXX))
            for i, line in enumerate(colsToXXX):

                # make sure each item is a list
                if not isinstance(line, list):
                    line = list(line)
                    colsToXXX[i] = line

                # if first item -elsetColumns- is string then convert to list
                # of strings
                elsetColumns = line[0]
                if isinstance(elsetColumns, basestring):
                    elsetColumns = [elsetColumns,]
                # convert to upper case and store back
                line[0] = map(str.upper, elsetColumns)

        # iteration over mining steps
        for rowCnt, row in enumerate(self.getLoadHistTableIter()):

            evalNamespaceDict = dict(extraLocals)
            evalNamespaceDict.update(row)
            errorLocTmp =("%%s value in the %dth data row of the load history"
                          % (rowCnt+1))

            # iteration over elsetColumns (like "Cave", "Stopes", "Drives")
            # ... for LoadStiffen
            for elsetColumns, amplist in colsToLoadStiffen:
                elsetList = [elsetName
                             for col in elsetColumns
                             for elsetName in row.get(col,[])]
                if not(elsetList):
                    continue

                # evaluation of amplitude values
                amplist = [
                    (self.evaluateString(
                        ti, evalNamespaceDict, errorLocTmp % "time"),
                     self.evaluateString(
                         amp, evalNamespaceDict, errorLocTmp % "amplitude"))
                    for ti, amp in amplist]

                # split off parts and create set_part
                elems = set(e for elsetName in elsetList
                            for e in elsetsDict[elsetName])
                newName = newNameTemplate.format(
                    oldName="{oldName}",   # will be assigned by splitParts
                    elsetColumns="_".join(col.upper() for col in elsetColumns),
                    miningStepName=row["miningStepName"])
                parts = model.splitParts(elems, name=newName)
                setPartId = model.setPart.append(parts)

                # define curve and check time points increasing monotonically
                tableId = model.table.append("curve", amplist)
                if not all(amplist[i+1][0] - amplist[i][0]
                           for i in range(len(amplist)-1)):
                    msg("ERROR: The tabular amplitude curve %d is defined"
                        " with time not monotonically increasing. %dth data"
                        " row of the load history."
                        "\nProceeding now, but the data may not make any"
                        " sense." % (tableId, rowCnt+1))

                # add load stiffen cmd to model
                model.loadStiffenPartSet[setPartId] = tableId

            # iteration over elsetColumns (like "Pit", Drives")
            # ... for LoadRemove
            for elsetColumns, t0, t1 in colsToLoadRemove:
                elsetList = [elsetName
                             for col in elsetColumns
                             for elsetName in row.get(col,[])]
                if not(elsetList):
                    continue

                # evaluation of start and end time values
                t0 = self.evaluateString(
                    t0, evalNamespaceDict, errorLocTmp % "removal start time")
                t1 = self.evaluateString(
                    t1, evalNamespaceDict, errorLocTmp % "removal end time")

                # split off parts and create set_part
                elems = set(e for elsetName in elsetList
                            for e in elsetsDict[elsetName])
                newName = newNameTemplate.format(
                    oldName="{oldName}",   # will be assigned by splitParts
                    elsetColumns="_".join(col.upper() for col in elsetColumns),
                    miningStepName=row["miningStepName"])
                parts = model.splitParts(elems, name=newName)
                setPartId = model.setPart.append(parts)

                # add load remove cmd to model
                model.loadRemovePartSet[setPartId] = (t0, t1)

            # iteration over elsetColumns (like "Pit", Drives")
            # ... for colsToInitStressSolid
            for (elsetColumns,
                 nint, nhisv, nthint, nthhsv,
                 initStress, initEps, initHisv) in colsToInitStressSolid:
                elsetList = [elsetName
                             for col in elsetColumns
                             for elsetName in row.get(col,[]) ]
                if not(elsetList):
                    continue

                # evaluation of initStress (initial stress)
                initStress = self.evaluateString(
                    initStress, evalNamespaceDict, errorLocTmp%"initStress")

                # evaluation of initHisv (initial history variables)
                initHisv = self.evaluateString(
                    initHisv, evalNamespaceDict, errorLocTmp%"initHisv")

                # get element set
                elems = set(e for elsetName in elsetList
                            for e in elsetsDict[elsetName])

                model.updateInitStressSolidData(
                    elems,
                    numbers=(nint, nhisv, nthint, nthhsv),
                    initStress=initStress,  # stress components
                    initEps=initEps,   # EPS, effective plastic strain
                    initHisv=initHisv,  # initial history variables
                    )
        return

    #} end of LS-Dyna keyword file commands

    #------------------------------------------------------------------
    #{ post data file commands

    def writePostDataFile(
            self, outputFile,
            postTypeCols=[('Excav', None),
                          ('Support',['Drives']),
                          ('Cave', ['Cave'])],
            restartFrom=None,
            useStepName=False,
            colNameATN="ATN"):
        r"""
        Write the python script file that will contain all post set
        definitions, frame names (and so on).

        @param outputFile: An already open file or a file name for the
          resulting file.

        @param postTypeCols: A list of (post set type, list of column labels)
          tupels specifying elsets (identified by the column label) for post
          set types.

          In place of the list of column labels you may also specify None. In
          this case all elsets that occur in any of the columns will be
          addressed.

          E.g: [('Support', ['DRIVES','FILL']), ('Cave', ['Cave']),
          ('Excavation',None),]

          Note: The post set type strings should be capitalized words,
          like "Support". Neither "support" nor "SUPPORT" form nice CamelCase
          identifiers in the resulting python script file.

          Note: All column labels are converted to upper case. It is essential
          to have had all column labels in the elsetColumns argument of
          self.__init__. Otherwise this will raise a ValueError.

        @param restartFrom: optional. If specified then create a post data file
          for a restart run that starts with the specified step. Typically set
          to 3 for the restart analysis in case of a two-step base run. Note:
          Omit this argument for the base run post data file. It will then also
          contain data for step 3 which will usually be ignored. (But can
          become useful if you join the odbs or for other special usecases.)

        @param useStepName: If True then write a string like "Step-2" as step
          identifier rather than the number we used in ancient times.
          If useStepName is set to True then restartFrom must not be given.

        @param colNameATN: optional. Column label of the Abaqus Transfer Number
           (ATN) values for the coresponding odb-frame.

        @Note: The elsets listed in the resulting post data file lack the
          'PART-1-1.' prefix that is needed in the abaqusMacros-case. It has to
          be prepended by that script! (Other scripts might not have to take it
          away though. Removing is usually harder than adding.)
        """
        # check parameter
        if not isinstance(postTypeCols, list):
            raise ValueError(
                "LoadHistory.writePostDataFile(): expect list as postTypeCols"
                " argument, got %s." % type(postTypeCols))

        if isinstance(outputFile, basestring):
            output = open(outputFile, 'w')
        else:
            output = outputFile

        if useStepName and restartFrom:
            raise ValueError(
                "LoadHistory.writePostDataFile(): Arguments useStepName and"
                " restartFrom may not be given at the same time!")

        if useStepName:
            def stepString(i):
                return '"Step-%d"' % i
        else:
            def stepString(i):
                return '%d' % i

        separator = ("# " + "-"*75 + "\n")


        #--- write sequence elsets
        output.write(separator)
        output.write("# postsets\n")

        for setType, colNameList in postTypeCols:
            msg("Processing post set type %s with columns %s"
                %(setType, colNameList), debugLevel=10)

            # if the list of column labels is None: take all elset columns
            if colNameList is None:
                colNameList = self.elsetColumns

            output.write("seqElsets%s = {\n" % setType)
            lastOdbStep = 0
            lastOdbFrame = 1
            elsetList = list()

            for rowCnt, row in enumerate(self.getLoadHistTableIter()):
                odbStep = row["odbStep"]
                odbFrame = row["odbFrame"]

                if restartFrom and odbStep<restartFrom:
                    # suppress output for steps before restart-start-step
                    pass

                elif (restartFrom and odbStep==restartFrom
                      and odbStep>lastOdbStep):
                    # write first step data
                    output.write(
                        '  %s: [\n' % stepString(odbStep-restartFrom+1))
                    output.write(
                        "    (%s),    # frame 0, totaltime %g\n"
                        % ("".join(('"%s",' % els) for els in elsetList),
                           row["totaltime"]))
                    # initialize new frame
                    elsetList = list()

                else:
                    if odbFrame>lastOdbFrame or odbStep>lastOdbStep>0:
                        # write data for last frame
                        output.write(
                            "    (%s),    # frame %d, totaltime %g\n"
                            % ("".join(('"%s",' % els) for els in elsetList),
                               lastOdbFrame, row["totaltime"]))
                        # initialize new frame
                        elsetList = list()

                    if odbStep>lastOdbStep:
                        if lastOdbStep>0:
                            # write step end for last step
                            output.write("    ],\n")
                        # initialize new step
                        if restartFrom:
                            output.write(
                                '  %s: [\n' % stepString(odbStep-restartFrom+1))
                        else:
                            output.write('  %s: [\n' % stepString(odbStep))
                        output.write("    (),     # frame 0\n")

                # iteration over elsetColumns (like "Cave", "Stopes", "Drives")
                for elsetColumn in colNameList:
                    elsetColumn = elsetColumn.upper()
                    if elsetColumn not in self.elsetColumns:
                        raise ValueError(
                            "Elset column %s not recognized as elsetColum (as"
                            " given by the corresponding argument to"
                            " LoadHistory.__init__). Please update the"
                            " elsetColumns argument." % elsetColumn)

                    # store elsets for this column in elsetList for this
                    # post-type and mining step
                    try:
                        elsetList.extend(row[elsetColumn])
                        msg("In column %s found %s"
                            %(elsetColumn, row[elsetColumn]), debugLevel=15)
                    except KeyError:
                        raise KeyError(
                            "LoadHistory.writePostDataFile(): Could not find"
                            " a column '%s' in the loadhistory file %s. This"
                            " column has been specified in the postTypeCols"
                            " parameter for post set type '%s'."
                            % (elsetColumn, self.loadHistTable.fileName,
                               setType))
                    except TypeError:
                        raise TypeError(
                            "Expected lists as contents of column '%s',"
                            " instead got %s. Maybe the column '%s' is missing"
                            " in the elsetColumns parameter of"
                            " LoadHistory.__init__()."
                            % (elsetColumn, type(row[elsetColumn]).__name__,
                               elsetColumn))

                    # potentially add split elsets to elsetList for this
                    # post-type and mining step
                    try:
                        splitElsets = self.stepColToSplitElsets[
                            (rowCnt, elsetColumn)]
                    except (AttributeError, KeyError):
                        # self.stepColToSplitElsets is not defined if
                        # intersectColumns has not been called before
                        # => AttributeError
                        # self.stepColToSplitElsets may not contain an item
                        # for current column and mining step (row)
                        # => KeyError
                        pass
                    else:
                        elsetList.extend(splitElsets)
                        msg("In column %s found split elset %s"
                            %(elsetColumn, splitElsets), debugLevel=15)

                # for next iteration
                lastOdbStep = odbStep
                lastOdbFrame = odbFrame

            # add last frame in the load history
            output.write(
                "    (%s),    # frame %d, totaltime %g\n"
                % ("".join('"%s",' % els for els in elsetList),
                   lastOdbFrame, row["totaltime"]+row["miningStepLength"]))
            output.write("    ],\n  }\n\n")


        #--- write frame names
        output.write(separator)
        output.write("# frameNames\n")
        output.write("frameNames = {\n")
        frameNames = dict()
        lastOdbStep = 0
        lastOdbFrame = 0

        for rowCnt, row in enumerate(self.getLoadHistTableIter()):
            odbStep = row["odbStep"]
            odbFrame = row["odbFrame"]

            # suppress output for steps before restart-start-step
            if not(restartFrom and odbStep<restartFrom):

                if odbStep>lastOdbStep:
                    if lastOdbStep>0 and (
                            not(restartFrom) or lastOdbStep>=restartFrom):
                        # write data for last frame
                        frameNames[lastOdbStep].append(lastName)
                        output.write(
                            '    "%s",    # frame %d, totaltime %g\n'
                            % (lastName, lastOdbFrame, row["totaltime"]))
                        # write step end for last step
                        output.write("    ],\n")
                    # initialize new step
                    frameNames[odbStep] = ["Start"]
                    if restartFrom:
                        output.write(
                            '  %s: ['
                            '\n    "Start",    # frame 0, totaltime %g\n'
                            % (stepString(odbStep-restartFrom+1),
                               row["totaltime"]))
                    else:
                        output.write(
                            '  %s: ['
                            '\n    "Start",    # frame 0, totaltime %g\n'
                            % (stepString(odbStep), row["totaltime"]))
                elif odbFrame>lastOdbFrame:
                    # write data for last frame
                    frameNames[lastOdbStep].append(lastName)
                    output.write(
                        '    "%s",    # frame %d, totaltime %g\n'
                        % (lastName, lastOdbFrame, row["totaltime"]))

            # for next iteration
            lastOdbStep = odbStep
            lastOdbFrame = odbFrame
            lastName = row["miningStepName"]

        # add last frame in the load history
        frameNames[lastOdbStep].append(lastName)
        output.write(
            '    "%s",    # frame %d, totaltime %g\n'
            %(lastName, lastOdbFrame, row["totaltime"]+row["miningStepLength"]))
        output.write("    ],\n  }\n")


        #--- write frame lists
        def getFrameListOutput(frameListName, frameList):
            """
            @param frameListName: something like "frameListAll",
               "frameListHalfYearly"
            @param frameList: a list with odbStep-number, odbFrame-number
               tuples: e.g. [(2,1), (2,2), ...]

            @returns: a string like "frameListAll = [
               ("Step-2", 1), ("Step-2", 2), ...]"
            """
            output = "%s = [\n    " % frameListName
            output += ",\n    ".join(
                ", ".join("(%s,%d)" % (stepString(odbStep), odbFrame)
                          for odbStep, odbFrame in frames)
                for frames in groupIter(frameList, maxcounts=4))
            output += ",]"
            return output

        # add frame lists: frameListAll
        frameListAll = list()
        for odbStep in sorted(self.odbStepLastFrameDict):
            lastOdbFrame = self.odbStepLastFrameDict[odbStep]
            if restartFrom:
                odbStep = odbStep-restartFrom+1
                if odbStep<1:
                    continue
            frameListAll.extend(
                (odbStep, i) for i in range(1, lastOdbFrame+1))
        output.write(
            "\n# frame lists"
            "\n%s\n" % getFrameListOutput("frameListAll", frameListAll))

        # add frame lists: frameListYearly
        dateListAll = [frameNames[odbStep][odbFrame]
                       for odbStep, odbFrame in frameListAll]
        frameList = [frameListAll[i] for i in filterDateListYearly(
            dateListAll, returnDates=False, returnIndexes=True)]
        output.write(
            "\n%s\n" % getFrameListOutput("frameListYearly", frameList))

        # add frame lists: frameListHalfYearly
        frameList = [frameListAll[i] for i in filterDateListInterval(
            dateListAll, interval="H", returnDates=False, returnIndexes=True)]
        output.write(
            "\n%s\n" % getFrameListOutput("frameListHalfYearly", frameList))

        # add frame lists: frameListQuarterly
        frameList = [frameListAll[i] for i in filterDateListInterval(
            dateListAll, interval="Q", returnDates=False, returnIndexes=True)]
        output.write(
            "\n%s\n" % getFrameListOutput("frameListQuarterly", frameList))


        #--- write Abaqus frame to cavesim step relation
        #... from ATN column
        if colNameATN in self.loadHistTable:
            stepFrameToATN = self.getStepFrameToATN(colNameATN=colNameATN)
            if useStepName:
                output.write(
                    "\n# Abaqus frame to cavesim step"
                    "\nstepFrameToCsStep = {\n")
                for odbStep in sorted(stepFrameToATN):
                    frames = stepFrameToATN[odbStep]
                    output.write(
                        '    "Step-%d": [%s],\n' % (odbStep, ",\n        ".join(
                            ",".join(line) for line in groupIter(
                                frames, maxcounts=10, format="%3d")
                        )))
                    continue
                # end of dict after all steps
                output.write("    }\n")

            else:
                if restartFrom:
                    stepFrameToATN = dict(
                        (odbStep-restartFrom+1, atnList)
                        for odbStep, atnList in stepFrameToATN.iteritems()
                        if odbStep>=restartFrom)
                output.write(
                    "\n# Abaqus frame to cavesim step"
                    "\nstepFrameToCsStep = %s\n" % repr(stepFrameToATN))

#} end of post data file commands

#{ functions to accces the postData database
def readPostData(postDataFileName):
    """import the data from the specified postData database.

    @param postDataFileName: Path to the postdata file. This file usually has
        the same name as the corresponding odb with ".odb" replaced by
        "_postData.py". It should also be found in the same directory.

    @returns: a dictionary containing the data
    """
    postData = dict()
    execfile(postDataFileName, postData)
    return postData

def sequenceElsetsIter(postData, frameList=None, seqElsets="seqElsetsExcav",
                       incremental=True):
    """Generator function that yields lists of elset names that have to be
    removed or added from the last frame to the current.

    @param postData: Either a postData dictionary as supplied by readPostData
        or the path to the postdata file. This file usually has
        the same name as the corresponding odb with ".odb" replaced by
        "_postData.py". It should also be found in the same directory.

    @param frameList: (stepnumber, framenumber)-tuples specifying the frames
        for which you want output. Specifiy None to take all frames that are
        listed in the particular sequence elsets object (see seqElsets
        argument) in the postData database. None means all frames.
        Example: frameList = [ (2, i) for i in range(5, 194, 5) ]

    @param seqElsets: A string specifying the sequence elsets object in the
        postData file.
        This string must match the variable name in the postData file. I.e.
        to reference the sequence from the example postdata file sketched
        below you would specify seqElsets="seqElsets".

        The corresponding variable in the postData file must be a dictionary
        {step nb : list of elset names tuples}. step nb 1 corresponds to
        "Step-1" in the odb, step nb 2 corresponds to "Step-2" and so on.

        The list of elset names tuples contains one tuple of elsets per odb
        frame that have to be removed or added for each odb frame. The index
        in the list equals the corresponding frame number. The first item in
        the list (index 0) correponding to frame 0 of this step is an empty
        list or tuple. The second item (index 1) is a list (or tuple) of
        elset names (excluding the 'PART-1-1.'-prefix) that shall be added or
        removed from the sequence in the first frame. Third item (index 2)
        for the second frame and so on. I.e.:
          >>> seqElsets = {
          >>>     2: [    # --- step-2
          >>>         (),                                        # frame 0
          >>>         ('TUNNEL_01', 'PIT_01'), # frame 1
          >>>         ...

    @param incremental: If True the generator yields just the new elset names
        in each iteration. Otherwise all elsets from the very beginning up to
        the current frame.
    """

    # read postData database
    if isinstance(postData, basestring):
        postData = readPostData(postData)
    sequenceSets = postData[seqElsets]

    # main loop
    if frameList is None:
        frameList = [
            (step, frame)
            for step in sorted(sequenceSets)
            for frame in range(len(sequenceSets[step]))]

    firstSeqOdbStep = sorted(sequenceSets)[0]
    lastFrame = (firstSeqOdbStep, -1)
    del firstSeqOdbStep
    frameList = iter(frameList)
    resultElsets = list()
    for stepNumber, frameNumber in frameList:

        if incremental:
            resultElsets = list()
        else:
            # make a copy to not change the results already delivered
            resultElsets = list(resultElsets)

        if frameNumber>len(sequenceSets[stepNumber]):

            # iterate trough frameList to the end of this step
            # issue a warning if the last framenumber asked for is below 1000
            # (because if you ask for more than 1000 you want all)
            lastQueriedFrame = frameNumber
            while True:
                step, fr = frameList.next()
                if step != stepNumber:
                    break
                lastQueriedFrame = fr
                continue
            if lastQueriedFrame<1000:
                msg("WARNING: No data for frame %d." % frameNumber)

            # if there is another step go on with this
            # Note: this is not working properly if the first frame number for
            # the next step is higher than the actual existing frame numbers!
            if step != stepNumber:
                stepNumber = step
                frameNumber = frame
            else:  # else fini
                break

        # get sequence sets
        newFrames = list()
        # complete previous steps
        while stepNumber>lastFrame[0]:
            newFrames.extend(
                (lastFrame[0], fr)
                for fr in range(lastFrame[1]+1,
                                len(sequenceSets[lastFrame[0]])))
            lastFrame = (lastFrame[0]+1, 0)
        # add all that passed in current step
        newFrames.extend(
            (stepNumber, fr)
            for fr in range(lastFrame[1]+1, frameNumber+1))

        for step, fr in newFrames:
            # get all elsets up to current frame and remove 'PART-1-1.'
            newElsets = sequenceSets[step][fr]
            resultElsets.extend(newElsets)

        yield resultElsets
        lastFrame = (stepNumber, frameNumber)
#} end of  functions to accces the postData database

#{ other functions to parse loadhistory csv files
def sequenceElsetsIterFromCsv(
        loadhistoryCsv, frameNameColumn, elsetColumns):
    """Deliver an iterator that yields (frame name, elsetContainer)-tuples,
    one of those for each frame (identified by the frame name).

    elsetContainer is a dict with a column header of the csv file as key and
    a list of corresponding elset names as value.

    Rows are accumulated until a row is found that has some content in the
    column identified by frameNameColumn.
    """

    # make sure each elset-column is mentioned only once
    elsetColumns = sorted(set(elsetColumns))

    tab = RowIterFromCsvFile(
        loadhistoryCsv,
        columns=[frameNameColumn,]+[[col,] for col in elsetColumns],
        )

    def resetColContainer():
        return dict((col, list()) for col in elsetColumns)

    columnContainer = resetColContainer()
    for row in tab.getDictIter():

        # update columnContainer with data from current row
        for col in elsetColumns:
            # extend by anything but empty strings
            columnContainer[col].extend(e for e in row[col] if e)

        if row[frameNameColumn]:
            # return this frame
            yield row[frameNameColumn], columnContainer

            # reset columnContainer
            columnContainer = resetColContainer()
#} end of other functions to parse loadhistory csv files

#{ exportSeqSurf
def exportSeqSurfOnePerStep(
        model, postData, fileNamePrefix=None,
        frameList=None,
        sequenceSetDictName="seqElsetsExcav",
        destDir=".",
        incremental=False, outputType="iv", singleOutputFile=True,
        objectNameTemplate="F%(stepNumber)d%(frameNumber)03d_%(frameName)s",
        outputFileName=None,
        ignoreElsets=set(), ignoreElems=set(), extraSurfaces=dict()):
    """
    Create iv (or vtk or stl or raw) files representing sequenced excavations.
    One surface/mesh object will be created for each output frame (sequence
    step).

    Example:
     >>> model = Model()
     >>> model.read("../1_mesh/G01/MyProj_G01_tetL.inp")
     >>> model.read("../../2_run/R01/MyProj_G01_Q01_IncElsets.inp")
     >>> for postType in ['Drives', 'Excav', 'ExcavCave']:
     >>>     exportSeqSurfOnePerStep(
     ...         model=model,
     ...         postDataFileName="../../2_run/R01/MyProj_R01_postData.py",
     ...         fileNamePrefix=(
     ...             "%(projectPrefix)s_%(meshVersion)s_%(seqVersion)s"
     ...             % vars(config) ),
     ...         sequenceSetDictName="seqElsets%s" % postType)

    See discussion of argument outputFileName for the naming scheme.

    @param model: list of input file names or a single input file name
       or an L{bae.abq_model_02.Model}-object.

    @param postData: Either a postData dictionary as supplied by readPostData
       or the path to the postdata file. This file usually has
       the same name as the corresponding odb with ".odb" replaced by
       "_postData.py". It should also be found in the same directory.

       This file / dictionary is expected to contain frameNames and a variable
       as specified by the argument sequenceSetDictName.

    @param fileNamePrefix: Beginning of the file name for the resulting output
       files and folder name. Ignored if outputFileName is given.

       Example:
        >>> fileNamePrefix=(
        ...     "%(projectPrefix)s_%(meshVersion)s_%(seqVersion)s"
        ...     % vars(config) )

    @param frameList: (stepnumber, framenumber)-tuples specifying the frames
       for which you want output. Specifiy None to take all frames that got
       names in the postData database.

       Example: frameList = [ (2, i) for i in range(5, 194, 5) ]

       Note: It is save to specify non consecutive frames even if you
       choose incremental=True

    @param sequenceSetDictName: variable name of the dictionary in the
       postData file that contains the sequence elsets for each frame.

       Preferably those strings start with "seqElsets". If so, only the
       remainder will be used in the output filename.

    @param destDir: destination directory. Note that this directory will
       contain either the single result file or another subdirectory
       with all the files for each individual frame. destDir will *not*
       receive more than one file or directory.

       Ignored if outputFileName is given.

    @param incremental: Bool: incremental or cumulative? I.e. do you want all
       sequence up to the current frame in the corresponding output file
       (incremental=False) or do you want only those elements that were
       added between the last frame and the current (incremental=True).

       David says (11. Mar 2012): No-one in the company is allowed to make
       sequence iv's incremental (for Voxler) anymore!

    @param outputType: can be "iv" for Voxler, "vtk", "raw" or "stl"

    @param singleOutputFile: a boolean. Write all frames into one file?
       Otherwise write one file per frame to a specific subdirectory.
       Defaults to True (all in one file).

       This only applies for outputType raw, otherwise one file per frame
       will be written regardless of this parameter.

    @param objectNameTemplate: template for the object name in raw and
       vtk output and for the filename-ending in case of multiple output
       files. Note that in case of vtk output a trailing "_%(frameName)s"
       is being removed from this string for the filename template. This is for
       paraview still automatically recognizing the correct sequence order.

       Can (should?) contain the following template strings:
       "%(stepNumber)d", "%(frameNumber)03d", "%(frameName)s"

       Alternatively can be a function taking the arguments stepNumber,
       frameNumber, frameName and returns the object name to be used in raw and
       vtk output and for the filename-ending in case of multiple output files.

    @param outputFileName: optional template for the output path. If given then
       the arguments destDir and fileNamePrefix are being ignored, i.e.
       outputFileName takes precedence over those.

       This string states the full path including possible subdirectories.
       Missing directories are automatically created. The string must generally
       contain a "%s" wildcard for the time step a.k.a. object-name (see
       objectNameTemplate argument). Only if singleOutputFile==True and
       outputType=="raw" then the string must not contain any wildcards.

       Note that the given output type is governed by the outputType argument.
       It's up to the caller to supply consistent outputType and outputFileName.

       If outputFileName is not given then the default goes like this:

       In case of multiple output files a directory will be used (and created
       if necessary) with the name assembled like that:
       <fileNamePrefix>_Seq<outputType>XXX, XXX is INC for incremental or CUM
       for cumulative output.

       Individual result files will be named like this:
       <fileNamePrefix>_<sequenceSetDictName>_<frameId>.<ext>.

       If sequenceSetDictName starts with "seqElsets" (conforming to the usual
       convention) then this prefix is being removed for that purpose.
       FrameId is generated according to the objectNameTemplate argument.

       In case of a single output file (type .raw) this file will be named
       <fileNamePrefix>_<sequenceSetDictName>_RawINC.raw or ..._RawCUM.raw.

    @param ignoreElsets: set of elset names to be ignored, i.e. removed from
       the sequence for this procedure.

    @param ignoreElems: set of element numbers to be ignored, i.e. removed from
       the sequence for this procedure.

    @param extraSurfaces: A dict {(step number, frame number) : list of
       L{bae.surface_03.TriangleSurface} objects to be displayed additionally }
       Those surfaces can come from arbitrary sources like .stl files. See
       the various fromXXX class methods of L{bae.surface_03.TriangleSurface}
       for options.

    @Note: Should add optional parameters to overwrite file and directory
       names, later.
    """

    # read postData database
    if isinstance(postData, basestring):
        msg("Importing post-data from %s." % postData)
        postData = readPostData(postData)

    sequenceFrameNames = postData["frameNames"]

    # read mesh and elsets
    if isinstance(model, basestring):
        model = Model().read(model)
    elif not isinstance(model, Model):
        inputFiles = model
        model = Model()
        for fileName in inputFiles:
            model.read(fileName)
        del fileName, inputFiles

    # prepare frameList
    if frameList is None:
        frameList = [
            (step, frame)
            for step in sorted(sequenceFrameNames)
            for frame in range(len(sequenceFrameNames[step]))]

    # preprocess singleOutputFile flag
    singleOutputFile = singleOutputFile and (outputType=="raw")

    if outputFileName is None:
        # create default outputFileName if not specified as argument

        # outputTypeString will be something like "RawINC"
        outputTypeString = "%s%s" % (outputType[0].upper(), outputType[1:])
        if incremental:
            outputTypeString += "INC"
        else:
            outputTypeString += "CUM"

        # seqSetsDescr: description for file names
        # from sequenceSetDictName, e.g. "seqElsetsExcav" -> "Excav"
        seqSetsDescr = re.sub(
            r"(?i)seqElsets(.*)", r"\1", sequenceSetDictName)

        if singleOutputFile:
            outputFileName = os.path.join(
                destDir, "%s_%s_%s.%s" % (
                    fileNamePrefix, seqSetsDescr, outputTypeString, outputType))
        else:
            outputFileName = os.path.join(
                destDir, "%s_Seq%s", "%s_%s_%%s.%s") % (
                    fileNamePrefix, outputTypeString,
                    fileNamePrefix, seqSetsDescr, outputType)
            if (outputType=="vtk"
                and isinstance(objectNameTemplate, basestring)
                and objectNameTemplate.endswith("_%(frameName)s")):
                # remove trailing framenames from vtk filenames so that paraview
                # retains the sequence order
                objectNameTemplate = objectNameTemplate.rsplit("_%(frameName)s")[0]

    # check if output dir exists and eventually create it
    outputDirName = os.path.dirname(outputFileName)
    if not(os.path.isdir(outputDirName)):
        os.makedirs(outputDirName)
        msg("Created output directory %s" % outputDirName)
    else:
        msg("Output directory %s exists already." % outputDirName)

    if singleOutputFile:
        # open single output file
        output = open(outputFileName, "w")
        # in the main loop we use currentOutputFileName...
        currentOutputFileName = outputFileName

    # main loop
    for (stepNumber, frameNumber), elsetNames in izip(
            frameList, sequenceElsetsIter(
                postData, frameList=frameList,
                seqElsets=sequenceSetDictName,
                incremental=incremental)):

        # get frame name
        frameName = sequenceFrameNames[stepNumber][frameNumber]

        # filter elsets ignoreElsets and single elements ignoreElems
        if ignoreElsets:
            lenBeforeElsetFilter = len(elsetNames)
            elsetNames = set(elsetNames).difference(ignoreElsets)

        displayedElems = model.getUnionSet("elset", elsetNames)
        if ignoreElems:
            lenBeforeElemsFilter = len(displayedElems)
            displayedElems.difference_update(ignoreElems)

        # some diagnostic output
        msg("Found %d sequence elements for frame %d (%s)"
            % (len(displayedElems), frameNumber, frameName))
        if ignoreElsets and len(elsetNames) < lenBeforeElsetFilter:
            msg("... filtered %d of the %d elsets."
                % (lenBeforeElsetFilter-len(elsetNames), len(elsetNames)))
        if ignoreElems and len(displayedElems) < lenBeforeElemsFilter:
            msg("... filtered %d of the %d elements."
                % (lenBeforeElemsFilter-len(displayedElems),
                   len(displayedElems)))

        # get exteriour surface of currently displayed sets
        log.sleep()
        vol = VolumeFromMesh(model, elset=displayedElems)
        log.wakeup()
        surf = vol.getExteriorSurface()

        msg("Created surface from sequence elsets for frame %d: %s"
            % (frameNumber, surf))
        try:
            extraSurfs = extraSurfaces[(stepNumber, frameNumber)]
        except KeyError:
            pass
        else:
            for extraSurf in extraSurfs:
                surf.insertSurface(extraSurf)

        argsDict = {"stepNumber":stepNumber,
                    "frameNumber":frameNumber,
                    "frameName":frameName}
        if callable(objectNameTemplate):
            objectName = objectNameTemplate(**argsDict)
        else:
            objectName = objectNameTemplate % argsDict

        if not singleOutputFile:
            currentOutputFileName = outputFileName % objectName

        if outputType=="iv":
            surf.exportAsIv(currentOutputFileName)

        elif outputType=="stl":
            surf.exportAsSTL(currentOutputFileName)

        elif outputType=="raw":
            if singleOutputFile:
                # if singleOutputFile: write to the already open file output
                surf.exportAsRaw(output, objectName=objectName)
            else:
                surf.exportAsRaw(currentOutputFileName, objectName=objectName)

        elif outputType=="vtk":
            surf.exportAsVtk(
                currentOutputFileName,
                objectName="Sequence %s" % objectName, format='binary')

        if singleOutputFile:
            msg("Wrote object %s for frame %d" % (objectName, frameNumber))
        else:
            msg("Wrote %s file for frame %d" % (outputType, frameNumber))

    if singleOutputFile:
        output.close()
        msg("Wrote results to file %s" % outputFileName)
    else:
        msg("Finished writing output files to %s." % outputDirName)
#} end of exportSeqSurf
