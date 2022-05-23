"""Miscellaneous functions and classes
"""

__version__ = "1.65"

__bugs__ = r"""
 - none (of course)
 - way too large
"""

__todo__ = """
Some functions should be split off into other/new packages or modules:

deprecated, to be removed:
 - Tools for diagnostic output: LogFile, QuietLog, MsgTicker
 - RemoteFile (not working)
 - ftpGet, ftpPut (deprecated?)

can stay in misc_02:
 - Config (or belonging to OS tasks / job control?)
 - groupIter
 - Container
 - AutoKeyDict (?, or remove?)
 - selectattrversion (camel-case?)

related to names / numbers, can stay?
 - incrementName
 - findNextRoundNumber
 - getDateStringSeq
 - convertDateMonthToNum
 - filterDateListYearly

Spreadsheets: data from csv files
 - "tabular data" ?
 - connection to topography (old mesh_01), fieldValues (old field_01), fields

Operating system tasks: file transfer, job control
 - Config (or staying in misc_02 ?)
 - withPickle
 - needsUpdate
 - newFile
 - waitForMarker, waitForProcess
 - StaFileParser, getCurrentOdbFrame, waitForFrame, waitForStepFrame
 - ftpGet, ftpPut (deprecated?)

geometry related
 - BoundingBox
 - KDTree
 - InvDistTree
 - getShortestDistKey
 - mapSamePoints
 - findDuplicatePoints
"""

_version_history_ = r"""
1.00 GP: new
1.01 GP: added DictTableFromCsvFile, Container, BoundingBox, selectattrversion
1.02 Sven: added a function
1.03 GP: added getShortestDistKey, waitForMarker in the meantime, increased
      version number just in case I forgot before.
1.04 GP: added ftpGet
1.05 GP added KDTree, InvDistTree and BoundingBox.getDiameter()
1.06 GP modif MsgTicker accepts period=0, ftpPut accepts list of files and
      directories
1.07 GP added getCurrentOdbFrame, waitForFrame
1.08 GP added BoundingBox.getCentroid
1.09 GP added BoundingBox.getVolume
1.10 GP added DictTableFromCsv warnings if columns without header.
1.11 GP modif DictTableFromCsvFile.mergeCols() returns the merged column, no
      need to store it in self.
1.12 GP added LogFile class.
1.13 GP added DictTableFromCsvFile.writeCsv(), .insertCol(), __delitem__ and
      removeCols
1.14 GP modif: quietlog is now a subclass of LogFile;
        added RowIterFromCsvFile
1.15 GP added RowIterFromAbqReport, convertAbqReport
1.16 FR, GP added mode argument to waitForMarker
1.17 GP added renameCol to DictTableFromCsvFile;
        modif LogFile from open file, only conditionally write initial header
1.18 GP added debugLevel functionality to LogFile class
        hopefully fixed logfile="stdout" option with IDLE bug (not tested yet)
1.19 GP added incrementName(), DictTableFromCsv.findInRow()
1.20 GP added mapSamePoints()
1.21 GP modif changed DictTableFromCsvFile, RowIterFromCsvFile and
      waitForMarker to bae.log_01 functionality
1.22 GP added KDTree's ability to take a {point label : point coords} dict as
      points argument and then return point labels instead of indices.
1.23 GP added KDTree's ability to tackle arbitrary dimensional spaces
1.24 GP fixed: getCurrentOdbFrame and waitForFrame
1.25 GP added: withPickle()
1.26 GP added: needsUpdate()
1.27 GP added: waitForStepFrame() as a replacement for waitForFrame
1.28 GP added: getDateStringSeq() for a sequence of date strings
1.29 GP added: getDateStringSeq() now accepts more sorts of input,
               mapSamePoints function and docs improved
1.30 GP added: new function convertDateMonthToNum()
1.31 GP added: findNextRoundNumber()
1.32 GP added: findDuplicatePoints(), fixed mapSamePoints issue around zero
1.33 FR, GP added: Config-container
1.34 GP changed: Config runs config scripts in a top down sequence now
1.35 GP changed: Config ignores backup files like "configBE.py~"
1.36 GP: needsUpdate: missing-prerequ-exception is optional, automatically
     create missing target directories
1.37 GP: added: getDateStringSeq() can now also produce quarters and so on
1.38 GP added: now Config works without config-file, some minor cleanup
1.39 GP fixed one of the previous "some minor cleanup" was to rename the
     groupIter() format argument to format_. Reversed that for compatibility.
1.40 GP added KDTree.boxSearch
1.41 GP added newFile, BoundingBox.getSpan()
1.42 GP added keepInMem option to withPickle()
1.43 GP added AutoKeyDict
1.44 GP fixed Windows problems with newFile if bak-file already there
1.45 GP added StaFileParser.getTimeSeries()
1.46 GP added BoundingBox.pointIsInside
1.47 GP changed: DictTableFromCsvFile.insertCol() now takes multi value
        columns as well. And appends to existing columns (same column key)
        fixed: BoundingBox recognizes the dim argument correctly
1.48 AF added PVFContainer class
1.49 AF added functionality for PVFContainer
1.50 GP added DictTableFromCsvFile accepts multiline header
1.51 GP added BoundingBox.scale()
1.52 GP fixed InvDistTree weights, added radius argument to interpolate()
1.53 TR added UnstructuredField class and functions readCsv, readCsvHeader,
        dTypeMaper, toArray, arrToStructArr
1.54 TR fixed: if not hasPandas: rowCondition = None
1.55 TR renamed and updated NpContainer (formerly UnstructuredField)
1.56 GP added configDir to class Config
1.57 GP added waitForProcess
1.58 GP added getSeqTypeFromDateString
1.59 GP added InvDistTree can take vector values to interpolate
1.60 TR added isIronPython
1.61 GP added CheckFileFinished
1.62 GP added make getDateStringSeq accept negative numbers
1.63 GP added postData option to Config class
1.64 GP added Config_02
1.65 GP added added filterDateListYearly and filterDateListInterval
"""

try:
    import numpy as np
except ImportError:
    np = None
    from bae.vecmath_01 import dist

try:
    import pandas as pd
    hasPandas = True
except ImportError:
    hasPandas = False
except AttributeError as e:
    '''
    bae_utils/setup.py causes an AttributeError because there is a module named
    csv.py in bae_utils/
    '''
    if str(e).find('excel'):
        pass
    else:
        raise e


import time, os, subprocess, csv, sys, ftplib, warnings, re
import platform
import cPickle as pickle
import functools

from itertools import count as icount, izip
from math import sqrt, floor, ceil
from getpass import getuser
from socket import gethostname
from collections import OrderedDict, Callable

from bae.vecmath_01 import vector, vector_scale, vector_plus
from bae.future_01 import izip_longest, any, all, defaultdict
from bae.log_01 import msg, MsgTicker as MsgTicker_

#check module is loaded by epydoc-build
#--> epyDocBuild=True to disable decorators
if 'epydoc' in sys.modules:
    epyDocBuild = True
else:
    epyDocBuild = False


###############################################################################
###############################################################################
###############################################################################
#{ other
#

class Container(object):
    """
    This is just a class to hold arbitrary data.
    It may be initialized with named arguments.

    example:
     >>> from bae.misc_01 import Container
     >>> c = Container()
     >>> c.name = "this is c's name"
     >>> c.type = "type of c"

    would be equivalent to:
     >>> from bae.misc_01 import Container
     >>> c = Container(name="this is c's name", type="type of c")
    """
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def groupIter(input, maxcounts=8, format=None):
    """yields lists of max maxcounts elements from the "input"-iterable object

    It's a generator function. You may use in a for clause:

    >>> from bae.misc_01 import groupIter
    >>> for group in groupIter(mylist, 3):
    ...     print ", ".join(group)

    @param format: string conversion of the items of input
     - may be a formating string in that case a list of strings is returned
     - may also be a function accepting a single item of input and returning
       the formated string for this item

    @Note: There are simpler and faster solutions to (part of) this problem as
    well.
     >>> izip(*[iter(input)]*maxcounts)

     >>> from itertools import islice
     >>> flat = iter(input)
     >>> result = [list(islice(flat, maxcounts))
     ...           for i in range(len(input)/maxcounts)]

    If len(input) % maxcounts == 0 then those two chunks will do the same
    (except for the formatting). Otherwise both alternatives will
    skip incomplete groups at the end:
     - groupIter([1,2,3,4,5],2) -> [ [1,2], [3,4], [5,] ]
     - list(izip(*[iter([1,2,3,4,5])]*2)) -> [ [1,2], [3,4] ]

    If you need incomple groups at the end here is one alternative:
     >>> flat = iter(input)
     >>> n = int((len(input)-1)/maxcounts) + 1
     >>> result = [list(islice(flat, maxcounts)) for i in range(n)]

    @Note: groupIter seems to occasionally cause problems in RhinoPython with
    the debugger. If possible use one of the alternatives above.
    """
    input = iter(input)
    if not(format is None
           or isinstance(format, basestring)
           or callable(format)):
        raise TypeError(
            "groupIter does not recognize the format argument %s of type %s."
            " Only strings and functions are valid."
            % (format, type(format)))
    while 1:
        output = list()
        for i in range(maxcounts):
            try:
                val = input.next()
            except StopIteration:
                break
            if isinstance(format, basestring):
                try:
                    val = (format % val)
                except TypeError:
                    raise TypeError(
                        "groupIter could not convert %s of type %s with"
                        " the format string '%s'"
                        % (str(val), str(type(val)), format))
            elif callable(format):
                try:
                    val = format(val)
                except TypeError:
                    raise TypeError(
                        "groupIter could not convert %s of type %s with"
                        " the given format function."
                        % (str(val), str(type(val))))
            # output val
            output.append(val)

        if output:  # no empty lists at the end
            yield output
        else:
            break


def accumulate(values):
    """Cumsum-generator function
        >>> print list(accumulate([4,5,6]))
        [4,9,15]        
    """
    total = 0
    for x in values:
        total += x
        yield total

def chunkList(aList, nChunks, withStartIndex=True):
    """Splits a list in nChunks, where the last will get the rest.
    @param aList: list to be chunked
    @param nChunks: number of chunks
    @param withStartIndex: if True, a list of tuples will be returned where
        each tulpe holds the startIndex and chunk.
    """
    nPerChunk = len(aList)//nChunks
    if not nPerChunk:
        if withStartIndex:
            return [(ii,[aa,]) for ii, aa in enumerate(aList)]
        else:
            return [[aa,] for aa in aList]

    nLists = [aList[i:i+nPerChunk] for i in xrange(0, len(aList), nPerChunk)]
    if len(nLists) > 1:
        last = nLists.pop()
        nLists[-1].extend(last)
    if not withStartIndex:
        return nLists
    else:
        startIdxs = accumulate([0,] + [len(nL) for nL in nLists][:-1])
        pieces = [(sIdx, nL) for sIdx, nL in zip(startIdxs, nLists)]
        return pieces

#------------------------------------------------------------------------------

def incrementName(name, numberpattern=re.compile(r"_(\d+)$"),
                  defaultnumber="_1"):
    """Create a modified name (e.g. a file name or any named item in an Abaqus
    model) by incrementing a number contained in this name.

    I.e. "Myfile_005" -> "Myfile_006"; "BEAMORIENT" -> "BEAMORIENT_1"

    @param name: The string to be modified.
    @param numberpattern: A string or a regexp object (as returned by
       re.compile) that serves as a pattern for the search (and replace) of
       the number to be increased. In case it's a string it will be fed to
       re.compile(). The number will be extracted as the first group in the
       regexp provided, i.e.:
        >>> oldnumber = int(re.search(numberpattern, name).group(1))

    @param defaultnumber: A string. If no number has been found in name, i.e.
       numberpattern does not match name, this string will be appended to name.
    @returns: The modified string.
    """
    if isinstance(numberpattern, basestring):
        numberpattern=re.compile(numberpattern)

    match = numberpattern.search(name)
    if not match:
        newname = name+defaultnumber
    else:
        oldNb = match.group(1)
        nbLen = len(oldNb)
        try:
            newNb = int(oldNb)+1
        except ValueError:
            newNb = 1
        newname = "%s%0*d%s" % (
            name[:match.start(1)],
            nbLen, newNb,
            name[match.end(1):])

    return newname


#------------------------------------------------------------------------------

def withPickle(pickleFileName, doOnce, protocol=pickle.HIGHEST_PROTOCOL,
               keepInMem=False, dataContainer={}):
    """Service function to simplify reusing data from a previous script
    invocation. Optionally also stores a copy in memory for potential
    subsequent re-use.

    Tries to load data from a pickle file. If this is not possible (file not
    there) then create the data and store it in a newly created pickle file.

    So on first invocation withPickle() will typically just call doOnce()
    then store its result in a new pickle file and return the same result.
    On subsequent invocations if will find this pickle file and instead of
    calling doOnce again will now load the data from the pickle file.

    Simple Example:
     >>> from bae.misc_01 import withPickle
     >>> def doOnce():
     >>>     return getSomeData()
     >>> data = withPickle("myDataStore.pickle", doOnce)

    Other example:
     >>> someParameters = ...
     >>> def doOnce():
     >>>     someMorePara = ...
     >>>     result = dosomething(someParameters, someMorePara)
     >>>     result2 = ...
     >>>     return (result, result2)
     >>> result, result2 = withPickle("myDataStore.pickle", doOnce)

    Example with keepInMem option and L{needsUpdate} function and
    L{class MeshToVolMapper<bae.utils.mapping_01.MeshToVolMapper>}:
     >>> from bae.utils.mapping_01 import MeshToVolMapper
     >>> from bae.misc_01 import withPickle, needsUpdate, Config
     >>>
     >>> meshName = "MyProj_G01"
     >>> def getMapper():
     >>>     def doOnce():
     >>>         return MeshToVolMapper("%s_tetL.inp"%meshName)
     >>>     return withPickle(
     >>>         "%s_mapper.pickle"%meshName, doOnce, keepInMem=True)
     >>>
     >>> volFileNames = ["volumeMesh.inp", "volumeElsets.inp"]
     >>> resultName = "mappedElsets.inp"
     >>> if needsUpdate(targets=[resultName,], prerequisites=volFileNames):
     >>>     mapper = getMapper()
     >>>     mout = Model()
     >>>     mout.elset.update( mapper.mapElsets(volModel=volFileNames) )
     >>>     mout.write(resultName)
     >>>     msg("Wrote %d elsets to %s" % (len(mout.elset), resultName)
     >>>
     >>> if needsUpdate(targets=[...], prerequisites=[...]):
     >>>     mapper = getMapper()
     >>>     ... some more mapping ...
     >>>
     >>> if needsUpdate(targets=[...], prerequisites=[...]):
     >>>     mapper = getMapper()
     >>>     ... some more mapping ...
     >>>

    @param pickleFileName: File name for the pickle file.
    @param doOnce: A function to be called without arguments to supply
       the actual data. Will only be called if the pickle file does not yet
       exist.

       It's suggested to have this function defined immediately before the
       call to withPickle and to name it "doOnce()".

       There is intentionally no possibilty to pass arguments to doOnce()
       to promote exactly this way of using withPickle().
    @param protocol: will be passed on to pickle.dump, see there for
       explanation or take the default.
    @param keepInMem: If True then store the data (no matter if it's been taken
       from a pickle file or freshly created) for possible later reuse in
       memory. A subsequent call will then neither freshly create the object
       nor load it from the pickle file but will simply take the data stored
       in memory from the previous call.

       This is very useful in conjunction with the L{needsUpdate} function:
       If several if-needsUpdate-branches need the same expensive object
       but it's not clear which of the braches actually executes (if any at
       all) then you can use withPickle(..., keepInMem=True) to only create
       or load the required data at the first occasion it's actually needed.
       And if it's needed in a second branch again it will be there already
       at no extra cost.

       Supply keepInMem=True with every invocation. If you don't then the
       requested data will not be looked up for in memory.
    @param dataContainer: Never assign anything to this argument. It's a
       static dictionary to store all items supplied with the keepInMem
       option. DON'T TOUCH.
    """
    if keepInMem:
        try:
            return dataContainer[pickleFileName]
        except KeyError:
            pass
    try:
        pickleFile = open(pickleFileName, "rb")
    except IOError:
        data = doOnce()
        pickleFile = open(pickleFileName, "wb")
        pickle.dump(data, pickleFile, protocol)
    else:
        data = pickle.load(pickleFile)
    if keepInMem:
        dataContainer[pickleFileName] = data
    return data


#------------------------------------------------------------------------------

def needsUpdate(targets, prerequisites=None,
                raiseMissingPres=False, ignoreMissingDir=False):
    """Checks wether all files listed in targets are newer than the newest
    file listed in prerequisites. The check is done by mtime which is the
    time of the last modification of the file.

    If prerequisites is not specified then return True if any of the target
    files is missing.

    Usage example:
     >>> targets = (
     >>>     "PERSE2013_G02_Q22_elsets_MAIN.inp",
     >>>     "notFound.csv",
     >>>     )
     >>> prerequisites = (
     >>>     '../../1_mesh/G02/PERSE2013_G02_tetL_elemCentroids.csv',
     >>>     'PERSE2013_G02_Q22_Opt4LvlDrop_ElCentroids_MAIN.csv',
     >>>     )
     >>> if needsUpdate(targets, prerequisites):
     >>>     createElset_findSameElemForSeqPt_takeElCentroids(
     >>>         elCentroidsFilename=prerequisites[0],
     >>>         elemCentrElsetFname=prerequisites[1],
     >>>         seqElemNotFound=targets[1],
     >>>         outputFileName=targets[0],    )

    @param targets: list or other iterable of filenames
    @param prerequisites: optional list or other iterable of filenames
    @param raiseMissingPres: If True then raise OSError on missing
       prerequisites. Otherwise just issue an message and return False.
    @param ignoreMissingDir: If False (the default) check target directories
       and create any of them that is missing.

    @Note: In case you have more than one if-needsUpdate-branches that
       potentially need the same expensive data then have a look at
       L{withPickle}(..., keepInMem=True) to easily reuse previous
       computations.
    """
    assert hasattr(targets, "__iter__"), "targets argument must be iterable."

    # check: prereq is a list or such
    assert prerequisites is None or hasattr(prerequisites, "__iter__"), \
        "prerequisites argument (if given) must be iterable."

    if prerequisites is None or len(prerequisites)==0:
        res = not all(os.access(fn, os.F_OK) for fn in targets)
        if res:
            reason = "target missing"
        else:
            reason = "no prerequisites, target already there"
    else:

        # check: all prereqs exist
        missingPres = [fn for fn in prerequisites
                       if not(os.access(fn, os.F_OK))]
        if missingPres:
            res = False
            reason = "prerequisite(s) missing: %s" % ", ".join(missingPres)
            if all(os.access(fn, os.F_OK) for fn in targets):
                reason = "target already there, %s" % reason
            if raiseMissingPres:
                raise OSError(reason)

        # check: targets exist (otherwise: needs update)
        elif not all(os.access(fn, os.F_OK) for fn in targets):
            res = True
            reason = "target missing"

        else:
            # check: mdates
            newestPre = max(os.path.getmtime(fn) for fn in prerequisites)
            oldestTarget = min(os.path.getmtime(fn) for fn in targets)
            res = oldestTarget < newestPre
            if res:
                reason = "new prerequisite"
            else:
                reason = "target up-to-date"

    # check target directories
    if res:
        for ttt in targets:
            dirname = os.path.dirname(ttt)
            if dirname and not os.path.isdir(dirname):
                msg("Creating missing target directory %s" % dirname)
                os.makedirs(dirname)

    # write message
    if len(targets)==0:
        targetStr = ""
    elif len(targets)==1:
        targetStr = targets[0]
    else:
        targetStr = "%s and more" % targets[0]
    if res:
        msg("Creating %s (%s)" % (targetStr, reason))
    else:
        msg("Skipping creation of %s (%s)" % (targetStr, reason))

    # fini
    return res


#------------------------------------------------------------------------------

# number of items (months, quarters, ...) per year for each sequence type.
# I.e for seq type="M" returns 12, for "Q" it's 4, ...
# Primarily intended for L{getSeqTypeFromDateString}
seqTypeToMaxMonth = {
    "H": 2,
    "Q": 4,
    "B": 6,
    "M": 12}
def getSeqTypeFromDateString(dateStr, rex=re.compile(r"\d+_?(\D*?)\d+")):
    """Determine sequence type and maxMonth for given dateStr. sequence type is
    one of M, B, Q or H extracted from date strings like 2017_M12 or Y2013Q1.
    For dates not containing a valid type like Y2012_11 it's "M".

    maxMonth is the number of items (months, quarters, ...) per year. I.e for
    seqType="M" maxMonth is 12, for "Q" it's 4, ...

    See also L{getDateStringSeq} often used in conjunction with this function.

    @returns: (sequence type, maxMonth)-tuple.
    @raises ValueError: if the type could not be determined
    """
    res = rex.search(dateStr)
    if not res:
        msg("ERROR: Could not determine date type (monthly, bi-monthly,"
            " quarterly, half-yearly).")
        raise ValueError(
            "Could not determine date type (monthly, bi-monthly,"
            " quarterly, half-yearly).")
    seqType = res.group(1).upper()

    if not(seqType and (seqType in "MBQH")):
        seqType = "M"

    maxMonth = seqTypeToMaxMonth[seqType]

    return seqType, maxMonth


def getDateStringSeq(startDate, endDate, maxMonth=12):
    """Create list of dates in a format similar to "YR2012_M10".
    The format of the resulting strings mimics the format of the input
    parameters. Years are always formatted as four digit number and months
    as two-digit number with possible leading zero like "03". There must be
    a non-numeric character between year and month.

    See also L{getSeqTypeFromDateString} to automatically determine the
    maxMonth argument.

    Examples:
     >>> print getDateStringSeq("YR2012_M12", "YR2013_M03")
     ["YR2012_M12","YR2013_M01","YR2013_M02","YR2013_M03"]
     >>> print getDateStringSeq("Y2012M12", "Y2013M02")
     ["Y2012M12","Y2013M01","Y2013M02"]
     >>> print getDateStringSeq("Y2012M12", 3)
     ["Y2012M12","Y2013M01","Y2013M02"]
     >>> print getDateStringSeq("Y2012M11", [0,1,2])
     ["Y2012M11","Y2012M12","Y2013M01"]
     >>> print getDateStringSeq("Y2012M11", [2, 12])
     ["Y2013M01","Y2013M11"]
     >>> print getDateStringSeq("Y2012Q3", 3, maxMonth=4)
     ["Y2012Q3","Y2012Q4","Y2013Q1"]

    More flexible use case, with the help of L{getSeqTypeFromDateString}:
     >>> from bae.misc_01 import getSeqTypeFromDateString, getDateStringSeq
     >>> startDate = "Y2025_Q3"
     >>> seqType, maxMonth = getSeqTypeFromDateString(startDate)
     >>> newNames = getDateStringSeq(
     >>>     startDate, 5, maxMonth=maxMonth)
     >>> print seqType, newNames
     Q ["Y2025_Q3", "Y2025_Q4", "Y2026_Q1", "Y2026_Q2", "Y2026_Q3"]

    @param startDate: A string specifying the first item in the resulting
      sequence and at the same time determining the format of the output.
    @param endDate: Can be:
      - A string (formatted like startDate) specifying the last item in the
        resulting sequence.
      - An integer >1 specifying the total number of months in the sequence.
      - An integer <-1 for a sequence going backwards from startDate. The
        total number of months in the sequence is determined by the absolute
        value of endDate.
      - A sequence (iterable) of integers specifying the offset in months to
        the given startDate.
    @param maxMonth: highest possible value of the "month" number.
      Default is 12 as there are 12 months in a year. If you want quarters
      specify 4. For half-years state 2.
    @raises ValueError: If the given arguments are inconsistent.
    @returns: A list of date strings.
    """

    rex = re.compile(r"(\D*?)(\d+)(\D+?)(\d+)")

    # evaluate start date
    res = rex.match(startDate)
    if not(res):
        raise ValueError("Start date not valid %s" % startDate)
    prefYr, startYr, prefMon, startMon = res.groups()
    nbDigitYr = len(startYr)
    startYr = int(startYr)
    nbDigitMon = len(startMon)
    startMon = int(startMon)

    # prepare result
    dateFmt = "%s%%0%dd%s%%0%dd" % (prefYr, nbDigitYr, prefMon, nbDigitMon)
    dateseq = list()
    def addDate(yr, mon):
        dateseq.append(dateFmt % (yr, mon))

    # evaluate end date
    if isinstance(endDate, basestring):
        # string as end date

        res = rex.match(endDate)
        if not(res):
            raise ValueError("End date not valid %s" % endDate)
        prefYr2, endYr, prefMon2, endMon = res.groups()
        if prefYr!=prefYr2 or prefMon!=prefMon2:
            raise ValueError("Date formats inconsistent")
        endYr = int(endYr)
        endMon = int(endMon)

        # create sequence
        yr = startYr
        mon = startMon
        while yr<=endYr:
            while (yr<endYr and mon<=maxMonth) or (mon<=endMon):
                addDate(yr, mon)
                mon += 1
            mon = 1
            yr += 1

    elif isinstance(endDate, (int, float)):
        if endDate>1:
            # number of months to create
            cnt = int(endDate)  # just to rename (and make sure it's integer)

            # create sequence
            yr = startYr
            mon = startMon
            while cnt>0:
                while mon<=maxMonth and cnt>0:
                    addDate(yr, mon)
                    mon += 1
                    cnt -= 1
                mon = 1
                yr += 1
        elif endDate<-1:
            # number of months to create
            cnt = -int(endDate)  # just to rename (and make sure it's integer)

            # create sequence
            yr = startYr
            mon = startMon
            while cnt>0:
                while mon>0 and cnt>0:
                    addDate(yr, mon)
                    mon -= 1
                    cnt -= 1
                mon = maxMonth
                yr -= 1
            dateseq.reverse()
        else:
            raise ValueError(
                "ERROR in getDateStringSeq: Number of dates must be smaller"
                " than -1 or larger than 1. We got: %s" % endDate)
    else:
        # list of offsets
        for offs in endDate:
            mon = startMon+offs
            yr = startYr + int((mon-1)/maxMonth)
            mon = ((mon-1) % maxMonth) + 1
            addDate(yr, mon)

    return dateseq


def filterDateListYearly(dates, returnDates=True, returnIndexes=False):
    """Filter date strings, return only the last for each year.

    This is a very simple algorithm: It simply determines the last date of each
    year in the given list. The list of dates is assumed to be sorted
    chronologically and ascending from early to late dates.

    The date strings are assumed to conform to the following convention:
    Everything up to the first underscore "_" identifies the year. If there is
    no underscore then the whole string identifies the year.

    Arguments returnDates and returnIndexes determine the kind of output:
    If only returnDates is True then returns the filtered list of date strings.
    If only returnIndexes is True then returns a list of zero based indexes in
    the input list (dates) of the filtered dates. If both returnDates and
    returnIndexes are True then return a list of (date-string, index)-tuples.

    @param dates: input: list of date strings
    @param returnDates: Determine kind of output, see above.
    @param returnIndexes: Determine kind of output, see above.
    """
    result = list()
    lastYear = None
    lastItem = None
    if returnDates and returnIndexes:
        def store(item):
            result.append(item)
            return
    elif returnDates:
        def store(item):
            result.append(item[0])
            return
    elif returnIndexes:
        def store(item):
            result.append(item[1])
            return
    else:
        raise ValueError(
            "filterDateListYearly: At least one of the arguments returnDates"
            " and returnIndexes must be True.")

    for idx, date in enumerate(dates):
        year = date.split("_", 1)[0]
        if lastItem is not None and year != lastYear:
            # store previous date --before year changes
            store(lastItem)
        lastYear = year
        lastItem = (date, idx)
    # store last date in any case
    if dates:
        store(lastItem)
    return result


def filterDateListInterval(dates, interval,
                           returnDates=True, returnIndexes=False):
    """Filter date strings, return only steps according to the specified
    interval.

    This is a very simple algorithm: It simply determines the last date of each
    interval in the given list. The list of dates is assumed to be sorted
    chronologically and ascending from early to late dates.

    The date strings are assumed to conform to the following convention:
    Everything up to the first underscore "_" identifies the year. If there is
    no underscore then the whole string identifies the year. After the
    underscore there is a seq-type identifier "H" for half yearly, "Q" for
    quarterly, "B" for bi-monthly, "M" for monthly. If the seq-type identifier
    is missing then monthly is assumed. E.g. "Y2012_M03", "Y2014_Q2", ...

    Arguments returnDates and returnIndexes determine the kind of output:
    If only returnDates is True then returns the filtered list of date strings.
    If only returnIndexes is True then returns a list of zero based indexes in
    the input list (dates) of the filtered dates. If both returnDates and
    returnIndexes are True then return a list of (date-string, index)-tuples.

    @param dates: input: list of date strings
    @param interval: "H" for half yearly, "Q" for quarterly,
        "B" for bi-monthly, "M" for monthly
    @param returnDates: Determine kind of output, see above.
    @param returnIndexes: Determine kind of output, see above.
    """

    interval = seqTypeToMaxMonth[interval]

    result = list()
    if returnDates and returnIndexes:
        def store(item):
            result.append(item)
            return
    elif returnDates:
        def store(item):
            result.append(item[0])
            return
    elif returnIndexes:
        def store(item):
            result.append(item[1])
            return
    else:
        raise ValueError(
            "filterDateListInterval: At least one of the arguments"
            " returnDates and returnIndexes must be True.")

    lastStep = None
    lastItem = None
    rex = re.compile(r"^(.*?)_(\D?)(\d+)")

    for idx, date in enumerate(dates):
        # year is everything up to the first underscore "_"
        res = rex.search(date)
        if res:
            year = res.group(1)
            t = res.group(2)
            n = int(res.group(3))

            substep = int(ceil(float(n * interval) / seqTypeToMaxMonth[t]))
            step = (year, substep)
        else:
            year = date
            step = (year, interval)

        if lastStep is not None and step != lastStep:
            # store previous frame --before step changes
            store(lastItem)
        lastItem = (date, idx)
        lastStep = step

    # store last step in any case
    if dates:
        store(lastItem)
    return result


def convertDateMonthToNum(olddate):
    """convert spelled out months in date strings to month numbers
    eg: '2013 Feb' -> 'Y2013_M02', 'Feb13' -> 'Y2013_M02'

    Anything that looks like a number is treated as year, anything looking like
    a three letter month is taken as the month. The output is converted to
    something like Y2013_M02.

    Takes a string and returns a string.
    """
    months = ["JAN","FEB","MAR","APR","MAY","JUN",
              "JUL","AUG","SEP","OCT","NOV","DEC"]
    rex = re.compile(r"(?i)(\d*)\W*(" + "|".join(months) + r")\W*(\d*)")

    res = rex.match(olddate)
    if res:
        year1 = res.group(1)
        month = res.group(2).upper()
        year2 = res.group(3)
        month = months.index(month) + 1
        if len(year1)==4:
            newdate = "Y%s_M%02d" % (year1, month)
        elif len(year2)==4:
            newdate = "Y%s_M%02d" % (year2, month)
        elif year1:
            year = int(year1)
            if year>60:
                year += 1900
            else:
                year += 2000
            newdate = "Y%04d_M%02d" % (year, month)
        elif year2:
            year = int(year2)
            if year>60:
                year += 1900
            else:
                year += 2000
            newdate = "Y%04d_M%02d" % (year, month)
        else:
            newdate = olddate
    else:
        newdate = olddate
    return newdate


#------------------------------------------------------------------------------

def findNextRoundNumber(lastNumber, offset):
    """Find the next round number greater than lastNumber. Only integers!

    Examples:
     >>> findNextRoundNumber(304, 100)
     400
     >>> findNextRoundNumber(19800, 1000)
     20000
     >>> findNextRoundNumber(max(model.elNodes), 100000) + 1
     150001

    @param lastNumber: Last number. The result will be greater than this.
    @param offset: The result will be an integer multiple of this number.
    @return: A round number that is a multiple of offset.
    """
    return int((floor(lastNumber/float(offset))+1)*offset)


#------------------------------------------------------------------------------

class BoundingBox(list):
    """
    This is essentially a list of two points: The point [xmin, ymin, zmin]
    and the point [xmax, ymax, zmax].

    This class delivers the functionality to gather this data from some
    elements or points and to write it in a easily readable format.

    Example:
     >>> from bae.misc_01 import BoundingBox
     >>>
     >>> # bounding box object by corner points
     >>> bb = BoundingBox([xyzmin, xyzmax])
     >>> print "bounding box", bb  # using the self.__str__() method
     >>>
     >>> # bounding box of a list of points
     >>> bb = BoundingBox()
     >>> bb.update([[1,0,0], [0,10,4], [4,6,2]])
     >>> print "bounding box %s" % bb  # using the self.__str__() method
    """

    coordFormat = "%.1f"

    def __init__(self, box=None, dim=None):
        """
        @param box: initialize the new BoundingBox object with this box
        @param dim: The dimensionality of the coordinate space, usually = 3.
           If not given will be derived from the box argument and if that's
           not given it defaults to 3.
        """
        # initialize bounding box
        if box is None:
            list.__init__(self, [[],[]])
            if dim is None:
                self.dim = 3
            else:
                self.dim = dim
        else:
            list.__init__(self, box)
            # check arguments (in case box is an iterator only check after init)
            if len(self) != 2:
                raise ValueError(
                    "ERROR: BoundingBox-box-argument must be an iterable of two"
                    " point coordinate lists. One for the min values and one"
                    " for the max values. Instead we got %s" % list(box))
            if any(x1>x2 for x1, x2 in izip(*self)):
                msg("WARNING: BoundingBox-box-argument doesn't make sense,"
                    " lower bound larger than upper bound: %s" % self)
            if len(self[0]) != len(self[1]):
                msg("WARNING: BoundingBox-box-argument doesn't make sense,"
                    " different dimensions on lower vs upper corner: %s" % self)
            if dim is None:
                self.dim = len(self[0])
            else:
                if len(self[0]) != dim:
                    msg("WARNING: BoundingBox-box dim-argument (=%d) is not"
                        " consistent with box argument %s. Resulting box: %s"
                        % (dim, box, self))
                self.dim = dim

    def update(self, nodeCoords, elNodes=None):
        """
        Updates the bounding box with the given nodes/points.

        If elNodes is given it is either a list of node numbers or a dictionary
        {element id: list of node numbers}. In this case only nodes listed in
        elNodes are considered.

        nodeCoords is a dictionary {node number: coordinates tuple} defining
        all nodes. If elNodes is not specified then nodeCoords may also be
        simply a list or other iterable of coordinates tuples.
        """

        if isinstance(nodeCoords, dict):
            # if nodeCoords is a dict: consider elNodes argument
            if isinstance(elNodes, dict):
                # if elNodes is a dict: find all node numbers
                nodes = set()
                for elemNodes in elNodes.itervalues():
                    nodes.update(elemNodes)
                elNodes = nodes
                pointIter = (nodeCoords[node] for node in elNodes)
            elif elNodes is None:
                pointIter = nodeCoords.itervalues()
            else:
                pointIter = (nodeCoords[node] for node in elNodes)
        else:
            # if nodeCoords is not a dict: ignore elNodes argument
            pointIter = nodeCoords

        # iterate over all points
        dim = self.dim
        for coords in pointIter:
            try:
                for i in xrange(dim):  # x,y,z
                    if coords[i]<self[0][i]:  # min
                        self[0][i]=coords[i]
                    if coords[i]>self[1][i]:  # max
                        self[1][i]=coords[i]
            except IndexError:
                self[0]=list(coords[:dim])
                self[1]=list(coords[:dim])

    def __str__(self):
        return ("min: [%s], max: [%s]"
                % (",".join(self.coordFormat % x for x in self[0]),
                   ",".join(self.coordFormat % x for x in self[1])))

    def getDiameter(self):
        """Diameter of the smallest circle or sphere enclosing the box. Equals
        the largest possible diagonal.
        """
        return sqrt(sum((self[1][i]-self[0][i])**2 for i in xrange(self.dim)))

    def getCentroid(self):
        """Centroid of the box.
        Returns the coordinate tuple (rather list) of the centre of this box
        """
        return [0.5*(xmin+xmax) for xmin, xmax in izip(*self)]

    def getVolume(self):
        """Returns the volume of the 3D box. Or the area in 2D.
        """
        res = 1.0
        for i in xrange(self.dim):
            res *= (self[1][i]-self[0][i])
        return res

    def getSpan(self):
        """dimension of the box in each coordinate direction.
        Returns the coordinate tuple (rather list) of
        [delta_x, delta_y, delta_z] of this box.
        """
        return [(xmax-xmin) for xmin, xmax in izip(*self)]

    def scale(self, scale_factor, scale_origin=None):
        """Scale the box relative to scale_origin by scale_factor.

        Example:
          >>> # scale this box by 2 in each direction
          >>> v.scale(scale_factor=2, scale_origin=None)

        @param scale_factor: May be a single number in which case it is applied
        to all three dimensions or a list of three numbers with one factor for
        each dimension.
        @param scale_origin: The scaling in done relative to this point.
        Defaults to the centre of the box.
        """
        if scale_origin is None:
            scale_origin = vector_scale(vector_plus(*self), 0.5)

        if type(scale_factor) in (float, int):
            scale_factor = [scale_factor,scale_factor,scale_factor]

        for i in xrange(2):
            orig2old = vector(scale_origin, self[i])
            orig2new = [orig2old[0]*scale_factor[0],
                        orig2old[1]*scale_factor[1],
                        orig2old[2]*scale_factor[2]]
            self[i] = vector_plus(scale_origin, orig2new)

    def pointIsInside(self, point):
        """test if point is inside the box
        point is tuple of coordinates
        """
        return (all(p>=b for p, b in izip(point, self[0]))
                and all(p<=b for p, b in izip(point, self[1])))


#------------------------------------------------------------------------------

class DefaultOrderedDict(OrderedDict):
    """
    Joins functionality of collections defaultdict and OrderedDict
    Found here U{http://stackoverflow.com/a/6190500/562769} .
    """
    
    def __init__(self, default_factory=None, *a, **kw):
        if (default_factory is not None and
           not isinstance(default_factory, Callable)):
            raise TypeError('first argument must be callable')
        OrderedDict.__init__(self, *a, **kw)
        self.default_factory = default_factory

    def __getitem__(self, key):
        try:
            return OrderedDict.__getitem__(self, key)
        except KeyError:
            return self.__missing__(key)

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value

    def __reduce__(self):
        if self.default_factory is None:
            args = tuple()
        else:
            args = self.default_factory,
        return type(self), args, None, None, self.items()

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        return type(self)(self.default_factory, self)

    def __deepcopy__(self, memo):
        import copy
        return type(self)(self.default_factory,
                          copy.deepcopy(self.items()))

    def __repr__(self):
        return 'OrderedDefaultDict(%s, %s)' % (self.default_factory,
                                               OrderedDict.__repr__(self))



class AutoKeyDict(dict):
    """
    Dictionary with integer keys and a stored factory function

    Usage:
    ======

    >>> d = AutoKeyDict()
    >>> d.append(MyObj(3,5))
    >>> print d
    {1: MyObj(3,5)}
    >>> d.append("hallo")
    >>> print d
    {1: MyObj(3,5), 2: "hallo"}

    >>> d = AutoKeyDict(factory=list)
    >>> d.append(5,4,3)
    >>> print d
    {1: [5, 4, 3]}

    >>> d = AutoKeyDict([(3, "hallo"), (5, [1,2]) ], factory=dict)
    >>> print d
    {3: "hallo", 5: [1,2]}
    >>> d.append(a=4, b=1)
    >>> print d
    { 3: "hallo", 5: [1,2], 6: {"a":4, "b":1} }
    """
    def __init__(self, *args, **kwargs):
        try:
            self.factory = kwargs["factory"]
            del kwargs["factory"]
        except KeyError:
            self.factory = None
        dict.__init__(self, *args, **kwargs)

    def append(self, *args, **kwargs):
        "Add an entry with the last key+1"
        if len(self):
            key = max(self)+1
        else:
            key = 1

        if self.factory:
            self[key] = self.factory(*args, **kwargs)
        else:
            self[key] = args[0]

        return key

#} end of other

###############################################################################
###############################################################################
###############################################################################
#{ Tools for diagnostic output
#

# ignore RuntimeWarning concerning os.tempnam in this module!
warnings.filterwarnings(
    action="ignore", message="tmpnam", category=RuntimeWarning,
    module=__name__)


class LogFile(file):
    r"""
    DEPRECATED! Use bae.log_01 instead

    A file object for diagnostic output.
    Automatically writes some header and footer upon initialization and
    deletion.

    @ivar debugLevelMax: Only write messages to the logfile whose
       debugLevel argument is less or equal to this value. See the
       corresponding argument to self.__init__().

    @Note: Usually an instance of this class should write a final message when
       it's being deleted. However this is not guaranteed (and appears to
       usually not happen) if the object still exists when the python
       interpreter exits. So it is most save to call self.close(). (del self
       is not that sensible because there might still be another object holding
       a reference self.
    """
    doubleSep = "="*75+"\n"
    singleSep = "-"*75+"\n"

    def __init__(self, logfile="stdout", mode="w", scriptName=None,
                 writeInitString="default", debugLevelMax=0):
        """Constructor, can be used like open()

        @param logfile: If "stdout" this object will write to sys.stdout. If
        "auto" opens a file named as the current skript. If this argument is
        an already opend file object, write to this. Otherwise opens a file
        with the specified name.

        @param mode: like the corresponding argument to open(). Default "w",
        i.e. open for write / overwrite. Use "a" to append to an existing file.

        @param scriptName: For use in the header of the logfile. Defaults
        to sys.argv[0] if not specified or None. For use with abaqus cae/viewer
        use cae_01.openDefaultLogFile() to extract the correct script name.

        @param writeInitString: Force initial head lines to be written (if True)
        or skipped (if False). If not specified then some heuristics decide.

        @param debugLevelMax: Only write messages to the logfile whose
        debugLevel argument is less or equal to this value. I.e. if
        debugLevelMax==0 then no debug messages (with a debugLevel>0) will be
        printed. If debugLevelMax==5 then all debugging messages up to level 5
        will be written to this logfile. Debugging messages of level 5 are
        messages with debugLevel=5 as argument to the msg() function.
        Note: The debugLevelMax functionallity only works with the msg()
        function. self.write() writes to the logfile unconditionally.

        This argument is written to the instance variable debugLevelMax which
        may be changed at a later stage.
        """

        # name of the running script
        # for Default header in logfile and "auto" logfile name
        if scriptName is None:
            scriptName = sys.argv[0]

        if logfile=="stdout":
            logfile=sys.stdout

        # use an already existing file-like object
        if (not(isinstance(logfile, basestring))
                and hasattr(logfile, "write")):

            # write method
            self.write = logfile.write

            # flush method
            try:
                self.flush = logfile.flush
            except AttributeError:
                def doNothing():
                    pass
                self.flush = doNothing

            if writeInitString == "default":
                try:
                    # write header only to the beginning of the logfile
                    writeInitString = (logfile.tell()==0)
                except (AttributeError, IOError, ValueError):
                    # if in doubt then do write the header
                    writeInitString = True

        # create a new logfile with a filename based on the script name
        elif logfile=="auto":
            scriptBaseName = scriptName.rsplit(".",1)[0]  # w/o .py extension
            logfile = scriptBaseName+".log"
            file.__init__(self, logfile, mode)
            if writeInitString == "default":
                writeInitString = True

        # create a new logfile with the given name
        else:
            file.__init__(self, logfile, mode)
            if writeInitString == "default":
                writeInitString = True

        # initialization of the setWritePattern mechanism
        self._defaultWrite = self.write
        self._pattern = ""

        # write initial head lines to file
        if writeInitString:
            self._defaultWrite(self.doubleSep)
            self._defaultWrite(
                "script '%s' executed\n  by %s\n  in %s:%s\n  on %s\n"
                % (scriptName, getuser(), gethostname(), os.getcwd()+os.sep,
                   time.strftime("%a, %Y-%m-%d at %H:%M:%S")
                   ))
            self._defaultWrite(self.singleSep)

        # save argument
        self.debugLevelMax = debugLevelMax
        self.writeInitString = writeInitString

    def close(self):
        if self.writeInitString:
            self._defaultWrite(self.singleSep)
            self._defaultWrite("%s : Fini.\n" % time.asctime())
            self._defaultWrite(self.doubleSep)
            self.flush()
            self.writeInitString = False  # don't write anything after this

        # does nothing in case the logfile argument to __init__ is an open file
        file.close(self)

    def __del__(self):
        self.close()

    def msg(self, text, end="\n", debugLevel=0):
        r"""A convienence function to write messages to the log file.

        The command
         >>> log.msg(message)

        is identical to:
         >>> log.write("%s : %s\n" % (time.asctime(), message))
         >>> log.flush()

        @param text: message to be written, will be converted to string as
           str() does
        @param end: append this string to the text. Defaults to newline so text
           does not need to have newline included. If you don't want to finish
           the current line then specify end="".

           I.e. the following two commands result in two lines of output.
            >>> log.msg("This and that happens.")
            >>> log.msg("Now something else happens.")

           Whereas the following code will write the first message then do
           something and append the finishing message to the very same line.
            >>> log.msg("Something is going to happen ...", end="")
            >>> do_something()
            >>> log.msg(" finished. Something happened.")

        @param debugLevel: Only write this message to file if the current
           debugLevelMax setting of self is equal or larger than this argument.
           See debugLevelMax argument to self.__init__() and the instance
           variable self.debugLevelMax.
        """
        if debugLevel <= self.debugLevelMax:
            self.write("%s : %s%s" % (time.asctime(), text, end))
            self.flush()

    def _writeWithPattern(self, text):
        """this is the write function that is used, after setWritePattern() has
        been used"""
        self._defaultWrite(self._pattern % text)

    def setWritePattern(self, pattern=None):
        r"""Switch on a write pattern. Subsequent calls of self.write(text)
        effectively yield the same as if self.write(pattern % text) would have
        been called.

        This is useful to prepend a certain prefex to all output of a function.
        Suppose there is a function myworker defined in some module that you
        can't or don't want to modify:
         >>> # module myModule
         >>> def myworker(logfile):
         ...     logfile.write("worker is working.\n")
         ...     return result

        Now you want to use this function but want to prepend some dots to the
        output of myworker to your logfile:
         >>> from myModule import myworker
         >>> log = LogFile("mylog.log")
         >>> log.write("Worker start.\n")
         >>> log.setWritePattern("... %s")
         >>> result = myworker(log)
         >>> log.setWritePattern()
         >>> log.write("Worker finished.\n")
         >>> del log

        This should yield the following log output (besides the default header
        and footer of LogFile)::

        Worker start.
        ... worker is working.
        Worker finished.

        @param pattern: The pattern for the output must contain a "%s". This
        will be replaced by the argument to the write method. If pattern is not
        specified or None, the original write function is reinstalled.
        """
        if pattern is None:
            self.write = self._defaultWrite
        else:
            self._pattern = pattern
            self.write = self._writeWithPattern


#------------------------------------------------------------------------------

class QuietLog(LogFile):
    """
    DEPRECATED! Use bae.log_01 instead

    A file like object that just sends every message into the black hole
    at the centre of the galaxy whose name was also swallowed by that very
    black hole.

    Use the instance quietlog instead of the class QuietLog in script.

    >>> from bae.misc_01 import quietlog
    >>> mymodel = Model(logfile=quietlog)
    """
    def __init__(self):
        pass

    def __del__(self):
        pass

    def flush(self):
        pass

    def write(self, s):
        pass

    def close(self):
        pass

# An instance of the QuietLog class, use this in your script.
#
# >>> from bae.misc_01 import quietlog
# >>> mymodel = Model(logfile=quietlog)
quietlog = QuietLog()


#------------------------------------------------------------------------------

class MsgTicker(object):
    r"""
    DEPRECATED! Use bae.log_01 instead

    Takes frequent messages (e.g. from a loop) as input and passes them on
    to a given file like object for diagnostic output only at a specified
    intervall.

    Usage:
     >>> from bae.misc_01 import MsgTicker
     >>> logfile = sys.stdout
     >>> jobList = [...]
     >>> logfile.write("Starting to process %d jobs\n" % len(jobList))
     >>> msgTemplate = "processed %d/"+str(len(jobList))+"\n"
     >>> ticker = MsgTicker(logfile, msgTemplate)
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
     >>> msgTemplate = "processing job %d of "+str(len(jobList))+"\n"
     >>> ticker = MsgTicker(logfile, msgTemplate)
     >>> for cnt, job in enumerate(jobList):
     ...     ticker.msg(cnt+1)
     ...     job.process() # this does the real work
     >>> ticker.msgTemplate = "Finished the last job.\n"
     >>> ticker.msg()
     >>> del ticker

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

    def __init__(self, output, msgTemplate="%s",
                 firstPeriod=None, period=None,
                 printLast=False):
        """
        @param output: output stream all passing messages will be sent to
        @param msgTemplate: message template with %- conversion specifiers
          matching the arguments to the subsequent self.msg() function calls.
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
        self.output = output
        self.msgTemplate = msgTemplate
        self.printLast = printLast

        if firstPeriod is None:
            firstPeriod = MsgTicker.defaultFirstPeriod
        self.nextOutputTime = time.time()+firstPeriod

        if period is None:
            self.period = MsgTicker.defaultPeriod
        else:
            self.period = period

        self.lastArgs = None

    def msg(self, *args, **kwargs):
        """
        Specify a set of state parameters that may be used to generate
        output if enough time has passed since the last output.

        There are two possible ways to call this function depending on the
        msgTemplate argument to __init__():

        If the msgTemplate argument includes %(name)-style format options,
        the arguments to this msg() function must be keyword arguments:

        >>> ticker=MsgTicker(log, "finished %(cnt)d items. State: %(state)s\\n")
        ...    ...
        ...    ticker.msg(cnt=1, state="go")

        If the msgTemplate argument includes format options without names,
        the arguments to this msg() function must be positional arguments:

        >>> ticker=MsgTicker(log, "finished %d items. State: %s\\n")
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

            self.output.write(self.msgTemplate % args)
            self.output.flush()
            self.lastArgs = None
            self.printLast = True

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
        if self.printLast and self.lastArgs:
            self.output.write(self.msgTemplate % self.lastArgs)
            self.output.flush()

    def __del__(self):
        self.close()
#} end of Tools for diagnostic output


###############################################################################
###############################################################################
###############################################################################
#{   Spreadsheets: Data from csv files
#

class PVFContainer(object):
    """Reads data from standardized CSV file to this container. We expect
    the CSV file that is read to have the following format:

    x, y, z, var, frmId1, frmId2, ...

    with many rows, also repeating x,y,z when multiple variables exist.

    @ivar frameIds: list of frameIds = [frmId1, frmId2, ...]

    @ivar ptData: the data dictionary, with (varId,listOfValues) dictionaries
    for each point.
     >>> self.ptData = dict{
     >>>     'x_y_z': dict(  # data for pt1
     >>>         'var1': [value1, value2, ... for all frames],
     >>>         'var2': [value1, value2, ... for all frames],
     >>>         ...,
     >>>         ),
     >>>     'x_y_z': dict(  # data for pt2
     >>>         'var1': [...],
     >>>         'var2': [...],
     >>>         ...,
     >>>         ),
     >>>     ...,
     >>>     }

    @ivar ptInfo: the info dictonary, with Containers for each point.

    >>> self.ptInfo = dict{
    >>>     'x_y_z': Container(  # info for pt1
    >>>         ptCoords = (x,y,z) tuple of ptCoords,
    >>>         ptId = the id of the point as string (optional),
    >>>         matName = the matName as string (optional),
    >>>         matpropsTable = dict{ statusVal: [ (pst,ucs,e,a,s,mb)-tuples,
    >>>                                            ...], ...} (optional),
    >>>         ),
    >>>     'x_y_z': Container(  # info for pt2
    >>>         ...
    >>>         ),
    >>>     ...
    >>>     }

    @Note: The name PVFContainer stands for -> PointsVariablesFramesContainer

    @Note: It is an alternative container class compared to
    L{bae.vtk_02.FramesFieldsOnUnstructuredPoints}.
    """

    def __init__(self, dataFilename=None, infoFilename=None, crdKeyFmt="%.0f",
                 crdKeySep="_"):
        """creates a PVFContainer instance

        @param dataFilename: optional, filename of the CSV datafile.
           Defaults to None.
        @param infoFilename: optional, filename of the PY infofile.
           Defaults to None.
        @param crdKeyFmt: optional, formatting string of a single coordinate
           component.
        @param crdKeySep: optional, separating string for coordinate components.
        """
        ### store as instance variables
        self.dataFilename = dataFilename
        self.infoFilename = infoFilename
        self.crdKeyFmt = crdKeyFmt
        self.crdKeySep = crdKeySep

        ### initialize further instance variables
        self.frameIds = list()
        self.ptData = dict()
        self.ptInfo = dict()

        ### if present, then already read data
        if dataFilename is not None:
            self.read(dataFilename,infoFilename)

        return

    def read(self, dataFilename, infoFilename=None):
        """reads data from CSV datafile, and (optionally) from PY infofile.

        the @ivar ptData dictonary and @ivar ptInfo dictonary are updated.
        """
        ### read data
        fr = open(dataFilename,"rb")
        cnt = 0
        for line in fr:
            row = line.split(",")
            if cnt==0:
                self.frameIds = [id.strip() for id in row[4:]]
            else:
                ptCoords = map(float,row[:3])
                #msg("DEBUG: crds = %s"%str(ptCoords)
                xyzKey = self.getCoordKey(
                    ptCoords, self.crdKeyFmt, self.crdKeySep)
                #msg("DEBUG: xyzKey = '%s'"%xyzKey)

                try:
                    ptDataDict = self.ptData[xyzKey]
                except KeyError:
                    self.ptData[xyzKey] = dict()
                    ptDataDict = self.ptData[xyzKey]
                varKey = str(row[3])
                try:
                    varList = ptDataDict[varKey]
                except KeyError:
                    ptDataDict[varKey] = map(float,row[4:])
                else:
                    ### dont allow for overwrite, never, multiple data def
                    ### for same point and var
                    #raise Exception(
                    #    "variable '%s' for point '%s' already exists!"
                    #    % (varKey,xyzKey))

                    ### or, allow to keep existing values if differences are
                    ### small
                    oldVals = ptDataDict[varKey]
                    newVals = map(float,row[4:])
                    chk = []
                    for vn,vo in zip(newVals,oldVals):
                        try:
                            ### relative error
                            flg = abs((vn-vo)/vo)<0.01
                        except ZeroDivisionError:
                            ### absolute error
                            flg = abs(vn-vo)<0.10
                        chk.append(flg)
                    if all(chk):
                        #msg("WARNING: variable '%s' for point '%s' already"
                        #    " exists! but at least, values differ not much or"
                        #    " not at all"%(varKey,xyzKey))
                        #msg(" -> was [%s]" % (", ".join(
                        #       [ "%.1f"%v for v in ptDataDict[varKey] ])))
                        #msg(" -> now [%s]"%(", ".join(
                        #        [ "%.1f"%v for v in map(float,row[4:]) ])))
                        pass
                    else:
                        raise Exception(
                            "variable '%s' for point '%s' already exists!"
                            " differences are too much!" % (varKey,xyzKey))
                    #varList = map(float,row[4:])

                try:
                    ptInfoContainer = self.ptInfo[xyzKey]
                except KeyError:
                    self.ptInfo[xyzKey] = Container(ptCoords = ptCoords)
                    ptInfoContainer = self.ptInfo[xyzKey]
                    #msg("DEBUG: created ptInfo['%s'] Container containing"
                    #    " ptCoords = %s"%(xyzKey,str(ptCoords)))
                else:
                    ### silently overwrite ptCoords, alternatively check and
                    ### warn if not identical
                    #msg("DEBUG: overwriting ptInfo['%s'].ptCoords=%s with"
                    #    " ptCoords = %s"
                    #    %(xyzKey,str(ptInfoContainer.ptCoords),str(ptCoords)))
                    ptInfoContainer.ptCoords = ptCoords

            cnt += 1

        ### read info (optional)
        if infoFilename is not None:
            readVars = dict()
            execfile(infoFilename, readVars)
            #msg("DEBUG: readVars = '%s'"%str(readVars))
            #assert(len(readVars)==1)
            #varName = readVars.keys()[0]
            varName = "ptsInfo"
            ptInfo = readVars[varName]
            assert(isinstance(ptInfo, dict))
            for k,v in ptInfo.iteritems():
                assert(isinstance(v,Container))
                try:
                    ptInfoContainer = self.ptInfo[k]
                except KeyError:
                    ### create new container
                    self.ptInfo[k] = v
                    #msg("DEBUG: created new Container()")
                else:
                    ### update variables of existing container
                    #msg("DEBUG: updating existing Container()")
                    for na, val in vars(v).iteritems():
                        #msg("DEBUG: name='%s', value='%s'"%(str(na),str(val)))
                        setattr(ptInfoContainer, na, val)
        return

    def updatePtInfoWithPtCoordsDict(self, ptCoordsDict):
        """
        @param ptCoordsDict: a dictionary with ptId as key and
            ptCoords (tuple of three coordinates) as value.
        """
        for ptId,ptCoords in ptCoordsDict.iteritems():
            ptKey = self.getCoordKey(ptCoords)
            try:
                ptInfo = self.ptInfo[ptKey]
            except KeyError:
                self.ptInfo[ptKey] = Container(ptId=ptId, ptCoords=ptCoords)
            else:
                ptInfo.ptId = ptId
                ptInfo.ptCoords = ptCoords
        return

    def getCoordKey(self, ptCoords, crdKeyFmt=None, crdKeySep=None):
        """returns a string from ptCoords to be used as dictionary-key

        @param ptCoords: a list of point coordinates
        @param crdKeyFmt: optional, point coordinate formatting string.
            defaults to @ivar self.crdKeyFmt
        @param crdKeySep: optional, point coordinates seperating string.
            defaults to @ivar self.crdKeySep
        """
        if crdKeyFmt is None:
            crdKeyFmt = self.crdKeyFmt
        if crdKeySep is None:
            crdKeySep = self.crdKeySep
        return crdKeySep.join(crdKeyFmt % c for c in ptCoords)

    def write(self, outDataFilename, outInfoFilename=None):
        """writes the data for all points to a CSV datafile, and also the
        points information for all points to a PY infofile.

        @param outDataFilename: the filename for the CSV datafile
        @param outInfoFilename: optional, the filename for the PY infofile.
            defaults to None, then, derived from outDataFilename
        """

        ### ----- def writeData(self, outDataFilename): -----
        fw = open(outDataFilename, "wb")

        ### write header
        fw.write("x,y,z,var,"+",".join(self.frameIds)+"\n")

        ### loop over all points
        for ptKey,ptDataDict in self.ptData.iteritems():
            try:
                ptCoords = self.ptInfo[ptKey].ptCoords
            except KeyError:
                ptCoords = map(float, ptKey.split(self.crdKeySep))
            ptXYZ = ",".join(str(c) for c in ptCoords)

            ### loop over all variables for this point
            for varId in self._defaultVarsList(ptDataDict.keys()):
                vals = ptDataDict[varId]
                fw.write(
                    ptXYZ+","+varId+","
                    + ",".join(str(v) for v in vals) + "\n")

        fw.close()
        del fw

        ### ----- def writeInfo(self, outInfoFilename): -----
        if outInfoFilename is None:
            outInfoFilename = outDataFilename.replace("_data.csv","")+"_info.py"

        fw = open(outInfoFilename, "wb")
        fw.write("\n".join([
            "from bae.misc_01 import Container",
            "",
            "ptsInfo = {",
            "    # ptCoordsKey: Container( ",
            "    #                  ptCoords = <tuple of ptCoords (x,y,z)>,",
            "    #                  ptId = <ptId as string> (optional),",
            "    #                  matName = <matName string>,",
            "    #                  matpropsTable = dict(",
            "    #                      <status as int>: <list of tuples of"
            " (pst, ucs, e, a, s, mb)>,",
            "    #                      ...,",
            "    #                      ),",
            "    #                  ),",
            #"    # ...",
            ])+"\n")

        ### sorted by ptKey
        #sortedPtKeys = sorted(self.ptInfo.keys())
        ### sorted by ptId
        sortedPtIds,sortedPtKeys = zip(*sorted(
            [(self.ptInfo[k].ptId,k)
             for k in self.ptInfo.iterkeys()]))

        for ptKey in sortedPtKeys:
            ptCont = self.ptInfo[ptKey]
            fw.write("\n".join([
                "    '%s': Container(" % ptKey,
                ]+[
                "        %s = '%s'," % (k,v)
                for k,v in vars(ptCont).iteritems()
                if isinstance(v, basestring)
                ]+[
                "        %s = %s," % (k,str(v))
                for k,v in vars(ptCont).iteritems()
                if not isinstance(v, basestring)
                ]+[
                "        ),"
                ])+"\n")
        fw.write("   }\n")
        fw.close()
        del fw

        return

    def writeSinglePoint(self, ptKey, outDataFilename, outInfoFilename=None,
                         varsList=None):
        """writes the data for a particaular point to a CSV datafile,
        which has the format of a L{bae.misc_01.DictTableFromCsvFile}.

        @param ptKey: the key of the point that should be written.
        @param outDataFilename: the filename for the CSV datafile.
        @param outInfoFilename: optional, the filename for the PY infofile.
            defaults to None, then, derived from outDataFilename
        @param varsList: optional, list of varIds to be written in given
            soring of this list. defaults to None, then all varIds are
            written in a standardized (non-alphabetically) sorting order.
        """

        ptDataDict = self.ptData[ptKey]
        if varsList is None:
            varsList = self._defaultVarsList(ptDataDict.keys())

        ### ----- def writeDataSinglePoint(self, outDataFilename): -----
        ### NOTE: This is classic write to file. Alternatively, one could
        ###       simply create an instance of class DictTableFromCsvFile
        ###       and write this to file.
        fw = open(outDataFilename, "wb")

        ### write header
        fw.write("frameName,"+",".join(varsList)+"\n")

        ### loop over all frames
        for frmIdx,frmName in enumerate(self.frameIds):

            ### values for this frame for all (requested) variables
            frmVals = [ptDataDict[varId][frmIdx] for varId in varsList]
            fw.write(frmName+"," + ",".join(str(v) for v in frmVals) + "\n")

        fw.close()
        del fw

        ### ----- def writeInfoSinglePoint(self, outInfoFilename): -----
        if outInfoFilename is None:
            outInfoFilename = outDataFilename.replace("_data.csv","")+"_info.py"

        ptInfoContainer = self.ptInfo[ptKey]

        fw = open(outInfoFilename, "wb")
        fw.write("\n".join(
            "%s = '%s'" % (k,v)
            for k,v in vars(ptInfoContainer).iteritems()
            if isinstance(v, basestring))+"\n")
        fw.write("\n".join(
            "%s = %s" % (k,str(v))
            for k,v in vars(ptInfoContainer).iteritems()
            if not isinstance(v, basestring))+"\n")
        fw.close()
        del fw
        return

    def _defaultVarsList(self, availableVars):
        varsList = []
        remainingVars = set(availableVars)
        idextns = ["", "_ENVMIN", "_ENVMAX"]
        ### two versions (for copy paste):
        ### 1) iterate over ids first, then over idextns:
        ###    --> U1,U2,U3, U1_ENVMIN,U2_ENVMIN,U3_ENVIM,
        ###                  U1_ENVMAX,U2_ENVMAX,U3_ENVMAX
        #varsList.extend([ id+idextn for idextn in idextns
        #                  for id in ids if id+idextn in remainingVars ])
        ### 2) iterate over idextns first, then over ids:
        ###    --> U1,U1_ENVMIN,U1_ENVMAX, U2,U2_ENVMIN,U2_ENVMAX,
        ###        U3,U3_ENVMIN,U3_ENVMAX,
        #varsList.extend([id+idextn for id in ids
        #                 for idextn in idextns if id+idextn in remainingVars])

        ### displacements 'U[R][_]1,2,3'
        #-- abaqus notation
        ids = ['U_1', 'U_2', 'U_3',]
        varsList.extend(id_+idextn for idextn in idextns
                        for id_ in ids if id_+idextn in remainingVars)
        remainingVars.difference_update(set(ids))
        ids = ['UR_1', 'UR_2', 'UR_3',]
        varsList.extend(id_+idextn for idextn in idextns
                        for id_ in ids if id_+idextn in remainingVars)
        remainingVars.difference_update(set(ids))
        #-- user notation
        ids = ['U1', 'U2', 'U3',]
        varsList.extend(id_+idextn for idextn in idextns
                        for id_ in ids if id_+idextn in remainingVars)
        remainingVars.difference_update(set(ids))
        ids = ['UR1', 'UR2', 'UR3',]
        varsList.extend(id_+idextn for idextn in idextns
                        for id_ in ids if id_+idextn in remainingVars)
        remainingVars.difference_update(set(ids))

        ### stress components 'S11', 'S22', 'S33', 'S12', 'S23', 'S13',
        # abaqus notation
        ids = ['S_11', 'S_22', 'S_33', 'S_12', 'S_23', 'S_13',]
        varsList.extend(id_+idextn for idextn in idextns
                        for id_ in ids if id_+idextn in remainingVars)
        remainingVars.difference_update(set(ids))
        # user notation
        ids = ['S11', 'S22', 'S33', 'S12', 'S23', 'S13',]
        varsList.extend([id_+idextn for idextn in idextns
                         for id_ in ids if id_+idextn in remainingVars])
        remainingVars.difference_update(set(ids))

        ### principal stresses 'S1', 'S2', 'S3',
        # abaqus notation
        ids = ['S_%s_PRINCIPAL' % v for v in ['MIN','MID','MAX']]
        varsList.extend(id_+idextn for idextn in idextns
                        for id_ in ids if id_+idextn in remainingVars)
        remainingVars.difference_update(set(ids))
        # user notation
        ids = ['S1', 'S2', 'S3',]
        varsList.extend(id_+idextn for idextn in idextns
                        for id_ in ids if id_+idextn in remainingVars)
        remainingVars.difference_update(set(ids))

        ### computed principal stresses,
        ### for S1,S2,S3: '_MAG', '_DIR1','_DIR2','_DIR3', '_BRG','_PLG'
        ids = ['S%d%s' % (i,c) for i in range(1,4)
               for c in ['_MAG', '_DIR1','_DIR2','_DIR3', '_BRG','_PLG']]
        varsList.extend(id_ for id_ in ids if id_ in remainingVars)
        remainingVars.difference_update(set(ids))

        ### computed stress invariants: p,q,th
        ids = ['INV1_p', 'INV2_q', 'INV3_th',]
        varsList.extend(id_ for id_ in ids if id_ in remainingVars)
        remainingVars.difference_update(set(ids))

        ### status, plastic strain: STATUS, PST,LOGP
        ids = ['STATUS',]
        varsList.extend(id_ for id_ in ids if id_ in remainingVars)
        remainingVars.difference_update(set(ids))
        ids = ['PST', 'LOGP']
        varsList.extend(id_+idextn for id_ in ids
                        for idextn in idextns
                        if id_+idextn in remainingVars)
        remainingVars.difference_update(set(ids))

        ### any material parameters: mpXXX
        ids = sorted(id_ for id_ in remainingVars if id_.startswith("mp"))
        varsList.extend(ids)
        remainingVars.difference_update(set(ids))

        ### any other
        ids = sorted(id_ for id_ in remainingVars)
        varsList.extend(ids)
        #remainingVars.difference_update(set(ids))

        return varsList


class IgnoredColumns(Exception):
    """Exception indicating that data columns without corresponding header
    field have been found.
    """
    pass


class DictTableFromCsvFile(dict):
    r"""
    An object representing tabular data with named columns and numbered rows.

    Reads a csv file and treats the first row containing something as
    column headers. Columns without value in the header row are ignored. If the
    same header value appears more than once those columns are treated as one
    and the corresponding values list will contain one list per row. The latter
    consists of as many values as this particular headline values appeared in
    the headline.

    Empty rows before the header line and between the header line and the
    first data row (the first row containing something after the header
    line) are ignored. Trailing empty lines are ignored as well.
    Intermediate empty rows are not ignored! They result in a row of None
    values.

    Objects of this class are basically dictionaries. Keys are the column
    headers from the first (nonempty) row. Values are lists with the values
    from the csv file or lists of lists of values in the case of multiple
    identical column headers. Empty cells result in a None value. Numbers
    are converted to floats even if the number was surrounded by quotes.

    To retrieve a column use the column header (a string value) as key. To
    retrieve a row use the index (an integer) as key. In the latter case you
    will get a dictionary with the column headers as keys. As values you will
    get the values of the row corresponding to the integer key, i.e. strings or
    floats or list in case of a column that appeared more than once in the
    header line.

    Suppose we've got the file mytable.csv:
     >>> # just list the contents of the input csv file:
     >>> print "".join(open("mytable.csv").readlines()) # "col01" appears twice!
     "col01","col02","col03","col01"
     "lala",1,2.5,"blub"
     10,"5",,"oops"
     ,-4,,3
     1e6,1.0E001,0.2,"bah"
     ...

    What we can do with this class is:
     >>> from bae.misc_01 import DictTableFromCsvFile
     >>> tab = DictTableFromCsvFile("mytable.csv")
     >>> print tab.keys()
     ["col01", "col02", "col03"]
     >>> print tab["col02"]     # returns a list of all items of that column
     [1.0, 5.0, -4.0, 10.0, None, -1.0, -5.0, 4.0]
     >>> print tab["col01"]     # if the same header occurs more than once:
     ...                        # returns a list of lists
     [['lala', 'blub'], [10.0, 'oops'], [None, 3.0], [1.0e6, 'bah'], ... ]
     >>> print tab["col02"][5]
     -1.0
     >>> print tab[5]["col02"]  # generally not recommended because it's slower
     -1.0
     >>> row = tab[0]     # query row with the row index: returns a dictionary
     >>> print row["col02"]
     1.0
     >>> print row["col01"]
     ['lala', 'blub']
     >>> for row in tab:  # can be used as an iterator over the rows
     ...     print row["col02"]
     1.0
     5.0
     -4.0
     ...

    We can also add, modify and delete columns and export to a csv file. For
    adding use only insertCol! If modifying/replacing a column make sure not to
    change the length of the list representing this column. Don't change a
    list of lists ("multicol") to a list of values ("singlecol") or vice versa.
    Don't change the length of multicol items (the length of the lists being
    items in the column list).
    ... And it goes like this:
     >>> from bae.misc_01 import DictTableFromCsvFile
     >>> tab = DictTableFromCsvFile("mytable.csv")
     >>> tab.insertCol("extracol", ["blublub", "oopoops", 6, "babah", ...])
     >>> tab["col01"][3]=[100,"buh"]  # modify a single value, was [1e6,'bah']
     >>> col2 = tab["col02"]          # just an abreviation
     >>> for i in range(len(col2)):   # modify values one by one in place
     >>>     if col2[i] is None:
     >>>         col2[i] = 0.0
     >>> # replace whole column (don't change length of multicol columns!)
     >>> tab["col02"] = [x*2 for x in tab["col02"]]
     >>> del tab["col03"]   # remove entire column
     >>> del tab[1]         # remove entire row from all columns
     >>> tab.writeCsv("myNewTable.csv")  # now save the whole stuff.

    @note: If you want to access a specific item with column id and row number
    use the form self[col][row]. This is faster than the other possibilty
    self[row][col].

    @ivar headline: headline as supplied by the csv reader. I.e. a list of
       strings stating the items in the first (non empty) line of the file.
    @ivar columnKeys: headline with items converted according to the
       replaceSpaceBy argument of the constructor. I.e. This lists the
       keys in self corresponding to each column in the csv file.
    @ivar colValListIndex: a list; for each column in the csv file
       (corresponding to the entries in self.headline and self.columnKeys)
       contains the index in the corresponding values-list if this value is
       a multi column value (if the corresponding header value appeared
       multiple times). If the corresponding value is a single column value
       (the corresponding header value appeared only once) then the
       corresponding item in this list is None.
    @ivar fileName: file name (mostly for diagnostic output)
    """

    def __init__(self, csvFile, replaceSpaceBy=None, logfile=None,
                 warnIgnoreCols=False, myOwnHeader=None,
                 headLines="auto", replaceHeaderFeedBy="_", **kwargs):
        r"""
        Constructor

        Besides those listed below you may specify additional keyword
        arguments. They will be passed to the csv.reader() function. See
        U{http://docs.python.org/library/csv.html#csv-fmt-params} for possible
        format parameters or use the parameter dialect to specify a certain
        flavour of csv file.

        Example, columns separated by tabs:
         >>> tab = DictTableFromCsvFile("mytable.csv", delimiter="\t")

        @param csvFile: either a filename or an already open file that
          contains the table data
        @param replaceSpaceBy: If not None, all spaces in each column header
          will be replaced by this string.
        @param logfile: DEPRECATED argument being ignored. use the bae.log_01
          functionality for diagnostic output
        @param warnIgnoreCols: If False data in columns without corresponding
          header field are silently ignored. Otherwise you get some kind of
          warning: If True or "message" there will be a warning message.
          If "exception" a IgnoredColumns Exception will be raised
        @param myOwnHeader: For reading csv files without a header line you can
          specify a would-be-header. I.e. a string containing comma-separated
          column names. Note: If specified then all lines (including the very
          first) will be parsed as data lines.
        @param headLines: If None or zero or False, no head line expected. All
          lines will be treated as data lines. If a positive integer N then the
          first N lines will be treated as head lines. If headLines=="auto" then
          the first non empty line will be used, following empty lines will be
          skipped.
        @param replaceHeaderFeedBy: Seperator between the contents of a
          header field on multiple rows. This conversion will be performed to
          get the keys for the data in the specified column.
        """

        # initialize the dict
        dict.__init__(self)

        # initialize the csv file reader
        if isinstance(csvFile, basestring):
            self.fileName = csvFile  # file name for diagnostic output
            csvFile = open(csvFile, "rb")
        else:
            # find a file name for diagnostic output
            try:
                self.fileName = csvFile.name
            except AttributeError:
                self.fileName = "input of type %s" % type(csvFile).__name__

        # initialize the table
        csvTable = csv.reader(csvFile, **kwargs)

        # find the first line that is not empty
        headline = None
        if myOwnHeader:
            if headLines and headLines!="auto":
                try:
                    N = int(headLines)
                except (TypeError, ValueError):
                    msg("WARNING: DictTableFromCsvFile argument myOwnHeader"
                        " specified. Ignoring argument headLines='%s'!"
                        % headLines)
                else:
                    # skip first headLines lines
                    msg("WARNING: DictTableFromCsvFile arguments myOwnHeader"
                        " and headLines specified. Skipping the first %d lines"
                        " ignoring its content." % N)
                    for i in N:
                        csvTable.next()
            headline = myOwnHeader
        else:
            if headLines=="auto":
                # find the first line that is not empty
                headline = None
                for headline in csvTable:
                    if len(headline) and any(map(len, headline)):
                        break
                # check if file is empty
                if headline is None or len(headline)==0:
                    # file empty, no headline
                    raise ValueError(
                        "DictTableFromCsvFile.__init__(): The specified csv"
                        " file %s is empty." % self.fileName)
            elif headLines>0:
                try:
                    lines = [csvTable.next() for cnt in range(headLines)]
                except StopIteration:
                    raise ValueError(
                        "DictTableFromCsvFile.__init__(): The specified csv"
                        " file %s does not contain the specified number (%d)"
                        " of header lines." % (self.fileName, headLines))
                # in python 2.4 you may not specify additional arguments after
                # a *args construct. The following is valid in python 2.6:
                # headline = ["\n".join(x)
                #             for x in izip_longest(*lines, fillvalue="")]
                # ... and here is the replacement for earlier versions
                headline = ["\n".join(x)
                            for x in izip_longest(*lines, **{"fillvalue":""})]
            else:
                raise ValueError(
                    "Illegal argument headLines to DictTableFromCsvFile"
                    " constructor: %s" % headLines)

        # permanently store headline (in its original form)
        self.headline = headline

        # convert headline according to replaceSpaceBy and replaceHeaderFeedBy
        if replaceSpaceBy is not None:
            headline = [item.replace(" ", replaceSpaceBy)
                        for item in headline]
        if replaceHeaderFeedBy is not None:
            headline = [item.replace("\n", replaceHeaderFeedBy)
                        for item in headline]

        # evaluate headline
        keyColNb = dict()
        recognizedColsList = set()
        self.colValListIndex = [None]*len(headline)
        self.columnKeys = [None]*len(headline)
        for cnt, item in enumerate(headline):
            if item!="":
                self.columnKeys[cnt] = item
                try:
                    # this will only work for the third and subsequent
                    # identical column headers
                    keyColNb[item].append(cnt)
                    self.colValListIndex[cnt] = len(keyColNb[item])+1
                except KeyError:
                    # first item of that name: just store the column number
                    keyColNb[item] = cnt
                except AttributeError:
                    # second item of that name: store a list of column numbers
                    lastColOfThisKey = keyColNb[item]
                    keyColNb[item] = [lastColOfThisKey, cnt]
                    self.colValListIndex[lastColOfThisKey] = 0
                    self.colValListIndex[cnt] = 1

                recognizedColsList.add(cnt)

        # initialize resulting dict
        for key in keyColNb:
            self[key] = []
        if warnIgnoreCols:
            ignoredCols = set()

        # loop over each row / line
        first = True
        nbEmptyRowsAtEnd = 0
        for row in csvTable:

            if len(row)==0 or not(any(map(len, row))):
                # ignore empty lines at the beginning
                # and count trailing empty lines (so they can be removed later)
                nbEmptyRowsAtEnd += 1
                if first:
                    continue
            else:
                # no trailing empty lines so far (current line is not empty)
                nbEmptyRowsAtEnd = 0
                # the following line is definitely not the first
                first = False

            for key, colNbList in keyColNb.iteritems():
                if isinstance(colNbList, list):
                    # multi column value: append a new list to the column list
                    # in the following loop append values to this new list
                    valueList = list()
                    dict.__getitem__(self, key).append(valueList)
                else:
                    # in the following loop append values directly to the
                    # column list
                    valueList = dict.__getitem__(self, key)
                    # only execute the loop once (and make colNbList a list)
                    colNbList = [colNbList]
                for colNb in colNbList:
                    try:
                        val = row[colNb]
                    except IndexError:
                        # if that column does not exist on this row
                        val = None
                    else:
                        # a column colNb exists
                        if val=="":
                            # ... but it's empty
                            val = None
                        else:
                            # there is actually a value in this column
                            try:
                                val = float(val)
                            except ValueError:
                                pass
                    valueList.append(val)

            # check wether there are values in not recognized columns
            if warnIgnoreCols:
                nonEmptyCols = set(
                    colNb for colNb, val in enumerate(row) if val!="")
                ignoredCols.update(nonEmptyCols.difference(recognizedColsList))

        # ignored columns, warning?
        if warnIgnoreCols and len(ignoredCols):
            msgText=("The following %d columns have data but no header:\n%s"
                     % (len(ignoredCols), sorted(ignoredCols)))
            if warnIgnoreCols=="exception":
                raise IgnoredColumns(msgText)
            else:
                msg("WARNING: %s" % msgText)

        # at last remove trailing empty rows from all columns
        if nbEmptyRowsAtEnd>0:
            for values in self.itervalues():
                del values[-nbEmptyRowsAtEnd:]

        return

    def __getitem__(self, key):
        """If key is a string then retrieve the corresponding column. If key is
        an integer return a dictionary containing the same keys as self with
        the values of the row corresponding to row number key.

        @note: If you want to access a specific item with column id and row
        number use the form self[col][row]. This is faster than the other
        possibilty self[row][col].
        """
        if isinstance(key, basestring):
            return dict.__getitem__(self,key)
        elif isinstance(key, int):
            return dict(
                (colKey, column[key])
                for colKey, column in self.iteritems())
        else:
            raise TypeError("Key must be either int or string.")

    def __iter__(self):
        for i in icount():
            try:
                yield self[i]
            except IndexError:
                break

    def __delitem__(self, key):
        """remove specified column or row
        @param key: the key (a string) of a column to be removed or a row index
           (an integer) of the row to be removed.

        @Note: Consider using self.removeCols if you want to remove many
           columns.
        """
        if isinstance(key, basestring):
            dict.__delitem__(self, key)
            removeIds = (idx for idx,key2 in enumerate(self.columnKeys)
                         if key == key2)
            for idx in removeIds:
                del self.columnKeys[idx]
                del self.headline[idx]
                del self.colValListIndex[idx]
        elif isinstance(key, int):
            for col in dict.itervalues(self):
                del col[key]
        else:
            raise TypeError("Key must be either int or string.")

    def getColumnLength(self):
        """Return the length of (supposedly all) columns of the table."""
        return len(dict.itervalues(self).next())

    def getRowColToList(self, rowIdx, listColumns=[],
                        removeNones=False):
        r"""
        Returns a dictionary {column header: column value of the specified row}
        exactly as the expression self[i] with an integer i does. Different
        from the latter this function makes sure, all values for the keys
        listed in the argument listColumns are lists.

        @param listColumns: identifies certain columns as possibly occuring
        multiple times and therefore making sure the resulting item will always
        be a list, no matter how often that column actually occurred in the
        csv file.
        If a column listed in listColumns occured multiple times in the columns
        header line the result is exactly the same as self[i], listColumns has
        no effect on that column in that case. If that column occured only
        once, a list with a single value is returned for that particular item.
        If there was no column with that name the resulting dictionary will
        contain a corresponding item with an empty list.

        @param removeNones:
        If the argument removeNones evaluates to True, the resulting lists for
        columns in listColumns will not contain None values. That might result
        in different lengths of the lists for different rows and the same
        column. Columns not in listColumns are not touched, they still may
        contain None values.

        E.g. suppose the following table

        >>> # just list contents:
        >>> print "".join(open("mytable.csv").readlines())
        "col01","col02","col03","col01"
        "lala",1,2.5,"blub"
        10,"5",,"oops"
        ,-4,,3
        >>> tab = DictTableFromCsvFile("mytable.csv")

        When retrieving a specific row usually columns occuring once yield a
        single value (a string or a float) whereas columns occuring multiple
        times yield a list values:

        >>> print tab[0]
        {'col01': ['lala', 'blub'], 'col02': 1.0, 'col03': 2.5}

        With getRowColToList() you can make sure you allways get lists of
        values for certain columns:

        >>> print tab.getRowColToList(rowIdx=0, listColumns=["col01","col02"])
        {'col01': ['lala', 'blub'], 'col02': [1.0], 'col03': 2.5}

        For the removeNones argument: Note that the last data row ",-4,,3"
        has some empty fields which would usually be imported as None values.
        The optional argument removeNones may be used to suppress Nones in
        the columns listed in the listColumns argument. Note that columns
        not listed in listColumns are not affected by this option.

        >>> print tab.getRowColToList(rowIdx=2, listColumns=["col01","col03"])
        {'col01': [None, 3], 'col02': -4.0, 'col03': [None]}
        >>> print tab.getRowColToList(rowIdx=2, listColumns=["col01","col03"],
        ...     removeNones=True)
        {'col01': [3], 'col02': -4.0, 'col03': []}
        >>> print tab.getRowColToList(rowIdx=2, removeNones=True)
        {'col01': [None, 3], 'col02': -4.0, 'col03': None}
        """
        if len(self)==0:
            # Note: IndexError is what you'd get from the line
            # "item = column[rowIdx]" below if you'd ask for a row that does
            # not exist. So it makes sense to also raise an IndexError if there
            # is no data at all in self.
            raise IndexError(
                "DictTableFromCsvFile.getRowColToList: There is no row %d"
                " because this DictTableFromCsvFile-object is empty, i.e. has"
                " no columns." % rowIdx)
        listColumns = set(listColumns)  # this is meant to speed up the search
        resultRow = dict()
        for colKey, column in self.iteritems():
            if colKey in listColumns:
                item = column[rowIdx]
                if isinstance(item, list):
                    resultItem = item
                else:
                    resultItem = [item]
                if removeNones:
                    resultItem = [i for i in resultItem if i is not None]
            else:
                resultItem = column[rowIdx]
            resultRow[colKey] = resultItem
        return resultRow

    def rowIterColToList(self, listColumns=[], removeNones=False):
        r"""
        Returns an iterator that returns a dictionary {column header: column
        value of the specified row} as from self.getRowColToList() for all rows
        one after the other. Makes sure that all column values for the keys
        listed in the argument listColumns are lists.

        See the documentation for self.getRowColToList().
        """
        for i in icount():
            try:
                yield self.getRowColToList(i, listColumns, removeNones)
            except IndexError:
                break

    def mergeCols(self, sourceCols, newColName=None):
        r"""
        Add a new column named newColName to self. Its items are lists
        containing the corresponding items of the sourceCols.

        E.g.:
         >>> tab = DictTableFromCsvFile("mydata.csv")
         >>> print tab.keys()
         ['x', 'y', 'z', 'value']
         >>> print tab["x"], tab["y"], tab["z"]
         [1.0, 5.0, -4.0] [None, -1.0, -3.0] [1.0, -5.0, 4.0]
         >>> tab.mergeCols(["x", "y", "z"], "coords")
         >>> print tab.keys()
         ['x', 'y', 'z', 'value', 'coords']
         >>> print tab["coords"]
         [[1.0, None, 1.0], [5.0, -1.0, -5.0], [-4.0, -3.0, 4.0]]

        And with numpy (possibly it doesn't like None values in the table):
         >>> from numpy import array
         >>> tab = DictTableFromCsvFile("mydata.csv")
         >>> tab.mergeCols(["x", "y", "z"], "coords")
         >>> coords = array(tab["coords"])
         >>> print coords[1]
         array([ 5., -1., -5.])

        Or shorter without storing the result in self:
         >>> from numpy import array
         >>> tab = DictTableFromCsvFile("mydata.csv")
         >>> coords = array(tab.mergeCols(["x", "y", "z"]))
         >>> print coords[1]
         array([ 5., -1., -5.])

        @param sourceCols: An iterable stating the names of the columns to be
        merged. Like ["u1", "u2", "u3"] or ["x", "y", "z"]. The last is
        equivalent to sourceCols="xyz".
        @param newColName: name of the newly created merged column. Must not
        exist beforehand. If None (the default) or an empty string then the
        resulting list of lists is *not* stored in self.

        @Returns: The new merged column (a list of lists).

        @raise ValueError: if there is already a column named newColName.
        @raise KeyError: if at least one of the specified sourceCols does not
           exist.

        @note: The original columns still exist independently.
        """
        if newColName and newColName in self:
            raise ValueError(
                "Cannot merge columns because there is already a column named"
                " %s in the DictTableFromCsvFile." % newColName)

        # The following command does something like:
        # self[newColName] = transpose( [ self[col] for col in sourceCols ] )
        # Remember: self is (derived from) a dictionary.
        #
        # Instead of self[colName] we use dict.__getitem__(self, colName) and
        # dict.__setitem__(self, colName) because we want to make sure that we
        # use the lookup function of the underlying dict class, not any
        # possibly overloaded method of the derived class DictTableFromCsvFile.
        #
        # transpose(matrix) can be implemented like that:
        # [list(row) for row in izip(*matrix)]
        # (The list comprehension [list(row) for row in ...] is necessary to
        # get a list of lists, zip(*matrix) would yield a list of tuples. Using
        # izip instead of zip is just for speed and low memory consumption.)
        #
        # With those two tipps you should be able to work out how the following
        # works.
        #
        # And here it is: The beautiful and amazing in any respect but
        # unfortunately completely incomprehensible one liner doing all the
        # trick: (Well it's really been a one liner until newColName became
        # optional and mergeCols() returned something, before ver 1.11.)
        newCol = [list(row) for row in izip(*[
            dict.__getitem__(self, col) for col in sourceCols])]
        if newColName:
            dict.__setitem__(self, newColName, newCol)
        return newCol

    def insertCol(self, key, data, columnIdx=None, exportName=None):
        """Insert the column data with the specified key into self.
        If a column with the same key already exists in self then the data will
        be appended to the existing columns, turning this column into a multi-
        value column (if it hasn't been one before).

        @param data: list of values. Will be truncated or padded with empty
           cells to match the length of the existing columns. Each item in the
           data-list might be a list in itself with possibly different numbers
           if items in it.
        @param columnIdx: Column index for file export. If not specified data
           will be appended to the end.
        @param exportName: Column header string. Defaults to key.
        """

        if exportName is None:
            exportName = key

        if columnIdx is None:
            columnIdx = len(self.headline)

        # correct length of rows
        data = [((not(isinstance(row, list)) and [row]) or row)
                for row in data]
        dataRowLen = max(len(x) for x in data)
        data = [(row + [None]*(dataRowLen-len(row))) for row in data]

        # correct length of data (nb of rows)
        columnLen = self.getColumnLength()
        data = data[:columnLen]
        for i in range(columnLen-len(data)):
            data.append([None]*dataRowLen)

        # append data
        try:
            column = dict.__getitem__(self, key)
        except KeyError:
            # this will be a new column
            if dataRowLen==1:
                # new single-value-column
                dict.__setitem__(self, key, [row[0] for row in data])
                self.headline.insert(columnIdx, exportName)
                self.columnKeys.insert(columnIdx, key)
                self.colValListIndex.insert(columnIdx, None)
            else:
                # new multi-value-column
                dict.__setitem__(self, key, data)
                self.headline[columnIdx:columnIdx] = [exportName]*dataRowLen
                self.columnKeys[columnIdx:columnIdx] = [key]*dataRowLen
                self.colValListIndex[columnIdx:columnIdx] = range(dataRowLen)
        else:
            if isinstance(column[0], list):
                # already is a multi-column
                oldRowLen = len(column[0])  # all rows have the same length
                for row, newItems in izip(column, data):
                    row.extend(newItems)
            else:
                # single column, turn into multi-column
                # update self.colValListIndex: old column is first item
                oldColIdx = self.columnKeys.index(key)
                self.colValListIndex[oldColIdx] = 0
                oldRowLen = 1
                for i, newItems in enumerate(data):
                    column[i] = [column[i]]+newItems
            self.headline[columnIdx:columnIdx] = [exportName]*dataRowLen
            self.columnKeys[columnIdx:columnIdx] = [key]*dataRowLen
            self.colValListIndex[columnIdx:columnIdx] = range(
                oldRowLen, oldRowLen+dataRowLen)

    def removeCols(self, columns):
        """remove specified columns
        @param columns: iterable of keys of columns (strings) to be removed.

        @Note: Use del self[column] to delete individual columns.
        @Note: This method should be faster than multiple del statements if
           removing many columns.
        @Note: Any column key not found in self is silently ignored.
        """
        try:
            removeIds = []
            columns = set(columns)
            for column in columns:
                dict.__delitem__(self, column)
            removeIds.extend((idx for idx,key in enumerate(self.columnKeys)
                              if key in column))
        except TypeError:
            raise ValueError(
                "DictTableFromCsvFile: Column argument to removeCol() must"
                " be a string or iterable of strings.")

        for idx in removeIds:
            del self.columnKeys[idx]
            del self.headline[idx]
            del self.colValListIndex[idx]

    def renameCol(self, oldname, newname):
        """renames a column, i.e. create a new column newname as a shallow
        copy of column oldname, then remove column oldname

        If there is already a column newname in self then it will be silently
        replaced by the former column oldname.
        """
        # move data
        data = dict.__getitem__(self, oldname)
        dict.__setitem__(self, newname, data)
        dict.__delitem__(self, oldname)

        # update self.headline (for output) and self.columnKeys
        # self.colValListIndex stays untouched
        for i in range(len(self.headline)):
            if self.columnKeys[i] == oldname:
                self.columnKeys[i] = newname
                self.headline[i] = newname

    def findInRow(self, value, row, inColumns=None, singleColValIndex=None):
        """Find value in a specific row of self. Find all occurences of value.

        @param value: value to be found

        @param row: Might be either a row index (an integer) referencing the
           row in self.
           Or it might be a row object (i.e. a dictionary) as returned by
           self[rowIndex] or self.getRowColToList() or self.rowIterColToList().

        @param inColumns: list of column names to be searched.

        @param singleColValIndex: Value to return for the list index (second
           item in each tuple of the returned list) if the value has been found
           in a singlular column (not a multi column value).

           The default None might be more suitable to easily identify the type
           of the column (singular column vs. multi column value). On the other
           hand a zero might make it more suitable to use in a
           getRowColToList(listColumns=[all columns]...) context.

        @Returns: A list of (column name, list index)-tuples. The column name
           is suitable to address the column in self. If the corresponding
           value is a list that is, if value is part of a multi column value
           then list index will be the index in that list. Otherwise list
           index will be singleColValIndex (defaulting to None).

           If value is not found an empty list will be returned.

        @Note: For searches in whole table (or over multiple rows) a similar
           but separate method should be written. For efficiency reasons.
        """
        result = list()

        # prepare iterator over all fields/columns in the specified row
        if isinstance(row, int):
            if inColumns is None:
                inColumns = dict.iterkeys(self)
            storedValIter = (
                (colKey, dict.__getitem__(self, colKey)[row])
                for colKey in inColumns)
        elif isinstance(row, dict):
            if inColumns is None:
                storedValIter = row.iteritems()
            else:
                storedValIter = ((colKey, row[colKey])
                                 for colKey in inColumns)

        for colKey, storedVal in storedValIter:
            if isinstance(storedVal, list):
                # multi columns value
                foundpos = 0
                while True:
                    try:
                        foundpos = storedVal.index(value, foundpos)
                    except ValueError:
                        break
                    else:
                        # found, store in result
                        result.append((colKey, foundpos))
                        foundpos += 1
            else:
                # single columns value
                if storedVal == value:
                    result.append((colKey, singleColValIndex))
        return result

    def writeCsv(self, output):
        """write self as csv file.

        @param output: open file object (or something suitable for csv.writer)
            to write the table to. The file should be opened for writing
            binary data, i.e. mode argument to the open function = "wb")
            In case it's a string it is used as filename.

        @note: Columns created by mergeCols are not being exported!
        @note: This is a slow implementation, should be possible to speed that
        up significantly.
        """
        # to speed this up avoid the "for row in self" loop!
        # instead use something similar to the following (needs to account
        # for multi column values and empty columns
        # columns = ([tabBaseSched[col] for col in oldcols]
        #            + [outLevel, outIX, outIY, drawPointsStartDate]
        #            + [tabBaseSched[col] for col in baseSchedColsDrawTons])
        # outtab.writerows(izip(*columns))

        if isinstance(output, basestring):
            output = open(output, "wb")

        outtab = csv.writer(output)
        outtab.writerow(self.headline)

        for row in self:
            line = []
            for key, colValIdx in izip(self.columnKeys, self.colValListIndex):
                if key is None:
                    line.append(None)
                elif colValIdx is None:
                    # single column value
                    line.append(row[key])
                else:
                    # multi column value
                    line.append(row[key][colValIdx])
            outtab.writerow(line)
        return


#------------------------------------------------------------------------------

class RowIterFromCsvFile(object):
    r"""Class for reading a csv file once. Faciliates retrieving data.

    See also L{DictTableFromCsvFile}. In contrast to the latter this class
    does not store the information (other than the header). RowIterFromCsvFile
    is somewhat leaner and meant to iterate through the data only once.

    An Example. Note that columns with identical column names are being merged
    and yield a list of values. First and fourth column labelled "col01" in the
    following example:
     >>> mycsv = CStringIO('''\ ...
     ... "col01","col02","col03","col01","with space"
     ... "lala",1,2.5,"blub","gaga"
     ... ''')
     >>> tab = RowIterFromCsvFile(mycsv)
     >>> for row in tab.getContainerIter():
     ...     print row.col01, row.col02, type(row.col02), row.with_space
     ['lala', 'blub'] 1 <type 'str'>

    Choose some columns and apply a conversion function:
     >>> mycsv = ...
     >>> tab = RowIterFromCsvFile(
     ...     mycsv, columns=["col01", ("col02", float), "with_space"])
     >>> for row in tab.getContainerIter():
     ...     print row.col01, row.col02, type(row.col02)
     ['lala', 'blub'] 1.0 <type 'float'>

    You may also get a dictionary instead of a Container object and specify
    column indices (zero based) instead of column names:
     >>> mycsv = ...
     >>> tab = RowIterFromCsvFile(
     ...     mycsv, columns=["col01", (1, float), "with_space"])
     >>> for row in tab.getDictIter():
     ...     print row["col01"], row["col02"]
     ['lala', 'blub'] 1.0

    The result -i.e. row- might also be a list:
     >>> mycsv = ...
     >>> tab = RowIterFromCsvFile(
     ...     mycsv, columns=["col01", (1, float), "with_space"])
     >>> for row in tab.getListIter():
     ...     print row
     [['lala', 'blub'], 1.0, "gaga"]

    Columns might be identified through a slice object as well and be assigned
    arbitrary names independent of the column name in the csv file:
    Note: list[slice(None, 4)] == list[:4] and list[slice(3,None)] == list[3:]
     >>> mycsv = ...
     >>> tab = RowIterFromCsvFile(
     ...     mycsv, columns=["col01", (slice(1,-1), float, "mycol")])
     >>> for row in tab.getContainerIter():
     ...     print row.col01, row.mycol
     ['lala', 'blub'] [1.0, 2.5]

    There might be no header in the csv file (or multi line headers, see the
    headLines argument):
     >>> mycsv = ...
     >>> tab = RowIterFromCsvFile(mycsv, columns=[[0,3], 1],
     ...                          headLines=0)
     >>> for row in tab.getListIter():
     ...     print row[0], row[1]
     ['col01', 'col01'], 'col02'
     ['lala', 'blub'] 1

    @ivar headline: headline (as is, no conversions)
    @ivar columnKeys: list of column keys for each column (converted headline)
    @ivar keyColNb: {key: column number or list of column numbers}
    @ivar recognizedColsList: set of indices of the recognized columns
    @ivar attrNameToColumnId: {attrName: list of column }
    @ivar columns: converted columns description. Like the corresponding
       argument to the constructor. All items are (key, converter, attrName)
       -tuples with sensible defaults for what has not been given. All column
       keys as integer indices.
    @ivar fileName: file name (mostly for diagnostic output)

    @Note: You can retrieve data from the same column multiple times, for
    example in a multi column field and again as single column field, i.e.
     >>> ... columns=[(["x", "y", "z"], float, "coords"),
     ...              ("x",             int,   "ix"    ),  ...

    @Note: replaceHeaderSpaceBy and replaceHeaderFeedBy are applied before
    the columns argument is evaluated. That means if you specify column names
    in the columns argument they must much the already converted form of the
    column key. E.g. suppose in the csv file you have "M3C2 distance" in the
    header line you'd address the corresponding column as "M3C2_distance" in
    the columns argument. And then after row = tab.getDictIter().next() you'd
    find a value with row["M3C2_distance"].
    Or after row = tab.getContainerIter().next() you'd do: row.M3C2_distance.
    """

    def __init__(self, csvFile, columns=None, logfile=None,
                 replaceHeaderSpaceBy="_", warnIgnoreCols=False,
                 headLines="auto", replaceHeaderFeedBy="_",
                 **kwargs):
        r"""
        Constructor

        Besides those listed below you may specify additional keyword
        arguments. They will be passed to the csv.reader() function. See
        http://docs.python.org/library/csv.html#csv-fmt-params for possible
        format parameters or use the parameter dialect to specify a certain
        flavour of csv file.

        Example, columns separated by tabs:
         >>> tab = RowIterFromCsvFile("mytable.csv", delimiter="\t")

        @param csvFile: either a filename or an already open file that
          contains the table data

        @param columns: A list to identify the columns to be read from the csv
          file one by one. Each item in this list represents one attribute of
          the Container object returned by self.next().

          An item of the columns list might be a single key value. In the
          simplest case this is an integer identifying a column by index or a
          string identifying a column by the key extracted from the head line.
           >>> ... columns = [ 0, "col02", "col03"] ...

          Instead, an item of the columns list might as well be a
          (key, conversion, attrName)-tuple (really a tuple not a list, to
          distinguish from key-lists, see below). In this case conversion would
          be a function applied to the value of each individual column (like
          int, float) or None to indicate no conversion. And attrName is the
          name of the corresponding attribute of the Container object returned
          by self.getContainerIter() or the key in the dictionary returned by
          self.getDictIter().
           >>> ... columns = [
           >>>         (0,       None,  "NameColumn"),
           >>>         ("col02", int,   "SomeNumber"),
           >>>         ("col03", float, "AnotherNumber") ] ...

          As said above a key might be an integer identifying a column by
          index. Or a string identifying a column by the key extracted from the
          head line(s). In those cases the resulting fields returned by
          self.getContainerIter() and the like would be single values (possibly
          converted to numbers). Or if there are multiple columns with the
          specified (string valued) key then the resulting field would be a
          list of those values. See examples in the description of the class.

          Now, additionally to what's been said before the key value specified
          for a particular column might also be a list of integers and/or
          strings identifying multiple columns in the csv file. Or it might be
          a slice object identifying columns by a range of indices. In those two
          cases (key being a list or a slice) the corresponding result will
          always be a list of values from the corresponding columns.
           >>> ... columns = [ ["x", "y", "z"], "col2" ] ...

          or:
           >>> ... columns = [
           >>>         (["x", "y", "z"], float, "coords"),
           >>>         (["Ux","Uy","Uz"], float, "displ")  ] ...

          mixed forms are possible as well:
           >>> ... columns = [
           >>>         ([0,1,2], float, "coords"),
           >>>         "col02", "col03  ] ...

          This is an optional argument. If not specified all columns will be
          read.

        @param replaceHeaderSpaceBy: If not None, all spaces in each column
          header will be replaced by this string.

        @param logfile: DEPRECATED argument being ignored. use the bae.log_01
          functionality for diagnostic output

        @param warnIgnoreCols: If False data in columns without corresponding
          header field are silently ignored. Otherwise you get some kind of
          warning: If True or "message" there will be a warning message.
          If "exception" a IgnoredColumns Exception will be raised

        @param headLines: If None or zero or False, no head line expected. All
          lines will be treated as data lines. If a positive integer N then the
          first N lines will be treated as head lines. If headLines=="auto" then
          the first non empty line will be used, following empty lines will be
          skipped.

        @param replaceHeaderFeedBy: Seperator between the contents of a
          header field on multiple rows. This conversion will be performed to
          get the keys for the data in the specified column.
        """

        # initialize the csv file reader
        #  ... and find a file name for diagnostic output
        if isinstance(csvFile, basestring):
            self.fileName = csvFile
            csvFile = open(csvFile, "rb")
        else:
            try:
                self.fileName = csvFile.name
            except AttributeError:
                self.fileName = "input of type %s" % type(csvFile).__name__

        csvTable = csv.reader(csvFile, **kwargs)

        self.ignoreEmptyRows = False
        if headLines=="auto":
            # find the first line that is not empty
            headline = None
            for headline in csvTable:
                if len(headline) and any(map(len, headline)):
                    break
            # check if file is empty
            if headline is None or len(headline)==0:
                # file empty, no headline
                raise ValueError(
                    "RowIterFromCsvFile.__init__(): The specified csv file %s"
                    " is empty." % self.fileName)
            self.ignoreEmptyRows = True
        elif not(headLines):
            headline=[]
        elif headLines>0:
            try:
                lines = [csvTable.next() for cnt in range(headLines)]
            except StopIteration:
                raise ValueError(
                    "RowIterFromCsvFile.__init__(): The specified csv file %s"
                    " does not contain the specified number (%d) of header"
                    " lines." % (self.fileName, headLines))
            # in python 2.4 you may not specify additional arguments after
            # a *args construct. The following is valid in python 2.6:
            # headline = ["\n".join(x)
            #             for x in izip_longest(*lines, fillvalue="")]
            # ... and here is the replacement for earlier versions
            headline = ["\n".join(x)
                        for x in izip_longest(*lines, **{"fillvalue":""})]
        else:
            raise ValueError("Illegal argument headLines to RowIterFromCsvFile"
                             " constructor: %s" % headLines)

        #--- evaluate headline ...

        # self.columnKeys ... list of column keys for each column (converted
        #                     headline)
        self.columnKeys = list(headline)
        if replaceHeaderSpaceBy is not None:
            self.columnKeys = [item.replace(" ", replaceHeaderSpaceBy)
                               for item in self.columnKeys]
        if replaceHeaderFeedBy is not None:
            self.columnKeys = [item.replace("\n", replaceHeaderFeedBy)
                               for item in self.columnKeys]

        # self.keyColNb ... {key: column number or list of column numbers}
        # self.recognizedColsList ... set of indices of the recognized columns
        self.keyColNb = dict()
        self.recognizedColsList = set()
        for cnt, item in enumerate(self.columnKeys):
            if item:
                try:
                    # this will only work for the third and subsequent
                    # identical column headers
                    self.keyColNb[item].append(cnt)
                except KeyError:
                    # first item of that name: just store the column number
                    self.keyColNb[item] = cnt
                except AttributeError:
                    # second item of that name: store a list of column numbers
                    lastColOfThisKey = self.keyColNb[item]
                    self.keyColNb[item] = [lastColOfThisKey, cnt]

                self.recognizedColsList.add(cnt)

        #-- create columns from header line if not specified
        if columns is None:
            if not self.columnKeys:
                raise ValueError(
                    "RowIterFromCsvFile.__init__(): if there is no header line"
                    " columns must be specified!")
                # this is a limitation: Otherwise we would have to read at
                # least the first row to learn how many columns we've got.

            # create columns
            columns = list()
            processedKeys = set()
            for col in self.columnKeys:
                if col in processedKeys:
                    continue
                processedKeys.add(col)
                columns.append(col)

        #-- reformat columns for easier access
        newcolumns = list()
        doNothing = lambda x:x     # a converter function that does nothing

        # also store minimum number of columns expected in each row
        self.minColNb = 0

        for col in columns:

            # make sure it's a tuple with three items
            if not isinstance(col, tuple):
                col = (col, None, None)
            else:
                col += (None,)*(3-len(col))

            # key...
            key = col[0]
            if isinstance(key, basestring):
                try:
                    key = self.keyColNb[key]
                except KeyError:
                    raise KeyError("No column key %s in csv file %s."
                                   % (key, self.fileName))
            elif isinstance(key, list):
                # convert all string keys to the corresponding indices
                newkey = list()
                for kkk in key:
                    if isinstance(kkk, basestring):
                        try:
                            kkk = self.keyColNb[kkk]
                        except KeyError:
                            raise KeyError("No column key %s in csv file"
                                           " %s." % (kkk, self.fileName))
                    if isinstance(kkk, int):
                        newkey.append(kkk)
                    else:
                        newkey.extend(kkk)
                key = newkey
            if not isinstance(key, (int, list, slice)):
                raise ValueError("RowIterFromCsvFile.__init__(): Could not"
                                 " handle column key <%s>." % key)

            # if no converter then use doNothing
            if col[1]:
                converter = col[1]
            else:
                converter = doNothing

            # attrName
            if isinstance(col[2], basestring) and col[2]:
                # attrName given in tuple
                attrName = col[2]
            elif isinstance(col[0], basestring):
                # attrName from type-string column key in argument
                attrName = col[0]
            elif self.columnKeys:
                # attrName from column header in csv file
                if isinstance(key, list):
                    attrName = self.columnKeys[key[0]]
                else:
                    attrName = self.columnKeys[key]
            elif isinstance(col[0], int):
                # generate column name from column index
                attrName = "col%d" % col[0]
            else:
                raise ValueError(
                    "RowIterFromCsvFile.__init__(): If attrName (third"
                    " item in the key tuple of the columns argument) is"
                    " not specified, the key (first item in the tuple)"
                    " must be a string or an int. Found key: %s,"
                    " attrName: %s" % (col[0], col[2]))

            newcolumns.append((key, converter, attrName))

            # update self.minColNb (minimum number of columns expected)
            # this number is later used to pad columns from the table to this
            # minimum length
            if isinstance(key, list):
                self.minColNb = max([self.minColNb]+[i+1 for i in key])
            elif isinstance(key, slice):
                if key.stop>0:
                    self.minColNb = max(self.minColNb, key.stop)
                elif self.columnKeys:
                    actualStopOfSlice = key.indices(len(self.columnKeys))[1]
                    self.minColNb = max(self.minColNb, actualStopOfSlice)
                else:
                    raise ValueError(
                        "RowIterFromCsvFile.__init__(): A slice stop relative"
                        " to the end can only be used if there is a header.")
                    # this is a limitation: Otherwise we would have to read at
                    # least the first row to learn how many columns we've got.
            else:
                self.minColNb = max(self.minColNb, key+1)

        # rename duplicate attrNames
        # (may occur if attrName is taken from the column headers in the csv)
        attrNameCnt = defaultdict(int)
        for key, converter, attrName in newcolumns:
            attrNameCnt[attrName] += 1
        # remove entries for singular columns
        newAttrNameCnt = dict()
        for attrName, cnt in attrNameCnt.iteritems():
            if cnt>1:
                newAttrNameCnt[attrName] = cnt
        attrNameCnt = newAttrNameCnt
        # rename multi column attrNames
        # generate attrNameToColumnId: {attrName: list of column keys}
        self.attrNameToColumnId = dict()
        columns = newcolumns
        newcolumns = list()
        for key, converter, attrName in columns[::-1]:
            try:
                cnt = attrNameCnt[attrName]
            except KeyError:
                pass
            else:
                attrNameCnt[attrName] -= 1
                attrName = "%s_%d" % (attrName, cnt)
            newcolumns.append((key, converter, attrName))
            self.attrNameToColumnId[attrName] = key
        newcolumns.reverse()

        # permanently store ...
        self.headline = headline   # headline (as is, no conversions)
        self.csvTable = csvTable
        self.columns = newcolumns  # converted columns description
        self.warnIgnoreCols = warnIgnoreCols

    def getContainerIter(self):
        """Return an iterator that yields Container objects.

        Example: Suppose the file mycsv contains column header items "col01",
        "col02" and "extracol"

         >>> tab = RowIterFromCsvFile(mycsv)
         >>> for row in tab.getContainerIter():
         >>>     print row.col01, row.col02, row.extracol
        """
        return self._getGenerator(newResult=Container, assignResult=setattr)

    def getDictIter(self):
        """Return an iterator that yields dictionaries.

        Example: Suppose the file mycsv contains column header items "col01",
        "col02" and "extracol"


         >>> tab = RowIterFromCsvFile(mycsv)
         >>> for row in tab.getDictIter():
         >>>     print row["col01"], row["col02"], row["extracol"]
        """
        return self._getGenerator(newResult=dict,
                                  assignResult=dict.__setitem__)

    def getListIter(self):
        """Return an iterator that yields lists.

        Example: Suppose the file mycsv contains 10 columns and only the first,
        fourth and 11th are of interest. The row objects returned by the
        iterator self.getListIter() are lists with three items:

         >>> tab = RowIterFromCsvFile(mycsv, columns=[0, 3, 9])
         >>> for row in tab.getListIter():
         >>>     print row[0], row[1], row[2]
        """
        def assignresult(res, attrName, val):
            res.append(val)
        return self._getGenerator(newResult=list,
                                  assignResult=assignresult)

    def _getGenerator(self, newResult, assignResult):
        for row in self.csvTable:

            # possibly ignore empty lines at the beginning
            if self.ignoreEmptyRows:
                if len(row)==0 or not(any(map(len, row))):
                    continue
                else:
                    self.ignoreEmptyRows = False

            # pad row with ""
            row += [""]*(self.minColNb-len(row))

            # check wether there are values in not recognized columns
            if self.warnIgnoreCols and self.recognizedColsList:
                nonEmptyCols = set(
                    colNb for colNb, val in enumerate(row) if val!="")
                ignoredCols = nonEmptyCols.difference(self.recognizedColsList)
                if len(ignoredCols):
                    msgText = (
                        "The following %d columns have data but no header:"
                        "\n%s" % (len(ignoredCols), sorted(ignoredCols)))
                    if self.warnIgnoreCols=="exception":
                        raise IgnoredColumns(msgText)
                    else:
                        msg("WARNING: %s" % msgText)

            # retrieve values
            res = newResult()
            for key, converter, attrName in self.columns:
                if isinstance(key, list):
                    val = [converter(row[i]) for i in key]
                elif isinstance(key, slice):
                    val = [converter(x) for x in row[key]]
                else:
                    val = converter(row[key])
                assignResult(res, attrName, val)

            yield res

    def getColumnKeysForAttrName(self, attrName):
        """returns the column keys (as in self.columnKeys) for the specified
        attrName.

        @param attrName: One of the attribute names of the Container
        objects returned by self.getContainerIter() or one of the keys in the
        dict returned by self.getDictIter().

        @returns: A string: the column key. Or a list of such strings in case
        of a multi column field.
        """
        key = self.attrNameToColumnId[attrName]
        if isinstance(key, list):
            columnKeys = [self.columnKeys[i] for i in key]
        else:  # key might be of type int or slice
            columnKeys = self.columnKeys[key]
        return columnKeys


#------------------------------------------------------------------------------

class RowIterFromAbqReport(object):
    r"""Class for reading an Abaqus report file. Faciliates retrieving data.
    Shall be used as an iterator yielding one row after the other.

    Example: read all data from last report (abaqus.rpt in current work dir)
     >>> from bae.misc_01 import RowIterFromAbqReport
     >>> report = RowIterFromAbqReport()
     >>> print "header:", report.columnKeys
     >>> tab = list(report)
     >>> print "second row:", tab[1]
     >>> print "second row, third column:", tab[1][2]

    Example: iterate over data rows:
     >>> report = RowIterFromAbqReport()
     >>> for row in report:
     >>>    print "X:", row[0], "Y:", row[1]

    @ivar columnKeys: list of column headers
    @note: Can read the file only once.
    @note: IMPORTANT! The abaqus report file format is exceptionally difficult
        to parse automatically. This class worked for some particular files.
        You should always check at least once for each use case whether the
        results are correct.
    """
    def __init__(self, abqReport="abaqus.rpt"):
        r"""
        @param abqReport: either a filename or an already open file that
          contains the table data
        """

        # open input file and find a file name for diagnostic output
        if isinstance(abqReport, basestring):
            self.fileName = abqReport
            self.inpFile = open(abqReport, "r")
        else:
            try:
                self.fileName = abqReport.name
            except AttributeError:
                self.fileName = "input of type %s" % type(abqReport).__name__
            self.inpFile = abqReport

        # find the first line that is not empty
        headline = None
        for headline in self.inpFile:
            if len(headline.strip()):
                break
        # check if file is empty
        if len(headline.strip())==0:
            # file empty, no headline
            raise ValueError(
                "RowIterFromAbqReport.__init__(): The specified file %s"
                " is empty." % self.fileName)

        if headline.strip()[0]=="X":
            # single line head line
            self.columnKeys = headline.strip().split()

            # skip following empty line
            headline = self.inpFile.next()
            assert len(headline.strip())==0, "Unrecognized format variant."

        else:
            # multi line head line. In this case the first of the head lines
            # has the full column width.
            # the first column "X" does not have anything on the first
            # head line. "X" stands on the second or last, I don't know.

            # columns: split at double space following a non space
            columnStarts = [m.start()
                            for m in re.finditer(r"(?<=[^ ])  ", headline)]
            # first col not detected by this,
            # prepend 0 and pos of first non space
            columnStarts[0:0] = (0, re.search(r"\S", headline).start())

            # column slices
            columnSlices = [slice(start, next) for start, next in izip(
                columnStarts, columnStarts[1:])]
            columnSlices.append(slice(columnStarts[-1], None))

            # create self.columnKeys
            self.columnKeys = [headline[sli].strip() for sli in columnSlices]

            # add next rows contents to columnkeys until a blank line occurs
            for headline in self.inpFile:
                if len(headline.strip())==0:
                    break
                else:
                    for i, sli in enumerate(columnSlices):
                        self.columnKeys[i] += headline[sli].strip()

    def __iter__(self):
        for line in self.inpFile:
            if len(line.strip())==0:
                break
            row = map(float, line.split())
            if len(row) != len(self.columnKeys):
                raise ValueError(
                    "RowIterFromAbqReport: number of column in %s is"
                    " inconsistent with header. Current line:\n%s\n%s\n"
                    "...header:\n%s"
                    % (self.fileName, line, row, self.columnKeys))
            yield row


#------------------------------------------------------------------------------

def convertAbqReport(inpName="abaqus.rpt", outName=None):
    """converts an Abaqus XY report to a csv file.

    @param inpName: name of the abaqus report file to be read
    @param outName: name of the resulting csv file. If None (the default),
      inpName is converted: ".rpt" is replaced by ".csv"

    @Note: Silently overwrites the file outName.
    """

    if outName is None:
        outName = inpName.rsplit(".rpt", 1)[0]+".csv"

    tab = RowIterFromAbqReport(inpName)
    output = csv.writer(open(outName, "wb"))
    output.writerow(tab.columnKeys)
    output.writerows(tab)

#} end of Spreadsheets: Data from csv files


###############################################################################
###############################################################################
###############################################################################
#{ customize to different run-time environments (e.g. with or without numpy)
#

def selectattrversion(flag, winp, wonp, doctext):
    r"""
    Choose adequate method (or other attribute) in class definition depending
    on a flag indicating whether numpy can be used or not (or any other
    condition).

    Usage:
     >>> from bae.misc_01 import selectattrversion
     >>> try:
     ...     import numpy as np
     ... except ImportError:
     ...     np = None
     >>>
     >>> class Myclass(object):
     ...     def _np_func(self, x, y):
     ...         dosomething(x,y)
     ...     def _nonp_func(self, x, y):
     ...         dosomethingelse(x,y)
     ...     func = selectattrversion(np, _np_func, _nonp_func, '''
     ...         The method self.func either does something or something else.
     ...         @param x: that's for singing
     ...         @param y: and that's for dancing
     ...         ''')
     >>> m = Myclass()
     >>> result = m.func(forsinging, fordancing)
     >>>
     >>> def _np_dot(x,y):
     ...     return np.dot(x,y)
     >>> def _nonp_dot(x,y):
     ...     return sum(xx*yy for xx,yy in zip(x,y))
     >>> dot = selectattrversion(np, _np_dot, _nonp_dot, '''
     ...     Calculate the dot (scalar) product of the two vectors x and y
     ...     ''')
     >>> a = [1,2,3]; b = [5,6,7]
     >>> print "a dot b = ", dot(a,b)

    The class Myclass has a func() method either being _np_func() if numpy has
    been imported, or otherwise the func() method will be the _nonp_func()
    function.

    The function dot() will either be _np_dot() which uses numpy or _nonp_dot()
    which does not use numpy.

    You will usually not call the _np_dot() and _nonp_dot() functions directly.

    @Note: The arguments to the selectattrversion() function must be defined
    before this function is called, therefore the order in which e.g.
    _np_dot(), _nonp_dot() and dot() are defined does matter: dot() must be
    last.

    @Note: There is a widespread convention that "secret", "hidden" or
    "private" objects are named with a single leading underscore. E.g. Epydoc
    does not list them by default. Therefore it is suggested to start the names
    of the special variants like _np_dot() with an underscore because then they
    are somehow hidden. Only the general function dot() should be visible.
    """
    if flag:
        res = winp
    else:
        res = wonp
    res.__doc__ = doctext
    return res
#} end of customize to different run-time environments (e.g. with or w/o numpy)


###############################################################################
###############################################################################
###############################################################################
#{ point search: find closest points, KDTree, inverse distance interpolation
#

# ------------------------------------------------------------------
# Svens getShortestDistKey
def _getShortestDistKeywNumpy(pt,dct):
    # numpy is imported already
    if isinstance(dct, dict):
        # keys and values of given dict:
        keys,vals = zip(*(dct.iteritems()))
    elif isinstance(dct, (list, tuple)):
        keys=None
        vals=dct
    else:
        raise ValueError("getShortestDistKey only accepts a list, tuple or"
                         " dict as second argument. (got %s)" % type(dct))

    # vectors from each centroid to the given point:
    # Say:
    # As input we have pt of dimension m x 3 and dct of dimension n x 3
    # we compute the difference between the values of dct and each of m
    # entries of pt
    # The result is an array of shape m x n x 3:
    # Or:
    # As input we have just one point, pt has dimension 3 and dct is of
    # dimension n x 3. We compute the difference between the values of dct
    # and pt. The result is an array of shape n x 3:
    pt = np.asanyarray(pt)
    d = np.asanyarray(vals)[np.newaxis,:,:] - pt[...,np.newaxis,:]
    # element wise squared:
    np.square(d,d)
    # sum of components = square of lenghts:
    # the components are addressed with the last index:
    lq=d.sum(axis=-1)
    # find index, at which norm is minimal:
    index=np.argmin(lq,axis=1)
    # return key (element number s)
    if keys is None:
        return index
    elif len(pt.shape)>1:
        # return an array of indexes
        return np.array(keys)[index]
    else:
        # if pt is just one point, i.e. pt.shape = (3,)
        # return just an index
        return keys[index[0]]


def _getShortestDistKeywoNumpy(pt,dct):
    # numpy could not be imported
    # vecmath used instead
    # dict of distances

    # not implemented for list of points yet
    if not isinstance(pt[0],(int,float)):
        raise Exception(
            "getShortestDistKey without numpy only accepts "
            "exactly one point as first argument (got %s)" % type(pt))

    if isinstance(dct, dict):
        dLst=[[dist(pt,otherpt),k] for k,otherpt in dct.iteritems()]
    elif isinstance(dct, (list, tuple)):
        dLst=[[dist(pt,otherpt),k] for k,otherpt in enumerate(dct)]
    else:
        raise ValueError("getShortestDistKey only accepts a list, tuple or"
                         " dict as second argument. (got %s)"
                         % type(dct))

    dLst.sort()
    # return key (element number):
    return dLst[0][1]

getShortestDistKey = selectattrversion(
    np, _getShortestDistKeywNumpy, _getShortestDistKeywoNumpy,
    r"""
compute the (key of the) shortest distance to a given point

Usage:
 >>> elCentDct={}
 >>> point=[0.,0.,0.]
 >>> elCentDct[1]=[1.,0., 0.]
 >>> elCentDct[2]=[3.,4., 0]
 >>> elCentDct[3]=[0.1,0.1,0.1]
 >>> nearest = getShortestDistKey(point,elCentDct)
 >>> print "Point %d is closest to %s." % (nearest, point)

Or:
 >>> pointList=[[1.,0., 0.], [3.,4., 0], [0.1,0.1,0.1]]
 >>> nearest = getShortestDistKey(point,pointList)
 >>> print "The %d. point is closest to %s." % (nearest+1, point)

@param pt: A point given as list or so. Or a list/array of points (this works
only for the numpy version so far).
@param dct: A dict or a list. A dict would consist of keys (e.g. element
number) and corresponding points (e.g. centroids). A list would simply
contain the points, the positions in the list would be treated as corresponding
keys.

@Return: The key or index of the closest point if the pt argument is just one
point. Or an array of keys or indexes in case the pt argument is a list of
points.

@Note: See also the KDTree class. It's practical for many (like 1E6) points.
""")


#------------------------------------------------------------------------------

def mapSamePoints(points1, points2, tolerance=1E-2):
    """The arguments are two lists of point coordinates. All points in points1
    are to be found in points2 as well within a certain tolerance.

    For each pt in points1 ("source") find the index in points2 ("target").

    Example:
     >>> points1 = [ [1e-4, -1e-4, 0.0], [1, 1e-4, -1e-4] ]
     >>> points2 = [ [1,0,0], [0,0,1], [0,0,0], [1,0,1] ]
     >>> result = mapSamePoints(points1, points2)
     >>> print result
     [ 2, 0 ]
     >>> len(result) == len(points1)
     True
     >>> zip(points1, [points2[i] for i in result])
     [ [ [1e-4,-1e-4,0.0], [0,0,0] ],
     [ [1,1e-4,-1e-4]], [1,0,0] ] ]

    @param points1: List of points that are to be assigned an index to.

    @param points2: List of points with a certain order. If two points in this
       list are the same (see parameter tolerance for the meaning of "same")
       then one of them hides the other. (The point occuring later in points2
       hides the other.)

    @param tolerance: Two points are always identified as same if their
       distance is at most tolerance. Two points are never being identified as
       same if their distance is larger than appr. 3.5*tolerance (exactly
       2*sqrt(3)*tolerance). If the distance is between those bounds the
       outcome depends on the actual positions of the points with respect to
       the grid used internally, i.e. it is random.

    @returns: A list of indexes for each point in points1. The index points to
       the corresponding point in points2. I.e. use this index to address
       points in points2.

       The length of the resulting list is the same as points1.

       If a certain point from point1 did not match any of the points in
       point2 the corresponding item in the result list is None.

       Many points from points1 may be found corresponding to the same point
       in points2.
    """
    # Do the calculations relative to the smallest coordinate
    # . no neg values! otherwise problem: int(-0.9) == int(0.9) == 0
    # . in case all points are close together compared to their coordinates
    #   this might help avoiding round-off errors.
    # . For small tolerances it might prevent the necessity to use arbitrary
    #   sized integer, which is done automatically by python but suffers from a
    #   performance penalty
    origin = [min(x[0] for x in points1+points2),
              min(x[1] for x in points1+points2),
              min(x[2] for x in points1+points2)]
    origin = [int(x/tolerance-0.5)*tolerance for x in origin]

    # note on variable names: ac is approximate coordinate
    # { approx coord: index in points } or all points
    # reason for +0.5: assuming that it's more likely that point coordinates
    # are a multiple of tolerance (i.e. coordinates are already rounded values)
    # then this +0.5 avoids direct-cell-misses and iterating through neighbour
    # cells.
    acId2Dict = dict(
        ((int((pt[0]-origin[0])/tolerance+0.5),
          int((pt[1]-origin[1])/tolerance+0.5),
          int((pt[2]-origin[2])/tolerance+0.5)), idx)
        for idx, pt in enumerate(points2))

    # list of approx coords of all points in points1
    points1Acs = [
        (int((pt[0]-origin[0])/tolerance+0.5),
         int((pt[1]-origin[1])/tolerance+0.5),
         int((pt[2]-origin[2])/tolerance+0.5))
        for pt in points1]

    # find all points directly in the cell (same approx coordinate)
    result = [acId2Dict.get(ac) for ac in points1Acs]

    # for all points not yet associated: list of indexes in points1
    notFoundId1s = [id1 for id1, ptId2 in enumerate(result)
                    if ptId2 is None]

    if len(notFoundId1s):
        searchOffsets = [
            # face neighbours
            ( 1, 0, 0),
            (-1, 0, 0),
            ( 0, 1, 0),
            ( 0,-1, 0),
            ( 0, 0, 1),
            ( 0, 0,-1),
            # edge neighbours
            ( 1, 1, 0),
            ( 1,-1, 0),
            (-1, 1, 0),
            (-1,-1, 0),
            ( 0, 1, 1),
            ( 0, 1,-1),
            ( 0,-1, 1),
            ( 0,-1,-1),
            ( 1, 0, 1),
            (-1, 0, 1),
            ( 1, 0,-1),
            (-1, 0,-1),
            # corner neighbours
            ( 1, 1, 1),
            ( 1, 1,-1),
            ( 1,-1, 1),
            ( 1,-1,-1),
            (-1, 1, 1),
            (-1, 1,-1),
            (-1,-1, 1),
            (-1,-1,-1),
            ]

    for id1 in notFoundId1s:
        ac = points1Acs[id1]
        for offs in searchOffsets:
            acoffs = tuple(x1+x2 for x1, x2 in zip(ac, offs))
            try:
                id2 = acId2Dict[acoffs]
            except KeyError:
                continue
            else:
                result[id1] = id2
                break

    return result


#------------------------------------------------------------------------------

def findDuplicatePoints(points, tolerance=1E-2):
    """Find identical points. Identical means: approximately closer than
    tolerance to each other.

    Example:
     >>> points = [ [1,0,0], [5,0,0], [4.9999,0.0,0.0], [1e-4, -1e-4, 0.0],
     >>>            [1, 1e-4, -1e-4] ]
     >>> result = findDuplicatePoints(points)
     >>> print result
     [[0,4],[1,2]]

    @param points: List of points.

    @param tolerance: Two points are always identified as same if their
       distance is at most tolerance. Two points are never being identified as
       same if their distance is larger than appr. 3.5*tolerance (exactly
       2*sqrt(3)*tolerance). If the distance is between those bounds the
       outcome depends on the actual positions of the points with respect to
       the grid used internally, i.e. it is random.

    @returns: A list of list of indexes of identical points.

    @Note: This seems not to work with Rhino.
       Test-scenario: a large ground-support-mesh (50000 lines) shall be
       exported with exportAsAbaqusModel. exportAsAbaqusModel uses
       findDuplicatePoints for identifying connections of the lines.
       Using this function Rhino crashes.

       There is a replacement version in bert.exim.exportAsAbaqusModel that
       does not contain the try dictionary access phrase:
       "try: other = acIdDict[ac]". This version does not crash.
    """
    # Do the calculations relative to the smallest coordinate
    # . no neg values! otherwise problem: int(-0.9) == int(0.9) == 0
    # . in case all points are close together compared to their coordinates
    #   this might help avoiding round-off errors.
    # . For small tolerances it might prevent the necessity to use arbitrary
    #   sized integer, which is done automatically by python but suffers from a
    #   performance penalty
    origin = [min(x[0] for x in points),
              min(x[1] for x in points),
              min(x[2] for x in points)]
    origin = [int(x/tolerance-0.5)*tolerance for x in origin]

    # linkedList always contains the index of the successor in the chain
    # singulars point at themselves, i.e. linkedList[i]==i
    # to traverse one chain loop until you encounter the initial item again
    linkedList = range(len(points))

    # note on variable names: ac is approximate coordinate
    # { approx coord: index in points/linkedList } for all points
    # reason for +0.5: assuming that it's more likely that point coordinates
    # are a multiple of tolerance (i.e. coordinates are already rounded values)
    # then this +0.5 avoids direct-cell-misses and iterating through neighbour
    # cells.
    acIdDict = dict()
    for idx, pt in enumerate(points):
        ac = (int((pt[0]-origin[0])/tolerance+0.5),
              int((pt[1]-origin[1])/tolerance+0.5),
              int((pt[2]-origin[2])/tolerance+0.5))

        try:
            other = acIdDict[ac]
        except KeyError:
            acIdDict[ac] = idx
        else:
            # here still: linkedList[idx] == idx
            linkedList[idx] = linkedList[other]
            linkedList[other] = idx
    msg("Created acIdDict with %d entries from %d points."
        % (len(acIdDict), len(points)), debugLevel=2)

    # # for debugging purposes only:
    # def getCurrentChain(linkedList):
    #     idxNext = dict(enumerate(linkedList))
    #     chainList = list()
    #     while idxNext:
    #         idx, next = idxNext.popitem()
    #         if next != idx:
    #             chain = []
    #             try:
    #                 while 1:
    #                     chain.append(next)
    #                     next = idxNext.pop(next)
    #             except KeyError:
    #                 pass
    #             chainList.append(chain)
    #     return chainList

    # check neighbouring cells, join corresponding chains
    for ac in sorted(acIdDict):
        idx = acIdDict[ac]
        predecessor = None
        for d in [(0,0,1),
                  (0,1,-1),(0,1,0),(0,1,1),
                  (1,0,-1),(1,0,0),(1,0,1),
                  (1,1,-1),(1,1,0),(1,1,1)]:
            dac = (ac[0]+d[0], ac[1]+d[1], ac[2]+d[2])
            try:
                other = acIdDict[dac]
            except KeyError:
                continue
            else:
                # merge in
                if predecessor is None:
                    # find predecessor -> idx -> linkedList[idx]
                    predecessor = idx
                    chain = set((idx,))
                    next = linkedList[idx]
                    while next!=idx:
                        chain.add(next)
                        predecessor = next
                        next = linkedList[predecessor]

                if other not in chain:
                    # merge chains
                    next = linkedList[other]
                    linkedList[other] = idx
                    linkedList[predecessor] = next
                    predecessor = other
                    while next!=idx:
                        chain.add(next)
                        next = linkedList[next]

    # extract chains from linkedList
    idxNext = dict(enumerate(linkedList))
    chainList = list()
    while idxNext:
        idx, next = idxNext.popitem()
        if next != idx:
            chain = []
            try:
                while 1:
                    chain.append(next)
                    next = idxNext.pop(next)
            except KeyError:
                pass
            chainList.append(chain)

    return chainList


#------------------------------------------------------------------------------

class KDTree(object):
    """KDTree for nearest neighbour search in 3D

    Usage:
     >>> from bae.misc_01 import KDTree
     >>> # create a list of 3D points
     >>> from random import uniform
     >>> from math import sqrt
     >>> box = [[0.0, 0.0, 0.0], [20.0, 10.0, 5.0]]; N=1000
     >>> points = [[uniform(box[0][i], box[1][i]) for i in xrange(3)]
     ...           for j in xrange(N)]
     >>> # now find closest points to centre, min and max corner
     >>> tree = KDTree(points)
     >>> targetPoints = [
     ...     [0.5*(box[0][i]+box[1][i]) for i in xrange(3)],
     ...     map(min, zip(*box)),
     ...     map(max, zip(*box))]
     >>> for targetpoint in targetPoints:
     >>>     closest = tree.knnSearch(targetpoint, K=1)
     >>>     print "closest to %s:" % targetpoint
     >>>     print ("  at distance %g: point %s"
     ...            % (sqrt(closest[0][0]), points[closest[0][1]]))

    @ivar pointList: list of point coordinates: either a reference to the
       points argument (if it's a list) or an ordered list of the values of the
       points argument (if it's a dict).
    @ivar pointLabels: (only if the points parameter is a dict) keys of the
       points argument in the same order as its values in pointList
    @ivar tree: list of nodes of the KD-tree.
       Each node is being represented by a tuple. A node is either a leaf or
       a non-leaf node. Each node is a list of five values:

       leaf node: [ [up to <leafsize> point ids], None, None, 0, 0], where
       point ids are indices in the self.pointList (and possibly
       self.pointLables) - lists

       non-leaf node: [None, left child hrect, right child hrect, left child
       node idx, right child node idx], hrect means "hyper rectangle" a synonym
       for bounding box.

    @author:
    Modified a recipe from this source:
    http://www.scipy.org/Cookbook/KDTree
    Copyleft 2008 Sturla Molden
    University of Oslo
    """

    def __init__(self, points, leafsize=10, dim=3):
        """
        build a kd-tree for O(n log n) nearest neighbour search

        @param points: list of points, being float-3tuples. Might also be a
           dictionary {point label: [point coords]}, e.g. like
           abq_model_02.Model.nodeCoords. In this case it's always this
           label that is being returned by the query functions
           (knnSearch, radiusSearch) instead of the point index.
        @param leafsize: max. number of data points to leave in a leaf
        @param dim: number of coords per point

        @Note: In case of points being just a list of points the object keeps
           a reference of the provided points list. Don't modify while this
           object is in use. In case of points being a dict a new list is
           created so nothing to worry about in the first place. But the
           contents of this list (the individual points) are the points
           supplied by the dict. So their individual coordinates must not be
           modified!

        @Note: You should check beforehand that the points list is not empty.
           Otherwise an assert command will fail.
        """

        # in case points is a dict store the point labels and convert points
        # to a list of coordinates
        if isinstance(points, dict):
            self.pointLabels = sorted(points)
            points = [points[i] for i in self.pointLabels]

        if len(points)<=0:
            raise ValueError("KDTree.__init__ needs at least one point!")

        ndata = len(points)

        # bounding box
        hrect = BoundingBox(dim=dim)
        hrect.update(points)

        #-- create root of kd-tree
        # stored in the list tree which later becomes instance attribute, see
        # class doc string
        idx = range(ndata)
        idx.sort(key=lambda point_id: points[point_id][0])

        # split along x-axis into equally large chunks
        splitval = points[idx[ndata/2]][0]

        left_hrect = [hrect[0][:], hrect[1][:]]
        right_hrect = [hrect[0][:], hrect[1][:]]
        left_hrect[1][0] = splitval
        right_hrect[0][0] = splitval

        tree = [[None, left_hrect, right_hrect, None, None]]

        if ndata<2:
            # in case of just one point we have to cheat and add the same point
            # twice in order to have two branches on the root
            stack = [(idx, 1, 0, True),
                     (idx, 1, 0, False)]
        else:
            stack = [(idx[:ndata/2], 1, 0, True),
                     (idx[ndata/2:], 1, 0, False)]

        #-- recursively split data in halves using hyper-rectangles:
        while stack:

            #-- pop data off stack
            didx, depth, parent, leftbranch = stack.pop()
            ndata = len(didx)
            nodeptr = len(tree)

            #-- update parent node
            parentNode = tree[parent]
            if leftbranch:
                parentNode[3] = nodeptr  # left branch := nodeptr
            else:
                parentNode[4] = nodeptr  # right branch := nodeptr

            #-- insert node in kd-tree

            # leaf node?
            if ndata <= leafsize:
                # append a leaf to the tree with a copy (!) of didx
                tree.append([didx[:], None, None, 0, 0])

            # not a leaf, split the data in two
            else:
                # split along a different axis each turn
                splitdim = depth % dim

                # in order to split now sort along split axis
                idx = sorted(didx,
                             key=lambda point_id: points[point_id][splitdim])
                nodeptr = len(tree)
                stack.append((idx[:ndata/2], depth+1, nodeptr, True))
                stack.append((idx[ndata/2:], depth+1, nodeptr, False))

                # split at splitval into equally large chunks
                splitval = points[idx[ndata/2]][splitdim]

                # hrect == hyper rectangle == bounding box
                # == [lower left front corner, upper right back corner]
                if leftbranch:
                    # initialize left and right hrect with values from
                    # the parents left hrect
                    parent_hrect = parentNode[1]
                else:
                    # initialize left and right hrect with values from
                    # the parents right hrect
                    parent_hrect = parentNode[2]

                # first initialize new left and right hrects with parent hrect
                left_hrect = [parent_hrect[0][:], parent_hrect[1][:]]
                right_hrect = [parent_hrect[0][:], parent_hrect[1][:]]
                # then update one dimension of each side with the splitval
                left_hrect[1][splitdim] = splitval
                right_hrect[0][splitdim] = splitval
                # append (a non-leaf) node to tree
                tree.append([None, left_hrect, right_hrect, None, None])

        self.pointList = points  # keep a reference of the points coordinates
        self.tree = tree
        self.dim = dim

    @staticmethod
    def _intersect(hrect, r2, centroid):
        """
        checks if the hyperrectangle hrect intersects with the
        hypersphere defined by centroid and r2

        r2 is the square of the radius of the sphere
        """
        p = map(max, hrect[0], centroid)
        p = map(min, hrect[1], p)
        return sum(map(lambda x,y:(x-y)**2, p, centroid)) < r2

    @staticmethod
    def _intersect_rect(hrect1, hrect2):
        """
        checks if the two hyperrectangles hrect1 and hrect2 intersect::

         (a2-b1)*(a1-b2)<0      (i.e. a2-b1 and a1-b2 have different sign)
         for (a1, a2) := intervall from hrect1
             (b1, b2) := intervall from hrect2
         in each coordinate x,y,z

        @Note: Not needed for the initial purpose (self.boxSearch)...
        """
        return all(
            ((max1-min2)*(min1-max2))<0
            for min1,max1,min2,max2
            in izip(hrect1[0],hrect1[1],hrect2[0],hrect2[1]))

    def _quadraticKnnSearch(self, targetpoint, lidx, K):
        """ find K nearest neighbours of targetpoint
        among [self.pointList[idx] for idx in lidx]

        This function could be integrated into knnSearch for speed. It's only
        called from there at one occasion within the main loop.
        """
        ndata = len(lidx)
        if K>=ndata:
            K=ndata
        points = self.pointList
        sqd = [(sum(map(lambda x,y:(x-y)**2, points[idx], targetpoint)), idx)
               for idx in lidx]
        sqd.sort()
        return sqd[:K]

    def knnSearch(self, targetpoint, K):
        """Find the K nearest neighbours of targetpoint in this kdtree.

        @param targetpoint: coordinates triple of the points to find neighbours
        for
        @param K: number of neighbours to find.

        @Returns: a list of K (square distance, point index)-tuples for the K
        nearest points in the KDTree to the given targetpoint. The point index
        (second item in each tuple) is the index of the point in the pointlist
        delivered to the KDTree-constructor.

        If the points argument to self.__init__() was a dict then return the
        corresponding point labels instead of the point indices.

        @Note: Always returns a list of tuples even if K=1. No distinction made
        for performance reasons.
        """
        tree = self.tree
        stack = [tree[0]]
        knn = [(1E300, None)]*K
        while stack:

            leaf_idx, left_hrect, right_hrect, left, right = stack.pop()
            # leaf
            if leaf_idx is not None:
                _knn = self._quadraticKnnSearch(targetpoint, leaf_idx, K)
                if _knn[0][0] < knn[-1][0]:
                    knn = sorted(knn + _knn)[:K]

            # not a leaf
            else:

                # check left branch
                if self._intersect(left_hrect, knn[-1][0], targetpoint):
                    stack.append(tree[left])

                # check right branch
                if self._intersect(right_hrect, knn[-1][0], targetpoint):
                    stack.append(tree[right])

        if hasattr(self, "pointLabels"):
            knn = [(sqd, self.pointLabels[idx])
                   for sqd, idx in knn]
        return knn

    def radiusSearch(self, targetpoint, radius):
        """find all points within radius of targetpoint

        @Returns: a list of (square distance, point index)-tuples for all
        points of the KDTree within the given sphere. The point index
        (second item in each tuple) is the index of the point in the pointlist
        delivered to the KDTree-constructor. The list not sorted.

        If the points argument to self.__init__() was a dict then return the
        corresponding point labels instead of the point indices.
        """

        tree = self.tree
        points = self.pointList
        r2 = radius**2
        stack = [tree[0]]
        inside = []
        while stack:

            leaf_idx, left_hrect, right_hrect, left, right = stack.pop()

            # leaf
            if leaf_idx is not None:

                sqd = [(
                       sum(map(lambda x,y:(x-y)**2, points[idx], targetpoint)),
                       idx
                       ) for idx in leaf_idx]
                inside.extend(
                    (d2, idx) for d2, idx in sqd if d2<=r2)

            # not a leaf
            else:

                # check left branch
                if self._intersect(left_hrect, r2, targetpoint):
                    stack.append(tree[left])

                # check right branch
                if self._intersect(right_hrect, r2, targetpoint):
                    stack.append(tree[right])

        if hasattr(self, "pointLabels"):
            inside = [(sqd, self.pointLabels[idx])
                      for sqd, idx in inside]

        return inside

    def boxSearch(self, hrect):
        """find all points within the given hyper-rectangle (i.e. box)

        @param hrect: tuple of two lists. The first states the minimum
        coordinates of the hyper-rectangle and the second states the maximum
        coordinates. I.e. for three dimensions it's
        [ [x_min,y_min,z_min], [x_max,y_max,z_max]]

        @Returns: a list of point indexes for all points of the KDTree within
        the given hyper-rectangle. The point index is the index of the point in
        the pointlist delivered to the KDTree-constructor. The list is not
        sorted.

        If the points argument to self.__init__() was a dict then return the
        corresponding point labels instead of the point indices.
        """
        tree = self.tree         # abbrev
        points = self.pointList  # abbrev
        stack = [(tree[0], 0)]   # tree-node, depth
        inside = []
        while stack:

            (leaf_idx, left_hrect,right_hrect,left,right), depth = stack.pop()
            splitdim = depth % self.dim

            # leaf
            if leaf_idx is not None:
                inside.extend(
                    idx for idx in leaf_idx if all(
                        (px>=minx and px<=maxx)
                        for px, minx, maxx in izip(
                            points[idx], hrect[0], hrect[1])))

            # not a leaf
            else:

                # note: left_hrect[1][splitdim] == right_hrect[0][splitdim]
                # == split value!

                # check left branch
                if left_hrect[1][splitdim] >= hrect[0][splitdim]:
                    stack.append((tree[left], depth+1))

                # check right branch
                if right_hrect[0][splitdim] <= hrect[1][splitdim]:
                    stack.append((tree[right], depth+1))

        if hasattr(self, "pointLabels"):
            inside = [self.pointLabels[idx] for idx in inside]

        return inside

    def searchBelow(self, targetpoint):
        """Find all points having smaller or equal coordinates in each
        dimension.

        @Returns: a list of point indexes for all points of the KDTree below
        targetpoint. A point index is the index of the point in the pointlist
        delivered to the KDTree-constructor. The list not sorted.

        If the points argument to self.__init__() was a dict then return the
        corresponding point labels instead of the point indices.
        """

        tree = self.tree
        points = self.pointList
        dim = self.dim

        hibound = tree[0][2][1]
        lobound = tree[0][1][0]
        if any(b>t for b, t in izip(lobound, targetpoint)):
            # whole box is above targetpoint
            stack = False
            found = []
        else:
            takeAll = [b<=t for b, t in izip(hibound, targetpoint)]
            if all(takeAll):
                # whole box is below targetpoint
                stack = False
                found = range(len(points))
            else:
                # box intersects
                # store tree node and take-all-flags, depth
                stack = [(tree[0], takeAll, 1)]
                found = []

        while stack:

            (leaf_idx,left_hrect,right_hrect,left,right), takeAll, depth \
                = stack.pop()

            # leaf
            if leaf_idx is not None:
                if takeAll is True:
                    found.extend(leaf_idx)
                else:
                    found.extend(
                        idx for idx in leaf_idx
                        if all(x<y for x,y in izip(points[idx], targetpoint)))

            # not a leaf
            elif takeAll is True:
                stack.append((tree[left], True, depth+1))
                stack.append((tree[right], True, depth+1))
            else:
                # splitted along a different axis each turn
                splitdim = depth % dim

                # separation along this axis is irrelevant already
                if takeAll[splitdim]:
                    stack.append((tree[left], takeAll, depth+1))
                    stack.append((tree[right], takeAll, depth+1))

                addedTakeAll = takeAll[:]
                addedTakeAll[splitdim] = True
                addedTakeAll = all(addedTakeAll) or addedTakeAll

                # check left branch
                if left_hrect[1][splitdim] <= targetpoint[splitdim]:
                    stack.append((tree[left], addedTakeAll, depth+1))
                else:
                    stack.append((tree[left], takeAll, depth+1))

                # check right branch
                if right_hrect[0][splitdim] <= targetpoint[splitdim]:
                    stack.append((tree[right], takeAll, depth+1))
                # if lower bound of right > targetpoint then branch is above

        if hasattr(self, "pointLabels"):
            found = [self.pointLabels[idx] for idx in found]

        return found


#------------------------------------------------------------------------------

class InvDistTree(object):
    """inverse-distance-weighted interpolation using KDTree:

    >>> from bae.misc_01 import InvDistTree
    >>> invdisttree = InvDistTree( X, z )  -- data points, values
    >>> interpol = invdisttree.interpolate( q, nnear=3, p=1)

    interpolates z from the 3 points nearest each query point q;

    The example above finds the 3 data points nearest q, at distances d1 d2 d3
    and returns the IDW average of the values z1 z2 z3::
        (z1/d1 + z2/d2 + z3/d3)
        / (1/d1 + 1/d2 + 1/d3)
        = .55 z1 + .27 z2 + .18 z3  for distances 1 2 3

    p: use 1 / distance**p as weighting factor for the values at different
    points

    How many nearest neighbors should one take ?
     - start with 8 11 14 .. 28 in 2d 3d 4d .. 10d; see Wendel's formula
     - make 3 runs with nnear= e.g. 6 8 10, and look at the results --
     |interpol 6 - interpol 8| etc., or |f - interpol*| if you have f(q).
     I find that runtimes don't increase much at all with nnear -- ymmv.

    p=1, p=2 ?
     - p=2 weights nearer points more, farther points less.
     - In 2d, the circles around query points have areas ~ distance**2,
     so p=2 is inverse-area weighting. For example::
        (z1/area1 + z2/area2 + z3/area3)
        / (1/area1 + 1/area2 + 1/area3)
        = .74 z1 + .18 z2 + .08 z3  for distances 1 2 3
    Similarly, in 3d, p=3 is inverse-volume weighting.

    Scaling:

    If different X coordinates measure different things, Euclidean distance
    can be way off.  For example, if X0 is in the range 0 to 1
    but X1 0 to 1000, the X1 distances will swamp X0;
    rescale the data, i.e. make X0.std() ~= X1.std() .

    A nice property of IDW is that it's scale-free around query points:
    If I have values z1 z2 z3 from 3 points at distances d1 d2 d3,
    the IDW average::
       (z1/d1 + z2/d2 + z3/d3)
       / (1/d1 + 1/d2 + 1/d3)
    is the same for distances 1 2 3, or 10 20 30 -- only the ratios matter.
    In contrast, the commonly-used Gaussian kernel exp( - (distance/h)**2 )
    is exceedingly sensitive to distance and to h.


    Taken from
    http://stackoverflow.com/questions/3104781/inverse-distance-weighted-idw-interpolation-with-python
    provided by http://stackoverflow.com/users/86643/denis

    ... and then modified a bit.
    """
# anykernel( dj / av dj ) is also scale-free
# error analysis, |f(x) - idw(x)| ? todo: regular grid, nnear ndim+1, 2*ndim

    def __init__(self, X, z):
        """
        @param X: point list, list of coordinate triples
        @param z: list of values: floats or vectors (lists of floats)
        """
        assert len(X) == len(z), "len(X) %d != len(z) %d" % (len(X), len(z))
        self.tree = KDTree(X)  # build the tree
        self.z = z

    def interpolateSingle(self, point, nnear=6, p=3):
        """Return the interpolated value at point, considering the nnear
        nearest neighbours. The corresponding values are weighted by
        weight = 1/(distance^p)
        """
        # nnear nearest neighbours of each query point --

        # note "sqd" means "square distance"
        sqdExp = 0.5*p
        sqdPtidList = self.tree.knnSearch(point, K=nnear)
        if sqdPtidList[0][0] < 1e-10 or nnear==1:
            # query point equals given data point
            # or consider only nearest point
            # => no interpolation
            wz = self.z[sqdPtidList[0][1]]
        else:
            w = [1 / sqd**sqdExp for sqd, pointId in sqdPtidList]
            ws = 1.0 / sum(w)
            wz = sum(wi*self.z[sqdPtidList[i][1]] for i, wi in enumerate(w)
                     ) * ws
#         print "*** Pt %s: sqdPtidList=%s" % (point, sqdPtidList)
#         print "*** pts: %s" % [self.tree.pointList[pd[1]] for pd in sqdPtidList]
#         print "*** w=%s, ws=%s, wz=%s" % (w, ws, wz)

        return wz

    def interpolate(self, points, nnear=6, p=3, msgTicker=None, radius=None):
        """Return the interpolated value at each of the given points,
        considering the nnear nearest neighbours. The corresponding values
        are weighted by weight = 1/(distance^p)

        For many points this should be faster then calling
        self.interpolateSingle() multiple times.

        @param radius: optional search radius. If specified disregard points
        further away.

        @param msgTicker: may be an instance of the MsgTicker class to provide a
        progress indicator. Its msg() method will be called with the number of
        values already interpolated. Does not work for nnear=1.
        """
        # nnear nearest neighbours of each query point --

        if radius:
            searchF = functools.partial(self.tree.radiusSearch, radius=radius)
        else:
            searchF = functools.partial(self.tree.knnSearch, K=nnear)

        # abbreviations
        # note "sqd" means "square distance"
        nan = float("nan")
        ptVals = self.z
        sqdExp = 0.5*p  # exponent to weight the square distance values

        # point search: iterator of (sqare distance, point id)-tuples
        sqdIdsList = (searchF(point) for point in points)
        if nnear==1:
            # consider only nearest point => no interpolation
            # if no near point then sqdIds is empty, and then return nan,
            # otherwise value of nearest point
            values = [((len(sqdIds) and ptVals[sqdIds[0][1]]) or nan)
                      for sqdIds in sqdIdsList]
        else:
            if radius:
                # filter points too far away
                sqRad = radius**2
                sqdIdsList = (
                    [sqdId for sqdId in sqdIds if sqdId[0]<=sqRad]
                    for sqdIds in sqdIdsList)

            dummyVal = ptVals[0]
            if isinstance(dummyVal, float):
                N = 0
            else:
                N = len(dummyVal)
            values = []
            for cnt, sqdIds in enumerate(sqdIdsList):
                if not sqdIds:
                    # no point in range as specified by radius argument
                    values.append(nan)
                elif sqdIds[0][0] < 1e-10 or len(sqdIds)==1:
                    # query point equals given data point => no interpolation
                    # to prevent division by zero in the weights
                    values.append(ptVals[sqdIds[0][1]])
                else:
                    weights = [1.0 / sqd**sqdExp
                               for sqd, pointId in sqdIds]
                    ws = 1.0 / sum(weights)
                    if N:
                        # averaging vector
                        res = [0.0]*N
                        for wi, sqdPtid in izip(weights, sqdIds):
                            for i, x in enumerate(ptVals[sqdPtid[1]]):
                                res[i] += ws*wi*x
                        # store vector result
                        values.append(res)
                    else:
                        # averaging scalars
                        values.append(ws * sum(
                            wi*ptVals[sqdPtid[1]]
                            for wi, sqdPtid in izip(weights, sqdIds)))

                if msgTicker is not None:
                    msgTicker.msg(cnt+1)
        return values

#} end of point search: find closest points, KDTree, inverse distance interpol

###############################################################################
###############################################################################
###############################################################################
#   ... ?
#


###############################################################################
###############################################################################
###############################################################################
#{ operating system tasks: file transfer, job control
#

def RemoteFile(filename, mode='r'):
    """Not working properly under all circumstances yet. Sorry...
    """
    hostname, filename = filename.split(":")

    if mode=='r':
        retcode = subprocess.call(["ssh", hostname, "test", "-f", filename])
        if retcode != 0:
            raise IOError("No such file or directory '%s:%s'"
                          % (hostname, filename))
        p = subprocess.Popen(["ssh", hostname, "cat", filename],
                             stdout=subprocess.PIPE)
        return p.stdout

    else:
        raise Exception("ERROR RemoteFile: Only mode='r' implemented so far.")


#------------------------------------------------------------------------------

def waitForMarker(markerFile, appear, waitPeriod=60, cntMsgLaunch=10,
                  logfile=None, mode=os.F_OK, addFullTime=False):
    """Stop the current program until a certain file appears or disappears.

    >>> from bae.misc_01 import waitForMarker
    >>> os.system("abaqus j=test")
    >>> print "Starting the job."
    >>> waitForMarker("test.lck", appear=True, cntMsgLaunch=None)
    >>> print "Job now running, waiting for it to finish..."
    >>> waitForMarker("test.lck", "disappear")
    >>> print "Job finished."

    @param markerFile: file name of the file to check
    @param appear: True or "appear": wait for markerFile to appear otherwise
    wait for it disappearing
    @param cntMsgLaunch: Every cntMsgLaunch lookups of the marker file a msg
    will be printed to logfile. If 0 or None no message will be printed.
    @param waitPeriod: every waitPeriod seconds the markerFile will be tested
    @param logfile: DEPRECATED argument being ignored. use the bae.log_01
    functionality for diagnostic output
    @param mode: mode argument to os.access or a corresponding string,
    i.e. "F_OK" or "R_OK"
    """

    if appear is True or appear=="appear":
        waitformarker = True
        waitoptionstr = 'appear'
    else:
        waitformarker = False
        waitoptionstr = 'disappear'

    if isinstance(mode, basestring):
        try:
            mode = getattr(os, mode)
        except AttributeError:
            raise ValueError("Illegal mode argument string to waitForMarker.")

    cnt = cntMsgLaunch
    if cntMsgLaunch is None or cntMsgLaunch==0:
        cnt = 0
        cntMsgLaunch = -1

    while 1:
        markerthere = os.access(markerFile, mode)
        if markerthere == waitformarker:
            break
        if cnt==cntMsgLaunch:
            if addFullTime:
                msg('%s, waiting for file %s to %s.'
                    % (time.asctime(),markerFile,waitoptionstr))
            else:
                msg("waiting for file %s to %s."
                    % (markerFile,waitoptionstr))
            cnt = 0
        cnt += 1
        time.sleep(waitPeriod)
    return


def waitForProcess(pid, waitPeriod=60, cntMsgLaunch=10):
    """Pause the current program until another process with the given pid
    finishes. The other process does not need to be a child of this process.

    Use command line tools like "top -u gero" or "ps -eaf | grep abaqus" to
    get the pid of the process you're after.

    >>> from bae.misc_01 import waitForProcess
    >>> os.system("abaqus j=test")
    >>> print "Starting the job."
    >>> pid = int(raw_input("Enter PID of the abaqus job:"))
    >>> waitForProcess(pid)
    >>> print "Job finished."

    @param pid: process id
    @param cntMsgLaunch: Every cntMsgLaunch lookups of the marker file a msg
    will be printed to logfile. If 0 or None no message will be printed.
    @param waitPeriod: every waitPeriod seconds the markerFile will be tested
    """

    cnt = cntMsgLaunch
    if cntMsgLaunch is None or cntMsgLaunch==0:
        cnt = 0
        cntMsgLaunch = -1

    procMarker = "/proc/%d" % pid
    while 1:
        if not os.path.isdir(procMarker):
            break
        if cnt==cntMsgLaunch:
            msg("Waiting for process %d to finish." % pid)
            cnt = 0
        cnt += 1
        time.sleep(waitPeriod)
    return


#------------------------------------------------------------------------------

def getCurrentOdbFrame(staFileName):
    """Read the frame number of the last frame written to the odb from the
    given .sta file.

    @returns: The frame number as an integer. Or zero if no appropriate line
       could be found in the sta file.
    """
    sedCmd = r"s/.*Number[ \t]*([0-9]+)[ \t].*/\1/"
    cmd = ("(cat %s | grep 'ODB Field Frame Number' | tail -n 1 | sed -r '%s')"
           " 2> /dev/null" % (staFileName, sedCmd))
    output = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE).communicate()[0]
    if not output:
        return 0
    else:
        return int(output)


def waitForFrame(staFileName, frameNb, waitPeriod=60):
    r"""Wait until the given frame number is written to the odb.

    DEPRECATED: Better use waitForStepFrame!

     >>> waitForFrameNb = 100
     >>> staFileName = odbName.replace(".odb", ".sta")
     >>>
     >>> print "Waiting for frame %d in %s" % (waitForFrameNb, staFileName)
     >>> sys.stdout.flush()
     >>>
     >>> waitForFrame(staFileName, waitForFrameNb)
     >>>
     >>> print "Frame %d finished. Now proceeding." % waitForFrameNb
     >>> sys.stdout.flush()

    @param staFileName: .sta-file to check which frames have been written so
    far.
    @param frameNb: wait until this frame appears in the odb (according to the
    .sta-file)
    @param waitPeriod: time interval to check the .sta file
    @note: This function does not check the actual odb. It does not check for
    the correct step either. If we are already past the given frame number,
    we won't wait.
    @note: At an odb step change this does not work well: Suppose you want to
    wait for frame 1 of step 3 and the odb currently has got the last frame of
    step 2. Then getCurrentOdbFrame returns the number of this last frame which
    is very likely to be greater than 1 and this function then assumes that
    frame one is already in the odb and continues without waiting for the
    actual step 3 frame 1 that you are waiting for.
    """
    while True:
        currentFrameNb = getCurrentOdbFrame(staFileName)
        if currentFrameNb>=frameNb:
            break
        time.sleep(waitPeriod)
    return


class StaFileParser(object):
    """Parse an sta file. Used to wait for specific frames to appear.

    Usage:
     >>> # note: this will block the program until the sta file appears
     >>> staFile = StaFileParser( odbName.replace(".odb", ".sta") )
     >>>
     >>> for stepNb, frameNb in stepFrameIter(framesList):
     >>>     msg("Waiting for step %d frame %d" % (stepNb, frameNb))
     >>>     staFile.waitForStepFrame(stepNb, frameNb)
     >>>     msg("Processing step %d frame %d" % (stepNb, frameNb))
     >>>     do_something()

    @note: Alternatively you can use the convenience function
    L{bae.misc_01.waitForStepFrame}.
    """
    def __init__(self, staFileName, waitPeriod=60, blocking=True):
        r"""
        @param staFileName: .sta-file to check which frames have been written so
        far.
        @param waitPeriod: time interval to check the .sta file
        @param blocking: if True then block the program until the given sta
        file appears and is readable. Otherwise just try to open the sta file
        which might raise an OSError.
        """
        while blocking and not os.access(staFileName, os.R_OK):
            time.sleep(waitPeriod)

        self.staFile = open(staFileName, "r")
        self.waitPeriod = waitPeriod
        self.currentStep = 0
        self.currentFrame = None

    def readline(self):
        "Read one line from the sta file. Wait until a new line pops up."
        while True:
            x = self.staFile.readline()
            if x!="":
                return x
            time.sleep(self.waitPeriod)

    def findNextFrame(self):
        "Read lines until you find the next frame."
        while True:
            line = self.readline()
            if line.startswith("Abaqus/Explicit"):
                msg("Found end step mark.", debugLevel=10)
                while not line.startswith(" STEP "):
                    line = self.readline()
                    continue
                self.currentStep = int(line.split(None,2)[1])
                self.currentFrame = None
                msg("Found step %d" % self.currentStep, debugLevel=10)
            elif line.startswith("ODB Field Frame Number"):
                line = line[24:]
                self.currentFrame = int(line.split(None,1)[0])
                msg("Found frame %d" % self.currentFrame, debugLevel=10)
                break

    def waitForStepFrame(self, stepNb, frameNb):
        "Read lines until you find the specified step frame. Or a later frame."
        while (
                (self.currentStep < stepNb)
                or (self.currentStep==stepNb and self.currentFrame<frameNb)):
            lastFrameNb = self.currentFrame
            self.findNextFrame()
            if self.currentStep > stepNb:
                msg("WARNING: When waiting for step %s frame %s found step %s."
                    " Last reported frame for the previous step is %s."
                    % (stepNb, frameNb, self.currentStep, lastFrameNb))

    def getTimeSeries(self):
        """Parse the Abaqus Explicit .sta-file and return lists for increment
        numbers, totaltime and CPU time (in seconds):
         >>> import csv
         >>> from itertools import izip
         >>> from bae.misc_01 import StaFileParser
         >>> from bae.log_01 import msg
         >>>
         >>> inputname = "NPM2016SLC_R02_G01_S01_Q02_M10_D01.sta"
         >>> outputName = "times.csv"
         >>> msg("Reading file %s." % inputname)
         >>> parser = StaFileParser(inputname, blocking=False)
         >>> (increment, totalTime, cpuTime) = parser.getTimeSeries()
         >>> msg("...found %d time points." % len(increment))
         >>>
         >>> # calculate deltas
         >>> cpuTime = [(x-cpuTime[0]) for x in cpuTime]
         >>> deltaIncr = [0]
         >>> deltaIncr.extend((x-y) for x, y in izip(increment[1:], increment[:-1]))
         >>> deltaTT = [0.0]
         >>> deltaTT.extend((x-y) for x, y in izip(totalTime[1:], totalTime[:-1]))
         >>> deltaCPU = [0]
         >>> deltaCPU.extend((x-y) for x, y in izip(cpuTime[1:], cpuTime[:-1]))
         >>>
         >>> fo = csv.writer(open(outputName, "wb"))
         >>> fo.writerow(["CPU/s", "Increment", "Totaltime", "d Incr/ d CPU", "d TT/ d CPU"])
         >>> for t, incr, tt, dInc, dTt, dCpu in izip(
         >>>         cpuTime, increment, totalTime, deltaIncr, deltaTT, deltaCPU):
         >>>     fo.writerow([t, incr, tt,
         >>>                  float(dCpu>0) and float(dInc)/dCpu,
         >>>                  float(dCpu>0) and dTt/dCpu])
         >>> del fo
         >>> msg("Wrote results to %s" % outputName)

        """
        # skip beginning
        for line in self.staFile:
            if line.startswith(" SOLUTION PROGRESS"):
                break

        # initialize result lists
        increment = list()
        totalTime = list()
        cpuTime = list()

        # now parse each line
        rex = re.compile(
            # increment number
            r"^\s+(?P<incr>\d+)"
            # step time
            r"\s+(?P<steptime>[0-9.E+-]+)"
            # total time
            r"\s+(?P<totaltime>[0-9.E+-]+)"
            # cpu time
            r"\s+(?P<cputime>[0-9:]+)")
        for line in self.staFile:
            res = rex.match(line)
            if not res:
                continue
            increment.append(int(res.group("incr")))
            totalTime.append(float(res.group("totaltime")))
            tTup = [int(x) for x in res.group("cputime").split(":")]
            cpuTime.append(3600*tTup[0] + 60*tTup[1] + tTup[2])
        return (increment, totalTime, cpuTime)


_staFileParserDict = dict()

def waitForStepFrame(staFileName, stepNb, frameNb, waitPeriod=None):
    r"""Wait until the given frame number is written to the odb.

     >>> waitForStepNb = 3
     >>> waitForFrameNb = 10
     >>> staFileName = odbName.replace(".odb", ".sta")
     >>>
     >>> msg("Waiting for step %d frame %d in %s"
     >>>     % (waitForStepNb, waitForFrameNb, staFileName))
     >>>
     >>> waitForStepFrame(staFileName, waitForStepNb, waitForFrameNb)
     >>>
     >>> msg("Frame %d finished. Now proceeding." % waitForFrameNb)

    @param staFileName: .sta-file to check which frames have been written so
    far.
    @param stepNb: step number to wait for
    @param frameNb: wait until this frame (of step stepNb) appears in the odb
    (according to the .sta-file)
    @param waitPeriod: time interval to check the .sta file. Defaults to 60sec.
    @note: This function does not check the actual odb but the given sta file.
    If we are already past the given frame number, we won't wait.
    @note: This function involves some magic in how it treats the sta file:
    It keeps the sta file open between subsequent calls to this function and
    later only reads newly appended lines. You should not do dirty tricks with
    the sta file on your side, this may not go well with this approach.
    @note: On subsequent calls to this function the sta file is identified by
    staFileName. If you change the current working directory you should usually
    not update staFileName to reflect the now different relative path unless
    you want to instatiate a new file descriptor which will then reread the sta
    file.
    @note: All sta files stay open until the end of the program. There is
    currently no simple means to close them.
    """

    # first find the correct StaFileParser object
    global _staFileParserDict
    try:
        parser = _staFileParserDict[staFileName]
    except KeyError:
        if waitPeriod is None:
            waitPeriod=60
        parser = StaFileParser(staFileName, waitPeriod)
        _staFileParserDict[staFileName] = parser
    else:
        if waitPeriod is not None:
            parser.waitPeriod = waitPeriod

    # actually wait (or not)
    parser.waitForStepFrame(stepNb, frameNb)


#------------------------------------------------------------------------------

def ftpGet(remotePath, login, passwd,
           localDestFolder=".",
           waitPeriod=300, cntMsgLaunch=12, logfile=sys.stdout,
           ftpServer="ftp.beckengineering.com.au"):
    """Try to download the file remotePath from the ftp server. Wait if it's
    not there (yet).

    >>> from bae.misc_01 import ftpGet
    >>> print "Waiting for new data on the ftp server."
    >>> fileName = "AIRGAP_YR2012_M04.txt"
    >>> ftpGet("Perseverance/2010/Coupling/%s" % fileName,
    ...        "gero", "AlleMeineEntchen")
    >>> f = open(fileName, "r")
    >>> print "Loaded %s to the current working directory." % fileName

    @param remotePath: complete path and filename on the ftp server
    @param login: ftp login name
    @param passwd: password for this login
    @param localDestFolder: where to copy the file from the ftp server.
    "." means into the current working directory (os.getcwd()). None or "" means
    don't get it at all, just check if it's there.
    @param waitPeriod: every waitPeriod seconds the ftp server will be checked.
    @param cntMsgLaunch: Every cntMsgLaunch lookups of the requested file a msg
    will be printed to logfile. If 0 or None no message will be printed.
    @param logfile: an open file object where to print the message to
    @param ftpServer: name of the ftp server

    @Note: The destination file will be overwritten silently.
    @Note: The remotePath argument contains the complete path, i.e. folder and
    file name. Only the file name without the folder is appended to the
    localDestFolder argument (with the proper "/" path separator in between)
    to yield the path name of the destination, i.e. where to write this file
    to on the local maschine.
    """

    cnt = cntMsgLaunch
    if cntMsgLaunch is None or cntMsgLaunch==0:
        cnt = 0
        cntMsgLaunch = -1

    remoteDir, fileName = os.path.split(remotePath)
    if localDestFolder:
        localDestPath = os.path.join(localDestFolder, fileName)

    while 1:

        # login to ftp server
        ftp = ftplib.FTP(ftpServer)
        ftp.login(login, passwd)

        # get list of files with at most one item: the file we are looking for
        fileList = list()
        ftp.retrlines("NLST %s" % remotePath, fileList.append)
        if len(fileList)>0:

            # file found on server, now get it
            fout = open(localDestPath, "wb")

            try:
                ftp.retrbinary("RETR %s" % remotePath, fout.write)
            except ftplib.error_perm:
                # do nothing and try again in the next loop iteration
                logfile.write(
                    "%s, ftpGet: Oops, file %s found but then it disappeared\n"
                    % (time.asctime(),remotePath))
                logfile.flush()
                del fout
                os.remove(localDestPath)
                pass
            else:
                # all went well, file fetch to local destination

                del fout  # close the file
                logfile.write("Read %s from the ftp server\n"
                              % localDestPath)
                logfile.flush()
                # only way out of the loop (except exceptions)
                break

        # in any case: quit ftp connection
        ftp.quit()

        # possibly write a message to the logfile
        if cnt==cntMsgLaunch:
            logfile.write('%s, ftpGet: waiting for file %s to appear on the'
                          ' ftp server.\n' % (time.asctime(),remotePath))
            logfile.flush()
            cnt = 0
        cnt += 1
        time.sleep(waitPeriod)
    return


#------------------------------------------------------------------------------

def _ftpPut_chkpath(connection, path):
    """function that will be called recursively by ftpPut to check and
    where applicable create all directories on the ftp server.
    """
    if not path or path==".":
        return

    head, tail = os.path.split(path)
    _ftpPut_chkpath(connection, head)
    _ftpPut_chkdir(connection, path)


def _ftpPut_chkdir(connection, path):
    """check directory path on the ftp connection for existence and create
    if applicable
    """
    currentdir = connection.pwd()
    try:
        connection.cwd(path)
    except ftplib.error_perm:
        connection.mkd(path)
    else:
        connection.cwd(currentdir)


def ftpPut(localPath, login, passwd,
           remoteDestFolder=".",
           ftpServer="ftp.beckengineering.com.au",
           msgTicker=None
           ):
    r"""Put one ore many files or directories on the ftp server.

    >>> from bae.misc_01 import ftpPut
    >>> fileName = "CaveShape_YR2012_M05.vtk"
    >>> ftpPut(fileName, "gero", "gerospassword",
    ...        "Perseverance/2010/Coupling")

    >>> from bae.misc_01 import ftpPut, MsgTicker
    >>> from sys import stdout as logfile
    >>> fileNames = ["README_mypngs.txt", "myPNGfolder"]
    >>> ticker = MsgTicker(logfile, "sending %s\n",
    ...                    firstPeriod=0, period=0)
    >>> ftpPut(fileNames, "gero", "gerospassword",
    ...        "Perseverance/2010/Coupling", msgTicker=ticker)
    >>> ticker.msgTemplate = "Ftp transfer finished.\n"
    >>> ticker.msg()
    >>> del ticker

    @param localPath: complete path and filename of the file to copy to the
    ftp server, or a list or other iterable of such path names. Must not start
    with a slash, only relative paths permitted!
    @param login: ftp login name
    @param passwd: password for this login
    @param remoteDestFolder: Directory on the ftp server to copy the file(s)
    to. "." or "" means into the root (home-) directory of the specified ftp
    account. This directory must already exist in order to write to it.
    (Whereas the subfolders as listed in localPAth will be created if
    necassary.) To create an ftp folder with python use module ftplib:
    ftp = ftplib.FTP(...), ftp.login(...), ftp.mkd(...), ftp.quit()
    @param ftpServer: name of the ftp server
    @param msgTicker: may be an instance of the MsgTicker class to provide a
    progress indicator. Its msg() method will be called with the local filename
    of the file currently processed. Note that a msgTicker with period>0 (the
    default) will not necessarily write a message for each file, so the
    messageTemplate should be something like "Ftp transfer in progress,
    currently processing %s\n". Or specify the period to 0 to get
    all messages. Since the messages are printed before the transfer actually
    starts, it might be advisable to post a final message.

    @note: The destination file will be overwritten silently.

    @note: The localPath argument contains the complete (relative) path(s),
    i.e. folder and file name. This path will be appended to the
    remoteDestFolder argument (with the proper "/" path separator in between)
    to yield the path name of the destination, i.e. where to write this file
    to on the ftp server.

    @note: Symbolic links in subfolders are ignored. Symbolic links directly
    stated in the localPath argument are transfered to the ftp server.
    """

    if isinstance(localPath, basestring):
        localPath = (localPath, )

    # login to ftp server
    connection = ftplib.FTP(ftpServer)
    connection.login(login, passwd)

    # change dir
    if remoteDestFolder and remoteDestFolder!=".":
        connection.cwd(remoteDestFolder)

    for path in localPath:

        if path.startswith(os.path.sep):
            print ("WARNING: Not putting %s on the ftp, it begins with a %s."
                   " Paths must be relative paths." % (path, os.path.sep))
            continue

        head, tail = os.path.split(path)
        _ftpPut_chkpath(connection, head)

        if os.path.isdir(path):
            for root, dirs, files in os.walk(path):
                # check remote directory
                _ftpPut_chkdir(connection, root)

                # store individual files
                for fileName in files:
                    fileName = os.path.join(root, fileName)
                    if not os.path.islink(fileName):
                        if msgTicker is not None:
                            msgTicker.msg(fileName)
                        connection.storbinary(
                            "STOR %s" % fileName, open(fileName, "rb"))

        else:  # path is a file
            if msgTicker is not None:
                msgTicker.msg(path)
            connection.storbinary(
                "STOR %s" % path, open(path, "rb"))

    connection.quit()
    return

#------------------------------------------------------------------------------

class Config(object):
    r"""Object for abaqus.env like variable configuration.
    It will search for a config file and execute the python code in that
    file.
    All files found upwards from the current folder up to the root-folder
    will be considered.

    Usage:
     >>> from bae.misc_01 import Config
     >>> config = Config("configBE.py")
     >>>
     >>> print config.thisVariableIsSetInsideConfigBE

    This code will search for a "configBE.py" file and execute the python code
    in that file. All files found upwards from the current folder up to
    the root-folder will be considered.

    I.e. if you start the script in /work/abaqus/MyProj/1_pre/9_doall/
    there must be a file configBE.py in this forlder or in
    /work/abaqus/MyProj/1_pre/ or in /work/abaqus/MyProj/ or in ...
    containing for example this code:
     >>> # this is the file configBE.py
     >>> thisVariableIsSetInsideConfigBE = "a value"

    If you don't want to use a config file you can also do:
     >>> from bae.misc_01 import Config
     >>> config = Config()
     >>> config.thisVariableIsSetInsideConfigBE = "a value"
     >>>
     >>> print config.thisVariableIsSetInsideConfigBE

    That is equivalent to:
     >>> from bae.misc_01 import Config
     >>> config = Config(thisVariableIsSetInsideConfigBE="a value")

    For path manipulation purposes there is an automatic variable available in
    the config file: configDir. So for example configBE.py could read:
     >>> import os.path
     >>> fileNameMeshSource = os.path.join(configDir, "10_MESHING/Konst.inp")

    """

    def __init__(self, *args, **kwargs):
        """
        If you specify a positional argument it will be the name of the
        config file.

        In this case there are two optional keyword arguments: printMessage
        and postData. All other keyword arguments are ignored in this case.
        I.e. if you specify a positional argument you can't pass additional
        variables as keyword arguments at the same time. You can however
        add additional attributes after initialization.

        If there is no positional argument then no file is being read.
        In this case all (optional) keyword arguments will be treated as
        initializer to self's attributes. (In this case printMessage and
        postData have no special meaning but will result in corresponding
        config-attributes.) I.e.:
         >>> config = Config(meshVersion="G01")
         >>> print config.meshVersion
         G01

        @kwarg printMessage: If True then write a message when reading a config
           file.
        @kwarg postData: if True then add contents of postdata file. The path
           for the post-data-file will be derived from self.odbPath. I.e.
           odbPath must be defined in the config file otherwise this option has
           no effect.
        """
        if len(args)==0:
            # optional arbitrary keyword args: Config behaves like Container
            self.__dict__.update(kwargs)
            return

        if len(args)>2:
            raise ValueError(
                "ERROR: %s.%s() accepts at most two positional arguments."
                " Instead %d were given:\n%s"
                % (__name__, self.__class__.__name__, len(args), args))

        # else: we do have one (or --deprecated-- two) positional arguments
        configFile = args[0]

        # deprecated feature: a second positional argument is treated as
        # printMessage argument
        try:
            printMessage = args[1]
        except IndexError:
            try:
                printMessage = kwargs["printMessage"]
            except KeyError:
                printMessage = False

        if configFile:
            # look for a file "configFile" above the current directory
            # exec every file found
            for c,d,f in self.walk_up(os.curdir):
                if configFile in f:
                    # execute commands in the execfile
                    if printMessage:
                        msg("Reading config file %s"
                            % os.path.join(c, configFile))
                    self.configDir = c
                    execfile(os.path.join(c, configFile), self.__dict__)
                    del self.configDir

        # add contents of postdata file if argument postData=True
        if kwargs.get("postData") and hasattr(self, "odbPath"):
            from bae.odb_04 import getOdbPostData
            if printMessage:
                msg("Trying to read post-data.")
            postData = getOdbPostData(self.odbPath)
            self.__dict__.update(postData.__dict__)
            if printMessage and not postData.__dict__:
                msg("Failed to read post-data for %s." % self.odbPath)

    def walk_up(self, bottom):
        """
        mimic os.walk, but walk 'up' instead of down the directory tree
        @param bottom: path to start the search
        """
        # Note: this could presumably be converted to an iteration instead of
        # the recursion.

        bottom = os.path.realpath(bottom)

        #get files in current dir
        try:
            names = os.listdir(bottom)
        except OSError as e:
            msg("Could not access %s.\n%s" % (bottom, e))
            return

        dirs = []
        nondirs = []
        for name in names:
            if os.path.isdir(os.path.join(bottom, name)):
                dirs.append(name)
            else:
                nondirs.append(name)

        new_path = os.path.realpath(os.path.join(bottom, '..'))

        # if we are not yet at the top then one more recursion step
        if new_path != bottom:
            for x in self.walk_up(new_path):
                yield x

        # give the results in current directory after the results further "up"
        yield bottom, dirs, nondirs


class Config_02(object):
    """This implementation is in Testing stage. Not production ready, yet.

    Configuration parameter container.

    The following attributes will be determined automatically:
    projectPrefix, jobName, runVersion, baseRunVersion,
    meshTypeSplitVersion, meshVersion, initStressVersion, seqVersion,
    matVersion, geoDomVersion,
    odbPath, projectBaseDir, runBaseDir (if not specified in config files).

    Then configBE.py files in  the current
    branch of the project tree will be
    located and parsed. All variables defined in those files will become
    attributes of self.

    Finallythe postData.py fileaccompanying the odb will be parsed in the
    same way, adding its attributes to self as well. You can specify runBaseDir in the config files: the folder under which
    we'll find the folder $(runVersion) and in it the odb.
    E.g. runBaseDir = "/mnt/boab2/data6/abaqus/Cadia2021028/2_RUN"This will be used to find the odb or postData or main input file
    and then further deduce the xxxVersion parameters and the job-name.

    runBaseDir will be searched for in the config files in a preliminary
    pass that only recognizes this single attribute. The config files will
    finally be evaluated a second time assuming that the automatic version
    parameters might be required to evaluate some of the expressions in the
    config files.

    If runBaseDir has not been specified in any config file then the odb is
    searched for in $(projectBaseDir)/2_RUN/$(runVersion).
    For path manipulation purposes there is an automatic variable available in
    the config file: configDir. So for example configBE.py could read:
     >>> import os.path
     >>> fileNameMeshSource = os.path.join(configDir, "10_MESHING/Konst.inp")
    """

    configFileName = "configBE.py"

    def __init__(self, printMessage=False):

        """
        @param printMessage: If True then write a message when reading a config
           file.
        """

        # rex_ prefix => "regular expression"
        rex_runfolder = re.compile(r"^R\d\d")

        # Walk up the tree to:
        # - collect configBE.py files
        # - determine self.projectBaseDir as parent of 1_PRE, 2_RUN or 3_POST
        # - determine self.runVersion="RXX" if cwd is subfolder of 1_PRE/RXX,
        #   2_RUN/RXX or 3_POST/RXX
        curpath = os.getcwd()  # doesn't contain a trailing slash
        configFiles = list()
        nextIsProjectBase = False
        while curpath:
            confFile = os.path.join(curpath, self.configFileName)
            if os.path.isfile(confFile):
                configFiles.append(confFile)

            if nextIsProjectBase:
                self.projectBaseDir = curpath
                break

            # next curpath and current dir and parent
            curpath, dirname = os.path.split(curpath)
            parentDir = os.path.basename(curpath)

            # identify runVersion by a folder name starting "R##"
            if rex_runfolder.match(dirname) and (
                    parentDir in ("1_PRE", "2_RUN", "3_POST")):
                self.runVersion = dirname

            if dirname in ("1_PRE", "2_RUN", "3_POST"):
                nextIsProjectBase = True

        # cleanup results of walking the folder tree
        configFiles.reverse()
        del nextIsProjectBase

        # check if runBaseDir is defined in config files:
        # parse config files once without all xxxVersions defined
        localsDict = dict()
        localsDict.update(self.__dict__)
        # As a first guess take the parent folder of 2_RUN as projectPrefix
        localsDict["projectPrefix"] = os.path.basename(self.projectBaseDir)
        # store dummys for several vaiables to make the config-files
        # executable even if they contain those variables that we'll only
        # define later after finding the odb / jobname
        for attr in [
                "jobName", "odbPath", "meshTypeSplitVersion", "meshVersion",
                "initStressVersion", "seqVersion", "matVersion",
                "geoDomVersion", "baseRunVersion"]:
            localsDict[attr] = "<%s>" % attr

        # parse config files only to find "runBaseDir"
        self.parseConfigFile(configFiles, localsDict, printMessage)
        runBaseDir = localsDict.get("runBaseDir")
        del localsDict

        # identify the base folder where we'd expect to find the odb
        if runBaseDir is None:
            # default location: $(projectBaseDir)/2_RUN
            self.runBaseDir = os.path.join(self.projectBaseDir, "2_RUN")
        else:
            # specified alternative location
            # ... normpath cleans trailing slashes
            self.runBaseDir = os.path.normpath(runBaseDir)

        odbDir = os.path.join(self.runBaseDir, self.runVersion)

        # try to locate the odb
        rex_odb = re.compile("([^_]+?)_%s_(.+)\\.odb" % self.runVersion)
        for fn in os.listdir(odbDir):
            res = rex_odb.match(fn)
            if res:
                break

        # if no odb then try to locate the postData file
        if not res:
            rex_postData = re.compile(
                "([^_]+?)_%s_(.+)_postData\\.py" % self.runVersion)
            for fn in os.listdir(odbDir):
                res = rex_postData.match(fn)
                if res:
                    break

        # if no odb then try to locate the postData folder
        if not res:
            rex_postData = re.compile(
                "([^_]+?)_%s_(.+)_postData" % self.runVersion)
            for fn in os.listdir(odbDir):
                res = rex_postData.match(fn)
                if res:
                    break

        # if no odb and no postdata then try to find the main input file
        if not res:
            rex_mainInp = re.compile(
                "([^_]+?)_%s_((?:.+_)*G\\d\\d.*_Q\\d\\d.*)\\.inp"
                % self.runVersion)
            for fn in os.listdir(odbDir):
                res = rex_mainInp.match(fn)
                if res:
                    break

        if res:
            # process job name from filename
            self.projectPrefix = res.group(1)
            self.jobName = "_".join((
                self.projectPrefix, self.runVersion, res.group(2)))
            self.odbPath = os.path.join(odbDir, "%s.odb" % self.jobName)
            otherVersionList = res.group(2).split("_")
        else:
            msg("WARNING: Couldn't find odb, postData or main input file."
                " Trying to determine further parameters from configBE.py"
                " files.")
            self.projectPrefix = os.path.basename(self.projectBaseDir)
            otherVersionList = []
        del res

        # determine other xxxVersions
        rex_meshTypeSplitVersion = re.compile(r"(M\d+)[MHT]v\d+")
        for word in otherVersionList:
            res_meshTypeSplitVersion = rex_meshTypeSplitVersion.match(word)
            if res_meshTypeSplitVersion:
                self.meshTypeSplitVersion = word
                self.meshVersion = res_meshTypeSplitVersion.group(1)
            elif word[0]=="G":
                # simple meshVersion
                self.meshTypeSplitVersion = word
                self.meshVersion = word
            elif word[0]=="S":
                self.initStressVersion = word
            elif word[0]=="Q":
                self.seqVersion = word
            elif word[0]=="M":
                self.matVersion = word
            elif word[0]=="D":
                self.geoDomVersion = word

        # if restart run then determine base run
        if "_" in self.runVersion:
            res_baseRunVersion = re.match(
                r"R(.+)_RST(.+)", self.runVersion)
            if res_baseRunVersion:
                self.baseRunVersion = "R%s" % res_baseRunVersion.group(2)

        # finally parse config files a second time
        self.parseConfigFile(configFiles, self.__dict__, printMessage)


        # add contents of postdata file if possible
        # deduct postData-filename from odb name
        self.postData = getOdbPostData(self.odbPath)


    def parseConfigFile(self, configFiles, localsDict, printMessage=False):
        """Parse the given config files in the order found in the argument
        using localsDict as locals name space.
        Add the directory of the current config file as local variable
        "configDir" to locals.
        """
        for confFile in configFiles:
            localsDict["configDir"] = os.path.dirname(confFile)
            localsDict["os"] = os
            if printMessage:
                msg("Reading config file %s" % confFile)
            execfile(confFile, localsDict)
            del localsDict["configDir"]

    
#------------------------------------------------------------------------------

def newFile(filename, mode="w"):
    """open a new file for writing but make sure to not overwrite an existing
    file with the same name by first renaming the existing.

    Note: May fail if the required permissions are missing.
    """
    if os.path.exists(filename):
        basename, ext = os.path.splitext(filename)
        basename += "_bak"
        bakFilename = basename + ext
        # This next chunk could be used if you want to keep all arbitrarily old
        # versions of this file. Don't overwrite the backup with a new backup.
        # But we don't want to drown in backups...
        # while os.path.exists(bakFilename):
        #     basename = incrementName(basename)
        #     bakFilename = basename + ext
        if os.path.isfile(bakFilename):
            os.remove(bakFilename)
        os.rename(filename, bakFilename)
        msg("WARNING: File %s already existed. Renamed it to %s"
            % (filename, bakFilename))

    if "w" not in mode:
        mode = "w"+mode

    return open(filename, mode)

#------------------------------------------------------------------------------

class CheckFileFinished(object):
    """Use this class to operate on files that may still be modified by other
    parties (like Cavesim, Dropbox) while you already try to access them.
    Makes sure that the file is not empty and that the file size does not
    change anymore.

    The functionality is taken from the cave coupling controller.
    Uses L{bae.log_01} to provide status information.

    Usage:
     >>> from bae.misc_01 import CheckFileFinished
     >>> c = CheckFileFinished(inputPath)
     >>> while True:
     >>>     process(inputPath)
     >>>     if c.finished():
     >>>         break

    Example: Wait until input file apears and has more than zero length. Then
    open the file and process its content. If the file size has changed after
    processing plus some delay then do it again.
     >>> from bae.misc_01 import CheckFileFinished
     >>> inputPath = "OUTPUT/VELOCITY_VTK/CAVESIM_VEL_4.vtk"
     >>> while not os.path.isfile(inputPath):
     >>>     time.sleep(1)
     >>> c = CheckFileFinished(inputPath)
     >>> while True:
     >>>     time.sleep(safetyDelaySeconds)
     >>>     with open(inputPath) as inp:
     >>>         process(inp)
     >>>     if c.finished():
     >>>         break
     >>> del c  # not necessary
    """
    def __init__(self, path, sizeZeroCheckIntervall=1.0, maxWaitPeriod=None):
        """
        @param path: The full path of the file.
        @param sizeZeroCheckIntervall: time in seconds to wait before checking
          again if the file has ztill zero length.
        @param maxWaitPeriod: Wait at most that much seconds for the file size
          to become >0. Specify 0 or None to wait forever.
        """
        self.path = path

        # check that file is larger than zero bytes
        starttime = time.clock()
        ticker = MsgTicker_("File size of %s still zero." % path)
        while os.path.getsize(path)==0:
            time.sleep(sizeZeroCheckIntervall)
            ticker.msg()
            if maxWaitPeriod and (time.clock()-starttime>maxWaitPeriod):
                raise Exception(
                    "ERROR: Waited %s seconds for %s to get non-zero file"
                    " size." % (maxWaitPeriod, path))
        del ticker

        self.lastFileSize = os.path.getsize(path)
        self.testIterations = 0

    def finished(self, maxIterations=20):
        """Return True if the file size of self.path has not changed since
        last call to self.__init__() or self.stillChanging().

        @param maxIterations: Stop after maxIterations and raise Exception
          If set to zero or None then check forever.
        """
        lfs = self.lastFileSize
        self.lastFileSize = os.path.getsize(self.path)
        finished = (lfs == self.lastFileSize)
        if finished:
            msg("Operation on %s has finished without file size having"
                " changed. Assuming success." % self.path)
        else:
            self.testIterations += 1
            if maxIterations and (self.testIterations > maxIterations):
                msg("File size has changed while operating on %s."
                    % self.path)
                raise Exception(
                    "ERROR: Stopping after trying the operation for %d times."
                    % maxIterations)
            else:
                msg("File size has changed while operating on %s. Repeating"
                    " the operation." % self.path)
                
        return finished

#} end of operating system tasks: file transfer, job control


#{Spreadsheets: Data from csv files
#
###############################################################################
### npContainer-Class and helper-functions
###############################################################################

### Exceptions
class FieldRequestError(Exception):
    pass


class DuplicateDataError(Exception):
    pass


### decorators
def _stringInputAsLists(func):
    ''' varNames, frameNames and columns (append var if needed) will be passed
    as lists of strings
    '''
    if epyDocBuild:
        return func
    else:
        def func_wrapper(*args,**kwargs):
            if 'varNames' in kwargs:
                if isinstance(kwargs['varNames'], basestring):
                    kwargs['varNames'] = [kwargs['varNames']]
            if 'frameNames' in kwargs:
                if isinstance(kwargs['frameNames'], basestring):
                    kwargs['frameNames'] = [kwargs['frameNames']]
            if 'columns' in kwargs:
                if isinstance(kwargs['columns'], basestring):
                    kwargs['columns'] = [kwargs['columns']]
            return func(*args,**kwargs)
        return func_wrapper


def _defaultNoneReplacement(func):
    ''' Unspecified varNames and frameNames will be replaced by existend values
    in self.info['varNames'] or self.info['frameNames']
    '''
    if epyDocBuild:
        return func
    else:
        def func_wrapper(*args,**kwargs):
            if ('varNames' not in kwargs) or (kwargs['varNames'] is None):
                kwargs['varNames'] = args[0].info['vars']
            if ('frameNames' not in kwargs) or (kwargs['frameNames'] is None):
                kwargs['frameNames'] = args[0].info['frames']
            return func(*args,**kwargs)
        return func_wrapper


def _defaultNoneWarning(func):
    ''' Raise warning if varNames is not specified
    '''
    if epyDocBuild:
        return func
    else:
        def func_wrapper(*args,**kwargs):
            if ('varNames' not in kwargs) or (kwargs['varNames'] is None):
                msg("WARNING! The order of variables will be the same as "+
                    "found in self.info['var']. Please specify varNames "+
                    "to enshure correct order.")
            return func(*args,**kwargs)
        return func_wrapper


def _defaultReaderWriterArgs(func):
    '''sets reader/writer keyword arguments (e.g. delimiter) if they are not
    given
    '''
    if epyDocBuild:
        return func
    else:
        def func_wrapper(*args,**kwargs):
            if ('delimiter' not in kwargs):
                kwargs['delimiter'] = ','
            return func(*args,**kwargs)
        return func_wrapper


def _itemString(items):
    '''returns nicely formated stringlist
    '''
    string = ', '.join(items)
    strings = string.rsplit(', ',1)
    return ' and '.join(strings)


### helper-functions npContainer ##############################################

def dTypeMapper(dTypeKeys, dType=None):
    """ Creates dtype list from string list

    Example:
    ========
        >>> dTypeMapper(['a', 'b'], dType=bool)
        >>> [('a', bool), ('b', bool)]
        >>> dTypeMapper(['a', 'b'], ['var', 'data'])
        >>> [('a', [('var', 'S8'), ('data', 'f4')]),
        >>> ('b', [('var', 'S8'), ('data', 'f4')])]

    @param dTypeKeys: list of strings or list of lists of strings. Later will
        create an nested dtype

    @param dType: sets dtype for inner dtype level. None, each dtype will be
        set by innerDType() according to is key
    """
    if (not(isinstance(dTypeKeys[0], list))
            and not(isinstance(dTypeKeys[0],np.ndarray))):
        dTypeKeys = [dTypeKeys, ]

    def innerDType(varName):
        'rule to assign dtype for inner dTypeKeys[-1]'
        if varName == 'var' or varName == 'ptId':
            return (varName,'S32')
        elif 'ptIdx' in varName:
            return (varName, 'i8')
        else:
            return (varName, float)

    inner = True
    while not dTypeKeys == []:
        keys = dTypeKeys.pop()
        if inner:
            if dType:
                dT = [(k,dType) for k in keys]
            else:
                dT = map(innerDType,keys)
            inner = False
        else:
            dT = [(key,dT) for key in keys]

    return dT


@_stringInputAsLists
def toArray(sArr, fieldNames=None, unTuple=True):
    '''Converts structured array to a regular numpy-array. By now only 2D-like
    data is implemented.

    Example:
    ========
        >>> theType = [('a',float),('b',bool),('c',float)]
        >>> structArr = np.asarray(map(tuple,np.arange(6).reshape(2,3)),theType)
        >>> structArr
        >>> array([(0.,  True, 2.), (3.,  True, 5.)],
        >>>           dtype=[('a', '<f8'), ('b', '?'), ('c', '<f8')])
        >>> toArray(structArr)
        >>> array([[0.0, True, 2.0],[3.0, True, 5.0]], dtype=object)
        >>> toArray(structArr,['a','c'])
        >>> array([[0., 2.],[3., 5.]])


    @param sArr: structured numpyArray of shape nSets holding nDTypes
        dtypes.

    @param fieldNames: specify fieldNames here if not allready done before
        calling toArray.
    '''
    #!!! varNames and/or frameNames will be converted to lists by
    #!!! decorator @_stringInputAsLists

    if fieldNames is None:
        fieldNames = [fN for fN in sArr.dtype.fields.keys()]

    outDType = dict((name, sArr.dtype.fields[name]) for name in fieldNames)

    outArr = sArr.getfield(np.dtype(outDType))
    if unTuple:
        dTypes = [outDType[key][0] for key in outDType]
        allSameDTypes = np.all([dTypes[0] == dT for dT in dTypes])
        if allSameDTypes:
            #The following code worked fine befor numpy 1.14. It did not create
            #a copy, so finally np.shares_memory(sArr,outArr) was True. Since
            #numpy 1.14 this works not reliable anymore.
            #tmpArr = np.ascontiguousarray(outArr).view(dTypes[0])
            #outArr = tmpArr.reshape(outArr.shape+(len(dTypes),))
            outArr = np.asarray(map(list,outArr), dtype=dTypes[0])
        else:
            outArr = np.asarray(map(list,outArr),dtype=object)
    return outArr


@_stringInputAsLists
def arrToStructArr(arr,varNames=None,frameNames=None,tensorSym=False):
    '''Transforms a regular numpy-array arr of shape
        1. (nPoints X nFrames X nVars) (vector-like) or
        2. (nPoints X nFrames X 3 X 3) (tensor)

        to a strucutured numpy-array with field-labels according to varNames
        and frameNames. A tensor (2.) will be reshaped to a vector-like array
        first - in this case, tensorSym can be set True to include only the
        independend tensor-components.

        @param arr: numpy array (ndim = 3 or 4) holding field-data

        @param varNames: labels for variables in data-set. Must match the
            last/third dim of arr (a, vector-like) or 9/6 (b, tensor). In case,
            a the data names will be set in order of appearence, in case b,
            the names will be interpreted in row-wise order.

        @param frameNames: labels for frames in data-set. Must match the
            second dim of arr

        @param tensorSym: if arr is tensor, only the independend components
            are returned

        @return: structured numpy-array

    '''
    #!!! varNames and/or frameNames will be converted to lists by
    #!!! decorator @_stringInputAsLists

    if (np.ndim(arr) == 4) and (arr.shape[-2:]==(3,3)):
        #rearange to vector-like field
        arr = arr.reshape(arr.shape[0],arr.shape[1],9)
        if tensorSym:
            indIdx = [0,1,2, 4,5, 8]
            arr = arr[:,:,indIdx]

    if not np.ndim(arr)==3:
        raise ValueError('Array has to be shape of Npoints X Nframes X Nvars'+
                         'or Npoints X Nframes X 3 X 3')

    if frameNames is not None:
        if not len(frameNames)== arr.shape[1]:
            raise IndexError('Length of frameNames-List must be the same as'+
                             ' second dim of array (arr.shape[1])')
    else:
        frameNames = ['F%00d' % d for d in range(arr.shape[1])]

    if varNames is not None:
        if not len(varNames)== arr.shape[2]:
            raise IndexError('Length of varsNames-List must be the same as'+
                             ' third dim of array (arr.shape[2])')
    else:
        varNames = ['var%d' % d for d in range(arr.shape[2])]

    dType = dTypeMapper([varNames,frameNames])
    #haven't found a way yet to convert an nD-array (nD>2) directly to
    #a structured array by using something like:
    #               structArray = nDArray.astype(dType)
    # --> create new array and assign data in a loop :-(
    structArr = np.zeros(arr.shape[0],dType)
    for ii,var in enumerate(varNames):
        for jj,frame in enumerate(frameNames):
            structArr[var][frame] = arr[:,jj,ii]

    return structArr


def readCsvHeader(inp, delimiter=','):
    '''Read headerline from CSV-file

    @param inp: csv-file

    @param delimiter: delimiter (default ',')
    '''
#    with open(fileName, 'r') as f:
#        headerStrings = f.readline().strip().split(delimiter)

    if isinstance(inp, basestring):
        fid = open(inp, 'r')
    else:
        fid = inp
    headerStrings = fid.readline().strip().split(delimiter)

    fid.seek(0)
    return headerStrings


@_defaultReaderWriterArgs
@_stringInputAsLists
def npReadCsv(fileName, columns=None, rowCondition=None, **readerArgs):
    '''Reads a regular csv file to structured numpy array. Each requested
    column can be indexed by its header string.

    @param fileName: Name of csvfile. Headerstring is used for indexing.

    @param readerArgs: depends on csvreader: eg. delimiter=','
        see the Documentation of  U{pandas read_csv <https://
        pandas.pydata.org/pandas-docs/stable/generated/pandas.read_csv.html>}
        or/and U{numpys genfromtxt<https://docs.scipy.org/
        doc/numpy-1.13.0/reference/generated/numpy.genfromtxt.html>}

    @param columns: list of columns to read (eg. ['x','y','var','F012'] )
        If columns==None --> read all columns in file

    @param rowCondition: This is only supported in pandas. For numpy
        implementaion see: U{https://stackoverflow.com/questions/14645789/
        numpy-reading-file-with-filtering-lines-on-the-fly}
        Pass a function that specifies under which condition each row will be
        stacked to output data.

        Example:
         >>> #return only rows where column 'var' equals 'var1'
         >>> rCondi = lambda data: data['var']=='var1'
         >>> readCsv(csvName, rowCondition=rCondi)

    @return: structured numpy array where each requested row can be indexed by
        its header string.
    '''
    #!!! ReaderArgs will be set automatically by decorator
    #!!! @_defaultReaderWriterArgs
    #!!! varNames and/or frameNames will be converted to lists by
    #!!! decorator @_stringInputAsLists

    #hasPandas = False #testSwitch

    ### 1) read and analyse header of csv-file
    cNames = readCsvHeader(fileName,delimiter=',')
    if columns is None:
        columns = cNames  # read all

    cols = [col if col in cNames else None for col in columns]
    if not all(cols):
        raise KeyError("Can't find row name(s) '%s' in header of %s "
                       % (','.join([d[1] for d in zip(cols,columns)
                                    if d[0] is None]), fileName))
    ### 2) create dtype for columns
    dType = dTypeMapper(cols)

    if not hasPandas:
        ### 3a) reading data to structured numpy array via numpy.genfromtxt
        if rowCondition:
            msg('Conditional reading by row is only supported for pandas')

        #msg('Start reading %s (using numpy) ...'%os.path.basename(fileName))

        data = np.genfromtxt(fileName,usecols=cols,
                             names=True, dtype=dType, **readerArgs)
        #msg('reading finished')
    else:
        ### 3b) use pandas an transform to structured numpy array
        if rowCondition is None:
            #msg('Start reading %s (using pandas) ...'%\
            #    os.path.basename(fileName))
            indata = pd.read_csv(fileName,
                                 usecols=cols,
                                 **readerArgs)
        else:
            #msg('Start reading %s (pandas, row-conditional) ...'%\
            #    os.path.basename(fileName))
            def valid(chunks):
                for i,chunk in enumerate(chunks):
                    mask = rowCondition(chunk)
                    yield chunk.loc[mask][cols]

            chunksize = 10 ** 5
            #to resolve the 'row-dependencies' caused by rowCondition, the
            #raised 'KeyErrors' will be parsed and the missing keys are
            #appended to usecols while they are not included yet. If they are
            #allready in usecols, than the fileheader does not contain this
            #key and the reading procedure will be canceled
            dependCols = []
            cnt = 0
            while True and cnt<100:
                cnt=cnt+1
                try:
                    chunks = pd.read_csv(fileName,
                                         chunksize=chunksize,
                                         usecols=cols + dependCols,
                                         **readerArgs)
                    indata = pd.concat(valid(chunks))
                    break  # finished
                except KeyError as e:
                    # recursively resolve column 'dependencies' of rowCondition
                    colName = re.search(r"'([A-Za-z]*)'",
                                        str(e)).group(0).replace("'","")
                    if colName not in dependCols:
                        dependCols.append(colName)
                    else:
                        raise KeyError(e)
                except ValueError as e:
                    if str(e) == "Usecols do not match names.":
                        raise ValueError((
                            "Check if the following columns can be found in " +
                            "header of %s:\n" +
                            "requested by colName: '%s'\n" +
                            "requested by rowCondition: '%s'\n" )
                            % (fileName, "', '".join(cols) +
                               "', '".join(dependCols)))
                    else:
                        raise ValueError(e)
        #msg('reading finished')
        # transform to structured numpy array
        # if all dTypes are float use as_matrix + ravel.view:
        dTs = [dt[1] for dt in dType]
        if all([dTs[0]==dt for dt in dTs]):
            #!!! as_matrix is deprecated since pandas 0.23.0
            #data = indata.as_matrix().ravel().view(dtype=dTs)
            data = indata.values().ravel().view(dtype=dTs)
        # for mixed dtypes preallocate output data:
        else:
            if indata.values.nbytes*1E-6 > 100:
                msg('The dataset extends the size of 100Mb. To speed up' +
                    ' read columns with different dtypes (e.g. float and ' +
                    'string) separately .')
            data = np.zeros(indata.shape[0], dtype=dType)
            for name in cols:
                data[name] = indata[name]
    return data


########## class npContainer ##################################################

class npContainer(object):
    '''Reader-, writer- and "management"-class for regular csv-files according
    to BE-standard::
        x, y, z, var, frmId1, frmId2, ...

    Point and field-data will be stored independently in self.points and
    self.data respectively. self.data is a structured numpy array with the
    same shape and order as self.points. The field-data can be indexed by
    (multible) variable and frame names:
     >>> self.data[['var1','var2',...]][['frame1','frame2',...]]

    self.points contains only a unique and sorted set of points. Because of the
    fact that the data is not necessarily is defined on all points, invalid
    subsets will be represented by np.nan. The mapping to valid points can be
    done by the boolean-indices in self.validPoints:
     >>> validVarData = self.data['var'][self.validPoints['var']]
    '''
    #register ErrorClass
    FieldRequestError = FieldRequestError
    DuplicateDataError = DuplicateDataError
    #Usage:
    #uF = npContainer
    #try:
    #   uF.readFromCsv(...)
    #except npContainer.FieldRequestError:
    #   print('Ups')

    def __init__(self):
        if np is None:
            raise NotImplementedError('This class requires numpy')
        if not hasPandas:
            msg('To speed up csv reading, please install pandas')

        self.data       = None
        self.points     = None
        self.info       = {}
        self.validPoints= None
        self.pointTol   = 1E-2
        self.defaultHead= ['x','y','z','var']

    def copy(self):
        import copy
        return copy.deepcopy(self)

    @_stringInputAsLists
    def fromArray(self,arr,xyz,varNames=None,frameNames=None):
        '''
        Creates a data-set from regular numpy array arr of shape
        (nPoints X nFrames X nVars) (vector-like) or
        (nPoints X nFrames X 3 X 3) (tensor) and its coordinate-array xyz
        (nPoints X 3). A tensor (b) will be reshaped to a vector-like array
        first - in this case, tensorSym can be specified to include only the
        independend tensor-components.

        Example:
        ========
         >>> sLabels = ['S11','S22','S33','S23','S13','S12']
         >>> fLabels = ['F000','F001']
         >>> points  = np.random.randn(100,3)
         >>> array   = np.random.randn(len(points),len(fLabels),len(sLabels))
         >>> uF = npContainer().fromArray(array,points,varLabels=sLabels,
         >>>                       frameLabels=fLabels)

        @param arr: numpy array (ndim = 3 or 4) holding field-data

        @param xyz: numpy array holding coordinates

        @param varNames: labels for variables in data-set. Must match the
            last/third dim of arr (a, vector-like) or 9/6 (b, tensor). In
            case (a), the data names will be set in order of appearence, in
            case (b), the names will be interpreted in row-wise order.
            See L{npContainer.reshapeAsTensor}

        @param frameNames: labels for frames in data-set. Must match the
            second dim of arr

        @return: self, so you can use:

            >>> dataSet = npContainer().fromArray(myArr,myXYZ)

        '''

        #!!! varNames and/or frameNames will be converted to lists by
        #!!! decorator @_stringInputAsLists

        arr = np.asarray(arr)
        xyz = np.asarray(xyz)

        #get type of arr and check dimensions
        dimErr = False
        if arr.ndim == 3:    # vector-like- has to be nPt X nFr X nVars
            isTensorType = False
            nPt,nFr,nVar = arr.shape
        elif arr.ndim == 4:  # tensor - has to be nPt X nFr X 3 X 3
            isTensorType = True
            nPt,nFr,nA,nB = arr.shape
            if not ((nA==3) and (nB==3)):
                dimErr = True
            nVar = nA*nB
        else:
            dimErr = True

        if dimErr:
            raise ValueError('DataArray has to be shape of '+
                             'nPoints X nFrames X nVars or '+
                             'nPoints X nFrames X 3 X 3')

        #check dimension of xyz
        ptErr = False
        if (np.obj2sctype(xyz) == np.void):  # is structured array
            if not (sorted(xyz.dtype.fields.keys()) == ['x','y','z']):
                ptErr is True
        else:
            if not xyz.shape == (nPt,3):
                ptErr is True
            else:
                xyz = np.array(map(tuple,xyz), dTypeMapper(['x','y','z'],float))
        if ptErr:
            raise ValueError(
                "PointArray has to be an array of shape"
                " (nPt,3) or a sturctured numpy array with"
                " dtype-keys ['x','y','z'] and shape (nPt,)."
                " The DataArray you passed has %d points"
                " (see DataArray.shape[0])" % nPt)

        #check consistency of varNames
        if varNames is None:
            varNames = ['var%d' % d for d in range(nVar)]

        tensorSym = False
        if isTensorType:
            if len(varNames) == 6:
                tensorSym = True  # assume symmetric tensor
            elif len(varNames) == 9:
                tensorSym = False
            else:
                raise IndexError(
                    'Number of varNames (%d) has to be 6 or 9 for tensor'
                    ' arrays' % len(varNames))
        elif not (len(varNames)==nVar):
            raise IndexError(
                'Number of varNames (%d) does not match the third dim of array'
                ' (%d)' % (len(varNames), nVar))

        #check consistency of frameNames
        if frameNames is None:
            frameNames = ['F%03d' % d for d in range(nFr)]
        elif not (len(frameNames)==nFr):
            raise IndexError(
                'Number of varNames (%d) does not match the third dim of array'
                ' (%d)' % (len(frameNames), nFr))

        #check for duplicate points
        xyz,ptIdx,ptCnt = self.uniquePtsTol(xyz)

        if np.any(ptCnt>1):
            raise DuplicateDataError(
                'There are %d duplicate points in xyz.' % np.sum(ptCnt>1))

        self.data = np.zeros(len(ptIdx), dTypeMapper([varNames,frameNames]))
        self.data[ptIdx] = arrToStructArr(arr,varNames=varNames,
                                          frameNames=frameNames,
                                          tensorSym=tensorSym)
        self.points = xyz
        self.validPoints = np.ones(nPt, dTypeMapper(varNames,bool))
        self.info['vars'] = varNames
        self.info['frames'] = frameNames

        return self

    def fromPVFContainer(self,pvf):
        '''Creates NP-Container from a PVF-Container
        (see L{bae.misc_01.PVFContainer}).

        @param pvf: PVF-Container

        @return: self, so you can do:
            >>> myNpC = NpContainer().fromPVFContainer(pvf)

        '''

        if not hasattr(pvf, 'ptInfo'):
            raise ValueError("The PVF-Container has no 'ptInfo'-attribute. "+
                             "This is necessary to convert to NP-Container.")

        ptKeys = pvf.ptData.keys()
        nPts = len(ptKeys)

        varKeys = []
        for ptVals in pvf.ptData.values():
            varKeys.append(ptVals.keys())
        varNames = np.unique(varKeys)
        frameNames = pvf.frameIds

        msg("PVF-Container contains  points:%d,  frames:%d, variables:%d"
            % (nPts,len(frameNames),len(varNames)))

        data        = np.empty(nPts, dTypeMapper([varNames,frameNames]))
        data[:]     = np.nan
        validPoints = np.zeros(nPts, dTypeMapper(varNames,dType=bool))
        xyz         = []

        for ptIdx,(ptKey,ptData) in enumerate(pvf.ptData.iteritems()):
            xyz.append(pvf.ptInfo[ptKey].ptCoords)
            for varKey,varData in ptData.iteritems():
                data[ptIdx][varKey] = np.array(varData)
                validPoints[ptIdx][varKey] = True

        xyz = np.array(map(tuple,xyz),dtype=dTypeMapper(['x','y','z']))
        xyz,ptIdx,ptCnt = self.uniquePtsTol(xyz)

        if len(ptIdx) < nPts:
            raise DuplicateDataError(
                "The PVF-Container seems to have %d duplicated points."
                % (nPts-len(ptIdx)))

        #copy or a new allocation (np.empty(len,dtype)) is necessary here!
        #Otherwise the view will be overwriten by itself when it is reordered
        self.data        = data.copy()
        self.validPoints = validPoints.copy()
        self.data[ptIdx] = data
        self.validPoints[ptIdx] = validPoints
        self.points = xyz
        self.info['vars'] = varNames
        self.info['frames'] = frameNames

        return self

    @_defaultReaderWriterArgs
    @_stringInputAsLists
    def readFromCsv(self,fileName, varNames=None, frameNames=None,
                    **readerArgs):
        '''Reads data from  PointVariableFrame-ordered csv

        Example:
        ========
        >>> sLabels = ['S11','S22','S33','S23','S13','S12']
        >>> uF = npContainer().readFromCsv(csvName,varLabels=slabels,
        >>>                         frameLabels=['F000'],delimiter=';')

        @param fileName: full path to csv-file

        @param varNames: list of variable names to read (None: read all)

        @param frameNames: list of frame names to read (None: read all)

        @param readerArgs: additional reader arguments. See L{npReadCsv}

        '''
        #!!! ReaderArgs will be set automatically by decorator
        #!!! @_defaultReaderWriterArgs
        #!!! varNames and/or frameNames will be converted to lists by
        #!!! decorator @_stringInputAsLists

        if hasPandas and not (varNames is None):
            rowCondition = lambda data: np.isin(data['var'],varNames)
        else:
            rowCondition = None

        if frameNames is not None:
            readLabels = self.defaultHead+[f for f in frameNames
                                           if f not in self.defaultHead]
        else:
            readLabels = None

        data = npReadCsv(fileName, columns=readLabels,
                         rowCondition=rowCondition,**readerArgs)
        if not data.size:
            raise FieldRequestError("Dataset is empty. Maybe the "
                                    "rowContition isn't defined correctly.")

        varNamesFound = np.unique(data['var']).tolist()
        if varNames is None:
            msg('Found variables %s' % _itemString(varNamesFound))
        else:
            foundAll = all([d in varNamesFound for d in varNames])
            if not foundAll:
                raise FieldRequestError("Can't find the requested variables "+
                                        "%s" % _itemString(varNames))
        varNames = varNamesFound

        frameNames = list(data.dtype.names)
        for noFrame in self.defaultHead:
            frameNames.remove(noFrame)

        #create lookuptable (ptIdx) for pointlist
        points = data[['x','y','z']]

        #get unique pointset and the index lookup table of data
        #Keep in mind that the original order of points as they appear in csv
        #migth be disrupted. It can be restored by using the 'return_index'.
        points,ptIdx = np.unique(points,return_inverse=True,axis=0)

        dType = dTypeMapper([varNames,frameNames])

        self.data = np.empty(len(points),dType)
        self.data[:] = np.nan

        self.validPoints = np.zeros(len(points),
                                    dTypeMapper(varNames,dType=bool))

        for var in varNames:
            varIdx = data['var']==var
            dat = data[varIdx]
            self.data[var][ptIdx[varIdx]] = dat[frameNames]
            self.validPoints[var][ptIdx[varIdx]] = True

        self.points = points
        self.info['frames'] = frameNames
        self.info['vars'] = varNames

        return self

    def pointArray(self):
        '''returns points as regular numpy array '''
        return toArray(self.points)

    @_defaultNoneReplacement
    @_stringInputAsLists
    def reshapeAsVector(self,varNames=None,frameNames=None,allPoints=False):
        ''' Reshape data from stuctured numpy-array to a regular numpy-array.

        Example:
        ========
        >>> myVars = ['var1','var2',...]
        >>> uF = NPC().readFromCsv(varData.csv,varNames=myVars,order='V')
        >>> myVector,ptIdx = uF.reshapeAsVector(varNames=myVars[:4])
        >>> myPoints = uF.points[ptIdx]

        @param varNames: list of related tensor arrays. If None the varialbles
            will appeare in vector according to their position in
            self.info['varNames']

        @param frameNames: List of frames to return

        @param allPoints: If False, only points with non-NaN-values will be
            returned. This will fail if not all requested variables are valid
            on the same points. Try allPoints = True, than all points,
            including those with NaN-values will be returned

        @return: numpy-array ( shape = (nPoints,nFrames,len(varLabels)) )
                 pointIndex  ( shape = (npoints,) )

        '''
        #!!! Note that None-inputvalues for varNames and frameNames will be
        #!!! replaced with 'all existing ones' by decorator
        #!!! varNames and/or frameNames will be converted to lists by
        #!!! decorator @_stringInputAsLists
        nFrames = len(frameNames)
        ptIdx = toArray(self.validPoints,varNames)
        if not allPoints:
            if not np.array_equiv(ptIdx[:,0].reshape(-1,1),ptIdx):
                raise FieldRequestError(
                    "The variables %s are defined on different points. Use"
                    " allPoints=True instead." % _itemString(varNames))
            ptIdx = ptIdx[:,0]
        else:
            ptIdx = np.ones(ptIdx.shape[0],dtype=bool)
        nPts = ptIdx.sum()

        vector = np.zeros((nPts,nFrames,len(varNames)))
        for ii,v in enumerate(varNames):
            for kk,frame in enumerate(frameNames):
                vector[:,kk,ii] = self.data[v][frame][ptIdx]

        return vector,ptIdx

    @_defaultNoneWarning
    @_defaultNoneReplacement
    @_stringInputAsLists
    def reshapeAsTensor(self, varNames=None, frameNames=None,
                        allPoints=False, order='R'):
        '''Reshape data from stuctured numpy array to a regular numpy-array

        Example:
        ========
        >>> myVars = ['var1','var2',...]
        >>> uF = NPC().readFromCsv(varData.csv,varNames=myVars,order='V')
        >>> myTensor,ptIdx = uF.reshapeAsTensor(varNames=myVars[:6])
        >>> myPoints = uF.points[ptIdx]

        @param varNames: list of related tensor arrays.
            If None (NOT RECOMMENDED) the stored dataset has to contain exactly
            6 or 9 variables, which will be mapped according to their position
            in self.info['varNames']

        @param frameNames: List of frames to return

        @param allPoints: If False, only points with non-NaN-values will be
            returned. This will fail if not all requested variables are valid
            on the same points. Try allPoints = True, than all points,
            including those with NaN-values will be returned

        @param order:
            The order of varLabels determines its position in 3x3-tensor
             1. symetric ordering
              - 'RS' [0,1,2,3,4,5] --> [[0 1 2] [1 3 4] [2 4 5]]
              - 'CS' [0,1,2,3,4,5] --> [[0 1 3] [1 2 4] [3 4 5]]
              - 'V'  [0,1,2,3,4,5] --> [[0 5 4] [5 1 3] [4 3 2]]
             2. non-symetric ordering:
              - 'R'  [0,1,2,3,4,5,6,7,8] --> [[0 1 2] [3 4 5] [6 7 8]]
              - 'C'  [0,1,2,3,4,5,6,7,8] --> [[0 3 6] [1 4 7] [2 5 8]]

        @return: numpy-array ( shape = (nPoints,nFrames,3,3)) )
                 pointIndex  ( shape = (npoints,) )

        '''
        #!!! Note that None-inputvalues for varNames and frameNames will be
        #!!! replaced with 'all existing ones' by decorator
        #!!! varNames and/or frameNames will be converted to lists by
        #!!! decorator @_stringInputAsLists

        nFrames = len(frameNames)

        ptIdx = toArray(self.validPoints,varNames)

        if not allPoints:
            if not np.array_equiv(ptIdx[:,0].reshape(-1,1),ptIdx):
                raise FieldRequestError(
                    "The variables %s are defined on different points. Use"
                    " allPoints=True instead." % _itemString(varNames))
            ptIdx = ptIdx[:,0]
        else:
            ptIdx = np.ones(ptIdx.shape[0],dtype=bool)
        nPts = ptIdx.sum()

        order = order.upper()
        if len(varNames) == 6:
            if order in ['R','RS','C','CS']:
                order = 'RS'
            elif order in ['C','CS']:
                order = 'CS'
            elif not order == 'V':
                e = "%s  no valid type for 'order' of nonsym tensor" % order
                raise ValueError(e)
        elif (len(varNames) == 9):
            if not (order in ['R','C']):
                e = "%s is no valid type for 'order' of nonsym tensor" % order
                raise ValueError(e)

        #rearrange varLabels with respect to their position in tensor
        if order == 'RS':                   # rowwise symmetric
            srcIdx = [0,1,2,1,3,4,2,4,5]
        elif order == 'CS':                 # columnwise symmetric
            srcIdx = [0,1,3,1,2,4,3,4,5]
        elif order == 'V':                  # Voigt symmetric
            srcIdx = [0,5,4,5,1,3,4,3,2]
        elif order == 'R':                  # rowwise not symmetric
            srcIdx = [0,1,2,3,4,5,6,7,8]
        elif order == 'C':                  # columnwise not symmetric
            srcIdx = [0,3,6,1,4,7,2,5,8]
        varNames = np.array(varNames)[srcIdx]  # ordered varNames

        #stacking data to nPointsXnFramesX3X3 numpy array
        arrayTensor = np.zeros((nPts,nFrames,3,3))
        i = 0
        for iR in range(3):
            for iC in range(3):
                Da = self.data[varNames[i]][frameNames]  # just creates a view
                for kk,frame in enumerate(frameNames):
                    arrayTensor[:,kk,iR,iC] = Da[frame][ptIdx]
                i = i+1

        return arrayTensor,ptIdx

    def uniquePtsTol(self,points,returnRounded=False):
        """Rounds pointData to self.pointTol before calculating unique points

        @param points: points as structured numpy-array with fields 'x','y'
            and 'z'

        @param returnRounded: return rounded unique pointset

        @note: it's rumoring that np.unique has a tolerance-method implemented
            so e.g. 0.5000000001 and 0.499999999 are the same, but this feature
            is not well documented. It seems that an additional tol-keywordArg
            will come with the next numpy-versions. Till than, the following
            code might not be 100% bullet-proof

        @return: uPoints uniquePoints
                 ptIdx inverseIndex (points = uniquePoints(inverseIndex))
                 ptCnt count unique points
        """
        rPoints = np.zeros(points.shape,points.dtype)
        for comp in ['x','y','z']:
            rPoints[comp] = np.floor(points[comp]/self.pointTol)
        uPoints,firstIdx,ptIdx,ptCnt = np.unique(rPoints,
                                                 return_index=True,
                                                 return_inverse=True,
                                                 return_counts=True,
                                                 axis=0)
        if returnRounded:
            for comp in ['x','y','z']:
                rPoints[comp] =rPoints[comp]*self.pointTol
            return uPoints,ptIdx,ptCnt
        else:
            return points[firstIdx],ptIdx,ptCnt

    def getPointsIndex(self, points):
        '''Returns index for a subSet of points

        @param points: list or numpy-array of  shape nPt X 3 or structured
            numpy-array of shape nPt and fields 'x','y','z'

        @return: indices for requested points in dataset ordered as requested
            by points

        '''

        # check dimension of points
        ptErr = False
        if isinstance(points,np.ndarray) and \
           (np.obj2sctype(points) == np.void):
            # is structured array
            if not (sorted(points.dtype.fields.keys()) == ['x','y','z']):
                ptErr = True
        else:
            points = np.atleast_2d(points)
            if not points.shape[1] == 3:
                ptErr = True
            else:
                points = np.array(map(tuple,points),
                                  dTypeMapper(['x','y','z'],float))
        if ptErr:
            raise ValueError("PointArray has to be an array of shape"
                             " (nPt,3) or a sturctured numpy array with"
                             " dtype-keys ['x','y','z'] and shape (nPt,).")

        if not points.size:
            raise ValueError("Requested PointArray is empty")

        nOldPts = len(self.points)
        bothPts = np.hstack([self.points,points])
        uPoints,ptIdx,ptCnt = self.uniquePtsTol(bothPts)
        #oldPtsIdx = ptIdx[:nOldPts]
        newPtsIdx = ptIdx[nOldPts:]
        if np.any(newPtsIdx > nOldPts-1):
            nNotFound = np.sum(newPtsIdx > nOldPts-1)
            raise IndexError(
                "%d of the requested points can't be found in dataSet within"
                " the given point tolerance (%.3f)"
                % (nNotFound, self.pointTol))

        if np.any(ptCnt[newPtsIdx] > 2):
            nDupFound = np.sum(ptCnt[newPtsIdx]>2)
            msg("WARNING! %d of the requested points were found multible"
                " times in dataSet. Reduce point tolerance (%.3f)."
                % (nDupFound, self.pointTol))

        return newPtsIdx

    @_defaultNoneReplacement
    @_stringInputAsLists
    def getPointSubset(self,points,varNames=None,frameNames=None):
        '''Returns a subSet of points

       Example:
       ========
        >>> arr,pts = fullData.getPointSubset(subPoints).reshapeAsArray()


        @param points: list or numpy-array of  shape nPt X 3 or structured
            numpy-array of shape nPt and fields 'x','y','z'

        @param varNames: List of variable names to be returned

        @param frameNames: List of frame names to be returned

        @return: npContainer holding subset of requested points

        '''
        #!!! Note that None-inputvalues for varNames and frameNames will be
        #!!! replaced with 'all existing ones' by decorator
        #!!! varNames and/or frameNames will be converted to lists by
        #!!! decorator @_stringInputAsLists

        #get pointIndex of requested points
        ptIdx = self.getPointsIndex(points)
        ptIdxSort = ptIdx[np.argsort(self.points[ptIdx])]
        #create npContainer for output
        outData = npContainer()
        outDType               = dTypeMapper([varNames,frameNames])
        outData.data           = np.zeros(len(ptIdxSort),outDType)
        for var in varNames:
            outData.data[var]  = self.data[ptIdxSort][var]
        outData.points         = self.points[ptIdxSort]
        outData.validPoints    = self.validPoints[ptIdxSort][varNames]
        outData.info['vars']   = varNames
        outData.info['frames'] = frameNames
        outData.pointTol       = self.pointTol
        return outData

    def appendFrames(self,newData,xyz=None,frameNames=None):
        '''Append frames to self.data

        @param newData: structured numpy-array or npContainer holding the data
            for frames. Must have same points and frames (dtypes) as
            self.data.

        @param xyz: List or structured numpy-array holding point coordinates
            Compulsory if newData is not type of NpContainer

        @param frameNames: List of frame names from newData to append

        '''
        if isinstance(newData,npContainer) and any([xyz,frameNames]):
            msg('Ignoring additional points or frameNames you passed.')
        elif isinstance(newData,np.ndarray) and isinstance(xyz,np.ndarray):
            try:
                newData =npContainer().fromArray(
                    newData, xyz,
                    frameNames=frameNames,
                    varNames=self.info['vars'])
            except:
                msg("While appending new data to npContainer a temporary"
                    " npContainer holding the new data is created."
                    " This creation raises the following Error:")
                raise

        #check that new dataset is present on a subset of existing one
        bothPts = np.hstack([self.points,newData.points])
        xyz,ptIdx,ptCnt = self.uniquePtsTol(bothPts)

        if not np.all(ptCnt == 2):
            raise IndexError('Points in existent and new dataSet differ.')

        fNamesA, fNamesB = self.info['frames'], newData.info['frames']
        frameNames = fNamesA + fNamesB
        dType = dTypeMapper([self.info['vars'],frameNames])
        data = np.zeros(len(self.data)).astype(dType)
        data[:] = np.nan
        for var in self.info['vars']:
            for frame in fNamesA:
                data[var][frame] = self.data[var][frame]
            for frame in fNamesB:
                data[var][frame] = newData.data[var][frame]

        self.data = data
        self.info['frames'] = frameNames

    def appendVars(self,newData,xyz=None,varNames=None):
        '''Append variables to self.data

        @param newData: structured numpy-array or npContainer holding the data
            for variables. Must have same points and frames (dtypes) as
            self.data.

        @param xyz: List or structured numpy-array holding point coordinates
            Compulsory if newData is not type of NpContainer

        @param varNames: List of variable names from newData to append

        '''
        if isinstance(newData,npContainer) and any([xyz,varNames]):
            msg('Ignoring additional points or varNames you passed.')
        elif isinstance(newData,np.ndarray) and isinstance(xyz,np.ndarray):
            try:
                newData =npContainer().fromArray(
                    newData,xyz,
                    frameNames=self.info['frames'],
                    varNames=varNames)
            except:
                msg("While appending new data to npContainer a temporary"
                    " npContainer holding the new data is created."
                    " This creation raises the following Error:")
                raise

        #check that new dataset is present on a subset of existing one
        bothPts = np.hstack([self.points,newData.points])
        xyz,ptIdx,ptCnt = self.uniquePtsTol(bothPts)

        ptIdxA = ptIdx[:len(self.points)]
        ptIdxB = ptIdx[len(self.points):]
        if any(ptIdx > len(self.points)-1):
            raise IndexError(
                'There are points which are not included in existing dataSet.')

        #merge dataset
        #1. create empty arrays
        varNames = self.info['vars'] + newData.info['vars']
        dType = dTypeMapper([varNames,self.info['frames']])
        data = np.zeros(len(xyz),dType)
        data[:] = np.nan
        valid = np.zeros(len(xyz), dTypeMapper(varNames,bool))
        #2. looping existent and new dataset
        for vN in self.info['vars']:
            data[vN][ptIdxA] = self.data[vN]
            valid[vN][ptIdxA] = self.validPoints[vN]
        for vN in newData.info['vars']:
            data[vN][ptIdxB] = newData.data[vN]
            valid[vN][ptIdxB] = True
        #3. overwrite existent with new
        self.data = data
        self.points = xyz
        self.validPoints = valid
        self.info['vars'] = varNames

    def appendPoints(self,newData,xyz=None):
        '''Append new point data to self.data

        @param newData: structured numpy-array or npContainer holding the data
            on new points. Must have same variables and frames (dtypes) as
            self.data.

        @param xyz: List or structured numpy-array holding point coordinates
            Compulsory if newData is not type of NpContainer

        '''
        if isinstance(newData,npContainer) and not(xyz is None):
            msg('Ignoring xyz-Data you passed.')
        elif isinstance(newData,np.ndarray) and isinstance(xyz,np.ndarray):
            try:
                newData =npContainer().fromArray(
                    newData,xyz,
                    frameNames=self.info['frames'],
                    varNames=self.info['vars'])
            except:
                msg("While appending new data to npContainer a temporary"
                    " npContainer holding the new data is created."
                    " This creation raises the following Error:")
                raise

        bothPts = np.hstack([self.points, newData.points])
        points,ptIdx,ptCnt = self.uniquePtsTol(bothPts)
        if np.any(ptCnt>1):
            raise DuplicateDataError(
                '%d of the points to append allready exist in dataSet'
                % np.sum(ptCnt>1))
        try:
            data = np.hstack([self.data,newData.data])
            validPoints = np.hstack([self.validPoints,newData.validPoints])
        except TypeError as e:
            if 'invalid type promotion' in str(e):
                raise TypeError(
                    "Existing and new dataSet seem to have different dTypes."
                    " existing.dtype=%s, new.dtype=%s"
                    % (str(self.data.dtype), str(self.newData.dtype)))
            else:
                raise TypeError(e)

        self.points = points
        self.data   = np.zeros(len(data),data.dtype)
        self.data[ptIdx] = data
        self.validPoints = np.zeros(len(validPoints),validPoints.dtype)
        self.validPoints[ptIdx] = validPoints

    @_defaultReaderWriterArgs
    @_defaultNoneReplacement
    @_stringInputAsLists
    def saveAsCsv(self,fileName,varNames=None,frameNames=None,
                  noHeader=False, mode='fragile',**writerArgs):
        '''Save dataset in following format::
            x, y, z, var, frmId1, frmId2, ...
        to csv file

        @param fileName: filename of csv-file

        @param varNames: specifies variables to be saved

        @param frameNames: specifies frames to be saved

        @param noHeader: True to skip creation of header

        @param mode:  'fragile'--> raise Error if file exists
                      'force'  --> overwrite file if exists
                      'append' --> append to existent file

        @param writerArgs: additional keyword-arguments for used
            csvWriter-method. See L{npReadCsv}
        '''
        #!!! Note that None-inputvalues for varNames and frameNames will be
        #!!! replaced with 'all existing ones' by decorator
        #!!! ReaderArgs will be set automatically by decorator
        #!!! @_defaultReaderWriterArgs
        #!!! varNames and/or frameNames will be converted to lists by
        #!!! decorator @_stringInputAsLists

        if self.data is None:
            raise ValueError('Nothing to save yet. Load data first.')

        if os.path.isfile(fileName):
            if mode.lower() == 'fragile':
                raise IOError('File %s already exists' % fileName)
            elif mode.lower() == 'force':
                msg("Overwriting file %s" % fileName)
                os.remove(fileName)
            elif mode.lower() in ['a','append']:
                msg("Append to file %s" % fileName)

        msg('Writing variables %s to file %s'
            % (_itemString(varNames), fileName))
        fid = open(fileName,'a')
        if not noHeader:
            fid.write(','.join(['x','y','z','var']+frameNames)+'\n')
        for varName in varNames:
            outXYZ  = toArray(self.points[self.validPoints[varName]])
            outData = self.data[varName][frameNames][self.validPoints[varName]]
            outData = toArray(outData)
            outVar  = np.repeat(varName,len(outXYZ))
            #appending data to csv-file

            np.savetxt(fid,np.c_[outXYZ,outVar,outData],fmt="%s",
                       **writerArgs)
        fid.close()
        msg('... done')

    @_defaultReaderWriterArgs
    @_defaultNoneReplacement
    @_stringInputAsLists
    def appendToCsv(self,fileName,varNames=None,frameNames=None,
                    noHeader=False, checkForDuplicates=True,**writerArgs):
        """Appends the dataset in following format::
            x, y, z, var, frmId1, frmId2, ...
        to an existent csv file.

        @param fileName: filename of csv-file

        @param varNames: specifies variables to be saved

        @param frameNames: specifies frames to be saved

        @param noHeader: True to skip creation of header

        @param  checkForDuplicates: True to read content of csv first and
            check if any variable to be stored allready exists on point

        @param writerArgs: additional keyword-arguments for used
            csvWriter-method. See L{npReadCsv}.

        """
        #!!! Note that None-inputvalues for varNames and frameNames will be
        #!!! replaced with 'all existing ones' by decorator
        #!!! ReaderArgs will be set automatically by decorator
        #!!! @_defaultReaderWriterArgs
        #!!! varNames and/or frameNames will be converted to lists by
        #!!! decorator @_stringInputAsLists

        # check for duplicates (same variable on same points)
        if checkForDuplicates:
            for varName in varNames:
                #reading each fieldname is expensive for large csvFiles
                #(>~50Mb). In this case it might be more reasonable to read
                #the whole file at once
                if hasPandas:
                    rowC = lambda data: data['var'] == varName
                else:
                    rowC = None

                fHead = readCsvHeader(fileName)
                if not all([t in fHead for t in self.defaultHead]):
                    raise ValueError(
                        "Can't columns 'x','y','z' and 'var' in %s" % fileName)

                fData = npReadCsv(
                    fileName,
                    columns=self.defaultHead,
                    rowCondition=rowC,
                    **writerArgs)

                fxyz = fData[['x','y','z']][fData['var']==varName]
                if fxyz.size:  # fieldName is in csv-->check equal points
                    fxyz = self.uniquePtsTol(fxyz)
                    xyz = self.points[self.validPoints[varName]]
                    if np.any(np.isin(xyz,fxyz)):
                        raise DuplicateDataError(
                            'Found duplicate for variable %s in file %s'
                            % (_itemString(varNames), fileName))

        dHead = self.defaultHead+frameNames
        cHead = np.array(list(set(fHead+dHead)))
        hTest = np.array([[h in fHead,h in dHead] for h in cHead])
        if np.any(~hTest[:,0]):
            raise FieldRequestError(
                "Can't find columns %s in file %s"
                % (_itemString(cHead[~hTest[:,0]]),fileName))
        if np.any(~hTest[:,1]):
            raise FieldRequestError(
                "Can't find %s in dataset %s"
                % (_itemString(cHead[~hTest[:,1]]),fileName))
        self.saveAsCsv(fileName,varNames=varNames,frameNames=frameNames,
                       noHeader=noHeader, mode='append',**writerArgs)

    @_defaultNoneReplacement
    @_stringInputAsLists
    def toPVFContainer(self,varNames=None,frameNames=None,
                       crdKeyFmt="%.0f", crdKeySep="_"):
        '''Stores points and fielddata to a PVFContainer
        (see L{bae.misc_01.PVFContainer})

        @param varNames: List of existing variable names to store in PVF-C

        @param varNames: List of existing frame names to store in PVF-C

        @param crdKeyFmt: Format string for point keys

        @param crdKeySep: Separator string for point keys
        '''
        #!!! Note that None-inputvalues for varNames and frameNames will be
        #!!! replaced with 'all existing ones' by decorator
        #!!! varNames and/or frameNames will be converted to lists by
        #!!! decorator @_stringInputAsLists

        from bae.misc_01 import PVFContainer

        pvf = PVFContainer()
        pvf.frameIds = frameNames
        ptDict = {}
        ptDataDict = {}
        for ipt,pt in enumerate(self.points):
            ptKey = pvf.getCoordKey(pt,crdKeyFmt=crdKeyFmt,crdKeySep=crdKeySep)
            ptDict[ptKey] = tuple(pt)
            ptDataDict[ptKey] = {}
            for varName in varNames:
                ptVar = [self.data[varName][ipt][fr] for fr in frameNames]
                ptDataDict[ptKey][varName] = ptVar
        pvf.ptData = ptDataDict
        pvf.updatePtInfoWithPtCoordsDict(ptDict)
        return pvf
### END npContainer-Class #####################################################


### test-function npContainer #################################################
def _testNpContainer():
    '''Testfunction for NPContainer
    allready tested:
    creation from array, creation from csv, creation from PVFContainer
    appending points, frames and variables
    requesting points by pointlist
    save to file, append vars/points to file

    @note: the data is created randomly but with a 'fixed' randomseed, so the
        results are replicable
    '''
    #create testdata
    nPts,nFrs,nVars =5,10,10
    varNames = ['var_%d' % d for d in range(nVars)]
    frameNames = ['F%03d' % d for d in range(nFrs)]
    np.random.seed(0)
    rawDataA = np.random.randn(nPts,nFrs,nVars)
    np.random.seed(1)
    xyzA     = 100*np.random.randn(nPts,3)
    np.random.seed(2)
    rawDataB = np.random.randn(int(3),nFrs,nVars)
    np.random.seed(3)
    xyzB     = 100*np.random.randn(int(3),3)
    np.random.seed(4)
    rawDataC = np.random.randn(nPts-1,nFrs,3)
    xyzC     = xyzA[:-1,:]
    np.random.seed(5)
    rawDataD = np.random.randn(nPts,3,nVars)
    xyzD     = xyzA

    def sortRaw(data,xyz):
        #sorting rawData for comparison
        xyz = np.array(map(tuple,xyz), dTypeMapper(['x','y','z']))
        sortIdx = np.argsort(xyz)
        return data[sortIdx],xyz[sortIdx]

    #test creation of npContainer
    dataA = npContainer()
    dataA.pointTol = 1E-6
    dataA.fromArray(rawDataA,xyzA,
                    varNames=varNames,
                    frameNames=frameNames)

    #test appending points
    dataB = dataA.copy()
    dataB.appendPoints(rawDataB,xyzB)

    #test requesting Points
    rawDataASort,_ = sortRaw(rawDataA,xyzA)
    tmpDataB,_ = dataB.getPointSubset(xyzA).reshapeAsVector(varNames=varNames)
    if np.array_equal(tmpDataB,rawDataASort):
        print 'Appending and requesting pointsets passed'
    else:
        print '!!! Appending and requesting pointsets failed !!!'

    #test appending vars
    dataC = dataA.copy()
    nVarNames = ['nVar_%d' % d for d in range(rawDataC.shape[-1])]
    dataC.appendVars(rawDataC, xyzC, nVarNames)

    #test append requesting vars
    rawDataCSort,_ = sortRaw(rawDataC,xyzC)
    rawDataASort,_ = sortRaw(rawDataA,xyzA)
    tmpDataA,_ = dataC.reshapeAsVector(varNames=varNames)
    tmpDataC,_ = dataC.reshapeAsVector(varNames=nVarNames)
    if np.array_equal(tmpDataA,rawDataASort) and\
       np.array_equal(tmpDataC,rawDataCSort):
        print 'Appending and requesting variables passed'
    else:
        print '!!! Appending and requesting variables failed !!!'

    #test appending frames
    dataD = dataA.copy()
    nFrameNames = ['nFrame_%d' % d for d in range(rawDataC.shape[-1])]
    dataD.appendFrames(rawDataD,xyzD,nFrameNames)

    #test requesting vars
    rawDataDSort,_ = sortRaw(rawDataD,xyzD)
    rawDataASort,_ = sortRaw(rawDataA,xyzA)
    tmpDataA,_ = dataD.reshapeAsVector(frameNames=frameNames,
                                       varNames=varNames)
    tmpDataD,_ = dataD.reshapeAsVector(frameNames=nFrameNames,
                                       varNames=varNames)
    if np.array_equal(tmpDataD,rawDataDSort) and\
       np.array_equal(tmpDataA,rawDataASort):
        print 'Appending and requesting frames passed'
    else:
        print '!!! Appending and requesting frames failed !!!'

    #test writing and reading
    testCsv = './testOut.csv'
    if os.path.exists(testCsv):
        os.remove(testCsv)

    dataA.saveAsCsv(testCsv, varNames=varNames[:1])
    dataA.appendToCsv(testCsv, varNames=varNames[1:])
    dataE = npContainer().readFromCsv(testCsv)

    rawDataASort,xyzASort = sortRaw(rawDataA,xyzA)
    tmpDataE,ptE = dataE.reshapeAsVector(varNames=varNames)
    if np.allclose(tmpDataE,rawDataASort):
        print 'Writing and reading datasets passed'
    else:
        print '!!! Writing and reading datasets failed !!!'
    os.remove(testCsv)

    #test data from pvf
    dataF = npContainer().fromPVFContainer(dataA.toPVFContainer())
    tmpDataF,ptF = dataF.reshapeAsVector(varNames=varNames)
    if np.allclose(tmpDataE,tmpDataF):
        print 'PVF conversation passed'
    else:
        print '!!! PVF conversation failed !!!'

###############################################################################
### END npContainer-Class and helper-functions
###############################################################################

def isIronPython():
    return platform.python_implementation().lower().find("ironpython")!=-1
