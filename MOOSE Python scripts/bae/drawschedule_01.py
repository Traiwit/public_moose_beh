"""Read draw schedule from cavesim input file or PCBC like csv.
The data is stored can be processed and exported to PCBC like output
or a cavesim input file.

Note: This interface is not fixed until we reach major-version-number 1!
Earlier versions should allways be copied to the working directory.

Implementation note: xyz coordinates are integer. Draw points are identified
by the xyz-tuple.

Example, create Cavesim input from PCBC:
 >>> from bae.drawschedule_01 import DrawPtDict,writeCavesimGradefileHomogeneous
 >>> from bae.misc_01 import BoundingBox
 >>> from bae.log_01 import msg
 >>>
 >>> #----------------------------------------------------------------
 >>> modelBox = BoundingBox([[14547, 21367, 4457], [16152, 22252, 5942]])
 >>> preambleData = dict(
 >>>     modelTitle="Cadia2012exPost_R51",
 >>>     dropboxPath="D:/NCA_Cadia2012exPost_R51/DROPBOX/COMPUTER_A",
 >>>     box=modelBox,
 >>>     spacing=3,
 >>>     lscale=5,)
 >>> #----------------------------------------------------------------
 >>>
 >>> drawPointDict = DrawPtDict().readPCBCfromCSV(
 ...     "Cadia2012_Q4_drawschedule_Q4.csv", tonnesFirstCol=7)
 >>> drawSchedule = drawPointDict.createDrawSchedule()
 >>> drawSchedule.filter(lambda event: event[2])
 >>> msg("We have %d mining steps, %d draw events, total tonnes: %g ..."
 >>>     % (len(drawSchedule.getMiningSteps()), len(drawSchedule),
 >>>        drawSchedule.getTotalTonnes()))
 >>> drawSchedule.createPreamble(**preambleData)
 >>> drawSchedule.writeCavesimInp()
 >>> msg("Wrote Cavesim input file.")
 >>>
 >>> writeCavesimGradefileHomogeneous(
 >>>     density=2700,
 >>>     initPorosity=0.0,
 >>>     maxPorosity=0.13,
 >>>     box=modelBox )
 >>> msg("Wrote Grade file.")
"""

__todo__= """
- new method to create (populate) DrawSchedule from four lists (dp-coords,
  tonnage, start date, end date), seq-type (weekly, fortnightly, monthly,
  two-monthly, quarterly, half-yearly, yearly), first-date-string.
  Idea:
  1. convert all dates to integer (sequential day number or seconds since
     the epoch or so)
  2. calculate those date-integers for starting points of each seq step
  3. for each draw event in the input data calculate intersection of the
     draw period (start date - end date) with each of the seq steps.
  4. calculate the tonnage-portion for each seq step and add this to the
     DrawSChedule.
  5. merge same date same draw point data points.
- new method to convert to quarterly, to half yearly, ... (accumulating steps)
- reverse of the above: split sequence e.g. from quarterly to monthly. This
  has been done in the 2017 Cadia project I think, maybe Cadia2017010?
- import of v5 CaveSim.inp files create error messages from lines after
  DRAWEVENTS... instead this data should be read and stored in self.footer,
  see __bugs__.
- check that v5 files are being recognized on import and versionCavesim is
  set properly!

for a version 2:
 - remove DrawSchedule.preamble and only use preambleData
 - automatically create preamble from self.preambleData when exporting
   cavesim.inp
 - think big:
   . new data structure: draw points, mining steps and draw events
     with references to draw points and mining steps
   . could be represented by a fieldsCollection object. Or convertible into one
   . think even bigger: it's just one kind of sequence (i.e. prescribed)
     events. Unite with makeHistory / loadHistory?
   . separate Cavesim job -with grid, solver and other data- from the actual
     production schedule. A Cavesim job contains a schedule. But a schedule
     can be stored as PCBC-csv for example without Cavesim data.
"""

__bugs__ = r"""
 - none (of course)
 - ehm, yes: reading from Cavesim.inp and writing back doesn't reproduce the
   footer yet. Reading doesn't update self.footer...
"""

__version__ = "1.29"

_version_history_ = r"""
0.01 GP: new
0.02 GP: changed recognition of first drawschedule line by "DSCHEDULE"
       at the beginning of the line, tonnes are floats now
0.03 GP added file.close statements for Rhino to work properly
        accept float coordinates on input, converting them into integers
0.04 GP modified: DrawPtDict.readPCBCfromCSV() returns self
        added DrawPtDict.createDrawPtDictAccumulated()
        modified: DrawPtDict.createDrawSchedule() preserves numeric draw dates
        fixed: DrawPtDict.readPCBCfromCSV() ignores zero-tonnes fields
        note: version 0.04 has at first been introduced as version 0.02 based
        on an older version!
0.05 GP added grade file creation for CaveSim: writeCavesimGradefile
        added createPreamble
        name change read/writeCAVESIMinp -> ...CavesimInp
0.06 GP added DrawPtDict.changeDates
0.07 GP added DrawPtDict.createDrawPtDictAccumulated(dateIndex=None) for
        accumulated up to current mining step
        incompatible change to default argument (was: dateIndex=-1)
1.00 GP nothing changed, just fixed the interface to make it usable as
        bae.drawschedule_01
1.01 GP fixed: DrawSchedule.readCavesimInp() skipped some early draw events
1.02 GP added DrawPtDict.drawDates and check in DrawPtDict.getAllDrawDates()
1.03 GP added DrawPtSchedule.add()
        modified: DrawSchedule.createPreamble preserves trailing comment-only
        block.
        added optional dateToStep-argument to DrawPtDict.createDrawSchedule()
1.04 GP added: storing date strings in DrawSchedule.stepToDate and converting
        mining step to date string.
1.05 GP added: strip leading and trailing white space from date strings on
        input from PCBC csv file
1.06 GP added DrawSchedule.sort(), DrawPtDict.createDrawSchedule() now creates
        positions as lists, not tuples, give sensible error message on
        PCBC file parsing error and skip this line.
1.07 GP DrawSchedule.readCavesimInp now returns self so you can do:
        schedule = DrawSchedule().readCavesimInp(...)
1.08 GP added DrawPtDict.joinSameDates()
1.09 GP added Cavesim format v4
1.10 GP added: recognize and store versionCavesim when reading input file
1.11 GP added: new class CaveSimSimple
1.12 GP added handling and persistence of optional eventId attribute on
        DrawEvent and DrawPtScheduleItem objects
1.13 GP added draw point radius attributes
1.14 TR fixed: collecting Rhino-GUIDs in DrawEvent.__init__
1.15 TR added v5.11 to DrawSchedule.writeCavesimInp
1.16 GP added DrawPtSchedule.addEvent
1.17 GP added DrawSchedule.addEventsWithDates
1.18 GP changed: DrawSchedule.createFooter function merged into createPreamble
        (this is an incompatible change but I don't expect problems.)
1.19 GP added writeCavesimGradefile, taken from Tobias'
        interpolateNearestToCSBlockModel.py
1.20 GP added CAVESIM.INP-export for Cavesim version 6.3
1.21 GP added DrawSchedule.readPosDateTonnesCsv()
1.22 TR added fs4-export type and zRotation
1.23 GP added DrawSchedule.getMiningDates(),
        added toStep argument to DrawSchedule.getTotalTonnes,
        added DrawPtDict.getHodFieldOnGrid()
1.24 GP changed default DrawSchedule export type (versionCavesim) to "fs4".
1.25 GP added radius, zRotation and drawPtAreaForHod arguments to
        DrawPtDict.getHodFieldOnGrid()
1.26 GP added heightOffset argument to DrawPtDict.getHodFieldOnGrid()
1.27 GP added DrawSchedule.fromPosDateTonnesList(), made
        DrawSchedule.addEventsWithDates() accept another DrawSchedule object.
1.28 GP fixed sumupMultipleDraweventsAtSamePointAndFrame to also work with
        gaps in the sequence of mining steps
1.29 GP added takePartialCells option to DrawPtDict.getHodFieldOnGrid()
"""

import csv
import re
from itertools import izip
from math import ceil, floor, sin, cos, pi
from copy import deepcopy

try:
    import pandas as pd
except ImportError:
    pd = None

from collections import defaultdict, OrderedDict
from bae.log_01 import msg, MsgTicker
from bae.misc_01 import newFile, BoundingBox
from bae.vecmath_01 import vector_plus, vector_minus, vector_scale
from bae.field_01 import Field


class DrawEvent(list):
    """One single draw event. It's basically a list of three items:
     1. position, a list [x,y,z] of three integers
     2. mining step, integer
     3. tonnes bogged, integer (or float?)

    @ivar eventId: (optional) When created by
    L{bert.lib.drawaro_01.getDrawScheduleFromRhinoObjs} then this attribute
    is present and contains the guid of the correponding Rhino-mesh-object.

    @ivar radius: (optional) an x,y,z-tuple of half-widths of the draw point.
    L{DrawSchedule.readCavesimInp}() sets this attribute.
    The radius attribute will transparently be transfered from und to
    L{DrawPtScheduleItem}s of a L{DrawPtDict}-object.
    The attribute can be added manually as well.
    @ivar zRotation: (optional) z-rotation of draw region in deg
    """
    def __init__(self, *args, **kwargs):
        """Accept position, mining step, tonnes as arguments
        or DrawEvent object (copy constructor)

        @kwarg radius: a triple of half-width of the draw point in
          x,y,z-direction

        @kwarg zRotation: z-rotation of draw region in deg
        """

        # single argument of type DrawEvent -> copy constructor
        if len(args)==1 and isinstance(args[0], DrawEvent):
            list.__init__(self, args[0])
            try:
                self.radius = args[0].radius
            except AttributeError:
                pass
            try:
                self.zRotation = args[0].zRotation
            except AttributeError:
                pass
            try:
                self.eventId = args[0].eventId
            except AttributeError:
                pass


        # three positional arguments + radius kwarg
        else:
            list.__init__(self, args)
            try:
                self.radius = tuple(kwargs["radius"])
            except (KeyError, TypeError):
                # KeyError if there is no keyword argument "radius"
                # TypeError if a keyword argument is no tuple, but e.g. None
                pass
            try:
                self.zRotation = kwargs["zRotation"]
            except (KeyError):
                # KeyError if there is no keyword argument "zRotation"
                pass


class DrawSchedule(list):
    """for processing draw schedules

    This is a chronologically ordered list of individual draw events, objects
    of type L{DrawEvent}.
    This view represents the CaveSim input file.

    The CaveSim input file contains a preamble with some more data than just
    the draw sequence. See the description of the instance variable preamble
    below and the method self.L{createPreamble}().

    @ivar preamble: First lines of the CaveSim input file stored as a list of
       strings. Will be updated by self.L{readCavesimInp}(). You can use the
       method self.L{createPreamble}() to generate it or supply it manually.

    @ivar preambleData: This dict contains all arguments to
       self.L{createPreamble}(). That are additional Cavesim parameters.

    @ivar stepToDate: Dictionary {mining step number: date string}.
       Mining step numbers are integer starting at 1. This dictionary will be
       filled by readCavesimInp if it finds lines like "** 1 : Y2012_M04" or
       so.
       Note: It's guaranteed that this disctionary exists. But not that it
       contains any entries, let alone entries for each mining step. It's
       suggested to access with the dict.get method instead of "[]"
       -- self.__getitem__():  sched.stepToDate.get(stepNb, "Step %d"%stepNb)

    @ivar xwid: half width of Cavesim draw point in x-direction in the
       CAVESIM.inp file
    @ivar ywid: half width of Cavesim draw point in y-direction in the
       CAVESIM.inp file
    @ivar zwid: width of Cavesim draw point in z-direction in the
       CAVESIM.inp file
    @ivar versionCavesim: "v2", "v4", "v5", "v6", "fs4"; will be determined by
       self.L{readCavesimInp}, otherwise defaults to "fs4".
    """

    def __init__(self, *args):
        """Takes initializer like a list object or takes another DrawSchedule
        object an then serves as copy-constructor.

        Makes sure that the actual items of self are of type L{DrawEvent}.

        Initializes mandatory instance attributes.
        """
        initFlag = False
        if len(args)==1:
            # convert items in first (and only) argument to DrawPtScheduleItem
            list.__init__(self, (DrawEvent(x) for x in args[0]))

            if isinstance(args[0], DrawSchedule):
                # copy constructor
                initFlag = True
                for name in [   # deep copy
                        "preamble","preambleData","stepToDate","keywordDict"]:
                    setattr(self, name, deepcopy(getattr(args[0], name)))
                for name in [  # shallow copy
                        "stepToDateHeader", "stepToDateTemplate",
                        "rexStepToDate", "rexModelTitle",
                        "xwid", "ywid", "zwid", "versionCavesim"]:
                    setattr(self, name, getattr(args[0], name))

        if not initFlag:

            # lines including the newline characters
            self.preamble = list()
            self.footer = list()

            # this dict contains all arguments to self.createPreamble()
            self.preambleData = dict()

            # dict {mining step number: date string}
            self.stepToDate = dict()

            # format of comment block describing the mining steps
            self.stepToDateHeader = "** Includes following mining steps."
            self.stepToDateTemplate = "**   %d : %s"
            self.rexStepToDate = re.compile(r"^\*\*\s+(\d+)\s*:\s*(\S+)")

            # model title in the Cavesim input file
            self.rexModelTitle = re.compile(
                r"^\s*\*\s*MODEL\s*TITLE\s*:\s*(.+?)\s*$")

            # Cavesim input file keyword to self.preambleData-(key, type)-tuple
            self.keywordDict = {
                "DROPBOXPATH": ("dropboxPath", str),
                "SPACING": ("spacing", int),
                "LSCALE": ("lscale", int),
                "NGRADEBLOCKS": ("nbGradeBlocks", int),
                "LDIM": ("ldim",
                         lambda x: [int(float(xx)) for xx in x.split()]),
                "LORIGIN": ("lorigin",
                            lambda x: [int(float(xx)) for xx in x.split()]),
                "FLOWSOLVER": ("flowSolver", int),}

            # some defaults
            self.xwid = 5
            self.ywid = 5
            self.zwid = 5
            self.versionCavesim = "fs4"


    @classmethod
    def fromPosDateTonnesList(cls, events):
        """Creates a new DrawSchedule from a list of events as (position,
        date-string, tonnes)-tuples.

        @param events: Iterable of event tuples consisting of (position,
           date, tonnes). position must be a xyz floats-triple, date must be
           a string that sorts chronologically and tonnes must be a float.
        """

        self = cls()

        if not isinstance(events, (tuple, list)):
            # we need to iterate twice, so an iterator is not suitable
            events = list(events)
        newDates = sorted(set(x[1] for x in events))

        self.stepToDate.update(
            (i+1, x) for i, x in enumerate(newDates))
        dateToStep = dict((s,i) for i,s in self.stepToDate.iteritems())
        self.extend(
            # position, mining step, tonnes
            DrawEvent(xyz, dateToStep[date], tonnes)
            for xyz, date, tonnes in events)
        return self


    def addEventsWithDates(self, events):
        """Adds events to the schedule according to their date rather than
        their mining step.

        @param events: Iterable of event tuples consisting of (position,
           date, tonnes). position must be a xyz floats-triple, date must be
           a string that sorts chronologically and tonnes must be a float.

           Or another DrawSchedule object.
        """
        if not isinstance(events, DrawSchedule):
            events = DrawSchedule.fromPosDateTonnesList(events)

        otherSteps = [(i, events.stepToDate[i])
                      for i in events.getMiningSteps()]
        otherSteps.append((None, None))  # add stop marker
        otherSteps = iter(otherSteps)

        oldSteps = [(i, self.stepToDate[i]) for i in self.getMiningSteps()]
        oldSteps.append((None, None))  # add stop marker
        oldSteps = iter(oldSteps)

        # reschedule, creating gaps for extra steps
        newToOldSteps = dict()
        newStepToDate = dict()
        dateToStep = dict()  # only dates in events

        newstep = 1
        otherstep, otherdate = otherSteps.next()
        oldstep, olddate = oldSteps.next()
        while otherstep or oldstep:
            if otherdate == olddate:
                newToOldSteps[newstep] = [oldstep,]
                dateToStep[otherdate] = newstep
                newStepToDate[newstep] = olddate
                # advance both
                newstep += 1
                otherstep, otherdate = otherSteps.next()
                oldstep, olddate = oldSteps.next()

            elif otherdate and (otherdate < olddate or olddate is None):
                dateToStep[otherdate] = newstep
                newStepToDate[newstep] = otherdate
                # advance other
                newstep += 1
                otherstep, otherdate = otherSteps.next()

            elif olddate and (olddate < otherdate or otherdate is None):
                newToOldSteps[newstep] = [oldstep,]
                newStepToDate[newstep] = olddate
                # advance old (self)
                newstep += 1
                oldstep, olddate = oldSteps.next()

        self.reschedule(newToOldSteps)
        self.stepToDate = newStepToDate
        for ev in events:
            # print "######"
            # print ev[1]
            # print events.stepToDate[ev[1]]
            # print dateToStep[events.stepToDate[ev[1]]]
            ev[1] = dateToStep[events.stepToDate[ev[1]]]
        self.extend(events)
        self.sort()


    def readCavesimInp(self, filename="CAVESIM.inp", ignoreZeroTonnes=True):
        """Import draw schedule and the preamble from a CAVESIM.inp file.

        @param ignoreZeroTonnes: if true then don't store draw events with zero
           tonnes.

        @Return: self so you can do:
          schedule = DrawSchedule().readCavesimInp(...)
        """

        # reset schedule
        del self[:]

        # open file, init
        msg("Reading data from %s." % filename)
        inp = open(filename, "r")
        cnt = 0

        # ignore beginning until draw schedule (incl. header line)
        msg("Skipping until draw schedule.")
        pos = inp.tell()
        inStepToDateBlock = False
        # Note: must use file.readline for file.tell to work
        # Here again: Shitty IronPython does not work:
        # ... instead of: for line in iter(inp.readline, ""):
        #     have to do (for IronPython):
        while 1:
            line = inp.readline()
            if not line: break
            # end of IronPython-workaround for "for line in iter(....):"
            if line.startswith("DSCHEDULE") or line.startswith("DRAWEVENT"):
                break
            elif line.startswith(self.stepToDateHeader):
                inStepToDateBlock = True
            elif inStepToDateBlock:
                # possibly read mining step to date relation in comment
                res = self.rexStepToDate.match(line)
                if res:
                    miningStep, dateStr = res.groups()
                    try:
                        miningStep = int(miningStep)
                    except ValueError:
                        pass
                    else:
                        self.stepToDate[miningStep] = dateStr
                else:
                    # not a step-to-date line, block ends here
                    inStepToDateBlock = False

            if not inStepToDateBlock:
                # store line as preamble
                self.preamble.append(line)

            # possibly update self.preambleData
            line = line.strip()
            while not inStepToDateBlock and line:

                # model title
                res = self.rexModelTitle.match(line)
                if res:
                    self.preambleData["modelTitle"] = res.group(1)
                    break

                # other parameters
                line = line.split(None, 1)
                if len(line)!=2 or line[0] not in self.keywordDict:
                    break
                key, type_ = self.keywordDict[line[0].upper()]
                self.preambleData[key] = type_(line[1])

                # no second iteration
                break

            # next line
            pos = inp.tell()
            cnt+=1

        # create self.preambleData["box"]
        try:
            lscale = self.preambleData["lscale"]
            ldim = self.preambleData["ldim"]
            origin = self.preambleData["lorigin"]
            nigiro = [x0 + lscale*dx for x0, dx in izip(origin, ldim)]
            self.preambleData["box"] = [origin, nigiro]
            del self.preambleData["ldim"], self.preambleData["lorigin"]
        except KeyError:
            pass

        # start reading the actual draw schedule
        msg("Now reading draw schedule starting at line %d, byte position %d"
            % ((cnt+1), pos))
        # rewind to beginning of this line
        # note: according to python-docs it's ok to use file.seek in these
        # circumstances as file.seek flushes the read-ahead buffer that
        # file.next() uses (i.e. in in for line in inp)
        inp.seek(pos)

        ticker = MsgTicker("processing line %d")
        totalTonnes = 0
        firstLoop = True
        warnDifferentRadius = False
        for line in inp:
            cnt += 1
            ticker.msg(cnt)
            row = line.split()

            # check nb of columns, might be 11 or 12 or 13
            if firstLoop:
                if len(row)==11:
                    versionCavesim = "v2"
                    colNbExpected = 11
                elif len(row)==12:
                    versionCavesim = "v4"
                    colNbExpected = 12
                elif len(row)==13:
                    versionCavesim = "v6"
                    colNbExpected = 13
                else:
                    msg("ERROR: Expected 11, 12 or 13 rows, found %d in line"
                        " %d. Ignoring line." % (len(row), cnt+1))
                    continue
                self.versionCavesim = versionCavesim
                msg("Found %d columns in the draw schedule, this is for"
                    " Cavesim %s." % (colNbExpected, versionCavesim))
            else:
                if len(row)!=colNbExpected:
                    msg("ERROR: Expected %d rows, found %d in line %d. Ignoring"
                        " line." % (colNbExpected, len(row), cnt+1))
                    continue

            # type of draw point dimensions
            radiusStr = row[5:8]
            if any("." in x for x in radiusStr):
                radiusType = float
            else:
                radiusType = int


            # convert all to numbers
            try:
                drawId = int(row[1])
                xyz = tuple(int(float(x)) for x in row[2:5])
                radius = tuple(radiusType(x) for x in radiusStr)
                miningStep = int(row[8])
                tonnes = float(row[9])
                if int(row[10])==1:
                    # draw zone type: 1==rotated
                    zRotation = float(row[11])
                else:
                    zRotation = 0
            except ValueError:
                msg("ERROR: Could not convert all numbers in line %d: %s"
                    % (cnt, row))
                continue

            if firstLoop:
                firstLoop = False
                # get draw point dimensions from first line
                self.xwid, self.ywid, self.zwid = radius
            elif ((self.xwid-radius[0])**2
                  + (self.ywid-radius[1])**2
                  + (self.zwid-radius[2])**2) > 0.9:
                warnDifferentRadius = True

            self.append(DrawEvent(xyz, miningStep, tonnes, radius=radius,
                                  zRotation=zRotation))

            totalTonnes += tonnes

        del ticker
        inp.close()   # required for Rhino Python
        del inp

        if warnDifferentRadius:
            msg("WARNING: We have different draw point half widths in this"
                " draw schedule. Some checks are disabled therefore.")

        msg("Finished processing %d lines." % cnt)
        msg("... total tonnes %d in %d draw events" % (totalTonnes, len(self)))
        return self

    def readPosDateTonnesCsv(self, input_):
        """Reads a csv file with columns x, y, z, date, tonnes.
        Currently a header line is expected but essentially ignored.
        A warning is printed if the first three columns are not labelled
        x,y,z.

        (We might want to add later that the header line is evaluated to
        identify the actual column order.)

        @param input_: can be a file name string or anything that csv.reader
         accecpts as first argument. I.e. any object which supports the
         iterator protocol and returns a string each time its next() method
         is called --file objects and list objects are both suitable.

        @returns: self, so this construct is possible:
          sched = DrawSchedule().readPosDateTonnesCsv("sched.csv")
        """

        # open input (file)
        if isinstance(input_, basestring):
            input_ = open(input_, "rb")
        tab = csv.reader(input_)

        # ignore header line
        header = tab.next()
        if tuple(x.lower() for x in header[:3]) != ("x","y","z"):
            msg("WARNING: Expected header line with x,y,z coordinates, then"
                " date and tonnes. Got: %s" % header)

        # reset schedule
        del self[:]

        # read data ignoring empty lines
        events = [
            (map(float, x[0:3]), x[3], float(x[4]))
            for x in tab if x and any(x)]
        self.addEventsWithDates(events)

        # close input file
        del tab
        return self

    def writeCavesimInp(self, filename="CAVESIM.inp", versionCavesim=None):
        r"""export draw schedule

        Usage:
         >>> sched = DrawSchedule()...  # some data to export
         >>> sched.createPreamble(**sched.preambleData, versionCavesim="v5")
         >>> sched.writeCavesimInp(filename, versionCavesim="v5")

        @param versionCavesim: "v2", "v4", "v6" or "fs4". If not given then
           self.L{versionCavesim} will be used.

        @Note: It's advisable to call
           self.L{createPreamble}(**self.L{preambleData})
           beforehand to update nb of draw events and so on. See example code
           above.
        @Note: For fs4 only quad-shaped drawzone geometries are implemented yet
           (see legacy- and rotated non-legacy format in
           U{http://192.168.31.24/wordpress/archives/3009#more-3009<BE-wordpress>})
        """
        output = newFile(filename)

        if versionCavesim is None:
            versionCavesim = self.versionCavesim

        if versionCavesim not in ("v2","v4","v5","v6","fs4"):
            msg("WARNING: Unrecognized Cavesim version %s. Instead assuming"
                " %s." % (versionCavesim, self.versionCavesim))
            versionCavesim = self.versionCavesim

        if versionCavesim=="v2":
            columnNames = (
                "SCHEDULE","ID","X","Y","Z","XWID","YWID","ZWID","MS_NUM",
                "TONNES_BOGGED","DRAWZONETYPE")
            # event-id, x, y, z, xwid, ywid, zwid, mining-step, tonnes
            dataLineTemplate = "DSCHEDULE"+"\t%d"*10+"\n"
        elif versionCavesim=="v4":
            columnNames = (
                "SCHEDULE","ID","X","Y","Z","XWID","YWID","ZWID","MS_NUM",
                "TONNES_BOGGED","SOG","DRAWZONETYPE")
            # event-id, x, y, z, xwid, ywid, zwid, mining-step, tonnes
            dataLineTemplate = "DRAWEVENT"+"\t%d"*10+"\t2\n"
        elif versionCavesim=="v5":
            columnNames = (
                "SCHEDULE","ID","X","Y","Z","XWID","YWID","ZWID","PERIOD",
                "TONNES","SOG","TYPE")
            # event-id, x, y, z, xwid, ywid, zwid, mining-step, tonnes
            dataLineTemplate = "DRAWEVENT"+"\t%d"*10+"\t2\n"
        elif versionCavesim=="v6":
            columnNames = (
                "SCHEDULE","ID","X","Y","Z","XWID","YWID","ZWID","PERIOD",
                "TONNES","SOG","TYPE","TEMPLATE")
            # event-id, x, y, z, xwid, ywid, zwid, mining-step, tonnes,
            # ... blast template
            dataLineTemplate = "DRAWEVENT"+"\t%d"*10+"\t2\t0\n"
        elif versionCavesim=="fs4":
            columnNames = (
                "SCHEDULE","ID","X","Y","Z","XWID","YWID","ZWID","PERIOD",
                "TONNES","DRAWZONETYPE","ROTATION")
            # event-id, x, y, z, xwid, ywid, zwid, mining-step, tonnes
            dataLineTemplate = "DRAWEVENT"+"\t%d"*11+"\n"
        else:
            # should never get here, choose a default before...
            raise NotImplementedError("Cavesim version <%s> not implemented"
                                      % versionCavesim)

        output.writelines(self.preamble)
        output.write(self.stepToDateHeader+"\n")
        output.writelines((self.stepToDateTemplate+"\n") % x
                          for x in sorted(self.stepToDate.iteritems()))
        output.write("** DRAW SCHEDULE **\n")
        output.write("* %s\n" % "\t".join(columnNames))

        for cnt, drawEvent in enumerate(self):
            xyz, miningStep, tonnes = drawEvent
            # drawEvent.radius may not exist, may be None or have sensible data
            try:
                radius = drawEvent.radius
            except AttributeError:
                radius = None
            if radius is None:
                radius = (self.xwid, self.ywid, self.zwid)

            try:
                zRotation = drawEvent.zRotation
                drawZoneType = 1  # has rotation
            except AttributeError:
                zRotation = None
                drawZoneType = 0  # legacy
            if zRotation is None:
                zRotation = 0

            if versionCavesim == "fs4":
                output.write(dataLineTemplate % (
                    cnt+1, xyz[0], xyz[1], xyz[2],
                    radius[0], radius[1], radius[2],
                    miningStep, round(tonnes),
                    drawZoneType, int(round(zRotation))))
            else:
                output.write(dataLineTemplate % (
                    cnt+1, xyz[0], xyz[1], xyz[2],
                    radius[0], radius[1], radius[2],
                    miningStep, round(tonnes), 0))

        output.writelines(self.footer)
        output.close()   # required for Rhino Python
        del output


    def createPreamble(
            self, modelTitle, dropboxPath, box, spacing=2, lscale=5,
            nbGradeBlocks=1, xwid=None, ywid=None, zwid=None,
            versionCavesim="v2", flowSolver=2, horizProbability=90):
        """Create / overwrite the preamble and footer of the CAVESIM.inp file.
        If there is already a preamble in self then preserve comments after the
        last non-comment line, i.e. preserve trailing comment-only-block.

        @param modelTitle: something like Cadia2015_R02
        @param dropboxPath: e.g. D:/NCA_Cadia2015_R01/DROPBOX/COMPUTER_A
        @param spacing: Number of Cavesim lattice cells per coupling grid
            interval on each axis. I.e. 2 means each coupling grid cell
            consists of 2^3=8 CaveSim cells.
        @param box: location of the model box, should correspond to the
            box defined in the GRADE.txt file. A L{bae.misc_01.BoundingBox}
            or equivalent pair of points. Should be larger than the
            coupling box by at least one coupling-box-cell.
        @param lscale: grid size of the Cavesim lattice in meters;
            "never use values larger 5 m", why?
        @param nbGradeBlocks: number of entries in the grade block model
        @param xwid: (and ywid and zwid) set x-, y- half width and z-width of
            drawpoints
        @param versionCavesim: "v2" or "v4" or "v5"
        @param flowSolver: 1 or 2. 1...slow Newton solver for grades evaluation
            2...fast diffusion solver for caveabilty
        @param horizProbability: probability of horizontal movements.
            Valid range: 45 to 90. Default from Glenn: 60. Fred recommends: 90.
        """
        # calculate number of cells from box and grid size "lscale"
        ldim = list()
        for i in range(3):
            dx = box[1][i] - box[0][i]
            nx = float(dx) / lscale
            if int(nx) < nx:
                nx = int(nx)+1
            else:
                nx = int(nx)
            ldim.append(nx)
        nbCells = ldim[0]*ldim[1]*ldim[2]

        # number of mining steps
        nbLastStep = self.getMiningSteps()[-1]

        # number of draw events
        nbDrawEvents = len(self)

        if xwid:
            self.xwid = xwid
        if ywid:
            self.ywid = ywid
        if zwid:
            self.zwid = zwid

        # Adjust drawpoint widths to match lscale = cavesim Mesh size
        if self.xwid*2<lscale:
            msg("Adjusted xwid (drawpt half x-width) from %s to %g"
                % (self.xwid, 0.5*lscale))
            self.xwid = 0.5*lscale
        if self.ywid*2<lscale:
            msg("Adjusted ywid (drawpt half y-width) from %s to %g"
                % (self.ywid, 0.5*lscale))
            self.ywid = 0.5*lscale
        if self.zwid<lscale:
            msg("Adjusted zwid (drawpt z-width) from %s to %g"
                % (self.zwid, lscale))
            self.zwid = lscale

        # Check: drawpoint size should usually be at least same size as
        # coupling grid size
        if self.xwid*2<lscale*spacing:
            msg("WARNING: xwid (drawpt half x-width) %s smaller than half"
                " coupling grid size %g" % (self.xwid, 0.5*lscale*spacing))
        if self.ywid*2<lscale*spacing:
            msg("WARNING: ywid (drawpt half y-width) %s smaller than half"
                " coupling grid size %g" % (self.ywid, 0.5*lscale*spacing))
        if self.zwid<lscale*spacing:
            msg("WARNING: zwid (drawpt z-width) %s smaller than coupling"
                " grid size %g" % (self.zwid, lscale*spacing))

        # from already existing preamble preserve comments after last
        # non-comment line
        oldComments = []
        for line in self.preamble:
            if not line.startswith("*"):
                oldComments = []
                continue
            oldComments.append(line)

        if versionCavesim=="v2":
            ngradefields = 15
            nLatticeFields = 19
            strNminsNdp = (
                "* NUMBER OF MINING STEPS\nNMINS %d"
                "\n\n* NUMBER OF DRAWPOINTS\nNDRAWZONES %d"
                % (nbLastStep, nbDrawEvents))
            strAutoStaAndSpacing = "SPACING %.1f" % spacing
            strGradeFile = (
                "NGRADEBLOCKS %d\nNGRADEFIELDS %d"
                % (nbGradeBlocks, ngradefields))

        elif versionCavesim=="v4":
            ngradefields = 19
            nLatticeFields = 19
            strNminsNdp = (
                "* NUMBER OF MINING STEPS\nNMINS %d"
                "\n\n* NUMBER OF DRAWPOINTS\nNDRAWZONES %d"
                % (nbLastStep, nbDrawEvents))
            strAutoStaAndSpacing = "AUTOSTART 1\nSPACING %.1f" % spacing
            strGradeFile = (
                "NGRADEBLOCKS %d\nNGRADEFIELDS %d"
                % (nbGradeBlocks, ngradefields))

        elif versionCavesim == "v5":
            nLatticeFields = 21
            ngradefields = 27
            strAutoStaAndSpacing = (
                "AUTOSTART 1\nSPACING %d\nSPACINGVEL %d" % (spacing, spacing))
            strNminsNdp = (
                # nb of draw points moves to footer in v5.11
                "* NUMBER OF MINING STEPS\nNMINS %d" % nbLastStep)
            strGradeFile = (
                "NGRADEBLOCKS %d\nNGRADEFIELDS %d"
                % (nbGradeBlocks, ngradefields))
            strFileIO = "READVTKS 1"

            strCaveability = "\n".join((
                "FREEFIELD 0",
                "PROPTYPE 3", "CRATE 0", "HRCRIT 9999",
                "PERIODDAYS 30", "CAVEFREQ 1", "BFACTORDEF 1.18",
                "CAVEANGLE 999", "ZSURFACEMIN 9999", "CRATEDEF 999",
                "DWIDTH 20", "DXWIDTH 16", "DYWIDTH 10",
                "XCAVED 1", "YCAVED 1", "ZCAVED 1",
                "CUSTOM 0", "SOCKET 0",
            ))
            strDrawZoneControl = ""
        elif versionCavesim == "v6":
            nLatticeFields = 21
            ngradefields = 27
            strAutoStaAndSpacing = (
                "AUTOSTART 1\nSPACING %d\nSPACINGVEL %d" % (spacing, spacing))
            strNminsNdp = (
                # nb of draw points moved to footer in v5.11
                "* NUMBER OF MINING STEPS\nNMINS %d" % nbLastStep)
            strGradeFile = "NGRADEFIELDS %d" % ngradefields
            strFileIO = "READVTKS 1\nIEZOUT 0"

            strCaveability = "\n".join((
                "READMRMR 1",
                "DRAWCOLUMNS 0",
                "PERIODDAYS 30", "CAVEFREQ 1",
                "FREE FIELD 0",
                "BFACTORDEF 1.2",
                "CAVEANGLE 999", "ZSURFACEMIN 9999", "CRATEDEF 999",
                "DWIDTH 20", "DXWIDTH 20", "DYWIDTH 20",
                "XCAVED 1", "YCAVED 1", "ZCAVED 1",
                "CUSTOM 0", "SOCKET 0",
                ))

            strDrawZoneControl = "\n".join((
                "SUPPRESSZERO 1",
                "WRITEBUCKETS 0",  # new in v6.3
                "FMR 3",
                "FMAX 0.2",
                "FRAGDIST 20",
                "FRAGFDIST 40",
                "FINEGENRATE 0.001",
                "READBLASTED 0",
                ))
        elif versionCavesim == "fs4":
            # create header
            from datetime import datetime
            now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
            self.preamble = list()
            self.preamble.extend('** FS4 DRAWSCHEDULE\n')
            self.preamble.extend(
                "** MODEL TITLE: {modelTitle}\n".format(**vars()))
            self.preamble.extend('** created %s\n' % now)
            self.preamble.extend("** number of drawevents: %d\n" % len(self))
            self.preamble.extend("** note: XWID, YWID, ZWID are half width\n")
            return
        else:
            raise ValueError(
                "versionCavesim %s not implemented yet." % versionCavesim)

        # Only Cavesim - versions go here
        # (versionCavesim == "fs4" finished earlier)
        self.preamble = list()
        self.preamble.extend(("""
* CAVESIM (3D CAVING SIMULATOR)
* COPYRIGHT (2010) Dr Glenn Sharrock
* INPUT PARAMETERS

* MODEL TITLE: {modelTitle}

DROPBOXPATH {dropboxPath}
{strAutoStaAndSpacing}

NINTERATIONS 1000000

* Number of Cores or CPU's
NPROC 20

* LX, LY, LZ; ... {nbCells} cells
LDIM {ldim[0]} {ldim[1]} {ldim[2]}

* LXO, LYO, LZO
LORIGIN {box[0][0]} {box[0][1]} {box[0][2]}

LSCALE {lscale}

* NUMBER OF LATTICE FIELDS (+1)
NLATTICEFIELDS {nLatticeFields}
ZSMIN 5
ZSMAX 95

{strNminsNdp}

* OUTPUT INTERVAL FOR GRAPHICS & PARTICLE DELETION ROUTINE
GFREQ 20

DRAWFREQ 5
ITRX 100000

* GRADE FIELDS
{strGradeFile}
"""
            ).format(**vars()).splitlines(True))

        if versionCavesim[0]=="v" and versionCavesim >= "v5":
            self.preamble.extend(("""
* SPECIFIC FOR v5.11
PTYPEDEF 1
MMAX 99

* 1=Newtonian (required for grade calculation)
* 2=Particle Diffusion (fast)
* 3=Void Diffusion (for grade calculations, new in v6.0)
FLOWSOLVER {flowSolver}
MAXVDP 1
MAXPORCHANGE 0.5

{strFileIO}
READCAVEPTS 0
POINTSOUT 0
CAVESHAPEMODE 1
DXFOUT 0
READMARKERS 0
MARKERSOUT 0
BMODELOUT 0
BMODELSTAGE 999
BMODELSPACING 1

CAVEABILITY 0
{strCaveability}
CRATER 0

DL 0
{strDrawZoneControl}

"""
            ).format(**vars()).splitlines(True))

            ABC = "ABCDEFGHIJKLMNOP"
            self.preamble.extend(
                ('%sCAL %d\n' % (a, horizProbability)) for a in ABC)

            self.preamble.extend('\n* INITIAL VIEW SETTING\n')
            for ii, dd in enumerate(['X','Y','Z']):
                c = (box[1][ii] + box[0][ii])/2
                self.preamble.extend(['%sCENTRE %d\n' % (dd,c),])
            self.preamble.extend("PLANESTRIKE 90\n\n")

        # append old comments (presumably containing the step-to-date table)
        # right before the draw schedule at the end of the parameter section
        self.preamble.extend(oldComments)

        # create footer
        self.footer = list()
        if versionCavesim[0]=="v" and versionCavesim >= "v5":
            nbDrawEvents = len(self)
            self.footer.extend([
                "\n* NUMBER OF DRAWPOINTS\n",
                "NDRAWZONES {nbDrawEvents}\n\n".format(**vars())])

            # This doesn't seem to be necessary,
            # does anybody know what it's for?
            # Yes, it's for the cave propagation algorithm that we don't use.
            # ABC = "ABCDEFGHI"
            # self.footer.extend(['span%s%s=10000\n' % (a,b)
            #                for a in ABC for b in ABC])


    def sort(self):
        """Sort the drawschedule chronologically. For example after you
        appended draw events.

         >>> mainSchedule.extend(extraProductionSchedule)
         >>> mainSchedule.sort()

        @Note: Remember that you may have to L{reschedule} the DrawSchedule
        before you merge it into another one to have consistent mining
        step numbers.
        """
        list.sort(self, key=lambda x:x[1])

    def reschedule(self, newToOldSteps, joinedDate=-1):
        """reassign mining steps in the draw schedule. Updates self.stepToDate
        as well.
        @param newToOldSteps: dict {new mining step nb: list of old mining
            step nbs}. mining steps that never occur in any list of old mining
            steps will be removed from the schedule.

        @param joinedDate: reschedule also updates self.stepToDate.
            In case many steps are joined into one we have to decide which
            date to assign to the merged step. joinedDate=0 takes the first,
            joinedDate=-1 takes the last.
            Other values are currently not implemented.

        @Note: You can savely merge multiple old mining steps into one single
        new mining step. That's one purpose of this function and the reason why
        the values of the newToOldSteps-dict are lists.
        """
        # debugLevel=1 because diagnostic output (other than progress
        # indication) should generally be done on the top level
        msg("Rescheduling drawSchedule: Originally we have %d mining events."
            % len(self), debugLevel=1)
        miningStepToEvents = defaultdict(list)
        for drawEvent in self:
            xyz, miningStep, tonnes = drawEvent
            try:
                radius = drawEvent.radius
            except AttributeError:
                radius = None

            try:
                zRotation = drawEvent.zRotation
            except AttributeError:
                zRotation = None

            miningStepToEvents[miningStep].append([xyz, tonnes, radius,
                                                   zRotation])

        msg("... occuring in %d mining steps."
            % len(miningStepToEvents), debugLevel=1)
        miningStepToEvents = dict(
            (newStep,
             sum((miningStepToEvents[oldStep] for oldStep in oldSteps), []))
            for newStep, oldSteps in newToOldSteps.iteritems())
        msg("After rescheduling we'll have %d mining steps."
            % len(miningStepToEvents), debugLevel=1)

        msg("Merging events on same draw point and updating the schedule.")
        ticker = MsgTicker("...processing %%d/%d mining steps."
                           % len(miningStepToEvents))
        del self[:]
        for miningStep in sorted(miningStepToEvents):
            ticker.tick()
            posTonsRad = OrderedDict()
            for xyz, tonnes, radius, zRot in miningStepToEvents[miningStep]:
                try:
                    tonsRad = posTonsRad[tuple(xyz)]
                except KeyError:
                    posTonsRad[tuple(xyz)] = [tonnes, radius, zRot]
                else:
                    tonsRad[0] += tonnes
                    tonsRad[1] = radius if tonsRad[1] is None else tonsRad[1]
                    tonsRad[2] = zRot if tonsRad[2] is None else tonsRad[2]
            self.extend(
                DrawEvent(list(xyz), miningStep, tonnes, radius=radius,
                          zRotation=zRot)
                for xyz, (tonnes, radius, zRot) in posTonsRad.iteritems())
        del ticker

        # updating stepToDate attribute as well
        if self.stepToDate:
            newStepToDate = dict(
                (newStep, self.stepToDate[oldSteps[joinedDate]])
                for newStep, oldSteps in newToOldSteps.iteritems())
            self.stepToDate = newStepToDate

    def sumupMultipleDraweventsAtSamePointAndFrame(self):
        miningSteps = self.getMiningSteps()
        # debugLevel=1 because diagnostic output (other than progress
        # indication) should generally be done on the top level
        msg("before, there are %d drawevents in %d miningsteps"
            % (len(self), len(miningSteps)), debugLevel=1)
        self.reschedule(newToOldSteps=dict(
            (i,[i,]) for i in miningSteps))
        msg("after,  there are %d drawevents in %d miningsteps"
            % (len(self), len(self.getMiningSteps())), debugLevel=1)

    def filter(self, filterFunction):
        """keep only draw events for which filterfunction evaluates to True
        when called with a particular L{DrawEvent} object

        To remove zero-tonnes draw events do the following:
         >>> myDrawSched.filter(lambda event: event[2])

        @param filterFunction: A function taking one DrawEvent object as single
        argument. Draw events for which filterFunction evaluates to False will
        be removed from the sequence.
        """
        self[:] = filter(filterFunction, self)

    def getTotalTonnes(self, toStep=None):
        """Sum up all tonnes up to the given step
        @param toStep: int mining step number or date string as in
           self.stepToDate
        """
        if toStep is None:
            #...remember: DrawEvent(xyz, miningStep, tonnes)
            return sum(event[2] for event in self)
        else:
            if isinstance(toStep, basestring):
                dateToStep = dict((d,s) for s,d in self.stepToDate.iteritems())
                stepNb = dateToStep.get(toStep, None)
                if stepNb is None:
                    raise ValueError(
                        "ERROR: getTotalTonnes got toStep argument %s. This"
                        " date is not stored in stepToDate." % toStep)
                toStep = stepNb
            #...remember: DrawEvent(xyz, miningStep, tonnes)
            return sum(event[2] for event in self if event[1]<=toStep)

    def getMiningSteps(self):
        """Return a sorted list of mining step numbers.

        To get the dates use L{getMiningDates}.
        """
        return sorted(set(miningStep for xyz, miningStep, tonnes in self))

    def getMiningDates(self):
        """Return a sorted list of the dates of all mining steps.

        To get the mining step numbers use L{getMiningSteps}.
        """
        return [self.stepToDate.get(step, "Step %d" % step)
                for step in self.getMiningSteps()]

    def createDrawPtDict(self):
        """Creates a L{DrawPtDict} object from self.
        Uses self.stepToDate to replace mining steps by meaningful date
        strings. If not all mining steps are included in self.stepToDate
        will keep the mining step number as date.
        """
        drawPtDict = defaultdict(DrawPtSchedule)

        # check self.stepToDate contains all mining steps
        allMiningSteps = self.getMiningSteps()
        if all((x in self.stepToDate) for x in allMiningSteps):
            stepToDate = self.stepToDate
        else:
            if len(self.stepToDate):
                msg("WARNING: Mining-step to date-string relation <stepToDate>"
                    " does not contain relations for all mining steps.")
                msg("... We have mining steps from %d to %d."
                    % (min(allMiningSteps), max(allMiningSteps)))
                minStepToDate = min(self.stepToDate.iteritems())
                maxStepToDate = max(self.stepToDate.iteritems())
                msg("... stepToDate contains %d->%s to %d->%s."
                    % (minStepToDate[0], minStepToDate[1],
                       maxStepToDate[0], maxStepToDate[1]))
                msg("... Therefore DrawSchedule.createDrawPtDict() will"
                    " deliver mining step numbers as dates!")
            # dummy-dict for no conversion
            stepToDate = dict((i, i) for i in allMiningSteps)

        # transfer mining steps and convert dates
        for drawEvent in self:
            xyz, miningStep, tonnes = drawEvent
            miningStep = stepToDate[miningStep]

            kwargs = dict()
            try:
                kwargs["eventId"] = drawEvent.eventId
            except AttributeError:
                pass
            try:
                kwargs["radius"] = drawEvent.radius
            except AttributeError:
                pass
            try:
                kwargs["zRotation"] = drawEvent.zRotation
            except AttributeError:
                pass

            item = DrawPtScheduleItem(miningStep, tonnes, **kwargs)
            drawPtDict[tuple(xyz)].append(item)

        # should not be necessary, just in case: sort chronologically
        for msTonnes in drawPtDict.itervalues():
            msTonnes.sort()

        return DrawPtDict(drawPtDict)


#----------------------------------------------------------------------

class DrawPtScheduleItem(list):
    """Basically a [mining step, tonnes] -list. Serves as items in
    L{DrawPtSchedule} - objects.

    @ivar eventId: (optional) identifier of a draw-event.
    E.g. can contain the guid of the correponding Rhino-mesh-object if it
    originates from such. See
    L{bert.lib.drawaro_01.getDrawScheduleFromRhinoObjs}

    @ivar radius: (optional) triple of half-width of the draw point in
       x,y,z-direction

    @ivar zRotation: (optional) z-rotation of draw region in deg

    @Note: Mining step values (first item of a L{DrawPtScheduleItem}-object)
    must be in a chronological order when sorted alphabetically. Otherwise
    some functions will not work properly.
    """
    def __init__(self, *args, **kwargs):
        """accept mining (step, tonnes,  as arguments
        or DrawEvent object (copy constructor)

        Possible usage with kwargs, optionally add eventId:
         >>> it = DrawPtScheduleItem(step=3, tonnes=5245.5, radius=(5,5,5))

        Possible usage with two positional args: step, tonnes:
         >>> it = DrawPtScheduleItem(3, 5245.5, radius=(5,5,5), tonnes=guid)

        Possible usage as copy constructor:
         >>> it = DrawPtScheduleItem(...)
         >>> it2 = DrawPtScheduleItem(it)

        @kwarg step: mining step, integer or string (e.g. date like Y2021_M05)
        @kwarg tonnes: tonnes bogged, integer or float
        @kwarg eventId: (optional) guid of the correponding Rhino-mesh-object
        @kwarg radius: (optional) a triple of half-width of the draw point in
            x,y,z-direction
        @kwarg zRotation: (optional) z-rotation of draw region in deg
        """

        # single argument of type DrawPtScheduleItem -> copy constructor
        if (len(args)==1 and isinstance(args[0], (list, tuple))
                and len(args[0])==2):
            list.__init__(self, args[0])
            try:
                self.radius = args[0].radius
            except AttributeError:
                pass
            try:
                self.zRotation = args[0].zRotation
            except AttributeError:
                pass
            try:
                self.eventId = args[0].eventId
            except AttributeError:
                pass
        else:
            # two positional args: step, tonnes
            if len(args)==2:
                list.__init__(self, args[:2])

            # only kwargs
            elif len(args)==0:
                self.extend((kwargs["step"], kwargs["tonnes"]))
            else:
                raise ValueError(
                    "Unexpected set of arguments to constructor of"
                    " DrawPtScheduleItem: args=%s, kwargs=%s" % (args, kwargs))

            # make sure radius is a tuple or None
            try:
                radius = kwargs["radius"]
            except KeyError:
                radius = None
            else:
                radius = tuple(radius) if radius else None
            self.radius = radius

            try:
                self.zRotation = kwargs["zRotation"]
            except KeyError:
                self.zRotation = None
            try:
                self.eventId = kwargs["eventId"]
            except KeyError:
                pass


class DrawPtSchedule(list):
    """Chronologically sorted list of L{DrawPtScheduleItem} objects.
    For compatibility (mining step, tonnes) tuples work as well.
    Objects of this class are the values in a L{DrawPtDict} - dict.

    @Note: Mining step values (first item of each L{DrawPtScheduleItem}-object)
    must be unique (occur only once in self) and must be in a chronological
    order when sorted alphabetically. Otherwise some functions will not work
    properly.
    """
    def getTotalTonnes(self):
        return sum(tonnes for (ms, tonnes) in self)

    def addEvent(self, eventToAdd, idxInSelf=0):
        """Add a single draw event.
        If the same mining step isn't in the schedule yet it's being added.
        If this mining step already is in the schedule then self will be
        modified such that the new tonnage is added to the previously existing.

         >>> drawPtDict = DrawPtDict()
         >>> for pt, data in drawPtData:
         >>>     sched = DrawPtSchedule()
         >>>     for step, tonnes in data:
         >>>         sched.addEvent(DrawPtScheduleItem(step, tonnes))
         >>>     drawPtDict[tuple(map(int, pt))] = sched

        @param eventToAdd: a L{DrawPtScheduleItem}-object to be added.
        @param idxInSelf: not for general use. Start index of the search
           in self. self.add() makes use of it.
        @returns: idxInSelf for internal use by self.add().
        """

        miningStep, tonnes = eventToAdd
        added = False

        while idxInSelf<len(self):
            oldEvent = self[idxInSelf]
            if oldEvent[0]==miningStep:
                if isinstance(oldEvent, DrawPtScheduleItem):
                    oldEvent[1] += tonnes
                else:
                    # for compatibility with older versions: accept event
                    # as (step, tonnes)-tuple
                    tonnes += oldEvent[1]
                    self[idxInSelf] = (miningStep, tonnes)
                idxInSelf += 1
                added = True
                break
            elif oldEvent[0] > miningStep:
                self.insert(idxInSelf, eventToAdd)
                idxInSelf += 1
                added = True
                break
            idxInSelf += 1
            continue
        if not added:  # in this case: idxInSelf==len(self)
            self.append(eventToAdd)
            idxInSelf += 1
        return idxInSelf

    def add(self, scheduleToAdd):
        """Add the production given as L{DrawPtSchedule}-object.
        Mining steps that have not been in the schedule before are being added.
        Mining steps that are already in the schedule self will be modified
        such that the new tonnage is added to the previously existing.

         >>> drawPtDict = DrawPtDict().readPCBCfromCSV(...)
         >>> for pt, data in drawPtData:
         >>>     pt = map(int, pt)
         >>>     try:
         >>>         sched = drawPtDict[pt]
         >>>     except KeyError:
         >>>         sched = DrawPtSchedule()
         >>>         drawPtDict[tuple(pt)] = sched
         >>>     sched.add(DrawPtSchedule(
         >>>         DrawPtScheduleItem(step, tonnes)
         >>>         for step, tonnes in data))

        @param scheduleToAdd: a L{DrawPtSchedule}-object to be added.
        """
        idxInSelf = 0
        for event in scheduleToAdd:
            idxInSelf = self.addEvent(event, idxInSelf)
        return


class DrawPtDict(dict):
    """
    This is the "PCBC"-like view of a draw schedule.

    It's a dict: The key is a draw point position, a tuple of integers (x,y,z).
    The value is a L{DrawPtSchedule}-object: A chronologically sorted list of
    L{DrawPtScheduleItem}-objects.

    Suppose you've got a file "myPCBC.csv" looking like this (CaveSim step
    numbers in the first line, then tonnes/step in the corresponding columns):
     >>> X,Y,Z,1,2,3
     >>> 10134,5124,252,500,1000,2000
     >>> 10152,5131,252,0,600,1000
     >>> 10170,5138,252,0,200,1500

    Then you could create a CaveSim input file like this:
    Note: We don't check for duplicate positions.
    Second note: There is a method DrawPtDict.readPCBCfromCSV() to accomplish
    this task.
     >>> drawPointDict = DrawPtDict()
     >>> input = csv.reader(open("myPCBC.csv", "rb"))
     >>> steps = input.next()[3:]
     >>> for row in input:
     >>>     drawPtPos = map(int, row[0:3])
     >>>     tonnes = row[3:]  # note: strings need to be converted later!
     >>>     drawPointDict[tuple(drawPtPos)] = DrawPtSchedule(
     >>>          (step, int(ton))
     >>>          for step, ton in zip(steps, tonnes) if ton )
     >>> del input  # close csv input file
     >>> # export to CaveSim
     >>> drawSced = drawPointDict.createDrawSchedule()
     >>> drawSced.preamble = " .... "  # all the other CaveSim data
     >>> drawSced.writeCavesimInp()

    @ivar drawDates: Chronological list of draw dates. Optional, attribute
    not guaranteed to exist at all. Expected to contain the same values as
    the mining steps in all the L{DrawPtScheduleItem}s.
    Filled by self.readPCBCfromCSV() from the header line. And by
    L{getAllDrawDates}() from the L{DrawPtScheduleItem}s.

    @Note: Mining step values (first item of a L{DrawPtScheduleItem}-object)
    must be in a chronological order when sorted alphabetically. Otherwise
    some functions will not work properly.
    """

    def getAllDrawDates(self):
        """returns a list containing all draw dates in chronological
        (actually: alphabetical) order.
        """
        allDrawDatesSet = set()
        for dateTonnesList in self.itervalues():
            allDrawDatesSet.update(dateTonnes[0]
                                   for dateTonnes in dateTonnesList)
        newDrawDates = sorted(allDrawDatesSet)
        if (hasattr(self, "drawDates")
            and newDrawDates != self.drawDates):
            msg("WARNING! Draw dates changed inexpectedly. One possible reason"
                " might be that the alphanumeric order of draw dates is not"
                " chronological!\n"
                "Current dates: %s\n"
                "Previously stored dates: %s"
                % (newDrawDates, self.drawDates))
            self.drawDates = newDrawDates
        return newDrawDates

    def readPCBCfromCSV(self, filename="PCBC_input.csv", coordsCol=None,
                        tonnesFirstCol=3, transVector=(0,0,0)):
        """Import draw schedule from a PCBC csv file. Overwriting existing
        drawpoints at the same location (integer coordinates).

        Format: It's a csv file. The header line gives the mining steps, the
        first three fields are being ignored. First three columns show the draw
        point coordinates followed by the tonnes per mining step.

        Rows and columns (below is a conceptual table)::
         x       | y       | z       | 1       | 2       | 3       | 4
         4321.2  | 1423.7  | -534.3  | 513.4   | 742.7   | 0       | 0
         4327    | 1425    | -534.3  | 0       | 513.4   | 742.7   | 0
         4334    | 1427    | -534    | 0       | 313     | 513     | 242
         ...

        @param filename: of the PCBC file to be read
        @param coordsCol: An integer identifies the first (x-)column by
          its zero-based index. If None (not specified) and tonnesFirstCol==3
          then use the first three columns for the point coordinates; i.e. in
          this case we expect the simplest format possible:
          x, y, z, tonnes, tonnes,...
          If None and tonnesFirstCol!=3 then try to find by column titles,
          i.e. search the header line for x, y, z (capital letters as well).
        @param tonnesFirstCol: Zero based index of the first column that
          holds tonnes to draw. And draw date in the head line.
        @param transVector: Translation vector for draw point coordinates.
          Default is (0,0,0)

        @Note: Coordinates are rounded to integers.
        @Returns: self so that drawPointDict = DrawPtDict().readPCBCfromCSV(...)
         is possible.
        """
        inputFile = open(filename, "rb")
        input = csv.reader(inputFile)
        headline = input.next()

        # if coordsCol is a number then this is the column index for
        if isinstance(coordsCol, int):
            coordsCols = [coordsCol, coordsCol+1, coordsCol+2]

        # if coordsCol not given
        # 1. heuristics: if tonnes/dates start in col 3 then coords just before
        elif coordsCol is None and tonnesFirstCol==3:
            coordsCols = [0, 1, 2]

        # if coordsCol not given
        # 2. heuristics: find x,y,z in head line
        elif coordsCol is None:
            coordsCols = list()
            for c in "xyz":
                try:
                    coordsCols.append(headline.index(c))
                except ValueError:
                    try:
                        coordsCols.append(headline.index(c.upper()))
                    except ValueError:
                        break
            # end of loop to find x,y,z
            if len(coordsCols) < 3:
                msg("ERROR: Could not determine drawpoint coordinates from"
                    " header line in file %s" % inputFile)
                inputFile.close()   # required for Rhino Python
                del input  # close csv input file
                return self
        msg("Loading coordinates from columns with index %s."
            % (", ".join(str(x) for x in coordsCols)))

        self.drawDates = [x.strip() for x in headline[tonnesFirstCol:]]
        msg("Starting with column index %d, found %d draw dates in the header"
            " line, first is %s, last is %s."
            % (tonnesFirstCol, len(self.drawDates),
               self.drawDates[0], self.drawDates[-1]))
        for cnt, row in enumerate(input):
            try:
                drawPtPos = tuple(int(round(float(row[i]))+transVector[i])
                                  for i in coordsCols)
            except ValueError as exc:
                msg("ERROR on line %d of input file %s:\n%s"
                    % (cnt+1, filename, exc.args))
                msg("Ignoring this line: %s" % row)
                continue
            # note: strings need to be converted later!
            tonnes = row[tonnesFirstCol:]
            self[drawPtPos] = DrawPtSchedule(
                # note: tonnes be converted from strings here!
                DrawPtScheduleItem(date, float(ton))
                for date, ton in zip(self.drawDates,tonnes)
                # ton may be None, "", or 0 but float("") raises Exception
                if ton and float(ton))

        inputFile.close()   # required for Rhino Python
        del inputFile  # close csv input file
        msg("Finished reading %s with %d lines after the header. There are %d"
            " distinct draw point positions" % (filename, (cnt+1), len(self)))
        return self

    def writePCBCtoCSV(self, filename="PCBC_output.csv", zeroValue=""):
        """export draw schedule as a PCBC csv file
        @param zeroValue: Value for dates without draw.
        """
        allDrawDates = self.getAllDrawDates()
        outputFile = open(filename, "wb")
        output = csv.writer(outputFile)
        output.writerow(["X","Y","Z"]+allDrawDates)
        for drawPtCoords,dateTonnesList in self.iteritems():
            raw = list(drawPtCoords)
            dateTonnes = dict(dateTonnesList)
            for date in allDrawDates:
                raw.append(dateTonnes.get(date, zeroValue))
            output.writerow(raw)
        outputFile.close()

    def createDrawSchedule(self, dateToStep=None):
        """Returns a L{DrawSchedule}-object that can then be exported for
        CaveSim.

        Draw dates are being converted to integer numbers: First simple type
        conversion is tried, i.e. newdate = int(float(olddate)).
        If that fails then secondly drawdates are enumerated in alphabetically
        sorted order starting with one.

        You can call self.L{changeDates}() beforehand for more sophisticated
        conversions.

        Example to have mining steps for each month regardless whether there
        is mining or not. I.e. there will be gaps in the sequence of mining
        steps if there is no activity in partical months:
         >>> from bae.misc_01 import getDateStringSeq
         >>> drawPtDict = DrawPtDict().readPCBCfromCSV(...)
         >>> dates = drawPtDict.getAllDrawDates()
         >>> dateToStep = dict(
         >>>     (date, i+1)
         >>>     for i, date in enumerate(
         >>>         getDateStringSeq(dates[0], dates[-1])))
         >>> drawSchedule = drawPtDict.createDrawSchedule(dateToStep)

        @param dateToStep: optional dict {date: integer step number} for the
        conversion of dates in self (might be strings) to the numeric mining
        step number in the resulting L{DrawSchedule}-object. The dict must
        contain all dates as can be determined by self.getAllDrawDates().

        It's save to provide more dates than necessary, dates that have no
        mining activity are simply ignored.
        """
        if not dateToStep:
            dates = self.getAllDrawDates()
            try:
                numDates = [int(float(d)) for d in dates]
            except (ValueError, TypeError):
                dateToStep = dict((date, i+1) for i,date in enumerate(dates))
            else:
                dateToStep = dict((date, n)
                                  for date, n in izip(dates, numDates))

        # collect all mining events in groups for each mining step
        # in dict miningStepToEvents {mining step number: list of events}
        miningStepToEvents = defaultdict(list)
        for pos, schedule in self.iteritems():
            for item in schedule:
                date, tonnes = item
                miningStep = dateToStep[date]
                event = DrawEvent(list(pos), miningStep, tonnes)

                # if present transfer eventId
                try:
                    event.eventId = item.eventId
                except AttributeError:
                    pass

                # if present transfer radius
                try:
                    event.radius = item.radius
                except AttributeError:
                    pass

                # if present transfer zRotation
                try:
                    event.zRotation = item.zRotation
                except AttributeError:
                    pass

                # store event
                miningStepToEvents[miningStep].append(event)

        # build the DrawSchedule object
        res = DrawSchedule()
        for miningStep in sorted(miningStepToEvents):
            res.extend(miningStepToEvents[miningStep])

        # add comments in the preamble
        res.stepToDate.update(
            (dateToStep[date], date)
            for date in dateToStep
            if dateToStep[date] in miningStepToEvents)
        return res

    def createDrawPtDictAccumulated(self, dateIndex=None):
        """Create a new DrawPtDict-object with accumulated draw events.
        Either all drawevents get the accumulated tonnes up to the current
        mining step. Or all drawevents for one drawpoint will be accumulated
        into one draw event.

        @param dateIndex: If not specified or None then all draw events for
            a particular draw point will get accumulated tonnage up to the
            current mining step.

            If specified and not None then all draw events for a single
            draw point will be accumulated into one single draw event.

            For the accumulated draw event take the mining step
            identified by this index in the draw-point-schedule-list (of the
            particular draw point).

            I.e. specify dateIndex=0 if you want the cumulative draw event
            happen when the first of the old individual draw events of a
            particular draw point occurred.

            Specify dateIndex=-1 if you want the cumulative draw event happen
            at the last draw event.
        """
        if dateIndex is None:
            res = DrawPtDict()
            for pos, schedule in self.iteritems():
                accTonnes = 0.0
                accSched = DrawPtSchedule()
                for event in schedule:
                    step, tonnes = event
                    accTonnes += tonnes
                    accSched.append(DrawPtScheduleItem(
                        step, accTonnes,
                        radius=event.__dict__.get("radius"),
                        eventId=event.__dict__.get("eventId"),
                        zRotation=event.__dict__.get("zRotation")
                        ))

                res[pos] = accSched
        else:
            res = DrawPtDict()
            for pos, schedule in self.iteritems():
                selectedOldEvent = schedule[dateIndex]
                newEvent = DrawPtScheduleItem(
                    step=selectedOldEvent[0],
                    tonnes=schedule.getTotalTonnes(),
                    radius=selectedOldEvent.__dict__.get("radius"),
                    eventId=selectedOldEvent.__dict__.get("eventId"),
                    zRotation=selectedOldEvent.__dict__.get("zRotation"))
                res[pos] = DrawPtSchedule((newEvent,))
        return res

    def joinSameDates(self):
        """Merge drawpoint events in the same mining step on the same
        draw point.
        Does not modify self but returns a new DrawPtDict with the draw events
        merged.
        """
        res = DrawPtDict()
        radiusMismatch,zRotMismatch = list(), list()
        for pos, schedule in self.iteritems():
            stepEventList = defaultdict(list)
            for event in schedule:
                stepEventList[event[0]].append(event)

            # join events for the same step
            newschedule = DrawPtSchedule()
            radiusMismatchFlag, zRotMismatchFlag = False, False
            for step in sorted(stepEventList):
                eventList = stepEventList[step]
                tonnes = sum(x[1] for x in eventList)
                radiusSet = set(getattr(x, "radius", None) for x in eventList)
                zRotSet   = set(getattr(x, "zRotation", None)
                                for x in eventList)
                radiusMismatchFlag = radiusMismatchFlag or (len(radiusSet)>1)
                zRotMismatchFlag   = zRotMismatchFlag or (len(zRotSet)>1)

                newschedule.append(DrawPtScheduleItem(
                    step, tonnes, radius=radiusSet.pop(),
                    zRotation=zRotSet.pop(),
                    # copy eventId if there is only one event
                    # ... and if it has eventId
                    eventId=(len(eventList)==1
                             and eventList[0].__dict__.get("eventId"))))
            if radiusMismatchFlag:
                radiusMismatch.append(pos)
            if zRotMismatchFlag:
                zRotMismatch.append(pos)
            res[pos] = newschedule

        if radiusMismatch:
            maxDisp = 20
            if len(radiusMismatch)>maxDisp:
                msgText = ("Some affected drawpoints: %s" % ", ".join(
                    (str(x) for x in radiusMismatch[:maxDisp])))
            else:
                msgText = ("Affected drawpoint(s): %s" % ", ".join(
                    (str(x) for x in radiusMismatch)))
            msg("WARNING: Draw point radius / half widths mismatch for %d draw"
                " points.\n%s" % (len(radiusMismatch), msgText))

        if zRotMismatch:
            maxDisp = 20
            if len(zRotMismatch)>maxDisp:
                msgText = ("Some affected drawpoints: %s" % ", ".join(
                    (str(x) for x in zRotMismatch[:maxDisp])))
            else:
                msgText = ("Affected drawpoint(s): %s" % ", ".join(
                    (str(x) for x in zRotMismatch)))
            msg("WARNING: Draw point zRotation mismatch for %d draw"
                " points.\n%s" % (len(zRotMismatch), msgText))

        return res

    def changeDates(self, oldToNewTuples, errorOnMissing=True):
        """Perform a 1 to 1 replacement of mining step labels.

        @param oldToNewTuples: Iterable of (old mining step label, new mining
          step label) tuples. (It will be converted to a dict on entry.)
        @param errorOnMissing: if True then raise a KeyError-exception in
          case you find a mining step in the schedule that's not listed as old
          mining step in oldToNewTuples. Otherwise such dates will be left
          untouched.

        @Note: See also the L{DrawSchedule.reschedule}-function that can be
        used to merge multiple old mining steps into one single
        new mining step.
        """
        oldToNewTuples = dict(oldToNewTuples)
        for pos, schedule in self.iteritems():
            for i, item in enumerate(schedule):
                # ( miningStep, tonnes )
                oldDate, tonnes = item
                if errorOnMissing:
                    newDate = oldToNewTuples[oldDate]
                else:
                    newDate = oldToNewTuples.get(oldDate, oldDate)
                try:
                    # in case item is a DrawPtScheduleItem object
                    item[0] = newDate
                except TypeError:
                    # otherwise assign updated tuple
                    schedule[i] = (newDate, tonnes)
        return

    def getHodFieldOnGrid(self, grid, miningStep, swellFactor, density,
                          fieldName="height_of_draw",
                          radius=None, zRotation=0.0, drawPtAreaForHod=None,
                          heightOffset=0.0, takePartialCells=False):
        """Create a L{bae.field_01} field object for the given miningStep with
        1 inside the Height Of Draw (HOD) and 0 elsewhere.

        Example: Get a field indicating HOD:
         >>> from bae.mesh_01 import MeshStructuredPoints
         >>> from bae.vtk_02 import FieldsOnStructPointGrid as Vtk
         >>>
         >>> drawPtDict = DrawSchedule().readCavesimInp(...).createDrawPtDict()
         >>> grid = MeshStructuredPoints(....)
         >>> hodField = drawPtDict.getHodFieldOnGrid(
         >>>     grid, "Y2021_M05", 0.2, 2700)
         >>>
         >>> vtk = Vtk(mesh=grid)
         >>> vtk.updateFields(hodField)
         >>> vtk.toVtk("hod_Y2021_M05.vtk")

        The cross-sectional area (as seen in plan view) of the draw zone is
        given by the drawpoint dimensions. It's a rectangular box, possibly
        rotated around the z-axis.

        @param grid: a L{MeshStructuredPoints} object
        @param miningStep: mining step, integer or string (e.g. date like
          Y2021_M05). All L{DrawPtScheduleItem}-items are considered with
          item[0] <= miningStep.
        @param density: in kg / m^3
        @param swellFactor: e.g. 20% -> 0.2
          height of draw per tonnes
          = (1.0 + swell factor) / (rho * draw point area * swell factor)
        @param fieldName: fieldName attribute of the field to be returned
        @param radius: Overrules draw region dimensions for all draw points.
          If not given then the radii from the draw schedule (stored in the
          radius attribute of each L{DrawPtScheduleItem} object are taken.
          If specified then this value is taken for all draw points and the
          values stored in each L{DrawPtScheduleItem} object are ignored.
          Note that the HOD column starts at the bottom of the draw region
          i.e. at z = z_drawpoint - radius[2].
        @param zRotation: rotation angle in deg of the rectangluar drawpoint
          area in the x-y-plane.
          Only considered if the radius argument is specified. If the radius
          (drawpoint dimensions) is taken from the draw schedule (i.e. the
          L{DrawPtScheduleItem} objects) then this argument is ignored.
        @param drawPtAreaForHod: if specified use this as cross sectional area
          to calculate the HOD. If not given then take the radii --specified
          or from the draw schedule-- to calculate this area. In square metres.
        @param heightOffset: This height is added to each HOD column. Can be
          used to account for the undercut heigth.
        @param takePartialCells: If True then consider cells that partially
          intersect the draw region as "in". If False then only consider cells
          "in" if the cell centroid is inside the draw region. Defaults to
          False.

          Note: This option together with rotated draw regions might be flawed!
          Further checks needed!

        @returns: L{bae.field_01.Field}-object (position="structPt",
            dataType="scalar") containing 1 for "in", 0 for "out"
        """

        warnings = set()  # collect warnings then write them once in the end

        # abbreviations
        orig = grid.origin
        gridPtNb = grid.gridPtNb
        dx = grid.spacing
        st = grid.strides
        eps = 1E-4

        # height = tonnesToHeightFactor * tonnes / draw_point_area
        tonnesToHeightFactor = (
            1000.0*(1.0+swellFactor)/(float(density)*swellFactor))

        # initialize resulting field with all zeros
        result = Field.classFromPosType(fieldName, "structPt", "scalar")(
            [0]*len(grid))

        # flag indicating if the size of the draw region has to be determined
        # for each draw point from the draw schedule
        drawRegionDimFromSchedule = (radius is None)

        # where do we get the draw point area from for the HOD-calculation?
        drawPtAreaFromSchedule = False
        if drawPtAreaForHod is None:  # area not given, calc form radius
            if drawRegionDimFromSchedule:
                drawPtAreaFromSchedule = True
            else:
                drawPtAreaForHod = 4.0*radius[0]*radius[1]

        # iterate over all draw points
        for pos, schedule in self.iteritems():
            tonnes = sum(tonnes for step, tonnes in schedule
                         if step <= miningStep)
            if not tonnes:
                continue

            if drawRegionDimFromSchedule:
                # determine radius (draw point half widths)
                try:
                    radii = set(x.radius for x in schedule)
                    radius = radii.pop()
                except AttributeError:
                    raise ValueError(
                        "Draw point width not given for each draw event on"
                        " drawpoint %s." % str(pos))
                except KeyError:
                    # Gero thinks this should never happen...
                    raise ValueError(
                        "No draw point width given on draw point %s."
                        % str(pos))
                if radii:
                    raise ValueError(
                        "Different radii given in draw point %s." % str(pos))
                del radii

                # determine zRotation
                try:
                    zRotations = set(
                        y for y in (
                            getattr(x, "zRotation", None) for x in schedule)
                        if y is not None)
                    zRotation = zRotations.pop()
                except KeyError:
                    # no zRotation given
                    zRotation = None
                if zRotations:
                    raise ValueError(
                        "Different z-rotations given in draw point %s." % pos)
                del zRotations

            if drawPtAreaFromSchedule:
                # determine draw point area for HOD from current radius
                drawPtAreaForHod = 4.0*radius[0]*radius[1]

            # determine z extent in xyz and ijk-grid coords
            z0 = float(pos[2]) - radius[2]
            z1 = (z0 + heightOffset
                  + tonnesToHeightFactor*tonnes/drawPtAreaForHod)
            if not takePartialCells:
                k0 = int(ceil((z0-orig[2])/dx[2]))
                k1 = min(int(ceil((z1-orig[2])/dx[2])), gridPtNb[2])
            else:
                # takePartialCells==True
                k0 = int(ceil(((z0-orig[2])/dx[2])-0.5))
                k1 = min(int(ceil(((z1-orig[2])/dx[2])+0.5)), gridPtNb[2])
            if k0<0:
                warnings.add("Lower draw point boundary outside grid.")
                k0 = 0

            if zRotation:
                # abbreviations
                s = sin(pi*zRotation/180)
                c = cos(pi*zRotation/180)

                #--- determine corner coordinates of rotated rectangle
                # corners in xyz coords: C[i][j] = j-th comp of i-th corner
                C = [[radius[0], radius[1]],
                     [-radius[0], radius[1]],
                     [-radius[0], -radius[1]],
                     [radius[0], -radius[1]]]
                # rotation through rot-mat [[c, s], [-s, c]]
                C = [[C[0][0]*c - C[0][1]*s, C[0][0]*s + C[0][1]*c],
                     [C[1][0]*c - C[1][1]*s, C[1][0]*s + C[1][1]*c],
                     [C[2][0]*c - C[2][1]*s, C[2][0]*s + C[2][1]*c],
                     [C[3][0]*c - C[3][1]*s, C[3][0]*s + C[3][1]*c]]
                # add centre of draw point: C now in global xyz coords
                C = [[x[0]+pos[0], x[1]+pos[1]] for x in C]
                # corners in grid coords: C[i][j] = j-th comp of i-th corner
                C = [[(x[0]-orig[0])/dx[0], (x[1]-orig[1])/dx[1]]
                     for x in C]

                #--- rotate indexes of corner points
                # after that corners in matrix C are in this order:
                # [max y, min x, min y, max x]
                indexOffset = int((zRotation % 360) / 90)
                C = [C[(i-indexOffset) % 4] for i in range(4)]

                #---
                # max range in i- direction
                if not takePartialCells:
                    i0 = int(ceil(C[1][0]))
                    i1 = 1+int(floor(C[3][0]))
                else:
                    i0 = int(ceil(C[1][0]-0.5))
                    i1 = 1+int(floor(C[3][0]+0.5))

                # check i- / x- boundaries
                if i0<0 or i1>gridPtNb[0]:
                    warnings.add("Draw point boundary outside grid in x-dir.")
                    i0 = max(i0, 0)
                    i1 = min(i1, gridPtNb[0])

                for i in range(i0, i1):

                    # j0: bottom j
                    if (abs(C[3][0]-C[2][0])<=eps or abs(C[2][0]-C[1][0])<=eps):
                        # edges vertical or horizontal
                        # then take C[2][1] as y_min
                        if not takePartialCells:
                            j0 = int(ceil(C[2][1]))
                        else:
                            j0 = int(ceil(C[2][1] - 0.5))
                    elif i<C[2][0]:
                        # "left" of C[2] -> boundary C[1]--C[2]
                        xi = (i-C[1][0])/(C[2][0]-C[1][0])
                        if not takePartialCells:
                            j0 = int(ceil(C[1][1]*(1.0-xi) + C[2][1]*xi))
                        else:
                            j0 = int(ceil(C[1][1]*(1.0-xi) + C[2][1]*xi - 0.5))
                    else:
                        # "right" of C[2] -> boundary C[2]--C[3]
                        xi = (i-C[2][0])/(C[3][0]-C[2][0])
                        if not takePartialCells:
                            j0 = int(ceil(C[2][1]*(1.0-xi) + C[3][1]*xi))
                        else:
                            j0 = int(ceil(C[2][1]*(1.0-xi) + C[3][1]*xi - 0.5))

                    # j1: top j
                    if (abs(C[3][0]-C[0][0])<=eps or abs(C[0][0]-C[1][0])<=eps):
                        # edges vertical or horizontal
                        # then take C[0][1] as y_max
                        if not takePartialCells:
                            j1 = 1+int(floor(C[0][1]))
                        else:
                            j1 = 1+int(floor(C[0][1] + 0.5))
                    elif i<C[0][0]:
                        # "left" of C[0] -> boundary C[1]--C[0]
                        xi = (i-C[1][0])/(C[0][0]-C[1][0])
                        if not takePartialCells:
                            j1 = 1+int(floor(C[1][1]*(1.0-xi) + C[0][1]*xi))
                        else:
                            j1 = 1+int(floor(C[1][1]*(1.0-xi) + C[0][1]*xi+0.5))
                    else:
                        # "right" of C[0] -> boundary C[0]--C[3]
                        xi = (i-C[0][0])/(C[3][0]-C[0][0])
                        if not takePartialCells:
                            j1 = 1+int(floor(C[0][1]*(1.0-xi) + C[3][1]*xi))
                        else:
                            j1 = 1+int(floor(C[0][1]*(1.0-xi) + C[3][1]*xi+0.5))

                    # check j- / y- boundaries
                    if j0<0 or j1>gridPtNb[1]:
                        warnings.add(
                            "Draw point boundary outside grid in y-dir.")
                        j0 = max(j0, 0)
                        j1 = min(j1, gridPtNb[1])

                    # set points inside bounds to 1
                    for j in range(j0, j1):
                        offs = i*st[0] + j*st[1]
                        result[offs+k0*st[2]:offs+k1*st[2]:st[2]] = [1]*(k1-k0)

            else:  # no zRotation
                # determine x,y ranges in grid coordinates
                if not takePartialCells:
                    i0 = int(ceil(float(pos[0]-radius[0]-orig[0])/dx[0]))
                    i1 = 1+int(floor(float(pos[0]+radius[0]-orig[0])/dx[0]))
                    j0 = int(ceil(float(pos[1]-radius[1]-orig[1])/dx[1]))
                    j1 = 1+int(floor(float(pos[1]+radius[1]-orig[1])/dx[1]))
                else:
                    i0 = int(ceil(float(pos[0]-radius[0]-orig[0])/dx[0]-0.5))
                    i1 = 1+int(floor(float(pos[0]+radius[0]-orig[0])/dx[0]+0.5))
                    j0 = int(ceil(float(pos[1]-radius[1]-orig[1])/dx[1]-0.5))
                    j1 = 1+int(floor(float(pos[1]+radius[1]-orig[1])/dx[1]+0.5))

                # check horizontal boundaries
                if i0<0 or i1>gridPtNb[0] or j0<0 or j1>gridPtNb[1]:
                    warnings.add(
                        "Draw point boundary outside grid horizontally.")
                    i0 = max(i0, 0)
                    i1 = min(i1, gridPtNb[0])
                    j0 = max(j0, 0)
                    j1 = min(j1, gridPtNb[1])

                # set points inside bounds to 1
                for i in range(i0, i1):
                    for j in range(j0, j1):
                        offs = i*st[0] + j*st[1]
                        result[offs+k0*st[2]:offs+k1*st[2]:st[2]] = [1]*(k1-k0)

        # write warnings only once in the end
        for txt in warnings:
            msg("WARNING: %s Result has been trimmed to the grid." % txt)

        return result


def writeCavesimGradefileHomogeneous(
        filename=None,
        density=2700,
        initPorosity=0.0,
        maxPorosity=0.13,
        particleDiameter=0.5,
        friction=32.0,
        cavePeriod=99999,       # never switch to cave automatically
        geoDomain=1,
        box=[[],[]],
        versionCavesim="v2",
        grades=[1.0, 0.0, 0.0, 0.0, 0.0],
        cavingRate=1,
        ):
    """Write Cavesim input file GRADE.TXT for homogeneous cave properties,
    i.e. just one block in the block model, one line in the grade file.

    @param filename: defaults to GRADE.txt for v2 and CAVESIM_BM.txt for v4
    @param density: ... in SI units kg/m^3, will be converted to tonnes/m^3 for
        Cavesim
    @param initPorosity: initial porosity
    @param maxPorosity: porosity of moved material (initPorosity=0.0 and
        maxPorosity=0.13 result in ~ 20% swell)
    @param particleDiameter: always use 0.5
    @param friction: friction zero means 30deg
    @param cavePeriod: Cavesim step in which to switch (this particular)
        block to cave.
        For coupling with Abaqus use 99999 (= never switch).
    @param geoDomain:  number of the material domain
    @param box: location of the model box, should correspond to the
        box defined in the CAVESIM.inp file. A L{bae.misc_01.BoundingBox}
        or equivalent pair of points.
    @param versionCavesim: "v2" or "v4"
    @param grades: (v4 only) grades for 5 different types of ore
    @param cavingRate: (v4 only) see Cavesim manual, last column in "grade"
        file CAVESIM_BM.txt
    """
    if filename is None:
        if versionCavesim=="v2":
            filename = "GRADE.txt"
        else:
            filename = "CAVESIM_BM.txt"
    fout = newFile(filename)

    # this is how Cavesim wants density [t/m3]
    density = density / 1000.0

    box = BoundingBox(box)   # make sure it's really a BoundingBox
    centre = box.getCentroid()
    span = box.getSpan()

    if versionCavesim=="v2":
        fout.write(
            "%6.0f\t%6.0f\t%6.0f\t%d\t%d\t%d\t%4.2f\t%4.2f\t%4.2f\t"
            "%d\t0\t0\t%4.2f\t%f\t%d\n" % (
                centre[0], centre[1], centre[2],
                span[0], span[1], span[2],
                density, maxPorosity, initPorosity, geoDomain,
                particleDiameter, friction, cavePeriod))
    elif versionCavesim=="v4":
        fout.write("\t".join(
            "\t".join("%6.0f" % x for x in centre),
            "\t".join("%d" % x for x in span),
            "\t".join("%4.2f"%x for x in (density, maxPorosity, initPorosity)),
            "%d" % geoDomain,
            "\t".join(str(x) for x in grades),
            "%4.2f" % particleDiameter,
            "%f" % friction,
            "%d" % cavePeriod,
            )+"\n")
    else:
        raise NotImplementedError(
            "Cavesim version %s not implemented in"
            " writeCavesimGradefileHomogeneous." % versionCavesim)
    fout.close()


def updateGradesFromMaterial(fieldData):
    """derive material-grades from matCode data/column

    Creates grade_01, grade_02, ... fields in L{bae.field_02.FieldsCollection}
    object fieldData according to the matCode field.
    """
    matCode = fieldData['matCode'].array
    matIds = np.unique(matCode)
    if not np.equal(np.mod(matIds, 1), 0).all():
        raise ValueError('MatCodes have to be int(like) values')

    if max(matIds) >= 100.:
        raise ValueError('Only matCodes up to 99 are allowed')

    for matId in sorted(matIds):
        name = 'grade_%02d' % matId
        fieldInput[name] = (matCode == matId).astype(int)


def writeCavesimGradefile(
        fieldData,
        gradeFieldNames=None,
        outputPath=None,
        versionCavesim="v5",
        **kwargs):
    """Write Cavesim input file GRADE.TXT or CAVESIM_BM.txt for inhomogeneous
    cave properties.
    The data is supplied by a (potentially inhomogeneous) field or constant
    values.

    @param fieldData: a L{bae.field_02.FieldsCollection} object with fields
        to be written to the Cavesim input "grade" file. For field names and
        values see the keyword arguments to this function. Field names for
        the grades fields can be specified by the gradeFieldNames argument.

        Note that density must be given in SI unit kg/m^3 with will be
        converted to the Cavesim unit tonnes per m^3.
    @param gradeFieldNames: list of keys in fieldData identifying the grade
        fields. Defaults to "grade_01", "grade_02", ...

    @param outputPath: defaults to GRADE.txt for v2 and CAVESIM_BM.txt for v4
        and v5
    @param versionCavesim: "v2", "v4" or "v5"

    @kwarg density: ... in SI units kg/m^3, will be converted to tonnes/m^3 for
        Cavesim
    @kwarg initPorosity: initial porosity
    @kwarg maxPorosity: porosity of moved material (initPorosity=0.0 and
        maxPorosity=0.13 result in ~ 20% swell)
    @kwarg matCode: number of the material domain
    @kwarg grade_01: density ratio for material 1; between 0 and 1
    @kwarg grade_02: ... for material 2, ...up to grade_10 (for Cavesim v5)
    @kwarg particleDiameter: always use 0.5
    @kwarg friction: friction zero means 30deg
    @kwarg cavePeriod: Cavesim step in which to switch (this particular)
        block to cave.
        For coupling with Abaqus use 99999 (= never switch).
    @kwarg UCS: not relevant?
    @kwarg fines: not relevant?
    @kwarg mrmr: not relevant?
    @kwarg caveAngle: not relevant?
    """

    if not pd:
        raise NotImplementedError(
            "Current implementation of writeCavesimGradefile relies on Pandas.")

    # preprocess arguments
    if versionCavesim=="v5":
        nMaxGrades = 10

        # caveSim-header used for naming in pandas dataframe
        # it will not appear in Cavesim-Blockmodel, but can be used for
        # debugging
        header = ['X', 'Y', 'Z', 'dX', 'dY', 'dZ',
                  'density', 'maxPorosity', 'initPorosity', 'matCode'] \
                  + ['grade_%02d' % (i+1) for i in range(nMaxGrades)] \
                  + ['partDiameter', 'frictionAngle', 'cavePeriod', 'UCS',
                     'fines', 'mrmr', 'caveAngle']

        defaults = dict(
            density=2700.0,
            initPorosity=0.0,
            maxPorosity=0.13,
            matCode=1,
            grade_01=1.0,
            particleDiameter=0.5,
            frictionAngle=42.0,
            cavePeriod=99999,       # never switch to cave automatically
            UCS=999,
            fines=999,
            mrmr=38,
            caveAngle=90,
            )
    else:
        raise NotImplementedError(
            "Currently writeCavesimGradefile is implemented for v5 only.")

    # prepare output
    cavesim_bm = pd.DataFrame(
        data=np.zeros((len(fieldInput.topo), len(header))),
        columns=header)

    # centre coordinates
    xyz = fieldInput.topo.getPoints()
    for ii, key in enumerate(['X', 'Y', 'Z']):
        cavesim_bm[key] = xyz[:,ii].ravel()

    # cell sizes / spacing
    spacing = list(fieldInput.topo.spacing)
    dxyz= np.ones(xyz.shape) * spacing
    for ii, key in enumerate(['dX', 'dY', 'dZ']):
        cavesim_bm[key] = dxyz[:,ii].ravel()

    # store grades
    # ... at the same time write CaveSimGradeIndex to gradeName-lookup
    gradeInfo = open(os.path.splitext(outputPath)[0] + '_gradeInfo.csv', 'w')
    gradeInfo.write('#CSgrade, gradeName\n')

    if gradeFieldNames is None:
        gradeFieldNames = ["grade_%02d" % (i+1) for i in range(nMaxGrades)]

    for ii, fieldName in enumerate(gradeFieldNames):
        if ii>=nMaxGrades:
            msg("WARNING: There are %d grade fields requested for output."
                " Cavesim version %s only accepts %d. Ignoring those extra"
                " fields." % (len(gradeFieldNames), versionCavesim, nMaxGrades))
            break

        gradeInfo.write('grade_%02d, %s\n' % (ii+1, fieldName))

        try:
            inpField = fieldInput[fieldName].array
        except KeyError:
            pass
        else:
            cavesim_bm['grade_%02d' % (ii+1)] = inpField

    gradeInfo.close()

    # take other fields from fieldData
    # all other fields from header after x, y, z, dx, dy, dz; from nb seven
    otherFieldNames = [fn for fn in header[6:]
                       if fn not in set(gradeFieldNames)
                       and fn not in kwargs]
    for fieldName in otherFieldNames:
        val = fieldInput[fieldName].array
        if fieldName=="density":
            val = val/1000.0
        cavesim_bm[fieldName] = val
        # remove current, so we can assign the remaining defaults in the end
        del defaults[fieldName]

    # take constant argument values
    for fieldName, val in kwargs.iteritems():
        if fieldName=="density":
            val = val/1000.0
        cavesim_bm[fieldName] = val
        del defaults[fieldName]

    # assign remaining defaults in the end
    for fieldName, val in defaults.iteritems():
        if fieldName=="density":
            val = val/1000.0
        cavesim_bm[fieldName] = val

    # writing results
    cavesim_bm.to_csv(
        outputPath,
        index=False,    # turn off index as 1st coloumn
        header=False,   # turn of header
        sep='\t',       # Fred says CS gets confused by ',' sep
        float_format='%.5f',  # CS struggles with long floats or exp
        )
    msg("Wrote Cavesim blockmodel data to %s" % outputPath)


decomposeDateStringRex = re.compile(r"(\D*?)(\d+)(?:(\D+?)(\d+))?")
def decomposeDateString(dateString):
    """evaluate date string

    Usage:
     >>> prefYr, year, prefMon, month = decomposeDateString(dateString)

    @returns: a (prefixYr, year, prefixMon, month)-tuple
    prefixMon and month will be None in cases like "Y22".
    Note that year and month are strings and might require conversion to
    integers.
    Note that the month -if present at all- does not necessarily need
    to be a month it may also be the quarter. Look at prefixMon.
    Note that prefixMon may contain something like "_M", i.e. not just
    the "M" for month.

    @raises ValueError: if dateString cannot be interpreted
    """

    res = decomposeDateStringRex.match(dateString)
    if not(res):
        raise ValueError("Date string %s not valid." % dateString)
    if len(res.groups())==2:
        prefYr, year = res.groups()
        return prefYr, year, None, None
    else:
        prefYr, year, prefMon, month = res.groups()
        return prefYr, year, prefMon, month


class CaveSimSimple(object):
    r"""Simplified CaveSim replacement with pure vertical movement.

    Basic usage:
     >>> from bae.drawschedule_01 import CaveSimSimple
     >>> cavesim = CaveSimSimple(grid)
     >>> cavesim.initRockFracFromGradeFile(gradeFileName)
     >>> for stepNb in ...:
     >>>     cavesim.processMiningStep(stepNb, caveableNew, drawEvents)

    Basic usage without grade file:
     >>> from bae.drawschedule_01 import CaveSimSimple
     >>> cavesim = CaveSimSimple(grid)
     >>> cavesim.density = 2.8  # in tonnes/m^3
     >>> cavesim.porosityFinal = 0.1
     >>> for stepNb in ...:
     >>>     cavesim.processMiningStep(stepNb, caveableNew, drawEvents)

    Usage (extended example):
     >>> from bae.drawschedule_01 import DrawSchedule, CaveSimSimple
     >>> from bae.vtk_02 import FieldsOnStructPointGrid as Vtk
     >>>
     >>> # get grid and new caveable zone from XFER1.vtk
     >>> abqInp = Vtk().fromVtk("XFER1.vtk")
     >>> cavesim = CaveSimSimple(grid=abqInp.mesh)
     >>> caveableNew = abqInp.data.itervalues().next()
     >>>
     >>> # apply data from GRADE.txt / CAVESIM_BM.txt
     >>> cavesim.initRockFracFromGradeFile("CAVESIM_BM.txt")
     >>>
     >>> # get drawEvents from CAVESIM.inp
     >>> drawschedule = DrawSchedule().readCavesimInp("drawschedule.info.txt")
     >>> dpHalfWidth = [drawschedule.xwid, drawschedule.ywid, drawschedule.zwid]
     >>> stepNb = 1
     >>> iDrawEvent = 0
     >>> drawEvents = list()
     >>> while drawschedule[iDrawEvent][1] <= stepNb:
     >>>     pos, step, tonnes = drawschedule[iDrawEvent]
     >>>     drawEvents.append((pos, dpHalfWidth, tonnes))
     >>>     iDrawEvent += 1
     >>>
     >>> # process mining step
     >>> fldDraw, fldVel, tonnesScheduled, tonnesDrawn = \
     >>>     cavesim.processMiningStep(stepNb, caveableNew, drawEvents)
     >>>
     >>> # Cavesim output XNC1.txt
     >>> outputName = os.path.join(outputDir, "XNC%d.txt" % stepNb)
     >>> output = open(tempOutputPath, "w")
     >>> for idx, (pt, rockFr, c, vel) in enumerate(izip(
     >>>         cavesim.grid.getPointsIter(), cavesim.rockFrac,
     >>>         cavesim.caveable, fldVel)):
     >>>     if c == 0:
     >>>         val = "0.0 0.0 0.0 0.0"
     >>>     else:
     >>>         if rockFr < 1E-5:
     >>>             val = "0.0 0.0 0.0 999999.0"
     >>>         else:
     >>>             val = "0.0 0.0 %.1f %.1f" % (-vel, vel)
     >>>     output.write("%d %.1f %.1f %.1f %s\n"
     >>>                  % (stepNb, pt[0], pt[1], pt[2],val))
     >>> del output
     >>> os.rename(tempOutputPath, outputName)
     >>> msg("Wrote result to %s" % outputName)
     >>>
     >>> # write diagnostic output
     >>> outputName = "results_%03d.vtk" % stepNb
     >>> output = Vtk(mesh=cavesim.grid)
     >>> output.updateFields(
     >>>     cavesim.rockFrac, fldDraw, cavesim.caveable, fldVel)
     >>> output.toVtk(outputName)
     >>> msg("Wrote diagnostic output to %s" % outputName)

    @ivar grid: L{MeshStructuredPoints<bae.mesh_01.MeshStructuredPoints>}-object
    @ivar rockFrac: grid field for fraction of rock in grid cell. It's a
       L{Field<bae.field_01.Field>}-object and can be easily appended to a
       L{FieldsOnStructPointGrid<bae.vtk_02.FieldsOnStructPointGrid>}-object,
       i.e. exported to a vtk-file.
    @ivar caveable: grid field of caveable region for last step.
       0=intact rock (not caveable), 1=caveable rock, 2=air. It's an ordinary
       list.
    @ivar density: in tonnes / m^3; constant in space, just a float scalar.
    @ivar porosityFinal: porosity of moved material; constant in space, just a
       float scalar.
    @ivar cellVol: grid cell volume (just an abbreviation)
    @ivar areaColumn: area of the horizontal cross section of one column (just
       an abbreviation)
    """

    def __init__(self, grid):
        """
        Initializes the following values:
         - self.L{rockFrac} ... constant 1.0 (no porosity)
         - self.L{caveable} ... constant 0 (nothing caveable)
         - self.L{density} ... 2.7
         - self.L{porosityFinal} ... 0.0 (no swell)
         - self.L{cellVol} and self.L{areaColumn}
        """
        self.grid = grid

        # check: this algorithm does not work if grid.strides[0] != 1
        if grid.strides[0]!=1:
            raise ValueError(
                "This Program does not work if grid.strides[0]!=1!")

        # init field: fraction of rock (vs void) in each cell
        FldRockFrac = Field.classFromPosType(
            fieldName="rockFrac", position="structPt", dataType="scalar")
        self.rockFrac = FldRockFrac([1.0]*len(self.grid))

        # old caveable: nothing in the first step
        self.caveable = [0]*len(grid)

        # some default values
        self.density = 2.7
        self.porosityFinal = 0.0

        # abbreviations
        self.cellVol = grid.spacing[0]*grid.spacing[1]*grid.spacing[2]
        self.areaColumn = grid.spacing[0]*grid.spacing[1]

        # temporary variables initialized only once (for speed)
        self.rockTonnesInCol = [None]*(grid.gridPtNb[0]*grid.gridPtNb[1])
        self.caveableFrom = [None]*(grid.gridPtNb[0]*grid.gridPtNb[1])
        self.caveableToPO = [None]*(grid.gridPtNb[0]*grid.gridPtNb[1])

    def initRockFracFromGradeFile(self, filename="CAVESIM_BM.txt"):
        """Initialize rock content fraction self.L{rockFrac} from the initial
        porosity field of the CaveSim-grade-file. Also updates self.L{density}
        and self.L{porosityFinal}.

        Currently ignores varying density and varying final porosity. But
        writes a warning if those values are varying.

        @param filename: of the CaveSim-GRADE file (usually GRADE.txt for
            CaveSim v2 and CAVESIM_BM.txt for v4) The file may be tab, space or
            comma separated. The following columns will be recognized:
             - nb 1,2,3 as centroid of the current block
             - nb 4,5,6 as width of the block
             - nb 7 as density in t/m^3
             - nb 8 "PMAX" as maximum (i.e. swelled, final) porosity
             - nb 9 "PINI" as initial porosity.

            Can also be an open file object. (Or StringIO...)

        @returns: number of lines being read.
        """
        if isinstance(filename, basestring):
            csvfile = open(filename, 'rb')
        else:
            csvfile= filename

        # detect format (tab, space or comma separated)
        dialect = csv.Sniffer().sniff(csvfile.read(1024))
        csvfile.seek(0)
        tab = csv.reader(csvfile, dialect)

        # initialize some containers for diagnostics
        otherPorMax = set()
        porIniValues = set()
        otherDensities = set()
        initMissing = [1]*len(self.grid)

        # process file contents ...
        density = None
        porosityFinal = None
        for cnt, line in enumerate(tab):
            # ignore blank lines
            if not line:
                continue

            line = map(float, line)  # convert to float

            # determine block / box
            centroid = line[:3]
            width = line[3:6]
            origin = vector_minus(centroid, vector_scale(width, 0.5))

            # determine density and porosity values
            if density is None:
                density = line[6]
            elif density != line[6]:
                otherDensities.add(density)
            if porosityFinal is None:
                porosityFinal = line[7]
            elif porosityFinal != line[7]:
                otherPorMax.add(porosityFinal)
            currentRockFrac = (1.0-line[8])
            porIniValues.add(line[8])

            # assign initial porosity
            ptIds = self.grid.getIdsInBox(
                [origin, vector_plus(origin, width)])
            for idx in ptIds:
                self.rockFrac[idx] = currentRockFrac
                initMissing[idx] = 0

        # diagnostic output
        msg("Found values for initial porosity from %g to %g. Final porosity:"
            " %g" % (min(porIniValues), max(porIniValues), porosityFinal))

        if sum(initMissing):
            msg("WARNING: %d grid points have not been initialized."
                % sum(initMissing))

        if otherDensities:
            msg("WARNING: Density varying from %g to %g. Assuming"
                " first value of %g."
                % (min(otherDensities),
                   max(otherDensities), porosityFinal))

        if otherPorMax:
            msg("WARNING: Maximum porosity varying from %g to %g. Assuming"
                " first value of %g."
                % (min(otherPorMax), max(otherPorMax), porosityFinal))

        if density is not None:
            self.density = density
        if porosityFinal is not None:
            self.porosityFinal = porosityFinal
        return cnt+1

    def pullFromColumn(
            self, draw, caveableFrom, caveableToPO, fldVel, idx0, didx):
        """Pull from column, with swell.
        Will be called by self.L{processMiningStep}. Updates self.L{rockFrac}
        and the velocity field being supplied as the fldVel argument.

        It must be guaranteed beforehand that the scheduled draw (as supplied
        by the corresponding argument) can actually be drawn from the current
        column. I.e. draw <= sum(rockFrac[idx(k)]*cellVol*density ...
        ... for k in range(caveableFrom, caveableToPO))

        @param draw: tonnes to be drawn from column
        @param caveableFrom: caveable region, vertical start index
        @param caveableToPO: caveable region, vertical end index Plus One
        @param fldVel: grid field for (vertical) velocity, will be updated by
            this method
        @param idx0: index in rockFrac and fldVel of the lowest (min z) grid
            cell
        @param didx: index offset going to the cell above
        """

        # abbreviation
        cellHeight = self.grid.spacing[2]

        # convert draw from unit tonnes to unit self.cellVol*self.density
        # (being the mass of a cell filled with rock of zero porosity)
        draw /= self.cellVol*self.density

        # loop bottom up and pull draw
        idx = idx0 + caveableFrom*didx
        idxEnd = idx0 + caveableToPO*didx
        while draw>1E-3 and idx<idxEnd:
            if draw>self.rockFrac[idx]:
                draw -= self.rockFrac[idx]
                self.rockFrac[idx] = 0.0
                idx += didx
            else:
                self.rockFrac[idx] -= draw
                draw = 0.0
                break  # actually not necessary, loop ends anyway

        # loop again to let fall rock
        idxSource = idx  # Where does the rock come from? All empty below idx.
        idx = idx0 + caveableFrom*didx  # where does it go?
        rockFracSwelled = 1.0-self.porosityFinal
        # msg("pullFromColumn: after draw before rock fall: idx %d,"
        #     " idxSource %d, idxEnd %d" % (idx, idxSource, idxEnd),
        #     debugLevel=10)

        # main loop primarily over idx and secondarily over idxSource
        # possibly moving rock from idxSource to idx
        while (idx<idxEnd and idxSource<idxEnd):

            # accumulating the flux into the current cell
            # in unit self.cellVol*self.density*cellHeight/timeStep
            # (it's actually a momentum: mass * velocity, in special units)
            influx = 0.0

            # start value for idxSource: pull from cell above (and higher)
            idxSource = max(idxSource, idx+didx)

            # pull completely empty the cell(s) at idxSource, moving rock to idx
            # loop over idxSource (moving up)
            while (idxSource<idxEnd and
                   # possible intake of idx-cell > rock in idxSource-cell
                   rockFracSwelled-self.rockFrac[idx]>self.rockFrac[idxSource]):

                # fraction of cell contents moving:
                rockFracMoving = self.rockFrac[idxSource]

                # update velocity:
                # distance travelled by the rock = (idxSource-idx)/didx
                influx += rockFracMoving * (idxSource-idx)/didx

                # move the rock
                self.rockFrac[idx] += rockFracMoving
                self.rockFrac[idxSource] = 0.0
                # msg("pulling empty idxSource %d to idx %d" %(idxSource, idx),
                #     debugLevel=10)
                idxSource += didx

            # fill up current cell (idx) with part of cell at idxSource
            # ... if idx-cell not yet full
            if idxSource<idxEnd and rockFracSwelled>self.rockFrac[idx]:

                # fraction of cell contents moving:
                rockFracMoving = rockFracSwelled-self.rockFrac[idx]

                # update velocity:
                # distance travelled by the rock = (idxSource-idx)/didx
                influx += rockFracMoving * (idxSource-idx)/didx

                self.rockFrac[idxSource] -= rockFracMoving
                self.rockFrac[idx] = rockFracSwelled
                # msg("filled up cell idx %d with rock from idxSource %d"
                #     % (idx, idxSource), debugLevel=10)

            # update velocity
            # influx in mass*velocity is being divided by total mass
            # ... mass in unit self.cellVol*self.density
            # then the resulting velocity is converted from unit cellHeight
            # to metres:
            if self.rockFrac[idx]:
                fldVel[idx] += cellHeight * influx/self.rockFrac[idx]

            # next cell
            idx += didx

    def processMiningStep(self, stepNb, caveableNew, drawEvents):
        """Process one mining step.

        @param stepNb: for diagnostic output
        @param caveableNew: grid field of the new caveable region for this
           step. 0=intact rock (not caveable), 1=caveable rock, 2=air
        @param drawEvents: list of (pos, dpHalfWidth, tonnes)-tuples

        @returns: (fldDraw, fldVel, tonnesScheduled, tonnesDrawn)
           fldDraw: grid field for accumulated production
           fldVel: grid field for (vertical) velocity
           tonnesScheduled: ... according to drawEvents-argument
           tonnesDrawn: actually pulled
           The returned grid fields are L{Field<bae.field_01.Field>}-objects
           and can be easily appended to a
           L{FieldsOnStructPointGrid<bae.vtk_02.FieldsOnStructPointGrid>}
           -object, i.e. exported to a vtk-file.
        """
        cellMass = self.cellVol*self.density  # abbreviation

        FldDraw = Field.classFromPosType(
            fieldName="draw", position="structPt", dataType="scalar")
        fldDraw = FldDraw([0.0]*len(self.grid))

        FldVelocity = Field.classFromPosType(
            fieldName="velocity", position="structPt", dataType="scalar")
        fldVel = FldVelocity([0.0]*len(self.grid))

        # find columns with new caveable rock
        newCaveableRock = list()
        for i in range(self.grid.gridPtNb[0]):
            for j in range(self.grid.gridPtNb[1]):
                # loop from top down
                ktop = None  # top index of new caveable column
                foundNewCaveable = False
                for k in range(self.grid.gridPtNb[2]-1, -1, -1):
                    idx = (i*self.grid.strides[0] + j*self.grid.strides[1]
                           + k*self.grid.strides[2])
                    if ktop is None:
                        if (0.5<caveableNew[idx]<1.5):
                            # was not caveable, now is caveable rock (not air)
                            ktop = k
                            foundNewCaveable = (self.caveable[idx]<0.5)
                    else:  # if ktop is not None:
                        if self.caveable[idx]<0.5:
                            foundNewCaveable = True
                        if caveableNew[idx]<0.5 and self.caveable[idx]<0.5:
                            # reached bottom of caveable (old or new!)
                            if foundNewCaveable:
                                newCaveableRock.append((i, j, k+1, ktop+1))
                            ktop = None
                            foundNewCaveable = False

                if ktop is not None and foundNewCaveable:
                    newCaveableRock.append((i, j, 0, ktop+1))


        msg("Found %d pieces, %d cells of new caveable rock."
            % (len(newCaveableRock), sum((x[3]-x[2] for x in newCaveableRock))))
        self.caveable = caveableNew
        del caveableNew, ktop

        # apply air: set rockFrac to 0 where XFER file is 2
        self.rockFrac = type(self.rockFrac)(
            int(a<2)*r for r, a in izip(self.rockFrac, self.caveable))
        msg("Applied air from Abaqus input. Now we have %d cells air (not from"
            " Abaqus alone)."
            % sum(1 for r in self.rockFrac if r==0.0), debugLevel=1)

        # rock fall before production from draw points
        for i, j, k0, k1 in newCaveableRock:
            # assuming strides[0] == 1 !!!!!
            flatIdx = i*self.grid.strides[0]+j*self.grid.strides[1]
            self.pullFromColumn(
                0.0, k0, k1, fldVel, flatIdx, self.grid.strides[2])

        # process draw
        tonnesScheduled = 0.0
        tonnesDrawn = 0.0
        for iDrawEvent, (pos, dpHalfWidth, tonnes) in enumerate(drawEvents):

            # process draw event
            tonnesScheduled += tonnes

            # region to draw from, according to draw schedule
            xmin = [x-dx for x, dx in izip(pos, dpHalfWidth)]
            ximin = [
                max(0, int(ceil(float(x-o)/s)))
                for x, o, s in izip(xmin, self.grid.origin, self.grid.spacing)]
            xmax = [x+dx for x, dx in izip(pos, dpHalfWidth)]
            ximax = [
                min(n-1, int(floor(float(x-o)/s)))
                for n, x, o, s in izip(
                    self.grid.gridPtNb,xmax,self.grid.origin,self.grid.spacing)]
            colsFlatIds = [  # assuming strides[0] == 1 !!!!!
                i*self.grid.strides[0]+j*self.grid.strides[1]
                for i in range(ximin[0], ximax[0]+1)
                for j in range(ximin[1], ximax[1]+1)]
            nbColsToDrawFrom = len(colsFlatIds)
            drawArea = nbColsToDrawFrom*self.areaColumn

            # commented out, see same block below
            #msg("Draw event %d for step %d, pull %.2f tonnes from %s."
            #    " Pulling from %d grid cells, an area of %g."
            #    % (iDrawEvent+1, stepNb, tonnes, ",".join(map(str, pos)),
            #       nbColsToDrawFrom, drawArea), debugLevel=2)

            # loop over columns, write max possible draw to rockTonnesInCol
            for flatIdx in colsFlatIds:
                alreadyInCave = False
                rockInCurrentCol = 0.0
                for k in range(ximin[2], self.grid.gridPtNb[2]):
                    idx = flatIdx + k*self.grid.strides[2]
                    if alreadyInCave and not self.caveable[idx]:
                        # reached the cave back, "PO" = Plus One
                        self.caveableToPO[flatIdx] = k
                        break

                    if not alreadyInCave:
                        # check if still/now in draw region
                        if k > ximax[2]:
                            # this column is not caveable at all
                            self.caveableFrom[flatIdx] = 0
                            self.caveableToPO[flatIdx] = 0
                            break

                        # check if caveable from current cell upwards
                        if self.caveable[idx]:
                            alreadyInCave = True
                            self.caveableFrom[flatIdx] = k
                            # assume top open, as long as we don't find
                            # ... the cave back
                            self.caveableToPO[flatIdx] = self.grid.gridPtNb[2]

                    if alreadyInCave:
                        rockInCurrentCol += self.rockFrac[idx]*cellMass

                # store max possible draw for current column
                self.rockTonnesInCol[flatIdx] = rockInCurrentCol
                del rockInCurrentCol

            # again loop over columns, smallest first, pull from columns
            supposedDrawPerCol = tonnes / nbColsToDrawFrom
            sortedColIdxTonnesRock = sorted(
                ((flatIdx, self.rockTonnesInCol[flatIdx])
                 for flatIdx in colsFlatIds),
                key=lambda (x,y):y)

            msg("Draw event %d for step %d, pull %.2f tonnes from %s."
                " Pulling from %d/%d grid cells, a total area of %g."
                % (iDrawEvent+1, stepNb, tonnes, ",".join(map(str, pos)),
                   sum(1 for (i,t) in sortedColIdxTonnesRock if t>0.0),
                   nbColsToDrawFrom, drawArea), debugLevel=2)

            for i, (flatIdx, rockInCurrentCol) in enumerate(
                    sortedColIdxTonnesRock):

                # decide how much to actually pull: currentColDraw
                # (ignoring anything less than 1kg)
                if (rockInCurrentCol-supposedDrawPerCol) > -1E-3:
                    # there is enough rock in this column
                    # ...to comply with supposedDrawPerCol
                    currentColDraw = supposedDrawPerCol
                else:
                    # column over-drawn
                    currentColDraw = rockInCurrentCol
                    if i<(nbColsToDrawFrom-1):
                        # if not in last column then distribute overhang
                        supposedDrawPerCol += (
                            (supposedDrawPerCol-rockInCurrentCol)
                            / (nbColsToDrawFrom-i-1))
                    else:
                        msg("WARNING: draw event %d in step %d could only"
                            " draw %g of the requested %g tonnes"
                            % (iDrawEvent+1, stepNb, rockInCurrentCol,
                               supposedDrawPerCol), debugLevel=1)

                # actually pull from column
                self.pullFromColumn(
                    currentColDraw,
                    self.caveableFrom[flatIdx], self.caveableToPO[flatIdx],
                    fldVel, flatIdx, self.grid.strides[2])

                # update fldDraw
                ii = flatIdx + self.caveableFrom[flatIdx]*self.grid.strides[2]
                fldDraw[ii] += currentColDraw

                tonnesDrawn += currentColDraw

        return (fldDraw, fldVel, tonnesScheduled, tonnesDrawn)
