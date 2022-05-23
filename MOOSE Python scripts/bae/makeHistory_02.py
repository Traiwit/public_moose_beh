r"""DEPRECATED, use L{bae.makeHistory_03}
It should be save to import the contents of this module with a wildcard as in:

>>> from bae.makeHistory_02 import *

Then you would also have the symbolic constants ES_... at hand.

See also the wiki http://cassowary:8080/TutorialBaeUtils/ModuleMakeHistory
"""


## ---------------------------------------------------------------------------
## version control
__version__ = '2.05'

_version_history_ = """\
2.01 GP new - created new from scratch
2.02 GP fixed - some errors
        changed - SDVINI ordered according to step number and extra comment
2.03 GP fixed - really ignore elsets not defined in the input file (no SDVINI)
        changed - DLOAD ordered according to step number (not yet finished)
2.03b FR changed - "STD_v2011" was added to the SDVINI output pattern        
2.04 GP write a warning if there are columns without headline in the csv file
2.05 FR adding the remark that updateDisjointElsets now takes a list of filenames
"""
## ---------------------------------------------------------------------------



import sys, itertools
from cStringIO import StringIO
from time import localtime, strftime

from bae.misc_01 import DictTableFromCsvFile, quietlog
from bae.future_01 import *
from bae.abq_model_02 import Model

class _LoadHistoryStep(object):
    r"""
    object that represents one time step in the load history

    usually you don't have to access this class directly.
    """

    def __init__(self, stepLength=1.0, stepName="",
                 elsetOffsetGravAmp=None,
                 elsetOffsetStatus=None,
                 timePointListOffset=None):
        self.stepLength = stepLength
        self.stepName = stepName

        self.elsetOffsetGravAmp = dict()
        if elsetOffsetGravAmp!=None:
            self.elsetOffsetGravAmp.update(elsetOffsetGravAmp)

        self.elsetOffsetStatus = dict()
        if elsetOffsetStatus!=None:
            self.elsetOffsetStatus.update(elsetOffsetStatus)

        self.timePointListOffset = defaultdict(list)
        if timePointListOffset!=None:
            for tpList, offset in zip(*timePointListOffset):
                self.timePointListOffset[tpList].append(offset)


# symbolic values for the element status value (ES = element status)
ES_delete, ES_intact, ES_startexc, ES_endexc, ES_fill = range(5)

class LoadHistory(object):
    r"""
    """

    def __init__(self, startTime=0.0, logfile=sys.stdout):
        r"""
        @param startTime: the first step in this LoadHistory object starts
          at startTime. This will be used as time offset for all amplitudes
          as well as for the timing values in the INITIAL CONDITIONS for the
          SDV variables.
        @param logfile: open file object to write diagnostic output to. If None
          don't write at all, defaults to sys.stdout.
        """
        self.startTime = startTime

        if logfile == None:
            self.logfile = quietlog
        else:
            self.logfile = logfile

        self.stepList = list()
        self.postTypeSets = dict()

    #-- get the data into the LoadHistory object

    defaultTPointsName = "TPOINTS"

    def addStep(self, stepLength=1.0,
                stepName="",
                elsetOffsetGravAmp=None,
                elsetOffsetStatus=None,
                timePointListOffset=True):
        r"""
        Insert a new step into the load history.

        Updates self.stepList, a list of _LoadHistoryStep objects.

        @param elsetOffsetGravAmp: A dict or tuples suitable for dict.update().
          The keys are elset names, the values are lists of (offset gravAmp)
          tuples. gravAmp is the new gravity DLOAD amplitude and offset is the
          delay when this values is to be applied relative to the start of this
          new step.

        @param elsetOffsetStatus: A dict or tuples suitable for dict.update().
          The keys are elset names, the values are lists of (offset, status)
          tuples. offset is the delay when this values is to be applied
          relative to the start of this new step. status is the new status, an
          integer meaning:

           - 1: intact rock mass (this is the value at the beginning of the
             analysis and will probably never appear in this argument)
           - 2: start of degradation due to excavation
           - 3: excavation ends
           - 4: refill
           - 0: deletion of the elements

           The values 1 to 4 correspond to the SDV9 internal state variable of
           the Levkovitch Reusch user material.

        @param timePointListOffset: In the simplest case if this is True, a non
          zero number, "yes" or "on" then the end of this load history
          step is put into the Abaqus *TIME POINTS list named
          LoadHistory.defaultTPointsName. If False it is not added to
          this list. This list is supposed to be used to determine
          field output to the odb.

          It may also be a tuple (tpName, offset). If tpName and
          offset both are lists they have to be of the same length and
          zip(tpName, offset) yields a list of (time points list name,
          time offset)-tuples. This would lead to several time points
          for this single load history step, one for each item in this
          list. The tpName value(s) state the name of the *TIME POINTS
          list and the time offset is the duration between the end of
          this load history step and the actual time point in the list
          (negative values refer to time points before the end of the
          step). If only one is a list this leads to several time
          points as well with the single value (not being a list) used 
          constantly for all of them. If both tpName and offset are
          single values that leads to only one time point for this single
          load history step.

          By this means you can define arbitrary time points relative
          to the step in arbitrary Abaqus *TIME POINTS lists.

        @returns: The new _LoadHistoryStep object
        """
        # TODO:
        # @param nsetBCOffset: not implemented yet. will be similar to the
        # former arguments to allow the application of BC to a certain nset.

        if (timePointListOffset==True
            or (isinstance(timePointListOffset,int) and timePointListOffset!=0)
            or (isinstance(timePointListOffset, float)
                and not (-1E-6 < timePointListOffset < 1E-6))
            or (isinstance(timePointListOffset, basestring)
                and timePointListOffset.lower() in ("yes", "on", str(1)))
            ):
            timePointListOffset = ([self.defaultTPointsName,], [0.0,])
        elif not(isinstance(timePointListOffset, (tuple, list))):
            timePointListOffset = None

        newStep = _LoadHistoryStep(stepLength=stepLength, stepName=stepName,
                                   elsetOffsetGravAmp=elsetOffsetGravAmp,
                                   elsetOffsetStatus=elsetOffsetStatus,
                                   timePointListOffset=timePointListOffset)

        self.stepList.append(newStep)

        return newStep

    @staticmethod
    def _evaluateString(value, evalNamespaceDict, description):
        r"""
        for internal use only:
        if the value is a string it is tried to evaluate that as a python
        expression with the evalNamespaceDict as local name space.

        @param description: description of the value for error messages.
        """
        if isinstance(value, basestring):
            try:
                value = eval(value,globals(), evalNamespaceDict)
            except Exception, exc:
                raise ValueError(
                    "Could not evaluate the expression '%s' for %s:\n%s"
                    % (value, description, exc.args[0]))
        return value


    def readFromCsv(self, loadhistoryCsv, colToGravAmp, colToStatus,
                    postTypeCols={},
                    stepLength="Steplength", stepName="Name",
                    timePointList=True,
                    extraLocals={}):
        r"""
        Reads the data for this loadhistory from the given csv file. Multiple
        calls to this function should append the data to the already existing.

        @param loadhistoryCsv: file name or open file containing the data

        @param colToGravAmp: A dictionary {column label:gravity amplitude list}
          Each dictionary key aka column label identifies the column(s) under
          which the elsets this timing applies to are listed. E.g. 'Pit',
          'Cave', 'Fill'.
          The corresponding value is the gravity amplitude list: A list of
          (time, amplitude) tuples. amplitude is the new amplitude value at the
          given time. time may be a float, in that case this same time value
          applies to all loadhistory steps. Or it may be a string evaluating to
          the time value. See the note about value evaluation.

        @param colToStatus: A dictionary {column label:status timing list}
          Each dictionary key aka column label identifies the column(s) under
          which the elsets this timing applies to are listed. E.g. 'Pit',
          'Cave', 'Fill'.
          The corresponding value is the status timing list: A list of
          (time, status) tuples. status is the new status value at the
          given time, one of ES_delete, ES_intact, ES_startexc, ES_endexc,
          ES_fill or the correspondig integer 0..4. time may be a float, in
          that case this same time value applies to all loadhistory steps. Or
          it may be a string evaluating to the time value. See the note about
          value evaluation below.

        @param postTypeCols: A dictionary {post set type:list of column labels}
          specifying elsets (identified by the column label) for post set
          groups e.g. {'Support':['DRIVES','FILL'], 'Cave':['Cave']}.

          Note: there will be an additional group containing all elsets that
          occur in any column that appears in either colToGravAmp, colToStatus
          or postTypeCols. See self.writeCaePostSets4abaqusMacros() for
          details.

        @param stepLength: Either the step length (a float), in this case this
          is the length of all steps.  Or it may be a string evaluating to the
          step length. See the note on value evaluation below.

        @param stepName: A string evaluating to the step name. See the note on
          value evaluation below. If you want to take it from a column supply
          that columns name. If you want the same name for all steps supply
          the name inclosed in quotation marks, i.e.:

          >>> stepName = "'OZM'"

        @param timePointList: This may be or evaluate to (see the note about
          value evaluation) True, a non zero number, "yes" or "on", in that
          case the corresponding load history step will be added to a
          time points list named LoadHistory.defaultTPointsName. If it
          evaluates to False, zero, "no" or "off" this step is not
          included in the list.

          It may also be a tuple (tpName, offset). If tpName and/or
          offset are strings it is tried to evaluate them (see note on
          value evaluation below) and then passed on as the
          timePointListOffset argument to self.addStep().

          This means you may have multiple columns in the csv file
          labeled say "tpName" and either one or an equal number of
          columns labeled "tpOffset" and then pass ("tpName",
          "tpOffset") as timePointList argument and it will create
          many time points in the Abaqus *TIME POINTS list named as
          the contents of the "tpName" columns with time offsets
          either all the same from the one "tpOffset" column or from
          the corresponding "tpOffset" column.

          tpName and offset may as well be lists by themselves (both
          or either of them) in that case the evaluation with the csv
          data is skipped and the tpName and offset vaues will be the
          same for all load history steps.

          tpName and offset may be single values (in the case of
          tpName that's the case if the evaluation with the csv data
          fails, in the case of offset this is the case if offset is a
          number and not a string). Then there will be one time point
          in the time points list tpName with offset being the
          duration between the end of this load history step and the
          actual time point in the list (negative values refer to time
          points before the end of the step).

        @param extraLocals: A namespace dictionary that contains additional
          variables that should be accessible during value evaluation, see
          below. Specify extraLocals=locals() in the argument list to make
          available all current variables.

        @note: Value evaluation
          Some strings (e.g. the stepLength argument if it's a string) are
          evaluated as python expression. In that expression all the values of
          the current row are accessible as variables named like the column.
          I.e. if you have a column labeled 'steptime' with a value 5 in the
          current row then the variable steptime evaluates to 5 and e.g. an
          expression '0.5*steptime' would be a valid expression for the
          stepLength argument yielding 2.5 in this case.

          Additionally you can access other varibales contained in the
          extraLocals dict. Specify extraLocals=locals() in the argument list
          to make available all current variables.

          Probably the most frequent use will be just to specify a string
          identifying the column in the csv file from which to read the value.

          Note that spaces in column names will be replaced by an equal number
          of underscores ("_") to make them valid python identifiers.

        @note: The order of columns with different names is not significant.
        """
        # initialize the data table form the csv file
        loadHistTable = DictTableFromCsvFile(
            loadhistoryCsv,replaceSpaceBy="_",
            logfile=self.logfile, warnIgnoreCols="message")
        if isinstance(loadhistoryCsv, basestring):
            loadhistoryCsvName = loadhistoryCsv
        else:
            try:
                loadhistoryCsvName = loadhistoryCsv.name
            except AttributeError:
                loadhistoryCsvName = "<generic csv input>"

        # initialize self.postTypeSets
        for postType in postTypeCols:
            if postType not in self.postTypeSets:
                # in case this is not the first call there already are some
                # steps in self.stepList. All new self.postTypeSets values must
                # be filled with as many empty sets of elset names as there are
                # already steps in the loadhistory
                self.postTypeSets[postType] = [frozenset(),]*len(self.stepList)

        # all elset columns
        elsetColumns = set(colToGravAmp)
        elsetColumns.update(colToStatus)
        for elsetNames in postTypeCols.itervalues():
            elsetColumns.update(elsetNames)

        # main loop over all lines in the loadhistory table csv file
        rowIterator = loadHistTable.rowIterColToList(
            listColumns=elsetColumns,removeNones=True)
        for rowCnt, row in enumerate(rowIterator):

            # the locals argument to the eval() call to evaluate non numeric
            # values in the stepLength argument and in the time offset fields
            # of colToGravAmp and colToStatus
            evalNamespaceDict = dict(extraLocals)
            evalNamespaceDict.update(row)

            # evaluate stepLength
            stepLengthVal = self._evaluateString(
                stepLength, evalNamespaceDict,
                "the step length in the %dth data row of the loadhistory table"
                " %s" % (rowCnt+1,  loadhistoryCsvName))

            # evaluate stepName
            stepNameVal = self._evaluateString(
                stepName, evalNamespaceDict,
                "the step name in the %dth data row of the loadhistory table"
                " %s" % (rowCnt+1,  loadhistoryCsvName))

            # fill dictionaries elsetOffsetGravAmp and elsetOffsetStatus
            # with data from the csv file
            elsetOffsetGravAmp = dict()
            elsetOffsetStatus = dict()
            for resultDict, colToDataDict in [
                (elsetOffsetGravAmp, colToGravAmp),
                (elsetOffsetStatus, colToStatus)]:

                for col, timeDataList in colToDataDict.iteritems():
                    try:
                        elsetNames = row[col]
                    except KeyError:
                        elsetNames = []
                    if len(elsetNames)==0:
                        continue

                    newTimeDataList = list()
                    for ti, data in timeDataList:
                        ti = self._evaluateString(
                            ti, evalNamespaceDict,
                            "the %dth data row in column %s of the loadhistory"
                            " table %s" % (rowCnt+1, col, loadhistoryCsvName))
                        newTimeDataList.append((ti, data))

                    for elsetName in elsetNames:
                        resultDict[elsetName] = newTimeDataList

            # evaluate timePointList
            timePointListVal = timePointList  # default, if nothing else works

            if isinstance(timePointList, (tuple, list)):
                # if its a tuple evaluate both items seperately
                (tpName, offset) = timePointList
                try:
                    tpName = self._evaluateString(
                        tpName, evalNamespaceDict, "")
                except ValueError:
                    pass
                try:
                    offset = self._evaluateString(
                        offset, evalNamespaceDict, "")
                except ValueError:
                    pass
                timePointListVal = (tpName, offset)
            else:
                # if its a string evaluate it, otherwise just pass it on
                try:
                    timePointListVal = self._evaluateString(
                        timePointList, evalNamespaceDict, "")
                except ValueError:
                    pass

            # add a new step to the load history
            self.addStep(stepLength=stepLengthVal, stepName=stepNameVal,
                         elsetOffsetGravAmp=elsetOffsetGravAmp,
                         elsetOffsetStatus=elsetOffsetStatus,
                         timePointListOffset=timePointListVal)

            # add the elset names to the corresponding post set groups
            postTypeList = list(postTypeCols)
            postTypeList.sort()
            for postType in postTypeList:
                newElsetsForStepInPostType = set()
                colsList = postTypeCols[postType]
                for colName in colsList:
                    try:
                        elsetNames = row[colName]
                    except KeyError:
                        elsetNames = []
                    newElsetsForStepInPostType.update(elsetNames)
                # add a new step to this self.postTypeSets item
                self.postTypeSets[postType].append(newElsetsForStepInPostType)


    def readFromCsvOldPitFillCaveDt(self, loadhistoryCsv):
        r"""
        Reads the data for this loadhistory from the given csv file (file name
        or open file).

        Expects columns labeled "Pit", "Cave", "Fill" (each column has to be
        labeled, the same label may appear multiple times, once or not at all)

        Timing: ...

        @Note: Using this function is discouraged. It will most probably be
        removed in the next version.
        """
        loadHistTable = DictTableFromCsvFile(loadhistoryCsv)

        elsetColumns = ["Pit", "Cave", "Fill"]
        for row in loadHistTable.rowIterColToList(listColumns=elsetColumns,
                                                  removeNones=True):
            stepLength = row["Time total"]

            # - t1 excavation starts: startTime + dt1
            # - t2 excavation ends: t1 + dt2
            # - t3 refill: Pit, Fill: t2 + dt3; Cave: never
            # - t4 element is deleted: Pit: t3 + dt4; Cave, Fill: never
            t1 = row["Time dt1"]
            t2 = t1 + row["Time dt2"]
            t3 = t2 + row["Time dt3"]
            t4 = t3 + row["Time dt4"]

            elsetOffsetGravAmp = dict()
            # Status: ES_delete, ES_intact, ES_startexc, ES_endexc, ES_fill
            elsetOffsetStatus = dict()

            for elsetName in row["Pit"]:
                elsetOffsetGravAmp[elsetName] = [(t1, 1.0), (t2, 0.0)]
                elsetOffsetStatus[elsetName] = [
                    (t1, ES_startexc), (t2, ES_endexc), (t3, ES_fill), (t4, ES_delete)]
            for elsetName in row["Fill"]:
                elsetOffsetGravAmp[elsetName] = [
                    (t1, 1.0), (t2, 0.0), (t3, 0.0), (t4, 1.0)]
                elsetOffsetStatus[elsetName] = [
                    (t1, ES_startexc), (t2, ES_endexc), (t3, ES_fill)]
            for elsetName in row["Cave"]:
                elsetOffsetGravAmp[elsetName] = [(t1, 1.0), (t2, 0.8)]
                elsetOffsetStatus[elsetName] = [
                    (t1, ES_startexc), (t2, ES_endexc)]
            
            self.addstep(stepLength, elsetOffsetGravAmp, elsetOffsetStatus)

    #------------------------------------------------------------------
    #-- get data out of the LoadHistory object

    def getGravAmps(self, applyStartOffset=True):
        r"""
        Returns a dictionary {elset name: list of (time, value) tuples}.

        The (time, value)-lists define the gravity dload amplitude for this
        elset.

        @param applyStartOffset: The first load history step starts at
          self.startTime. If applyStartOffset=False then the time values are
          relative to this start point, otherwise self.startTime is added to
          each time value.
        """
        result = defaultdict(list)

        # add all gravity amplitude items off all steps to the result dict
        if applyStartOffset:
            stepStart = self.startTime
        else:
            stepStart = 0.0
        for step in self.stepList:
            
            for elset, timesAmplList in step.elsetOffsetGravAmp.iteritems():
                for ti, amp in timesAmplList:
                    result[elset].append((ti+stepStart, amp))
            stepStart += step.stepLength

        # 1. sort all amplitudes (Suppose e.g. step 1: length 1.0, time offset
        # for a certain amplitude value for a certain elset is 10.0. Then in
        # the second step for the same elset there is an amplitude value with a
        # time offset of 5.0. Obviously the amplitude change associated to step
        # two happens before that of step 1.)
        # 2. make sure there are no two identical items
        for elsetName, ampList in result.iteritems():
            # task 1.
            ampList.sort(key=lambda x:x[0])
            # task 2. 
            lastTime = None
            lastAmp = None
            for ti, amp in list(ampList): # make a copy so we can modify it
                if lastTime==ti and lastAmp==amp:
                    ampList.remove((lastTime, lastAmp))
                lastTime = ti
                lastAmp = amp

        return result


    def getStrGravAmpDef(self):
        r"""
        Return a string containing the *amplitude commands needed to apply the
        gravity loads in an abaqus input file.
        """
        model = Model(logfile=None)
        elsetAmpList = self.getGravAmps(applyStartOffset=False)
        for elsetName, ampList in elsetAmpList.iteritems():
            ampName = "GRAVAMP_%s" % elsetName
            model.amplitude.updateAmp(
            ampName, ampList, smooth=None, time="STEP TIME")
        output = StringIO()
        model.write(output, "AMPLITUDE")
        return output.getvalue()


    def getStrGravDload(self):
        r"""
        Return a string containing the *dload commands needed to apply the
        gravity loads in an abaqus input file.
        """
        elsetAmpList = self.getGravAmps()
        elsetList = elsetAmpList.keys()
        elsetList.sort()
        output = StringIO()
        for elsetName in elsetList:
            ampList = elsetAmpList[elsetName]
            ampName = "GRAVAMP_%s" % elsetName
            output.write(
                "*DLOAD, OP=NEW, AMPLITUDE=%(ampName)s\n"
                "%(elsetName)s, GRAV, 9.81, 0., 0., -1.\n"
                % locals())
        return output.getvalue()


    # default values for the time t1 .. t4 in the INITIAL CONDITIONS command
    defaultTimes = [100000.0, 100001.0, 100002.0, 100003.0]

    # A dictionary { umat version : data line template for the
    # *INITIAL CONDITIONS command}
    umatVersionToSdvIniString = {
        "STD_v2008" : (
            "%(elsetName)s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0\n"
            "0.0, 1.0, 1.0, %(t1)6.3f, %(t2)6.3f, %(t3)6.3f, %(t4)6.3f\n"),
        "STD_v2008_wFredsExtraState" : (
            "%(elsetName)s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0\n"
            "0.0, 1.0, 1.0, %(t1)6.3f, %(t2)6.3f, %(t3)6.3f, %(t4)6.3f, 0.0\n"),
        "STD_v2010_wTempCouplingExtraState" : (
            "%(elsetName)s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0\n"
            "0.0, 1.0, 1.0, %(t1)6.3f, %(t2)6.3f, %(t3)6.3f, %(t4)6.3f, 0\n"
            "0.0, 0.0, 0.0\n"),
        "STD_v2011" : (
            "%(elsetName)s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0\n"
            " 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0\n"
            " 1.0,%(t1)6.3f, %(t2)6.3f, %(t3)6.3f, %(t4)6.3f\n"),       
        }

    def getElsetTimes(self):
        r"""
        Returns a dictionary {elset name: [t1, t2, t3, t4, stepIndx]}. The time
        values in the list have the following meaning:
         - t1 ... start excavation
         - t2 ... end excavation
         - t3 ... refill
         - t4 ... element delete

        stepIndx is the index in self.stepList where this elsets comesfrom
        """
        result = dict()
        statToInd = {ES_delete:3, ES_startexc:0, ES_endexc:1, ES_fill:2}

        # add all elsetOffsetStatus items of all steps to the result dict
        stepStart = self.startTime
        for stepCnt, step in enumerate(self.stepList):
            
            for elsetName, offsetStatus in step.elsetOffsetStatus.iteritems():
                timingList = list(self.defaultTimes)
                for ti, status in offsetStatus:
                    try:
                        index = statToInd[status]
                    except KeyError:
                        raise ValueError(
                            "A status value of %d is not implemented so far in"
                            " makeHistory_02.LoadHistory.getElsetTimes()."
                            % status)
                    timingList[index] = ti+stepStart
                if elsetName in result:
                    self.logfile.write(
                        "WARNING: In step %d a status timing for elset %s"
                        " occurs which is not the first time. This later"
                        " status timing will be ignored!"
                        % (stepCnt+1, elsetName))
                else:
                    result[elsetName] = timingList + [stepCnt]
            stepStart += step.stepLength

        return result


    def getStrSdvIni(self, umatVersion, elsetAll="ALL"):
        r"""
        Return a string containing the *INITIAL CONDITIONS,TYPE=SOLUTION
        command
        
        @param umatVersion: A string specifying the version of the umat
          subroutine: "LR_DP_v2009_12", ...

          Implemented so far: "STD_v2008", "STD_v2008_wFredsExtraState"

        @param elsetAll: All elements of this elset are initialized with a
          standard set of initial values before any of the elsets listed in
          self.stepList[###].elsetOffsetStatus are initialized. This elset
          should contain at least all elements not listed in elsetOffsetStatus
          that use the umat. It must not contain elements that are not assigned
          to this user material.
        """
        try:
            sdvIniString = self.umatVersionToSdvIniString[umatVersion]
        except KeyError:
            raise ValueError(
                "umatVersion %s not listed in"
                " makeHistory_02.LoadHistory.umatVersionToSdvIniString."
                % umatVersion)
        
        output = StringIO()
        output.write("*INITIAL CONDITIONS,TYPE=SOLUTION\n")
        elsetName = elsetAll
        t1,t2,t3,t4 = self.defaultTimes
        output.write(sdvIniString % locals())
        
        elsetTimes = self.getElsetTimes()
        stepElset = defaultdict(list)
        for elset, times in elsetTimes.iteritems():
            stepElset[times[4]].append(elset)
        for stepIndx, step in enumerate(self.stepList):
            output.write("**\n** %03d: %s\n"
                         % (stepIndx+1, step.stepName))
            elsetList = stepElset[stepIndx] # ok, since stepElset is a defaultdict!
            elsetList.sort()
            for elsetName in elsetList:
                t1,t2,t3,t4,stepIndx = elsetTimes[elsetName]
                output.write(sdvIniString % locals())
        return output.getvalue()


    def getTimePointsDict(self, addZeroPoints=False, applyStartOffset=False):
        r"""
        Returns a dictionary {time points list name: list of time points}.

        @param addZeroPoints: If True, add the start of the step to all time
          point lists.
        @param applyStartOffset: The first load history step starts at
          self.startTime. If applyStartOffset=False then the time values are
          relative to this start point, otherwise self.startTime is added to
          each time value.
        """
        result = defaultdict(list)

        if applyStartOffset:
            stepEnd = self.startTime
        else:
            stepEnd = 0.0

        # add all timePointListOffset items of all steps to the result dict
        for step in self.stepList:
            stepEnd += step.stepLength

            for tpName, offsets in step.timePointListOffset.iteritems():
                result[tpName].extend([ offs+stepEnd for offs in offsets ])

        # sort all time point lists
        for tlist in result.itervalues():
            if addZeroPoints:
                tlist.append(self.startTime)
            tlist.sort()

        return result


    def getStrTimePointsDef(self):
        r"""
        Return a string containing the *TIME POINTS commands for the field
        output requests in an abaqus input file.
        """
        model = Model()
        tpDict = self.getTimePointsDict(
            addZeroPoints=False, applyStartOffset=False)
        model.timepoints.update(tpDict)
        output = StringIO()
        model.write(output, "TIMEPOINTS")
        return output.getvalue()


    def getTotalLength(self):
        """
        return the sum of all stepLength values of all steps in self.stepList
        """
        return sum([step.stepLength for step in self.stepList])

    #------------------------------------------------------------------
    #-- other service

    def disjointSeqElsets(self, elsetsDict, updateStepList=False,
                          additionalElsets=[]):
        r"""
        Does what the boolElsets script did: Modify the given elsets in such a
        way that they are disjoint according to the order those elsets appear
        in self.stepList and after that in additionalElsets. Elsets that don't
        appear in the load history (or in the additionalElsets list) or that
        are empty after this procedure are not considered in the result. At the
        same time it removes empty elsets from the load history steps (given
        the corresponding optional argument).

        @param elsetsDict: A dictionary like abq_model_02.Model.elset. This is
          not touched. Instead a modified copy is returned

        @param updateStepList: If True remove all elsets not found in the
          resulting elsets dictionary from the self.stepList items. I.e. elsets
          that are empty after this process or have not been there from the
          beginning won't be assigned a gravity amplitude anymore nor will they
          appear in the initial conditions lists.

        @param additionalElsets: List of elset names that appear in elsetsDict
          and are supposed to appear in the output (if not empty). They are
          being made disjoint from the all others as well. This list acts as if
          there where additional steps at the end of the loadhistory, one for
          each elset in the order of this list. (This should be used for the
          GRAV-01 elsets.)

        @returns: A modified copy of the input elsetsDict with all the elsets
          made disjoint according to the order those elsets appear in
          self.stepList and after that in additionalElsets.

        @note: The order of elsets that occur in the same step for the first
          time is random. Neither does the column order in the underlying csv
          file influence the order of the elsets (if that was the data source)
          nor do different timing offsets for the different elsets. The latter
          may be considered a bug that could be fixed when necessary.
        """

        # check, if all additionalElsets are found in elsetsDict
        additionalNotFound = set(additionalElsets).difference(elsetsDict)
        if len(additionalNotFound):
            if len(additionalNotFound)<10:
                msgElsets = ":\n%s" % (", ".join(additionalNotFound))
            else:
                msgElsets = ""
            raise ValueError(
                "ERROR: %d elsets from the additionalElsets argument to"
                " disjointSeqElsets() not found in elsetsDict%s"
                % (len(additionalNotFound), msgElsets))

        # elsetNamesToStep is a dict {elset name: set of step indexes}. A step
        # index is the index in self.stepList of the step this elset occurs
        # in either the elsetOffsetGravAmp or elsetOffsetStatus argument.
        elsetNamesToStep = dict()
        # elsetsOrder: elset names in the order they appear in self.stepList
        elsetsOrder = list()
        for thisStepIndex, thisStep in enumerate(self.stepList):
            elsetsInCurrentStep = set(thisStep.elsetOffsetGravAmp).union(
                thisStep.elsetOffsetStatus)

            # remove elsets already mentioned earlier
            elsetDoubles = elsetsInCurrentStep.intersection(elsetNamesToStep)
            elsetsInCurrentStep.difference_update(elsetDoubles)
            if len(elsetDoubles):
                self.logfile.write(
                    "WARNING: %d elset(s) appear for the second time in step"
                    " %d of the loadhistory. This second occurence will be"
                    " disregarded!"
                        % (len(elsetDoubles), thisStepIndex+1))
            if updateStepList:
                for elsetName in elsetDoubles:
                    # remove this doubled elset from the steps
                    try:
                        del thisStep.elsetOffsetGravAmp[elsetName]
                    except KeyError:
                        pass
                    try:
                        del thisStep.elsetOffsetStatus[elsetName]
                    except KeyError:
                        pass

            # update elsetsOrder and elsetNamesToStep
            for elsetName in elsetsInCurrentStep:
                elsetsOrder.append(elsetName)
                try:
                    elsetNamesToStep[elsetName].add(thisStepIndex)
                except KeyError:
                    elsetNamesToStep[elsetName] = set((thisStepIndex,))
        elsetsOrder.extend(additionalElsets)

        # create: confirmedList, list of elset names that exist in the load
        # history and still contain elements after removing the intersection
        # with all previosly touched elsets.
        # disjointElsets: dictionary of disjoint elsets
        # allElements: set of all elements in any of the elsets
        # confirmedList = list()   not used at the moment
        disjointElsets = dict()
        allElements = set()
        emptyElsets = list()
        newElset = []
        for elsetName in elsetsOrder:
            removeThisElset = False
            try:
                newElset = elsetsDict[elsetName].difference(allElements)
            except KeyError:
                self.logfile.write(
                    "WARNING: Elset %s listed in the loadhistory table"
                    " not found. Will be ignored.\n" % elsetName)
                removeThisElset = True
            if len(newElset):
                allElements.update(newElset)
                disjointElsets[elsetName] = newElset
                # confirmedList.append(elsetName)
            else:
                emptyElsets.append(elsetName)
                removeThisElset = True

            if (removeThisElset and updateStepList
                and elsetName in elsetNamesToStep):
                # remove this empty elset from the steps
                for stepIdx in elsetNamesToStep[elsetName]:
                    thisStep = self.stepList[stepIdx]
                    try:
                        del thisStep.elsetOffsetGravAmp[elsetName]
                    except KeyError:
                        pass
                    try:
                        del thisStep.elsetOffsetStatus[elsetName]
                    except KeyError:
                        pass
                    
        if len(emptyElsets):
            self.logfile.write("WARNING: %d elsets have been removed because"
                               " they are empty"
                               % len(emptyElsets))
            self.logfile.write(":\n" + ", ".join(emptyElsets[:20]))
            if len(emptyElsets)>20:
                self.logfile.write("...\n")
            else:
                self.logfile.write("\n")

        # element sets that don't appear in the load history
        notUsedElsets = set(elsetsDict).difference(emptyElsets) \
            .difference(disjointElsets)
        if len(notUsedElsets):
            self.logfile.write(
                "WARNING: %d elsets do not appear in the load history. They"
                " won't appear as one of the resulting disjoint elsets either"
                % len(notUsedElsets))
            nbToD = 20
            self.logfile.write(":\n" + ", ".join(sorted(notUsedElsets)[:nbToD]))
            if len(notUsedElsets)>nbToD:
                self.logfile.write("... and %d more ...\n"
                                   % (len(notUsedElsets)>nbToD))
            else:
                self.logfile.write("\n")
                
        return disjointElsets

    def updateDisjointElsets(self, setsSequenceFileName, disjointSetsFileName,
                             additionalElsets=[]):
        """
        Removes empty elsets from the load history and (re)creates the Abaqus
        input file that contains the disjoint elsets if it appears to be
        necessary.

        This is basically a wrapper around self.disjointSeqElsets() that reads
        and writes the necessary files.

        @param setsSequenceFileName: file name or list of file names of the
          Abaqus input file(s) containing the elsets before they have been made
          disjoint
        @param disjointSetsFileName: file name of the Abaqus input file
          to be written on output containing the modified elsets
        @param additionalElsets: May be a list or a dictionary of additional
          elsets not in the load history steps that shall be made disjoint from
          the rest as if they appeared at the end of the load history. (This
          should be used for the GRAV-01 elsets.)

          If it is a list of elset names then those elsets must appear in the
          Abaqus input file mentioned as setsSequenceFileName.

          If the elsets don't appear in this Abaqus input file, you have to
          supply a dictionary {elset name: set of elements} e.g. like
          abq_model_02.Model.elset.
        """

        if type(setsSequenceFileName)==list:
            setsSequenceModel = Model(logfile=None)
            for setSequenceFile in setsSequenceFileName:
                setsSequenceModel.read(setSequenceFile)
        else:
            setsSequenceModel = Model(logfile=None).read(setsSequenceFileName)
        oldSeqElsets = setsSequenceModel.elset

        if isinstance(additionalElsets, dict):
            oldSeqElsets.update(additionalElsets)
            additionalElsetNames = additionalElsets.keys()
            additionalElsetNames.sort()
        elif isinstance(additionalElsets, list):
            additionalElsetNames = additionalElsets
        else:
            raise ValueError(
                "The additionalElsets argument must be either a dictionary or"
                " a list. (It's a %s.)" % type(additionalElsets))

        disjointSets = self.disjointSeqElsets(
            elsetsDict = oldSeqElsets,
            updateStepList = True,
            additionalElsets = additionalElsetNames)
        disjointSetsModel = Model(logfile=None)
        disjointSetsModel.elset.update(disjointSets)
        disjointSetsModel.write(disjointSetsFileName)


    def writeCaePostSets4abaqusMacros(self, outputFile, frameZeroName=None):
        r"""
        @param outputFile: An already open file or a file name for the python
          script file that will contain all post set definitions.
        @param frameZeroName: If it is a string it is appended to the list
          sequenceFrameNames in the output. Being the last item of this list
          it should be used as the name for the initial frame zero.
        """
        if isinstance(outputFile, basestring):
            output = open(outputFile, 'w')
        else:
            output = outputFile

        stringSeqIn = ("# " + "-"*75 + "\n"
                       "# Sequence %(setType)s\n\n"
                       "modelHistory%(setType)s = [\n")
        stringSeqOut = "  ]\n"

        #-- write the post elsets group for all elsets
        output.write(stringSeqIn % {'setType':''})
        allElsets = set()
        elsetsInThisFrame = list()
        for step in self.stepList:

            # elsets not found in any previous steps
            newElsets = set(step.elsetOffsetGravAmp) \
                .union(step.elsetOffsetStatus)
            newElsets.difference_update(allElsets)
            allElsets.update(newElsets)

            # add those to the elsets to assign to the next (possibly this)
            # output frame
            elsetsInThisFrame.extend(newElsets)

            # if this step corresponds to a frame: flush elsetsInThisFrame
            if len(step.timePointListOffset):
                output.write('    (')
                for setName in elsetsInThisFrame:
                    output.write("'PART-1-1."+setName+"', ")
                output.write('),\n' )
                elsetsInThisFrame = list()
        output.write(stringSeqOut)

        #-- write the post elsets groups for the groups in self.postTypeSets
        for postType, elsetLists in self.postTypeSets.iteritems():
            output.write(stringSeqIn % {'setType':postType})
            allElsets = set()
            elsetsInThisFrame = list()
            for step, elsetNames in itertools.izip(self.stepList, elsetLists):

                # elsets not found in any previous steps
                newElsets = elsetNames.difference(allElsets)
                allElsets.update(newElsets)

                # add those to the elsets to assign to the next (possibly this)
                # output frame
                newElsets = list(newElsets)
                newElsets.sort()
                elsetsInThisFrame.extend(newElsets)

                # if this step corresponds to a frame: flush elsetsInThisFrame
                if len(step.timePointListOffset):
                    output.write('    (')
                    for setName in elsetsInThisFrame:
                        output.write("'PART-1-1."+setName+"', ")
                    output.write('),\n' )
                    elsetsInThisFrame = list()
            output.write(stringSeqOut)

        #-- write the sequenceFrameNames
        output.write('\nsequenceFrameNames=[')
        for step in self.stepList:
            if len(step.timePointListOffset):
                # this step corresponds to a frame in the odb
                output.write("'%s',"%step.stepName)
        if isinstance(frameZeroName, basestring):
            output.write("'"+frameZeroName+"',")
        output.write(']\n')

        #-- fini
        return


#------------------------------------------------------------------------

class AbqMainInput(file):
    r"""
    This file class is meant to create the main Abaqus input file for an
    analysis.

    It is not meant to read, it will always overwrite an existing file of that
    same name.
    """
    def __init__(self, filename):
        file.__init__(self, filename, "w")

    # some default values
    stringSeparator="** "+"-"*75+"\n"

    def writeFileHeader(self):
        """
        write a default file header with time stamp, script name and module
        version information, an empty *heading command and a default *PREPRINT
        command.
        """
        self.write(self.stringSeparator)
        self.write("** created %s\n"
                   % strftime("%a, %d %b %Y %H:%M:%S", localtime()))
        self.write("** by %s - using makeHistory_02, version %s\n"
                   % (sys.argv[0], __version__))
        self.write(self.stringSeparator)
        self.write("*HEADING\n")
        self.write("*PREPRINT, ECHO=NO, MODEL=NO, HISTORY=NO, CONTACT=NO\n")
        self.write(self.stringSeparator)
