r'''makeHistory_01.py DEPRECATED, use L{bae.makeHistory_03}

used in a makeHistory script. This can look like this:
======================================================

>>> """makeHistory.py
>>> writes input file for VUMAT analysis
>>> 
>>> usage:
>>> ------
>>> 
>>> 1. copy sequence from the excel sheet
>>> 2. remove all <tabs>, all <'',>, all <,''> , all <''> in that order.
>>> 3. python makeHistory.py 
>>> 
>>> 
>>> makeHistory VUMAT 2.5
>>> 
>>> 1. Equilibrium step: 10s stress ramp
>>> 2. Sequence step
>>> """
>>> _version_ = '3.01'
>>> 
>>> from time import localtime, strftime
>>> from bae.makeHistory_01 import makeHistory
>>> 
>>> # ----------------------------------------------------------------------------
>>> # Format 
>>> sequence = [ \
>>>     {'Name':'Equilibrium',                 # comment line at begin of step
>>>      'Time':['total,dt1,dt2,dt3,dt4'],     # can overlap with next intervall
>>>      'Fill':['Setname1','Setname2','...'], # Set names
>>>      'Pit':['Setname1','Setname2','...'], # Set names
>>>      'Cave':['Setname1','Setname2','...'], # Set names
>>>      'Skip':'YES/NO'}]                     # optional, no ODB frame written.
>>> 
>>> sequence = [ \
>>> ... copy the contents of the excel sheet here ...
>>> ]
>>> 
>>> 
>>> # ----------------------------------------------------------------------------
>>> # other parameters
>>> 
>>> makeHistory(
>>>     sequence,
>>>     projectname='PRM805L_L1',
>>>     projectversion='01',
>>>     materialversion='LR01',
>>>     elsetsInputFileName='PRM805L_SEQ1',
>>>     elsetAllName='GRAV-01', 
>>> 
>>>     gravTimeFactor=0.5,
>>>     stepTimeEq = 10.,
>>> 
>>>     stringHeading="""\
>>> ** ----------------------------------------------------------------
>>> *HEADING
>>> *PREPRINT, ECHO=NO, MODEL=NO, HISTORY=NO, CONTACT=NO
>>> ** ----------------------------------------------------------------
>>> ** 
>>> ** Mesh
>>> ** 
>>> *INCLUDE,INPUT=RBC_SEQ7-meshQ1.inp
>>> *INCLUDE,INPUT=RBC_SEQ7-meshQ1_Regions.inp
>>> *INCLUDE,INPUT=PRM805L_SEQ1_IncElsets-%(version)s.inp
>>> **
>>> ** ----------------------------------------------------------------
>>> ** 
>>> ** Materials
>>> ** 
>>> *INCLUDE,INPUT=mat_RWAY_%(matName)s.inp
>>> **
>>> **
>>> ** SECTIONS
>>> **
>>> **
>>> *Solid Section, elset=Volcanics-s, material=VOLCANICS, controls=EC-1
>>> 1.,
>>> *Solid Section, elset=Volcanics, material=VOLCANICS, controls=EC-1
>>> 1.,
>>> *Solid Section, elset=Sediment-s, material=SEDIMENT, controls=EC-1
>>> 1.,
>>> *Solid Section, elset=Sediment, material=SEDIMENT, controls=EC-1
>>> 1.,
>>> *Solid Section, elset=Porphyry-s, material=VOLCANICS, controls=EC-1
>>> 1.,
>>> *Solid Section, elset=Porphyry, material=VOLCANICS, controls=EC-1
>>> 1.,
>>> *Solid Section, elset=Faults-s, material=FAULT, controls=EC-1
>>> 1.,
>>> *Solid Section, elset=Faults, material=FAULT, controls=EC-1
>>> 1.,
>>> *Solid Section, elset=Monzo-s, material=MONZONITE, controls=EC-1
>>> 1.,
>>> *Solid Section, elset=Monzo, material=MONZONITE, controls=EC-1
>>> 1.,
>>> **
>>> *Section Controls, name=EC-1, DISTORTION CONTROL=YES
>>> 1., 1., 1.
>>> **
>>> ** ----------------------------------------------------------------
>>> ** 
>>> ** Initial Conditions
>>> **
>>> """,
>>> #
>>> stringStepDef="""\
>>> **
>>> ** ----------------------------------------------------------------
>>> **
>>> ** Step
>>> **
>>> *STEP
>>> *DYNAMIC, EXPLICIT
>>> ,%(stepTimeEq)6.3f
>>> **
>>> *VARIABLE MASS SCALING, TYPE=BELOW MIN, DT=1.0e-2, NUMBER INTERVAL=10
>>> **
>>> *BULK VISCOSITY
>>> 0.1, 1.2
>>> **
>>> *BOUNDARY, TYPE=VEL
>>> FIX,  1,3, 0.
>>> ** 
>>> *DLOAD, OP=NEW, AMPLITUDE=ALL
>>> ALL, GRAV, 9.81, 0., 0., -1.
>>> ** 
>>> *OUTPUT, FIELD, NUMBER INTERVAL=1
>>> *Element Output
>>> S, LE, SDV, STATUS
>>> *Node Output
>>> U, V
>>> ** 
>>> *End Step
>>> ** ----------------------------------------------------------------
>>> **
>>> ** Step
>>> **
>>> *STEP
>>> *DYNAMIC, EXPLICIT
>>> ,%(stepTime)6.3f
>>> **
>>> *VARIABLE MASS SCALING, TYPE=BELOW MIN, DT=1.0e-2, NUMBER INTERVAL=50
>>> **
>>> *BULK VISCOSITY
>>> 0.1, 1.2
>>> ** 
>>> *BOUNDARY, OP=NEW
>>> FIX,  1,3
>>> ** 
>>> *DLOAD, OP=NEW, AMPLITUDE=GRAV
>>> GRAV-01, GRAV, 9.81, 0., 0., -1.
>>> ** 
>>> """
>>> )
'''

## ---------------------------------------------------------------------------
## version control
_version_ = '1.04'

_version_history_ = """\
1.01 GP new - extracted from boolElset1.08.py and some recent makeHistory.py
           (scriptVersion 2.5)
        changed - even faster boolElset algorithm
        fixed - correct baseList and emptysetList
1.02 GP changed - all %5.2f formats to %6.3f because small time steps produced
           error in amplitude commands.
1.03 GP changed to abq_model_02
1.04 GP corrected: Fill uses fill time 
added readSequenceFromCsv
"""

## ---------------------------------------------------------------------------
## known bugs/problems
_bugs_ = """
- if an elset appears twice in the sequence it seems to be deleted completely
- does not process ABAQUS parts/instances input files
  use command mdb.models[modelName].setValues(noPartsInputFile=ON)
"""

## ---------------------------------------------------------------------------
## import modules

from time import localtime, strftime
from sys import stdout
from bae.abq_model_02 import Model
from bae.misc_01 import groupIter, QuietLog
from bae.future_01 import *
import csv

############################################################################
############################################################################
############################################################################
############################################################################

def boolElsets(sourceFileName, setsequence, projectversion,
               renameSets = {'ALL': 'GRAV'}, verbose=1):
    r"""
    boolElsets

    usage:
    ======

    >>> from bae.makeHistory_01 import boolElsets, makeHistory
    >>> confirmedList, baseList, emptysetList = boolElsets(
    ...      sourceFileName, setsequence, projectversion,
    ...      renameSets, verbose)

    writes two input files as output: One ends with _IncElsets-NN.inp and
    contains the incremental elsets used in the sequence. The second ends with
    _BaseElsets-NN.inp and contains all elsets found in the model but not in
    setsequence (NN == projectversion)

    @param sourceFileName: is the name of the abaqus input file that contains
       all element set definitions mentioned in setsequence.
    @param setsequence: states the elset names in cronological order
    @param projectversion: is only added to the output file names, should be
       something like '01'
    @param renameSets: a dict stating elsets to be renamed in the output
    @param verbose: if set print diagnostic output (progress indicator) to
       stdout
    @returns: a tuple consisting of
       - confirmedList = chronological list of elset names used in the sequence
       - baseList = list of elset names found in the model but not in
         setsequence.
       - emptysetList = list of elsets that have no elements in them anymore
    @note: renameSets is applied to confirmedList and baseList, but not to
          emptysetList!
    """

    if verbose:
        logfile = stdout
    else:
        logfile = None

    # restore filename extension
    sourceFileName=sourceFileName.split(".inp")[0]+".inp"

    # see if files are readable, otherwise raise IOError
    sourceInputFile=open(sourceFileName)
    sourceInputFile.close()


    ## ------------------------------------------------------------------------
    ## read source input file

    m = Model(logfile=logfile)
    m.read(sourceFileName, ['ELSET'])

    ## ------------------------------------------------------------------------
    ## variables, check if all empty

    elsetList = list()
    for elset in setsequence:
        if elset in m.elset:
            elsetList.append(elset)
        elif verbose:
            print ('elset %s from the sequence not found in %s'
                   % (elset, sourceFileName))
    # shorter and faster but no diagnostic output:
    # elsetList = [elset for elset in setsequence if elset in m.elset]
    elsetN=len(elsetList)
    if elsetN==0:
        return list(), list(), list()

    ## ------------------------------------------------------------------------
    ## sort and split

    confirmedList=list()
    emptysetList=list()
    if len(elsetList[0])>0:
        confirmedList.append(elsetList[0])
    else:
        emptysetList.append(elsetList[0])

    if verbose:
        print ('( 0 / %d ) elset %s is the first and remains untouched.'
               % (elsetN, elsetList[0]))

    elsetCounter = 1
    for elsetName1 in elsetList[1:]:
        if verbose:
            print '(',elsetCounter,'/',elsetN-1,')',
        if elsetName1 in confirmedList:
            if verbose:
                print ("Warning: elset %s appears the second time in the"
                       " sequence. The second occurence is ignored."
                       % elsetName1)
            continue
            
        setN1 = m.elset[elsetName1]
        oldsetlen = len(setN1)
        for elsetName2 in elsetList[:elsetCounter]:
            setN1.difference_update(m.elset[elsetName2])
        if verbose:
            print ('elset %s shrinks from %d to %d elements.'
                   % (elsetName1, oldsetlen, len(setN1)))
        if len(setN1)>0:
            confirmedList.append(elsetName1)
        else:
            emptysetList.append(elsetName1)
        elsetCounter += 1

    baseList = list(set(m.elset).difference(setsequence))
    baseList.sort()

    ## ------------------------------------------------------------------------
    ## rename sets
    for oldName, newName in renameSets.iteritems():
        if oldName not in m.elset: continue
        if oldName.upper() == newName.upper(): continue

        if newName not in m.elset:
            m.elset[newName] = m.elset[oldName]
            try:
                ix = confirmedList.index(oldName)
                confirmedList[ix] = newName
            except ValueError:
                pass
            try:
                ix = baseList.index(oldName)
                baseList[ix] = newName
            except ValueError:
                pass
            if verbose:
                print ('elset %s is renamed to %s.'
                       % (oldName, newName))
        else:
            print "WARNING:"
            print ("When trying to rename elset %s to %s, the new elset has"
                   " been there already." % (oldName, newName))
            print ("Now adding the contents of elset %s to %s."
                   % (oldName, newName))
            m.elset[newName].update(m.elset[oldName])
            try:
                confirmedList.remove(oldName)
            except ValueError:
                pass
            try:
                baseList.remove(oldName)
            except ValueError:
                pass
        del m.elset[oldName]
    
    ## ------------------------------------------------------------------------
    ## write incremental elsets as input file

    outfile = open(sourceFileName.split(".inp")[0]
                   + "_IncElsets-%s.inp" % projectversion, 'w')
    for elsetName in confirmedList:
        m.elset.writeOne(outfile, elsetName)
    outfile.close()

    ## ---------------------------------------------------------------------------
    ## write base elsets as separate input file

    outfile = open(sourceFileName.split(".inp")[0]
                   + "_BaseElsets-%s.inp" % projectversion, 'w')
    for elsetName in baseList:
        m.elset.writeOne(outfile, elsetName)
    outfile.close()

    return confirmedList, baseList, emptysetList



#############################################################################
#############################################################################
#############################################################################
#############################################################################

def makeHistory(sequence, 
                projectname, projectversion, materialversion,
                elsetsInputFileName = None,
                elsetAllName='ALL', 

                gravTimeFactor=0.5,
                stepTimeEq = 10.,

                totalTimeAmplitudes = False,

                stringHeading=None,

# now some other strings you are more likely to leave untouched
stringStepDef="""\
**
** ----------------------------------------------------------------
**
** Step
**
*STEP
*DYNAMIC, EXPLICIT
,%(stepTimeEq)6.3f
**
*VARIABLE MASS SCALING, TYPE=BELOW MIN, DT=1.0e-2, NUMBER INTERVAL=10
**
*BULK VISCOSITY
0.1, 1.2
**
*BOUNDARY, TYPE=VEL
FIX,  1,3, 0.
** 
*DLOAD, OP=NEW, AMPLITUDE=ALL
ALL, GRAV, 9.81, 0., 0., -1.
** 
*OUTPUT, FIELD, NUMBER INTERVAL=1
*Element Output
S, LE, SDV, STATUS
*Node Output
U, V
** 
*End Step
** ----------------------------------------------------------------
**
** Step
**
*STEP
*DYNAMIC, EXPLICIT
,%(stepTime)6.3f
**
*VARIABLE MASS SCALING, TYPE=BELOW MIN, DT=1.0e-2, NUMBER INTERVAL=50
**
*BULK VISCOSITY
0.1, 1.2
** 
*BOUNDARY, OP=NEW
FIX,  1,3
** 
*DLOAD, OP=NEW, AMPLITUDE=GRAV
GRAV-01, GRAV, 9.81, 0., 0., -1.
** 
""",
#
stringStepEnd="""\
**
*OUTPUT, FIELD, TIME POINTS=TPOINTS
*Element Output
S, LE, SDV, STATUS
*Node Output
U, V
** 
*End Step
""",
#
stringTimePoints="""\
**
** Time Points
**
*TIME POINTS, NAME=TPOINTS
""",
#
stringIntroGravAmp=None,
stringGravAmp=None,
stringGravAmpCave=None,
#
stringIntroGravDload="""\
**
** Gravity Dload
""",
#
stringGravDload="""\
*DLOAD, OP=NEW, AMPLITUDE=%(ampName)s
 %(setName)s, GRAV, 9.81, 0., 0., -1.
""",
#
stringIntroSDVINI="""\
*INITIAL CONDITIONS,TYPE=SOLUTION
**
 ALL, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
 0.0, 1.0, 1.0, 100000.0, 100001., 100002.0, 100003.
""",
#
stringSDVINI="""\
 %(setName)s,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
 0.0, 1.0, 1.0, %(t1)6.3f, %(t2)6.3f,%(t3)6.3f,%(t4)6.3f
""",
#
stringCommentTemplate="""\
** 
** %(id)03d: %(comment)s
""",
#
stringStar="** \n",
stringSeparator="** "+"-"*75+"\n"
                ):
    """
    create the main input file for an abaqus run

    Timings: The values for the t1..t4 times are calculated from dt1..dt4 in
    the Times item of the sequence argument:

     - t1 excavation starts: startTime + dt1
     - t2 excavation ends: t1 + dt2
     - t3 refill: Pit, Fill: t2 + dt3; Cave: totalTime+1.0 = never
     - t4 element is deleted: Pit: t3 + dt4; Cave, Fill: totalTime+2.0 = never
    
    gravity amplitudes::

      CAVE -----t1                     (100%)
                   t2-------------t3t4 (80% / 99%, see stringGravAmpCave)
      
      PIT  -----t1                     (100%)
                 \                     
                  \                    
                   t2-------------t3t4 (0%)
      
      FILL -----t1             t4----- (100%)
                 \            /         
                  \          /           
                   t2------t3          (0%)

    @param sequence: all the info from the excel-sheet, it's a list of dict's
    Format:
    
    >>> sequence = [{
    >>>     'Name':'Equilibrium',             # comment line at begin of step
    >>>     'Time':['total,dt1,dt2,dt3,dt4'], # can overlap with next intervall
    >>>     'Fill':['Setname1','Setname2','...'], # Elset names
    >>>     'Pit':['Setname1','Setname2','...'],  # Elset names
    >>>     'Cave':['Setname1','Setname2','...'], # Elset names
    >>>     'Skip':'YES/NO'}]                 # optional, no ODB frame written.

    @param projectname: is what all the output file names start with

    @param projectversion: is only added to the output file names, should be
    something like '01'
    
    @param materialversion: is appended to the output file names, something like
    'LR01'

    @param elsetsInputFileName: is the name of the input file that contains all
    the elset definitions.
    Defaults to (projectname+'setsSequence'+projectversion+'.inp')

    @param elsetAllName: is the name of the elset containing all elements,
    defaults to 'ALL'

    @param totalTimeAmplitudes: If true write the amplitude commands with the
    TIME=TOTAL TIME option otherwise use the TIME=STEP TIME option. This flag
    also affects the default values for the arguments stringIntroGravAmp,
    stringGravAmp and stringGravAmpCave. Default: False.

    @param stringIntroGravAmp: Intro string for the gravity amplitudes. If
    totalTimeAmplitudes==False, defaults to::

        **
        ** Gravity Amplitudes
        **
        *AMPLITUDE, NAME=ALL, TIME=STEP TIME
        0.0, 0.0, %(stepTimeEq)6.3f, 1.0
        *AMPLITUDE, NAME=GRAV, TIME=STEP TIME
        0.0, 1.0, %(stepTime)6.3f, 1.0

    If totalTimeAmplitudes==True, defaults to::

        **
        ** Gravity Amplitudes
        **
        *AMPLITUDE, NAME=ALL, TIME=TOTAL TIME
        0.0, 0.0, %(stepTimeEq)6.3f, 1.0
        *AMPLITUDE, NAME=GRAV, TIME=TOTAL TIME
        0.0, 1.0, 100000.0, 1.0

    @param stringGravAmp: String for each of the gravity amplitudes. If
    totalTimeAmplitudes==False, defaults to::

        *AMPLITUDE, NAME=%(setName)s, TIME=STEP TIME
        %(t1)6.3f, 1., %(t2)6.3f, 0., %(t3)6.3f, 0., %(t4)6.3f, 1.

    If totalTimeAmplitudes==True, defaults to::

        *AMPLITUDE, NAME=%(setName)s, TIME=TOTAL TIME
        %(t1)6.3f, 1., %(t2)6.3f, 0., %(t3)6.3f, 0., %(t4)6.3f, 1.
    
    @param stringGravAmpCave: String for the gravity amplitudes for caves. It
    differs from the former that it only reduces the gravity during excavation
    by one percent. If totalTimeAmplitudes==False, defaults to::

        *AMPLITUDE, NAME=%(setName)s, TIME=STEP TIME
        %(t1)6.3f, 1., %(t2)6.3f, 0.99, %(t3)6.3f, 0.99, %(t4)6.3f, 1.

    If totalTimeAmplitudes==True, defaults to::

        *AMPLITUDE, NAME=%(setName)s, TIME=TOTAL TIME
        %(t1)6.3f, 1., %(t2)6.3f, 0.99, %(t3)6.3f, 0.99, %(t4)6.3f, 1.
    """

    # default arguments
    if elsetsInputFileName == None:
        elsetsInputFileName = ('%s_setsSequence-%s.inp'
                               % (projectname, projectversion))

    if stringHeading == None:
        raise ValueError, "stringHeading argument not supplied to makeHistory."

    if stringIntroGravAmp==None:
        if totalTimeAmplitudes:
            stringIntroGravAmp="""\
**
** Gravity Amplitudes
**
*AMPLITUDE, NAME=ALL, TIME=TOTAL TIME
0.0, 0.0, %(stepTimeEq)6.3f, 1.0
*AMPLITUDE, NAME=GRAV, TIME=TOTAL TIME
0.0, 1.0, 100000.0, 1.0
"""
        else:
            stringIntroGravAmp="""\
**
** Gravity Amplitudes
**
*AMPLITUDE, NAME=ALL, TIME=STEP TIME
0.0, 0.0, %(stepTimeEq)6.3f, 1.0
*AMPLITUDE, NAME=GRAV, TIME=STEP TIME
0.0, 1.0, %(stepTime)6.3f, 1.0
"""

    if stringGravAmp==None:
        if totalTimeAmplitudes:
            stringGravAmp="""\
*AMPLITUDE, NAME=%(setName)s, TIME=TOTAL TIME
%(t1)6.3f, 1., %(t2)6.3f, 0., %(t3)6.3f, 0., %(t4)6.3f, 1.
"""
        else:
            stringGravAmp="""\
*AMPLITUDE, NAME=%(setName)s, TIME=STEP TIME
%(t1)6.3f, 1., %(t2)6.3f, 0., %(t3)6.3f, 0., %(t4)6.3f, 1.
"""

    if stringGravAmpCave==None:
        if totalTimeAmplitudes:
            stringGravAmpCave="""\
*AMPLITUDE, NAME=%(setName)s, TIME=TOTAL TIME
%(t1)6.3f, 1., %(t2)6.3f, 0.99, %(t3)6.3f, 0.99, %(t4)6.3f, 1.
"""
        else:
            stringGravAmpCave="""\
*AMPLITUDE, NAME=%(setName)s, TIME=STEP TIME
%(t1)6.3f, 1., %(t2)6.3f, 0.99, %(t3)6.3f, 0.99, %(t4)6.3f, 1.
"""

    # constants
    setTypesList = ['Cave', 'Pit', 'Fill']

    #------------------------------------------------------------------

    # create elset sequence list
    set_sequence = list()
    for step in sequence:
        for setType in setTypesList:
            try:
                set_sequence.extend(step[setType])
            except KeyError:
                # if there is no setType - item in this step, just go on
                pass

    # add the elset containing all elements
    set_sequence.append(elsetAllName)

    # run boolElset to create incremental sets without intersection
    confirmedSets = boolElsets(
        elsetsInputFileName, set_sequence, projectversion,
        renameSets = {elsetAllName: 'GRAV-%s'%projectversion}, verbose=1)[0]

    ##-----------------------------------------------------------------------
    ## now here starts what the old makeHistory script does in its second run 
    ##-----------------------------------------------------------------------

    parameterdict = {
        'fileName': projectname,
        'matName' : materialversion,
        'version' : projectversion,
        }


    # -------------------------------------------------------------------------
    # eliminate sets not confirmed (i.e. empty sets after boolean operations)

    for step in sequence:
        for key in setTypesList:
            if key not in step:
                continue
            sequenceKeep=[]
            for setName in step[key]:
                if setName in confirmedSets:
                    sequenceKeep.append(setName)
                else:
                    print 'delete',setName
            step[key]=sequenceKeep

    # -------------------------------------------------------------------------
    # get step length adding up intervals and time points

    stepTime = 0
    for step in sequence:
        step['Start'] = stepTime
        stepTime += step['Time'][0]
    totalTime = stepTime+stepTimeEq



    # -------------------------------------------------------------------------

    outfile='%(fileName)s_%(matName)s-%(version)s.inp' % parameterdict
    output=open(outfile,'w')

    output.write(stringSeparator)
    output.write("** %s - makeHistory_01.makeHistory, version %s\n"
                 % (strftime("%a, %d %b %Y %H:%M:%S", localtime()),_version_))
    output.write(stringHeading % parameterdict)

    # -------------------------------------------------------------------------
    # initial conditions (TOTAL TIME)

    output.write(stringIntroSDVINI)

    n=0
    for step in sequence:
        n=n+1
        time = step['Start']+stepTimeEq
        output.write(stringCommentTemplate % {'id':n,'comment':step['Name']})

        # CAVE
        t1 = time+step['Time'][1]
        t2 = time+step['Time'][1]+step['Time'][2]
        t3 = totalTime+1.0
        t4 = totalTime+2.0
        if step.has_key('Cave'):
            for setName in step['Cave']:
                output.write(stringSDVINI % {'setName':setName,
                                             't1':t1,'t2':t2,'t3':t3,'t4':t4})

        # PIT
        t1 = time+step['Time'][1]
        t2 = time+step['Time'][1]+step['Time'][2]
        t3 = time+step['Time'][1]+step['Time'][2]+step['Time'][3]
        t4 = time+step['Time'][1]+step['Time'][2]+step['Time'][3]+step['Time'][4]
        if step.has_key('Pit'):
            for setName in step['Pit']:
                output.write(stringSDVINI % {'setName':setName,
                                             't1':t1,'t2':t2,'t3':t3,'t4':t4})

        # FILL
        t1 = time+step['Time'][1]
        t2 = time+step['Time'][1]+step['Time'][2]
        t3 = time+step['Time'][1]+step['Time'][2]+step['Time'][3]
        ## t4 = time+step['Time'][1]+step['Time'][2]+step['Time'][3]+step['Time'][4]
        t4 = totalTime+2.0
        if step.has_key('Fill'):
            for setName in step['Fill']:
                output.write(stringSDVINI % {'setName':setName,
                                             't1':t1,'t2':t2,'t3':t3,'t4':t4})

    output.write(stringStar)

    # -------------------------------------------------------------------------
    # gravity amplitudes (STEP TIME)

    output.write(stringSeparator)
    output.write(stringIntroGravAmp %
                 {'stepTimeEq':stepTimeEq, 'stepTime':stepTime})

    n=0
    for step in sequence:
        n=n+1
        time = step['Start']
        output.write(stringCommentTemplate % {'id':n,'comment':step['Name']})

        # CAVE -----t1
        #              t2-------------t3t4 (80%)
        t1 = time+step['Time'][1]
        t2 = time+step['Time'][1]+step['Time'][2]*gravTimeFactor
        (t3,t4) = (1000003,1000004)
        if totalTimeAmplitudes:
            (t1,t2,t3,t4) = [tt+stepTimeEq for tt in (t1,t2,t3,t4)]
        if step.has_key('Cave'):
            for setName in step['Cave']:
                output.write(stringGravAmpCave % {'setName':setName,
                                                  't1':t1,'t2':t2,'t3':t3,'t4':t4})

        # PIT  -----t1
        #            \
        #             \
        #              t2-------------t3t4
        t1 = time+step['Time'][1]
        t2 = time+step['Time'][1]+step['Time'][2]*gravTimeFactor
        (t3,t4) = (1000003,1000004)
        if totalTimeAmplitudes:
            (t1,t2,t3,t4) = [tt+stepTimeEq for tt in (t1,t2,t3,t4)]
        if step.has_key('Pit'):
            for setName in step['Pit']:
                output.write(stringGravAmp % {'setName':setName,'t1':t1,'t2':t2,'t3':t3,'t4':t4})

        # FILL -----t1              t4-----
        #            \            /
        #             \          /
        #              t2------t3
        t1 = time+step['Time'][1]
        t2 = time+step['Time'][1]+step['Time'][2]*gravTimeFactor
        t3 = time+step['Time'][1]+step['Time'][2]+step['Time'][3]+step['Time'][4]*(1.0-gravTimeFactor)
        t4 = time+step['Time'][1]+step['Time'][2]+step['Time'][3]+step['Time'][4]
        if totalTimeAmplitudes:
            (t1,t2,t3,t4) = [tt+stepTimeEq for tt in (t1,t2,t3,t4)]
        if step.has_key('Fill'):
            for setName in step['Fill']:
                output.write(stringGravAmp % {'setName':setName,'t1':t1,'t2':t2,'t3':t3,'t4':t4})

    output.write(stringStar)

    # -------------------------------------------------------------------------
    # time points 

    output.write(stringSeparator)
    output.write(stringTimePoints)
    timesList = [step['Start']+step['Time'][0]
                 for step in sequence
                 if (('Skip' not in step)
                     or (step['Skip'] in ['No','NO','Nada','Nope']))]

    for line in groupIter(timesList, maxcounts=8, format='%6.3f'):
        output.write(", ".join(line)+"\n")
    output.write(stringStar)

    # -------------------------------------------------------------------------
    # step 

    parameterdict.update({'stepTimeEq':stepTimeEq, 'stepTime':stepTime})
    output.write(stringStepDef % parameterdict)

    # -------------------------------------------------------------------------
    # dload 

    output.write(stringSeparator)
    output.write(stringIntroGravDload)

    n=0
    for step in sequence:
        n=n+1

        output.write(stringCommentTemplate % {'id':n,'comment':step['Name']})

        step['Sets']=[]
        if step.has_key('Cave'):
            for setName in step['Cave']:
                output.write(stringGravDload % {'ampName':setName,'setName':setName})
                step['Sets'].append(setName)
        if step.has_key('Pit'):
            for setName in step['Pit']:
                output.write(stringGravDload % {'ampName':setName,'setName':setName})
                step['Sets'].append(setName)
        if step.has_key('Fill'):
            for setName in step['Fill']:
                output.write(stringGravDload % {'ampName':setName,'setName':setName})
                step['Sets'].append(setName)

    output.write(stringStar)

    # -------------------------------------------------------------------------
    # step end

    output.write(stringSeparator)
    output.write(stringStepEnd)

    output.close()


    # -------------------------------------------------------------------------
    # post sets / adds sets from skipped ODB output in next line

    stringSeqIn="""
# -------------------------------------------------------------------------
# Sequence %(setType)s

modelHistory%(setType)s = [\\
"""
    stringSeqOut="  ]\n"

    setTypeString={'Sets':'','Fill':'Support','Cave':'Cave','Pit':'Pit'}

    outfile='%(fileName)s_%(version)s_postSets.inp' % parameterdict
    output=open(outfile,'w')

    output.write('# POST SETS %(version)s\n' % parameterdict)


    for setType in ['Sets','Fill','Cave','Pit']:

        output.write(stringSeqIn % {'setType':setTypeString[setType]})

        n=0
        skipSets=[]
        for step in sequence:
            if not step.has_key('Skip') or (step['Skip'] in ['No','NO','Nada','Nope']):
                output.write('    (')
                if skipSets:
                    for setName in skipSets:
                        output.write("'PART-1-1."+setName+"', ")
                    skipSets=[]
                if step.has_key(setType):
                    for setName in step[setType]:
                        output.write("'PART-1-1."+setName+"', ")
                output.write('),\\\n' )
            else:
                if step.has_key(setType):
                    for setName in step[setType]:
                        skipSets.append(setName)

        output.write(stringSeqOut)

    for setType in ['Name']:

        output.write('\nsequenceFrameNames=[')
        for step in sequence:
            if not step.has_key('Skip') or (step['Skip'] in ['No','NO','Nada','Nope']):
                if step.has_key(setType):
                    output.write("'"+step[setType]+"',")
        output.write(']\n')

    output.close()


#############################################################################
#############################################################################
#############################################################################
#############################################################################

def readSequenceFromCsv(sequenceCsvFile,
                        defaultTimeList = [1.0,0.0,0.0,0.0,0.0],
                        logfile=stdout):
    r"""

    @param sequenceCsvFile: either a filename or an already open file that
      contains the sequence data

    @param defaultTimeList: a list of timing values (floats) containing the
      default values if the corresponding column as a whole or the value in
      this particular row does not exist in the sequence csv file. The list
      contains the five values for

       1. sequence step total length or frame length,
       2. dt1: time period from the beginning of the step to the start of the
          excavation
       3. dt2: duration of the excavation
       4. dt3: time period between finished excavation and refill.
          Refill is the instantanous change to linear elastic usually softer
          material properties.
       5. dt4: time period between refill and deletion of the elements.

      Those values are replaced (if applicable) by the values found in columns
      labeled "Time total", "Time dt1", "Time dt2", "Time dt3" and "Time dt4".

    @param logfile: if provided it is an already opened log file
      for diagnostic output, use None to suppress diagnostic output,
      if not provided defaults to stdout 

    @returns: a sequence list suitable as first argument to the makeHistory
      function

    @note: the values for dt3 and dt4 are ignored for Cave and dt4 is ignored
      for Fill
    """

    if logfile == None:
        logfile = QuietLog()

    if isinstance(sequenceCsvFile, basestring):
        sequenceCsvFile = open(sequenceCsvFile, "rb")

    sequenceCsv = csv.reader(sequenceCsvFile)

    # find the first line that is not empty
    for headline in sequenceCsv:
        if len(headline) and any(map(len, headline)):
            break

    # valid fields
    validCategories = ['Cave','Pit','Fill']
    validFields = [('Skip', str.upper), ('Time', eval), ('Name', str.upper)]
    validTimeFields =["Time "+item for item in "total","dt1","dt2","dt3","dt4"]
    validFields.extend([(item, float) for item in validTimeFields])

    # evaluate headline
    colNbToKeyConverter = dict()
    class FieldFound(Exception): pass
    for cnt, item in enumerate(headline):
        try:
            for valid in validCategories:
                if item.upper().startswith(valid.upper()):
                    converter = None
                    raise FieldFound()
            for valid, converter in validFields:
                if item.upper().startswith(valid.upper()):
                    raise FieldFound()
        except FieldFound:
            colNbToKeyConverter[cnt] = (valid, converter)
        else:
            if len(item):
                colContent = " (heading '%s')"%item
            else:
                colContent = " (empty heading)"
            logfile.write("Column number %d%s does not match any of the valid"
                          " categories or fields, it will be ignored.\n"
                          % (cnt+1, colContent))

    # loop over each row / line
    sequence = list()
    for row in sequenceCsv:

        # ignore empty lines
        if len(row)==0 or not(any(map(len, row))):
            continue

        # add elsets and other fields
        sequenceItem = dict()
        for category in validCategories:
            sequenceItem[category] = list()

        for col, item in enumerate(row):

            # ignore empty fields
            if item=="":
                continue

            # look up to which category this column belongs
            try:
                (category, converter) = colNbToKeyConverter[col]
            except KeyError:
                # ignore columns to be ignored
                continue

            # process this item / column
            if converter==None:
                sequenceItem[category].append(item)
            else:
                sequenceItem[category] = converter(item)

        # special treatment of Time items
        timeList = list(defaultTimeList)
        for cnt, category in enumerate(validTimeFields):
            try:
                timeList[cnt] = sequenceItem[category]
                del sequenceItem[category]
            except KeyError:
                pass
        sequenceItem['Time'] = timeList

        # store sequenceItem
        sequence.append(sequenceItem)

    return sequence
