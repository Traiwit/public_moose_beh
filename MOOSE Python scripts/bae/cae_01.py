"""Some tools to use from within cae. They are not supposed to be run as
ordinary python scripts but rely on the cae environment, i.e. they have to be
launched from an interactiv cae session or by abaqus cae noGUI=myscript.py
"""

__version__ = "1.02"

_version_history_ = r"""
1.00 GP: no versioning so far
1.01 GP: added getFieldValueAtPath(), getScriptName()
        incompatible change: removed odbVarLabels["PST"] because we always use
        logP and sometimes call it PST.
1.02 GP added exportElsets()
"""

from bae.abq_model_02 import Model
from bae.misc_01 import quietlog, LogFile
from bae.colormap_01 import colormap_rainbow
from abaqus import session, mdb
from abaqusConstants import *
import sys, csv, itertools, os

_debug = False

def write_nset(setname):
    """
    - create an nset named 'myset' (or whatever)
    - write: 'write_nset("myset")' on the cae prompt.
    """
    m=Model()
    vpn = session.currentViewportName
    obj = session.viewports[vpn].displayedObject
    nodes = obj.sets[setname].nodes
    m.nset[setname] = set([n.label for n in nodes])
    m.write('nset_'+setname+'.inp', ['NSET'])

def write_elset(setname):
    """
    - create an elset named 'myset' (or whatever)
    - write: 'write_elset("myset")' on the cae prompt.
    """
    m=Model()
    vpn = session.currentViewportName
    obj = session.viewports[vpn].displayedObject
    if type(obj).__name__=='Odb':

        # try to find it in the rootassembly
        # in the rootassembly:
        # rootAssembly.elementSets[elsetName].elements is a tuple of
        # OdbMeshElementArray-objects (i.e. elements)
        elsetContainer = obj.rootAssembly.elementSets
        try:
            els1 = elsetContainer[setname]
        except KeyError:
            elsetContainer = None

        # try to find it in a part instance
        # in part instances:
        # instance.elementSets[elsetName].elements is an OdbMeshElementArray
        if elsetContainer == None:
            inst = None
            if "." in setname:
                partName, partsetname = setname.split(".")
                try:
                    inst = obj.rootAssembly.instances[partName]
                except KeyError:
                    pass
                else:
                    elsetContainer = inst.elementSets
                    try:
                        els1 = elsetContainer[partsetname]
                    except KeyError:
                        els1 = None     
            if inst == None:
                # try to find in any part
                for inst in obj.rootAssembly.instances.values():
                    elsetContainer = inst.elementSets
                    try:
                        els1 = elsetContainer[setname]
                    except KeyError:
                        els1 = None
            if els1==None:
                elsetContainer = None

        if elsetContainer == None:
            raise KeyError("Elset %s not found in the odb." % setname)

        # gather all elements
        elset = set()
        for els2 in els1.elements:
            if type(els2).__name__ == 'OdbMeshElementArray':
                elset.update([e.label for e in els2])
            else:
                raise Exception("Element set of type %s not implemented so far"
                                % type(els2).__name__)
    else:
        # in mdb:
        # obj.sets[elsetName].elements is an OdbMeshElementArray
        els = obj.sets[setname].elements
        elset = set([e.label for e in els])
    m.elset[setname] = elset
    m.write('elset_'+setname+'.inp', ['ELSET'])


def exportElsets(filename="elsets.inp"):
    """export all elsets to an input file

    Note: There is an odb_02.OdbReader.getAbqModel() function with an optional
    recognizedBlocks argument that does the same if the current displayed
    object is an odb:
     >>> from bae.cae_01 import initOdb
     >>> from bae.odb_02 import OdbReader
     >>> svp, odb = initOdb()
     >>> odb2 = OdbReader(odb)
     >>> model = odb2.getAbqModel(recognizedBlocks = "ELSET")
     >>> model.write("elsets.inp")
    """
    m=Model()
    vpn = session.currentViewportName
    obj = session.viewports[vpn].displayedObject

    def getElset(elementArray):
        if type(elementArray).__name__ in ('OdbMeshElementArray', 'MeshElementArray'):
            return set(e.label for e in elementArray)
        else:
            raise Exception("Element set of type %s not implemented so far"
                            % type(elementArray).__name__)
        
    if type(obj).__name__=='Odb':
        # in the rootassembly:
        # rootAssembly.elementSets[elsetName].elements is a tuple of
        # OdbMeshElementArray-objects (i.e. elements)
        elsetContainer = obj.rootAssembly.elementSets
        for elsetName in elsetContainer.keys():
            elset = set()
            for elementArray in obj.rootAssembly.elementSets[elsetName].elements:
                elset.update(getElset(elementArray))
            m.forceElset(elsetName).update(elset)
        # in part instances:
        # instance.elementSets[elsetName].elements is an OdbMeshElementArray
        for inst in obj.rootAssembly.instances.values():
            elsetContainer = inst.elementSets
            for elsetName in elsetContainer.keys():
                elset = getElset(elsetContainer[elsetName].elements)
                m.forceElset(elsetName).update(elset)
    else:
        # in mdb:
        # obj.sets[elsetName].elements is an OdbMeshElementArray
        elsetContainer = obj.sets
        for elsetName in elsetContainer.keys():
            elset = getElset(elsetContainer[elsetName].elements)
            m.forceElset(elsetName).update(elset)

    m.write(filename, ['ELSET'])


def write_surface(name):
    """
    - create an surface named 'mysurf' (or whatever)
    - write: 'write_surface("mysurf")' on the cae prompt.
    - only works for element based surfaces on tets
    """
    #
    faceSym2Num = {FACE1: 1, FACE2: 2, FACE3: 3, FACE4: 4, FACE5: 5, FACE6: 6}
    def setname(name, facenum):
        return "%s_S%1d" % (name.upper(),facenum)
    #
    m=Model()
    vpn = session.currentViewportName
    obj = session.viewports[vpn].displayedObject
    surf = obj.surfaces[name]
    els = surf.elements
    sides = surf.sides
    allfaces = set()
    for i in range(len(els)):
        try:
            facenum = faceSym2Num[sides[i]]
        except KeyError:
            raise Exception(
                "value for session.viewports['%s'].displayedObject.surfaces"
                "['%s'].sides[%d] unknown to write_surface: %s"
                % (vpn, name, i, str(sides[i])))
        allfaces.add(facenum)
        m.forceElset(setname(name,facenum)).add(els[i].label)

    newsurf = dict()
    m.surface[name] = newsurf
    for facenum in allfaces:
        eln = setname(name,facenum)
        if eln in m.elset:
            newsurf["S%1d"%facenum] = eln

    m.write('surface_'+name+'.inp', ['ELSET', 'SURFACE'])


##--------------------------------------------------------------------

def noGUIStart(odbName):
    """
    This function is obsolete. Use initOdb instead!

    if your script was started with the -noGUI option on the command line
    create a viewport, open the given odb and display this odb in the viewport.
    After that your script should work with the noGUI option as it did from the
    File -> Run Script menu option of the interactive cae session.
    """
    if '-noGUI' in sys.argv:
        session.Viewport(
            name='Viewport: 1', origin=(0.0, 0.0), width=320.,height=240.)
        session.viewports['Viewport: 1'].makeCurrent()
        session.viewports['Viewport: 1'].maximize()
    
        o1 = session.openOdb(name=odbName)
        session.viewports['Viewport: 1'].setValues(displayedObject=o1)

        logFile = open(sys.argv[3].replace(".py", ".log"), "w")
    else:
        logFile = sys.stdout

    return logFile

def initOdb(odbName = None):
    """
    Return the current or the specified odb and the current viewport.

    If the "displayed" object in the current viewport is not an open odb then
    this function opens the odb with the specified odbName and make it the
    displayed object in the current viewport. Otherwise the argument odbName is
    ignored and the currently displayed odb is used.

    This function should make sure that your script works with the noGUI option
    as it did from the File -> Run Script menu option of the interactive cae
    session.

    @returns: A tuple (svp, odb) with svp being the current viewport and odb
    the current odb object.
    """
    svp = session.viewports[session.currentViewportName]
    odb = svp.displayedObject
    if (type(svp.displayedObject).__name__ != "Odb"
        or svp.displayedObject.closed):
        import viewerModules
        # ... and there must be an odbName argument
        odb = session.openOdb(name=odbName)
        svp.setValues(displayedObject=odb)

    return (svp, odb)

def getScriptName():
    r"""extracts the name of the running script from sys.argv

    @returns: The script name or "<unknown>" if the detection failed.
    """
    scriptName = "<unknown>"
    
    try:
        scriptNamePos = sys.argv.index('-noGUI') + 1
    except ValueError:
        pass
    else:
        if len(sys.argv)>scriptNamePos and sys.argv[scriptNamePos][0]!="-":
            # next argument is the script name
            scriptName = sys.argv[scriptNamePos]

    return scriptName

def openDefaultLogFile():
    r"""
    Open a logfile for diagnostic output.

    The name of the logfile with a .log extension is derived from the name of
    the script if it can be found in sys.argv, otherwise from the open odb.

    @returns: the open logfile
    """
    scriptName = getScriptName()

    if scriptName!="<unknown>":
        logfileName = scriptName.rsplit(".py",1)[0] + ".log"
    else: # if scriptName=="<unknown>"
        # no script name found, so take the odb name as base
        svp = session.viewports[session.currentViewportName]
        odbName = svp.odbDisplay.name
        logfileName = odbName.rsplit(".odb",1)[0] + ".log"

    # open logfile
    logfile = LogFile(logfileName, scriptName=scriptName)

    return logfile


##--------------------------------------------------------------------

def annotatePoints(labels, coords, colors=None,
                   xOffset=0, yOffset=10,
                   font="-*-helvetica-medium-r-normal--12-*",
                   colormap=colormap_rainbow
                   ):
    """
    create annotations with a label, pointing to a point all from lists.

    >>> from bae.misc_01 import DictTableFromCsvFile
    >>> from bae.cae_01 import annotatePoints
    >>> tab = DictTableFromCsvFile("myannotions.csv")
    >>> labels = tab["labels"]
    >>> coords = zip(tab["X"], tab["Y"], tab["Z"])
    >>> colors = tab["colors"]
    >>> annotatePoints(labels, coords, colors)

    The labels, coords and colors lists must have the same length!

    @param labels: list of annotation labels
    @param coords: list of coordinate tuples
    @param colors: list of values to be mapped to a colormap (floats or ints)

    @param xOffset: x-offset of the text to the point
    @param yOffset: y-offset of the text to the point
    @param colormap: a Colormap object from module bae.colormap_01
    """

    svp = session.viewports[session.currentViewportName]

    if colors is None:
        for ID, coord in itertools.izip(labels, coords):
            coord = tuple(coord)
            arrow = mdb.Arrow(
                name='Arrow-'+ID, startPoint=(xOffset, yOffset),
                startAnchor=coord, endAnchor=coord)
            svp.plotAnnotation(annotation=arrow)
            text = mdb.Text(
                name=ID, text=ID, offset=(xOffset, yOffset),
                anchor=coord, font=font)
            svp.plotAnnotation(annotation=text)
    else:
        colormap.set_xrange(min(colors), max(colors))
        for ID, coord, color in itertools.izip(labels, coords, colors):
            coord = tuple(coord)
            colorText = "#%02X%02X%02X" % tuple(colormap.get_rgb_color(color))
            arrow = mdb.Arrow(
                name='Arrow-'+ID, startPoint=(xOffset, yOffset),
                startAnchor=coord, endAnchor=coord, color=colorText)
            svp.plotAnnotation(annotation=arrow)
            text = mdb.Text(
                name=ID, text=ID, offset=(xOffset, yOffset),
                anchor=coord, font=font, color=colorText)
            svp.plotAnnotation(annotation=text)

##--------------------------------------------------------------------

class OdbVarLabel(object):
    """
    object to store the data you need to access a certain variable in the odb
    """
    __slots__ = ("label", "position", "refinement")
    def __init__(self, label, position=None, refinement=()):
        """
        Constructor

        If you only specify the label, tries to get the data from the
        odbVarLabels dictionary.
        """
        try:
            if position!=None:
                raise KeyError("Not a real error, jump to the except clause.")
            self.label, self.position, self.refinement = odbVarLabels[label]
        except KeyError:
            self.label = label
            self.position = position
            self.refinement = refinement

# Predefined variables for the varLabels argument to getFieldValueAtPoints()
# At the moment contains:
# PST (=SDV3), logP (=SDV4), U_1 ..., UMAG, V_1 ..., VMAG, SMAJ, SMIN,
# SDV1 ... SDV15
odbVarLabels = {
    # 'PST': OdbVarLabel(label='SDV3', position=INTEGRATION_POINT, refinement=()),
    'logP':OdbVarLabel(label='SDV4', position=INTEGRATION_POINT, refinement=()),

    'U_1': OdbVarLabel(label='U', position=NODAL, refinement=(COMPONENT,'U1')),
    'U_2': OdbVarLabel(label='U', position=NODAL, refinement=(COMPONENT,'U2')),
    'U_3': OdbVarLabel(label='U', position=NODAL, refinement=(COMPONENT,'U3')),
    'UMAG':OdbVarLabel(label='U', position=NODAL,
                       refinement=(INVARIANT, 'Magnitude')),
    'V_1': OdbVarLabel(label='V', position=NODAL, refinement=(COMPONENT,'V1')),
    'V_2': OdbVarLabel(label='V', position=NODAL, refinement=(COMPONENT,'V2')),
    'V_3': OdbVarLabel(label='V', position=NODAL, refinement=(COMPONENT,'V3')),
    'VMAG':OdbVarLabel(label='V', position=NODAL,
                       refinement=(INVARIANT, 'Magnitude')),

    'S_11': OdbVarLabel('S', INTEGRATION_POINT, (COMPONENT,'S11')),
    'S_22': OdbVarLabel('S', INTEGRATION_POINT, (COMPONENT,'S22')),
    'S_33': OdbVarLabel('S', INTEGRATION_POINT, (COMPONENT,'S33')),
    'S_12': OdbVarLabel('S', INTEGRATION_POINT, (COMPONENT,'S12')),
    'S_23': OdbVarLabel('S', INTEGRATION_POINT, (COMPONENT,'S23')),
    'S_13': OdbVarLabel('S', INTEGRATION_POINT, (COMPONENT,'S13')),
    'S_MIN_PRINCIPAL': OdbVarLabel(label='S', position=INTEGRATION_POINT,
                        refinement=(INVARIANT,'Min. Principal')),
    'S_MAX_PRINCIPAL': OdbVarLabel(label='S', position=INTEGRATION_POINT,
                        refinement=(INVARIANT,'Max. Principal')),
    'S_MID_PRINCIPAL': OdbVarLabel(label='S', position=INTEGRATION_POINT,
                        refinement=(INVARIANT,'Mid. Principal')),
    # don't use the following: just for compatibility reasons
    'SMAJ': OdbVarLabel(label='S', position=INTEGRATION_POINT,
                        refinement=(INVARIANT,'Min. Principal')),
    'SMIN': OdbVarLabel(label='S', position=INTEGRATION_POINT,
                        refinement=(INVARIANT,'Max. Principal')),
    'SMID': OdbVarLabel(label='S', position=INTEGRATION_POINT,
                        refinement=(INVARIANT,'Mid. Principal')),
    }
for i in range(15):
    varName = "SDV%d" % i
    odbVarLabels[varName] = OdbVarLabel(
        label=varName, position=INTEGRATION_POINT, refinement=())


def getFieldValueAtPoints(
    varLabels, pointCoords, shape='UNDEFORMED',stepIndex=-1, frameList=range(10000),
    odbName=None, tempFileName="tempReport.rpt", logfile="CreateDefaultLog", averageElementOutput=True):
    """
    Get field values at specified points out of an odb using the cae path
    feature.

    @param varLabels: Something like this:

    >>> from abaqusConstants import *
    >>> varLabels = {
    ...     'PST':{'label':'SDV3', 'position':INTEGRATION_POINT,
    ...     'refinement':()},
    ...     'SMAJ':{'label':'S', 'position':INTEGRATION_POINT,
    ...             'refinement':(INVARIANT,'Min. Principal')},
    ...     'SMIN':{'label':'S', 'position':INTEGRATION_POINT,
    ...             'refinement':(INVARIANT,'Max. Principal')},
    ...     'Uxx':{'label':'U', 'position':NODAL,
    ...            'refinement':(COMPONENT, 'U1')},
    ...     'Uyy':{'label':'U', 'position':NODAL,
    ...            'refinement':(COMPONENT, 'U2')},
    ...     'Uzz':{'label':'U', 'position':NODAL,
    ...            'refinement':(COMPONENT, 'U3')},
    ...     'UHor':{'label':'UHor', 'position':NODAL,
    ...            'refinement':()},
    ...     }

    Or like this:

    >>> varLabels = {'UVERT' : 'U_3',
    ...              'UHor': {'label':'UHor', 'position':NODAL,
    ...                       'refinement':()},
    ...              'UNorth': OdbVarLabel(label='U', position=NODAL,
    ...                        refinement=(COMPONENT,'U2')),
    ...             }
    
    Or simply like this:

    >>> varLabels = ["PST", "SMAJ"]

    Any variable specified only by label (U_3, PST, SMAJ in the examples above)
    must be contained in the dictionary odbVarLabels.

    Ordering of the variables in the resulting csv file: If varLabels is a list
    the order in that list is retained. If varLabels is a dictionary its keys
    are ordered alphabetically determining the order in the result.

    @param pointCoords:

    @param shape: 'UNDEFORMED' or 'DEFORMED': Supposedly the path moves with the material 
      if you specify 'UNDEFORMED'. For 'UNDEFORMED', the value is extracted for a material
      point, that was sitting at the given position in the UNDEFORMED state.
      If you specify 'DEFORMED' the values are extracted for a particle sitting at the 
      given position in the 'UNDEFORMED' state. Needs to be checked. Whoever is sure, that 
      the above is true, should remove the questionmarks here...

    @param stepIndex: zero based step index, -1 is the last, -2 the second
      last ...
    @param frameList: list of zero based frame index numbers, frames that are
      not in the odb are silently ignored so frameList=range(10000) means all
      existing frames.
    @param odbName: If specified and if there is not already an odb open in the
      current viewport then open this odb and make it the displayed object in
      the current viewport.
    @param tempFileName:
    @param logfile: Where to write diagnostic output?
     - If not specified create a new logfile according to the initOdb()
       function.
     - If None: no output.
     - If it's a string write to a file with that name.
     - Anything else is treated as the logfile to write to.

    @returns: A tuple (varList, dict).
    
    varList: the ordered keys of the varLabels argument and       
    
    dict: is a dictionary {frame number: data list}
    where data list is a list of list-per-point lists, one per point in the
    pointCoords list. Each of those list-per-point lists finally contains
    the values in the order of varList.

    So result looks something like that:
    {frame number: [ pt1-[var1, var2, ...], pt2-[var1, var2, ...], ...]}
    """
    global odbVarLabels

    (svp, odb) = initOdb(odbName)
    if logfile == "CreateDefaultLog":
        logfile = openDefaultLogFile()
    elif logfile == None:
        logfile = quietlog
    elif isinstance(logfile, basestring):
        logfile = open(logfile, "w")

    # path from points
    pth = session.Path(name="pointlistpath", type=POINT_LIST,
                       expression=pointCoords)
    if shape.upper()=='UNDEFORMED':
        shape = UNDEFORMED
    elif shape.upper()=='DEFORMED':
        shape = DEFORMED
    else:
        raise ValueError("shape can only be 'DEFORMED' or 'UNDEFORMED'.")

    # correct the step index
    nbSteps = len(odb.steps)
    while stepIndex<0:
        stepIndex += nbSteps
    # filter frameList: only valid frames
    validFrames = set(range(len(odb.steps.values()[stepIndex].frames)))
    for idx, frame in enumerate(frameList):
        while frame<0:
            frame += len(validFrames)
            frameList[idx] = frame
    validFrames.intersection_update(frameList)
    frameList = list(validFrames)
    frameList.sort()

    nbFrames = len(frameList)
    logfile.write("will examine %d frames of step %s\n"
                  % (nbFrames, odb.steps.keys()[stepIndex]))

    if isinstance(varLabels, list):
        varList = varLabels
        varLabels = dict()
        for var in varList:
            try:
                varLabels[var] = odbVarLabels[var]
            except KeyError:
                raise ValueError(
                    "Odb variable %s not recognized (in odbVarLabels dict)"
                    % (var))
    elif isinstance(varLabels, dict):
        varList = varLabels.keys()
        varList.sort()
        for varName in varList:
            var = varLabels[varName]
            if isinstance(var, str):
                try:
                    varLabels[varName] = odbVarLabels[var]
                except KeyError:
                    raise ValueError(
                        "Odb variable %s not recognized (in odbVarLabels dict)"
                        % (var))
            elif isinstance(var, dict):
                try:
                    varLabels[varName] = OdbVarLabel(**var)
                except TypeError:
                    raise ValueError(
                        "Wrong items in varLabels dictionary argument, odb"
                        " variable description for variable %s."
                        % varName)
    else:
        raise ValueError("varLabels argument must be a dict or a list.")

    result = dict()

    for frameIdx in frameList:
        svp.odbDisplay.setFrame(step=stepIndex, frame=frameIdx)
        svp.odbDisplay.basicOptions.setValues(averageElementOutput=averageElementOutput)
        logfile.write ("processing frame number %d, extracting data to report"
                       " file %s\n"% (frameIdx, tempFileName))

        xyDataList = list()
        for var in varList:
            odbVar = varLabels[var]
            svp.odbDisplay.setPrimaryVariable(
                variableLabel=odbVar.label,
                outputPosition=odbVar.position,
                refinement=odbVar.refinement)

            # xyName = "%s_%03d" % (var, svp.odbDisplay.fieldFrame[1])
            xyName = var
            xyData = session.XYDataFromPath(
                name=xyName, path=pth, includeIntersections=FALSE,
                shape=shape, labelType=TRUE_DISTANCE)
            xyDataList.append(xyData)

        xyDataList = tuple(xyDataList)
        session.writeXYReport(
            fileName=tempFileName,
            appendMode=OFF, xyData=xyDataList)

        logfile.write("processing frame number %d, moving data from the"
                      " temporary report file %s to the output data"
                      " structures.\n"
                      % (frameIdx, tempFileName))

        tempFile = open(tempFileName)
        # skip first three lines
        tempFile.next(); tempFile.next(); tempFile.next()
        thisU = list() # result list for this one frame
        lastX = None
        for line in tempFile:
            vals = line.strip().split()
            if len(vals)==0 or lastX == vals[0]:
                continue
            lastX = vals[0]
            try:
                u = map(float, vals[1:])
            except ValueError:
                raise Exception(
                    "ERROR: Could not convert all values (%s) to floats.\n"
                    "The code could possibly be made more flexible, but at the"
                    " the moments it's not."
                    % vals[1:])
            thisU.append(u)
        tempFile.close()
        if _debug:
            logfile.write("frame %d, u: %s\n" % (frameIdx, thisU))

        if len(thisU)!=len(pointCoords):
            raise Exception(
                "ERROR: Found %d values in the report file %s for %d specified"
                " points." % (len(thisU), tempFileName, len(pointCoords)))

        result[frameIdx] = thisU
        logfile.flush()

    #-- fini
    try:
        os.remove(tempFileName)
    except OSError:
        pass
    return varList, result


def getFieldValueAtPath(
    varLabels, pointCoords, shape='UNDEFORMED',stepIndex=-1, frameList=range(1001),
    odbName=None, tempFileName="tempReport.rpt", logfile="CreateDefaultLog",
    averageElementOutput=True):
    """
    Get field values along a specified path out of an odb using the cae path
    feature. There are values at each of the given points and all element
    boundary intersections in between.

    @param varLabels:
      Something like this:
        >>> from abaqusConstants import *
        >>> varLabels = [
        ...     # [name, [label, position, refinement]]
        ...     ['PST', ['SDV3', INTEGRATION_POINT, ()]],
        ...     ['S_MIN_PRINCIPAL', ['S', INTEGRATION_POINT, (INVARIANT,'Min. Principal')]],
        ...     ['S_MAX_PRINCIPAL', ['S', INTEGRATION_POINT, (INVARIANT,'Max. Principal')]],
        ...     ['Uxx', ['U', NODAL, (COMPONENT, 'U1')]],
        ...     ['Uyy', ['U', NODAL, (COMPONENT, 'U2')]],
        ...     ['Uzz', ['U', NODAL, (COMPONENT, 'U3')]],
        ...     ['UHor', ['UHor', NODAL, ()]],
        ...     ]

      Or like this:
        >>> varLabels = [
        ...     # [name for output, name in odbVarLabels
        ...     ['UVERT', 'U_3'],
        ...     # [name, [label, position, refinement]]
        ...     ['UHor', ['UHor', NODAL, ()]],
        ...     ['UNorth', ['U', NODAL, (COMPONENT,'U2')]],
        ...     ]
    
      Or simply like this:
        >>> varLabels = ["PST", "S_MIN_PRINCIPAL"]

      Any variable specified only by label (U_3, PST, S_MIN_PRINCIPAL in the examples above)
      must be contained in the dictionary cae_01.odbVarLabels.

      Ordering of the variables in the resulting csv file is retained from the
      order of varLabels.

    @param pointCoords: List of point coordiante tuples defining the path.

    @param shape: 'UNDEFORMED' or 'DEFORMED': Supposedly the path moves with the material 
      if you specify 'UNDEFORMED'. For 'UNDEFORMED', the value is extracted for a material
      point, that was sitting at the given position in the UNDEFORMED state.
      If you specify 'DEFORMED' the values are extracted for a particle sitting at the 
      given position in the 'UNDEFORMED' state. Needs to be checked. Whoever is sure, that 
      the above is true, should remove the questionmarks here...

    @param stepIndex: zero based step index, -1 is the last, -2 the second
      last ...
    @param frameList: list of zero based frame index numbers, frames that are
      not in the odb are silently ignored so frameList=range(1001) means all
      existing frames.
    @param odbName: If specified and if there is not already an odb open in the
      current viewport then open this odb and make it the displayed object in
      the current viewport. If there is an odb displayed in the current
      viewport then this argument is ignored and the displayed odb is used.
    @param tempFileName: The file abaqus view writes its results to. This is
      then parsed.
    @param logfile: Where to write diagnostic output?
     - If not specified create a new logfile according to the initOdb()
       function.
     - If None: no output.
     - If it's a string write to a file with that name.
     - Anything else is treated as the logfile to write to.

    @returns: A tuple (varList, valDict). The first variable is called
      "TRUE_DISTANCE" and contains the position along the path of the
      corresponding data point.
    
      varList: list of variable names as from the varLabels argument and       
    
      valDict: is a dictionary {frame number: data list}
      where data list is a list of list-per-point-lists, one per point in the
      pointCoords list. Each of those list-per-point lists finally contains
      the values in the order of varList.

      So valDict looks something like that with td being the "TRUE_DISTANCE" of
      the specific point along the path:
      {frame number: [ [td, var1, var2, ...], [td, var1, var2, ...], ...]}

      Data values that cannot be converted to floats yield a None in the
      resulting list of values.
    """
    global odbVarLabels

    (svp, odb) = initOdb(odbName)
    if logfile == "CreateDefaultLog":
        logfile = openDefaultLogFile()
    elif logfile == None:
        logfile = quietlog
    elif isinstance(logfile, basestring):
        logfile = LogFile(logfile, "w", scriptName=getScriptName())

    # path from points
    pth = session.Path(name="pointlistpath", type=POINT_LIST,
                       expression=pointCoords)
    if shape.upper()=='UNDEFORMED':
        shape = UNDEFORMED
    elif shape.upper()=='DEFORMED':
        shape = DEFORMED
    else:
        raise ValueError("shape can only be 'DEFORMED' or 'UNDEFORMED'.")

    # correct the step index
    nbSteps = len(odb.steps)
    while stepIndex<0:
        stepIndex += nbSteps
    # filter frameList: only valid frames
    validFrames = set(range(len(odb.steps.values()[stepIndex].frames)))
    for idx, frame in enumerate(frameList):
        while frame<0:
            frame += len(validFrames)
            frameList[idx] = frame
    validFrames.intersection_update(frameList)
    frameList = sorted(validFrames)

    nbFrames = len(frameList)
    logfile.write("will examine %d frames of step %s\n"
                  % (nbFrames, odb.steps.keys()[stepIndex]))

    # create a local copy of the global odbVarLabels dict
    odbVarLabels = dict(odbVarLabels)
    varList = list()
    for var in varLabels:
        if isinstance(var, basestring):
            # just a variable name to look up in odbVarLabels
            if var not in odbVarLabels:
                raise ValueError(
                    "Odb variable %s not recognized (in odbVarLabels dict)"
                    % (var))
            varList.append(var)
        elif (len(var)==2 and
              isinstance(var[0], basestring) and len(var[1])==3):
            # (name, odbVarLabel)-tuple
            varList.append(var[0])
            odbVarLabels[var[0]] = OdbVarLabel(*var[1])
        elif (len(var)==2 and
              isinstance(var[0], basestring) and isinstance(var[1], basestring)):
            # (name, odbVarLabel)-tuple
            try:
                odbVarLabels[var[0]] = odbVarLabels[var[1]]
            except KeyError:
                raise ValueError(
                    "Odb variable %s (for variable %s) not recognized (in"
                    " odbVarLabels dict)" % (var[1], var[0]))
            else:
                varList.append(var[0])
        else:
            raise ValueError(
                "The items of the parameter varLabels must either be"
                " strings specifying a variable from the global"
                " odbVarLabels dict or (name, odbVarLabel)-tuples with"
                " odbVarLabel being a (label, position, refinement)-tuple."
                "\n"
                "Cannot interpret the following item, it's being ignored:"
                "\n%s" % (var))
    del varLabels
    result = dict()

    for frameIdx in frameList:
        svp.odbDisplay.setFrame(step=stepIndex, frame=frameIdx)
        svp.odbDisplay.basicOptions.setValues(averageElementOutput=averageElementOutput)
        logfile.write ("processing frame number %d, extracting data to report"
                       " file %s\n"% (frameIdx, tempFileName))

        xyDataList = list()
        for var in varList:
            odbVar = odbVarLabels[var]
            svp.odbDisplay.setPrimaryVariable(
                variableLabel=odbVar.label,
                outputPosition=odbVar.position,
                refinement=odbVar.refinement)

            # xyName = "%s_%03d" % (var, svp.odbDisplay.fieldFrame[1])
            xyName = var
            xyData = session.XYDataFromPath(
                name=xyName, path=pth, includeIntersections=True,
                shape=shape, labelType=TRUE_DISTANCE)

            xyDataList.append(xyData)

        xyDataList = tuple(xyDataList)
        session.writeXYReport(
            fileName=tempFileName,
            appendMode=OFF, xyData=xyDataList)

        logfile.write("processing frame number %d, moving data from the"
                      " temporary report file %s to the output data"
                      " structures.\n"
                      % (frameIdx, tempFileName))

        tempFile = open(tempFileName)
        # skip first three lines
        tempFile.next(); tempFile.next(); tempFile.next()
        thisU = list() # result list for this one frame
        for line in tempFile:
            vals = line.strip().split()
            if len(vals)==0:
                continue
            try:
                u = map(float, vals)
            except ValueError:
                u = []
                for x in vals:
                    try:
                        u.append(float(x))
                    except ValueError:
                        u.append(None)
            thisU.append(u)
        tempFile.close()
        if _debug:
            logfile.write("frame %d, u: %s\n" % (frameIdx, thisU))

        result[frameIdx] = thisU
        logfile.flush()

    #-- fini
    try:
        os.remove(tempFileName)
    except OSError:
        pass

    varList.insert(0, "TRUE_DISTANCE")
    return varList, result
