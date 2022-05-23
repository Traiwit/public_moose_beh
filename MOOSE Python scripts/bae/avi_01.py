r"""Module for abaqusMacros scripts.

Description of basic functions
==============================

Sample abaqusMacros.py
----------------------

>>> from bae.avi_01 import *
>>>
>>> # parameter section
>>> projectPrefix = "KBGS2012"
>>> runVersion = "R04"
>>> mainOdbName = "../../2_Run/%s/myOdbName.odb" % runVersion[1:]
>>> frameListStatic = [ (2, i) for i in range(25) ]
>>> frameTest = (2,8)
>>>
>>> #------------------------
>>>
>>> dgRock = DGCtrlRemoveSeq("seqElsetsExcav")
>>> dgBolts = DGCtrlAddSeqFilter("Support",
>>>                              filterElsetsByElset=["DCB_BOLT","DSB_BOLT"])
>>>
>>> viewLongSection = View(
>>>     name='User-1', nearPlane=190.14, farPlane=267.08, width=16.934,
>>>     height=12.372, projection=PERSPECTIVE, cameraPosition=(204.58, 220.8,
>>>     206.17), cameraUpVector=(0.012998, 0.0099758, 0.99987), cameraTarget=(
>>>     177.38, 456.72, 204.17), viewOffsetX=1.0325, viewOffsetY=-1.1135,
>>>     autoFit=OFF)
>>>
>>> class LayerTunnelLongSectionRock(LayerContour):
>>>     def __init__(self, varToDisplay, translucencyFactor=1.0):
>>>         LayerContour.__init__(
>>>             self, "TunnelLongSectionRock",
>>>             varToDisplay=varToDisplay,
>>>             viewCut=ViewCutPlane(
>>>                 "LongSect",
>>>                 origin=[0,0,0], normal=[0.3,1.0,0.0], axis2=[0,1,0],
>>>                 planeCutPosition=481.0, planeAboveOnBelow=[1,0,0]),
>>>             displayGroupController=dgRock,
>>>             translucencyFactor=translucencyFactor)
>>>
>>> class LayerBoltsGrey(LayerSolidColour):
>>>     def __init__(self, fillColor='#222222', translucencyFactor=1.0):
>>>         LayerSolidColour.__init__(
>>>             self, "Bolts",
>>>             displayGroupController=dgBolts,
>>>             fillColor=fillColor,
>>>             edgeLineThickness=THICK,
>>>             translucencyFactor=translucencyFactor)
>>>
>>> sceneList = [
>>>     Scene(
>>>         "LongSection_PST",
>>>         frameList=frameListStatic,
>>>         layerList=[
>>>             LayerTunnelLongSectionRock(VarLogP_SDV4(2.0)),
>>>             LayerBoltsGrey()
>>>             ],
>>>         view=viewLongSection,
>>>         ),
>>>     ]
>>>
>>> # file and folder name templates from jobInfo.py data
>>> import os.path as osp
>>> outputDir = "%s_PNG/%s" % (
>>>         osp.splitext(osp.basename(mainOdbName))[0],
>>>         "%(sceneName)s" )  # sceneName will be inserted later
>>> pictureName = "%s_%s" % (
>>>         "%(projectPrefix)s_%(runVersion)s" % locals(), # expands now
>>>         "%(sceneName)s_%(varInfo)s_%(frameInfo)s" )    # expands later
>>>
>>> postJob = PostJob(mainOdbName=mainOdbName,
>>>                   projectPrefix=projectPrefix,
>>>                   outputDir=outputDir,
>>>                   pictureName=pictureName,
>>>                   sceneList=sceneList)
>>>
>>> def AAA_Testit():
>>>     postJob.test(frameListTest=[frameTest])
>>>
>>> def AAA_Makeit():
>>>     postJob.run()
>>>
>>> msg("abaqusMacros script reloaded and initialized...")
>>> if (__name__ == '__main__' ):
>>>     postJob.autorun()


To plot a legend only
---------------------

E.g. for the variable on a second layer of a certain scene. By default you
would only get the legend of the displayed variable of the first layer of a
scene. Plot a scene with that layer only with empty frameList.
 >>> ...
 >>> sceneList = [
 >>>     ...
 >>>     Scene("Boltslip_legend",
 >>>         frameList=[], layerList=[ LayerBoltSlip() ]),
 >>> ...


Record parameter of an interactive session
==========================================

There are functions provided to store current parameters of the current
interactive session (from the active viewport) in a way to be easily used for
the classes of this module avi_01. They will automatically appear in the macros
if you import avi_01 like so:
 >>> from bae.avi_01 import *

The function names all start with recordPara. The parameters will be appended
to the file recordedSessionPara.py in the current working directory.

If you don't want those function to appear in the macro list of the A/Viewer
then hide them with the L{hidden} function:
 >>> recordParaView = hidden(recordParaView)

Currently implemented:
 - recordParaAddSeparator: Adds just a long line of hash-marks to visually
   separate sections in the recorded file.
 - recordParaView: Stores aspect ratio (imageSize) of the current viewport
   and view (camera angles and so forth).
 - recordParaLayer: Stores currently active layer with display variable and
   colour scale properties (for contour-plot), display colour (for
   LayerSolidColour), colour mapping (for LayerMultiColour, i.e. for material
   regions, elsets, ...), status variable, viewcut, translocency, type of edge
   display, ...

The commands stored in the file recordedSessionPara.py can be copied into the
abaqusMacros.py script.

It's advisable and partly required to rename some of the items.

In order to become fully functional the following pieces should be added to
your abaqusMacros.py file:
 >>> sceneList = [
 >>>     Scene(
 >>>        "my_scene_name",
 >>>         frameList=[(2,i) for i in range(30, 68)],
 >>>         layerList=[
 >>>             layerXXX,  # should be renamed...
 >>>             ],
 >>>         view=viewXXX,  # should be renamed...
 >>>         imageSize=(1280, 977),  # take from recordedSessionPara.py
 >>>         ),
 >>>     ]
 >>> postJob = PostJob(
 >>>     mainOdbName="path/to/my_projects.odb",
 >>>     sceneList=sceneList)
 >>> 
 >>> def AAA_Testit():
 >>>     frameTest = (2,40)
 >>>     postJob.test(frameListTest=[frameTest])
 >>> 
 >>> msg("abaqusMacros script reloaded and initialized...")
 >>> if (__name__ == '__main__' ):
 >>>     postJob.autorun()

Recording in multilayer mode
----------------------------

The recordPara... functions work in single layer mode as well as in multi layer
mode. In multi layer mode the function recordParaLayer stores the currently
active layer. To store all layers you have to activate each, then store it,
for all layers, one by one.

Note: Before creating the first layer in the A/Viewer session you should set
the view to "front view" otherwise the coordinate system will be inconsistent
between multilayer mode and single layer mode.
"""

__version__ = "1.20"

_version_history = r"""\
Versions:
=========

0.01 GP: start of version numbering
         - no PART-1-1. in sequence elset lists of ...postdata.py
         - ...postdata.py automatically being found near the odb
0.02 GP: using log_01 module,
         no "PART-1-1" in displaygroup controllers initialState param any more,
         typo "acis" -> "axis",
0.03 GP: added DGCtrlStaticSomeSurfaces
0.04 GP: added vpSize argument to Scene, hidden function decorator
0.05 GP: cleaned up file and directory name generation, removed vpSize
         argument, automatic viewport size, added Scene.plotLegend()
         DGCtrlAddSeqFilter filter by element and by elset
0.06 GP: merge the ViewCutPlane variants. Add getOpposite() method.
         fixed plotLegend: doesn't change view anymore.
0.07 GP: added postJob.test, postJob.autorun, added sample in docstring
         fixed: default initialState for DGCtrlAddSeqFilter
         added: statusVariable in LayerContour to suppress elements based on
         result values.
0.08 GP: plotLegend with other than mainlayer possible
0.08b GP: fixed some doc string issues
0.08c GP: added some more docs
0.08d GP: improved example abaqusMacros.py: file names, testFrame
0.09 GP: added ViewCutMulti
0.10 GP: added renderBeamProfile
0.11 GP: prevent cae crash by computeOrder=EXTRAPOLATE_COMPUTE_AVERAGE in
         LayerContour.setup()
0.12 GP: added renderShellThickness (this seems to crash Abaqus CAE sometimes)
0.12b GP: added doku, parameter check
0.13 GP: incompatible change: renamed all DisplayGroupController starting with
         DgCtrl... to be consitent with the rest: DGCtrl...:
         DGCtrlStaticSomeNodes, DGCtrlAddSeqFilter.
0.14 GP: Added default view None for Scene class. Handy for plotting the
         legend only. Added docu.
         Added followDeformation option to some viewcuts.
         Added DGCtrlRemoveSeqFilter class
0.15 AF: adding param edgeColorFillShade to class LayerSolidColour and class
         LayerContour change in Scene.setup() that layer.odbFileName is
         condensed to absolute path in case that odbFileName is defined as
         relative path from working directory
0.16 GP: automatically convert all sequence elsets to uppercase
0.17 GP: deformed argument to two more layers, more todo notes
0.18 GP: frameForSetup argument for Layer can be a (stepNo, frameNo)-tuple
         default frameForSetup more robust now
0.19 GP: changed some defaults in predefined VarToDisplay Smin, Smax.
0.20 GP: added LayerContour()-argument visibleEdges (and replaced instance
         variable self.lineStyle by self.visibleEdges)
0.21 GP: added conMin argument to VarSMin
0.22 GP: added ViewWideAngle class for awesome views
1.01 GP: Switched to version 1, API changes must now be upwards compatible.
         Improved and tested ViewWideAngle class, removed some VarToDisplay
         child classes. Frame idenitier in png-name (default) switched from
         SyFxxx to Fyxxx for multistep analyses
1.02 GP: Use basicOptions->sectionResults=USE_BOTTOM_AND_TOP for shells if
         renderShellThickness is on. (I.e. data from different section points
         for top and bottom side of the shell)
1.03 GP: Added: VarLogP, VarRER_SDV12
1.04 GP: added explicit file name extension because Abaqus seems to add one
         only if there is no dot in the filename.
1.05 GP: added synonym VarToDisplay.name for VarToDisplay.varDisplayName
         (because everything else has a name attribute instead)
1.06 GP: changed default frameInfo string from _042_ to _F42_
1.06 FR: LayerSolidColour has no attribute edgeColorWireHide that allows to
         set the edge color in wireframe modus. This may be needed if you want
         to display the wireframe of a 2D 4-node quad-mesh as a coordinate
         grid.
         Default value is black "#000000".
1.07 don't know anymore
1.08 GP: incompatible change! stepNo now refers to the position in the steps
         sequence in the odb. This is consistent with the current meaning in
         bae.odb_02. Earlier stepNo was the number part in the step name.
         Removed FrameList.getStepNameFrameNbIterator()
1.09 GP: added imageType "VRML"
1.10 GP added LayerMultiColour
1.11 GP added frameInfoNameTemplate-argument to PostJob or Scene objects can
         now be a function as well.
1.12 GP added edgeLineStyle attribute and renderBeamProfiles in LayerContour
1.13 GP added deformedVariable argument to Layer classes
1.14 GP added: frameList can contain negative frame numbers. Until now
        DGCtrlRemoveSeq did not work with frameList = [(2,-1)].
1.15 GP fixed: projection (PARALLEL or PERSPECTIVE) in view definition works
1.16 FR added: contourStyle=UNIFORM is now the standard contour style.
        varToDisplay now allows to set the contourStyle.
1.17 GP fixed: argument statusVariable=None did not switch off the use of the
        statusVariable but left the default: in deformed plotstate
        statusVariable STATUS was used. Fixed: argument statusVariable=None
        now means: don't switch off elements because of any value of any
        variable.
1.18 GP changed: default statusVariable=None, was statusVariable="STATUS"
        which caused Cae to crash if STATUS was not in the odb as is usually
        the case in our models. (This is Valds idea and wish.)
1.19 VL added contourType (QUILT) to LayerContour
1.20 GP fixed DGCtrlAddRemoveSeq to correctly add and remove stope and fill
        elsets even if extraction and fill is in one step,
"""

_todo_ = r"""\
BUGS:
 - filename extension .png seems to be missing if there is a dot in the first
   part of the filename.

OTHER:
 - for principal stress directions (and in fact all other vector plots) Sven
   implemented a trick incorporating an elset. Have a look at
   http://boab2/Wiki/AbaqusViewerPostProcessing%282f%29PlotPrincipalDirections.html
   http://boab2/Wiki/AbaqusPreProcessing%282f%29CreateVolumeSetsEquallySpaced.html
 - do we need Scene.setup() independently from Scene.plot()?
   . __init__ must be separate because we don't have all objects working yet
   . In class Layer we have setup once at the beginning and update for each
     frame but for Scene it's different.
   . In class Scene do we need separate setup for easy overloading? Where would
     be the distinction between setup and plot? Should setup be called by
     plot(), then?
   . There's been Scene_setup functions in the old abaqusMacros...
 - AF recommends to take out majority of subclasses from VarToDisplay because
   these are project-specific subclasses (example: LOGP is not always 'SDV4')
   but to put them in a project-specific 'abaqusMacros_PRJbase.py' that then can
   easily being imported to any other script
 - AF wishes for a method 'apply()' in class View(), such that an instance of
   View can be easily being applied to the currentLayer in currentViewport of an
   interactive Abq/Vi session, example in 'abaqusMacros.py':
    >>> myViewInstance = View(...)
    >>> def applyVIEW_myViewname(): myViewInstance.apply()

   the above defined function will be visible and executable in Abq/Vi GUI
   MacroManager
   coding therefore within class View():
    >>>    def applyView(self):
    >>>        ### apply this view to current layer in current viewport
    >>>        sv = session.viewports[session.currentViewportName]
    >>>        self.setupDisplayView()
    >>>        sv.view.setValues(session.views[self.name])
 - GP: do viewCut objects need a name? General arbitrarily defined planes need
   one and isosurfaces need one. But what about ViewCutCoordPlane and
   ViewCutMulti?
"""


# from abaqusConstants import *
# specify all constants to make possible syntax checking feature
from abaqusConstants import \
    ON, OFF, TRUE, FALSE, \
    SINGLE, OVERLAY, \
    SOLID, GRADIENT, \
    CURRENT, ODB_REGIONS, INVARIANT, COMPONENT, \
    WIREFRAME, FILLED, HIDDEN, SHADED, \
    PERSPECTIVE, PARALLEL, \
    DEFORMED, UNDEFORMED, CONTOURS_ON_DEF, CONTOURS_ON_UNDEF, \
    SYMBOLS_ON_DEF, SYMBOLS_ON_UNDEF, \
    ALL, NONE, FEATURE, EXTERIOR, FREE, \
    DASHED, DOT_DASH, DOTTED, \
    VERY_THIN, THIN, MEDIUM, THICK, \
    RESULTANT, VECTOR_COMPONENT, \
    MODEL_SIZE, SCREEN_SIZE, \
    EXTRAPOLATE_COMPUTE_AVERAGE, EXTRAPOLATE_AVERAGE_COMPUTE, \
    USE_BOTTOM, USE_BOTTOM_AND_TOP, \
    SPECIFY, SPECTRUM, UNIFORM, LOG, USER_DEFINED, DEFAULT_COLORS, \
    PNG, TIFF, SVG, \
    ALL_NODES, ALL_ELEMENTS, ALL_SURFACES, DEFAULT_MODEL, EMPTY_LEAF, \
    PLANE, ISOSURFACE, \
    NODAL, INTEGRATION_POINT, ELEMENT_FACE, ELEMENT_NODAL, ELEMENT_CENTROID, \
    WHOLE_ELEMENT, WHOLE_REGION, WHOLE_PART_INSTANCE, WHOLE_MODEL, \
    GENERAL_PARTICLE, \
    BANDED, QUILT, LINE, ISOSURFACE

from abaqus import session, RangeError
import os
from math import log10

#------------------------------------------------------------------------------
#--- diagnostic output, debugging

from bae.log_01 import log as logfile, msg
logfile.toFile("abaqusMacros.log")

#{ Miscellaneous
#------------------------------------------------------------------------------
#--- hide functions

class _HideHelper(object):
    """class used solely by the hidden decorator function to hide functions
    in the CAE Macro Manager.
    """
    def __init__(self, f):
        self.f = f
    def call(self, *args, **kwargs):
        return self.f(*args, **kwargs)
def hidden(f):
    """decorator function to hide functions in the CAE Macro Manager.

    In your abaqusMacros.py declare a hidden function not to appear in the
    CAE Macro Manager list like so:
     >>> @hidden
     >>> def myVerySecretFunction(...):
     >>>     ...
    """
    a = _HideHelper(f)
    return a.call
hidden = hidden(hidden)  # hide the hidden function as well

#------------------------------------------------------------------------------
#--- temporaliy change state (and later reset to old)
class SaveChangeState(object):
    """Save current values of certain attributes of an object, change those
    attributes and later reset to the old state.

    Works only for objects with a setValues method that takes keyword arguments
    exactly matching the corresponding attributes of the object itself.
    I.e. object.setValues(attr=value) is used to set object.attr to value.

    Usage:
     >>> saveVpAnnotationOptions = SaveChangeState(
     >>>     sv.viewportAnnotationOptions, triad=OFF,legend=ON )
     >>> ... do something with the new state ...
     >>> saveVpAnnotationOptions.restore()
     >>> ... now the former state of sv.viewportAnnotationOptions is reset

    Currently used by Scene.plotLegend() to temporarily change some plotting
    states based on what Scene.setup() prescribed. It is being reset to old in
    order to not interfer with a possibly following Scene.plot().
    """
    def __init__(self, obj, **kwargs):
        """
        @param obj: object of which the attributes are to be saved and later
            restored. obj must have a setValues method.
        @param kwargs: attributes and their (temporally limited) new values
        """
        self.saved = dict( (attrName, getattr(obj, attrName))
                           for attrName in kwargs )
        self.obj = obj
        obj.setValues(**kwargs)

    def restore(self):
        """restore the saved values of the object"""
        self.obj.setValues(**self.saved)

#} end of Miscellaneous

#------------------------------------------------------------------------------
#{ functions in the Macro Manager menu

def startup_ViewportOptions(backgroundColor='#B5B5B5', vpSize=(320.0, 240.0),
                            pngImageSize=(1280, 1024)):
    """Setting basic viewport options.
    This function has to be called once before anything can be plotted.

    @param backgroundColor:
    other popular values: '#999999', '#FFFFFF'
    @param vpSize: viewport size (width, height)-tuple

    @param pngImageSize:
    other popular values: (int(1280*1.5), int(1024*1.5))
    """
    sv  = session.viewports[session.currentViewportName]
    svo = sv.odbDisplay

    # be sure to operate in SINGLE mode
    sv.setValues(displayMode=SINGLE)

    #--- some all time options
    session.printOptions.setValues(vpDecorations=OFF)

    # setbasic viewportoptions
    sv.viewportAnnotationOptions.setValues(
        triad=OFF,legend=OFF, title=OFF, state=OFF, annotations=OFF,
        compass=OFF)
    session.graphicsOptions.setValues(
        backgroundStyle=SOLID,backgroundColor=backgroundColor)
    # session.graphicsOptions.setValues(backgroundColor='#FFFFFF')

    # this is for gradient
    # session.graphicsOptions.setValues(backgroundStyle=GRADIENT)
    # session.graphicsOptions.setValues(
    #     backgroundColor='#999999',backgroundBottomColor='#333333')

    # this is required since 6.10 to allow individual settings on each layer
    sv.setValues(fieldOutputLayers=CURRENT)
    sv.setValues(animationLayers=CURRENT)
    sv.setValues(plotOptionsLayers=CURRENT)
    sv.setValues(plotStateLayers=CURRENT)
    sv.setValues(viewManipLayers=CURRENT)

    #--- some defaults
    svo.commonOptions.setValues(translucency=OFF)
    svo.basicOptions.setValues(translucencySort=ON)

    svo.deformedShapeOptions.setValues(
        visibleEdges=NONE,deformationScaling=UNIFORM,
        uniformScaleFactor=1.0)
    svo.commonOptions.setValues(visibleEdges=FEATURE)

    # outside limits colour: cae default
    # svo.contourOptions.setValues(outsideLimitsMode=SPECTRUM)

    # outside limits colour: BE default for standard RAINBOW
    svo.contourOptions.setValues(
        outsideLimitsMode=SPECIFY, outsideLimitsAboveColor='#E10000',
        outsideLimitsBelowColor='#0000DB',
        contourStyle=UNIFORM)

    # outside limits colour: BE default for reverse
    # svo.contourOptions.setValues(
    #     outsideLimitsMode=SPECIFY, outsideLimitsBelowColor='#E10000',
    #     outsideLimitsAboveColor='#0000DB')

    svo.basicOptions.setValues(featureAngle=80)
    svo.basicOptions.setValues(noResultsColor='#CCCCCC')

    # Averaging element output at nodes
    svo.basicOptions.setValues(
        averageOnlyDisplayed=TRUE,
        averagingThreshold=75,
        # # needed for min and max principal values of S, otherwise: crash
        # # this computeOrder=EXTRAPOLATE_AVERAGE_COMPUTE is the cae default
        computeOrder=EXTRAPOLATE_COMPUTE_AVERAGE,
        regionBoundaries=ODB_REGIONS,
        useRegionBoundaries=True,
        )

    # setting the image size
    sv.setValues(width=vpSize[0], height=vpSize[1], origin=(0,-0.5*vpSize[1]))
    sv.restore()

    session.pngOptions.setValues(imageSize=pngImageSize)
    session.printOptions.setValues(reduceColors=FALSE)
    session.printOptions.setValues(vpBackground=OFF)

    # set anti aliasing of line off
    session.graphicsOptions.setValues(antiAlias=ON)

    return sv, svo
#} end of functions in the Macro Manager menu


#{ Top level: PostJob and Scene
#------------------------------------------------------------------------------
#--- PostJob
class PostJob(object):
    """Class to save all information about on job that is typically run at
    a particular invocation of abaqus viewer and may consist of a couple of
    scenes.

    @ivar projectPrefix: Argument from the constructor. May be used in file
       names and such.
    @ivar postData: A dictionary {odb name: postData}. For example:
       seqElsDict = self.postData["2_Run/02/supiDupiJob.odb"]["seqElsetsExcav"]

       There are postData atttributes in each L{layer<Layer>} object that get
       initialized from this dictionary, i.e.
       layer.postData = postJob.postData[layer.odbFileName]
    """
    def __init__(self, mainOdbName=None, projectPrefix=None, sceneList=[],
                 frameNames=None, frameInfoTemplate=None,
                 outputDir=r"%(projectPrefix)s_PNG/%(sceneName)s",
                 pictureName=(r"%(projectPrefix)s_%(sceneName)s_%(varInfo)s"
                              r"_%(frameInfo)s%(fileNameExtension)s"),
                 imageSize = (1280, 1024), imageType=PNG):
        """
        @param mainOdbName: Default odbFileName for all layers. It may as
            well be specified for each layer separately as an argument to
            Layer.__init__().
            Can be omitted if an odb is already loaded in A/Cae.
        @param projectPrefix: A name prefix for this project. By default will
            be used to construct folder and file names. See the name and
            pictureName arguments of Scene.__init__()
            If not given will be assumed from the the odb name.
        @param sceneList: List of L{Scene} objects.
        @param frameNames: Usually not specified! A dict {stepname: list of
            frame names}; The frame at index FrNb in the list states the frame
            name for frame FrNb in the odb. This name is added to the end of
            the file name of each frame image like FRAMENAME in the following
            image name:
            BE-M03-SUB02-N33000_STANDARD_UMAG_0_50m_0_00m_FrNb_FRAMENAME.png
            Optional. Can be overridden by a corresponding argument to the
            Scene-constructor. If not specified (neither for the PostJob nor
            for the Scene it's being taken from the external postData file
            corresponding to the odb of the Scene's main layer.

        @param frameInfoTemplate: Template for the frameInfo part of the
            picture name. If None some heuristics will be used to create a
            sensible frameInfo. May be overriden by a corresponding argument to
            Scene.__init__(). May contain the following tags:
             - "%(stepNo)d": step number, the first step in the odb has got
               stepNo=1.
               converted to integer.
             - "%(frameNo)02d": current frame number (an integer)
             - "%(frameName)s": current frame name taken from the frameNames
               dictionary.

            Can also be a function that takes three arguments stepNo, frameNo
            and frameName. E.g.:
             >>> def frameInfoTemplate(st, fr, na):
             >>>     return "%d%03d_%s" % (st, fr, na)

            For restart analyses use something like this (We only have step 2
            and 3, the latter sometimes referred to as 1):
             >>> def frameInfoTemplate(stepNb, frameNb, frameName):
             >>>     if stepNb==2:
             >>>         return "F2%03d_%s" % (frameNb, frameName)
             >>>     else:
             >>>         return "F3%03d_%s" % (frameNb, frameName)

        @param outputDir: Name template for the subdirectory of all files of
            one scene. May be overriden by a corresponding argument to
            Scene.__init__(). May contain the following tags:
             - "%(projectPrefix)s": argument to PostJob.__init__
             - "%(sceneName)s": argument to Scene.__init__
        @param pictureName: Name template for picture files (without the
            directory part). May be overriden by a corresponding argument to
            Scene.__init__(). May contain the following tags:
             - "%(projectPrefix)s": argument to PostJob.__init__
             - "%(sceneName)s": argument to Scene.__init__
             - "%(viewCutName)s": name attribute of the viewCut
             - "%(varDisplayName)s": varDisplayName attribute of the variable
               to be displayed in the first layer of the current scenes
               layerList.
             - "%(varInfo)s": variable name together with colour scale limits
             - "%(frameInfo)s": this is what results from frameInfoTemplate
             - "%(fileNameExtension)s": the file name extension corresponding to
               the given value of imageType.
               Note: Abaqus seems to append a sensible file name extension
               only if there is no dot in the filename.
        @param imageSize: canvas size in pixels for images. A tuple (width,
            height).

            To get an appropriate imageSize: Resize the CAE viewport window
            with the mouse and look up the logged width and height values in
            the file abaqus.rpy (default size: 320x240). Then get an imageSize
            tuple with the same aspect ratio. Default: (1280, 1024)

            Note: This parameter does not determine the displayed size on the
            screen. It does however change the aspect ratio of the displayed
            viewport.

        @param imageType: Usually not specified! TIFF, SVG or PNG. Or "VRML".
            Note that only "VRML" is a string!
        """
        if mainOdbName is None:
            try:
                mainOdbName = session.odbs.keys()[0]
            except IndexError:
                msgText = ("mainOdbName argument must be given if no odb is"
                           " loaded yet.")
                msg("ERROR: %s" % msgText)
                raise ValueError(msgText)
        self.mainOdbName = mainOdbName

        if projectPrefix is None:
            # assume projectPrefix from the odb name
            projectPrefix = os.path.basename(mainOdbName).split("_")[0]
        self.projectPrefix = projectPrefix
        self.sceneList = sceneList

        # default values for all Scenes:
        self.frameNames = frameNames
        self.frameInfoTemplate = frameInfoTemplate
        self.outputDir = outputDir
        self.pictureName = pictureName
        self.imageSize = imageSize
        self.imageType = imageType

        # initialize reference to this PostJob in each of the scene objects
        # and transfer postJob default attributes to all scenes that don't have
        # those attributes defined on their own.
        for scene in sceneList:
            scene.postJob = self
            # set self's defaults for some scene attributes
            for attrName in ["frameNames", "frameInfoTemplate",
                             "outputDir", "pictureName",
                             "imageSize", "imageType"]:
                if getattr(scene, attrName) is None:
                    setattr(scene, attrName, getattr(self, attrName))

        # odb name: postData dictionary
        # i.e. seqElsDict = self.postData["supiDupiJob.odb"]["seqElsetsExcav"]
        #
        # There are postData atttributes in each layer object that get
        # initialized from this dictionary here:
        # layer.postData = postJob.postData[layer.odbFileName]
        #
        # The dictionary will be filled in Scene.setup() from external file
        # named odbFileName.replace(".odb", "_postData.py")
        #
        self.postData = dict()

    def run(self):
        """Plot all images as defined for the post job.
        """
        for currentScene in self.sceneList:
            msg("Processing scene %s..." % currentScene.name, debugLevel=1)
            currentScene.setup()
            msg("... scene %s setup() done." % currentScene.name, debugLevel=1)
            currentScene.plot()
            msg("... scene %s plot() done." % currentScene.name, debugLevel=1)
            currentScene.plotLegend()
            msg("... scene %s plotting legend done."
                % currentScene.name, debugLevel=1)

    def test(self, frameListTest=None):
        """Plot just the one image for the last frame of the first scene
        in the scene list.
        """
        # only first scene
        saveSceneList = self.sceneList
        testScene = self.sceneList[0]
        self.sceneList = [testScene, ]

        # only frameListTest
        saveFrameList = testScene.frameList
        if frameListTest is None:
            testScene.frameList = testScene.frameList[-1:]
        else:
            testScene.frameList = frameListTest

        self.run()

        # restore original values
        testScene.frameList = saveFrameList
        self.sceneList = saveSceneList

    def autorun(self):
        """Convenience function to be called if abaqusMacros is invoked as
        stand-alone script. At the end of the abaqusMacros.py script you would
        typically have:
         >>> msg("abaqusMacros script reloaded and initialized...")
         >>> if (__name__ == '__main__' ):
         >>>     postJob.autorun()
        """

        import viewerModules  # don't know if this is needed
        msg("stand alone call...")

        session.Viewport(
            name='Viewport: 1', origin=(0.0, 0.0), width=320.,height=240.)
        session.viewports['Viewport: 1'].makeCurrent()
        session.viewports['Viewport: 1'].maximize()

        msg("loading odb...")
        o1 = session.openOdb(name=self.mainOdbName)
        session.viewports['Viewport: 1'].setValues(displayedObject=o1)

        startup_ViewportOptions()
        self.run()

#------------------------------------------------------------------------------
#--- Scene

class Scene(object):
    """Class to save information about a scene.
    contains all information about picture names, sequence, layers etc.

    Subclassing should be considered in the following cases:
     - overload createOutputDir() and/or createFileName() methods for more
       control on where results go.
     - to be continued...

    Here is how to get Freds awesome "Darkroom" Scene. You have to use the
    layer class LayerDark for the rock (showing PST) in a scene of type
    SceneDarkRoom:
     >>> varPSTDark = VarToDisplay(
     >>>     'SDV1',INTEGRATION_POINT,(),'PST',
     >>>     0.3, 1E-4, 10, 'DarkRoomSpectrum', '', '', "%2.2f", LOG)
     >>>
     >>> class LayerDark(LayerContour):
     >>>     def __init__(self, name, viewCut=None,
     >>>                  displayGroupController=None,
     >>>                  translucencyFactor=1.0):
     >>>         LayerContour.__init__(
     >>>             self, name, varToDisplay=varPSTDark, viewCut=viewCut,
     >>>             displayGroupController=displayGroupController,
     >>>             translucencyFactor=translucencyFactor)
     >>>     def setup(self, sv):
     >>>         LayerContour.setup(self, sv)
     >>>         # modify light
     >>>         self.sv.lightOptions.lights[1].setValues(enabled=True)
     >>>
     >>> class SceneDarkRoom(Scene):
     >>>     def plot(self):
     >>>         graphicsOptions = SaveChangeState(
     >>>             session.graphicsOptions,
     >>>             backgroundColor='#000000',
     >>>             backgroundBottomColor='#D0D0D0',
     >>>             translucencyMode=2)
     >>>         session.Spectrum(name="DarkRoomSpectrum",
     >>>                          colors =('#ADADAD', '#262626', ))
     >>>         Scene.plot(self)
     >>>         graphicsOptions.restore()

    @ivar currentFrameNo:   current frame displayed (internal use only)
    @ivar currentStepNo:   current step displayed (internal use only)
    @ivar layerNamesList: a list of the names of the layers in the layerList
       argument to the contructor
    @ivar mainLayer: The first layer in the layerList.
    @ivar sv: session.viewport, will be initialized by the setup() method.

    @note: There are additional instance variables like the arguments to the
       constructor.
    """

    def __init__(self, name, frameList, layerList, view=None,
                 frameNames=None, frameInfoTemplate=None,
                 outputDir=None, pictureName=None,
                 imageSize=None, imageType=None):
        """
        @param name: scene name to identify the layer and to generate a
            sub-folder containing the images
        @param frameList: frame numbers (a FrameList object)
        @param layerList: List of L{Layer} instances (actually of its subclasses
            like L{LayerContour}, L{LayerSolidColour} or such). Those layers
            are used for overlay plot. The first layer in the list is the main
            layer: That's the layer that is used to look up the variable name
            and frame names for the images. (So the variable frameNames only
            has to be defined in the postData.py file for the odb of the main
            layer.)
        @param view: a View object that is used for the view settings
        @param frameNames: Usually not specified (None). Otherwise a dict
            {stepname: list of frame names}; The frame at index FrNb in the
            list states the frame name for frame FrNb in the odb. This name can
            be used in constructing the file name of each frame image like
            FRAMENAME in the following image name:
            BE-M03-SUB02-N33000_STANDARD_UMAG_0_50m_0_00m_FrNb_FRAMENAME.png

            If not specified then take the value as specified by the frameNames
            argument of PostJob.__init__(). If not specified for the L{PostJob}
            object either then take the default from the external postData file
            associated with the odb of the main layer.

        @param frameInfoTemplate: Template for the frameInfo part of the
            picture name. If None will be taken from L{postjob<PostJob>}.
            If still None some heuristics will be used to create a sensible
            name. May contain the following tags:
             - "%(stepNo)d": step number, the first step in the odb has got
               stepNo=1.
               converted to integer.
             - "%(frameNo)02d": current frame number (an integer)
             - "%(frameName)s": current frame name taken from the frameNames
               dictionary.

            Can also be a function that takes three arguments stepNo, frameNo
            and frameName. E.g.:
             >>> def frameInfoTemplate(fr, st, na):
             >>>     return "%d%03d_%s" % (fr, st, na)

            For restart analyses use something like this (We only have step 2
            and 3, the latter sometimes referred to as 1):
             >>> def frameInfoTemplate(stepNb, frameNb, frameName):
             >>>     if stepNb==2:
             >>>         return "F2%03d_%s" % (frameNb, frameName)
             >>>     else:
             >>>         return "F3%03d_%s" % (frameNb, frameName)

        @param outputDir: Usually specified on the postjob level! Name template
            for the subdirectory of all files of one scene. If not specified
            defaults to that of the L{postJob<PostJob>}. May contain the
            following tags:
             - "%(projectPrefix)s": argument to PostJob.__init__
             - "%(sceneName)s": argument to Scene.__init__
        @param pictureName: Usually specified on the postjob level! Name
            template for picture files (without the directory part). If not
            specified defaults to that of the L{postJob<PostJob>}. May contain
            the following tags:
             - "%(projectPrefix)s": argument to L{PostJob.__init__}
             - "%(sceneName)s": argument to L{Scene.__init__}
             - "%(viewCutName)s": name attribute of the L{viewCut<ViewCut>}
             - "%(varDisplayName)s": varDisplayName attribute of the variable
               to be displayed in the first layer of the current scenes
               layerList.
             - "%(varInfo)s": variable name together with colour scale limits
             - "%(frameInfo)s": this is what results from frameInfoTemplate
             - "%(fileNameExtension)s": the file name extension corresponding
               to the given value of imageType.
               Note: Abaqus seems to append a sensible file name extension
               only if there is no dot in the filename.
        @param imageSize: Usually specified on the postjob level! Canvas size
            in pixels for images. A tuple (width, height). If not specified
            defaults to that of the L{postJob<PostJob>}.

            To get an appropriate imageSize: Resize the CAE viewport window
            with the mouse and look up the logged width and height values in
            the file abaqus.rpy (default size: 320x240). Then get an imageSize
            tuple with the same aspect ratio. Default: (1280, 1024)

            Note: This parameter does not determine the displayed size on the
            screen. It does however change the aspect ratio of the displayed
            viewport.

        @param imageType: Usually specified on the postjob level! TIFF, SVG or
            PNG or "VRML". Note that only "VRML" is a string!
            If not specified defaults to that of the L{postJob<PostJob>}.

        @Note: The viewport size will be automatically adapted to match the
            imageSize as closely as possible.
        """
        self.name                = name            
        self.frameList           = frameList

        assert len(layerList), \
            "There must be at least one layer in the layerList!"
        self.layerList           = layerList
        self.layerNamesList      = [layer.name for layer in layerList]
        self.mainLayer           = layerList[0]

        self.view                = view

        # The following attributes will later be overridden by PostJob.__init__
        # if they got the default value None here
        self.frameNames          = frameNames
        self.frameInfoTemplate   = frameInfoTemplate

        self.outputDir           = outputDir
        self.pictureName         = pictureName

        self.imageSize           = imageSize
        self.imageType           = imageType

        # height of the title bar of the viewport window
        # this value is important to determine the viewport dimensions matching
        # the desired image size
        # SHOULD BE A FLOAT! (there won't be explicit type conversion to float)
        self.titleBarHeigth = 5.0

        self.sv = None

        # store self as scene in all layers
        for layer in self.layerList:
            layer.scene = self

    def createOutputDir(self):
        """Return the path to the output files of this scene. Make sure the
        directory exists.

        Needs self.setup() to be run first. In order to initialize
        self.pictureNameInfo

        This is meant to be overloaded by derived classes for more control on
        the output directory.
        """
        pathForResults = self.outputDir % self.pictureNameInfo
        if not (os.access(pathForResults,os.F_OK)):
            os.makedirs(pathForResults)
        return pathForResults

    def createFileName(self, stepNo, frameNo):
        """Return a sensible file name for the picture for the current frame.

        Have a look at self.setupFileNameInfo() if you want to adapt filename
        creation. self.setupFileNameInfo() creates some parts of the filename
        during scene setup (only once, not for each frame again)

        Needs self.setup() to be run first. In order to initialize
        self.frameNames and self.pictureNameInfo

        This is meant to be overloaded by derived classes for more control on
        the output file names.
        """

        try:
            frameName = self.frameNames[stepNo][frameNo]
        except (KeyError, IndexError):
            frameName = ""
        frameInfoData = {
            "stepNo": stepNo,
            "frameNo": frameNo,
            "frameName": frameName,
            }

        self.pictureNameInfo["stepNo"] = stepNo
        self.pictureNameInfo["frameNo"] = frameNo
        self.pictureNameInfo["frameName"] = frameName
        if isinstance(self.frameInfoTemplate, basestring):
            self.pictureNameInfo["frameInfo"] = (
                self.frameInfoTemplate % frameInfoData)
        else:
            self.pictureNameInfo["frameInfo"] = self.frameInfoTemplate(
                *(frameInfoData[i] for i in ("stepNo", "frameNo", "frameName")))

        name = self.pictureName % self.pictureNameInfo

        return name

    def setup(self):
        """Does the initialization of the scene, the viewport and so on.

        Initializes all the layers and so on. Does what switchLayersOn()
        did in recent versions.
        """
        import visualization
        import xyPlot
        import displayGroupOdbToolset as dgo

        # initialize sv, self.sv
        sv = session.viewports[session.currentViewportName]
        self.sv = sv

        #----------------------------------------------------------------------
        #--- set viewport size
        width = sv.currentWidth
        height = sv.currentHeight

        newwidth = width
        newheight = (self.imageSize[1]*width/self.imageSize[0]
                     + self.titleBarHeigth )
        try:
            sv.setValues(width=newwidth, height=newheight)
        except RangeError, err:
            # Note: width did not change, so it must be the height that is wrong
            # maxheight: -1 to guarantee that it's not too large
            maxheight = float(str(err.args[0]).rsplit(None,1)[-1]) - 1
            msg("Changing viewport size from %gx%g to %gx%g failed."
                " Max height: %g"
                % (width, height, newwidth, newheight, maxheight),
                debugLevel=1)
            newwidth *= ( (maxheight-self.titleBarHeigth)
                          / (newheight-self.titleBarHeigth) )
            newheight = maxheight
            sv.setValues(width=newwidth, height=maxheight)

        msg("Changed viewport size from %gx%g to %gx%g."
            % (width, height, newwidth, newheight), debugLevel=1)

        #----------------------------------------------------------------------
        #--- check odb for each layer, initialize layer.postData
        for layer in self.layerList:

            # set default odbFileName to each layer if needed
            # this is the first time odbFileName is needed!
            if layer.odbFileName is None:
                layer.odbFileName = self.postJob.mainOdbName

            # check if all required odbs are already loaded
            # Note: therefore, reduce/condense relative paths to absolute paths
            layer.odbFileName = os.path.abspath(layer.odbFileName).replace(
                "\\","/")
            if not session.odbs.has_key(layer.odbFileName):
                session.openOdb(name=layer.odbFileName)

            # assign postData dict for current layer
            # either from the self.postJob.postData - repository
            # or read this data from external file
            try:
                layer.postData = self.postJob.postData[layer.odbFileName]
            except KeyError:
                # which postData file
                postDataFileName = layer.odbFileName.replace(
                    ".odb", "_postData.py")

                # load/parse postData from external file
                layer.postData = dict()
                execfile(postDataFileName, layer.postData)

                # store in the self.postJob.postData - repository
                self.postJob.postData[layer.odbFileName] = layer.postData

        #----------------------------------------------------------------------
        #--- create all layers

        # switch to SINGLE mode
        sv.setValues(displayMode=SINGLE)
        # delete all layers
        for name in sv.layers.keys():
            del sv.layers[name]

        msg('Switch layers on...')

        # set a defined view at the beginning of the setup
        sv.view.setValues(session.views['Front'])

        # create all layers
        msg("Create all layers ...")
        refLayer=self.layerList[0]

        for layer in self.layerList:

            # associate layer with odb
            odb = session.odbs[layer.odbFileName]
            sv.setValues(displayedObject=odb)

            ########## most probably not needed anymore:
            # # set a scenario if a mainPart is defined
            # # otherwise display the default model
            # svo = sv.odbDisplay
            # if layer.mainPart:
            #     setsToDisplay=(layer.mainPart)
            #     leaf = dgo.LeafFromElementSets(elementSets=(setsToDisplay))
            #     svo.displayGroup.replace(leaf=leaf)
            # else:
            #     leaf = dgo.Leaf(leafType=DEFAULT_MODEL)
            #     svo.displayGroup.replace(leaf=leaf)

            # create layer
            if layer==refLayer:
                sv.Layer(layer.name)

                # set layer transformation matrix
                LTM = (1.,0.,0.,0., 0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1.)
                # apply layer transformation matrix
                sv.layers[layer.name].view.setLayerTransform(layerTransform=LTM)

                refLayer=layer
            else:
                sv.Layer(layer.name,refLayer.name)

        # switch to OVERLAY mode (with all newly created layers)
        sv.setValues(visibleLayers=self.layerNamesList,
                     viewManipLayers=ALL, displayMode=OVERLAY)

        #----------------------------------------------------------------------
        #--- layer setup
        renderShellThickness = any(
            layer.renderShellThickness for layer in self.layerList)
        for layer in self.layerList:

            # some diagnostic output
            msg("Setup Layer %s ..." % layer.name, end="", debugLevel=0)
            # ...end this line here only for debugging output channels
            msg(" starting.", debugLevel=1)

            # call the layers setup function, also makes the layer active
            layer.setup(sv)

            # workaround for a CAE bug: renderShellThickness seems to act
            # globally therefor we switch it for all layers here
            sv.odbDisplay.basicOptions.setValues(
                renderShellThickness=renderShellThickness,
            )

            # ...add line prefix only to debugging output channels
            msg("Setup Layer %s " % layer.name, end="", debugLevel=1)
            msg("finished.", debugLevel=0)

        # make reference layer active
        refLayer.makeCurrent()

        #----------------------------------------------------------------------
        #--- update frameList: convert negative frame numbers
        newFrameList = list()
        for (stepNo,frameNo) in self.frameList:
            if frameNo<0:
                currentStepName = refLayer.odb.steps.keys()[stepNo-1]
                nbOfFrames = len(refLayer.odb.steps[currentStepName].frames)
                while frameNo<0:
                    frameNo -= nbOfFrames
            newFrameList.append((stepNo,frameNo))
        self.frameList = newFrameList

        #----------------------------------------------------------------------
        #--- create: self.frameNames
        if self.frameNames is None:
            # if not specified explicitly for the scene (usually the case)
            # it's being taken from external postData file
            try:
                self.frameNames = self.mainLayer.postData["frameNames"]
            except KeyError:
                msg("WARNING: frameNames not defined in postData for %s"
                    % self.mainLayer.odbFileName)
                self.frameNames = dict()

        #----------------------------------------------------------------------
        #--- create: self.frameInfoTemplate
        if self.frameInfoTemplate is None:
            # are there many steps or just one?

            # first try to guess from postData frameNames
            relevantSteps = self.frameNames.keys()
            if len(relevantSteps)==1:
                lastFrame = len(self.frameNames[relevantSteps[0]])

            # second try to guess from displaygroup controller
            if len(relevantSteps)==0:
                try:
                    seqElsetsDict = self.mainLayer.dgCtrl.seqElsetsDict
                except AttributeError:
                    pass
                else:
                    relevantSteps = seqElsetsDict.keys()
                    if len(relevantSteps)==1:
                        lastFrame = len(seqElsetsDict[relevantSteps[0]])

            # third try to guess from frameList
            if len(relevantSteps)==0:
                relevantSteps = set(s for s,f in self.frameList)
                if len(relevantSteps)==1:
                    lastFrame = max( f for s,f in self.frameList )

            # template for the frame string part of the file name
            # step and frame number or frame number only
            if len(relevantSteps) == 1:
                # supposedly the whole sequence is in one step
                nDigits = max(2, int(log10(lastFrame)+1))
                self.frameInfoTemplate = 'F%%(frameNo)0%dd' % nDigits
            else:
                self.frameInfoTemplate = 'F%(stepNo)s%(frameNo)03d'

            if len(self.frameNames):
                self.frameInfoTemplate += "_%(frameName)s"

        #----------------------------------------------------------------------
        #--- create: self.pictureNameInfo
        self.pictureNameInfo = dict()
        self.pictureNameInfo["projectPrefix"] = self.postJob.projectPrefix
        self.pictureNameInfo["sceneName"] = self.name

        try:
            viewCutName = self.mainLayer.viewCut.name
        except AttributeError:
            self.pictureNameInfo["viewCutName"] = ""
        else:
            self.pictureNameInfo["viewCutName"] = viewCutName

        try:
            varToDisplay = self.mainLayer.varToDisplay
        except AttributeError:
            self.pictureNameInfo["varDisplayName"] = ""
            self.pictureNameInfo["varInfo"] = ""
        else:
            self.pictureNameInfo["varDisplayName"] = varToDisplay.name
            self.pictureNameInfo["varInfo"] = varToDisplay.varInfo

        # note: Abaqus seems to append a sensible file name extension only if
        # there is no dot in the filename.
        if self.imageType==PNG:
            self.pictureNameInfo["fileNameExtension"] = ".png"
        elif self.imageType==TIFF:
            self.pictureNameInfo["fileNameExtension"] = ".tiff"
        elif self.imageType==SVG:
            self.pictureNameInfo["fileNameExtension"] = ".svg"
        elif self.imageType=="VRML":
            self.pictureNameInfo["fileNameExtension"] = ".wrl"
        else:
            self.pictureNameInfo["fileNameExtension"] = ""

    def plot(self):
        """taken from the old createSequence() function

        Needs self.setup() to be run first. To initialize everything properly.
        """
        sv = session.viewports[session.currentViewportName]

        msg('Starting picture sequence...')

        # select the user views
        if self.view is not None:
            self.view.setupDisplayView()

        # switch overlay on and switch to selected view using the 1st layer
        self.mainLayer.makeCurrent()
        if self.view is not None:
            # projection seems to be required separately
            sv.view.setProjection(projection=self.view.para["projection"])
            sv.view.setValues(session.views[self.view.name])

        # create/set directory according to the scene name
        pathForResults = self.createOutputDir()

        # file size info
        session.pngOptions.setValues(imageSize=self.imageSize)
        session.tiffOptions.setValues(imageSize=self.imageSize)

        # loop over all steps,frames in the odb
        msg("Start loop ....")
        for cnt, (stepNo,frameNo) in enumerate(self.frameList):

            #
            # define current options for the layers
            #
            for layer in self.layerList:
                layer.makeCurrent()
                layer.updateFrame(stepNo, frameNo)

            # reference main layer to set variable limit etc in file name
            name=self.createFileName(stepNo,frameNo)

            # setup a file name
            msg("Printing: %s Frame: %d Count: %d" % (name,frameNo,cnt))
            fileName = os.path.join(pathForResults, name)
            # create image from step
            # imageType may be TIFF, PNG, SVG
            if self.imageType=="VRML":
                session.writeVrmlFile(fileName=fileName, compression=0,
                                      canvasObjects=(sv, ))
            else:
                session.printToFile(fileName=fileName, format=self.imageType,
                                    canvasObjects=(sv, ))

    def plotLegend(self, imageSize=(240, 320), fileName=None,
                   legendTitle=False, legendBox=False,
                   layerIndex=0):
        """plot a legend for the variable on the current main layer

        Needs self.setup() to be run first. To initialize everything properly.
        plotLegend() can be called before or after plot().

        @param imageSize: number of pixels in the resulting image. The viewport
            size will be 0.25 this number in each direction.
        @param fileName: of the resulting legend image. Defaults to an
            ordinary image name with the frameinfo replaced by "legend"
        @param legendTitle: boolean whether to show the Abaqus variable name
            (No way to adapt that I'm afraid)
        @param legendBox: boolean whether to plot a frame around the legend
        @param layerIndex: index in the layer list of the layer for which the
            legend is to be plotted. If you specify a nondefault value here
            then you should also specify the filename argument, otherwise the
            filename of the legend will very likely be misleading.
        """
        if self.imageType=="VRML":
            return

        msg('Printing legend...')

        sv = session.viewports[session.currentViewportName]

        # switch off all but the main (first) layer
        saveVpLayers = SaveChangeState(
            sv, visibleLayers=self.layerNamesList[layerIndex:layerIndex+1],
            viewManipLayers=ALL, displayMode=OVERLAY)
        svo = self.layerList[layerIndex].makeCurrent()

        # create/set directory according to the scene name
        pathForResults = self.createOutputDir()

        # fileName
        if fileName is None:
            saveFrameInfoTemplate = self.frameInfoTemplate
            self.frameInfoTemplate = "legend"
            fileName=self.createFileName(stepNo=0,frameNo=0)
            self.frameInfoTemplate = saveFrameInfoTemplate
        fileName = os.path.join(pathForResults, fileName)

        # legend on, objects invisible
        legendTitle = {True: ON, False: OFF}[legendTitle]
        legendBox = {True: ON, False: OFF}[legendBox]
        saveVpAnnotationOptions = SaveChangeState(
            sv.viewportAnnotationOptions,
            triad=OFF,legend=ON, title=OFF, state=OFF, annotations=OFF,
            compass=OFF,
            legendTitle=legendTitle, legendBox=legendBox)
        saveVisibility = SaveChangeState(
            svo.commonOptions, translucencyFactor=0, visibleEdges=NONE)

        # change view further to the right
        # (should come before viewport size change)
        saveView = SaveChangeState(
            sv.view, viewOffsetX=sv.view.viewOffsetX+3.0*sv.view.width)

        # viewport size
        # (should come after view change)
        saveVpSize = (sv.currentWidth, sv.currentHeight)
        sv.setValues(width=imageSize[0]/4, height=imageSize[1]/4)

        # image size
        session.pngOptions.setValues(imageSize=imageSize)
        session.tiffOptions.setValues(imageSize=imageSize)

        # imageType may be TIFF, PNG, SVG
        msg("Printing: %s" % (fileName))
        session.printToFile(fileName=fileName, format=self.imageType,
                            canvasObjects=(sv, ))

        #--- clean up (restore order matters!)

        # restore viewport size
        sv.setValues(width=saveVpSize[0], height=saveVpSize[1])

        # restore view
        saveView.restore()

        # legend, visibility back to old
        saveVpAnnotationOptions.restore()
        saveVisibility.restore()

        # switch back on all layers
        saveVpLayers.restore()

#} end of Top level: PostJob and Scene


#------------------------------------------------------------------------------
#{ Layer

class Layer(object):
    """Base class for all kinds of layers.
    Holds information how to display an individual layer.

    Each subclass can provide a setup and an update method. The setup method
    is called once in the beginning and the update method once for each frame.

    @param name: Name to identify the layer. Does not occur in the output in
           any way.
    @param odbFileName: filename of the odb that is loaded into the layer
       defaults to the mainOdbName of the L{PostJob} object
    @param viewCut: a L{ViewCut} object
    @param displayGroupController: None means show all allways

    @ivar name: used to identify the layer in CAE's layer manager
    @ivar odbFileName: path to the odb to be displayed on this layer
    @ivar viewCut: might be a L{ViewCutNone} instance or of another subclass of
       L{ViewCut}
    @ivar dgCtrl: instance of (a subclass of) L{DisplayGroupController}, like
       L{DGCtrlStatic}, L{DGCtrlRemoveSeq} or such.

    @ivar postData: A dictionary with the contents of the ..._postData.py file
       corresponding to self.odbFileName. For example:
       seqElsDict = self.postData["seqElsetsExcav"]

       This atttribute is only created in the L{Scene.setup} function and not
       available earlier.
    """
    def __init__(self, name, odbFileName=None,
                 viewCut=None, displayGroupController=None):
        self.name = name
        self.odbFileName = odbFileName

        if viewCut is None:
            self.viewCut = ViewCutNone()
        else:
            self.viewCut = viewCut

        if displayGroupController is None:
            self.dgCtrl = DGCtrlStatic(DEFAULT_MODEL)
        else:
            self.dgCtrl = displayGroupController

    def setup(self, sv):
        """
        Make sure this function is called before any update() or makeCurrent
        method.

        This should also be invoked by all derived sub-classes.
        """

        # set viewport, odb-object, viewport.odbDisplay-object
        self.sv = sv
        self.odb = session.odbs[self.odbFileName]
        self.svo = self.makeCurrent()

        msgTemplate = "Layer %s setup(): %s finished." % (self.name, "%s")

        self.dgCtrl.setup(self)
        msg(msgTemplate % "self.dgCtrl.setup()", debugLevel=2)
        self.viewCut.setup(self.svo)
        msg(msgTemplate % "self.viewCut.setup()", debugLevel=2)

    def makeCurrent(self):
        self.sv.setValues(layerOffset=0, currentLayer=self.name,
                          viewManipLayers=ALL, displayMode=OVERLAY)
        svo = self.sv.odbDisplay
        return svo


class LayerSolidColour(Layer):
    """Plot the configuration with a specified (solid, fixed) colour,
    view cut activated.
    """
    def __init__(self, name, odbFileName=None,
                 viewCut=None, displayGroupController=None,
                 fillColor='#888888',
                 visibleEdges=NONE,
                 edgeLineStyle=SOLID,
                 edgeLineThickness=VERY_THIN,
                 edgeColorFillShade='#666666',
                 edgeColorWireHide='#000000',
                 renderBeamProfiles=False, beamScaleFactor=1.0,
                 renderShellThickness=False, shellScaleFactor=1.0,
                 translucencyFactor=1.0,
                 statusVariable=None,
                 statusMinimum=None, statusMaximum=None, statusInsideRange=None,
                 deformed=True, deformedVariable=None
                 ):
        """
        @param name: Name to identify the layer. Does not occur in the output
           in any way.
        @param viewCut: a ViewCut object
        @param displayGroupController: None means show all allways
        @param fillColor: an RGB value coded as three byte hexadecimal code,
           e.g. '#888888'
        @param edgeLineStyle: SOLID (default), DASHED, DOT_DASH, DOTTED
        @param edgeLineThickness: VERY_THIN (default), THIN, MEDIUM, THICK
        @param edgeColorFillShade: any color, defaults to grey
        @param visibleEdges: ALL, EXTERIOR, FEATURE, FREE, NONE.
           NONE can be used only when renderStyle=SHADED
        @param renderBeamProfiles: True, False (If True it's suggested to also
           set edgeLineStyle=DOTTED, edgeLineThickness=VERY_THIN.)
        @param beamScaleFactor: scale factor for the beam cross section.
           Applies only if renderBeamProfiles=True.
        @param renderShellThickness: True, False
        @param shellScaleFactor: scale factor for the shell thickness. Applies
           only if renderShellThickness=True.
        @param translucencyFactor: between 0 (transparent/invisible) and
           1 (opaque)
        @param statusVariable: field output variable for filtering element
           display based on a status criteria. None if not used. If "STATUS"
           then use the STATUS-variable, which is the A/cae-default in deformed
           plot state.
        @param statusMinimum: A Float specifying the minimum result value to
           be considered for element removal. I.e. this is the maximum allowed
           value for an element in order to be displayed.
        @param statusMaximum: A Float specifying the maximum result value to
           be considered for element removal. I.e. this is the minimum allowed
           value for an element in order to be displayed.
        @param statusInsideRange: A Boolean utilized when both statusMinimum
           and statusMaximum are given. Elements will be removed when they
           contain values between the minimum and maximum, inclusive, when
           true. Elements will be removed when they contain values outside
           of the minimum and maximum, exclusive, when false.
        @param deformed: True ...on deformed configuration, False ...on
           undeformed
        @param deformedVariable: deformation variable (default is U) of type
           L{VarToDisplay}.
        """

        Layer.__init__(self, name,odbFileName,viewCut,displayGroupController)
        self.fillColor = fillColor
        self.edgeColorFillShade=edgeColorFillShade
        self.edgeLineStyle = edgeLineStyle
        self.edgeLineThickness = edgeLineThickness
        self.edgeColorWireHide = edgeColorWireHide

        self.mainPart = ''   # mainPart is what is placed on the single view
                    # if this is empty the WHOLEMODEL default will be placed

        if deformed:
            self.plotState = DEFORMED
        else:
            self.plotState = UNDEFORMED

        # session.viewports[name].odbDisplay.commonOptions. ...
        # visibleEdges: ALL, EXTERIOR, FEATURE, FREE, NONE
        #   ... NONE can be used only when renderStyle=SHADED
        if visibleEdges is None:
            visibleEdges = NONE
        self.visibleEdges = visibleEdges
        self.renderBeamProfiles = renderBeamProfiles
        self.beamScaleFactor = beamScaleFactor
        self.renderShellThickness = renderShellThickness
        self.shellScaleFactor = shellScaleFactor

        # renderStyle: WIREFRAME, FILLED, HIDDEN, SHADED
        self.renderStyle = SHADED

        if translucencyFactor<1.0:
            self.translucency=ON             # ON, OFF
            self.translucencyFactor=translucencyFactor
        else:
            self.translucency=OFF

        if statusVariable=="STATUS":
            self.statusVariable = varStatus
            self.statusMinimum = 0
            self.statusMaximum = 0.5
            self.statusInsideRange = True
        else:
            self.statusVariable = statusVariable
            self.statusMinimum = statusMinimum
            self.statusMaximum = statusMaximum
            self.statusInsideRange = statusInsideRange

        self.deformedVariable = deformedVariable

    def setup(self, sv):
        import displayGroupOdbToolset as dgo
        Layer.setup(self, sv)

        self.svo.commonOptions.setValues(
            visibleEdges=self.visibleEdges,renderStyle=self.renderStyle,
            fillColor=self.fillColor,
            edgeColorFillShade=self.edgeColorFillShade,
            edgeLineStyle=self.edgeLineStyle,
            edgeLineThickness=self.edgeLineThickness,
            edgeColorWireHide=self.edgeColorWireHide,
            )
        self.svo.display.setValues(plotState=(self.plotState, ))

        # render beam profiles
        if self.renderBeamProfiles:
            self.svo.basicOptions.setValues(
                renderBeamProfiles=ON,
                beamScaleFactor=self.beamScaleFactor,
                )
        else:
            self.svo.basicOptions.setValues( renderBeamProfiles=OFF )

        # render shell thickness
        if self.renderShellThickness:
            self.svo.basicOptions.setValues(
                renderShellThickness=ON,
                shellScaleFactor=self.shellScaleFactor,
                )
        else:
            self.svo.basicOptions.setValues( renderShellThickness=OFF )

        # translucency
        if self.translucency==ON:
            self.svo.commonOptions.setValues(
                translucency=ON, translucencyFactor=self.translucencyFactor)
            self.svo.basicOptions.setValues(translucencySort=ON)

        # modify light
        self.sv.lightOptions.lights[1].setValues(enabled=False)

        # status variable to restrict plotted elements based on result values
        if self.statusVariable is None:
            self.svo.setStatusVariable(useStatus=False)
            msg("... setStatusVariable(useStatus=False).",
                debugLevel=4)
        else:
            var = self.statusVariable
            msg("... setStatusVariable(var=%s)." % var.name,
                debugLevel=4)
            kwargs = dict(
                variableLabel=var.varToDisplay,
                outputPosition=var.varPosition,
                refinement=var.varRefinement,
                useStatus=True,
                applyStatusToUndeformed=bool(self.plotState==UNDEFORMED) )
            for k in ["statusMinimum","statusMaximum","statusInsideRange"]:
                v = getattr(self, k)
                if v is not None:
                    kwargs[k] = v
            self.svo.setStatusVariable(**kwargs)

        # deformed variable
        if self.plotState==DEFORMED and self.deformedVariable is not None:
            self.svo.setDeformedVariable(
                variableLabel=self.deformedVariable.varToDisplay)

        # make sure we're only displaying one single colour
        self.sv.disableMultipleColors()

    def updateFrame(self, currentStepNo, currentFrameNo):
        self.dgCtrl.update(currentStepNo, currentFrameNo)
        if (self.plotState == DEFORMED
            or isinstance(self.viewCut, ViewCutIsoSurface)
            or self.statusVariable is not None):
            # display the correct frame
            currentStepName = self.odb.steps.keys()[currentStepNo-1]
            try:
                f = self.odb.steps[currentStepName].frames[currentFrameNo]
            except IndexError:
                raise IndexError(
                    "%s, frame %d does not exist in the current odb %s"
                    % (currentStepName, currentFrameNo, self.odb.name))
            self.svo.setFrame(frame=f)


class LayerMultiColour(LayerSolidColour):
    """Plot the configuration with a specified multi-colour colour mapping,
    view cut activated.

    Primarily used to display material regions:
     >>> layer = LayerMultiColour(
     >>>     "Rocktype", viewCut=viewCutLong, colorMapType='Material',
     >>>     colorMap={"Volcanics": '#555577', "Sediment": '#A09050'})
    """
    def __init__(self, name, odbFileName=None,
                 viewCut=None, displayGroupController=None,
                 colorMapType='Material',
                 colorMap={},
                 visibleEdges=NONE,
                 edgeLineStyle=SOLID,
                 edgeLineThickness=VERY_THIN,
                 edgeColorFillShade='#666666',
                 edgeColorWireHide='#000000',
                 renderBeamProfiles=False, beamScaleFactor=1.0,
                 renderShellThickness=False, shellScaleFactor=1.0,
                 translucencyFactor=1.0,
                 statusVariable=None,
                 statusMinimum=None, statusMaximum=None, statusInsideRange=None,
                 deformed=True, deformedVariable=None
                 ):
        """
        @param name: Name to identify the layer. Does not occur in the output
           in any way.
        @param viewCut: a ViewCut object
        @param displayGroupController: None means show all allways
        @param colorMapType: possible values: 'Material','Section',
           'Element set', 'Averaging region', 'Element type', 'Internal set',
           'Display group'. For further options see in A/CAE:
           session.viewports['Viewport: 1'].colorMappings.keys()
        @param colorMap: A dictionary defining the colour mapping. For
           colorMapType='Material' keys are material names. Values are
           RGB colour codes in the form of a three byte hex value like
           '#FFFF00' (this would be yellow). If the value is None or '' then
           this region shall not be shown.
           Keys not specified will receive a default value from A/CAE.
        @param edgeLineStyle: SOLID (default), DASHED, DOT_DASH, DOTTED
        @param edgeLineThickness: VERY_THIN (default), THIN, MEDIUM, THICK
        @param edgeColorFillShade: any color, defaults to grey
        @param visibleEdges: ALL, EXTERIOR, FEATURE, FREE, NONE.
           NONE can be used only when renderStyle=SHADED
        @param renderBeamProfiles: True, False (If True it's suggested to also
           set edgeLineStyle=DOTTED, edgeLineThickness=VERY_THIN.)
        @param beamScaleFactor: scale factor for the beam cross section.
           Applies only if renderBeamProfiles=True.
        @param renderShellThickness: True, False
        @param shellScaleFactor: scale factor for the shell thickness. Applies
           only if renderShellThickness=True.
        @param translucencyFactor: between 0 (transparent/invisible) and
           1 (opaque)
        @param statusVariable: field output variable for filtering element
           display based on a status criteria. None if not used
        @param statusMinimum: A Float specifying the minimum result value to
           be considered for element removal. I.e. this is the maximum allowed
           value for an element in order to be displayed.
        @param statusMaximum: A Float specifying the maximum result value to
           be considered for element removal. I.e. this is the minimum allowed
           value for an element in order to be displayed.
        @param statusInsideRange: A Boolean utilized when both statusMinimum
           and statusMaximum are given. Elements will be removed when they
           contain values between the minimum and maximum, inclusive, when
           true. Elements will be removed when they contain values outside
           of the minimum and maximum, exclusive, when false.
        @param deformed: True ...on deformed configuration,
           False ...on undeformed
        @param deformedVariable: deformation variable (default is U) of type
           L{VarToDisplay}.
        """

        LayerSolidColour.__init__(
            self, name, odbFileName=odbFileName,
            viewCut=viewCut,
            displayGroupController=displayGroupController,
            visibleEdges=visibleEdges,
            edgeLineStyle=edgeLineStyle,
            edgeLineThickness=edgeLineThickness,
            edgeColorFillShade=edgeColorFillShade,
            edgeColorWireHide=edgeColorWireHide,
            renderBeamProfiles=renderBeamProfiles,
            beamScaleFactor=beamScaleFactor,
            renderShellThickness=renderShellThickness,
            shellScaleFactor=shellScaleFactor,
            translucencyFactor=translucencyFactor,
            statusVariable=statusVariable,
            statusMinimum=statusMinimum, statusMaximum=statusMaximum,
            statusInsideRange=statusInsideRange,
            deformed=deformed, deformedVariable=deformedVariable)
        self.colorMapType = colorMapType
        self.colorMap = colorMap

    def setup(self, sv):
        LayerSolidColour.setup(self, sv)
        self.sv.enableMultipleColors()
        self.sv.setColor(initialColor='#BDBDBD')
        cmap = self.sv.colorMappings[self.colorMapType]
        if self.colorMap:
            cmapDict = dict()
            for key, colour in self.colorMap.iteritems():
                if colour:
                    cmapDict[key] = (True, colour, 'Default', colour)
                else:
                    cmapDict[key] = (False,)
            cmap.updateOverrides(overrides=cmapDict)
        self.sv.setColor(colorMapping=cmap)


class LayerContour(Layer):
    """Contourplot on deformed or undeformed configuration
    """
    def __init__(self, name, varToDisplay,
                 odbFileName=None,
                 viewCut=None, displayGroupController=None,
                 frameForSetup=None,
                 visibleEdges=FEATURE,
                 edgeLineStyle=SOLID,
                 edgeLineThickness=VERY_THIN,
                 edgeColorFillShade='#666666',
                 renderBeamProfiles=False, beamScaleFactor=1.0,
                 renderShellThickness=False, shellScaleFactor=1.0,
                 translucencyFactor=1.0,
                 statusVariable=None,
                 statusMinimum=None, statusMaximum=None, statusInsideRange=None,
                 deformed=True, deformedVariable=None,
                 contourType=BANDED,
                 ):
        """
        @param name: Name to identify the layer. Does not occur in the output in
           any way.
        @param varToDisplay: a VarToDisplay object
        @param odbFileName: filename of the odb to display data from. *Only*
           specify if different from the mainOdbName of the postJob object.
        @param viewCut: a ViewCut object
        @param visibleEdges: ALL, EXTERIOR, FEATURE, FREE, NONE (can be used
           only when renderStyle=SHADED)
        @param edgeLineStyle: SOLID (default), DASHED, DOT_DASH, DOTTED
        @param edgeLineThickness: VERY_THIN (default), THIN, MEDIUM, THICK
        @param edgeColorFillShade: any color, defaults to grey
        @param displayGroupController: None means show all allways
        @param frameForSetup: Optional frame number or (step number, frame
           number)-tuple. This frame will be used during initialization of the
           layer. If not specified, the last frame in the odb will be taken for
           this purpose. If the given frame number is negative then the last
           frame (of the specified step) is taken.
           This option is required if the last frame (=default frameForSetup)
           does not contain the requested output variable. In this case the
           setup()-function would fail to switch to contour plot state resulting
           in DEFORMED or UNDEF plot state (objects all in green).
        @param renderShellThickness: True, False (WARNING: This seems to cause
           Abaqus Viewer to crash sometimes.)
        @param shellScaleFactor: scale factor for the shell thickness. Applies
           only if renderShellThickness=True.
        @param translucencyFactor: between 0 (transparent/invisible) and
           1 (opaque)
        @param statusVariable: field output variable for filtering element
           display based on a status criteria. None if not used. If "STATUS"
           then use the STATUS-variable, which is the A/cae-default in deformed
           plot state.
        @param statusMinimum: A Float specifying the minimum result value to
           be considered for element removal. I.e. this is the maximum allowed
           value for an element in order to be displayed.
        @param statusMaximum: A Float specifying the maximum result value to
           be considered for element removal. I.e. this is the minimum allowed
           value for an element in order to be displayed.
        @param statusInsideRange: A Boolean utilized when both statusMinimum
           and statusMaximum are given. Elements will be removed when they
           contain values between the minimum and maximum, inclusive, when
           true. Elements will be removed when they contain values outside
           of the minimum and maximum, exclusive, when false.
        @param deformed: True ...on deformed configuration, False ...on
           undeformed
        @param deformedVariable: deformation variable (default is U) of type
           L{VarToDisplay}
        @param contourType: BANDED or QUILT are the most used values.
           LINE, ISOSURFACE are possible as well.
        """
        Layer.__init__(self, name,odbFileName,viewCut,displayGroupController)
        self.varToDisplay = varToDisplay

        self.frameForSetup = frameForSetup

        self.mainPart = ''  # mainPart is what is placed on the single view
                       # if this is empty the WHOLEMODEL default will be placed

        # deformed contour plot: CONTOURS_ON_DEF or CONTOURS_ON_UNDEF
        if deformed:
            self.plotState = CONTOURS_ON_DEF
        else:
            self.plotState = CONTOURS_ON_UNDEF

        self.renderBeamProfiles = renderBeamProfiles
        self.beamScaleFactor = beamScaleFactor
        self.renderShellThickness = renderShellThickness
        self.shellScaleFactor = shellScaleFactor

        # renderStyle: WIREFRAME, FILLED, HIDDEN, SHADED
        self.renderStyle = SHADED
        # session.viewports[name].odbDisplay.commonOptions. ...
        # visibleEdges: ALL, EXTERIOR, FEATURE, FREE, NONE
        #   ... NONE can be used only when renderStyle=SHADED
        if visibleEdges is None:
            visibleEdges = NONE
        self.visibleEdges = visibleEdges
        self.edgeLineStyle = edgeLineStyle
        self.edgeLineThickness = edgeLineThickness
        self.edgeColorFillShade=edgeColorFillShade

        if translucencyFactor<1.0:
            self.translucency=ON             # ON, OFF
            self.translucencyFactor=translucencyFactor
        else:
            self.translucency=OFF

        if statusVariable=="STATUS":
            self.statusVariable = varStatus
            self.statusMinimum = 0
            self.statusMaximum = 0.5
            self.statusInsideRange = True
        else:
            self.statusVariable = statusVariable
            self.statusMinimum = statusMinimum
            self.statusMaximum = statusMaximum
            self.statusInsideRange = statusInsideRange

        self.deformedVariable = deformedVariable
        self.contourType = contourType

    def switchVarOn(self):
        """Set primary variable for contour plot.

        Usually called from self.setup(). Only if the displayGroup is empty
        this action must be postponed until some elements become visible.

        Don't do anything if the variable to be displayed is not in the set
        of available variables for the current step
        """
        msg("LayerContour.switchVarOn() start.", debugLevel=4)
        var = self.varToDisplay
        msg("... determine variables existing in current frame.", debugLevel=4)
        variablesInCurrentFrame = set(
            x[0] for x in self.svo.fieldVariables.variableList)
        if var.varToDisplay in variablesInCurrentFrame:
            msg("... setPrimaryVariable(var=%s)." % var.name,
                debugLevel=4)
            self.svo.setPrimaryVariable(
                variableLabel=var.varToDisplay, outputPosition=var.varPosition,
                refinement=var.varRefinement)
            msg("... set plotState: %s." % self.plotState, debugLevel=4)
            self.svo.display.setValues(plotState=(self.plotState, ))
            self.isVarSwitchedOn = True
        else:
            msg("... variable to display (%s) does not exist in current frame."
                % var.varToDisplay, debugLevel=4)

        msg("LayerContour.switchVarOn() end.", debugLevel=4)

    def setup(self, sv):
        import displayGroupOdbToolset as dgo
        import visualization

        msg("LayerContour.setup() for %s" % self.name, debugLevel=2)
        Layer.setup(self, sv)
        msg(" ... Layer.setup() finished.", debugLevel=2)

        # needed for min and max principal values of S, otherwise: crash
        self.svo.basicOptions.setValues(
            computeOrder=EXTRAPOLATE_COMPUTE_AVERAGE)

        if len(self.odb.steps)>0:
            if self.frameForSetup is None:
                for stepName, lastStep in self.odb.steps.items()[::-1]:
                    frameNo = len(lastStep.frames)-1
                    if frameNo>=0:
                        break
                self.frameForSetup = None  # to trigger initalization
            elif isinstance(self.frameForSetup, int):
                for stepName, lastStep in self.odb.steps.items()[::-1]:
                    if self.frameForSetup>=0:
                        frameNo = self.frameForSetup
                        # only consider given framenumber in the last step
                        self.frameForSetup = -1
                        if 0<=frameNo<len(lastStep.frames):
                            # found a valid frame number
                            break
                    else:
                        frameNo = len(lastStep.frames)-1
                        if frameNo>=0:
                            break
                self.frameForSetup = None  # to trigger initalization
            elif isinstance(self.frameForSetup, (list, tuple)):
                stepName = self.odb.steps.keys()[self.frameForSetup[0]-1]
                frameNo = self.frameForSetup[1]
                self.frameForSetup = None  # to trigger initalization

            if self.frameForSetup is None:
                try:
                    self.frameForSetup=self.odb.steps[stepName].frames[frameNo]
                except IndexError:
                    raise IndexError(
                        "Suggested frameForSetup (%s, frame %d) does not exist"
                        " in the current odb %s"
                        % (stepName, frameNo, self.odb.name))
                else:
                    msg(" ... identified (%s,%d) as frameForSetup."
                        % (stepName, frameNo), debugLevel=1)

            self.svo.setFrame(frame=self.frameForSetup)
            msg(" ... setFrame(frameForSetup=%d) finished."
                % self.frameForSetup.incrementNumber, debugLevel=2)

        self.svo.commonOptions.setValues(
            visibleEdges=self.visibleEdges,renderStyle=self.renderStyle,
            edgeLineStyle=self.edgeLineStyle,
            edgeLineThickness=self.edgeLineThickness,
            )

        # set contour options
        self.varToDisplay.updateContourColourScale(self.svo)
        self.svo.contourOptions.setValues(contourType=self.contourType)

        # render beam profiles
        if self.renderBeamProfiles:
            self.svo.basicOptions.setValues(
                renderBeamProfiles=ON,
                beamScaleFactor=self.beamScaleFactor)
        else:
            self.svo.basicOptions.setValues(renderBeamProfiles=OFF)

        # render shell thickness, use different section points for both sides
        # bug in CAE: renderShellThickness seems to act globally, so
        # we don't switch it here. In Scene.setup it's switched globally for
        # each layer of the scene.
        if self.renderShellThickness:
            self.svo.basicOptions.setValues(
                # renderShellThickness=ON,
                shellScaleFactor=self.shellScaleFactor,
                sectionResults=USE_BOTTOM_AND_TOP)
        else:
            self.svo.basicOptions.setValues(
                # renderShellThickness=OFF,
                sectionResults=USE_BOTTOM)

        # translucency
        if self.translucency==ON:
            self.svo.commonOptions.setValues(
                translucency=ON, translucencyFactor=self.translucencyFactor)
            self.svo.basicOptions.setValues(translucencySort=ON)
        msg(" ... setting some options finished.", debugLevel=2)

        # flag if the variable and contour plot is already switched on
        # (may not happen during setup if the display group is empty)
        self.isVarSwitchedOn = False

        # set contour var and plot state
        if self.dgCtrl.elemsDisplayed:
            self.switchVarOn()
            msg(" ... self.switchVarOn() finished.", debugLevel=2)
        else:
            msg(" ... self.switchVarOn() skipped, dgCtrl empty.", debugLevel=2)

        # modify light
        self.sv.lightOptions.lights[1].setValues(enabled=False)

        # status variable to restrict plotted elements based on result values
        if self.statusVariable is None:
            self.svo.setStatusVariable(useStatus=False)
            msg("... setStatusVariable(useStatus=False).",
                debugLevel=4)
        else:
            var = self.statusVariable
            msg("... setStatusVariable(var=%s)." % var.name,
                debugLevel=4)
            kwargs = dict(
                variableLabel=var.varToDisplay,
                outputPosition=var.varPosition,
                refinement=var.varRefinement,
                useStatus=True,
                applyStatusToUndeformed=bool(
                    self.plotState==CONTOURS_ON_UNDEF))
            for k in ["statusMinimum","statusMaximum","statusInsideRange"]:
                v = getattr(self, k)
                if v is not None:
                    kwargs[k] = v
            self.svo.setStatusVariable(**kwargs)

        # deformed variable
        if self.plotState==DEFORMED and self.deformedVariable is not None:
            self.svo.setDeformedVariable(
                variableLabel=self.deformedVariable.varToDisplay)

        msg("LayerContour.setup() for %s finished." % self.name, debugLevel=2)

    def updateFrame(self, currentStepNo, currentFrameNo):

        self.dgCtrl.update(currentStepNo, currentFrameNo)

        # display the correct frame
        currentStepName = self.odb.steps.keys()[currentStepNo-1]
        try:
            f = self.odb.steps[currentStepName].frames[currentFrameNo]
        except IndexError:
            raise IndexError(
                "%s, frame %d does not exist in the current odb %s"
                % (currentStepName, currentFrameNo, self.odb.name))
        self.svo.setFrame(frame=f)

        # possibly switch on correct variable and plot state
        if not(self.isVarSwitchedOn) and self.dgCtrl.elemsDisplayed:
            self.switchVarOn()


class LayerSymbolVector(Layer):
    """vector symbols on deformed or undeformed configuration

    @param name:
    @param varToDisplay: a VarToDisplay object
    @param viewCut: a ViewCut object (optional)

    @param symbolPosition: optional, can be NODAL, INTEGRATION_POINT,
       ELEMENT_FACE, ELEMENT_NODAL, ELEMENT_CENTROID, WHOLE_ELEMENT,
       WHOLE_REGION, WHOLE_PART_INSTANCE, WHOLE_MODEL, and GENERAL_PARTICLE
    @param symbolMax: If None will be taken from varToDisplay
    @param symbolMin:

    @param vectorColor: If None vectors will be coloured according to the
       colourscale associated with varToDisplay. Otherwise this serves as the
       colour string, e.g. "#333333"

    @param symbolDensity: A Float specifying the factor for randomized
       sampling. (0 is high, 2 is about half way and 4 is lowest)
    @param arrowScaleMode: MODEL_SIZE or SCREEN_SIZE
    @param arrowSymbolSize:

    @param displayGroupController: None means show all allways
    @param frameForSetup: OdbFrame object to be displayed during self.setup().
       If None then take the first frame of the first step.

    @param deformed: True ...on deformed configuration, False ...on undeformed
    @param deformedVariable: deformation variable (default is U) of type
       L{VarToDisplay}.
    """
    def __init__(self, name, varToDisplay,
                 odbFileName=None, viewCut=None,
                 symbolPosition=None,
                 symbolMax=None, symbolMin=None,
                 vectorColor=None,
                 arrowSymbolSize=2.0, arrowScaleMode=MODEL_SIZE,
                 symbolDensity=1.0,
                 displayGroupController=None, frameForSetup=None,
                 deformed=True, deformedVariable=None
                 ):
        Layer.__init__(self, name, odbFileName, viewCut, displayGroupController)
        self.varToDisplay = varToDisplay

        self.frameForSetup = frameForSetup

        self.mainPart = ''  # mainPart is what is placed on the single view
                       # if this is empty the WHOLEMODEL default will be placed

        if deformed:
            self.plotState=SYMBOLS_ON_DEF
        else:
            self.plotState=SYMBOLS_ON_UNDEF

        # special for symbol plots
        if symbolPosition is None:
            self.symbolPosition = varToDisplay.varPosition
        else:
            self.symbolPosition = symbolPosition

        if symbolMax is None:
            self.symbolMax = self.varToDisplay.conMax
        else:
            self.symbolMax = symbolMax
        if symbolMin is None:
            self.symbolMin = self.varToDisplay.conMin
        else:
            self.symbolMin = symbolMin

        self.vectorColor = vectorColor

        self.arrowSymbolSize = arrowSymbolSize
        self.arrowScaleMode = arrowScaleMode
        self.symbolDensity = symbolDensity

        self.deformedVariable = deformedVariable

    def switchVarOn(self):
        """Set primary variable for symbol plot.

        Usually called from self.setup(). Only if the displayGroup is empty
        this action must be postponed until some elements become visible.

        Don't do anything if the variable to be displayed is not in the set
        of available variables for the current step
        """
        var = self.varToDisplay
        variablesInCurrentFrame = set(
            x[0] for x in self.svo.fieldVariables.variableList)
        if var.varToDisplay in variablesInCurrentFrame:

            self.svo.setSymbolVariable(
                variableLabel=self.varToDisplay.varToDisplay,
                outputPosition=self.symbolPosition,
                vectorQuantity=RESULTANT)
            self.svo.display.setValues(plotState=(self.plotState, ))
            self.isVarSwitchedOn = True

    def setup(self, sv):
        import displayGroupOdbToolset as dgo
        import visualization
        Layer.setup(self, sv)

        if len(self.odb.steps)>0:
            if self.frameForSetup is None:
                self.frameForSetup = self.odb.steps.values()[-1].frames[0]
            self.svo.setFrame(frame=self.frameForSetup)

        self.svo.symbolOptions.setValues(
            vectorMaxValueAutoCompute=OFF, vectorMaxValue=self.symbolMax,
            vectorMinValueAutoCompute=OFF, vectorMinValue=self.symbolMin)

        if self.vectorColor is None:
            # set colourscale options
            self.svo.symbolOptions.setValues(
                vectorColorMethod=SPECTRUM,
                vectorColorSpectrum=self.varToDisplay.conSpec,
                vectorIntervalNumber=self.varToDisplay.conLev)
        else:
            # set fixed colour for symbols
            self.svo.symbolOptions.setValues(
                vectorColorMethod=UNIFORM,vectorColor=self.vectorColor)

        # symbol size and density
        self.svo.symbolOptions.setValues(
            arrowSymbolSize=self.arrowSymbolSize,
            arrowScaleMode=self.arrowScaleMode,
            symbolDensity=self.symbolDensity)

        # translucency, hide all
        self.svo.commonOptions.setValues(
            translucency=ON, translucencyFactor=0, visibleEdges=NONE,
            renderStyle=FILLED)
        self.svo.basicOptions.setValues(translucencySort=ON)

        # flag if the variable and symbol plot is already switched on
        # (may not happen during setup if the display group is empty)
        self.isVarSwitchedOn = False

        # set symbol var and plot state
        if self.dgCtrl.elemsDisplayed:
            self.switchVarOn()

        # modify light
        self.sv.lightOptions.lights[1].setValues(enabled=False)

        # deformed variable
        if self.plotState==DEFORMED and self.deformedVariable is not None:
            self.svo.setDeformedVariable(
                variableLabel=self.deformedVariable.varToDisplay)

    def updateFrame(self, currentStepNo, currentFrameNo):

        self.dgCtrl.update(currentStepNo, currentFrameNo)

        # display the correct frame
        currentStepName = self.odb.steps.keys()[currentStepNo-1]
        try:
            f = self.odb.steps[currentStepName].frames[currentFrameNo]
        except IndexError:
            raise IndexError(
                "%s, frame %d does not exist in the current odb %s"
                % (currentStepName, currentFrameNo, self.odb.name))
        self.svo.setFrame(frame=f)

        # possibly switch on correct variable and plot state
        if not(self.isVarSwitchedOn) and self.dgCtrl.elemsDisplayed:
            self.switchVarOn()


class LayerSymbolTensor(Layer):
    """tensor symbols on deformed or undeformed configuration

    @param name:
    @param varToDisplay: a VarToDisplay object
    @param viewCut: a ViewCut object

    @param symbolPosition: optional, can be NODAL, INTEGRATION_POINT,
       ELEMENT_FACE, ELEMENT_NODAL, ELEMENT_CENTROID, WHOLE_ELEMENT,
       WHOLE_REGION, WHOLE_PART_INSTANCE, WHOLE_MODEL, and GENERAL_PARTICLE
    @param symbolMax: If None will be taken from varToDisplay
    @param symbolMin:

    @param tensorQuantity: Possible values are ALL_PRINCIPAL_COMPONENTS,
       PRINCIPAL_COMPONENT, ALL_DIRECT_COMPONENTS, and DIRECT_COMPONENT.
       Abaqus default is ALL_PRINCIPAL_COMPONENTS.
    @param tensorRefinement:

    @param vectorQuantity: can be RESULTANT or VECTOR_COMPONENT. Abaqus default
       is RESULTANT.
    @param symbolDensity: A Float specifying the factor for randomized
       sampling. (0 is high, 2 is about half way and 4 is lowest)
    @param arrowScaleMode: MODEL_SIZE or SCREEN_SIZE
    @param arrowSymbolSize:

    @param displayGroupController: None means show all allways
    @param frameForSetup: OdbFrame object to be displayed during self.setup().
    If None then take the first frame of the first step.

    @param deformed: True ...on deformed configuration, False ...on undeformed
    @param deformedVariable: deformation variable (default is U) of type
       L{VarToDisplay}.

    @note: construction work ahead. Not ready for use!
    """
    def __init__(self, name, odbFileName, varToDisplay, viewCut,
                 symbolPosition,
                 symbolMax=None, symbolMin=None,
                 tensorQuantity=None, tensorRefinement=None,
                 vectorQuantity=None, vectorColorMethod = SPECTRUM,
                 arrowSymbolSize=2.0, arrowScaleMode=MODEL_SIZE,
                 symbolDensity=1.0,
                 displayGroupController=None, frameForSetup=None,
                 deformed=True, deformedVariable=None
                 ):
        Layer.__init__(self, name, odbFileName, viewCut, displayGroupController)
        self.varToDisplay = varToDisplay

        self.frameForSetup = frameForSetup

        self.mainPart = ''  # mainPart is what is placed on the single view
                       # if this is empty the WHOLEMODEL default will be placed

        if deformed:
            self.plotState=SYMBOLS_ON_DEF
        else:
            self.plotState=SYMBOLS_ON_UNDEF

        # special for symbol plots
        if symbolPosition is None:
            self.symbolPosition = varToDisplay.varPosition
        else:
            self.symbolPosition = symbolPosition

        if symbolMax is None:
            self.symbolMax = self.varToDisplay.conMax
        else:
            self.symbolMax = symbolMax
        if symbolMin is None:
            self.symbolMin = self.varToDisplay.conMin
        else:
            self.symbolMin = symbolMin

        self.tensorQuantity = tensorQuantity
        self.tensorRefinement = tensorRefinement

        self.vectorQuantity = vectorQuantity
        self.vectorColorMethod = vectorColorMethod

        self.arrowSymbolSize = arrowSymbolSize
        self.arrowScaleMode = arrowScaleMode
        self.symbolDensity = symbolDensity

        self.deformedVariable = deformedVariable

    def setup(self, sv):
        import displayGroupOdbToolset as dgo
        import visualization
        Layer.setup(self, sv)

        if len(self.odb.steps)>0:
            if self.frameForSetup is None:
                self.frameForSetup = self.odb.steps.values()[-1].frames[0]
            self.svo.setFrame(frame=self.frameForSetup)

        self.svo.commonOptions.setValues(
            visibleEdges=self.visibleEdges,renderStyle=self.renderStyle)

        if self.tensorQuantity:
            self.svo.setSymbolVariable(
                variableLabel=self.varToDisplay.varToDisplay,
                outputPosition=self.symbolPosition,
                tensorQuantity=self.tensorQuantity,
                refinement=self.tensorRefinement)
            self.svo.symbolOptions.setValues(
                tensorMaxValueAutoCompute=OFF, tensorMaxValue=self.symbolMax,
                tensorMinValueAutoCompute=OFF, tensorMinValue=self.symbolMin)
            self.svo.symbolOptions.setValues(tensorColorSpectrum='Rainbow',
                                             tensorIntervalNumber=12)
            # self.svo.symbolOptions.setValues(
            #     tensorColorMethod=UNIFORM,tensorColor='#FFFFFF')
        else:
            self.svo.setSymbolVariable(
                variableLabel=self.varToDisplay.varToDisplay,
                outputPosition=self.symbolPosition,
                vectorQuantity=self.vectorQuantity)
            self.svo.symbolOptions.setValues(
                vectorMaxValueAutoCompute=OFF, vectorMaxValue=self.symbolMax,
                vectorMinValueAutoCompute=OFF, vectorMinValue=self.symbolMin)
            if self.vectorColorMethod==SPECTRUM:
                self.svo.symbolOptions.setValues(
                    vectorColorMethod=SPECTRUM,
                    vectorColorSpectrum='Rainbow', vectorIntervalNumber=10)
            else:
                self.svo.symbolOptions.setValues(
                    vectorColorMethod=UNIFORM,vectorColor='#333333')

        self.svo.symbolOptions.setValues(arrowSymbolSize=self.arrowSymbolSize,
                                         arrowScaleMode=self.arrowScaleMode,
                                         symbolDensity=self.symbolDensity)

        var = self.varToDisplay
        self.svo.setPrimaryVariable(
            variableLabel=var.varToDisplay, outputPosition=var.varPosition,
            refinement=var.varRefinement)

        self.svo.display.setValues(plotState=(self.plotState, ))

        # set contour options
        self.varToDisplay.updateContourColourScale(self.svo)

        # translucency
        self.svo.commonOptions.setValues(
            translucency=ON, translucencyFactor=0, visibleEdges=NONE,
            renderStyle=FILLED)
        self.svo.basicOptions.setValues(translucencySort=ON)

        # modify light
        self.sv.lightOptions.lights[1].setValues(enabled=False)

        # deformed variable
        if self.plotState==DEFORMED and self.deformedVariable is not None:
            self.svo.setDeformedVariable(
                variableLabel=self.deformedVariable.varToDisplay)

    def updateFrame(self, currentStepNo, currentFrameNo):

        self.dgCtrl.update(currentStepNo, currentFrameNo)

        # display the correct frame
        currentStepName = self.odb.steps.keys()[currentStepNo-1]
        try:
            f = self.odb.steps[currentStepName].frames[currentFrameNo]
        except IndexError:
            raise IndexError(
                "%s, frame %d does not exist in the current odb %s"
                % (currentStepName, currentFrameNo, self.odb.name))
        self.svo.setFrame(frame=f)
#} end of Layer


#------------------------------------------------------------------------------
#{ View cuts

class ViewCut(object):
    """base class for view cuts

    Usage
    =====

    Pass an object of one of ViewCut's child classes to the viewcut argument
    when creating the corresponding Layer.


    Notes on subclassing
    ====================

    ViewCut objects are handled by Layer objects, i.e. Layer.setup() calls
    the setup() method of its ViewCut object.

    In general each class derived from ViewCut should provide an initViewerVC()
    method in order for this viewcut class to be usable for the ViewCutMulti
    class. initViewerVC() must accept one argument svo and return a tuple of
    viewcut names of the Abaqus viewer viewcut objects.

    In general derived classes do not need an own setup method. The base
    class's ViewCut.setup() method calls the derived class's
    self.initViewerVC().

    Each subclass must provide a name attribute for output purposes.
    """

    def __init__(self, name, planeCutPosition, planeAboveOnBelow,
                 followDeformation=False):
        self.name = name
        self.planeCutPosition = planeCutPosition
        self.planeAboveOnBelow = planeAboveOnBelow
        self.followDeformation = followDeformation

    def setCutPositions(self, viewerVC):
        r"""Sets the showModel[Above|On|Below]Cut and the position arguments
        of the Abaqus viewer viewcut object viewerVC.

        viewerVC is required to be the Abaqus viewer viewcut-object associated
        with self.

        @Note: Sets the showFreeBodyCut option to False (Setting to True would
        show resultant force and moment vectors on this viewcut plane).

        Setting to False seems to be required for multiple viewcuts.

        Manipulating this option for other than plane viewcuts
        (i.e. isosurfaces) might not be feasible.
        """
        viewerVC.setValues(
            showModelBelowCut=self.planeAboveOnBelow[2],
            showModelOnCut=self.planeAboveOnBelow[1],
            showModelAboveCut=self.planeAboveOnBelow[0],
            position=self.planeCutPosition,
            followDeformation=self.followDeformation,
            showFreeBodyCut=False)

    def setup(self, svo):
        """This method (or the corresponding of a derived class) will be called
        by Layer.setup()
        """
        namesOfCut = self.initViewerVC(svo)
        svo.setValues(viewCutNames=namesOfCut, viewCut=ON)


class ViewCutNone(ViewCut):
    """view cut off, display all
    """
    def __init__(self):
        self.name = "All"

    def setup(self, svo):
        svo.setValues(viewCut=OFF)


class ViewCutCoordPlane(ViewCut):
    """view cut perpendicular to a certain coordinate axis. Uses the predefined
    viewcuts X-Plane, Y-Plane or Z-Plane.

    @param name: name of the view cut
    @param axis: "x"=="NS", "y"=="EW" or "z"=="PLAN" (case insensitive)
    @param planeAboveOnBelow: something like [FALSE,TRUE,TRUE]
    @param followDeformation: True or False

    @Note: Only one viewcut based on one of the three predefined viewcuts can be
    used *on the same layer*. I.e. you cannot create two "x"-viewcuts with
    different positions and use them in a multiple viewcuts assembly (of type
    ViewCutMulti). Or use an "x" viewcut and it's opposite returned by
    self.getOpposite().
    """
    def __init__(self, name, axis, planeCutPosition, planeAboveOnBelow,
                 followDeformation=False):
        ViewCut.__init__(self, name, planeCutPosition, planeAboveOnBelow,
                         followDeformation)
        axis = axis.upper()
        if axis in ("X", "NS", 'X-PLANE'):
            self.viewerVCName = 'X-Plane'
        elif axis in ("Y", "EW", 'Y-PLANE'):
            self.viewerVCName = 'Y-Plane'
        elif axis in ("Z", "PLAN", 'Z-PLANE'):
            self.viewerVCName = 'Z-Plane'
        else:
            raise ValueError("viewcut axis '%s' not implemented." % axis)

    def initViewerVC(self, svo):
        viewerVC = svo.viewCuts[self.viewerVCName]
        self.setCutPositions(viewerVC)
        nameOfCut = viewerVC.name
        return (nameOfCut,)

    def getOpposite(self, name=None, showOnPlane=None):
        """return the viewcut that shows exactly what self does not show.

        @param name: name for the new viewcut. If not specified defaults to
            "anti_%s"%self.name
        @param showOnPlane: If defined can override the middle item of the
            planeAboveOnBelow-triple. I.e. you can force the new viewcut to
            show the cutting plane even though self shows it as well. Or force
            the new viewcut to not show it although self does not show it
            either. If specified must be either TRUE or FALSE. The default
            behaviour for the new viewcut is to show the plane if self does
            not and to not show it if self does.
        """
        switch = {FALSE:TRUE, TRUE:FALSE}
        planeAboveOnBelow = [switch[i] for i in self.planeAboveOnBelow]
        if showOnPlane is not None:
            planeAboveOnBelow[1] = showOnPlane
        if name is None:
            name="anti_%s" % self.name

        return ViewCutCoordPlane(
            name=name,
            axis=self.viewerVCName,
            planeCutPosition=self.planeCutPosition,
            planeAboveOnBelow=planeAboveOnBelow,
            followDeformation=self.followDeformation)


class ViewCutPlane(ViewCut):
    """ViewCut defined by an arbitrary plane.

    Can be definied by either
     - origin, normal and axis2 vectors (as in the Abaqus Viewer dialog box),
     - or by three points.
    """

    def __init__(self, name,
                 planeCutPosition, planeAboveOnBelow, followDeformation=False,
                 **kwargs):
        """
        @param name: name of the view cut
        @param planeCutPosition: offset between the origin and the plane
        @param planeAboveOnBelow: something like [FALSE,TRUE,TRUE]
        @param followDeformation: True or False
        @kwarg origin: point on the plane
        @kwarg normal: vector perpendicular to the plane
        @kwarg axis2: vector in the plane identifying the 2-axis
        @kwarg pointNormal: point on the line from the origin perpendicular
            to the plane
        @kwarg pointOnAxis2: point on the line from the origin into the
            direction of axis2
        """

        ViewCut.__init__(self, name, planeCutPosition, planeAboveOnBelow,
                         followDeformation)

        if all(k in kwargs for k in ("origin", "normal", "axis2")):
            self.origin = kwargs["origin"]
            self.normal = kwargs["normal"]
            self.axis2 = kwargs["axis2"]
        elif all(k in kwargs for k in ("origin","pointNormal","pointOnAxis2")):
            def subPoint(a,b):
                c=(a[0]-b[0],a[1]-b[1],a[2]-b[2])
                return c
            self.origin = kwargs["origin"]
            self.normal = subPoint(kwargs["pointNormal"],self.origin)
            self.axis2 = subPoint(kwargs["pointOnAxis2"],self.origin)
        else:
            raise ValueError(
                "ViewCutPlane needs either of the following sets of keyword"
                " arguments: ('origin', 'normal', 'axis2') or ('origin',"
                " 'pointNormal', 'pointOnAxis2').")

    def initViewerVC(self, svo):
        viewerVC = svo.ViewCut(
            name=self.name, shape=PLANE,
            origin=self.origin, normal=self.normal, axis2=self.axis2)
        self.setCutPositions(viewerVC)
        nameOfCut = viewerVC.name
        return (nameOfCut,)

    def getOpposite(self, name=None, showOnPlane=None):
        """return the viewcut that shows exactly what self does not show.

        @param name: name for the new viewcut. If not specified defaults to
            "anti_%s"%self.name
        @param showOnPlane: If defined can override the middle item of the
            planeAboveOnBelow-triple. I.e. you can force the new viewcut to
            show the cutting plane even though self shows it as well. Or force
            the new viewcut to not show it although self does not show it
            either. If specified must be either TRUE or FALSE. The default
            behaviour for the new viewcut is to show the plane if self does
            not and to not show it if self does.
        """
        switch = {FALSE:TRUE, TRUE:FALSE}
        planeAboveOnBelow = [switch[i] for i in self.planeAboveOnBelow]
        if showOnPlane is not None:
            planeAboveOnBelow[1] = showOnPlane
        if name is None:
            name="anti_%s" % self.name
        return ViewCutPlane(
            name=name,
            planeCutPosition=self.planeCutPosition,
            planeAboveOnBelow=planeAboveOnBelow,
            origin=self.origin, normal=self.normal, axis2=self.axis2,
            followDeformation=self.followDeformation)


class ViewCutIsoSurface(ViewCut):
    """isosurface view cut
    """
    def __init__(self, name,
                 variable, value, showAboveOnBelow):
        self.name = name
        self.variable = variable
        self.value = value
        self.showAboveOnBelow = showAboveOnBelow

    def initViewerVC(self, svo):
        import visualization

        svo.setPrimaryVariable(
            variableLabel=self.variable.varToDisplay,
            outputPosition=self.variable.varPosition,
            refinement=self.variable.varRefinement)

        viewerVC = svo.ViewCut(name=self.name, shape=ISOSURFACE)
        viewerVC.setValues(
            showModelBelowCut=self.showAboveOnBelow[2],
            showModelOnCut=self.showAboveOnBelow[1],
            showModelAboveCut=self.showAboveOnBelow[0],
            overrideAveraging=TRUE,
            value=self.value)

        belowOptions = visualization.OptionArg(
            translucency=ON, translucencyFactor=0.0,
            renderStyle=SHADED, visibleEdges=NONE)
        svo.viewCutOptions.setValues(
            useBelowOptions=TRUE, belowOptions=belowOptions)

        return (viewerVC.name,)


class ViewCutMulti(ViewCut):
    """Viewcut container for the multiple viewcut option.

    @Note: Only one viewcut based on one of the three predefined viewcuts
    (i.e. of class ViewCutCoordPlane with the same axis) can be used in a
    ViewCutMulti instance. I.e. you must not use two
    ViewCutCoordPlane(axis="x", ...)-viewcuts in one ViewCutMulti instance.
    """
    def __init__(self, name, *args):
        """
        @param name: First argument is the name of the view cut
        @param args: Subsequent arguments are ViewCut instances. The
        superposition of everything shown by those instances will be show by
        the multi-viewcut.
        """
        assert isinstance(name, basestring), \
            "First argument of ViewCutMulti must be the name"
        assert all(isinstance(vc, ViewCut) for vc in args), \
            "Second to last argument of ViewCutMulti must be ViewCut instances."
        self.name = name
        self.viewCuts = args

    def initViewerVC(self, svo):
        namesOfViewCuts = list()
        for viewCut in self.viewCuts:
            msg("DEBUG-MultiCut: names for '%s': " % viewCut.name, debugLevel=1)
            currentNames = viewCut.initViewerVC(svo)
            msg("DEBUG-MultiCut: names are %s" % str(currentNames),
                debugLevel=1)
            namesOfViewCuts.extend(currentNames)
        return tuple(namesOfViewCuts)


#} end of View cuts

#------------------------------------------------------------------------------
#{ Display Group Control

class DisplayGroupController(object):
    """Base class for all display group controllers.
    Those classes control which sets / elements to display in a particular
    frame / time step.

    Production classes derived from this must implement a setup and update
    method and update the elemsDisplayed flag.

    @ivar elemsDisplayed: True if there are elements (at least one) displayed,
    False if not. Must be updated by the setup and update methods. This is
    important for layer to see wether they can set on the displayed variable.

    @Note: Do not forget to update self.elemsDisplayed in self.setup() of
    child classes!

    @Note: It may be problematic to use the same DisplayGroupController-instance
    in different layers in one scene. Better create a specific subclass and
    create a different instance for each layer.
    """
    def setup(self, layer, state=None):
        """Basic setup function should be invoked by setup functions of
        derived classes. Optionally sets an initial state
        (displayGroup.replace) and initializes self.elemsDisplayed

        @param state: If not None (if specified) set the initial state. Must
        be ALL_ELEMENTS, ALL_NODES, ALL_SURFACES, DEFAULT_MODEL, EMPTY_LEAF or
        something that dgo.LeafFromElementSets() accepts as input, i.e. an
        elset name or a sequence of elset names.

        Note that "PART-1-1." is being prepended automatically to each given
        elset name. And it's being converted to uppercase.
        """
        self.svo = layer.svo
        self.elemsDisplayed = False

        if state is not None:
            import displayGroupOdbToolset as dgo

            if state in (ALL_NODES, ALL_ELEMENTS, ALL_SURFACES,
                         DEFAULT_MODEL, EMPTY_LEAF):
                leaf = dgo.Leaf(leafType=state)
                self.elemsDisplayed = state in (
                    ALL_ELEMENTS, ALL_SURFACES, DEFAULT_MODEL)
            else:
                if isinstance(state, basestring):
                    self.elemsDisplayed = bool(state)
                    state = "PART-1-1.%s" % state.upper()
                else:
                    # if state is an iterable, anything in it?
                    self.elemsDisplayed = any(state)
                    state = ["PART-1-1.%s" % x.upper() for x in state]
                # display group from elset(s)
                leaf = dgo.LeafFromElementSets(elementSets=state)
            # switch on this display group
            self.svo.displayGroup.replace(leaf=leaf)

    def update(self, currentStepNo, currentFrameNo):
        pass


class DGCtrlStatic(DisplayGroupController):
    """Display group controller class for a certain static display group.
    The objects to be displayed are defined once and remain as they are.

    @param state: If not None (if specified) set the initial state. Must
    be ALL_ELEMENTS, ALL_NODES, ALL_SURFACES, DEFAULT_MODEL, EMPTY_LEAF or
    something that dgo.LeafFromElementSets() accepts as input, i.e. an
    elset name or a sequence of elset names.

    Note that "PART-1-1." is being prepended automatically to each given
    elset name. And it's being converted to uppercase.
    """
    def __init__(self, state):
        self.state = state

    def setup(self, layer):
        DisplayGroupController.setup(self, layer, self.state)

class DGCtrlStaticSomeElems(DisplayGroupController):
    """Display group controller class for a certain static display group.
    The elements to be displayed are defined once and remain as they are.

    @param initialElems: A sequence of element numbers. They are assumed to
       belong to part 'PART-1-1'
    """
    def __init__(self, initialElems):
        self.initialElems = initialElems

    def setup(self, layer):
        import displayGroupOdbToolset as dgo
        DisplayGroupController.setup(self, layer)
        leaf = dgo.LeafFromModelElemLabels(elementLabels=(
            ('PART-1-1', self.initialElems), ))
        layer.svo.displayGroup.replace(leaf=leaf)
        # non empty sequence?
        self.elemsDisplayed = bool(self.initialElems)

class DGCtrlStaticSomeNodes(DGCtrlStatic):
    """Display group controller class for a certain static display group.
    The nodes (no elements) to be displayed are defined once and remain as
    they are.

    @param nodes: A sequence of nodes numbers. They are assumed to belong to
       part 'PART-1-1'
    """
    def __init__(self, nodes):
        self.nodes = nodes

    def setup(self, layer):
        import displayGroupOdbToolset as dgo
        DisplayGroupController.setup(self, layer, EMPTY_LEAF)
        leaf = dgo.LeafFromModelNodeLabels(nodeLabels=(
                ('PART-1-1', self.nodes), ))
        layer.svo.displayGroup.replace(leaf=leaf)
        # non empty sequence?
        self.elemsDisplayed = bool(self.nodes)


class DGCtrlStaticSomeSurfaces(DisplayGroupController):
    """Display group controller class for a certain static display group.
    The surfaces to be displayed are defined once and remain as they are.

    @param surfaces: surface names each given as an individual argument.
       Potential 'PART-1-1.'-prefixes must be included. (Because not all
       surface names have the prefix.)
    """
    def __init__(self, *surfaces):
        self.surfaces = surfaces

    def setup(self, layer):
        import displayGroupOdbToolset as dgo
        DisplayGroupController.setup(self, layer)
        leaf = dgo.LeafFromSurfaceSets(surfaceSets=tuple(self.surfaces))
        layer.svo.displayGroup.replace(leaf=leaf)
        # non empty sequence?
        self.elemsDisplayed = bool(self.surfaces)


class DGCtrlFromSequence(DisplayGroupController):
    """Base class for display group controllers based on a certain sequence.
    Adds or removes elsets according to the current frame / time step.

    @Note: Don't use the same DisplayGroupController-instance in different
    layers in one scene. Better create a specific subclass and create a
    different instance for each layer.
    """
    def __init__(self, seqElsets, initialState):
        """
        @param seqElsets: A string specifying the sequence elsets object in the
          external postData file associated with the odb of the current layer.
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

        @param initialState: ALL_ELEMENTS, ALL_NODES, ALL_SURFACES,
          DEFAULT_MODEL, EMPTY_LEAF or something that dgo.LeafFromElementSets()
          accepts as input, i.e. an elset name or a sequence of elset names,
          e.g.: ('BOLT', 'INSTALLER',). This argument might
          be optional for derived classes like DGCtrlRemoveSeq and
          DGCtrlAddSeq.

        @note: Each elset name from the external postData file associated with
          the odb of the current layer is automatically being prepended by
          "PART-1-1." and is being converted to uppercase.
        """

        # IMPORTANT NOTE: (This might be a bit confusing.)
        # self.seqElsets is only a string identifying the corresponding
        # variable in the external postData file. The actual values are later
        # initialized as self.seqElsetsDict = layer.postData[self.seqElsets].
        # This happens in self.setup().
        self.seqElsets = seqElsets

        self.initialState = initialState

    def initFilter(self, filterElsetsByElement, filterElsetsByElset):
        """Common function of the __init__ method of the ...filter variants
        of ...Add and ...Remove display group controllers. To be called from
        the __init__ methods.

        Initializes self.filterElsetsByElement and self.filterElsetsByElset
        """
        # store filterElsetsByElement, add PART-1-1 - prefix
        if filterElsetsByElement is None:
            self.filterElsetsByElement = None
        elif isinstance(filterElsetsByElement, basestring):
            self.filterElsetsByElement = [
                "PART-1-1.%s" % filterElsetsByElement.upper()]
        else:
            # assume: filterElsetsByElement is an iterable
            self.filterElsetsByElement = [
                "PART-1-1.%s" % x.upper() for x in filterElsetsByElement]

        # convert filterElsetsByElset to a set
        if filterElsetsByElset is None:
            self.filterElsetsByElset = None
        elif isinstance(filterElsetsByElset, basestring):
            self.filterElsetsByElset = set(
                ("PART-1-1.%s" % filterElsetsByElset.upper(), ) )
        else:
            # assume: filterElsetsByElset is an iterable
            self.filterElsetsByElset = set(
                "PART-1-1.%s" % x.upper() for x in filterElsetsByElset)

    def setup(self, layer):
        DisplayGroupController.setup(self, layer, self.initialState)
        try:
            self.seqElsetsDict = layer.postData[self.seqElsets]
        except KeyError:
            m = ("ERROR: seqElsets-variable %s not defined in postData for"
                 " layer %s, odb %s." % (self.seqElsets, layer.name,
                                         layer.odbFileName))
            msg(m)
            raise KeyError(m)

    def setupFilterLeaf(self):
        """Common addition to the setup function for the ...filter variants
        of ...Add and ...Remove display group controllers.

        Creates self.filterLeaf.
        """
        if self.filterElsetsByElement is not None:
            import displayGroupOdbToolset as dgo
            self.filterLeaf=dgo.LeafFromElementSets(
                elementSets=self.filterElsetsByElement)

    def getSetsFromSequence(self, currentStepNo, currentFrameNo):
        """Returns all element set names from the beginning up to the given
        current frame. The sequence is taken from self.seqElsetsDict which is
        the object specified by the seqElsets argument to self.__init__ in the
        postdata file.

        All element sets are assumed to belong to part 'PART-1-1'. This
        function adds the prefix 'PART-1-1.' to each of the elset names.
        """
        setList=list()
        for stepNo in range(1, currentStepNo+1):
            try:
                activeSequence = self.seqElsetsDict[stepNo]
            except KeyError:
                continue

            if stepNo==currentStepNo:
                currentStepEndFrame = min(len(activeSequence),currentFrameNo+1)
            else:
                currentStepEndFrame = len(activeSequence)

            for frNb in range(0, currentStepEndFrame):
                setList.extend(activeSequence[frNb])

        return ['PART-1-1.%s' % x for x in setList]

    def getFilteredLeaf(self, setsToManipulate):
        """Common function of the update method of the ...filter variants
        of ...Add and ...Remove display group controllers.

        Creates the leaf object to be added or subtracted from the displayed
        objects. None if the filter(s) leave nothing.
        """
        leaf = None
        if setsToManipulate and self.filterElsetsByElset is not None:
            setsToManipulate = tuple(self.filterElsetsByElset.intersection(
                setsToManipulate))
        if setsToManipulate:
            import displayGroupOdbToolset as dgo
            if self.filterElsetsByElement is None:
                leaf = dgo.LeafFromElementSets(elementSets=setsToManipulate)
            else:
                # create new display group of all in filterElsetsByElement
                dgToAdd = session.DisplayGroup("ToAdd", self.filterLeaf)
                # strip display group to intersection with current new sets
                dgToAdd.intersect(
                    dgo.LeafFromElementSets(elementSets=setsToManipulate))
                leaf = dgo.LeafFromDisplayGroup(dgToAdd)
        return leaf


class DGCtrlRemoveSeq(DGCtrlFromSequence):
    """Starting from an initial state (usually displaying everything)
    step by step remove elsets as specified in the seqElsets.

    See the description of the base class. The parameter initialState defaults
    to DEFAULT_MODEL (the whole model is displayed initially).

    @Note: self.elemsDisplayed is not updated when elements are removed to the
      last. I.e. When after some frames all elements have been removed step by
      step so that nothing is displayed anymore, self.elemsDisplayed still
      stays True!

    @Note: Don't use the same DisplayGroupController-instance in different
    layers in one scene. Better create a specific subclass and create a
    different instance for each layer.
    """
    def __init__(self, seqElsets, initialState=DEFAULT_MODEL):
        DGCtrlFromSequence.__init__(self, seqElsets, initialState)

    def update(self, currentStepNo, currentFrameNo):
        import displayGroupOdbToolset as dgo
        setsToManipulate = self.getSetsFromSequence(
            currentStepNo, currentFrameNo)
        if setsToManipulate:
            leaf = dgo.LeafFromElementSets(elementSets=setsToManipulate)
            self.svo.displayGroup.remove(leaf=leaf)


class DGCtrlAddSeq(DGCtrlFromSequence):
    """Starting from an initial state (usually displaying nothing)
    step by step add elsets as specified in the seqElsets.

    See the description of the base class. The parameter initialState defaults
    to EMPTY_LEAF (nothing is displayed initially).

    @Note: Don't use the same DisplayGroupController-instance in different
    layers in one scene. Better create a specific subclass and create a
    different instance for each layer.
    """
    def __init__(self, seqElsets, initialState=EMPTY_LEAF):
        DGCtrlFromSequence.__init__(self, seqElsets, initialState)

    def update(self, currentStepNo, currentFrameNo):
        import displayGroupOdbToolset as dgo
        setsToManipulate = self.getSetsFromSequence(
            currentStepNo, currentFrameNo)
        if setsToManipulate:
            leaf = dgo.LeafFromElementSets(elementSets=setsToManipulate)
            self.svo.displayGroup.add(leaf=leaf)
            self.elemsDisplayed = True
            msg("DGCtrlAddSeq <%s> for step %d, frame %d added %d elsets:\n%s"
                % (self.seqElsets, currentStepNo, currentFrameNo,
                   len(setsToManipulate), setsToManipulate),
                debugLevel=5)
        else:
            msg("DGCtrlAddSeq <%s> for step %d, frame %d no change"
                % (self.seqElsets, currentStepNo, currentFrameNo),
                debugLevel=5)


class DGCtrlAddRemoveSeq(DGCtrlFromSequence):
    """Starting from an initial state (usually displaying nothing or whole
    model) step by step add elsets as specified in the seqElsetsAdd and
    remove elsets as specified in seqElsetsRemove.

    See the description of the L{base class<DGCtrlFromSequence>}. The parameter
    initialState defaults to EMPTY_LEAF (nothing is displayed initially).

    In case the same elset is due to be added and removed in the same mining
    step (odb frame) then the set will be added, because that's how it's done
    in our sequences: remove first and add later.

    @Note: Don't use the same DisplayGroupController-instance in different
    layers in one scene. Better create a specific subclass and create a
    different instance for each layer.
    """
    def __init__(self, seqElsetsAdd, seqElsetsRemove, initialState=EMPTY_LEAF):
        DGCtrlFromSequence.__init__(
            self, seqElsets=None, initialState=initialState)
        self.seqElsetsAdd = seqElsetsAdd
        self.seqElsetsRemove = seqElsetsRemove
        self.initialState = initialState

    def setup(self, layer):
        DisplayGroupController.setup(self, layer, self.initialState)
        try:
            self.seqElsetsAddDict = layer.postData[self.seqElsetsAdd]
        except KeyError:
            m = ("ERROR: seqElsets-variable %s not defined in postData for"
                 " layer %s, odb %s." % (self.seqElsetsAdd, layer.name,
                                         layer.odbFileName))
            msg(m)
            raise KeyError(m)
        try:
            self.seqElsetsRemDict = layer.postData[self.seqElsetsRemove]
        except KeyError:
            m = ("ERROR: seqElsets-variable %s not defined in postData for"
                 " layer %s, odb %s." % (self.seqElsetsRemove, layer.name,
                                         layer.odbFileName))
            msg(m)
            raise KeyError(m)

        self.lastStepFrame = (1, 0)

    def getSetsFromSequence(self, currentStepNo, currentFrameNo):
        setListAdd=set()
        setListRem=set()
        for stepNo in range(self.lastStepFrame[0], currentStepNo+1):
            try:
                activeSequenceAdd = self.seqElsetsAddDict[stepNo]
                activeSequenceRem = self.seqElsetsRemDict[stepNo]
            except KeyError:
                continue

            if stepNo==self.lastStepFrame[0]:
                currentStepStartFrame = self.lastStepFrame[1]
            else:
                currentStepStartFrame = 0

            if stepNo==currentStepNo:
                currentStepEndFrame = min(
                    len(activeSequenceAdd), len(activeSequenceRem),
                    currentFrameNo+1)
            else:
                currentStepEndFrame = min(
                    len(activeSequenceAdd), len(activeSequenceRem))

            for frNb in range(currentStepStartFrame, currentStepEndFrame):
                # Note: the order of the following steps matters!
                # Elsets are first removed then added in our sequences.

                # new elsets to be removed: remove those from old add-list
                setListAdd.difference_update(activeSequenceRem[frNb])

                # add elsets to be removed to remove-list
                setListRem.update(activeSequenceRem[frNb])

                # new elsets to be added: remove those from remove-list
                setListRem.difference_update(activeSequenceAdd[frNb])

                # add elsets to be added to add-list
                setListAdd.update(activeSequenceAdd[frNb])

        # update self.lastStepFrame
        self.lastStepFrame = (currentStepNo, currentFrameNo)

        ##### DEBUG
        msg("for step %d frame %d:\nadd: %s\nremove: %s" % (
            currentStepNo, currentFrameNo, sorted(setListAdd),
            sorted(setListRem)))

        return ( ['PART-1-1.%s' % x for x in setListAdd],
                 ['PART-1-1.%s' % x for x in setListRem] )

    def update(self, currentStepNo, currentFrameNo):
        import displayGroupOdbToolset as dgo
        setsToAdd, setsToRemove = self.getSetsFromSequence(
            currentStepNo, currentFrameNo)

        if setsToAdd:
            leaf = dgo.LeafFromElementSets(elementSets=setsToAdd)
            self.svo.displayGroup.add(leaf=leaf)
            self.elemsDisplayed = True
            msg("DGCtrlAddRemoveSeq <%s> for step %d, frame %d added %d elsets:"
                "\n%s" % (self.seqElsetsAdd, currentStepNo, currentFrameNo,
                          len(setsToAdd), setsToAdd), debugLevel=5)
        else:
            msg("DGCtrlAddRemoveSeq <%s> for step %d, frame %d nothing to Add"
                % (self.seqElsetsAdd, currentStepNo, currentFrameNo),
                debugLevel=5)
        if setsToRemove:
            leaf = dgo.LeafFromElementSets(elementSets=setsToRemove)
            self.svo.displayGroup.remove(leaf=leaf)
            self.elemsDisplayed = True
            msg("DGCtrlAddRemoveSeq <%s> for step %d, frame %d removed %d"
                " elsets:\n%s"
                %(self.seqElsetsRemove, currentStepNo, currentFrameNo,
                  len(setsToRemove), setsToRemove), debugLevel=5)
        else:
            msg("DGCtrlAddRemoveSeq <%s> for step %d, frame %d nothing to"
                " Remove"
                % (self.seqElsetsRemove, currentStepNo, currentFrameNo),
                debugLevel=5)


class DGCtrlRemoveSeqFilter(DGCtrlFromSequence):
    """Starting from an initial state (by default displaying everything)
    step by step remove elsets as specified in the seqElsets.

    From the sets in seqElsets remove only the intersection with a certain
    filter region. I.e. nothing will be removed that is not in this filter
    region.

    @Note: Don't use the same DisplayGroupController-instance in different
    layers in one scene. Better create a specific subclass and create a
    different instance for each layer.
    """
    def __init__(self, seqElsets, initialState=DEFAULT_MODEL,
                 filterElsetsByElement=None,
                 filterElsetsByElset=None):
        """
        There are two types of filter that can be used. Typically you would use
        either of them by supplying the corresponding parameter. If both are
        specified they are combined in an exclusive manner. I.e. the resulting
        filter is the intersection of both filters.

        See description of parameter filterElsetsByElement and
        filterElsetsByElset. Both parameters (if given and not None) state a
        tuple (list,...) of elset names as found in the odb. Each item in those
        tuples will be prepended by "PART-1-1.". Those elsets constitute a
        superset of all elements ever to be removed. Or the other way around:
        No element that is not in those elsets will ever be removed. Only the
        intersection of the applicable sequence elsets with this filter elset
        will be removed.

        The "applicable list of sequence elset names" is the list of elsets
        to be removed for a particular mining step as given in the postData
        module of the corresponding layer.

        @param seqElsets: see L{DGCtrlFromSequence}
        @param initialState: see L{DGCtrlFromSequence}
        @param filterElsetsByElement: Elsets listed in this parameter form a
        a filter region independent of the sets in the applicable list of
        sequence elset names. The intersection is computed element wise and
        therefore filterElsetsByElement does not need to contain elsets listed
        in the corresponding seqElsets. The filtering is done element by
        element. Each elset name will be prepended by "PART-1-1." automatically.
        @param filterElsetsByElset: Tuple of elset names as found in the
        applicable list of sequence elset names. From the applicable list of
        sequence elset names only those elsets are considered that are in this
        tuple as well. The elset names in this list must be included in the
        applicable list of sequence elset names as defined in the postData
        module of the corresponding layer. I.e. they must have the same name.
        Each elset name will be prepended by "PART-1-1." automatically.

        @Note: The functionallity provided by filterElsetsByElset is a part of
        what is provided by filterElsetsByElement already. The method used for
        parameter filterElsetsByElset is assumed to be considerably faster.
        """
        DGCtrlFromSequence.__init__(
            self, seqElsets=seqElsets, initialState=initialState)
        self.initFilter(filterElsetsByElement, filterElsetsByElset)

    def setup(self, layer):
        DGCtrlFromSequence.setup(self, layer)
        DGCtrlFromSequence.setupFilterLeaf(self)

    def update(self, currentStepNo, currentFrameNo):
        setsToManipulate = self.getSetsFromSequence(
            currentStepNo, currentFrameNo)
        if setsToManipulate:
            leaf = DGCtrlFromSequence.getFilteredLeaf(self, setsToManipulate)
            self.svo.displayGroup.remove(leaf=leaf)


class DGCtrlAddSeqFilter(DGCtrlFromSequence):
    """Starting from an initial state (usually displaying nothing)
    step by step add elsets as specified in the seqElsets.

    From the sets in seqElsets add only the intersection with a certain filter
    region. I.e. nothing will be added that is not in this filter region.

    @Note: Don't use the same DisplayGroupController-instance in different
    layers in one scene. Better create a specific subclass and create a
    different instance for each layer.
    """
    def __init__(self, seqElsets, initialState=EMPTY_LEAF,
                 filterElsetsByElement=None,
                 filterElsetsByElset=None):
        """
        There are two types of filter that can be used. Typically you would use
        either of them by supplying the corresponding parameter. If both are
        specified they are combined in an exclusive manner. I.e. the resulting
        filter is the intersection of both filters.

        See description of parameter filterElsetsByElement and
        filterElsetsByElset. Both parameters (if given and not None) state a
        tuple (list,...) of elset names as found in the odb. Each item in those
        tuples will be prepended by "PART-1-1.". Those elsets constitute a
        superset of all elements ever to be added. Or the other way around: No
        element that is not in those elsets will ever be added. Only the
        intersection of the applicable sequence elsets with this filter elset
        will be added.

        The "applicable list of sequence elset names" is the list of elsets
        to be added for a particular mining step as given in the postData
        module of the corresponding layer.

        @param seqElsets: see L{DGCtrlFromSequence}
        @param initialState: see L{DGCtrlFromSequence}
        @param filterElsetsByElement: Elsets listed in this parameter form a
        a filter region independent of the sets in the applicable list of
        sequence elset names. The intersection is computed element wise and
        therefore filterElsetsByElement does not need to contain elsets listed
        in the corresponding seqElsets. The filtering is done element by
        element. Each elset name will be prepended by "PART-1-1." automatically.
        @param filterElsetsByElset: Tuple of elset names as found in the
        applicable list of sequence elset names. From the applicable list of
        sequence elset names only those elsets are considered that are in this
        tuple as well. The elset names in this list must be included in the
        applicable list of sequence elset names as defined in the postData
        module of the corresponding layer. I.e. they must have the same name.
        Each elset name will be prepended by "PART-1-1." automatically.

        @Note: The functionallity provided by filterElsetsByElset is a part of
        what is provided by filterElsetsByElement already. The method used for
        parameter filterElsetsByElset is assumed to be considerably faster.
        """
        DGCtrlFromSequence.__init__(
            self, seqElsets=seqElsets, initialState=initialState)
        self.initFilter(filterElsetsByElement, filterElsetsByElset)

    def setup(self, layer):
        DGCtrlFromSequence.setup(self, layer)
        DGCtrlFromSequence.setupFilterLeaf(self)

    def update(self, currentStepNo, currentFrameNo):
        setsToManipulate = self.getSetsFromSequence(
            currentStepNo, currentFrameNo)
        if setsToManipulate:
            leaf = DGCtrlFromSequence.getFilteredLeaf(self, setsToManipulate)
            if leaf is not None:
                self.svo.displayGroup.add(leaf=leaf)
                self.elemsDisplayed = True


#} end of Display Group Control

#------------------------------------------------------------------------------
#{ View (viewing perspective)
#
class View(object):
    """Objects that describe the view. Those objects typically are passed to
    a scene-object. All keyword arguments can (should) be copied as is from
    the abaqus.rpy replay file.

    Usage:
     >>> viewNE = View(
     >>>    name='User-1', nearPlane=1725.3, farPlane=2415, width=19.298,
     >>>    height=14.112, projection=PERSPECTIVE,
     >>>    cameraPosition=(280.65, -2050.7, 48.489),
     >>>    cameraUpVector=(0.0053282, 0.024365, 0.99969),
     >>>    cameraTarget=(1.1104, -0.062256, 0),
     >>>    viewOffsetX=-0.12938, viewOffsetY=2.3477, autoFit=OFF)

    @keyword name: Name of this view in the viewer session. Usually (but not
    restricted to) 'User-1'. If not specified defaults to 'User-1'

    @keyword projection: PERSPECTIVE or PARALLEL

    @ivar name: initialized to the name argument of the constructor
    """
    def __init__(self, **kwargs):
        self.para = dict(**kwargs)
        if "name" not in self.para:
            self.para["name"] = 'User-1'
        self.name = self.para["name"]

    def setupDisplayView(self):
        session.View(**self.para)


class ViewWideAngle(View):
    """View object for "awesome" wide angle views.

    Usage:
     >>> viewPersp = ViewWideAngle(
     >>>    cameraPosition=(280.65, -2050.7, 48.489),
     >>>    cameraDirection=(1,0,0),  # look East
     >>>    )

    @param name: Name of this view in the viewer session. Usually (but not
    restricted to) 'User-1'. If not specified defaults to 'User-1'
    @param cameraPosition: position of the observer.
    @param cameraTarget: Specifies the view direction in terms of a point in
    from of the observer.
    @param cameraDirection: Specifies the view direction in terms of a vector.
    This value is ignored if cameraTarget is given.
    @param cameraUpVector: If not specified defaults to perpendicular to
    cameraDirection and approximately up. Or if cameraDirection is almost
    vertical than cameraUpVector defaults to North (0,1,0).
    @param nearPlane: The distance of the screen (where the 2D-image is being
    drawn) to the observer (cameraPosition). Everything in front of the
    nearPlane (closer to the observer) is not being drawn.
    @param viewAngle: In degrees. Between 0 and 180 degree. The smaller of
    width and height determines whether this is the view angle between left
    and right border or lower and upper border of the image.

    @ivar name: initialized to the name argument of the constructor
    """

    def __init__(self,
                 cameraPosition, cameraTarget=None, cameraDirection=None,
                 cameraUpVector=None,
                 nearPlane=0.5, viewAngle=90.0,
                 name='User-1'):
        from bae.vecmath_01 import vector, vector_plus, cross, norm, length
        from math import pi, tan

        self.name = name

        if cameraTarget is None and cameraDirection is None:
            raise ValueError(
                "cameraTarget or cameraDirection argument must be specified.")
        elif cameraTarget is None:
            cameraTarget = vector_plus(cameraPosition, cameraDirection)
        else:
            # ignore specified cameraDirection, recalculate
            cameraDirection = vector(cameraPosition, cameraTarget)

        if cameraUpVector is None:
            pointLeft = cross([0,0,1], norm(cameraDirection))
            if length(pointLeft) <= 1E-4:
                # looking approx. vertical, then up defaults to North
                cameraUpVector = (0,1,0)
            else:
                cameraUpVector = norm(cross(cameraDirection, pointLeft))

        if viewAngle<=0 or viewAngle>=180:
            raise ValueError(
                "cameraTarget or cameraDirection argument must be specified.")
        width = 2.0*nearPlane*tan(viewAngle*pi/360)

        self.para = dict(
            name=name,
            # Gero thinks: nearPlane is the distance of the screen (where the
            # 2D-image is being drawn) to the observer (cameraPosition)
            nearPlane=2,

            # if (width/height > ratio of image sizes) then only width is
            # relevant, otherwise only height is relvant
            # (only the value with the larger ratio is relevant)
            # Gero thinks: width (if it's relevant) is the size in model units
            # (i.e. metres) of the nearPlane. I.e. a horizontal line in the
            # model of length width at a distance of nearPlane in front of the
            # observer should be represented by a 2D line on the screen exactly
            # reaching from one side of the image on the screen to the other.
            width=width, height=width,

            farPlane=100,    # seems to not have any effect
            projection=PERSPECTIVE,
            cameraPosition=cameraPosition,
            cameraUpVector=cameraUpVector,
            cameraTarget=cameraTarget,
            viewOffsetX=0, viewOffsetY=0,  # moves the image left/right, up/down
            movieMode=ON, autoFit=OFF )
#} end of View


#------------------------------------------------------------------------------
#{ frame lists
#
class FrameList(list):
    """Base class for lists of odb frames.

    This is a python list of (step number, frame number) - tuples. Step numbers
    count from one according to the sequence in the odb.
    """
    pass
    # Removed this in v1.08 because now we'd need the actual odb to retrieve the
    # actual step name. ...
    # def getStepNameFrameNbIterator(self):
    #     for stepNo, frameNo in self:
    #         yield ("Step-%d"%stepNo, frameNo)


class FrListOneStep(FrameList):
    """A list of odb frames.

    This is a python list of (step number, frame number) - tuples.

    @param stepNo: An integer identifying the step number. Step numbers
    count from one according to the sequence in the odb.

    @param frameNoList: list of frame numbers. Negative values can be used
    to count from the end like usual list indices.
    """
    def __init__(self, stepNo, frameNoList):
        FrameList.__init__(
            self, ((stepNo, frameNo) for frameNo in frameNoList))


class FrListSomeSteps(FrameList):
    """A list of odb frames.

    This is a python list of (step number, frame number) - tuples.

    @param stepNoFrameNos: an iterable of (stepNo, frameNoList)-tuples,
    each frameNoList is the list of frameNumbers in the corresponding
    step. Negative frame number values can be used to count from the
    end like usual list indices.
    """
    def __init__(self, stepNoFrameNos):
        """
        @param stepNoFrameNos: an iterable of (stepNo, frameNoList)-tuples,
          each frameNoList is the list of frameNumbers in the corresponding
          step
        """
        FrameList.__init__(self, (
            (stepNo, frameNo)
            for stepNo, frameNoList in stepNoFrameNos
            for frameNo in frameNoList))

#} end of frame lists


##
## ---------------------------------------------------------------------------
##
#{ Variable to display and colourscale, Status variable
##
class VarToDisplay(object):
    r"""Holds the data needed to access a certain variable in the viewer and its
    name for output.

    Examples, integration point:
     >>> class VarLogP(VarToDisplay):
     >>>     "common values: 1.52, 1.71 (5%), 1.92 (8%), 2.0"
     >>>     def __init__(self, conMax=2.0):
     >>>         VarToDisplay.__init__(
     >>>             self, 'SDV4',INTEGRATION_POINT,(),'LogP',
     >>>             conMax, -0.1,10, 'Rainbow', '', '', "%2.2f", UNIFORM)
     >>>         self.varInfo = "LogP10pc"  # variable descriptor in file names
     >>> ..., varToDisplay=VarLogP_SDV4(1.8), ...

     >>> class VarSxx(VarToDisplay):
     >>>     def __init__(self, component, conMax=1E7, conMin=-1E7):
     >>>         "@param component: 'S11', 'S12' or such"
     >>>         VarToDisplay.__init__(
     >>>             self,"S",INTEGRATION_POINT,(COMPONENT,component),component,
     >>>             conMax, conMin, 20,'Rainbow', '', '', "%2.2ePa", UNIFORM)
     >>>         self.varInfo = "Sxx_10MPa"  # variable descriptor in file names

     >>> class VarSMax(VarToDisplay):
     >>>     def __init__(self, conMin=-4.0E7):
     >>>         VarToDisplay.__init__(
     >>>             self, 'S', INTEGRATION_POINT,
     >>>             (INVARIANT,'Min. Principal'),'SMAX',
     >>>             0.0,conMin,20,'Reversed rainbow','','',"%2.2ePa",UNIFORM)
     >>>         # self.varInfo = "SMAX_40MPa"  # description text
     >>>         self.varInfo = ("SMAX_%.0fMPa" % (-conMin/1E6)).replace(".","_")

     >>> class VarPWP(VarToDisplay):
     >>>     def __init__(self, conMax=5.E6):
     >>>         VarToDisplay.__init__(
     >>>             self, 'TEMP', INTEGRATION_POINT, (), 'PWP',
     >>>             conMax, 0.0, 20,'Rainbow', '', '', "%2.2ePa", UNIFORM)

    nodal values
     >>> class VarUR3(VarToDisplay):
     >>>     def __init__(self, conMin=-2.5):
     >>>         VarToDisplay.__init__(
     >>>             self, 'UR',NODAL,(COMPONENT, 'U3'),'SUB',
     >>>             0.0,conMin,10,'Reversed rainbow','','',"%2.2f",UNIFORM)
     >>>         self.varInfo = "SUB_2_5m"  # description text

     >>> class VarURH(VarToDisplay):
     >>>     def __init__(self, conMax=1.0):
     >>>         VarToDisplay.__init__(
     >>>             self, 'URH',NODAL,(),'UH',
     >>>             conMax, 0.0, 10,'Rainbow', '', '', "%2.2f", UNIFORM)
     >>>         self.varInfo = "UH_1m"  # description text

    magnitude of vector
     >>> # incremental displacements
     >>> class VarDUMag(VarToDisplay):
     >>>    def __init__(self, conMax=0.05):
     >>>       VarToDisplay.__init__(
     >>>          self, 'DU',NODAL,(INVARIANT,'Magnitude'),'DUMAG',
     >>>          conMax,0.0,10,'Rainbow','#E10000','#0000DB',"%2.2fm",UNIFORM)

    for cave coupling
     >>> class VarDU3(VarToDisplay):
     >>>     def __init__(self, conMin=-1.0):
     >>>         VarToDisplay.__init__(
     >>>             self, 'DU',NODAL,(COMPONENT,'U3'),'DU3',
     >>>             0.0,conMin,20,'Reversed rainbow',
     >>>             '#0000DB','#E10000',"%2.2fm",UNIFORM)

     >>> class VarCave(VarToDisplay):
     >>>     def __init__(self, abqName='SDV13'):
     >>>         VarToDisplay.__init__(
     >>>             self, abqName, INTEGRATION_POINT, (), 'CAVE',
     >>>             1.0,0.0,10, 'Black to white', '', '', "%2.2f", UNIFORM)

    for USER_DEFINED contours
     >>> class VarUR3USRDEF(VarToDisplay):
     >>>     def __init__(self, conIntervalValues=(-2.,-1.,-.5,-2.,-1.,0.,+.1,+.2,+.5,+1.,+2.) ):
     >>>         conMin = min(conIntervalValues)
     >>>         conMax = max(conIntervalValues)
     >>>         conLev = len(conIntervalValues)-1
     >>>         VarToDisplay.__init__(
     >>>             self, 'UR', NODAL, (COMPONENT, 'U3'), 'UR3',
     >>>             conMax,conMin,conLev, 'Rainbow', '', '', "%2.2f",
     >>>             conIntervalType=USER_DEFINED, conIntervalValues=conIntervalValues)

    @ivar name: Display name of variable for image (e.g. use "LogP"
       if you display SDV4)
    @ivar varToDisplay: Abaqus identifier (SDV1, U, S, E ...)
    @ivar conMax: Maximum value, None means autocompute
    @ivar conMin: Minimum value, None means autocompute
    @ivar formatString: Format string to display min and max values in the
       file name. Any "." will be replaced with "_" for the final image name.
    @ivar varInfo: variable description for file names. Is automatically
       predefined by the constructor and L{setColourScale}() but can be
       overridden manually.
    """
    def __init__(self,
                 varToDisplay, varPosition, varRefinement,
                 name=None,
                 conMax=None, conMin=None,
                 conLev=10, conSpec='Rainbow',
                 outsideLimitsAboveColor='', outsideLimitsBelowColor='',
                 formatString='%s', conIntervalType=UNIFORM,
                 varDisplayName=None,
                 contourStyle=UNIFORM,
                 conIntervalValues=None):
        """Constructor.

        @param varToDisplay: Abaqus identifier (SDV1, U, S, E ...)
        @param varPosition: Abaqus output positiion (INTEGRATION_POINT, NODAL,
           WHOLE_ELEMENT, ...)
        @param varRefinement: Abaqus variable refinement
           (e.g. (INVARIANT,'Magnitude') or (COMPONENT,'S11') or () ...)
        @param name: Display name of variable for image (e.g. use "LogP"
           if you display SDV4)

        @param conMax: Maximum value, None means autocompute
        @param conMin: Minimum value, None means autocompute
        @param conLev: Number of intervals
        @param conSpec: Contour type: 'Rainbow', 'Reversed rainbow',
            'Black to white'
        @param outsideLimitsAboveColor: Colour code for below minimum (or blank)
        @param outsideLimitsBelowColor: Colour code for above minimum (or blank)
        @param formatString: Format string to display min and max values in the
           file name. Any "." will be replaced with "_" for the final image
           name.
        @param conIntervalType: interval type: UNIFORM, LOG, USER_DEFINED
           (this must be understood by the layer used!, also note: for type
           USER_DEFINED, keyword argument conIntervalValues is mandatory)
        @param contourStyle: UNIFORM or CONTINUOUS style of the color legend
        @param varDisplayName: alias for name argument, deprecated, for
           compatibility only
        @param conIntervalValues: the USER_DEFINED contour values, must be a
           tuple of floats in ascending order
        """
        self.varToDisplay  = varToDisplay
        self.varPosition   = varPosition
        self.varRefinement = varRefinement
        self.contourStyle  = contourStyle

        # deprecated alias varDisplayName (if name not given)
        if name is None:
            name = varDisplayName
        self.name           = name
        self.varDisplayName = name  # deprecated alias

        self.setColourScale(
            conMax, conMin, conLev, conSpec,
            outsideLimitsAboveColor, outsideLimitsBelowColor,
            formatString, conIntervalType, contourStyle, conIntervalValues)

    def setColourScale(self, conMax, conMin,
                       conLev=10, conSpec='Rainbow',
                       outsideLimitsAboveColor='', outsideLimitsBelowColor='',
                       formatString='%s', conIntervalType=UNIFORM,
                       contourStyle=UNIFORM,
                       conIntervalValues=None):
        """To (re-)define the colourscale used for this variable.
        Initializes self.varInfo.

        @param conMax: Maximum value, None means autocompute
        @param conMin: Minimum value, None means autocompute
        @param conLev: Number of intervals
        @param conSpec: Contour type: 'Rainbow', 'Reversed rainbow',
           'Black to white'
        @param outsideLimitsAboveColor: Colour code for below minimum (or blank)
        @param outsideLimitsBelowColor: Colour code for above minimum (or blank)
        @param formatString: Format string to display min and max values in the
           file name. Any "." will be replaced by "_" for the final image name.
        @param conIntervalType: interval type: UNIFORM, LOG, USER_DEFINED
           (this must be understood by the layer used!)
        @param contourStyle: A SymbolicConstant specifying the interval style
           of the contour plot. Possible values are CONTINUOUS and UNIFORM.
        @param conIntervalValues: the contour values used for conIntervalType=
           USER_DEFINED, must be a tuple of floats in ascending order
        """

        self.conMax = conMax
        self.conMin = conMin
        self.conLev = conLev
        self.conSpec = conSpec
        self.outsideLimitsAboveColor = outsideLimitsAboveColor
        self.outsideLimitsBelowColor = outsideLimitsBelowColor
        self.conIntervalType = conIntervalType
        self.formatString = formatString
        self.contourStyle = contourStyle
        self.conIntervalValues = conIntervalValues

        if ( (self.conMax==0.0 and self.conMin==0.0)
             or (self.conMax is None or self.conMin is None) ):
            contourMinMaxInfo=""
        else:
            contourMinMaxInfo=(
                ("_"+ self.formatString + "_" + self.formatString)
                % (self.conMax,self.conMin))
        contourMinMaxInfo = contourMinMaxInfo.replace(".","_")
        if self.name:
            self.varInfo = self.name+contourMinMaxInfo
        else:
            self.varInfo = ""

        if self.conIntervalType==USER_DEFINED:
            #if self.conIntervalValues is None:
            if not self.conIntervalValues:
                raise Exception(
                    "VarToDisplay requires an iterable of ascending floats"
                    " for keyword 'conIntervalValues' when using"
                    " intervalType=USER_DEFINED, current value is"
                    " conIntervalValues=%s" % str(conIntervalValues))
            self.varInfo += "_USRDEF"
        # else:
        #     self.varInfo += "_%d" % conLev

    def updateContourColourScale(self, svo):
        """set values in the viewport"""

        # set values
        if self.conMax is not None:
            svo.contourOptions.setValues(
                maxAutoCompute=OFF, maxValue=self.conMax)
        else:
            svo.contourOptions.setValues(maxAutoCompute=ON)
        if self.conMin is not None:
            svo.contourOptions.setValues(
                minAutoCompute=OFF, minValue=self.conMin)
        else:
            svo.contourOptions.setValues(minAutoCompute=ON)

        svo.contourOptions.setValues(
            spectrum=self.conSpec,
            numIntervals=self.conLev, intervalType=self.conIntervalType)

        # set "out of spectrum"-colors according to contour definition
        if self.outsideLimitsAboveColor and self.outsideLimitsBelowColor:
            svo.contourOptions.setValues(
                outsideLimitsMode=SPECIFY,
                outsideLimitsAboveColor=self.outsideLimitsAboveColor,
                outsideLimitsBelowColor=self.outsideLimitsBelowColor)
        else:
            svo.contourOptions.setValues(outsideLimitsMode=SPECTRUM)

        # CONTINUOUS and UNIFORM
        svo.contourOptions.setValues(contourStyle=self.contourStyle)

        # setting contour values for USER_DEFINED contours
        if self.conIntervalType==USER_DEFINED:
            svo.contourOptions.setValues(
                #intervalType=USER_DEFINED,
                intervalValues=self.conIntervalValues,
                contourType=BANDED,
                )

# strain
class VarPST_SDV1(VarToDisplay):
    """use like this:
     >>> ... , varToDisplay=VarPST_SDV1(0.1), ...

    common values: 0.05 (5%), 0.08 (8%), 0.1
    """
    def __init__(self, conMax=0.1, conMin=1E-4):
        VarToDisplay.__init__(
            self, 'SDV1',INTEGRATION_POINT,(),'PST',
            conMax, conMin,10, 'Rainbow', '', '', "%2.2f", UNIFORM)
        # description text for file names
        self.varInfo = ("PST_%.1fpc" % (conMax*100)).replace(".","_")

class VarLogP_SDV1(VarToDisplay):
    """use like this:
     >>> ... , varToDisplay=VarLogP_SDV1(0.1), ...

    common values: 0.05 (5%), 0.08 (8%), 0.1
    """
    def __init__(self, conMax=0.1, conMin=1E-4):
        VarToDisplay.__init__(
            self, 'SDV1',INTEGRATION_POINT,(),'LogP',
            conMax, conMin,10, 'Rainbow', '', '', "%2.2f", LOG)
        # description text for file names
        self.varInfo = ("LogP_%.1fpc" % (conMax*100)).replace(".","_")

class VarLogP(VarToDisplay):
    """use like this:
     >>> ..., varToDisplay=VarLogP(1.8), ...

    common values: 1.52, 1.71 (5% PST), 1.92 (8% PST), 2.0
    """
    def __init__(self, conMax=2.0):
        VarToDisplay.__init__(
            self, 'SDV4',INTEGRATION_POINT,(),'LogP',
            conMax, -0.1,10, 'Rainbow', '', '', "%2.2f", UNIFORM)

class VarPEEQ(VarToDisplay):
    def __init__(self, conMax=0.001):
        VarToDisplay.__init__(
            self, 'PEEQ',INTEGRATION_POINT,(),'PEEQ',
            conMax, 0.0, 20, 'Rainbow', '', '', "%2.2f", UNIFORM)

# energy
class VarRER_SDV12(VarToDisplay):
    """use like this:
     >>> ... , varToDisplay=VarRER_SDV12(1E6), ...

    common values: 1E6
    """
    def __init__(self, conMax=1E6, conMin=1):
        VarToDisplay.__init__(
            self, 'SDV12',INTEGRATION_POINT,(),'RER',
            conMax, conMin,12, 'Rainbow', '', '', "%2.2e", LOG)

# stress
class VarSMax(VarToDisplay):
    def __init__(self, conMin=-4.0E7):
        VarToDisplay.__init__(
            self, 'S',INTEGRATION_POINT,(INVARIANT,'Min. Principal'),'SMAX',
            0.0, conMin, 20,'Reversed rainbow', '', '', "%2.2ePa", UNIFORM)
        # description text for file names
        self.varInfo = ("SMAX_%.0fMPa" % (-conMin/1E6)).replace(".","_")

class VarSMin(VarToDisplay):
    def __init__(self, conMax=0.0, conMin=-5E6):
        VarToDisplay.__init__(
            self, 'S',INTEGRATION_POINT,(INVARIANT,'Max. Principal'),'SMIN',
            conMax, conMin, 20,'Reversed rainbow', '', '', "%2.2ePa", UNIFORM)
        if conMax:
            self.varInfo = (
                "SMIN_%.2ePa_%.2ePa" % (-conMax, -conMin)).replace(".","_")
        else:
            self.varInfo = ("SMIN_%.2ePa" % (-conMin)).replace(".","_")

# displacements
class VarUMag(VarToDisplay):
    def __init__(self, conMax=1.0):
        VarToDisplay.__init__(
            self, 'U',NODAL,(INVARIANT,'Magnitude'),'UMAG',
            conMax,0.0,10,'Rainbow','#E10000','#0000DB',"%03.1fm",UNIFORM)

class VarURMag(VarToDisplay):
    def __init__(self, conMax=0.5):
        VarToDisplay.__init__(
            self, 'UR',NODAL,(INVARIANT,'Magnitude'),'UMAG',
            conMax,0.0,10,'Rainbow','#E10000','#0000DB',"%2.2fm",UNIFORM)

# velocities
class VarVMag(VarToDisplay):
    def __init__(self, conMax=1.0):
        VarToDisplay.__init__(
            self, 'V',NODAL,(INVARIANT,'Magnitude'),'VMAG',
            conMax, 0.0, 10,'Rainbow', '', '', "%1.1f", UNIFORM)

# state
varStatus = VarToDisplay(
    'STATUS',WHOLE_ELEMENT,(),'STATUS',
    1.0, 0.0, 5, 'Rainbow', '', '', "%2.2f", UNIFORM)

#} end of Variable to display and colourscale, Status variable

##
## ---------------------------------------------------------------------------
##
#{ query and store parameters of interactive session

#: File name where to store parameter from the current session
recordParaFileName = "recordedSessionPara.py"

def recordParaAddSeparator():
    """Add a separator to the recordPara file.
    """
    with open(recordParaFileName, "a") as f:
        f.write("\n%s\n\n" % ("#"*79))
    return

def recordParaView():
    """Store the current view and imageSize in the recordPara file."""
    sv = session.viewports[session.currentViewportName]

    # default size to scale: (1280, 1024)
    factor = 1280.0 / sv.width
    imageSize = (1280.0, sv.height*factor)
    imageSize = tuple(int(round(x)) for x in imageSize)
    with open(recordParaFileName, "a") as f:
        f.write("imageSize = %s\n" % str(imageSize))

    # current view para
    lines = ",\n    ".join(

        # each line consists of key=val tuples separated by ", "
        ", ".join(
            "%s=%s" % (attr, getattr(sv.view, attr))
            for attr in attrLine)

        for attrLine in [
                ["nearPlane", "farPlane"],
                ["width", "height", "projection"],
                ["cameraPosition"],
                ["cameraUpVector"],
                ["cameraTarget"],
                ["viewOffsetX", "viewOffsetY"],
                ])

    with open(recordParaFileName, "a") as f:
        f.write('viewXXX = View(\n    name="User-1", autoFit=OFF, %s)\n'
                % lines)

    return

@hidden
def getVarInfo(svoVar):
    varPosition = {
        2: INTEGRATION_POINT,
        1: NODAL,
        5: WHOLE_ELEMENT}[svoVar[1]]
    if svoVar[4]==0:
        refinement = ()
    elif svoVar[4]==1:
        refinement = (INVARIANT, svoVar[5])
    elif svoVar[4]==2:
        refinement = (COMPONENT, svoVar[5])
    else:
        raise ValueError(
            "Refinement %d not implemented." % svoVar[4])
    return (svoVar[0], varPosition, refinement)

def recordParaLayer():
    """Store the current viewport configuration as Layer-object in the
    recordPara file."""

    sv = session.viewports[session.currentViewportName]
    svo = sv.odbDisplay

    # determine viewcut
    viewCut = None
    if svo.viewCut==ON:
        viewCuts = list()
        for vcName in svo.viewCutNames:
            viewerVC = svo.viewCuts[vcName]
            planeAboveOnBelow = (
                viewerVC.showModelAboveCut,
                viewerVC.showModelOnCut,
                viewerVC.showModelBelowCut)

            if vcName in ('X-Plane', 'Y-Plane', 'Z-Plane'):
                vcStr = 'ViewCutCoordPlane("%s", "%s", %s, %s' % (
                    vcName, vcName[0], viewerVC.position, planeAboveOnBelow)
                if viewerVC.followDeformation:
                    vcStr += ", followDeformation=True"
                vcStr += ")"

            elif viewerVC.shape==PLANE:
                vcStr = (
                    'ViewCutPlane("%s", %s, %s,'
                    ' origin=%s, normal=%s, axis2=%s' % (
                        vcName, viewerVC.position, planeAboveOnBelow,
                        viewerVC.origin, viewerVC.normal, viewerVC.axis2))
                if viewerVC.followDeformation:
                    vcStr += ", followDeformation=True"
                vcStr += ")"

            elif viewerVC.shape==ISOSURFACE:
                vcStr = (
                    'ViewCutIsoSurface("%s", %s, %s, %s)'
                    % (vcName, "varIsoSurf",
                       viewerVC.position, planeAboveOnBelow))
            viewCuts.append(vcStr)

        if len(viewCuts) == 1:
            viewCut = viewCuts[0]
        else:
            viewCut = ('ViewCutMulti("viewCutXXX", [\n%s]\n' %
                       ",\n        ".join(viewCuts))

    # deformed display?
    deformed = svo.display.plotState in (CONTOURS_ON_DEF, DEFORMED)
    if deformed:
        deformedStr = ""
    else:
        deformedStr = "    deformed=False,\n"

    # determine other plot options
    if svo.commonOptions.translucency:
        translucencyStr = ("    translucencyFactor=%4.2f,\n"
                           % svo.commonOptions.translucencyFactor)
    else:
        translucencyStr = ""
    visibleEdges = svo.commonOptions.visibleEdges

    # determine status variable
    statusVarStr = ""
    if svo.useStatus and (deformed or svo.applyStatusToUndeformed):
        statusVarInfo = getVarInfo(svo.statusVariable)
        statusVarStr = (
            '    statusVariable=VarToDisplay("%s", %s, %s),\n'
            '    statusMinimum=%s,\n'
            '    statusMaximum=%s,\n'
            '    statusInsideRange=%s,\n'
            % (statusVarInfo[0], statusVarInfo[1], statusVarInfo[2],
               svo.statusMinimum, svo.statusMaximum, svo.statusInsideRange))

    # LayerContour
    if svo.display.plotState in [(CONTOURS_ON_UNDEF,), (CONTOURS_ON_DEF,)]:

        # deformed = svo.display.plotState==CONTOURS_ON_DEF
        # ... add argument "deformed=deformed," ...

        # determine varToDisplay
        varInfo = getVarInfo(svo.primaryVariable)
        conMin = None
        if svo.contourOptions.minAutoCompute==OFF:
            conMin = svo.contourOptions.minValue
        conMax = None
        if svo.contourOptions.maxAutoCompute==OFF:
            conMax = svo.contourOptions.maxValue
        if varInfo[0] == "SDV1":
            varToDisplay = "VarPST_SDV1(conMax=%s)" % conMax
        elif varInfo[0] == "S" and varInfo[2] == (INVARIANT, "Min. Principal"):
            varToDisplay = "VarSMax(conMin=%s)" % conMin
        elif varInfo[0] == "S" and varInfo[2] == (INVARIANT, "Max. Principal"):
            varToDisplay = "VarSMin(conMax=%s)" % conMax
        elif varInfo[0] == "U":
            varToDisplay = "VarUMag(conMax=%s)" % conMax
        elif varInfo[0] == "UR":
            varToDisplay = "VarURMag(conMax=%s)" % conMax
        else:
            varToDisplay = (
                'VarToDisplay(\n    %r, %s, %s,\n    name="varXXX",\n'
                % varInfo)
            newLine = False
            if conMax:
                varToDisplay += '    conMax=%s,' % conMax
                newLine = True
            if conMin:
                varToDisplay += '    conMin=%s,' % conMin
                newLine = True
            if newLine:
                varToDisplay += '\n'

            varToDisplay += '    conLev=%s,\n' % svo.contourOptions.numIntervals
            varToDisplay += '    conSpec=%r,\n' % svo.contourOptions.spectrum

            if svo.contourOptions.outsideLimitsMode==SPECIFY:
                varToDisplay += ('    outsideLimitsAboveColor=%r,\n'
                                 % svo.contourOptions.outsideLimitsAboveColor)
                varToDisplay += ('    outsideLimitsBelowColor=%r,\n'
                                 % svo.contourOptions.outsideLimitsBelowColor)
            if svo.contourOptions.intervalType != UNIFORM:
                varToDisplay += ('    conIntervalType=%s,\n'
                                 % svo.contourOptions.intervalType)
            if svo.contourOptions.contourStyle != UNIFORM:
                varToDisplay += ('    contourStyle=%s,\n'
                                 % svo.contourOptions.contourStyle)
            varToDisplay += '    )'

        with open(recordParaFileName, "a") as f:
            f.write('layerXXX = LayerContour(\n    "XXX",\n')
            f.write('    varToDisplay=%s,\n'
                    % "\n    ".join(varToDisplay.splitlines()))
            f.write(statusVarStr)
            if viewCut:
                f.write('    viewCut=%s,\n' % viewCut)
            f.write('    displayGroupController=None,\n')
            f.write(translucencyStr)
            f.write(deformedStr)
            if visibleEdges!="FEATURE":
                f.write('    visibleEdges=%s,\n' % visibleEdges)
            f.write('    )\n')

    # LayerSolidColour
    elif (svo.display.plotState in [(UNDEFORMED,), (DEFORMED,)]
            and sv.colorMode==DEFAULT_COLORS):

        with open(recordParaFileName, "a") as f:
            f.write('layerXXX = LayerSolidColour(\n    "XXX",\n')
            f.write('    fillColor="%s",\n' % svo.commonOptions.fillColor)
            f.write(statusVarStr)
            if viewCut:
                f.write('    viewCut=%s,\n' % viewCut)
            f.write('    displayGroupController=None,\n')
            f.write(translucencyStr)
            f.write(deformedStr)
            if visibleEdges!="NONE":
                f.write('    visibleEdges=%s,\n' % visibleEdges)
            f.write('    )\n')

        ### to be added... LayerSolidColour(...)
        pass

    # LayerMultiColour
    elif svo.display.plotState in [(UNDEFORMED,), (DEFORMED,)]:
        colorMapType = {
            "SET_MAP_COLORS": 'Element set',
            "SECTION_MAP_COLORS": 'Section',
            "MATERIAL_MAP_COLORS": 'Material',
            "ELTYPE_MAP_COLORS": 'Element type',
            "AVERAGING_REGION_MAP_COLORS": 'Averaging region',
            "INTERNAL_SET_MAP_COLORS": 'Internal set',
            }[str(sv.colorMode)]
        colorMap = dict(
            (it[0], it[1] and it[2] or None)
            for it in sv.colorMappings[colorMapType].attributeColors)

        with open(recordParaFileName, "a") as f:
            f.write('layerXXX = LayerMultiColour(\n    "XXX",\n')
            f.write('    colorMapType="%s",\n' % colorMapType)
            f.write('    colorMap=%s,\n' % colorMap)
            f.write(statusVarStr)
            if viewCut:
                f.write('    viewCut=%s,\n' % viewCut)
            f.write('    displayGroupController=None,\n')
            f.write(translucencyStr)
            f.write(deformedStr)
            if visibleEdges!="NONE":
                f.write('    visibleEdges=%s,\n' % visibleEdges)
            f.write('    )\n')

    return

##
#} end of query and store parameters of interactive session
