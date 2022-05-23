from bae.log_01 import msg
import os
import collections

#import numpy (required for matplotlib)
try:
    import numpy as np
except ImportError:
    raise NotImplementedError('This module requires numpy')

#import matplotlib
try:
    import matplotlib as mpl
except ImportError:
    raise NotImplementedError('This module requires matplotlib')

#import pyplot an set backend
useBackend = None  # "TkAgg"
try:
    import matplotlib.pyplot as plt  # this will fail on nukus (no x-server)
    mplBackend = mpl.get_backend()
    if useBackend and not mplBackend==useBackend:
        mpl.use(useBackend)
        msg("matplotlib backend is changed from %s (default) to %s"%
            (mplBackend,useBackend))
except ImportError:                  # force backend 'Agg'
    mpl.use('Agg')
    import matplotlib.pyplot as plt



### checking matplotlib version
mplVersion = mpl.__version__
msg("matplotlib version: %s" % mplVersion)

### setting desired matplotlib default-parameters (based on mpl version >= 2.0)
mpl.rcdefaults()
msg("matplotlib style settings: reset the defaults")

logoPath = os.path.join(os.path.dirname(__file__),"logoBE_300_white.png")

useStyleDir = os.path.dirname(__file__).replace(os.sep,"/")+"/"
# useStyle = "ggplot"
# useStyle = "dark_background"
useStyle = "mpl_styleBE.mplstyle"

if mplVersion.startswith("2."):
    msg("style to use: '%s'" % (useStyleDir+useStyle))
    plt.style.use(useStyleDir+useStyle)
    msg("matplotlib style settings: used from '%s'" % useStyle+"\n")
else:
    msg("WARNING: Default style %s required matplotlib version >= 2.0."
        " Current version is %s. Will not use the plot style."
        % (useStyle, mplVersion))
    useStyle = None

def truncateColormap(cmap, minval=0.0, maxval=1.0, n=100):
    #To create a truncated subset of a given colormap
    import matplotlib.colors as colors
    newCmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return newCmap

### ##########################################################################
### (MATPLOTLIB) M A I N P L O T
### ##########################################################################
class Mainplot(object):
    """base class for the main plot.

    This is the standard case, and can be directly be used for most cases.

    For special purposes (like for a standardized report plots containing a
    standardized arrangement of subplots), a particuar derived class can be
    outsourced to mpl_samples.py and be called from the driver script.

    Usage
    =====


    getOutFilename function
    =======================

    To determine output file name(s) for the plot and plotprogress methods
    you can either supply corresponding arguments or define an instance
    method getOutFilename(idx). This function is expected to return an
    appropriate filename to save the plot for the current index of the data.
    If None is returned, the plot won't be saved when created.
     >>> myMainplot = MainPlot()
     >>> sp = SubPlot(...)
     >>> sp.addData(...)
     >>> myMainplot.addSubplot(sp)
     >>> ### now before plotting and in order to save as PNGs:
     >>> def getFilename(idx=None):
     >>>     pngFiletmpl = "PrjXXXX_Rxx_Qxx_WYSIWYG_%s.png"
     >>>     if idx is None:
     >>>         frmId = "ALL"
     >>>     else:
     >>>         frmId = lkpFrameIdFromIdx[idx]
     >>>     return pngFiletmpl%frmId
     >>>
     >>> myMainplot.getOutFilename = getFilename
     >>> myMainplot.plot()
     >>> myMainplot.plotprogress(range(0,30,3))

    @ivar fig: the matplotlib figure as self.fig with figsize
        (width = nCols*gridScale, height = nRows*gridScale) and gridScale as
        gridSpacing dimension in cm.
    @ivar gridSpacingCm: the gridspacing in cm.
    @ivar gridSpacingIn: the gridspacing in inches.
    @ivar dpi: the dpi, defaults to 150., or when using None and style
        'mpl_styleBE.mplstyle' also defaults to 150. The method
        self.fig.get_dpi() should confirm this self.dpi value.
    @ivar subplots: a list of L{bae.plot_01.mpl_core.SubPlot} instances,
        which can be added using L{self.addSubplot} method
    @ivar subplotNames: a list of subplot.names (used to iterate over
        subplots in a consistent sort order)
    """

    def __init__(self, figNb=1, gridRowsCols=(1,1), gridScale=None, dpi=150.):
        """creates an instance of the main plot, containing

        @param figNb: optional, figure number of matplotlib.Figure as integer
        @param gridRowsCols: optional, tuple of (nRows,nCols) of the grid
        @param gridScale: optional, the grid spacing (in cm). defaults to None,
          then the gridScale is computed such that the figwidth is 17.34cm
          (this works well with dpi=150. and the mpl_styleBE.mplstyle settings).
        @param dpi: optional, the dots-per-inch, dpi=150. recommended, which is
          equivalent to dpi=None when using the mpl_styleBE.mplstyle settings.

        @Note: with figure width of 17.34cm (= nCols*gridScale) and dpi=150. the
        figure width is 1024px.
        """

        ### ----- OVERALL FIGURE AND LAYOUT -----
        self.figNb = figNb
        self.gridRowsCols = gridRowsCols
        if gridScale is None:
            gridScale = 17.34/gridRowsCols[1]
        self.gridScaleCm = gridScale          # cm
        self.gridScaleIn = gridScale/2.54     # cm to inch
        self.dpi = dpi
        self.fig = plt.figure(figNb,                            
            figsize=list(reversed([ self.gridScaleIn*v for v in gridRowsCols ])),
            dpi=dpi,
            facecolor='white',
            edgecolor=None,
            linewidth=0.0,
            frameon=None,
            subplotpars=None,
            tight_layout=True,
            )
        self.fig.clear()
        #fig = self.fig
        #msg("DEBUG: created a MainPlot figure with DPI=%s, sizeInches=%s"
        #    %(str(fig.get_dpi()),str(fig.get_size_inches())))

        ### ----- list of SUBPLOTS -----
        self.subplots = dict()
        self.subplotNames = list()

        return


    def addSubplot(self, subplot):
        """add a L{bae.plot_01.Subplot} instance to L{self.subplots}
        and its name to L{self.subplotNames}.
        """
        name = subplot.name
        if name not in self.subplotNames:
            self.subplots[name] = subplot
            self.subplotNames.append(name)
        else:
            raise Exception("A SubPlot instance with name '%s' already "
                            "exists on MainPlot(%d)."%(name, self.figNb))

        return


    def show(self, block=True):
        """show MainPlot.fig (matplotlib.Figure) in interactive mode.
        """
        plt.show(block=block)


    def plot(self, outFilename=None):
        """updates the MainPlot.fig (matplotlib.Figure) showing all
        subplots and their data, if any.

        @param outFilename: The file name for the resulting png file.

           Note: An optionally supplied self.getOutFilename function
           takes precedence. See class description.
        """
        for spName in self.subplotNames:
            sp = self.subplots[spName]
            sp.setup()
            sp.plot()

        if hasattr(self, "getOutFilename"):
            outFilename = self.getOutFilename(None)
            msg("taking self.getOutFilename(None) for outFilename -> '%s'"%outFilename)

        if outFilename is not None:
            self.fig.savefig(
                outFilename,
                transparent=False
                )
        return

    def plotprogress(self, idx=None, show_marker=True, show_rest=True,
                     verbose=True, outFileNameTemplate=None):
        """updates the MainPlot.fig (matplotlib.Figure) showing all
        subplots and their data, if any, upto a certain idx.

        @param idx: the range of indeces to plot (and save) the data.
            idx can either be an integer or a list of integers. if saved, this
            produces as many pictures as the length of the list is.
        @param show_marker: boolean, whether or not to show a marker at the
            position of the current index.
        @param show_rest: boolean, whether or not to show the rest of the
            curve(s) from the position of the current index up to the end of
            the datavalues.
        @param verbose: boolean, for diagnostic output using L{bae.log_01.msg}.
        @param outFileNameTemplate: How shall the file name for the resulting
           png file be generated? Can either be a template string with one
           integer placeholder (%d). Or a function to be called with one integer
           file index as argument to deliver the complete file name.

           Note: An optionally supplied self.getOutFilename function
           takes precedence. See class description.
        """
        if isinstance(idx, int):
            idx = [ idx, ]
        elif idx is None:
            N = max(len(curveX)
                    for sp in self.subplots.itervalues()
                    for curveX in sp.xValues)
            idx = range(N)

        if hasattr(self, "getOutFilename"):
            getOutFilename = self.getOutFilename
        elif outFileNameTemplate is None:
            getOutFilename = lambda i: None
        elif isinstance(outFileNameTemplate, basestring):
            getOutFilename = lambda i: (outFileNameTemplate % i)
        elif callable(outFileNameTemplate):
            getOutFilename = outFileNameTemplate
        else:
            raise ValueError(
                "ERROR: outFileNameTemplate argument of type %s not"
                " implemented." % type(outFileNameTemplate))

        for i in idx:
            if verbose: msg("  - plotprogress idx %d..."%i)
            for spName in self.subplotNames:
                sp = self.subplots[spName]
                sp.setup()
                sp.plotsingle(i, show_marker=show_marker, show_rest=show_rest)
            if verbose: msg("    plotted")

            outFilename = getOutFilename(i)
            if outFilename is not None:
                self.fig.savefig(
                    outFilename,
                    transparent=False
                    )
                if verbose: msg("    saved to '%s'"%outFilename)
        return


### ##########################################################################
### (MATPLOTLIB) S U B P L O T
### ##########################################################################
class Subplot(object):
    """base class for Subplots that are going to be added to Mainplot(s).

    This is the standard 2D_rectilinear case, and can be directly be used for
    these cases.

    For other cases (like polar plots, 3d plots, etc), or for special purposes
    of this case (like S1_S3 diagrams with or without yield curve, etc), new
    classes can be derived from this base class.
    """

    def __init__(self, name, mpInst, mpGridPos, mpGridSpan,
            bgcolor='white',
            projection=u'rectilinear'):
        """
        creating a subplot being positioned on a L{Mainplot} instance
        using 2D rectilinear coordinate system (typical x-y plots).
        """
        self.name = name
        self.mpInst = mpInst
        self.mpGridPos = mpGridPos
        self.mpGridSpan = mpGridSpan

        try:
            # for matplotlib 2.x
            self.ax = plt.subplot2grid(
                shape=mpInst.gridRowsCols,
                loc=mpGridPos,
                rowspan=mpGridSpan[0],
                colspan=mpGridSpan[1],
                #fig=mpInst.figNb,
                #fig=mpInst.fig,
                facecolor=bgcolor,
                projection=projection,
                )
        except AttributeError:
            # for matplotlib 1.5
            self.ax = plt.subplot2grid(
                shape=mpInst.gridRowsCols,
                loc=mpGridPos,
                rowspan=mpGridSpan[0],
                colspan=mpGridSpan[1],
                #fig=mpInst.figNb,
                #fig=mpInst.fig,
                axisbg=bgcolor,
                projection=projection,
                )

        self.xValues = []
        self.yValues = []
        self.axLimitsUser = None
        self.axLimitsAuto = None
        self.computeAxLimitsAuto()
        self.crvProps = []
        self.axTitle  = dict([ ('text',""), ('props',dict()), ])
        self.axXLabel = dict([ ('text',""), ('props',dict()), ])
        self.axYLabel = dict([ ('text',""), ('props',dict()), ])
        return


    def addData(self, ptDataDict, xId, yId, lc=None, lp=None):
        ### y-data
        try:
            yvals = ptDataDict[yId]
        except:
            crvNb = str(len(self.yValues))
            if crvNb.endswith("1"): crvNb+="st"
            if crvNb.endswith("2"): crvNb+="nd"
            if crvNb.endswith("3"): crvNb+="rd"
            else: crvNb+="th"
            raise Exception("could not set y-data for the %s curve (yId='%s')"
                            % (crvNb,yId))
        else:
            self.yValues.append(yvals)
            crvProps = dict()
            if lc is not None:
                crvProps['color'] = lc
            else:
                crvNb = len(self.yValues)-1
                crvProps['color'] = "C%d"%crvNb  # using the mpl.cycler
            if lp is not None:
                if 'c' in lp.keys():
                    lp['color'] = lp['c']
                    del lp['c']
                crvProps.update(lp)
            self.crvProps.append(crvProps)
        ### x-data
        try:
            xvals = ptDataDict[xId]
        except KeyError:
            xvals = map(float,range(len(yvals)))
        self.xValues.append(xvals)
        self.computeAxLimitsAuto()
        return


    def addDataXY(self, xvals, yvals, lc=None, lp=None):
        """Store data for a single line.

        @param xvals: vector of x-values
        @param yvals: vector of y-values
        @param lc: line colour ("blue", "red", "black",...?)
        @param lp: other line props: a dict {prop: val}

        @Note: self.addData(ptDataDict, xId, yId, lc, lp) is equivalent to
        self.addDataXY(ptDataDict[xId], ptDataDict[yId], lc, lp)
        """
        crvNb = len(self.yValues)

        ### y-data
        self.yValues.append(yvals)

        ### line colour, line props
        crvProps = dict()
        if lc is not None:
            crvProps['color'] = lc
        else:
            crvProps['color'] = "C%d"%crvNb  # using the mpl.cycler
        if lp is not None:
            if 'c' in lp:
                lp['color'] = lp['c']
                del lp['c']
            crvProps.update(lp)
        self.crvProps.append(crvProps)

        ### x-data
        self.xValues.append(xvals)
        self.computeAxLimitsAuto()
        return


    def computeAxLimitsAuto(self):
        if self.axLimitsAuto is None:
            self.axLimitsAuto = [ [ None, None, None ], [ None, None, None ] ]

        ### x-axis
        if self.xValues:
            xminList = [ min(xvals) for xvals in self.xValues if len(xvals) ]
            if xminList:
                xmin = min(xminList)
            else:
                xmin = 0.
            xmaxList = [ max(xvals) for xvals in self.xValues if len(xvals) ]
            if xmaxList:
                xmax = max(xmaxList)
            else:
                xmax = 1.
        else:
            xmin, xmax = self.ax.get_xlim() # only makes sense after any plot, otherwise [0.,1.]
        self.axLimitsAuto[0] = [xmin,xmax,None]

        ### y-axis
        if self.yValues:
            yminList = [ min(yvals) for yvals in self.yValues if len(yvals) ]
            if yminList:
                ymin = min(yminList)
            else:
                ymin = 0.
            ymaxList = [ max(yvals) for yvals in self.yValues if len(yvals) ]
            if ymaxList:
                ymax = max(ymaxList)
            else:
                ymax = 1.
        else:
            ymin, ymax = self.ax.get_ylim() # only makes sense after any plot, otherwise [0.,1.]
        self.axLimitsAuto[1] = [ymin,ymax,None]

        return

    def setAxLimitsUser(self, x=None, y=None):
        if self.axLimitsUser is None:
            self.axLimitsUser = [ [ None, None, None ], [ None, None, None ] ]

        ### x-axis
        if x is None:
            axx = [ None, None, None ]
        else:
            if len(x)==2:
                axx = [ x[0], x[1], None ]
            elif len(x)==3:
                axx = [ x[0], x[1], x[2] ]
        self.axLimitsUser[0] = axx

        ### y-axis
        if y is None:
            axy = [ None, None, None ]
        else:
            if len(y)==2:
                axy = [ y[0], y[1], None ]
            elif len(y)==3:
                axy = [ y[0], y[1], y[2] ]
        self.axLimitsUser[1] = axy

        return

    def applyAxLimits(self):
        if self.axLimitsUser is None:
            self.axLimitsUser = [ [ None, None, None ], [ None, None, None ] ]
        #msg("DEBUG: self.axLimitsAuto = %s"%str(self.axLimitsAuto))
        #msg("DEBUG: self.axLimitsUser = %s"%str(self.axLimitsUser))

        ### x-axis
        xmin,xmax,xticks = self.axLimitsUser[0]
        if xmin is None:
            xmin = self.axLimitsAuto[0][0]
        if xmax is None:
            xmax = self.axLimitsAuto[0][1]
        if xticks is None:
            xticks = self.axLimitsAuto[0][2]
            if xticks is not None:
                xticks = [ x for x in xticks if xmin<=x<=xmax ]

        ### y-axis:
        ymin,ymax,yticks = self.axLimitsUser[1]
        if ymin is None:
            ymin = self.axLimitsAuto[1][0]
        if ymax is None:
            ymax = self.axLimitsAuto[1][1]
        if yticks is None:
            yticks = self.axLimitsAuto[1][2]
            if yticks is not None:
                yticks = [ y for y in yticks if ymin<=y<=ymax ]

        ### apply
        #msg("DEBUG: appling %s"%str([[xmin,xmax,xticks],[ymin,ymax,yticks]]))

        self.ax.set_xlim(xmin,xmax)
        self.ax.set_ylim(ymin,ymax)
        if xticks is not None: self.ax.set_xticks(xticks)
        if yticks is not None: self.ax.set_yticks(yticks)

        return

    def set_title(self, titletext, **kwargs):
        self.axTitle['text'] = titletext
        self.axTitle['props'].update(**kwargs)
        return

    def setup(self):
        #tmp_title = self.ax.get_title()
        tmp_xlabel = self.ax.get_xlabel()
        tmp_ylabel = self.ax.get_ylabel()
        #msg("DEBUG: tmp_title = %s"%tmp_title)

        self.ax.clear()
        #self.ax.set_title(tmp_title)
        self.ax.set_title(self.axTitle['text'], self.axTitle['props'])
        self.ax.set_xlabel(tmp_xlabel)
        self.ax.set_ylabel(tmp_ylabel)

        self.applyAxLimits()
        return

    def plot(self):
        for xvals,yvals,props in zip(self.xValues,self.yValues,self.crvProps):
            self.ax.plot(xvals,yvals, **props)
        return

    def plotsingle(self, idx, show_marker=True, show_rest=True):
        for xvals,yvals,props in zip(self.xValues,self.yValues,self.crvProps):
            if show_rest:
                self.ax.plot(xvals[idx:],yvals[idx:], c='grey', ls='--', lw=0.6)
            self.ax.plot(xvals[:idx+1],yvals[:idx+1], **props)
            if show_marker:
                self.ax.plot(xvals[idx],yvals[idx], c=props['color'], marker='o')
        return


### ##########################################################################
### (MATPLOTLIB) S I N G L E P L O T
### ##########################################################################
class Singleplot(Mainplot):
    """Plot a simple diagram.

    Usage
    =====
     >>> p = Singleplot(10, 7.5)
     >>> p.addDataXY(tvec, vals)
     >>> p.addDataXY(tvec, vals)
     >>> p.plotprogress()
    """

    def __init__(self, width=10, height=7.5, dpi=80,
                 bgcolor='white', projection=u'rectilinear'):
        """
        @param width: in cm
        @param height: in cm
        @param dpi: optional, dots-per-inch.
        """
        Mainplot.__init__(self, figNb=1, gridRowsCols=(75,100),
                          gridScale=0.1, dpi=dpi)
        self.sp = Subplot(name="singleplot",
                     mpInst=self, mpGridPos=(0,0), mpGridSpan=(75,100),
                     bgcolor=bgcolor, projection=projection)
        self.addSubplot(self.sp)
        return

    def addDataXY(self, xvals, yvals, lc=None, lp=None):
        """Store data for a single line.

        @param xvals: vector of x-values
        @param yvals: vector of y-values
        @param lc: line colour ("blue", "red", "black",...?)
        @param lp: other line props: a dict {prop: val}

        @Note: self.addData(ptDataDict, xId, yId, lc, lp) is equivalent to
        self.addDataXY(ptDataDict[xId], ptDataDict[yId], lc, lp)
        """
        crvNb = len(self.sp.yValues)

        ### y-data
        self.sp.yValues.append(yvals)

        ### line colour, line props
        crvProps = dict()
        if lc is not None:
            crvProps['color'] = lc
        else:
            crvProps['color'] = "C%d"%crvNb  # using the mpl.cycler
        if lp is not None:
            if 'c' in lp:
                lp['color'] = lp['c']
                del lp['c']
            crvProps.update(lp)
        self.sp.crvProps.append(crvProps)

        ### x-data
        self.sp.xValues.append(xvals)
        self.sp.computeAxLimitsAuto()
        return


###############################################################################
# defaultPlots
###############################################################################
class DefaultPlots(object):
    '''The defaultPlots class provides some straight forward methods to create
    and save matplotlib-figures. Subplots will be created in a flexible way
    using specified parameters for size, position and padding/spacing.
    Compared to the Mainplot class defaultPlots uses the handles of figures,
    axes and plots so it preserves the accessibility for all these objects via
    maplotlib api.

    Example
    =======
     >>> #create 3x2 subplots, first 2 rows in 2nd column as polar plots
     >>> fig, axs = dP.createAxes(nRows=3,nCols=2,
     >>>                          colScale=[1.5,1],
     >>>                          polar=[1,3],
     >>>                          grid='all')

    see L{exampleDefaultPlots}() for a more comprehensive example
    '''
    def __init__(self,mode='onPaper'):
        toinch = 1/2.54
        self.layout = {'mode'     :  mode,#'standalone','onPaper'
                 'paperSize':  {'width':21.*toinch,
                                 'height':1.033*29.7*toinch},
                 'landscape': False,
                 'noHeader' : False,
                 'figMargins': #relative to paperSize
                                {'top':.01,'bottom':.05,
                                 'headerHeight':.075,'footerHeight':0,
                                 'headerPad':.05, 'footerPad':0.,
                                 'left':0.05,'right':0.07},
                 'subFigSpacing':#relative to paperSize
                                {'horz':.12,'vert':.075},
                 'output':      {'format':'png','dpi':150,},
                 }
        self.font = {'family' : ['sans-serif', 'monospace',
                                  'cursive', 'fantasy', 'serif'][0],
                    'weight' : ['light', 'normal', 'medium',
                                 'semibold', 'bold', 'heavy', 'black'][1],
                    'style'  : ['normal', 'italic', 'oblique'][0],
                    'variant': ['normal', 'small-caps'][0],
                    'size'   : 10,}

    def _calAxesDims(self, axLay):
        '''Calculate the dimensions and offsets for a simple axes-grid
        @param axLay: [numberOfSubplotRows, numberOfSubplotColums]
        '''
        ps    = self.layout['paperSize'].copy()
        fm    = dict(self.layout['figMargins'].copy())
        axSp  = self.layout['subFigSpacing']

        if self.layout['noHeader']:
            fm['headerHeight'] = 0.
            fm['headerPad']  = 0.

        if self.layout['landscape']:
            ps['width'],ps['height'] = ps['height'],ps['width']

        plotAreaW = (1.-(fm['left']+fm['right'])     \
                            -(axLay[1]-1.)*axSp['horz'] )
        plotAreaH = (1.-(fm['top']+fm['bottom'])                     \
                             -(axLay[0]-1.)*axSp['vert']             \
                             -(fm['footerHeight']+fm['headerHeight'])\
                             -(fm['footerPad']+fm['headerPad']))

        axDim, axSep, axOff = [None,None], [None,None], [None,None]
        axDim[0] = .95*plotAreaW/axLay[1]   #width of each subplot/axes
        axDim[1] = .95*plotAreaH/axLay[0]   #height of each subplot/axes
        axSep[0] = axSp['horz']   #horz. space betw. subplots
        axSep[1] = axSp['vert']   #vert. space betw. subplots
        axOff[0]  = fm['left']
        axOff[1]  = fm['footerHeight']+fm['footerPad']+fm['bottom']
        figDim    = [ps['width'], ps['height']]

        return figDim, axDim, axSep, axOff

    def createAxes(self, nRows=1, nCols=1, **kwargs):
        '''This creates a simple subplotgrid of nRows X nCols. It allows the
        positioning on a 'sheet of paper' as defined in self.layout. The size,
        position and spacing between subplots can be adjusted by changing the
        values in self.layout-dictionary. There are also build in functions in
        matplotlib (e.g.: plt.subplot or mpl.gridspec), but using those impedes
        the positioning of subplots on a 'page'.

        B{Example}:
        Create a 3x2-subfigure grid, where the second column is slightly wider 
        than the first. The upper-left subPlot will hold an mapPicture/axes, the
        first two plots in right column will have a southernHemisphere-projection
        and the lower-right plot will be a harrison chart.
         >>> fig, axs = createAxes(
         >>>                   nRows=3, nCols=2,
         >>>                   colScale=[1,1.5],
         >>>                   mapaxes=[0,],
         >>>                   shemisphere=[1,3,],
         >>>                   harrison=[5,])

        @param nRows: number of rows of subplots

        @param nCols: number of columns of subplots

        @kwarg polar: indexlist of axes for PolarAxes projection.

        @kwarg invpolar: indexlist for a invPolarAxes projection.
            See L{invPolarAxes} for details.

        @kwarg shemisphere: indexlist for a southernHemisphere.
            See L{southernHemisphere} for details.
        
        @kwarg mapaxes: indexlist for a mapaxes.
            See L{MapAxes} for details.
            
        @kwarg harrison: indexlist for a harrison charts.
            See L{HarrisonAxes} for details

        @kwarg grid: display grid? ('all'(default),'any',indexlist)

        @kwarg colScale: relative size for each column (e.g. [1,.2,2])

        @kwarg rowScale: relative size for each row (e.g. [1,.2,2])

        @return: figurehandle, axeshandles (np.array of dim nRows X nCols)
        '''

        ### resolve kwargs
        # kwargs for projections
        isProjection = np.array([None for i in range(nRows*nCols)])
        projections = ['polar', 'invpolar', 'shemisphere',
                       'mapaxes', 'harrisson']
        for projection in projections:
            isProjection[kwargs.pop(projection, [])] = projection


        # kwargs for non uniform scaling of subplots
        if 'rowScale' not in kwargs:
            rowScale = np.ones((nRows,))
        else:
            rowScale = np.asarray(kwargs['rowScale'])
        if 'colScale' not in kwargs:
            colScale = np.ones((nCols,))
        else:
            colScale = np.asarray(kwargs['colScale'])

        # simple check of dimentions
        assert colScale.shape == (nCols,), \
            "ERROR: length of colScale doesn't match nCols"
        assert rowScale.shape == (nRows,), \
            "ERROR: length of colScale doesn't match nCols"

        # normalize scales, so their total-dim is 1
        colScale = nCols*(1.*colScale/colScale.sum())
        rowScale = nRows*(1.*rowScale/rowScale.sum())

        # kwars for grid on/off / can also be changed later
        isGridOn = np.array([False for i in range(nRows*nCols)])
        if 'grid' in kwargs:
            if kwargs['grid'] == 'all':
                isGridOn[:] = True
            elif kwargs['grid'] == 'any':
                pass
            else:
                isGridOn[kwargs['grid']] = True

        ### calculate dimensions of figure and subplots
        figDim, axDim, axSep, axOff = self._calAxesDims([nRows,nCols])

        ### set specified font
        mpl.rc('font', **self.font)

        ### setup figure object
        fig = plt.figure(dpi=self.layout['output']['dpi'],
                         frameon=False,)
        fig.set_size_inches(figDim[0],figDim[1])
        if self.layout['mode'].lower() != 'standlone':
            #invisible 'helper'-axes to define sheet of paper
            axbackgrnd = plt.Axes(fig,[0.,0.,1.,1.])
            axbackgrnd.axis('off')
            fig.add_axes(axbackgrnd)

        axs = np.empty((nRows,nCols),dtype=object)  # collects axes-handles

        isProjection = isProjection.reshape(nRows, nCols)[::-1,:]
        isGridOn = isGridOn.reshape(nRows, nCols)[::-1,:]
        colScale = colScale[::-1]
        rowScale = rowScale[::-1]
        
        ### create axes-grid bottom up
        for irow in range(nRows)[::-1]:
            for icol in xrange(nCols):
                #position and size off the new 'subplot'
                xL = icol*axSep[0]+axDim[0]*colScale[:icol].sum()+axOff[0]
                xU = irow*axSep[1]+axDim[1]*rowScale[:irow].sum()+axOff[1]
                yL = axDim[0]*colScale[icol]
                yU = axDim[1]*rowScale[irow]
                place = [xL, xU, yL, yU]
                
                #create subplot-axes / if polartype: a projection is needed
                ax = fig.add_axes(place, projection=isProjection[irow,icol])

                #do some refinements
                ax.tick_params(axis='both', which='major',
                               # labelsize=self.font['size']-2
                               )
                
                try:
                    ax.adjustAxes()
                except AttributeError:
                    pass

                #turn grid on/off
                if isGridOn[irow,icol]: 
                    ax.grid(linewidth=.25)

                axs[irow,icol] = ax  # append axes
        axs=axs[::-1,:]
        return fig, axs

    def createInfoHeader(self, fig, infoCell, tabWidth=None, logo=True):
        '''Create a simple header with some information in tables and BE-Logo
        @param fig: figurehandle

        @param infoCell: list with celldata for tables
        Example:
         >>> infocell= [ [ ['a','b'],['1','2'] ],
         >>>             [ ['c',],['d'] ],
         >>>             [ ['x','y'] ]]

        will create three tables::

        a   | b     c
        ----+----  ----  x | y
        1   | 2     d

        @param tabWidth: list/numpy array with relative size of each single
        table defined in infocell (example [2,1,0.7])

        @param logo: set logo True/False

        @return: a tuple::
                    hTabs = list of table-handles
                    axTabs = list of axes-handles 'holding' tables
                    axLogo = axes-handle holding logo
        '''
        import matplotlib.image as mpimg
        fm = self.layout['figMargins']

        if self.layout['noHeader']:
            msg('The value for noHeader is set True. Header creation skipped!')
            axLogo, axTabs, hTabs = None, None, None
            logoDim = [0,0]
        else:
            if logo:
                # show BeckEngineering-logo in the upper-right corner
                beLogo  = mpimg.imread(logoPath)
                logoDim = [fm['headerHeight']*beLogo.shape[1]/beLogo.shape[0],
                           fm['headerHeight']]
                axLogo = fig.add_axes([1-logoDim[0]-fm['right'],
                                       1-logoDim[1]-fm['top']]+logoDim)
                axLogo.imshow(beLogo)
                axLogo.axis('off')
            else:
                logoDim = [0,fm['headerHeight'],]

            # subdivide the left space into (axes) fields and show tables for
            # info cells
            horzSep = .02
            bWidth = 1-logoDim[0]-fm['right']-fm['left']-horzSep
            if not tabWidth:  # apply equal tablewidth
                tabWidth = bWidth*np.ones((len(infoCell),))/len(infoCell)\
                           - horzSep
            else:
                assert len(tabWidth)==len(infoCell), \
                'Number of tables and dimension of tableWidth are different'
                tabWidth = 1.*np.array(tabWidth)
                tabWidth = tabWidth/tabWidth.sum()
            tabWidth = (bWidth*tabWidth - horzSep)

            #np.array/list to collect tablehandles
            axTabs = np.empty(len(tabWidth),dtype='object')
            hTabs = np.empty(len(tabWidth),dtype='object')
            for iTab,tabW in enumerate(tabWidth):
                #calculate position and size
                pos = [fm['left'] + tabWidth[0:iTab].sum() + iTab*horzSep,
                       1 - logoDim[1]- fm['top'],
                       tabWidth[iTab],logoDim[1]]
                #create axes where the table is set
                axhTab = fig.add_axes(pos)
                axhTab.axis('off')
                #create table and set celltext
                hTab = axhTab.table(cellText=infoCell[iTab],
                                    loc="center")
                #do some refinements
                hTab.set_fontsize(self.font['size']-4)
                #add some vertical space to each cell
                table_props = hTab.properties()
                table_cells = table_props['child_artists']
                for cell in table_cells:
                    cell.set_height(.275)

                axTabs[iTab] = axhTab
                hTabs[iTab]  = hTab
        return hTabs, axTabs, axLogo

    @staticmethod
    def changeHeaderTableCell(hTab, cellIds, newVals):
        ''' Updates values in table cellwise

        @param hTab: handle of table

        @param cellIds: list of tuples or single tuple (iRow,iCol) holding the
            cells to be updated

        @param newVals: list of values or single value to be set to cellIds
        '''
        if not type(hTab)==mpl.table.Table:
            raise TypeError('handleHeadTab needs to be a valid table-handle'
                            ' You passed: type(handleHeadTab) =%s' % type(hTab))
        if not isinstance(cellIds,list):
            cellIds = [cellIds]
        if not all((isinstance(cId,tuple) and len(cId)==2) for cId in cellIds):
            raise ValueError('The cellIds must come as a list of '+
                             'tuples, each of len==2. e.g. [(0,1),(2,2)]')

        if not isinstance(newVals, collections.Iterable) or \
           isinstance(newVals,str):
            newVals = [newVals]

        if not (len(cellIds) == len(newVals)):
            raise IndexError('newVals must be the same length as cellIds')

        for cellId,val in zip(cellIds,newVals):
            hTab._cells[cellId]._text.set_text(val)

    @staticmethod
    def setCellWidth(tabh, ratios):
        '''Sets the relative column-sizes of tabular

        Example:
         >>> setCellWidth(tabh,[1,2,3])

        @param tabh: reference of table-object
        @param ratios: unnormalized list of column ratios
        '''
        cells = tabh.get_celld()
        nCols = max([key[1] for key in cells.keys()])+1
        nRows = max([key[0] for key in cells.keys()])+1
        if not len(ratios) == nCols:
            raise IndexError('Len of ratio-list must match the number of columns')
        width = sum([cells[(0,ii)].get_width() for ii in range(nCols)])
        #normalize ratios
        ratios = np.asarray(ratios,dtype=float)
        ratios = ratios/ratios.sum()
        #set new cell width
        for ii,rat in enumerate(ratios):
            for kk in range(nRows):
                cells[(kk,ii)].set_width(width*rat)

                
    @staticmethod
    def autoscaleAxes(axhs):
        '''Rescale axes (eg. after plot data is updated)
        @param axhs: single axes handle or list/numpy array holding handles
        @warning: this may fail for polarlike axes, see:
            U{https://github.com/matplotlib/matplotlib/issues/7130} (jan 2018)
            Workaround: get all y-data of plot objects* in axh using getMaxR
            and set rmax manually. Different plot-types have different methods
            to get the data, so append type+method to getMaxR if needed. By now

            * matplotlib.lines.Line2D (created by 'plot')
            * matplotlib.collections.PathCollection (created by 'scatter')

            are implemented.
        '''
        if not isinstance(axhs, collections.Iterable): 
            axhs = [axhs]

        polarPrjs = ['PolarAxes','invPolarAxes','southernHemisphere']
        polarPrjs = polarPrjs + ['%sSubplot'%prj for prj in polarPrjs]
        def getMaxR(axh):
            '''Get all y-data in plot objects* of axes and return max value.
                Objects:
                * matplotlib.lines.Line2D (created by 'plot')
                * matplotlib.collections.PathCollection (created by 'scatter')
            '''
            rVals = []
            for child in axh.get_children():
                if type(child) == mpl.lines.Line2D:
                    rVals.append(child.get_ydata())
                elif type(child) == mpl.collections.PathCollection:
                    rVals.append(child.get_offsets()[:,1])
            if rVals:
                maxval = np.hstack(rVals).max()
            else:
                msg("Warning! Can't find plot objects in polar like axes."+
                    "Maxval for r can not be determined and will be set to 1")
                maxval = 1
            return maxval

        for axh in np.asarray(axhs).ravel():
            axesType = str(type(axh)).split('.')[-1][:-2]
            if axesType in ['AxesSubplot','Axes']:   #regular xy-axes
                axh.relim()
                axh.autoscale_view()
            elif axesType in polarPrjs: #polartyp axes
                axh.set_rmax(getMaxR(axh))

    
    @staticmethod
    def estimateNiceAxisLims(array, mod = 1.):
        '''Rounds the min and max values of an array to next nth of mod.
        estimateAxisLims([11,12,22],10) >> (10.0, 30.0)
        estimateAxisLims(np.pi,.01) >> (3.14, 3.15)
        '''
        arrayMin = np.floor(np.nanmin(1.*np.asarray(array))/mod)*mod
        arrayMax = np.ceil(np.nanmax(1.*np.asarray(array))/mod)*mod
        return arrayMin, arrayMax
    
    @staticmethod
    def getRangeString(array, allDims=False, fmt='%.0f'):
        def rangeString(the_range):
            if the_range[0]==the_range[1]:
                return fmt%the_range[0]
            else:
                return (fmt+' to '+fmt)%(the_range[0],the_range[1])

        array = np.atleast_2d(array)
        if allDims:
            rangeStr = [rangeString([array.min(),array.max()])]
        else:
            rangeStr = []
            for ii in range(array.shape[1]):
                ranges = ( array[:,ii].min(),array[:,ii].max() )
                rangeStr.append(rangeString(ranges))
        return rangeStr


    @staticmethod
    def setByIdx(axs, axesIdxs, objLabel,
                 pattern=None,valList=None,addargs=''):
        '''Set axes-parameters for a column of subplots (axs[:,colIdx]).
        It iterates all axes objects in column and applies an available set_*-
        method:
        ax.set_#objLabel(#arg,#addargs)
        where #arg can be interpreted independently for each row.
        Example:
         >>> setByIdx(axs, [0,1,4], 'title', '%s Spam',
         >>>             valList = ['','Eggs and', 'Spam and'],
         >>>             addargs = 'loc="left"')

        This will create 3 independent (sub)titles for the axes in second
        column which are aligned in the left corner of axes.

        @param axs: numpy array holding axes handles
        @param axesIdxs: flat axes-index where set_methode is applied to
        @param objLabel: available object label of axes which can be changed
        via set methode (see https://matplotlib.org/api/axes_api.html)
        @param pattern: argument of set_method. Can conatain '%'-replacement-
        pattern which will be interpreted for each n-th row according to the
        n-th value of valList
        @param valList: List of values to be for each row according to pattern
        @param addargs: additional constant arguments as string. In most cases
        this is redundant to an appropriate definition of pattern, but usually
        the usage of addargs makes it easier to debug the pattern string.
        '''
        def fstr(value):
            if isinstance(value, basestring):   out = '"%s"'%value
            else:                               out = str(value)
            return out

        axesIdxs = np.asarray(axesIdxs)
        for i,axesIdx in enumerate(axesIdxs):
            ax = axs.ravel()[axesIdx]
            if valList and pattern:
                evalStr = r'ax.set_%s'%(objLabel) + \
                          '( %s%%(%s),%s )'%(fstr(pattern),
                                          ','.join(map(fstr,valList[i])),
                                          addargs,)
            elif valList:
                evalStr = r'ax.set_%s( %s,%s )'%(objLabel,
                                               ','.join(map(fstr,valList[i])),
                                               addargs,)
            else:
                evalStr = r'ax.set_%s( %s,%s )'%(objLabel,pattern,addargs)
            try:
                eval(evalStr)
            except:
                msg('ERROR: The objectLabel and valList you passed gives: %s'%\
                    evalStr)

    @classmethod
    def setByColIdx(cls, axs, colIdx, objLabel, pattern, valList=None, addargs=''):
        '''Set axes-parameters for a column of subplots (axs[:,colIdx]).
        See setByIdx for detailed information.
        '''
        assert axs.shape[1] >=colIdx+1, \
            'Figure does not contain an axis-colunm with indx %d'%colIdx
        assert (valList == None) or (axs.shape[0]==len(valList)),\
            'Length of valuelist does not match the number of subplots columns'

        #flat index list of axes in column
        axesIdxs = np.arange(axs.size).reshape(axs.shape)[:,colIdx]
        cls.setByIdx(axs,axesIdxs,objLabel,pattern,valList,addargs)

        
    @classmethod   
    def setByRowIdx(cls, axs, rowIdx, objLabel, pattern, valList=None, addargs=''):
        '''Set axes-parameters for a row of subplots (axs[rowIdx,:]).
        See setByIdx for detailed information.
        '''
        assert axs.shape[0] >=rowIdx+1, \
            'Figure does not contain an axis-rows with indx %d'%rowIdx
        assert (valList == None) or (axs.shape[1]==len(valList)),\
            'Length of valuelist does not match the number of subplots in col'

        #flat index list of axes in row
        axesIdxs = np.arange(axs.size).reshape(axs.shape)[rowIdx,:]
        cls.setByIdx(axs,axesIdxs,objLabel,pattern,valList,addargs)

    @staticmethod
    def getAxesExtent(axhs, pad=0.0,includeLabels=True, includeTitles=True):
        from matplotlib.transforms import Bbox
        import collections
        """Get the full extent/bounding box of one ore more axes objects
        including their labels, tick labels and/ore titles.
        @param axhs: single axes handle or list of axes handles to be saved
        @param pad: extra spacing (rel to bonding box of axhs)
        @param includeLabels: if True, labels where used to determine boundig box
        @param includeTitles: if True, titles where used to determine boundig box
        """
        if not isinstance(axhs, collections.Iterable): axhs = [axhs,]

        items = []
        for axh in axhs:
            axh.figure.canvas.draw()
            items += axh.get_xticklabels() + axh.get_yticklabels()
            items += [axh,]
            if includeLabels:
                items += [axh.xaxis.label, axh.yaxis.label]
            if includeTitles:
                items += [axh.title,]
                bboxes = []
        for item in items:
            #Axes-objects which are not drawn will raise a RuntimeError
            try :
                bb = item.get_window_extent()
                #append only if bb has non zero extent
                if bb.height*bb.width != 0:
                    bboxes.append(bb)
            except RuntimeError:
                pass
        bbox = Bbox.union(bboxes)
        return bbox.expanded(1.0 + pad, 1.0 + pad)

    @classmethod
    def saveSingleAxes(cls, figh, axhs, fileName,
                       pad=0.0,
                       includeLabels=True,
                       includeTitles=True):
        ''' Saving specified set of axes oject(s) in figh into file.
        If multiple axes handles are given, then the saved area is determined by
        their (union) bounding box.
        @param figh: handle of existing figure
        @param axhs: single axes handle or list of axes handles to be saved
        @param pad: extra spacing (rel to bonding box of axhs)
        @param includeLabels: if True, labels where used to determine boundig box
        @param includeTitles: if True, titles where used to determine boundig box
        '''

        extent = cls.getAxesExtent(axhs,pad,
                              includeLabels = includeLabels,
                              includeTitles = includeTitles)
        extent = extent.transformed(figh.dpi_scale_trans.inverted())
        figh.savefig(fileName, bbox_inches=extent)

    @staticmethod
    def setDepthAxes(axh,z0,
                         zLabel = 'Depth [m]',
                         tickLabelFmt = '%d',
                         onlyPosTicks = True):
        ''' Takes a (left side) z-axis (axh) flips it, adds an offset z0 and
        creates a new right side z-axis
        @param axh: existing axes handle
        @param z0: offset where axh.zaxes has its zero depth
        @param zLabel: specifies label of the new z-axis
        @param tickLabelFmt: eg: %d, %.2f,...
        @param onlyPosTicks: if True, only ticks >=0 will be set
        @return: axes handle of new z-axis
        '''
        #create new axes on right side
        axhRight = axh.twinx()
        #get z-ticks of axh and calculate offset
        ylimLeft = axh.get_ylim()
        off = ylimLeft[1]-z0
        #reverse new axis
        axhRight.invert_yaxis()
        #set limits + offset to new axes
        axhRight.set_ylim(ylimLeft[1]-off,ylimLeft[0]-off)
        if onlyPosTicks:
            ticks = np.asarray(axhRight.get_yticks())
            axhRight.set_yticks(ticks[ticks >= 0])
        #set label
        axhRight.set_ylabel(zLabel)
        axhRight.tick_params(axis='both', which='major', labelsize=8)

        return axhRight


def exampleDefaultPlots(savePath='.'):
    import matplotlib.image as mpimg
    plt.close("all")
    dP = DefaultPlots()

    #1) prepare plots / layout
    dP.layout['mode']       = 'onPaper'
    dP.layout['landscape']  = False
    #create 3x2 subplots
    fig, axs = dP.createAxes(nRows=3,nCols=2,
                             colScale=[1.5,1],
                             polar=[1,3],
                             grid='all')

    #2) setup header
    infoCell = [np.array([['Project:','MYPRJ2018'],
                          ['Date:','12/31/2018'],
                          ['Info:','somthing usefull']]),
                np.array([['East [m]:','1000'],
                          ['North [m]:','1001'],
                          ['Depth (range) [m]:','100 - 3000']])]

    tabs,_,_ = dP.createInfoHeader(fig,infoCell,tabWidth=[1.5,1])

    #3) create some data to plot
    t = np.arange(0.0, 6.0, 0.1)*np.pi
    s = 1+np.sin(t)
    h = 0.1*np.random.randn(2,10000)
    h[1,:] = 0.25*h[1,:]+0.1

    #4) plot data using common matplotlib/pyplot commands
    #a) first subplot column
    ph = [None,None,None,None]
    ph[0], = axs[0,0].plot(t,np.sin(t),'og-')
    ph[1], = axs[0,0].plot(t,np.random.randn(t.shape[0])+1,'^:')

    axs[1,0].hist(h[0,:],50)
    axs[2,0].hist(h[1,:],50)
    #b) second subpolt column
    axs[0,1].imshow(mpimg.imread(logoPath))
    axs[0,1].axis('off')
    ph[2], = axs[1,1].plot(t,np.sin(t),'og')
    ph[3]  = axs[2,1].scatter(.5*t,s,s=15*s,c=t)

    #5) some stupid adjustments showing robust usage of setByIdx-methods
    #setting x labels for first column
    dP.setByColIdx(axs,0,'xlabel',
                   pattern='$\sigma_{%s}$ [MPa]',
                   valList=[['I'],['II'],[3]])#
    #setting titles for first row
    dP.setByRowIdx(axs,0,'title',
                   pattern='title: %s',
                   valList=[['a plot'],['png via imshow']])
    #unify axis limits of 2nd and 4th x-axis
    xlims = dP.estimateNiceAxisLims(h,.1)
    dP.setByIdx(axs,[2,4],'xlim',
                valList=[xlims,xlims])

    fig.savefig(os.path.join(savePath,'foo.png'), bbox_inches='tight')

    #6) example for incremental update for some plots
    for i in range(2):
        #create data to update 1st plot in 1st axes
        dataOld    = ph[0].get_data()
        dataNew    = [[],[]]
        dataNew[0] = t + dataOld[0][-1] + t[1] #last timestamp + one timesample
        dataNew[1] = np.sin(dataNew[0])+1+i
        #update data of plot via plothandle
        ph[0].set_data(np.c_[dataOld,dataNew])  #update data in plot object
        ph[3].set_offsets(np.c_[dataNew[0],dataNew[1]])  #append data to scatter plot

        dP.autoscaleAxes([axs[0,0],axs[2,1]])

        #update headerdata
        dP.changeHeaderTableCell(tabs[0],(2,1),'Run %d'%i)
        currentPosRangeStr = dP.getRangeString(100*np.random.randn(100,3))
        dP.changeHeaderTableCell(tabs[1],[(0,1),(1,1),(2,1)],
                                         currentPosRangeStr)

        fig.savefig(os.path.join(savePath,'foo%d.png'%i), bbox_inches='tight')
        dP.saveSingleAxes(fig,axs[2,1],os.path.join(savePath,'bar%d.png'%i))

###############################################################################
#'Pseudo'-projections
###############################################################################
#
#Matplotlib allows to create custom projections to plot abitrary non-cartesian
#data. This is defined on axes level and all plot methods can be used as usual.
#(see https://matplotlib.org/devel/add_new_projection.html)
#But for some reasons (singularities?) i wasn't able yet to create a PolarAxes
#where the radial axis having a reversed direction (90deg in center to 0deg on
#outer bound). Thus the following axes classes / projections just reuse the
#PolarAxes class and append some transformed plotting methods (plotT, scatterT)

from matplotlib.axes import Axes
from matplotlib.projections import PolarAxes
from matplotlib.projections import register_projection

class invPolarAxes(PolarAxes):
    ''' A new axes-type depending on PolarAxes. It just appends the plotmethods
    plotT and scatterT to plot the axis-direction for dip from inside out. By
    default the r-axes of PolarAxes is pointing inside-out. There is no simple
    way to change this behavior yet. Manually swiching the direction and
    renaming the axis-labels seems to be the easiest way by now.

    @warning: this works just for basic plot routines (plot,scatter,contour,
        contourf). If you want to use other polt-methods or set_data you need
        to transform the input data before:
        >>> myPlotHandle.set_data(soutPoleAxes.transformData(strike,dip))
    '''
    name = 'invpolar'

    def transformData(self,strike,dip):
        '''
        @param strike: Vector of strike-direction in radians [0,2pi]
        @param dip: Vector of dip-direction in radians [0,pi/2]
        '''
        dip      = np.abs(np.asarray(dip)-0.5*np.pi)
        return strike, dip


    def plotT(self, strike, dip, *args, **kwargs):
        """Identical to plot - for legacy"""
        self.plot(strike, dip, *args, **kwargs)

    def scatterT(self, strike, dip, *args, **kwargs):
        """Identical to scatter - for legacy"""
        self.scatter(strike, dip, *args, **kwargs)

    def plot(self, strike, dip, *args, **kwargs):
        """ Wraps plot-method and applies the transformation"""
        strike, dip = self.transformData(strike, dip)
        return PolarAxes.plot(self,strike, dip, *args, **kwargs)

    def scatter(self, strike, dip, *args, **kwargs):
        """ Wraps scatter-method and applies the transformation"""
        strike, dip = self.transformData(strike, dip)
        return PolarAxes.scatter(self, strike, dip, *args, **kwargs)



    def adjustAxes(self, aStepStrike=15, aStepDip=15):
        '''Sets a inverse direction of dip-axes and a nice look. This includes:
         - all angles are displayed in degree (see warning!)
         - 0degree points upwards
         - positive rotation of strike is clockwise
         - dip-axis (here r-axis) ranges from 90deg (inside) to 0deg (outside)
        @param aStepStrike: sets delta angle (deg) for strike-ticks/gridlines
        @param aStepDip: sets delta angle (deg) for dip-ticks/gridlines
        @warning: this axes is just a 'lookalike' and includes no trans-
        formation. This means that the data needs to be transformed (see Class
        southPoleAxes).
        !!! before plotting --> inpStrike=deg2rad(strike);inpDip=np.abs(Dip-90)
        '''
        self.set_theta_zero_location('N')  #North at 12:00
        self.set_theta_direction(-1)       #Rotation clockwise

        #Ticks and labels for strike-angle
        strikeLabels = np.array([str(val)+u"\u00b0" \
                               for val in range(0, 360, aStepStrike)] )
        for deg,direct in zip([0,90,180,270],['N','E','S','W']):
            strikeLabels[strikeLabels == str(deg)+u"\u00b0"] = direct

        self.set_thetagrids(np.arange(0, 360, aStepStrike),
                        labels=strikeLabels)

        #Limits, ticks and labels for dip-angle
        self.set_rlim(0, np.pi/2)
        self.set_yticks(np.arange(0, 1.05*np.pi/2, 2*np.pi*aStepDip/360))
        self.set_yticklabels([str(val)+u"\u00b0"\
                          for val in range(90, -1,-aStepDip)])
        self.yaxis.set_tick_params(labelsize=5)
        self.tick_params(axis='x', which='major', pad=0)

register_projection(invPolarAxes)


class southernHemisphere(PolarAxes):
    '''A new axes-type depending on PolarAxes. It creates a southern-hemi-
    sphere-projection-like look and appends the plotmethods plotT and scatterT
    to plot the axis-direction for dip from inside out. By default the r-axes
    of PolarAxes is pointing inside-out. There is no simple way to change this
    behavior yet. Manually swiching the direction and renaming the axis-labels
    seems to be the easiest way by now.

    @warning: this works just for basic plot routines (plot,scatter,contour,
        contourf). If you want to use other polt-methods or set_data you need
        to transform the input data before:
        >>> myPlotHandle.set_data(soutPoleAxes.transformData(strike,dip))
    '''
    name = 'shemisphere'

    def transformData(self,strike,dip):
        '''
        @param strike: Vector of strike-direction in radians [0,2pi]
        @param dip: Vector of dip-direction in radians [0,pi/2]
        '''
        dip = 2*np.sin( (0.5*np.pi - np.asarray(dip))/2 )
        return strike, dip

    def plotT(self, strike, dip, *args, **kwargs):
        """Identical to plot - for legacy"""
        self.plot(strike, dip, *args, **kwargs)

    def scatterT(self, strike, dip, *args, **kwargs):
        """Identical to scatter - for legacy"""
        self.scatter(strike, dip, *args, **kwargs)

    def plot(self, strike, dip, *args, **kwargs):
        """ Reimplemtation of plot-methode (includes transformation)"""
        strike, dip = self.transformData(strike, dip)
        return PolarAxes.plot(self, strike, dip, *args, **kwargs)

    def scatter(self, strike, dip, *args, **kwargs):
        """ Reimplemtation of scatter-plot (includes transformation)"""
        strike, dip = self.transformData(strike, dip)
        return PolarAxes.scatter(self, strike, dip, *args, **kwargs)

    def contourf(self, strike, dip, *args, **kwargs):
        """ Reimplemtation of contourf-plot (includes transformation)"""
        strike, dip = self.transformData(strike, dip)
        return PolarAxes.contourf(self, strike, dip, *args, **kwargs)


    def contour(self, strike, dip, *args, **kwargs):
        """ Reimplemtation of contour-plot (includes transformation)"""
        strike, dip = self.transformData(strike, dip)
        return PolarAxes.contour(self, strike, dip, *args, **kwargs)

    def adjustAxes(self, aStepStrike=15, ticksDip=range(0,90,15)):
        '''Sets a southern-hemisphere-projection-like look for a regular polar-
        axes. This includes:
         - all angles are displayed in degree (see warning!)
         - 0degree points upwards
         - positive rotation of strike is clockwise
         - dip-axis (here r-axis) ranges from 90deg (inside) to 0deg (outside)

        @param aStepStrike: sets delta angle (deg) for strike-ticks/gridlines
        @param ticksDip:    list of angles (deg) for dip-ticks/gridlines
        @warning: this axes is just a 'lookalike' and includes no transformation
        This means that the data needs to be transformed (see class
        L{southPoleAxes})
        before plotting --> inpStrike=deg2rad(strike), inpDip=np.abs(Dip-90)
        '''
        self.set_theta_zero_location('N')  # North at 12:00
        self.set_theta_direction(-1)       # Rotation clockwise

        #Ticks and labels for strike-angle
        strikeLabels = np.array(
            [str(val)+u"\u00b0" for val in range(0, 360, aStepStrike)])
        for deg,direct in zip([0,90,180,270],['N','E','S','W']):
            strikeLabels[strikeLabels == str(deg)+u"\u00b0"] = direct

        self.set_thetagrids(np.arange(0, 360, aStepStrike),
                            labels=strikeLabels)

        #Limits, ticks and labels for dip-angle
        _,angs = self.transformData(None,np.deg2rad(ticksDip))
        self.set_yticks(angs)
        self.set_yticklabels([str(val)+u"\u00b0" for val in ticksDip])
        self.set_rlim(0, .91*np.pi/2)
        self.yaxis.set_tick_params(labelsize=5)
        self.tick_params(axis='x', which='major', pad=0)

register_projection(southernHemisphere)


class HarrissonAxes(Axes):
    '''
    Creates the Background for Harrisson charts.

    Default-settings are stored in::
        cls._harrissonRegions -- quadratic parameters according to harrisson
        cls._harrissonLabels -- description strings for harrison regions
        cls._harrissonBackgroundCmap -- colormap for background (default grey)

    The created plot/annotation objects get stored in::
        self._harrissonBackgroundObjs

    Example:
     >>> import numpy as np
     >>> import matplotlib.pyplot as plt
     >>> from bae.plot_01.mpl_core import HarrisonAxes
     >>> ax = plt.subplot(111, projection='harrisson')
     >>> ax._harrissonLabels[-1]= 'BOOM!' #change last label to BOOM
     >>> ax._harrissonBackgroundCmap= plt.get_cmap('Reds')
     >>> ax.updateHarrisson()
     >>> ax.plot(np.arange(7), .25*np.arange(7)*(1 + .2*np.random.rand(7)))
     >>> plt.show()
    '''

    name = 'harrisson'

    _harrissonRegions = [
        [0, 0, 0],                  # none
        [-0.3000, -0.0500, 0.4500], # very slightly
        [-0.2654, -0.0404, 0.7466], # slightly
        [-0.1377, -0.0071, 1.5465], # moderate
        [-0.0669, -0.0009, 3.0026], # severe
#        [0, 0, inf],                # very severe
        ]

    _harrissonLabels = [
        'Negligible damage',
        'Very slight damage',
        'Slight damage',
        'Moderate to Severe damage',
        'Severe to Very severe damage',
        ]

    _harrissonBackgroundCmap = truncateColormap(
        plt.get_cmap('Greys'),
        minval=0.0, maxval=.35, n=100)

    def __init__(self, *args, **kwargs):
        super(Axes, self).__init__(*args,**kwargs)
        self.updateHarrisson()

    def updateHarrisson(self):

        pObjKeys = ['edges', 'filledAreas', 'annotations']

        # delete existing harrisionRegions
        for key in pObjKeys:
            try:
                for obj in self._harrissonBackgroundObjs[key]:
                    try:
                        obj.remove()
                    except:
                        pass
            except:
                pass

        self._harrissonBackgroundObjs = dict((p,[]) for p in pObjKeys)

        def xy(p):
            # to get a clean cut-off, calculate root for each harrisson region
            try:
                x0 = np.roots(p).max()
                x = np.linspace(0., x0, 100)
                y = p[0]*x**2 + p[1]*x + p[2]
            except ValueError:
                x, y = 0, p[-1]

            return x, y

        # plot in reversal order to simplify fill_between-call
        regs = self._harrissonRegions[::-1]
        labs = self._harrissonLabels[::-1]

        #last color will be set as axes background -> independent from limits
        cols = self._harrissonBackgroundCmap(np.linspace(0.0, 1, len(regs)))
        cols = cols.tolist()
        bcol = cols.pop()
        cols = cols[::-1]

        self.set_facecolor(bcol)

        # loop the harrisson-regions
        vals = zip(regs[1:], regs[:-1], cols, labs)
        for ii, (pl, pu, col, lab) in enumerate(vals):
            x, y = xy(pu)
            fill = self.fill_between(x, y, 0, color=col)
            edge = self.plot(x, y, 'k:')
            try:
                text = CurvedText(x=x, y=y, text=lab,
                                  fontsize='xx-small',
                                  color=(0.4, 0.4, 0.4, 1.0),
                                  va='bottom', axes=self)
            except Exception:
                text = None

            #append plotObjects to edit them later on
            self._harrissonBackgroundObjs['edges'].append(edge)
            self._harrissonBackgroundObjs['filledAreas'].append(fill)
            self._harrissonBackgroundObjs['annotations'].append(text)

        # preadjust lower plotBounds
        self.set_xlim(left=0)
        self.set_ylim(bottom=0)

register_projection(HarrissonAxes)


class MapAxes(Axes):
    '''
    Creates a map-Background from pngFile.

    Example:
        >>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> from bae.plot_01.mpl_core import MapAxes
        >>> ax = plt.subplot(111, projection='mapaxes')
        >>> ax.setMap('./myFunkyMap.png', [(-123, 2345), (1456, 3456)])
        >>> ax.scatter(xPos, yPos, s=10, c='r')
        >>> plt.show()

    '''

    name = 'mapaxes'

    def __init__(self, *args, **kwargs):

        self.mapFile = kwargs.pop('mapFile', None)
        self.mapExtend = kwargs.pop('mapExtend', None)

        super(Axes, self).__init__(*args, **kwargs)

        self.updateMap()

    def setMap(self, mapFile, mapExtend):
        self.mapFile = mapFile
        self.mapExtend = mapExtend
        self.updateMap()

    def updateMap(self):
        if self.mapFile and self.mapExtend:
            img = plt.imread(self.mapFile)
            self.mapPlot = self.imshow(img, zorder=0, extent=self.mapExtend)
            self.autoscale(enable=True, axis='x', tight=True)
            self.autoscale(enable=True, axis='y', tight=True)
            
register_projection(MapAxes)

###############################################################################
from matplotlib import text as mtext
import math


class CurvedText(mtext.Text):
    """
    A text object that follows an arbitrary curve defined by a series of
    x and y values.

    Usage:
     >>> import numpy as np
     >>> import matplotlib.pyplot as plt
     >>> from bae.plot_01.mpl_core import CurvedText
     >>> axh = plt.subplot(111)
     >>> x = np.linspace(0, 2*np.pi, 100),
     >>> y = np.sin(np.linspace(0, 2*np.pi, 100))
     >>> axh.plot(x,y)
     >>> someRegularTextProps = {color='b', va = 'bottom',}
     >>> ct = CurvedText(x, y, 'boring sin', axh, **someRegularTextProps )
     >>> axh.show()

    Found here:
    U{https://stackoverflow.com/questions/19353576/curved-text-rendering-in-matplotlib}
    """
    def __init__(self, x, y, text, axes, **kwargs):
        super(CurvedText, self).__init__(x[0],y[0],' ', **kwargs)

        axes.add_artist(self)

        ##saving the curve:
        self.__x = x
        self.__y = y
        self.__zorder = self.get_zorder()

        ##creating the text objects
        self.__Characters = []
        for c in text:
            if c == ' ':
                ##make this an invisible 'a':
                t = mtext.Text(0,0,'a')
                t.set_alpha(0.0)
            else:
                t = mtext.Text(0,0,c, **kwargs)

            #resetting unnecessary arguments
            t.set_ha('center')
            t.set_rotation(0)
            t.set_zorder(self.__zorder +1)

            self.__Characters.append((c,t))
            axes.add_artist(t)


    ##overloading some member functions, to assure correct functionality
    ##on update
    def set_zorder(self, zorder):
        super(CurvedText, self).set_zorder(zorder)
        self.__zorder = self.get_zorder()
        for c,t in self.__Characters:
            t.set_zorder(self.__zorder+1)

    def draw(self, renderer, *args, **kwargs):
        """
        Overload of the Text.draw() function. Do not do
        do any drawing, but update the positions and rotation
        angles of self.__Characters.
        """
        self.update_positions(renderer)

    def update_positions(self,renderer):
        """
        Update positions and rotations of the individual text elements.
        """

        #preparations

        ##determining the aspect ratio:
        ##from https://stackoverflow.com/a/42014041/2454357

        ##data limits
        xlim = self.axes.get_xlim()
        ylim = self.axes.get_ylim()
        ## Axis size on figure
        figW, figH = self.axes.get_figure().get_size_inches()
        ## Ratio of display units
        _, _, w, h = self.axes.get_position().bounds
        ##final aspect ratio
        aspect = ((figW * w)/(figH * h))*(ylim[1]-ylim[0])/(xlim[1]-xlim[0])

        #points of the curve in figure coordinates:
        x_fig,y_fig = (
            np.array(l) for l in zip(*self.axes.transData.transform([
            (i,j) for i,j in zip(self.__x,self.__y)
            ]))
        )

        #point distances in figure coordinates
        x_fig_dist = (x_fig[1:]-x_fig[:-1])
        y_fig_dist = (y_fig[1:]-y_fig[:-1])
        r_fig_dist = np.sqrt(x_fig_dist**2+y_fig_dist**2)

        #arc length in figure coordinates
        l_fig = np.insert(np.cumsum(r_fig_dist),0,0)

        #angles in figure coordinates
        rads = np.arctan2((y_fig[1:] - y_fig[:-1]),(x_fig[1:] - x_fig[:-1]))
        degs = np.rad2deg(rads)


        rel_pos = 10
        for c,t in self.__Characters:
            #finding the width of c:
            t.set_rotation(0)
            t.set_va('center')
            bbox1  = t.get_window_extent(renderer=renderer)
            w = bbox1.width
            h = bbox1.height

            #ignore all letters that don't fit:
            if rel_pos+w/2 > l_fig[-1]:
                t.set_alpha(0.0)
                rel_pos += w
                continue

            elif c != ' ':
                t.set_alpha(1.0)

            #finding the two data points between which the horizontal
            #center point of the character will be situated
            #left and right indices:
            il = np.where(rel_pos+w/2 >= l_fig)[0][-1]
            ir = np.where(rel_pos+w/2 <= l_fig)[0][0]

            #if we exactly hit a data point:
            if ir == il:
                ir += 1

            #how much of the letter width was needed to find il:
            used = l_fig[il]-rel_pos
            rel_pos = l_fig[il]

            #relative distance between il and ir where the center
            #of the character will be
            fraction = (w/2-used)/r_fig_dist[il]

            ##setting the character position in data coordinates:
            ##interpolate between the two points:
            x = self.__x[il]+fraction*(self.__x[ir]-self.__x[il])
            y = self.__y[il]+fraction*(self.__y[ir]-self.__y[il])

            #getting the offset when setting correct vertical alignment
            #in data coordinates
            t.set_va(self.get_va())
            bbox2  = t.get_window_extent(renderer=renderer)

            bbox1d = self.axes.transData.inverted().transform(bbox1)
            bbox2d = self.axes.transData.inverted().transform(bbox2)
            dr = np.array(bbox2d[0]-bbox1d[0])

            #the rotation/stretch matrix
            rad = rads[il]
            rot_mat = np.array([
                [math.cos(rad), math.sin(rad)*aspect],
                [-math.sin(rad)/aspect, math.cos(rad)]
            ])

            ##computing the offset vector of the rotated character
            drp = np.dot(dr,rot_mat)

            #setting final position and rotation:
            t.set_position(np.array([x,y])+drp)
            t.set_rotation(degs[il])

            t.set_va('center')
            t.set_ha('center')

            #updating rel_pos to right edge of character
            rel_pos += w-used


###############################################################################
def make_get_proj(self, rx, ry, rz):
    '''
    MatplotLib still dosn't provide a setAspectRatio in 3D.
    This workaround can be found here:
    U{https://frageit.de/questions/38149736/3d-plot-aspect-ratio-matplotlib}

    Returns a variation on :func: '~mpl_toolkit.mplot2d.axes3d.Axes3D.getproj'
    that makes the box aspect ratio equal to *rx:ry:rz*, using an axes object
    *self*.

    Usage:
     >>> from mpl_toolkits.mplot3d import Axes3D
     >>> import matplotlib.pyplot as plt
     >>> from bae.plot_01.mpl_core import make_get_proj

     >>> fig = plt.figure()
     >>> ax = fig.add_subplot(111, projection='3d')
     >>> data = np.random.randn(30,3)
     >>> ax.scatter(data[:,0],data[:,1],data[:,2],c=data[:,2],s=10)

     >>> ax.get_Proj = make_get_proj(ax, 1, 1, 1)
     >>> ax.set_aspect(1.0)
    '''
    from mpl_toolkits.mplot3d import proj3d
    rm = max(rx, ry, rz)
    kx = rm / rx
    ky = rm / ry
    kz = rm / rz

    # Copied directly from mpl_toolkit/mplot3d/axes3d.py. New or modified lines
    # are marked by ##
    def get_proj():
        relev, razim = np.pi * self.elev/180, np.pi * self.azim/180

        xmin, xmax = self.get_xlim3d()
        ymin, ymax = self.get_ylim3d()
        zmin, zmax = self.get_zlim3d()

        # transform to uniform world coordinates 0-1.0,0-1.0,0-1.0
        worldM = proj3d.world_transformation(xmin, xmax,
                                             ymin, ymax,
                                             zmin, zmax)

        # adjust the aspect ratio                          ##
        aspectM = proj3d.world_transformation(-kx + 1, kx, ##
                                              -ky + 1, ky, ##
                                              -kz + 1, kz) ##

        # look into the middle of the new coordinates
        R = np.array([0.5, 0.5, 0.5])

        xp = R[0] + np.cos(razim) * np.cos(relev) * self.dist
        yp = R[1] + np.sin(razim) * np.cos(relev) * self.dist
        zp = R[2] + np.sin(relev) * self.dist
        E = np.array((xp, yp, zp))

        self.eye = E
        self.vvec = R - E
        self.vvec = self.vvec / proj3d.mod(self.vvec)

        if abs(relev) > np.pi/2:
            # upside down
            V = np.array((0, 0, -1))
        else:
            V = np.array((0, 0, 1))
        zfront, zback = -self.dist, self.dist

        viewM = proj3d.view_transformation(E, R, V)
        perspM = proj3d.persp_transformation(zfront, zback)
        M0 = np.dot(viewM, np.dot(aspectM, worldM)) ##
        M = np.dot(perspM, M0)
        return M
    return get_proj



# class Singleplot(Subplot):
#     """Plot a simple diagram.

#     Usage
#     =====
#      >>> p = Singleplot(10, 7.5)
#      >>> p.addDataXY(tvec, vals)
#      >>> p.addDataXY(tvec, vals)
#      >>> p.plotprogress()
#     """

#     def __init__(self, width=10, height=7.5, dpi=80,
#                  bgcolor='white', projection=u'rectilinear'):
#         """
#         @param width: in cm
#         @param height: in cm
#         @param dpi: optional, dots-per-inch.
#         """
#         self.gridRowsCols = (1,1)
#         self.fig = plt.figure(
#             1, figsize=(width/2.54, height/2.54), dpi=dpi,
#             facecolor='white',
#             # edgecolor=None,
#             # linewidth=0.0,
#             # frameon=None,
#             # subplotpars=None,
#             # tight_layout=True,
#             )

#         Subplot.__init__(
#             self, name="singleplot",
#             mpInst=self, mpGridPos=(1,1), mpGridSpan=(0,0),
#             bgcolor=bgcolor, projection=projection)

#         return

#     def show(self, block=True):
#         """show figure (matplotlib.Figure) in interactive mode.
#         """
#         plt.show(block=block)


#     def plot(self, outFilename=None):
#         """updates the self.fig (matplotlib.Figure) showing the plot.

#         @param outFilename: The file name for the resulting png file.

#            Note: An optionally supplied self.getOutFilename function
#            takes precedence. See class description.
#         """
#         for spName in self.subplotNames:
#             sp = self.subplots[spName]
#             sp.setup()
#             sp.plot()

#         if hasattr(self, "getOutFilename"):
#             outFilename = self.getOutFilename(None)

#         if outFilename is not None:
#             self.fig.savefig(
#                 outFilename,
#                 transparent=False
#                 )
#         return

#     def plotprogress(self, idx=None, show_marker=True, show_rest=True,
#                      verbose=True, outFileNameTemplate=None):
#         """updates the MainPlot.fig (matplotlib.Figure) showing all
#         subplots and their data, if any, upto a certain idx.

#         @param idx: the range of indeces to plot (and save) the data.
#             idx can either be an integer or a list of integers. if saved, this
#             produces as many pictures as the length of the list is.
#             Pass None to get an image for each step.
#         @param show_marker: boolean, whether or not to show a marker at the
#             position of the current index.
#         @param show_rest: boolean, whether or not to show the rest of the
#             curve(s) from the position of the current index up to the end of
#             the datavalues.
#         @param verbose: boolean, for diagnostic output using L{bae.log_01.msg}.
#         @param outFileNameTemplate: How shall the file name for the resulting
#            png file be generated? Can either be a template string with one
#            integer placeholder (%d). Or a function to be called with one integer
#            file index as argument to deliver the complete file name.

#            Note: An optionally supplied self.getOutFilename function
#            takes precedence. See class description.
#         """
#         if isinstance(idx, int):
#             idx = [ idx, ]
#         elif idx is None:
#             idx = range(len(self.xValues))

#         if hasattr(self, "getOutFilename"):
#             getOutFilename = self.getOutFilename(None)
#         elif outFileNameTemplate is None:
#             getOutFilename = lambda i: None
#         elif isinstance(outFileNameTemplate, basestring):
#             getOutFilename = lambda i: (outFileNameTemplate % i)
#         elif callable(outFileNameTemplate):
#             getOutFilename = outFileNameTemplate
#         else:
#             raise ValueError(
#                 "ERROR: outFileNameTemplate argument of type %s not"
#                 " implemented." % type(outFileNameTemplate))

#         for i in idx:
#             if verbose: msg("  - plotprogress idx %d..."%i)
#             self.setup()
#             self.plotsingle(i, show_marker=show_marker, show_rest=show_rest)
#             if verbose: msg("    plotted")

#             outFilename = getOutFilename(i)
#             if outFilename is not None:
#                 self.fig.savefig(
#                     outFilename,
#                     transparent=False
#                     )
#                 if verbose: msg("    saved to '%s'"%outFilename)
#         return

if __name__=='main':
    msg('No syntax errors')
    exampleDefaultPlots()
