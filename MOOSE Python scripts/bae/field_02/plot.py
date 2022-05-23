"""WARNING: Interface still in alpha stage and subject to change without notice.

Collection of common plots often used for scalar-, vector- and tensor fields

If not opted-out a Scalar-, Vector- or TensorField will have an plot-object
which has a reference to Field-Values, Frames and Points and provides some
fieldtype-specific plotmethods::
    Usage:
        >>> # plot scalar values on some points over time/frame:
        >>> myScalarField[1:4,:].plot.plotOverFrames(currentIndex=5)
        >>> # quiver-plot for different frames to specified axes:
        >>> import matplotlib.pyplot as plt
        >>> fig,axhs = plt.subfigures(2,2)
        >>> for ii, axh in enumerate(axhs):
        >>>     axh.set_title(myVectorField.frames.names[ii][0])
        >>>     myVectorField[:,ii].plot.quiver(axh=axh)
"""

import sys

import numpy as np

from numpy_geom import vector2plunge_bearing

from bae.log_01 import msg

### check if grafic interface is available
import platform
if 'nuku' in platform.node().lower():
    import matplotlib #as mpl
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
else:
    try:
        import matplotlib as mpl
        import matplotlib.pyplot as plt
    except RuntimeError:
        # funny enough: if import failed once, you can not reset
        # matplotlib.use('Agg') afterwards
        msg('It seems that you are running this script on a server, ' +
            'because no graphic interface is available. ' +
            'So matplotlib-backend has to be swiched to "Agg". ' +
            'See the "nuku"-solution above.')


#check module is loaded by epydoc-build
#--> epyDocBuild=True to disable decorators
if 'epydoc' in sys.modules:
    _epyDocBuild = True
else:
    _epyDocBuild = False


########### HELPERFUNCTIONS #######################################################
#{ decorators and helper functions
def _checkDataConnection(func, passThrough=True):
    '''Decorator that checks if data.topo.coords and data.frames.ids
    is available
    '''
    if _epyDocBuild:
        return func
    else:
        def func_wrapper(self, *args, **kwargs):
            try:
                hasPoints = np.any(self.data.topo.coords)

            except:
                hasPoints = False

            try:
                hasFrames = np.any(self.data.frames.ids)
            except:
                hasFrames = False

            if hasPoints and hasFrames:
                return func(self, *args, **kwargs)

            else:
                if not hasPoints:
                    msg('WARNING! data has no valid points assigned!')
                if not hasFrames:
                    msg('WARNING! data has no valid frames assigned!')

                if passThrough:
                    return func(self, *args, **kwargs)
                else:
                    msg('Function %s skipped' % func.__name__)

        return func_wrapper


def _allowFieldAsInput(func):
    """ This decorator allows to pass ScalarFields or VectorFiedls instead of
    values + points + frames
    """
    if _epyDocBuild:
        return func
    else:
        def func_wrapper(*args, **kwargs):
            field = args[0]
            if type(field).__name__ in ('ScalarField', 'VectorField'):
                args = (field.values, field.topo, field.frames) + args[1:]
            return func(*args, **kwargs)
        return func_wrapper


def _initPlotHandles(**kwargs):
    ''' Creates or passes through an axes-object depending on presence of
    'axh' in kwargs. If it is specified it will be removed from kwargs.

    @kwarg axh: axes-handle (default = None)

    @kwarg projection: if axh is None, the created axes will have this
        projection-type

    @returns: axh, parent-figure-handle, kwargs

    @note: if axh and/or projection is specified, each will be removed from
        kwargs-dict
    '''
    # from mpl_toolkits.mplot3d import Axes3D
    # from bae.plot_01.mpl_core import southernHemisphere
    # get existing axes handle
    # if existent: this will remove axh from kwargs
    axh = kwargs.pop('axh', None)

    projection = kwargs.pop('projection', None)

    if axh is None:
        # create new fig and axes
        fig = plt.figure()
        axh = fig.add_subplot(1,1,1, projection=projection)

    else:
        # get fig/canvas from axh
        fig = axh.figure
        if projection == '3d':
            if not hasattr(axh, 'get_zlim'):
                msg('Warning!!! The axes you provided does not support ' +
                    'a %s projection.' % projection)

                axh = fig.add_axes(axh.get_position(),
                                   projection=projection)
        elif projection == 'shemisphere':
            if not axh.name == 'shemisphere':
                msg('Warning!!! The axes you provided does not support ' +
                    'a %s projection.' % projection)

                axh = fig.add_axes(axh.get_position(),
                                   projection=projection)

    return fig, axh, kwargs


def _uniqueMarkers(n):
    '''Returns n unique polygone-markers
    @param n: number of unique markers to be created
    '''
    markers = [(2 + ii, 1 + ii%2, ii/n * 90.0) for ii in range(1, n+1)]
    return markers

#} # End decorators and helper functions



###############################################################################
### ScalarPlots
#{Plotting scalar fields
@_allowFieldAsInput
def hist(values, points, frames, sortByTime = True, **kwargs):
    ''' Creates a HistogramPlot of values.
    @kwarg field: ScalarField (!instead of values, points and frames)

    @param values: numpyArray(fieldValues) of size Npoints x Nframes

    @param points: Points-Object holding pointsData

    @param frames: Frames-Object holding framesData

    @kwarg sortByTime: if True (default) each frame will have its own hist

    @kwarg axh: axes(handle) to plot on

    @returns: figureHandle, axesHandle, plotHandles

    @note: You can pass a ScalarField OR values, points and frames
    '''
    fig, axh, kwargs = _initPlotHandles(**kwargs)

    handles = []

    if sortByTime:
        for frame in xrange(values.shape[1]):
            _, _, pltH = axh.hist(values[:,frame].ravel(),
                                  **kwargs)
            handles.append(pltH)
    else:
        _, _, pltH = axh.hist(values.ravel(), **kwargs)
        handles.append(pltH)

    return fig, axh, handles


@_allowFieldAsInput
def scatter3D(values, points, frames, sortByTime = True, **kwargs):
    ''' Creates a 3D scatterPlot of values.

    @kwarg field: ScalarField (!instead of values, points and frames)

    @param values: numpyArray(fieldValues) of size Npoints x Nframes

    @param points: Points-Object holding pointsData

    @param frames: Frames-Object holding framesData

    @kwarg sortByTime: if True (default) each frame will have its own plot

    @kwarg axh: axes(handle) to plot on, needs to be projection='3d'

    @returns: figureHandle, axesHandle, plotHandles

    @note: You can pass a ScalarField OR values, points and frames

    @note: if you specified a none-'3d'-axes this will create a new axes
        overlaying the specified one
    '''
    fig, axh, kwargs = _initPlotHandles(projection='3d', **kwargs)

    handles = []
    if sortByTime:
        x,y,z = (points.coords[:,ii] for ii in range(3))
        markers = _uniqueMarkers(values.shape[1])

        for frame in xrange(values.shape[1]):

            try:
                label = 'frame: %s' % str(frames.ids[frame])
            except:
                label = 'frame: %d' % frame

            pltH = axh.scatter(x, y, z,
                               c=values[:,frame].ravel(),
                               marker=markers[frame],
                               label=label,
                               **kwargs)
            handles.append(pltH)

    else:
        coords = np.tile(points.coords, (values.shape[1],1))
        x,y,z = (coords[:,ii] for ii in range(3))
        pltH = axh.scatter(x, y, z, c=values.ravel(), **kwargs)

        handles.append(pltH)

    return fig, axh, handles


@_allowFieldAsInput
def plotOverFrames(values, points, frames, currentIndex=None, **kwargs):
    '''Creates a plot of values over time/frame

    @kwarg field: ScalarField (!instead of values, points and frames)

    @param values: numpyArray(fieldValues) of size Npoints x Nframes

    @param points: Points-Object holding pointsData

    @param frames: Frames-Object holding framesData

    @kwarg currentIndex: specifies the frameIndex up to which the plot will
        have a solid linestyle

    @kwarg axh: axes(handle) to plot on

    @kwarg seqTicks: if True it will try to set frames.seq as x-Ticks

    @kwarg seqTicksStartIndex: index of sequence-ticks to start with
        (only if seqTicks is set True, default = 0)

    @kwarg seqTicksDelta: number frames to skip for sequence-ticks
        (only if seqTicks is set True, default = 4)

    @kwarg color: can be
        - a single color for all ploted points given as a valid colorstring
        ('r', 'b', 'g', ...) or a (4 value) RGBA list
        - a list/np.array of point colors as RGBA-values (shape: nPts x 4)

        The defaultValue is None, so matplotlib will cycle default colors

    @returns: figureHandle, axesHandle, plotHandles

    @note: You can pass a ScalarField OR values, points and frames
    '''
    fig, axh, kwargs = _initPlotHandles(**kwargs)

    nFrames = values.shape[1]
    nPoints = values.shape[0]
    iFrames = np.arange(0,nFrames)

    if currentIndex is not None:
        currentIndex = frames._parseItem(currentIndex)
        if isinstance(currentIndex, slice):
            currentIndex = currentIndex.start

    else:
        currentIndex = nFrames

    seqTicks = kwargs.pop('seqTicks', False)
    seqTicksStartIndex = kwargs.pop('seqTicksStartIndex', 0)
    seqTicksDelta = kwargs.pop('seqTicksDelta', 4)

    handles = []
    cMarkerSize = kwargs.pop('cMarkerSize', 15)
    markers = _uniqueMarkers(nPoints)

    colors = kwargs.pop('color', None)
    try:
        if colors is None:
            colors = nPoints * [None,]
        elif isinstance(colors, basestring):
            colors = nPoints * [colors,]

        # convert into numpyArray if not already done
        colors = np.asarray(colors)

        if colors.shape in [(4,),(4,1)]:
            colors = nPoints * [colors,]
        elif (not colors.shape == (nPoints, 4) or   # RGBA-list
              not (colors.shape == (nPoints,) and
                   colors.dtype.kind in ['S', 'U'])): # string list
            pass
        else:
            raise ValueError

    except ValueError:
        print 'colors: ', colors
        raise ValueError('Pass a single color (valid colorString or' +
                         ' 4-Values-RGBA list)' +
                         ' Or pass a RGBA field of length of Points')

    for ptIdx in xrange(nPoints):
        phA = axh.plot(iFrames[:currentIndex+1],
                       values[ptIdx,:currentIndex+1],
                       color = colors[ptIdx],
                       **kwargs)

        if not currentIndex == nFrames:
            col, width = phA[0].get_color(), phA[0].get_linewidth()

            marker = kwargs.get('currentMarker', markers[ptIdx])

            phB = axh.scatter([currentIndex,],
                              [values[ptIdx,currentIndex],],
                              s=cMarkerSize, c=col, marker=marker)

            phC = axh.plot(iFrames[currentIndex:],
                           values[ptIdx,currentIndex:],
                           color=col, linewidth=width, linestyle=':')

        else:
            phB, phC = None, None

        handles.append({'current':phB,'past':phA,'future':phC})

    if seqTicks:
        try:
            seq = np.asarray(frames.seq.labels)
            ticks = np.arange(seqTicksStartIndex, seq.shape[0], seqTicksDelta)
            labels = seq[ticks]
            axh.set_xticks(ticks)
            axh.set_xticklabels(labels)
        except AttributeError:
            msg('WARNING! No sequence data stored in frames object yet. ' +
                'Use e.g. frames.regularSequence to create.')

    return fig, axh, handles


@_allowFieldAsInput
def plotAlongAxis(values, points, frames, axis='z', sortByTime=True,
                  flipAxis=False, **kwargs):
    ''' Creates a plot of values along specified axis

    @kwarg field: ScalarField (!instead of values, points and frames)

    @param values: numpyArray(fieldValues) of size Npoints x Nframes

    @param points: Points-Object holding pointsData

    @param frames: Frames-Object holding framesData

    @param sortByTime: if True (default) each frame will have its own hist

    @kwarg axis: axis to plot along (default : 'z')

    @kwarg flipAxis: if True (default : False) x and y - axis get flipped

    @kwarg axh: axes(handle) to plot on

    @returns: figureHandle, axesHandle, plotHandles

    @note: You can pass a ScalarField OR values, points and frames
    '''
    fig, axh, kwargs = _initPlotHandles(**kwargs)

    handles = []
    dirIdx = {'x':0, 'y':1, 'z':2}[axis.lower()]
    ii = np.argsort(points.coords[:,dirIdx])
    x = points.coords[ii,dirIdx]

    if sortByTime:
        for frame in xrange(values.shape[1]):
            if flipAxis:
                pltH = axh.plot(values[ii,frame].ravel(), x, **kwargs)
            else:
                pltH = axh.plot(x, values[ii,frame].ravel(), **kwargs)
            handles.append(pltH)
    else:
        if flipAxis:
            pltH = axh.plot(values[ii].ravel(), x, **kwargs)
        else:
            pltH = axh.plot(x, values[ii].ravel(), **kwargs)
        handles.append(pltH)

    return fig, axh, handles


class CommonScalarPlots(object):
    """ Provides a Collection of often used plots for ScalarFields

    Usage:
        >>> from bae.fields_00 import ScalarField
        >>> import matplotlib.pyplot as plt
        >>> sField = ScalarField('S', np.random.randn(128,4),
        ...                       points = np.random.randn(128,3),
        ...                       frames = [(2,d) for d in range(4)])
        >>> fig, axhs = plt.subplots(2,2)
        >>> for ii, axes in enumerate(axhs.flatten()):
        >>>     title = 'Frame (%d, %d)' % vecField.frames[ii][0]
        >>>     axes.set_title(title)
        >>>     vecField[5:10].plot.plotOverFrames(axh=axes, currentIndex=ii)
        >>> plt.show()

    @ivar data: reference of ScalarField
    """

    def __init__(self, field):
        self.data = field


    @_checkDataConnection
    def hist(self, sortByTime=True, **kwargs):
        """ Creates a HistogramPlot of values.

        @kwarg sortByTime: if True (default) each frame will have its own hist

        @kwarg axh: axes(handle) to plot on

        @returns: figureHandle, axesHandle, plotHandles
        """
        handles = hist(self.data, sortByTime=sortByTime, **kwargs)
        return handles


    @_checkDataConnection
    def scatter3D(self, sortByTime=True, **kwargs):
        ''' Creates a 3D scatterPlot of values.

        @kwarg sortByTime: if True (default) each frame will have its own plot

        @kwarg axh: axes(handle) to plot on, needs to be projection='3d'

        @returns: figureHandle, axesHandle, plotHandles

        @note: if you specified a none-'3d'-axes this will create a new axes
            overlaying the specified one
        '''
        handles = scatter3D(self.data, sortByTime=sortByTime, **kwargs)

        return handles


    @_checkDataConnection
    def plotOverFrames(self, currentIndex=None, **kwargs):
        ''' Creates a plot of values over time/frame

        @kwarg currentIndex: specifies the frameIndex up to which the plot will
            have a solid linestyle

        @kwarg axh: axes(handle) to plot on

        @returns: figureHandle, axesHandle, plotHandles
        '''
        handles = plotOverFrames(self.data,
                                 currentIndex=currentIndex, **kwargs)

        return handles

    def plotAlongAxis(self, axis='z', sortByTime=True, **kwargs):
        ''' Creates a plot of values along specified axis

        @param sortByTime: if True (default) each frame will have its own hist

        @kwarg axis: axis to plot along (default : 'z')

        @kwarg axh: axes(handle) to plot on

        @returns: figureHandle, axesHandle, plotHandles
        '''
        handles = plotAlongAxis(self.data, axis='z', sortByTime=True, **kwargs)

        return handles
#} # End Plotting scalar fields 


###############################################################################
### VectorPlots
#{Plotting vector fields
@_allowFieldAsInput
def quiver3D(values, points, frames, length=.1, **kwargs):
    '''Creates a 3D quiverPlot of VectorData

    @kwarg field: VectorField (!instead of values, points and frames)

    @param values: numpyArray(fieldValues) of size Npoints x Nframes x 3

    @param points: Points-Object holding pointsData

    @param frames: Frames-Object holding framesData

    @kwarg length: relative length of arrows (default .1)

    @kwarg axh: axes(handle) to plot on, needs to be a projection='3d'

    @returns: figureHandle, axesHandle, plotHandles

    @note: You can pass a VectorField OR values, points and frames

    @note: see U{https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.quiver.html}
        for further keywordarguments

    @note: if you specified a none-'3d'-axes this will create a new axes
            overlaying the specified one
    '''
    fig, axh, kwargs = _initPlotHandles(projection='3d', **kwargs)

    coords = points.coords
    handles = []

    for ii in xrange(values.shape[1]):
        h = axh.quiver(coords[:,0], coords[:,1], coords[:,2],
                       values[:,ii,0], values[:,ii,1], values[:,ii,2],
                       length=length, **kwargs)
        handles.append(h)

    return fig, axh, handles


def plungeBearingScatter(*args, **kwargs):
    """
    Plots a scatter distribution of irregular plunge-bearing-pairs
    to a southernHemisphere-projection. You can pass::
        a) a VectorField
        b) a Npoints x NFrames x 3 numpyArray or
        c) plunge and bearing as a numpyArray

    @kwarg field: a VectorField

    @kwarg values: Npoints x NFrames x 3 numpyArray

    @kwarg plg: plunge as numpyarray in radians

    @kwarg brg: bearing as numpyarray in radians

    @kwarg axh: axes(handle) to plot on, needs to be a projection='shemisphere'

    @returns: figureHandle, axisHandle, plotHandles

    @note: if you specified a none-'shemisphere'-axes this will create a new
        axes overlaying the specified one
    """


    errTxt = ('You can pass: a) a VectorField, ' +
              'b) a Npoints x NFrames x 3 numpyArray or ' +
              'c) plunge and bearing as a numpyArray')

    if len(args) == 1:
        if type(args[0]).__name__ == 'VectorField':
            values = args[0].values

        elif isinstance(args[0], np.ndarray) and args[0].shape[2] == 3:
            values = args[0]

        else:
            raise ValueError(errTxt)

        plg, brg = vector2plunge_bearing(values[:,:,0], values[:,:,1],
                                         values[:,:,2], asRad=True)

    elif (len(args) == 2
          and isinstance(args[0], np.ndarray)
          and isinstance(args[0], np.ndarray)):
        plg, brg = args[0], args[1]

    else:
        raise ValueError(errTxt)

    fig, axh, kwargs = _initPlotHandles(projection='shemisphere', **kwargs)

    handles = []

    # get scatter size/color from kwargs or set to default
    size = kwargs.pop('s',.5)
    color = kwargs.pop('c',None)

    # loop over frames
    for ii in xrange(plg.shape[1]):
        # set current scatterSize
        if isinstance(size, np.ndarray):
            s = size[:,ii,:]
        else:
            s = size

        # set current scatterColor
        if isinstance(color, np.ndarray):
            c = color[:,ii,:]
        else:
            c = color

        # plot plunge and bearing
        h = axh.scatter(brg[:,ii], plg[:,ii], s=s, c=c, **kwargs)
        handles.append(h)

    # set nice axes labels
    axh.adjustAxes()

    return fig, axh, handles


def plungeBearingDensity(*args, **kwargs):
    """
    Plots an estimated distribution of irregular plunge-bearing-pairs
    to a southernHemisphere-projection.
    This function requires the (great) sklearn-package, as this is the only
    implementation of a KernelDensityEstimation (see also scipy
    stats.gaussian_kde ) that allows to use a noneEuklidian metric.
    See a full example here:
    U{https://scikit-learn.org/stable/auto_examples/neighbors/plot_species_kde.html}

    You can pass::
        a) a VectorField
        b) a Npoints x NFrames x 3 numpyArray or
        c) plunge and bearing as a numpyArray

    @kwarg field: a VectorField

    @kwarg values: Npoints x NFrames x 3 numpyArray

    @kwarg plg: plunge as numpyarray in radians

    @kwarg brg: bearing as numpyarray in radians

    @kwarg colorMap: matplotlib-colormap (default plt.cm.Reds)

    @kwarg denseBrg: number of interpolationPoints for bearing (default 64)

    @kwarg densePlg: number of interpolationPoints for plunge (default 16)

    @kwarg bandWidth: bandwidth-parameter for KernelDensityEstimation

    @kwarg densityType: lin (default) or exp denstity

    @returns: figureHandle, axisHandle, contourfHandle

    @note: if you specified a none-'shemisphere'-axes this will create a new
        axes overlaying the specified one
    """
    # from transform import line, vector2plunge_bearing

    # from fields import VectorField

    # check for sklearn
    try:
        from sklearn.neighbors import KernelDensity
    except ImportError as e:
        raise ImportError('You need to install sklearn-package to use ' +
                          'a proper (angular) metric for KernelDensity' +
                          'Estimation. Use: ' +
                          '"pip install -U scikit-learn" or ' +
                          '"conda install scikit-learn"' +
                          '>>> %s' % str(e))


    ### process inputVariables
    errTxt = ('You can pass: a) a VectorField, ' +
              'b) a Npoints x NFrames x 3 numpyArray or ' +
              'c) plunge and bearing as a numpyArray')

    if len(args) == 1:
        if type(args[0]).__name__ == 'VectorField':
            values = args[0].values

        elif isinstance(args[0], np.ndarray) and args[0].shape[2] == 3:
            values = args[0]

        else:
            raise ValueError(errTxt)

        plg, brg = vector2plunge_bearing(values[:,:,0], values[:,:,1],
                                         values[:,:,2], asRad=True)

    elif (len(args) == 2
          and isinstance(args[0], np.ndarray)
          and isinstance(args[0], np.ndarray)):
        plg, brg = args[0], args[1]

    else:
        raise ValueError(errTxt)

    # set/create current axesHandle and set shemisphere-projection
    fig, axh, kwargs = _initPlotHandles(projection='shemisphere', **kwargs)

    # if not already done, flatten fields
    plg = plg.ravel()
    brg = brg.ravel()

    # use only not-nan-values for density-estimation
    ii = ~np.isnan(plg*brg)

    # get mesh-settings or set to default
    denseBrg = kwargs.pop('denseBrg', 64)
    densePlg = kwargs.pop('densePlg', 16)
    # get bandwidth-parameter or set to default
    bandWidth= kwargs.pop('bandWidth', .15)
    # get densityType or set default
    densityType = kwargs.pop('densityType', 'lin')

    # get color-map or set to default
    colorMap = kwargs.pop('cmap', plt.cm.Reds)

    ### estimate density
    # create regular plunge/bearing-grid
    P,B = np.meshgrid(np.linspace(0, .5*np.pi, densePlg),
                      np.linspace(0, 2*np.pi, denseBrg))
    pb = np.vstack([P.ravel(), B.ravel()]).T

    # estimate kernelDensity - use haversine-metric for correct angular
    # weightening
    kde = KernelDensity(bandwidth=bandWidth,
                        metric='haversine',
                        kernel='gaussian',
                        algorithm='ball_tree')

    kde.fit(np.vstack( (plg[ii], brg[ii]) ).T)

    D = kde.score_samples(pb)
    if densityType == 'exp':
        D = np.exp(D)

    D = D.reshape(P.shape)

    # plot density
    handles = axh.contourf(B, P, D, cmap=colorMap)

    # set nice axeslabels
    axh.adjustAxes()

    return fig, axh, handles


class CommonVectorPlots(object):
    """ Provides a Collection of often used plots for VectorFields.

    Usage:
        >>> from bae.fields_00 import VectorField
        >>> from bae.plot_01.mpl_core import southernHemisphere
        >>> import matplotlib.pyplot as plt
        >>> vecField = VectorField('V', np.random.randn(128,4,3),
        ...                       points = np.random.randn(128,3),
        ...                       frames = [(2,d) for d in range(4)])
        >>> fig, axhs = plt.subplots(2,2,
        ...                          subplot_kw=dict(projection='shemisphere'))
        >>> for ii, axes in enumerate(axhs.flatten()):
        >>>     title = 'Frame (%d, %d)' % vecField.frames[ii][0]
        >>>     axes.set_title(title)
        >>>     vecField[:,ii].plot.plungeBearingScatter(axh=axes)
        >>> plt.show()

    @ivar data: reference of VectorField
    """
    def __init__(self, field):
        self.data = field


    @_checkDataConnection
    def quiver3D(self, length=0.1, **kwargs):
        """Creates a 3D quiverPlot of VectorData

        @kwarg length: relative length of arrows (default .1)

        @kwarg axh: axes(handle) to plot on

        @returns: figureHandle, axesHandle, plotHandles

        @note: see U{https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.quiver.html}
            for further keyword arguments
        """
        handles = quiver3D(self.data, length=length, **kwargs)

        return handles


    @_checkDataConnection
    def plungeBearingScatter(self, absAsMeta='', **kwargs):
        """
        Transforms a VectorField to plunge and bearing and creates a
        scatterplot to a southernHemisphere-projection

        @param absAsMeta: sets abs(vector) as metaValue to color and/or
            size of scatterplot (default:''; opts: 'color', 'size',
            'sizeAndColor'))

        @kwarg size: scale parameter for scatter-point size (default: 10)

        @kwarg color: specifies color of scatter-points, if 'color' is set_title
            in absAsMeta this will be overriden

        @kwarg axh: axes(handle) to plot on

        @returns: figureHandle, axisHandle, contourfHandle

        @note: if you specified a none-'shemisphere'-axes this will create a
            new axes overlaying the specified one
        """
        # points = self.data.topo
        # frames = self.data.frames
        values = self.data.values

        # get scatter size/color from kwargs or set to default
        size = kwargs.pop('size', 10)
        color = kwargs.pop('color', None)

        # normalizing vector
        norm = np.linalg.norm(values, axis=2)
        ii = norm > 1E-12
        values[~ii,:] = 0
        values[ii,:] = values[ii,:] / norm[ii,np.newaxis]

        # force pointing downwards
        sign = np.sign(values[:,:,2])[:,:,np.newaxis]
        sign[sign==0] = 1
        values *= -sign

        # transform to plunge and bearing
        plg, brg = vector2plunge_bearing(values[:,:,0], values[:,:,1],
                                         values[:,:,2], asRad=True)

        # if requested apply vector-norm to scatter- size and/or color
        if absAsMeta.lower().find('size')>=0:
            size = size * norm / norm.max()
        if absAsMeta.lower().find('color')>=0:
            color = norm

        handles = plungeBearingScatter(plg, brg, s=size, c=color, **kwargs)

        return handles


    def plungeBearingDensity(self, **kwargs):
        """
        Plots an estimated distribution of irregular plunge-bearing-pairs
        to a southernHemisphere-projection.

        @kwarg colorMap: matplotlib-colormap (default plt.cm.Reds)

        @kwarg denseBrg: number of interpolationPoints for bearing (default 64)

        @kwarg densePlg: number of interpolationPoints for plunge (default 16)

        @kwarg bandWidth: bandwidth-parameter for KernelDensityEstimation

        @kwarg axh: axes(handle) to plot on

        @returns: figureHandle, axisHandle, contourfHandle

        @note: all frames/points stored in VectorField will be taken in account
            at once

        @note: if you specified a none-'shemisphere'-axes this will create a
            new axes overlaying the specified one

        @warning: This funtion requires the sklearn-package.
        """

        handles = plungeBearingDensity(self.data, **kwargs)

        return handles
#} # End Plotting vector fields
