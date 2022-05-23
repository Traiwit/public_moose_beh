# -*- coding: utf-8 -*-
"""
Some default plotting proceedures for materials and yield surfaces

Example:
========
    >>> from bae.material_01 import PlasticLevkovitchReuschAniso
    >>> from bae.material.plot import plotAniso, plotYieldSurf
    >>> plr = PlasticLevkovitchReuschAniso()
    >>> plr.fromComputator(50E6, 55, anisoS=1.2, anisoN=1.08)
    >>> figAniso,_,_ = plotAniso(plr)
    >>> figYield,_,_ = plotYieldSurf(plr)
    >>> figAniso.show()
    >>> figYield.show()
"""
from bae.log_01 import msg

import numpy as np
try:
    import matplotlib

    try: # test if graphical backend is available
        import Tkinter
    except ImportError:
#        msg('It seems that no graphic interface is available. '
#            'So matplotlib-backend has to be swiched to "Agg" '
#            'and you can not use matplotlib interactively.')
        matplotlib.use('Agg')

    import matplotlib.pyplot as plt
    import pylab
    import matplotlib.gridspec as gridspec

except ImportError:
    plt, gridspec, pylab = None, None, None


try:
    from scipy.interpolate import InterpolatedUnivariateSpline
except ImportError:
    InterpolatedUnivariateSpline = None

from material import (PlasticLevkovitchReusch,
                      PlasticLevkovitchReuschAniso,
                      SinglePlaneMohrCoulomb,)
from models import restoreComputatorInputs
from yieldsurfaces import (getPAndQFromThetaAndS,
                           pQThetaToEigs,
                           QuasiMohrCoulomb,)

from matpropcollection import characteristics
from bae.colormap_01 import ColorBrewer




#%%
# defaultOrder is used for output to sort (potential) material properties
# in a reasonable order
defaultOrder = ['type', 'pst', 'density', 'E', 'nu', 'K', 'G',
                'UCSi','GSI', 'D', 'cohesion', 'friction',
                's', 'mb', 'alpha', 'e',]

# the following values will be ommited for e.g. excel output
defaultHiddenValues = ['umatDelete', 'umatDepvarNum']
defaultDocHiddenValues = (defaultHiddenValues +
                          ['dampingAlpha','critDamage','voidStiffnessRatio'])



plotParams = {
    ### note: sizes can be an integer, or one of ['x-large','large',
    ### 'medium','small',x-small',...]

    ### axes and axis formatting:
    'axes.titlesize': 'medium',     #default: 'large'
    'axes.labelsize': 'x-small',      #default: 'medium'
    #'axes.labelcolor': 'black',     #default: 'black'
    'xtick.labelsize': 'x-small',   #default: 'medium'
    'xtick.direction': 'in',        #default: 'out'
    'ytick.labelsize': 'x-small',    #default: 'medium'
    'ytick.direction': 'in',        #default: 'out'
    'legend.loc': 'upper left',     #default: 'best'
    'legend.fontsize': 'x-small',     #default: 'medium
    #'legend.title_fontsize': 'x-small',
}

#%% Helper and Service functions

#%% getting minima and maxma via interpolated splines
def quadraticSplineRoots(spl):
    roots = []
    knots = spl.get_knots()
    for a, b in zip(knots[:-1], knots[1:]):
        u, v, w = spl(a), spl((a+b)/2), spl(b)
        t = np.roots([u+w-2*v, w-u, 2*v])
        t = t[np.isreal(t) & (np.abs(t) <= 1)]
        roots.extend(t*(b-a)/2 + (b+a)/2)
    return np.array(roots)

def getApproxMinFromSpline(phis, samples):
    if not InterpolatedUnivariateSpline:
        raise ImportError('You need scipy to estimate minimum by splineFit')
    f = InterpolatedUnivariateSpline(phis, samples, k=4)
    # get local minimum
    crPts = quadraticSplineRoots(f.derivative())
    # also check the endpoints of the interval
    crPts = np.append(crPts, (phis[0], phis[-1]))
    crVals = f(crPts)
    minIdx = np.argmin(crVals)
    return crPts[minIdx], crVals[minIdx]


def getApproxMaxFromSpline(phis, samples):
    if not InterpolatedUnivariateSpline:
        raise ImportError('You need scipy to estimate minimum by splineFit')
    f = InterpolatedUnivariateSpline(phis, samples, k=4)
    # get local minimum
    crPts = quadraticSplineRoots(f.derivative())
    # also check the endpoints of the interval
    crPts = np.append(crPts, (phis[0], phis[-1]))
    crVals = f(crPts)
    maxIdx = np.argmax(crVals)
    return crPts[maxIdx], crVals[maxIdx]

#%% Some UTF-8-Letters for Copy and Paste #####################################
# # Greek
# # α β γ δ ε ζ η θ ι κ λ μ ν ξ ο π ρ ς σ τ υ φ χ ψ ω
# #
# # Super- and Subscripts
# # numbers and common mathematical symbols
# # ⁰ ¹ ² ³ ⁴ ⁵ ⁶ ⁷ ⁸ ⁹ ⁺ ⁻ ⁼ ⁽ ⁾ ₀ ₁ ₂ ₃ ₄ ₅ ₆ ₇ ₈ ₉ ₊ ₋ ₌ ₍ ₎
# #
# # a full superscript Latin lowercase alphabet except q
# # ᵃ ᵇ ᶜ ᵈ ᵉ ᶠ ᵍ ʰ ⁱ ʲ ᵏ ˡ ᵐ ⁿ ᵒ ᵖ ʳ ˢ ᵗ ᵘ ᵛ ʷ ˣ ʸ ᶻ
# #
# # a limited uppercase Latin alphabet  - no C, F, Q, S, X, Y, Z
# # ᴬ ᴮ ᴰ ᴱ ᴳ ᴴ ᴵ ᴶ ᴷ ᴸ ᴹ ᴺ ᴼ ᴾ ᴿ ᵀ ᵁ ⱽ ᵂ
# #
# # a few subscripted lowercase letters
# # ₐ ₑ ₕ ᵢ ₖ ₗ ₘ ₙ ₒ ₚ ᵣ ₛ ₜ ᵤ ᵥ ₓ
# #
# # and some Greek letters
# # ᵅ ᵝ ᵞ ᵟ ᵋ ᶿ ᶥ ᶲ ᵠ ᵡ ᵦ ᵧ ᵨ ᵩ ᵪ
# =============================================================================

#%% outputConversions
def outputConversions(name, value, latex=False, valueNone='------'):
    """"Service fuction that changes a name and its value according to the
    specified rules.
    """
    from bae.future_01 import utf8Subs

    if name in ['E', 'K', 'G']:
        name = (name + ' [GPa]', r'$%s$ [GPa]' % name)
        value /= 1E9

    elif name in ['UCSi',]:
        name = (name + ' [MPa]', r'$UCS_i$ [MPa]')
        value /= 1E6

    elif name in ['GSI',]:
        name = (name, r'$GSI$')

    elif name == 'pst':
        #\epsilon_{plast}
        name =( u'\u03B5' + utf8Subs('plast',1) + ' [%]',
                r'$\varepsilon_{plast}$ [$\%$]')
        value *= 100.

    elif name == 'nu':
        name = (u'\u03BD', r'$\nu$')

    elif name == 'density':
        name = (u'\u03C1 [kg/m³]', r'$\rho$ '+u'[kg/m³]')

    elif name == 'mb':
        name = (u'm_b', r'$m_b$')

    elif name == 'phi':
        name = (u'φ [°]', r'$\phi$ [$\degree$]')

    elif name == 'c':
        name = ('cohesion [kPa]', r'$c$ [kPa]')
        value /= 1000

    elif name == 'anisoN':
        name = (u'n' + utf8Subs('aniso',1), r'$n_{aniso}$')

    elif name == 'anisoS':
        name = (u's' + utf8Subs('aniso',1), r'$s_{aniso}$')

    elif name == 'dilation':
        name = ('dilation', r'$d$')

    else:
        name = (name, r'$%s$' % name)

    if value is None:
        value = valueNone

    if isinstance(name, tuple):
        if latex:
            name = name[1]
        else:
            name = name[0]
    return name, value


# naming of pstSamples
def pstSampleNames(mat):
    """Returns the verbal degradation-states of a material depending on the
    number of stored pst-samples. If length is:
        1 'PEAK'
        2 'PEAK', 'RES'
        3 'PEAK','TRANS','RES'
        more than 3 'PEAK','TRANS1', 'TRANS2',...,'RES'
    is returned.
    @param mat: can be a enumerable(e.g. pst-sample-list) or integer value
    """
    if isinstance(mat, int):
        count = mat
    else:
        count = len(mat)
    if count == 1:
        return ['PEAK',]
    elif count == 3:
        return ['PEAK', 'TRANS', 'RES']
    else:
        names = ['TRANS%d'%ii for ii in range(len(mat))]
        names[-1] = 'RES'
        names[0] = 'PEAK'
        return names


def _getSampleNamesAndIndices(mat, requstedSamples):
    sampleNames = pstSampleNames(mat)

    if isinstance(requstedSamples, basestring):
        if requstedSamples.upper() == 'ALL':
            pstSamples = sampleNames
        else:
            pstSamples = [requstedSamples,]
    else:
        pstSamples = requstedSamples

    pstSamples = [s.upper() for s in pstSamples]

    sampleIdxs = [ii for ii,sampleName in enumerate(sampleNames)
                   if sampleName.upper() in pstSamples]

    return pstSamples, sampleIdxs

def valueFmt(val):
    if isinstance(val, float) and val > 0 and np.log10(val) < -2:
        return '%.2e'
    elif isinstance(val, float):
        return '%.2f'
    else:
        return '%s'

def valueToFmt(val):
    return valueFmt(val) % val



#%% matplotlib functions
def plotYieldSurfCompressionMeridian(mat, S3Min=None, S3Max=None, axh=None,
                  thetaRange=False, pstSamples='all'):
    """Plots the yieldsurfaces of a L{material<bae.material_01.material.DefaultPst>}
        or a single L{yield surface<bae.material_01.yieldsurfacesExtendedMenetreyWillam>}
        object on compression meridian (M{S{theta}=S{pi}/3}).

    Example:
        >>> import matplotlib.pyplot as plt
        >>> from bae.material_01 import PlasticLevkovitchReuschAniso
        >>> from bae.material_01.plot import plotYieldSurfCompressionMeridian
        >>> plr = PlasticLevkovitchReuschAniso()
        >>> plr.fromComputator(50E6, 55, version='V14')
        >>> plotYieldSurfCompressionMeridian(plr, S3Max=10E6, thetaRange=True,
        ...                                  pstSamples=['PEAK', 'RES'])
        >>> plt.show()

    @param mat: a L{material<bae.material_01.material.DefaultPst>} or a single
        L{yield surface<bae.material_01.yieldsurfacesExtendedMenetreyWillam>}
        object
    @param S3Min: lower bound of SMinor to plot. If None, the tensile
        strength of (peak) yield surface is used.
    @param S3Max: upper bound of SMinor to plot. If None, -15 of tensile
        strength of (peak) yield surface is used.
    @param axh: axes-handle to plot yieldsurfaces to. If None, a new figure and
        axes will be created.
    @param thetaRange: if True, the area between compression meridian and
        projection of tensile meridian (M{S{theta}=S{pi}/6}) gets plotted.
    @param pstSamples: if 'all', all samples will be plotted. Use verbal names
        of pstSamples or a list of those to filter (see L{pstSampleNames} for
        naming convention).
    @returns figHandle, axesHandle, (plotHandle,...), (legendHandle,...)
    """
    nSamples = 128
    scale = 1E-6
    # parse inputs an prepare
    if not hasattr(mat, '__iter__'): # in case you only passed a sample
        mat = [mat,]

    sampleNames = pstSampleNames(mat)
    pstSamples, sampleIdxs = _getSampleNamesAndIndices(mat, pstSamples)

    if not sampleIdxs:
        msg('No pstSamples left to plot when using %s as filter'
            % str(pstSamples))
        return 4*(None,)

    if not axh:
        fig, axh = plt.subplots()
        # axh.set_xlabel(r'$S_{Minor}$ in MPa')
        # axh.set_ylabel(r'$S_{Major}$ in MPa')
        axh.set_xlabel(r'$S_2=S_3$ in MPa')
        axh.set_ylabel(r'$S_1$ in MPa')
    else:
        fig = axh.get_figure() # to return a reasonalbe value for fig

    lineStyles = ['--' for i in sampleIdxs]
    lineStyles[-1] = ':'
    lineStyles[0] = '-'

    if S3Max is None:
        tensStrength = mat[0].yieldSurf.getIsoTensileStrength()
        if np.isclose(tensStrength,0):
            S3Max = 1e6
        else:
            S3Max = -15*mat[0].yieldSurf.getIsoTensileStrength()


    plotHandles = []
    for sIdx in sampleIdxs:
        plotHandles.append([])
        try:
            ys = mat[sIdx].yieldSurf
        except AttributeError:
            ys = mat[sIdx]

        # stress at compression meridian: theta = pi/3
        S1c,_,S3c = ys.sampleMeridianEigs(S3Max=S3Max, S3Min=S3Min,
                                    distribution='alpha',theta=np.pi/3.,
                                    samples=nSamples, miningNotation=True)
        S1c *= scale
        S3c *= scale

        if thetaRange:
            # stress at tensile meridian: theta = pi/6
            S1t,_,S3t = ys.sampleMeridianEigs(S3Max=S3Max, S3Min=S3Min,
                                    distribution='alpha',theta=np.pi/6.,
                                    samples=nSamples, miningNotation=True)
            S1t *= scale
            S3t *= scale
            S3cIp = np.interp(S1t, S1c, S3c)
            phf1 = axh.fill_betweenx(S1t, S3cIp, S3t, alpha=0.5)

            phu = axh.plot(S3t, S1t, ls='--', c='gray')
            phl = axh.plot(S3c, S1c, ls='-.', c='gray')
            plotHandles[-1].extend([phu[0], phl[0], phf1])

        else:
            ph = axh.plot(S3c, S1c, ls=lineStyles[sIdx],
                          label=sampleNames[sIdx])
            plotHandles[-1].append(ph)

    if thetaRange:
        lh1 = axh.add_artist(axh.legend(plotHandles[0][:2],
                            #  [r'$\theta=30\degree$',r'$\theta=60\degree$',],
                            [r'$S_1 \geq S_2=S_3$',r'$S_1=S_2 \geq S_3$',],
                             loc='lower right'))
        lh2 = axh.legend([p[-1] for p in plotHandles],
                         [sampleNames[ii] for ii in sampleIdxs])
        legendHandles = [lh1, lh2]
    else:
        legendHandles = [axh.legend(),]


    return fig, axh, plotHandles, legendHandles


def plotYieldSurfpq(mat, pMin=None, pMax=None, S3Min=None, S3Max=None,
                    axh=None, thetaRange=False, pstSamples='all'):
    """Plots hydrostatic pressure vs. vonMises stress of the yieldsurfaces of
        a L{material<bae.material_01.material.DefaultPst>} or a single
        L{yield surface<bae.material_01.yieldsurfacesExtendedMenetreyWillam>}
        object.

    Example:
        >>> import matplotlib.pyplot as plt
        >>> from bae.material_01 import PlasticLevkovitchReuschAniso
        >>> from bae.material_01.plot import plotYieldSurfpq
        >>> plr = PlasticLevkovitchReuschAniso()
        >>> plr.fromComputator(50E6, 55, version='V14')
        >>> plotYieldSurfpq(plr, S3Max=10E6, thetaRange=True,
        ...                 pstSamples=['PEAK', 'RES'])
        >>> plt.show()

    @param mat: a L{material<bae.material_01.material.DefaultPst>} or a single
        L{yield surface<bae.material_01.yieldsurfacesExtendedMenetreyWillam>}
        object
    @param pMin: lower bound of p to plot.
    @param pMax: upper bound of p to plot.
    @param S3Min: lower bound of SMinor to plot. Only used if pMin is None.
    @param S3Max: upper bound of SMinor to plot. Only used if pMax is None.
    @param axh: axes-handle to plot yieldsurfaces to. If None, a new figure and
        axes will be created.
    @param thetaRange: if True, the area between compression meridian and
        projection of tensile meridian gets plotted.
    @param pstSamples: if 'all', all samples will be plotted. Use verbal names
        of pstSamples or a list of those to filter (see L{pstSampleNames} for
        naming convention).
    @note: if no upper/lower bound (neither p nor S3) is set the
        tensileStrength/-15*tensileStrength is used
    @returns figHandle, axesHandle, (plotHandle,...), (legendHandle,...)
    """
    nSamples = 128
    scale = 1E-6

    # parse inputs an prepare
    if not hasattr(mat, '__iter__'): # in case you only passed a sample
        mat = [mat,]

    sampleNames = pstSampleNames(mat)
    pstSamples, sampleIdxs = _getSampleNamesAndIndices(mat, pstSamples)

    if not sampleIdxs:
        msg('No pstSamples left to plot when using %s as filter'
            % str(pstSamples))
        return 4*(None,)

    if not axh:
        fig, axh = plt.subplots()
        axh.set_xlabel(r'$p$ in MPa')
        axh.set_ylabel(r'$\sigma_{Mises}$ in MPa')
    else:
        fig = axh.get_figure() # to return a reasonalbe value for fig

    try:
        ys = mat[0].yieldSurf
    except AttributeError:
        ys = mat[0]

    if pMin is None and S3Min is not None:
        pMin = getPAndQFromThetaAndS(ys, theta=np.pi/3, S3=[S3Min,],
                                     miningNotation=True)[0][0]
    if pMax is None and S3Max is not None:
        pMax = getPAndQFromThetaAndS(ys, theta=np.pi/3, S3=[S3Max,],
                                     miningNotation=True)[0][0]

    if pMin is None:
        pMin = ys.getIsoTensileStrength()
    if pMax is None:
        pMax = -15*ys.getIsoTensileStrength()

    if pMin >= pMax:
        raise ValueError('pMin is larger than pMax. Check input.')

    lineStyles = ['--' for i in sampleIdxs]
    lineStyles[-1] = ':'
    lineStyles[0] = '-'

    plotHandles = []
    for sIdx in sampleIdxs:
        p = np.linspace(pMin, pMax, nSamples)
        plotHandles.append([])
        try:
            ys = mat[sIdx].yieldSurf
        except AttributeError:
            ys = mat[sIdx]

        # stress at compression meridian: theta = pi/3
        thetac = np.full_like(p, np.pi/3.)
        thetat = np.full_like(p, np.pi/6.)
        qc = ys.ptOnYieldSurf(p=p, theta=thetac)
        qt = ys.ptOnYieldSurf(p=p, theta=thetat)
        p *= scale
        qc *= scale
        qt *= scale

        if thetaRange:
            phf1 = axh.fill_between(p, qc, qt, alpha=0.5)
            phu = axh.plot(p, qt, ls='--', c='gray')
            phl = axh.plot(p, qc, ls=':', c='gray')
            plotHandles[-1].extend([phu[0], phl[0], phf1])

        else:
            ph = axh.plot(p, qc, ls=lineStyles[sIdx],
                          label=sampleNames[sIdx])
            plotHandles[-1].append(ph)

    if thetaRange:
        lh1 = axh.add_artist(axh.legend(plotHandles[0][:2],
                            #  [r'$\theta=30\degree$',r'$\theta=60\degree$',],
                            [r'$S_1 \geq S_2=S_3$',r'$S_1=S_2 \geq S_3$',],
                             loc='lower right'))
        lh2 = axh.legend([ph[-1] for ph in plotHandles],
                         [sampleNames[ii] for ii in sampleIdxs])
        legendHandles = [lh1, lh2]
    else:
        legendHandles = [axh.legend(),]

    return fig, axh, plotHandles, legendHandles



def plotYieldSurf(mat, view='compressionMeridian', **kwargs):
    """Convenience fuction that calls:
        - L{plotYieldSurfCompressionMeridian} if view is 'comp'
        - L{plotYieldSurfpq} if view is 'pq'
    Other keyword-arguments will be passed to these functions.
    @param mat: a L{material<bae.material_01.material.DefaultPst>} or a single
        L{yield surface<bae.material_01.yieldsurfacesExtendedMenetreyWillam>}
        object
    @param view: sets the ploting space
    @returns figHandle, axesHandle, (plotHandle,...), (legendHandle,...)
    """
    if 'comp' in view.lower():
        return plotYieldSurfCompressionMeridian(mat, **kwargs)
    elif 'pq' in view.lower():
        return plotYieldSurfpq(mat, **kwargs)
    else:
        raise NotImplementedError


def plotAniso(mat, confinements=None, ucsScaled=False, maxScaled=True,
               axh1=None, axh2=None, overwriteAxesProps=[False,False]):
    """Plots the anisotropic behaviour of a L{material<bae.material_01.material.DefaultPst>}
    (only PEAK-level) or a single
    L{yield surface<bae.material_01.yieldsurfacesExtendedMenetreyWillam>}.
        - first plot: compressive strength vs. orientation (S{phi}) for
            different confinements
        - second plot: compressive strength vs. confinement (in compression
            meridian) for M{S{phi}=0} and for the weakest orientation
            M{S{phi}_w(confinement)}

    Example:
        >>> import matplotlib.pyplot as plt
        >>> from bae.material_01 import PlasticLevkovitchReuschAniso
        >>> from bae.material_01.plot import plotYieldSurfpq
        >>> plr = PlasticLevkovitchReuschAniso()
        >>> plr = PlasticLevkovitchReuschAniso()
        >>> plr.fromComputator(50E6, 55, anisoS=1.2, anisoN=1.08)
        >>> plotAniso(plr)
        >>> plt.show()

    @param mat: a L{material<bae.material_01.material.DefaultPst>} or a single
        L{yield surface<bae.material_01.yieldsurfacesExtendedMenetreyWillam>}
        object
    @param confinements: list of confinements (tensile=positive). If None
        [-2*ii*UCS for ii in range(6)] is used.
    @param ucsScaled: if True, the strength axes of first plot gets scaled
        by UCS
    @param axh1: axesHandle for first plot. If None, a new fig and axes will
        be created.
    @param axh2: axesHandle for second plot. If None, a new fig and axes will
        be created.
    @returns (figHandle,...), (axesHandle,...), (plotHandle,...), (legendHandle,...)
    """
    import matplotlib.colors as colors

    try:
        ys = mat[0].yieldSurf
    except AttributeError:
        ys = mat
    UCS = ys.getUCS()

    # createting confinementSteps if not provided
    if confinements is None:
        confinements = [-ii*UCS for ii in range(6)]
        confinements[0] *= -1 # avoids '-0' in legend

    # default-colors: use bluesColormap and remove the first (close to white)
    colors1 = ColorBrewer().Blues.asHexStr(len(confinements)+1)[1:][::-1]
    zeroColor = None # will be set if confinement is zero

    # preparing plot axes if not present
    if not axh1:
        fig1, axh1 = plt.subplots()
        overwriteAxesProps[0] = True
    else:
        fig1 = axh1.get_figure() # to return a reasonalbe value for fig

    if not axh2:
        fig2, axh2 = plt.subplots()
        overwriteAxesProps[1] = True
    else:
        fig2 = axh2.get_figure() # to return a reasonalbe value for fig

    if overwriteAxesProps[0]:
        axh1.grid(True)
        axh1.set_xlabel(r'$\phi=\measuredangle$(foliation-normal, load) in $\degree$ ')
        if maxScaled:
            # axh1.set_ylabel('$S_{Major}$ / $S_{Minor}(\phi=0)$ in %')
            axh1.set_ylabel('$S_{1}$ / $S_{3}(\phi=0)$ in %')
        elif ucsScaled:
            # axh1.set_ylabel('$S_{Major}$ / $UCS$ in %')
            axh1.set_ylabel('$S_{1}$ / $UCS$ in %')
        else:
            #axh1.set_ylabel('$S_{Major}$ in MPa')
            axh1.set_ylabel('$S_{1}$ in MPa')
        axh1.set_xlim((0,90))

    if overwriteAxesProps[1]:
        axh2.grid(True)
        # axh2.set_xlabel(r'$S_{Minor}$ in MPa')
        # axh2.set_ylabel(r'$S_{Major}$ in MPa')
        axh2.set_xlabel(r'$S_2=S_3$ in MPa')
        axh2.set_ylabel(r'$S_1$ in MPa')

    phis = np.linspace(-0, 90, 45+1)
    plotHandles1 = []
    for ii,confinement in enumerate(confinements):
        strengths = ys.getOrientedStrength(phis, S11=confinement)
        valid = ~np.isnan(strengths)
        if not valid.all():
            msg("Warning! Can't find strength values for all angles."
                " Mostlikely the newton-iteration did not converge.")
        if not valid.sum() >= 3:
            msg('Even worse: less than 3 valid strength values found for %f!'
                % confinement)
            continue

        if ucsScaled:
            scale = -100./float(UCS)
        elif maxScaled:
            scale = 100./strengths[0]
        else:
            scale = -1E-6

        phiMin,sMin = getApproxMaxFromSpline(phis[valid], strengths[valid])

        label = '%.1f' % (-1*confinement * 1E-6)
        p1 = axh1.plot(phis, scale*strengths,
                       color=colors1[ii], label=label)
        p2 = axh1.scatter([phiMin,], [scale*sMin],
                          color=colors1[ii])
        plotHandles1.extend([p1,p2])
        if confinement == 0:
            zeroColor = colors1[ii]
            if not ucsScaled and not maxScaled:
                # annotate unconfined plot
                axh1.scatter([0,], [-UCS*scale,], color=colors1[ii])
                axh1.text(0, -UCS*scale,
                         r'$UCS_i\sqrt{s}=%.1f$MPa'%(UCS*1E-6),
                         fontsize='x-small')

        leg1 = axh1.legend(title=r'confinement in MPa', ncol=2,
                          bbox_to_anchor=(0.01,0.01), loc='lower left')
        leg1.get_title().set_fontsize(leg1.prop._size)

    phis = np.linspace(0, 90, 45+1)

    S3s = np.linspace(np.min(confinements), np.max(confinements), 128)

    S1best = ys.getOrientedStrength(0, S11=S3s)

    phiWorst, S1worst = [],[]
    for S3 in S3s:
        ss = ys.getOrientedStrength(phis, S11=S3)
        valid = ~np.isnan(ss)
        if not valid.any():
            phiWorst.append(np.nan)
            S1worst.append(np.nan)
        phiMin,sMin = getApproxMaxFromSpline(phis[valid], ss[valid])
        phiWorst.append(phiMin)
        S1worst.append(sMin)

    S1worst = np.asarray(S1worst)
    S3s, S1best, S1worst = -S3s*1E-6, -S1best*1E-6, -S1worst*1E-6

    phb1 = axh2.plot(S3s, S1best, label=r'$\phi=0.0\degree$',c=zeroColor)

    phb2 = axh2.scatter(S3s, S1worst, c=phiWorst,s=1,
                        norm=colors.Normalize(vmin=45, vmax=90),
                        label=r'worst $\phi$')
    #phb2 = axh2.plot(S3s, S1worst, c='lightgrey')
    cbar = fig2.colorbar(phb2)
    cbar.set_label('$\phi$ of lowest strength')
    leg2 = axh2.legend()

    leg1.get_title().set_fontsize(leg1.prop._size)

    return (fig1,fig2), (axh1,axh2), (plotHandles1, (phb1,phb2)),  (leg1, leg2)


def _setTableMatplotlib(axes, title, **kwargs):
    #makes some default settings for tabulars
    titleArgs = dict(loc='left',
                     fontsize='small')
    titleArgs.update( kwargs.pop('title', {}) )
    axes.set_title(title, **titleArgs)

    tableArgs = dict(loc='best',
                     colLoc='center',
                     )
    tableArgs.update(kwargs)
    tab = axes.table(**tableArgs)
    tab.auto_set_font_size(False)
    tab.set_fontsize(6)
    return tab

def _setExcelFormats(workBook):
    ## cellFormats
    formats = {}
    formats['title'] = workBook.add_format({'bold':True, 'font_size':14})
    formats['subTitle'] = workBook.add_format({'bold':True, })

    formats['header'] = workBook.add_format({'bottom':2, 'bg_color':'#C0C0C0',
                                     'align': 'center',})

    formats['std'] = workBook.add_format({})
    formats['exp'] = workBook.add_format({'num_format': '0.00E+#0'})
    formats['flt'] = workBook.add_format({'num_format': '##0.00'})

    formats['stdULine'] = workBook.add_format({'bottom':1,})
    formats['expULine'] = workBook.add_format({'num_format': '0.00E+#0',
                                               'bottom':1,})
    formats['fltULine'] = workBook.add_format({'num_format': '##0.00',
                                               'bottom':1,})


    def valueFmt(val, underLine=False):
        if isinstance(val, float) and val > 0 and np.log10(val) < -2:
            ff = 'exp'
        elif isinstance(val, float):
            ff = 'flt'
        else:
            ff = 'std'

        if underLine:
            ff += 'ULine'
        return formats[ff]

    formats['valFun'] = valueFmt

    return formats

def _setExcelSheetTable(sheet, cellData, formats, row=0, col=0, sepLines=None,
                        rowLabels=None, columnLabels=None, title=None):

    valueFmt = formats['valFun']
    if title:
        sheet.write(row, col, title, formats['title'])
        row += 1

    if rowLabels:
        sCol = col+1
    else:
        sCol = col

    if not sepLines:
        sepLines = [False for _ in range(len(cellData))]

    _colA, _colB = sCol, sCol
    _rowA, _rowB = row, row

    if columnLabels:
        for ii, cLabel in enumerate(columnLabels):
            sheet.write(row, col+ii, cLabel, formats['header'])
        row +=1
        _colA += ii

    if rowLabels:
        for ii, rLabel in enumerate(rowLabels):
            fmt = formats['stdULine'] if sepLines[ii] else formats['std']
            sheet.write(row+ii, col, rLabel, fmt)
        _rowA = ii+row

    for ii, dataRow in enumerate(cellData):
        for jj, cell in enumerate(dataRow):
            sheet.write(row+ii, sCol+jj, cell, valueFmt(cell, sepLines[ii]))
        _colB = max(_colB, sCol +jj)

    _rowB = ii+row
    row = max(_rowA, _rowB) + 1
    col = max(_colA, _colB)
    return row, col

def descriptiveCellData(matItem, pstSamples='all', latex=True):
    def conversion(name, value):
        return outputConversions(name, value, latex=latex)

    pstSamples, sampleIdxs = _getSampleNamesAndIndices(matItem, pstSamples)

    if not sampleIdxs:
        msg('No pstSamples left to plot when using %s as filter'
            % str(pstSamples))
        return

    descriptiveMatType, descriptiveMatTypeLong = characteristics(matItem)

    cellTextA = [['Type', descriptiveMatTypeLong]]

    rho = getattr(matItem, 'density', None)
    cellTextA.append(conversion('density', rho))

    if descriptiveMatType not in ['MC', 'SPMC', 'elastic', 'empty',]:
        cellTextA.append(conversion('UCSi', matItem.UCSi))
        GSI = getattr(matItem, 'GSI', None)
        cellTextA.append(conversion('GSI', GSI))

    if descriptiveMatType not in ['SPMC', 'elastic', 'empty',]:
        anisoN = getattr(matItem, 'anisoN', 1.)
        anisoS = getattr(matItem, 'anisoS', 1.)
        cellTextA.append(conversion('anisoN',anisoN))
        cellTextA.append(conversion('anisoS',anisoS))

    cellTextB = []
    def getDataForSamples(variable):
        name,_ = conversion(variable,0)
        vals = []
        for ii in sampleIdxs:
            try:
                vals.append(getattr(matItem[ii], variable))
                continue
            except AttributeError:
                pass

            try:
                vals.append(getattr(matItem[ii].elastic, variable))
                continue
            except AttributeError:
                pass

            try:
                vals.append(getattr(matItem[ii].yieldSurf, variable))
                continue
            except AttributeError:
                pass

        return [name,] + [conversion(variable,val)[1] for val in vals]

    if descriptiveMatType == 'elastic':
        cellTextA.append(getDataForSamples('E')[:2])
        cellTextA.append(getDataForSamples('nu')[:2])

    elif not descriptiveMatType == 'empty':
        cellTextB.append(getDataForSamples('pst'))
        cellTextB.append(getDataForSamples('E'))
        cellTextB.append(getDataForSamples('nu'))

    if descriptiveMatType == 'MC':
        fit = QuasiMohrCoulomb.cohsFrictionfromExtendedMenetreyWillamFit
        phis, cohs = [], []
        for ii in sampleIdxs:
            c,p = fit(matItem[ii].yieldSurf, nSamples=4, returnFit=False)
            if c < 10.: #set zero for very small cohesion values
                c = 0.
            phis.append(p)
            cohs.append(c)
        cellTextB.append([conversion('c',0)[0],] +
                         [conversion('c',v)[1] for v in cohs])
        cellTextB.append([conversion('phi',0)[0],] +
                         [conversion('phi',v)[1] for v in phis])
        cellTextB.append(getDataForSamples('dilation'))

    elif descriptiveMatType not in ['MC', 'SPMC', 'elastic', 'empty',]:
        cellTextB.append(getDataForSamples('s'))
        cellTextB.append(getDataForSamples('mb'))
        cellTextB.append(getDataForSamples('a'))
        cellTextB.append(getDataForSamples('e'))
        cellTextB.append(getDataForSamples('dilation'))

    return cellTextA, cellTextB


def plotPropertyTables(matItem, tax1=None, tax2=None, pstSamples='all',
                       latex = True):

    def conversion(name, value):
        return outputConversions(name, value, latex=latex)

    # preparing plot axes if not present
    if not tax1:
        fig1, tax1 = plt.subplots()
        tax1.axis('off')
    else:
        fig1 = tax1.get_figure() # to return a reasonalbe value for fig


    if not tax2:
        fig2, tax2 = plt.subplots()
        tax2.axis('off')
    else:
        fig2 = tax2.get_figure() # to return a reasonalbe value for fig

    taxHandles = [tax1, tax2,]
    figHandles = [fig1, fig2]

    sampleNames = pstSampleNames(matItem)
    pstSamples, sampleIdxs = _getSampleNamesAndIndices(matItem, pstSamples)

    cellTextA, cellTextB = descriptiveCellData(matItem, pstSamples=pstSamples,
                                               latex=latex)
    tabHandles = []
    if cellTextA:
        rowLabelsA = [ct[0] for ct in cellTextA]
        colLabelsA = None
        cellTextA = [map(valueToFmt, ct[1:]) for ct in cellTextA]
        tab1 = _setTableMatplotlib(tax1, 'General Properties',
                         cellText=cellTextA,
                         rowLabels=rowLabelsA,
                         colLabels=colLabelsA)
        tabHandles.append(tab1)

    if cellTextB:
        colLabelsB = [sampleNames[ii] for ii in sampleIdxs]
        rowLabelsB = [ct[0] for ct in cellTextB]
        cellTextB = [map(valueToFmt, ct[1:]) for ct in cellTextB]
        tab2 = _setTableMatplotlib(tax2, 'LR2 Parameters',
                        cellText=cellTextB,
                        rowLabels=rowLabelsB,
                        colLabels=colLabelsB)
        tabHandles.append(tab2)

    return figHandles, taxHandles, tabHandles


#%% Excel functions

def _xlsxColumnString(n):
    "Service fuction that returns excel column name for column-index n"
    string = ""
    while n > 0:
        n, remainder = divmod(n - 1, 26)
        string = chr(65 + remainder) + string
    return string



def matItemToXlsx(matItem, workBook, plot=True, pstSamples='all', **kwargs):
    """
    Writes a general and sample-dependent properties to a specified
    xlsx-file. It creates a worksheet named as the fullname of the
    material (e.g. M01.rock.MY_HOST).

    @param workBook: xlsx-filename or xlsxwriter.Workbook - object

    @param plot: creates a S3-S1-Plot for all pst-samples on worksheet
        default=True

    @kwarg plotMaxS3: upper Range for S3 in compression meridian plot.
        If not specified, 20 times the absolut value of minimal
        tensileStrength of all PST-samples is choosen

    @kwarg plotNSamples: number of samples over p-range

    @kwarg plotBisection: plots bisection line

    @note: this function requires xlsxwriter (and numpy if plot==True)
    """
    try:
        import xlsxwriter as xls
    except ImportError:
        raise ImportError("xlsxwriter is not installed. Use:\n" +
                          ">>> pip install XlsxWriter\n" +
                          "on regular python or\n" +
                          ">>> conda install -c anaconda xlsxwriter\n"
                          "if you are using anaconda")

    if (not hasattr(matItem, '_matType')
        or
        not matItem._matType in [PlasticLevkovitchReuschAniso,
                                PlasticLevkovitchReusch,
                                SinglePlaneMohrCoulomb]):
        raise ValueError("Can only deal with MatItemPLR and MatItemPLRAniso "
                         "(from matpropscollection) but got %s"
                         % type(matItem))

    if not isinstance(workBook, xls.Workbook):
        workBook = xls.Workbook(workBook, {'nan_inf_to_errors': True})
        closeIfDone = True
    else:
        closeIfDone = False

    sheetName = matItem.fullName
    matName = sheetName.split('.')[-1]
    if not matName == sheetName:
        originName = '.'.join(sheetName.split('.')[:-1])
    else:
        originName = ''

    sheetName = '.'.join(sheetName.split('.')[1:])
    if not sheetName:
        sheetName = matName
    sheetName = sheetName.replace('rock.', 'R.').replace('fill.', 'F.')
    sheet = workBook.add_worksheet(sheetName)

    try:
        computator,_ = matItem.isFromComputator()
    except AttributeError:
        computator = None
    description, descriptionLong = characteristics(matItem)

    ## cellFormats
    formats = _setExcelFormats(workBook)

    sheet.write('A1', matName, formats['title'])
    sheet.write('E1', originName, formats['subTitle'])

    sampleNames = pstSampleNames(matItem)
    pstSamples, sampleIdxs = _getSampleNamesAndIndices(matItem, pstSamples)

    cellTextA, cellTextB = descriptiveCellData(matItem, pstSamples=pstSamples,
                                               latex=False)

    row, column = 2, 0
    if cellTextA:
        columnLabelsA = [ct[0] for ct in cellTextA]
        cellDataA = [[ct[1] for ct in cellTextA],]
        row,_ =_setExcelSheetTable(sheet, cellDataA, formats,
                                   row=row, col=column,
                                   rowLabels=None, columnLabels=columnLabelsA,
                                   title='General Properties')

    if cellTextB:
        rowLabelsB = [sampleNames[ii] for ii in sampleIdxs]
        colLabelsB = ['',] + [ct[0] for ct in cellTextB]
        cellDataB = [ct[1:] for ct in cellTextB]
        cellDataB = map(list, zip(*cellDataB)) # transpose
        row,_ =_setExcelSheetTable(sheet, cellDataB, formats,
                                   row=row+2, col=column,
                                   rowLabels=rowLabelsB,
                                   columnLabels=colLabelsB,
                                   title='Plastic Strain Dependend Properties')


    if not computator and not description in ['MC', 'arbitrary',]:
        plot = False

    if plot:
        blues = ['#08519c', '#4292c6', '#9ecae1', '#c6dbef']

        defaultLStyle = [{'color': blues[0], 'width': 1.5,}, ] + \
                        [{'color': col, 'width': 1., 'dash_type': dt }
                           for dt in ['long_dash', 'dash', 'round_dot',
                                      'square_dot', 'dash_dot',
                                      'long_dash_dot','long_dash_dot_dot',]
                           for col in blues ]

        if matItem.UCSi == 0:
            matItem = matItem.copy()
            matItem.UCSi = 1.

        isoTens = [s.yieldSurf.getIsoTensileStrength() for s in matItem]
        minP = min(isoTens)
        maxP = kwargs.get('plotMaxS3', 20*abs(minP))
        nSamples = kwargs.get('plotNSamples', 32)
        plotBisection = kwargs.get('plotBisection', True)

        # setup chart
        chart = workBook.add_chart({'type': 'scatter',
                                    'subtype': 'smooth'})

        chart.set_title({'name':matName,
                         'name_font':{'size':11}})
        chart.set_x_axis({'name': u'S₃ / MPa'})
        chart.set_y_axis({'name': u'S₁ / MPa',
              'major_gridlines': {'visible': True}})

        limS1, limS3 = [np.nan, np.nan], [np.nan, np.nan]
        startCol = 26
        for ii, sample in enumerate(matItem):
            theta = np.deg2rad(60)
            ys = sample.yieldSurf

            pMax,qMax = getPAndQFromThetaAndS(ys, theta,
                                              S3=maxP,
                                              miningNotation=True)

            ps = np.linspace(ys.getIsoTensileStrength(), pMax, nSamples)

            # calculate stressValues
            thetaC = np.full_like(ps, theta)

            qC = ys.ptOnYieldSurf(p=ps, theta=thetaC)

            S1,S2,S3 = pQThetaToEigs(ps,qC,thetaC,miningNotation=True)

            S1 /= 1E6 # to MPa
            S3 /= 1E6

            limS1[0] = np.nanmin([np.nanmin(S1), limS1[0]])
            limS1[1] = np.nanmax([np.nanmax(S1), limS1[1]])
            limS3[0] = np.nanmin([np.nanmin(S3), limS3[0]])
            limS3[1] = np.nanmax([np.nanmax(S3), limS3[1]])

            # replacing NaNValues by empty cells
            nanIdx = np.isnan(S1)
            S1 = S1.astype(object)
            S3 = S3.astype(object)
            S1[nanIdx] = ''
            S3[nanIdx] = ''

            ## writing series to excelSheet

            # create label and headerData
            label = u'%d: εₚₗₐₛₜ=%.2f%%' % (ii, 100*sample.pst)
            sheet.write(0, startCol + 2*ii, label, formats['title'])
            sheet.write(1, startCol + 2*ii+1, 'S1', formats['header'])
            sheet.write(1, startCol + 2*ii, 'S3',formats['header'])

            c0 = _xlsxColumnString(startCol + 2*ii+1)
            c1 = _xlsxColumnString(startCol + 2*ii +2)

            #write data
            sheet.write_column(2, startCol + 2*ii+1, S1)
            sheet.write_column(2, startCol + 2*ii, S3)
            sheet.write_comment(2, startCol + 2*ii,
                                'isotropic tensile strength (S1=S3)')
            sheet.write_comment(2, startCol + 2*ii+1,
                                'isotropic tensile strength (S1=S3)')

            # create chart
            # def-Strings for chards
            categoriesStr = "=%s!$%s$%d:$%s$%d" % (sheetName,
                                               c0, 3, c0, 3 + len(S1))
            valuesStr = "=%s!$%s$%d:$%s$%d" % (sheetName,
                                               c1, 3, c1, 3 + len(S1))

            seriesDict = {'values': valuesStr,
                          'categories':categoriesStr,
                          'name': label,}

            try:  # fallback if ii exeeds defaultStyleList
                seriesDict.update({'line':defaultLStyle[ii]})
            except IndexError:
                pass

            chart.add_series(seriesDict)

        if plotBisection:
            sheet.write(5 + nSamples, startCol,
                        'Bisecting Line', formats['header'])
            sheet.write_column(6 + nSamples, startCol, limS3)

            c0 = _xlsxColumnString(startCol + 1)

            categoriesStr = "=%s!$%s$%d:$%s$%d" % (
                sheetName,
                c0, 7 + nSamples,
                c0, 8 + nSamples)

            chart.add_series({'values': categoriesStr,
                              'categories':categoriesStr,
                              'name':'Bisector',
                              'line':{'color': '#C0C0C0', 'width': .75}})

        sheet.insert_chart(row+3, 0, chart)

    if closeIfDone:
        workBook.close()


#%% Testfunctions

def _testPlotAniso():
    plr = PlasticLevkovitchReuschAniso()
    plr.fromComputator(50E6, 55, anisoS=1.2, anisoN=1.191)
    plotAniso(plr)
    plt.show()


def _testPlotYieldSurf():
    plr = PlasticLevkovitchReuschAniso()
    plr.fromComputator(50E6, 55, anisoS=1.2, anisoN=1.15, version='V14.1')

    plotYieldSurf(plr, S3Max=10E6, thetaRange=True,pstSamples=['PEAK', 'RES'])
    plotYieldSurf(plr, S3Max=10E6, thetaRange=True, view='pq')
    plt.show()

def _testPlotPropertyTables():
    plr = PlasticLevkovitchReuschAniso()
    plr.fromComputator(50E6, 55, anisoS=1.2, anisoN=1.08, version='V14.1')
#    print plr
#    plr.pop()
#    print plr
    # plr.UCSi =10000*plr.UCSi # test elastic

    plr.a = 1 # test MC
    plotPropertyTables(plr, latex=True, pstSamples=['PEAK','TRANS', 'RES'])
    plt.show()

if __name__ == '__main__':
    print('No Syntax Errors')
    plt.close('all')
#    _testPlotYieldSurf()
#    _testPlotAniso()
    _testPlotPropertyTables()
