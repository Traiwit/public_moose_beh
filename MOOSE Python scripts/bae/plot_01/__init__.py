"""package to plot data using matplotlib

Intended to be imported like a module, importing L{Mainplot}, L{Subplot},
L{Singleplot} and L{DefaultPlots}
from L{bae.plot_01.mpl_core.py}.

Usage Mainplot, Subplot
=====
 >>> from bae.plot_01 import Mainplot, Subplot
 >>>
 >>> # read data, for example for multiple points, vars, frames from CSV:
 >>> pvf = PVFContainer(...)
 >>>
 >>> # create Mainplot instance
 >>> mp = Mainplot(figNb=1, gridRowsCols=(3,1), gridScale=None, dpi=None)
 >>>
 >>> # add Subplot instance, data for one particular point only -> single curve
 >>> sp1 = Subplot(name="subplot1", mpInst=mp, mpGridPos=(0,0),mpGridSpan=(1,1))
 >>> sp1.addData(pvf.ptData[ptKey], xId='t', yId='LOGP')
 >>> mp.addSubplot(sp1)
 >>>
 >>> # add Subplot instance, data for one particular point only,
 >>> # but two vars -> two curves
 >>> sp2 = Subplot(name="subplot2", mpInst=mp, mpGridPos=(1,0),mpGridSpan=(1,1))
 >>> sp2.addData(pvf.ptData[ptKey], xId='t', yId='S1',
 >>>             lp=dict([ ('c','C0'), ('ls','-'), ]) )
 >>> sp2.addData(pvf.ptData[ptKey], xId='t', yId='S3',
 >>>             lp=dict([ ('c','C0'), ('ls',':'), ]) )
 >>> mp.addSubplot(sp2)
 >>>
 >>> # add Subplot instance, data for multiple points, single var
 >>> # -> as many curves as points
 >>> sp3 = Subplot(name="subplot3", mpInst=mp, mpGridPos=(2,0),mpGridSpan=(1,1))
 >>> ptKeys = pvf.sortPtKeysById()
 >>> for k in ptKeys:
 >>>     # by default using color-cycler
 >>>     sp3.addData(pvf.ptData[k], xId='t', yId='S1_MAG')
 >>> mp.addSubplot(sp3)
 >>>
 >>> # finally:
 >>> # to simply plot all data
 >>> mp.plot()
 >>> # and/or, to plot the data as progress for each frame with options to show
 >>> # the current frame with a 'o' marker (show_marker=True), and/or to how
 >>> # the rest of the data for subsequent frames as grey-dashed line
 >>> # (show_rest=True)
 >>> mp.plotprogress(range(0,100,10), show_rest=True, show_marker=True)

Usage DefaultPlots (see exampleDefaultPlots in bae.plot_01.mpl_core.py)
====
 >>> from bae.plot_01 import DefaultPlots
 >>> import numpy as np
 >>>
 >>> #create some random data
 >>> angleA = np.mod(np.random.randn(100,1),2*np.pi)
 >>> angleB = np.mod(np.random.randn(100,1),np.pi/2)
 >>> mag = 10*np.random.randn(100,1)
 >>> xyz = 100*np.random.randn(100,3)
 >>>
 >>> #setup plots
 >>> dP = DefaultPlots()
 >>> dP.layout['mode']='onpaper'
 >>> fig, axs = dP.createAxes(nRows=2,nCols=2,
 >>>                          colScale=[1,.8],
 >>>                          shemisphere=[0,2])
 >>>
 >>> #define cell data for tables in header
 >>> rangeStr = dP.getRangeString(xyz)
 >>> infoCell = [np.array([['Project:','MyProject'],
 >>>                       ['Info:','some random data']
 >>>                      ]),
 >>>             np.array([['East [m]:',rangeStr[0]],
 >>>                       ['North [m]:',rangeStr[1]],
 >>>                       ['Depth (range) [m]:',rangeStr[2]]
 >>>                      ])]
 >>> tabs,_,_ = dP.createInfoHeader(fig,infoCell,tabWidth=[.4,.6])
 >>>
 >>> #set titles and axes labels
 >>> dP.setByColIdx(axs,0,'title','$S_{%d}$',
 >>>                valList=[[1],[2]],
 >>>                addargs='loc="left"')
 >>> #plot data as usual
 >>> axs[0,0].scatter(angleA,angleB,c=mag,s=mag)
 >>> axs[0,1].plot(xyz[:,2],mag,'og-')
 >>> axs[1,0].scatter(angleA,angleB,c='r',s=xyz[:,2])
 >>> axs[1,1].plot(xyz[:,2],mag,'-')
 >>> fig.savefig('/my/save/path/myplot.png', bbox_inches='tight')
"""

__version__ = "1.5"

_version_history_ = """
plot_01 package is intended to make easy, consistent looking matpotlib plots.

1.0 AF new
1.1 GP accept matplotlib-versions prior to 2.0, though without plot styles;
   added arguments to plot and plotprogress as alternative to the
   Mainplot.getOutFilename() functionality. Added Singleplot class.
1.2 TR New class DefaultPlots to create and manage subplot (grids).
    It is less capsuled (and less comfortable) than the Mainplot-Class but it
    preserves the functionality to access all figure, axes and plot objects by
    their handle (as intended in matplotlib).
    Added (pseudo)projections for reversed polarplot and southernhemisphere.
1.3 TR new: DefaultPlots.changeHeaderTableCell and DefaultPlots.getRangeString
       fixed: DefaultPlots.setDepthAxes 'uncertain' behaviour of z0
1.4 TR southernHemisphere and invPolarAxes now can be used directly with the
    'raw'-PolarAxis-plot-methods plot, scatter, contour and contourf
1.5 TR added Harrisson-plot and MapAxes projection, CurvedText-class and
    truncatedColormap
"""

__todo__ = """
Geros suggestions for incompatible changes:
- remove getOutFilename option or rename to getOutFileName (?)
- rename plotprogress to plotProgress (or something else?)
- remove the verbose option in plotprogress, it's redundant to msg-debugLevel!
- Subplot.addData: Why take a dict and a key as argument instead of just the
  corrsponding value?
"""

from bae.plot_01.mpl_core import Mainplot, Subplot, Singleplot, DefaultPlots
