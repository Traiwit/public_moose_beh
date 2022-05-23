import os
import numpy as np  # for Harrisson curve calculations
from matplotlib.image import imread
from matplotlib.patches import Rectangle

from bae.plot_01 import Subplot
from bae.misc_01 import Container


### ##########################################################################
### ##### SPECIFIC  S U B P L O T  VARIATIONS ################################ 
### ##########################################################################

### --------------------------------------------------------------------------
### T I T L E
### --------------------------------------------------------------------------
class SubplotTitle(Subplot):

    def __init__(self, name, mpInst, mpGridPos, mpGridSpan):
        Subplot.__init__(self, name, mpInst, mpGridPos, mpGridSpan)
        ### 
        self.setTitleData()
        return 
        
    def setTitleData(self, 
            info=Container(
                ptId="TESTSITE",
                ptNb=0,
                ptCoords=(0.,0.,0.),
                coordsOffset=(0.,0.,0.),
                ),
            ):
        ### check for translationVector, defaults to (0.,0.,0.)
        try:
            siteCoordsOffset = info.coordsOffset
        except AttributeError:
            siteCoordsOffset = (0.,0.,0.)
        ### check for relevant info
        try:
            siteName = info.ptId
            siteNumber = info.ptNb
            siteCoords = info.ptCoords
        except AttributeError:
            raise Exception("Missing information for SubplotTitle!")
        else:
            self.titleData = Container(
                siteName=siteName,
                siteNumber=siteNumber,
                siteCoords=[ c-dc for c,dc in zip(siteCoords,siteCoordsOffset) ],
                )
        
    def setup(self):
        aspectratio = self.mpGridSpan[1]/self.mpGridSpan[0]
        h,w = 5.,5.*aspectratio
        self.ax.clear()
        self.ax.set_frame_on(False)
        self.ax.set_axis_off()
        self.ax.set_xlim(0,w)
        self.ax.set_ylim(0,h)
        ### text properties for descriptions
        props1 = dict(
            ha='left',
            va='center',
            #color=None,
            #style=None,
            fontsize='x-small',
            fontweight='bold',
            #bbox={
            #    'facecolor': 'white',
            #    'edgecolor': 'black',
            #    'pad': 2,
            #    'lw': 0.5,
            #    },
            )
        props2 = props1.copy()
        props2.update(dict(
            fontweight='normal',
            ))
        props2b = props2.copy()
        props2b.update(dict(
            ha='right',
            ))
        propsBox = dict(
            lw=0.5, 
            facecolor='white',
            edgecolor='black',
            )
        x = 0.
        dx = w/3.
        dy = 0.90
        self.ax.add_artist(Rectangle((x,4.05),dx,dy, **propsBox))
        self.ax.text(x+0.1,4.5, "SITE NAME", **props1)
        self.ax.add_artist(Rectangle((x,3.05),dx,dy, **propsBox))
        self.ax.text(x+0.1,3.5, "SITE NUMBER", **props1)
        self.ax.add_artist(Rectangle((x,2.05),dx,dy, **propsBox))
        self.ax.text(x+0.1,2.5, "EAST [m]", **props1)
        self.ax.add_artist(Rectangle((x,1.05),dx,dy, **propsBox))
        self.ax.text(x+0.1,1.5, "NORTH [m]", **props1)
        self.ax.add_artist(Rectangle((x,0.05),dx,dy, **propsBox))
        self.ax.text(x+0.1,0.5, "RL [m]", **props1)
        
        x  = w/3.+0.1*aspectratio
        x2 = w/3.+0.1*aspectratio+w/6.
        dx = w-x-0.1
        self.ax.add_artist(Rectangle((x,4.05),dx,dy, **propsBox))
        self.ax.text(x+0.1,4.5, str(self.titleData.siteName), **props2)
        self.ax.add_artist(Rectangle((x,3.05),dx,dy, **propsBox))
        self.ax.text(x+0.1,3.5, "%s"%self.titleData.siteNumber, **props2)
        self.ax.add_artist(Rectangle((x,2.05),dx,dy, **propsBox))
        self.ax.text(x2+0.1,2.5, "%.1f"%self.titleData.siteCoords[0], **props2b)
        self.ax.add_artist(Rectangle((x,1.05),dx,dy, **propsBox))
        self.ax.text(x2+0.1,1.5, "%.1f"%self.titleData.siteCoords[1], **props2b)
        self.ax.add_artist(Rectangle((x,0.05),dx,dy, **propsBox))
        self.ax.text(x2+0.1,0.5, "%.1f"%self.titleData.siteCoords[2], **props2b)

        return
        
### --------------------------------------------------------------------------
### L O G O
### --------------------------------------------------------------------------
class SubplotLogo(Subplot):

    def __init__(self, name, mpInst, mpGridPos, mpGridSpan):
        Subplot.__init__(self, name, mpInst, mpGridPos, mpGridSpan)
        self.setLogo(os.path.join(
            os.path.dirname(__file__),
            "logoBE_300_white.png"
            ))
        self.imgLogoSize = (300,170)
        return 
        
    def setLogo(self, imgFilename):
        try:
            self.imgLogoFilename = imgFilename
            self.imgLogoArray = imread(imgFilename)
        except Exception as e:
            raise Exception("Any Error reading the logo-image!\n%s"%e)
        return 
        
    def setup(self):
        self.ax.clear()
        self.ax.set_frame_on(False)
        self.ax.set_axis_off()
        figdpi = self.mpInst.fig.dpi
        figdpi = self.mpInst.fig.get_dpi()
        figsizeIn = self.mpInst.fig.get_size_inches()
        figsizePx = map(int, [figdpi*v for v in figsizeIn])
        gridShape = map(float, self.mpInst.gridRowsCols)
        gridShape.reverse()
        gridPos = map(float, self.mpGridPos)
        gridPos.reverse()
        gridSpan = map(float, self.mpGridSpan)
        gridSpan.reverse()
        logoOffsetPx = map(int, [ 
            (0.0+(gridPos[0]+0.0        )/gridShape[0])*figsizeIn[0]*figdpi,
            (1.0-(gridPos[1]+gridSpan[1])/gridShape[1])*figsizeIn[1]*figdpi
            ])
        #msg("logoOffsetPx = %s"%str(logoOffsetPx))
        
        if self.imgLogoArray is not None:
            #msg("imgLogoArray = %s"%str(self.imgLogoArray))
            self.ax.figure.figimage(
                self.imgLogoArray, 
                #xo=logoOffsetPx[0],
                xo=figsizePx[0]-self.imgLogoSize[0],
                #yo=logoOffsetPx[1], 
                yo=figsizePx[1]-self.imgLogoSize[1], 
                resize=False,
                #alpha=1.,
                #origin='lower',
                zorder=1
                )
        return
        
        
### --------------------------------------------------------------------------
### H A R R I S S O N   C H A R T 
### --------------------------------------------------------------------------
class SubplotHarrisson(Subplot):
    
    #def __init__(self, name, mpInst, mpGridPos, mpGridSpan):
    #    Subplot.__init__(self, name, mpInst, mpGridPos, mpGridSpan)
    #    return 
    
    def addData(self, ptDataDict, xId, yId, lc=None, lp=None):
        Subplot.addData(self, ptDataDict, xId, yId, lc, lp)
        ### difference: now only consider values with 'STATUS'=1
        try:
            statusMask = [ 
                True if int(v+0.1)==1 else False 
                for v in ptDataDict['STATUS']
                ]
        except KeyError:
            msg("WARNING: variable 'STATUS' not found, values not masked.")
            pass
        else:
            xValsMasked = []
            yValsMasked = []
            for xvals,yvals in zip(self.xValues,self.yValues):
                xValsMasked.append([ v for v,m in zip(xvals,statusMask) if m ])
                yValsMasked.append([ v for v,m in zip(yvals,statusMask) if m ])
            self.xValues = xValsMasked
            self.yValues = yValsMasked
            self.computeAxLimitsAuto()
        return

    def setup(self):
        Subplot.setup(self)
        ### differences to default: 
        ### 1. add harrisson curves as dotted black lines
        x = np.linspace(0.,7.,50)
        yDict = {
            "veryslightly":      -0.3000*x**2.-0.0500*x+0.4500,
            "slightly":          -0.2654*x**2.-0.0404*x+0.7466,
            "mod_severe":        -0.1377*x**2.-0.0071*x+1.5465,
            "severe_verysevere": -0.0669*x**2.-0.0009*x+3.0026,
            }
        for k,y in yDict.iteritems():
            self.ax.plot(x,y, c='black', ls=':')
            
        ### 2. add shadings between curves [color is  color=(R,G,B,alpha) ]
        self.ax.fill_between(x,yDict["veryslightly"],     0.,                         color=(0.9,0.9,0.9,1.0) )
        self.ax.fill_between(x,yDict["slightly"],         yDict["veryslightly"],      color=(0.8,0.8,0.8,1.0) )
        self.ax.fill_between(x,yDict["mod_severe"],       yDict["slightly"],          color=(0.7,0.7,0.7,1.0) )
        self.ax.fill_between(x,yDict["severe_verysevere"],yDict["mod_severe"],        color=(0.6,0.6,0.6,1.0) )
        self.ax.fill_between(x,3.5,                       yDict["severe_verysevere"], color=(0.5,0.5,0.5,1.0) )

        ### 3. add annotations
        self.ax.text(5.50,3.25, 'Severe damage to\nVery severe damage', ha='center', va='center', multialignment='center', fontsize='xx-small')
        self.ax.text(4.50,0.25, 'Moderate damage to\nSevere damage', ha='center', va='center', multialignment='center', fontsize='xx-small')
        self.ax.text(0.10,1.00, 'Slight\ndamage', ha='left', va='bottom', multialignment='left', fontsize='xx-small')
        #self.ax.text(-1.50,0.50, 'Very slight\ndamage', color='grey', ha='left', va='center', multialignment='left', fontsize='xx-small')
        self.ax.annotate('Very slight\ndamage', xy=(0.1,0.6), xytext=(-1.70,0.80), 
            ha='left', va='top', multialignment='left', fontsize='xx-small',
            #arrowprops=None,
            arrowprops=dict(fc='black',lw='.05',arrowstyle='->'),
            )
        #self.ax.text(0.50,-0.50, 'Negligible\ndamage', ha='center', va='center', multialignment='center', fontsize='xx-small')
        self.ax.annotate('Negligible\ndamage', xy=(0.5,0.1), xytext=(0.5,-0.5), 
            ha='center', va='center', multialignment='center', fontsize='xx-small',
            #arrowprops=None,
            #arrowprops=dict(fc='black',lw='.05',arrowstyle='simple'),
            arrowprops=dict(fc='black',lw=0.2, width=.20,headwidth=.20,headlength=.40),
            )

        return

        
### --------------------------------------------------------------------------
### M A P C H A R T 
### --------------------------------------------------------------------------
class SubplotMapchart(Subplot):

    def __init__(self, name, mpInst, mpGridPos, mpGridSpan):
        Subplot.__init__(self, name, mpInst, mpGridPos, mpGridSpan)
        ### differences to default:
        self.setMapImage()
        self.setMarkerPosition()
        return 
        
    def setMapImage(self, mapFilename="map.png", extend=None):
        """for further keywords, see 
        U{https://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.imshow}
        """
        self.mapFilename = mapFilename
        self.mapExtend = extend
        return
        
    def setMarkerPosition(self, ptCoords=(0.,0.,0.)):
        self.markerPosition = ptCoords
        return
    
    def setup(self):
        Subplot.setup(self)
        ### differences to default: 
        ### 1. read and show the image of the map
        #img = plt.imread(self.mapFilename)
        img = imread(self.mapFilename)
        self.ax.imshow(
            img,
            zorder=0,
            extent=self.mapExtend,
            )

        ### 2. show a marker at the ptCoord position
        x = self.markerPosition[0]
        y = self.markerPosition[1]
        self.ax.plot(x,y, c='black',marker='o')
        return
        
