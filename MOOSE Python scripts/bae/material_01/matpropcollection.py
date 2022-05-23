# -*- coding: utf-8 -*-
"""Classes for handeling multiple materials and fills.
Initially, theses classes where the API of the beRocker-tool which will
never get alive unfortunatly. Besides an easy access on material properties, an
easy to use plotting module and some other usefull functions, the main
intention of the beRocker tool was to log all changes made on materials and
material collections. That's why the classes L{MatItem}, L{MatItemCollection},
L{MatItemAssignments} and L{MatPropsCollection} do have special functions for
adding, changing, removing, ... which could be neglected if no logging is
required.
So there are two possible scenarios to use the following classes or functions:
    1 Uses the provided functions (only) to enable logging (as used in
     beRocker)
    2 Use regular python syntax or a mixture with scenario (1)

Assume you'll have to move a material from one collection to another:

>>> matCollectionA.moveItemToOtherCollection(matCollectionB, 'matFromA') # (1)
>>> matCollectionB['matFromA'] = matCollectionA.pop('matFromA') #scenario (2)

"""

import os
import datetime
import warnings
from collections import OrderedDict
from itertools import izip

from math import log10
from copy import deepcopy

from bae.misc_01 import Container
from material import (PlasticLevkovitchReusch,
                      PlasticLevkovitchReuschAniso,
                      SinglePlaneMohrCoulomb,
                      ElasticLinear,)
from yieldsurfaces import getPAndQFromThetaAndS
from models import restoreComputatorInputs

from bae.log_01 import log, msg



def characteristics(matItem):
    """Returns a major characteristic of material like elastic or
    computator if identified. The returned tuple holds a (sortDescrition,
    longDescription) of the identified characteristic.
    """
    try:
        if matItem._matType == SinglePlaneMohrCoulomb:
            return ('SPMC', 'SinglePlaneMohrCoulomb')
    except AttributeError:
        if type(matItem) == SinglePlaneMohrCoulomb:
            return ('SPMC', 'SinglePlaneMohrCoulomb')

    if not len(matItem):
        return ('empty', 'empty PLR')

    if matItem[0].yieldSurf.UCSi > .9E9:
        return ('elastic', 'elastic PLR')
    try:
        if matItem[0].yieldSurf.a == 1.:
            return ('MC', 'MOHR-COULOMB PLR')
    except (IndexError, AttributeError):
        pass

    computator, _ = matItem.isFromComputator(firstMismatch = True)
    if computator is None:
        return ('arbitrary', 'arbitrary PLR')
    else:
        return (computator, 'Computator %s PLR' % computator )


##############################################################################
#{Logging for Collection classes
class LogEvent(object):
    """A single log event which is described by:
        - an event(name/description)
        - changes
        - a logLevel (0 is highest)
        - a creation time.
    The changes will be passed as keyWordArguments where the key describes the
    change/changed value and the argument is an abitrary object (e.g. a tuple
    holding old an new value).
    """
    timeFmt = "%I:%M:%S%p on %B %d, %Y"

    def __init__(self, eventName, level = 0, time=None, **kwargs):
        self.event = eventName
        self.level = level
        if time is None:
            self.time  = datetime.datetime.now()
        else:
            self.time = time
        for key, val in kwargs.iteritems():
            setattr(self, key, val)

    def toDict(self):
        return dict(eventName=self.event, level=self.level,
                    time=self.fTime(), data=self.getLogDataDict())

    @classmethod
    def fromDict(cls, input):
        event = input['eventName']
        level = input['level']
        time = input['time']
        time = datetime.datetime.strptime(time, cls.timeFmt)
        data = input['data']

        return cls(event, level=level, time=time, **data)

    def fTime(self):
        """Returns a formated time string of log-events creation date/time"""
        return self.time.strftime(self.timeFmt)

    def getLogDataDict(self):
        poper = ['event', 'time', 'level']
        return dict((k,v) for k, v in vars(self).iteritems() if k not in poper)

    def logToFormatedArray(self):
        """Returns a tuple of (logEventName, changes, logLevel, creationTime)
        """
        dat = [k + ': %s'%str(v) for k, v in self.getLogDataDict().iteritems()]
        return self.event, '\t'.join(dat), self.level, self.fTime()


class MetaData(object):
    """Container to store single L{LogEvent}s in. Each LogEvent gets stored
    successivly in MetaData.logs (simple list)
    """
    def __init__(self, creationStr = 'new'):
        """Creates an emty MetaData-Container. Each newly created MetaData
        object will hold a creation-L{LogEvent} of highest Level
        """
        self.logs = []
        self.logState = True
        self.log(creationStr, level = 0)

        self.comment = ''

    def toDict(self):
        return dict(logs=[log.toDict() for log in self.logs],
                    logState=self.logState,
                    comment=self.comment)

    @classmethod
    def fromDict(cls, input):
        out = cls()
        out.logs = [LogEvent.fromDict(log) for log in input['logs']]
        out.logState = input['logState']
        out.comment = input['comment']
        return out

    def logsToFormatedArray(self, level = None):
        """Will return a list of formated logArrays (see L{LogEvent.logToFormatedArray})
        up to a specified logLevel.
        """
        if level is None:
            level = 999
        arr = []
        for log in self.logs:
            dd = log.logToFormatedArray()
            if log.level <= level:
                arr.append(dd)
        return arr


    def log(self, event, **kwargs):
        """Stores a single L{LogEvent} in MetaData.logs
        """
        if self.logState:
            level = kwargs.pop('level',1)
            logEvent = LogEvent(event, level=level, **kwargs)
            self.logs.append(logEvent)


    def clearExceptRootLevel(self):
        """Removes all logEvents except those on root-level (0)
        """
        self.logs = [log for log in self.logs if log.level == 0]


    def clearAllBeforeLastLevel(self, level = 0, keepLevels = None):
        levels = [log.level for log in self.logs]
        try:
            lastLevel = len(levels) - 1 - levels[::-1].index(level)
        except ValueError:
            msg('No logLevel %d in logs. Clearing skipped.' % level)
            return

        lastLogs = self.logs[lastLevel:]
        if type(keepLevels) == bool and not keepLevels or keepLevels is None:
            self.logs = lastLogs
            return

        elif type(keepLevels) == bool and keepLevels:
            keepLevels = 0

        else:
            raise ValueError("Can't understand keepLevels=%s" %
                             str(keepLevels))

        keepLogs = [log for log in self.logs if log.level <= keepLevels]
        try:
            idx = keepLogs.index(lastLogs[0])
            self.logs = keepLogs[:idx] + lastLogs[1:]

        except ValueError:
            self.logs = keepLogs + lastLogs

#}#End Logging for Collection classes
##############################################################################


##############################################################################
#{MatItems
class MatItem(object):
    """Represents a single rock or fill material of a L{MatItemCollection}.

    It will get 'mixed-in' with L{material}s (
    L{PlasticLevkovitchReuschAniso}, L{PlasticLevkovitchReusch},
    L{SinglePlaneMohrCoulomb}) to L{MatItemPLR}, L{MatItemPLRAniso} and
    L{MatItemSPMC}. See L{itemInit} to see the __init__ function.
    """
    # _defaultOrder is used for output to sort (potential) material properties
    # in a reasonable order
    _defaultOrder = ['type', 'pst',
                     'density', 'E', 'nu', 'K', 'G',
                     'UCSi','GSI', 'D',
                     'cohesion', 'friction',
                     's', 'mb', 'alpha', 'e',]

    # the following values will be ommited for e.g. excel output
    _defaultHiddenValues = ['umatDelete', 'umatDepvarNum']
    _defaultDocHiddenValues = _defaultHiddenValues + ['dampingAlpha',
                                                      'critDamage',
                                                      'voidStiffnessRatio']

    # PST-sample class that describes this Material
    # will be set dynamically from dynamicMatItem
    _matType = None

    def copy(self):
        """Creates a deepcopy of MatItem"""
        out = self.__class__(self)

        for key,value in vars(self).iteritems():
            setattr(out, key, value)

        return out


    def toDict(self):
        extraMatPropKeys = ['GSI',
                            'density',
                            'dampingAlpha',
                            'critDamage',
                            'voidStiffnessRatio',]

        matData = dict(samples=[s.toDict() for s in self])
        for key in extraMatPropKeys:
            matData[key] = getattr(self, key, None)

        parentName = getattr(self._parent, 'name', None)
        return dict(matType=self.__class__.__name__,
                    parent=parentName,
                    name=self.name,
                    frozen=self._frozen,
                    metaData=self._meta.toDict(),
                    matData=matData,)


    @classmethod
    def fromDict(cls, input):
        implementedTypes = {'MatItemPLRAniso'   : MatItemPLRAniso,
                            'MatItemPLR'        : MatItemPLR,
                            'MatItemSPMC'       : MatItemSPMC,}
        input = dict(**input)
        matType  = implementedTypes[input.pop('matType')]
        parent   = input.pop('parent')
        name     = input.pop('name')
        metaData = MetaData.fromDict(input.pop('metaData'))
        frozen   = input.pop('frozen')

        matData  = dict(**input.pop('matData'))
        samples  = matData.pop('samples')

        pstSampleList = matType._matType.fromDict(dict(samples=samples))
        matItem = matType(pstSampleList, name=name, parent=parent, **matData)
        matItem._meta = metaData
        matItem._frozen = frozen

        return matItem


    def freeze(self, doIt=True):
        """Freezes a MatItem so it can't be edited using matpropscollection
        -fuctions anymore. Use self.freeze(False) to unfreeze.
        """
        self._frozen = doIt
        if doIt:
            self._meta.log('Material got frozen', level=0)
        else:
            self._meta.log('Material got unfrozen', level=0)


    @property
    def fullName(self):
        """Returns the complete name of a matItem e.g. 'M01.fill.defaultFill'
        """
        try:
            return self._parent.fullName + '.' + self.name
        except Exception:
            return self.name


    @property
    def character(self):
        return characteristics(self)


    @property
    def isAnisotropic(self):
        try:
            anisoN, anisoS = self.anisoN, self.anisoS
        except AttributeError:
            return False
        if anisoN==1. and anisoS==1.:
            return False
        else:
            return True

    def getProperty(self, prop, idx=None):
        """Gets a single material property of MatItem.

        @param prop: name of property

        @param idx: index of pst-sample. If None the property is assumed to
            be a general one (independend of pst like G, UCSi, ...). If idx is
            not specified/None for sample values (e.g. pst, mb, s, ...) an
            IndexError is raised.
        """
        try:
            return vars(self)[prop]
        except KeyError:
            if idx is None:
                raise('Material %s has no general property named %s. ' %
                      (self.name, prop) +
                      'Maybe you need to specify an Index.')

        return self._sampleVal(idx, prop)


    def _sampleVal(self, *args):
        '''Service function that sets/gets single pstSample to Material
        '''
        idx, prop = args[0], args[1]

        if len(args) == 2:
            setGet = lambda t: getattr(t, prop)
        else:
            val  = args[2]
            def setGet(t):
                oldVal = getattr(t, prop)
                setattr(t, prop, val)
                return oldVal

        # check if index is valid
        try:
            sample = self[idx]
        except IndexError:
            raise IndexError('Index for pst-sample out of range')
        except TypeError:
            raise IndexError('The Index %s does not refer to a pstSample' %
                             str(idx))

        if hasattr(sample, prop):
            #try to get/set general sampleData
            return setGet(sample)

        elif hasattr(sample.elastic, prop):
            #try to get/set elastic sampleData
            return setGet(sample.elastic)

        elif hasattr(sample.yieldSurf, prop):
            #try to get/set yieldSurf sampleData
            return setGet(sample.yieldSurf)

        else:
            raise AttributeError('Can not find attribute %s in sample.' %
                                 prop)


    def setProperty(self, prop, val, idx=None):
        '''Sets a material property of an unfrozen MatItem.

        >>> matItem.setProperty('UCSi', 40E6) # will set UCSi for all samples
        >>> matItem.setProperty('mb', .1, idx=1) # sets mb for 2nd sample
        >>> matItem.setProperty('mb', .1) # raises an IndexError

        @param prop: name of property to set

        @param val: new value of material property

        @param idx: index of pst-sample. If None, general (independend of pst)
            will be set (e.g. G, UCSi, ...). If idx is not specified/None for
            sample values (e.g. pst, mb, s, ...) an IndexError is raised

        @note: changes will get stored as eventLogs
        '''
        if self._frozen:
            msg('Material %s is frozen. No changes allowed.' % self.name)
            return

        sampleFixed = self._sampleType._fixedSampleIndex.keys()
        sampleVars  = set(self._sampleType._varSampleIndex.keys() +
                          ['E', 'nu'])

        if prop in sampleFixed:
            if not idx is None:
                msg(('WARNING! %s is defined as general material value and ' +
                         'thus will be set to all samples and not only %d') %
                         (prop, idx))

            for idx in range(len(self)):
                old = self._sampleVal(idx, prop, val)

            self._meta.log('changed %s for all samples' % prop,
                           old = old, new = val)

        elif prop in sampleVars:
            if idx is None:
                raise IndexError(('Property %s is defined as sample-value' +
                                  'and not as general material value.') % prop)

            old = self._sampleVal(idx, prop, val)

            self._meta.log('changed %s for sample %d' % (prop, idx),
                           old = old, new = val)

        elif hasattr(self, prop):
            old = getattr(self, prop)
            vars(self)[prop] = val

            self._meta.log('changed %s' % (prop),
                           old = old, new = val)

        else:
            raise AttributeError('Material %s has no attribute %s' %
                                 (self.name, prop))


    def getPropertyDicts(self, idx=None, order=None):
        '''Returns material properties of a MatItem as (two) dictionaries.

        The first holds the general material properties like density, UCSi,
        GSI. The second holds (if idx is not None) the elastic and yieldsurface
        properties of the specified pst-sample index.

        @param idx: index of pst-sample. If None the sampleMatPropsDict is
         empty
        @param order: specifies the order of properties in dictionary. This is
         usefull for more strutured outputs. If order is None, the
         MatItem.E{_}defaultOrder is used.
         If properties are not contained in order, they will be by their key
         and appended to the ordered dictionaries.
        @returns: generalMatPropsDict, sampleMatPropsDict
        '''
        if idx is None:
            idx = 0 # to get pst-fixed yieldSurfData
            onlyGeneral = True
        else:
            onlyGeneral = False

        general = vars(self.copy())
        general.update(vars(self[idx]))

        general.update({'E' : self[idx].elastic.E,
                        'nu' : self[idx].elastic.nu,})
        general.update(vars(self[idx].yieldSurf))

        # all pst-variable propertie-keys
        # elastic is a crank here because it can be defined by nu+E or G+K
        sampleKeys = set( self._sampleType._varSampleIndex.keys() +
                          ['E','nu'] )

        sample = OrderedDict()
        for key, val in general.items():
            #collect just float/int attributes
            if (key[0] == '_' or
                not (isinstance(val,float) or isinstance(val,int))):
                general.pop(key,None)

            if key in sampleKeys:
                sample[key] = val
                general.pop(key, None)

        general['type'] = self._matType.__name__

        if onlyGeneral:
            sample = {}

        if order is None:
            order = self._defaultOrder

        gen, sam = OrderedDict(), OrderedDict()
        for name in order:
            if name in general.keys():
                gen[name] = general.pop(name)
            if name in sample.keys():
                sam[name] = sample.pop(name)

        gen.update(sorted(general.iteritems()))
        sam.update(sorted(sample.iteritems()))

        return gen, sam


    def checkValidPst(self):
        """Service fuction that checks the correct order of plastic strain
        values of PST-samples.
        """
        testPst = [sample.pst for sample in self]
        deltas = [(j-i)>0 for i, j in zip(testPst[:-1], testPst[1:])]

        if not len(testPst) == len(set(testPst)):
            sampleError = {'pstError': 'duplicated PST-Values'}
        elif not all(deltas):
            sampleError = {'pstError': 'unsorted PST-Values'}
        else:
            sampleError = {}

        return sampleError



    def toXlsx(self, workBook, plot=True, skipOutput = None, **kwargs):
        """
        Writes a general and sample-dependent properties to a specified
        xlsx-file. It creates a worksheet named as the fullname of the
        material (e.g. M01.rock.MY_HOST).

        @param workBook: xlsx-filename or xlsxwriter.Workbook - object

        @param plot: creates a S3-S1-Plot for all pst-samples on worksheet
            default=True

        @param skipOutput: list of values to be excluded for output. If None
            (default), all values listed in L{defaultDocHiddenValues<plot.defaultDocHiddenValues>}
            will be skiped. An empty list [] will output all parameters.

        @kwarg plotMaxS3: upper Range for S3 in compression meridian plot.
            If not specified, 20 times the absolut value of minimal
            tensileStrength of all PST-samples is choosen

        @kwarg plotNSamples: number of samples over p-range

        @kwarg plotBisection: plots bisection line

        @note: this function requires xlsxwriter (and numpy if plot==True)
        """
        from plot import matItemToXlsx
        return matItemToXlsx(self, workBook, plot=plot, skipOutput=skipOutput,
                             **kwargs)


    def plotYieldSurfpq(self, **kwargs):
        """Plots hydrostatic pressure vs. vonMises stress of the yieldsurfaces
        of a L{material<bae.material_01.material.DefaultPst>} or a single
        L{yield surface<bae.material_01.yieldsurfacesExtendedMenetreyWillam>}
        object. See L{plot.plotYieldSurfpq} for more details.
        """
        from plot import plotYieldSurfpq
        return plotYieldSurfpq(self, **kwargs)


    def plotYieldSurfCompressionMeridian(self, **kwargs):
        """Plots the yieldsurfaces of a L{material<bae.material_01.material.DefaultPst>}
        or a single L{yield surface<bae.material_01.yieldsurfacesExtendedMenetreyWillam>}
        object on compression meridian (M{S{theta}=S{pi}/3}). See
        L{plot.plotYieldSurfCompressionMeridian} for more details.
        """
        from plot import plotYieldSurfCompressionMeridian
        return plotYieldSurfCompressionMeridian(self, **kwargs)


    def plotAniso(self, **kwargs):
        """Plots the anisotropic behaviour of a L{material<bae.material_01.material.DefaultPst>}
        (only PEAK-level) or a single
        L{yield surface<bae.material_01.yieldsurfacesExtendedMenetreyWillam>}.
            - first plot: compressive strength vs. orientation (S{phi}) for
                different confinements
            - second plot: compressive strength vs. confinement (in compression
                meridian) for M{S{phi}=0} and for the weakest orientation
                M{S{phi}_w(confinement)}
        See L{plot.plotAniso} for more details.
        """
        from plot import plotAniso
        return plotAniso(self, **kwargs)


    def plotPropertyTables(self, **kwargs):
        from plot import plotPropertyTables
        return plotPropertyTables(self, **kwargs)


    def createYieldSurfPlot(self, saveName=None, interactive=False,
                            propertyTables=False, pstSamples='all',
                            thetaRange=[True,True],
                            postPlotFun=None,
                            **kwargs):
        # import matplotlibs plt, pylab and gridspec via plot.py
        from plot import plt, pylab, gridspec, plotParams

        pylab.rcParams.update(plotParams)
        fig = plt.figure(figsize=(17.34/2.54, 11.00/2.54), dpi=150,)
        if propertyTables:
            # create 3columns, where the first is divided into 4 rows
            gs = gridspec.GridSpec(3, 3, height_ratios=[.3, .55, .15,],
                                         width_ratios=[3,3,3])

            # first 'axis-column', one axis for each table
            tax1 = plt.subplot(gs[:1, 0])
            tax2 = plt.subplot(gs[1:2, 0])
            tax3 = plt.subplot(gs[2:, 0])
            taxs = [tax1, tax2, tax3]
            for tax in taxs:
                tax.axis('off')

            # ploting axis
            ax1 = plt.subplot(gs[:, 1])
            ax2 = plt.subplot(gs[:, 2])
        else:
            taxs = []
            ax1 = fig.add_subplot(1,2,1)
            ax2 = fig.add_subplot(1,2,2)

        # tweak: reduce axis extents slightly so they dosn't override main-title
        fig.tight_layout(rect=[0.03, 0.03, 0.95, 0.95])
        fig.tight_layout(pad=2.0)
        # ax1.set_xlabel(r'$S_{Minor}$ in MPa')
        # ax1.set_ylabel(r'$S_{Major}$ in MPa')
        ax1.set_xlabel(r'$S_2=S_3$ in MPa')
        ax1.set_ylabel(r'$S_1$ in MPa')
        ax1.grid(True)
        ax2.set_xlabel(r'$p$ in MPa')
        ax2.set_ylabel(r'$\sigma_{Mises}$ in MPa')
        ax2.grid(True)
        fig.suptitle(self.fullName)

        if pstSamples is None:
            if len(self) > 1:
                pstSamples = ['PEAK', 'RES']
            else:
                pstSamples = ['PEAK']

        if thetaRange is None:
            thetaRange = [False,  # plot range in SmajorVsSminor
                          False]  # plot range in vonMisesVsP
        if not hasattr(thetaRange, '__iter__'):
            thetaRange = [thetaRange, thetaRange]

        # initialize handleLists
        axHs, pltHs, legHs, tabHs = [ax1, ax2], [], [], []

        if self.character[0] not in ['empty', 'elastic', 'SPMC']:
            h0 = self.plotYieldSurfCompressionMeridian(axh=ax1, thetaRange=thetaRange[0],
                                                  pstSamples=pstSamples,
                                                  **kwargs)
            h1 = self.plotYieldSurfpq(axh=ax2, thetaRange=thetaRange[1],
                                 pstSamples=pstSamples, **kwargs)
            pltHs += [h0[2], h1[2]]
            legHs += [h0[3], h1[3]]
        else:
            msg('Skipping yieldsurface plot for material %s (%s)' %
                (self.fullName, self.character[1]))

        if propertyTables:
            hTab = self.plotPropertyTables(tax1=tax1, tax2=tax2,
                                    pstSamples=pstSamples)
            axHs += taxs
            tabHs += hTab[2]

        if postPlotFun is not None:
            postPlotFun(fig, axHs, pltHs, legHs, tabHs, self)

        if interactive:
            fig.show()

        if saveName:
            outDir = os.path.dirname(saveName)
            if not os.path.isdir(outDir):
                os.makedirs(outDir)
            fig.savefig(saveName, bbox_inches='tight')

        return fig

###############################################################################
# init function for different MatItem-types

def itemInit(self, pLR, **kwargs):
    """__init__-fuction for dynamically created MatItems.

    @param pLR: material PST-sample of type L{material.PlasticLevkovitchReusch},
        L{material.PlasticLevkovitchReuschAniso} or
        L{material.SinglePlaneMohrCoulomb}
    @kwarg parent: set MatItemCollection as parent if available, default is
        None
    @kwarg name: name of MatItem, default is an empty string
    @kwarg creationString: description of log-event for newly created MatItem,
        default is 'new'
    @kwarg others: other kwargs will be set as attributes of MatItem. Use
        these to set the for 'umatDelete', 'umatDepvarNum', 'dampingAlpha',
        'critDamage', 'voidStiffnessRatio' and 'density'
    """
    # copying all entries of provided material to new MatItem
    self._parent = kwargs.pop('parent', None)
    self.name = kwargs.pop('name', '')

    # update metaData
    self._frozen = False
    creationStr = kwargs.pop('creationStr', 'new')
    self._meta = MetaData(creationStr = creationStr)


    for key,val in kwargs.iteritems():
        setattr(self, key, val)

    # used by beRocker
    self._plotState = dict(general = False)
    for t in range(len(self)):
        self._plotState[t] = True

    super(pLR.__class__, self).__init__(pLR)



### different MatItemClasses
class MatItemPLR(MatItem, PlasticLevkovitchReusch):
    """MatItem of material type L{material.PlasticLevkovitchReusch}.

    Will be dynamically created by L{dynamicMatItem}.
    """
    _matType = PlasticLevkovitchReusch

    def __init__(self, pLR, **kwargs):
        itemInit(self, pLR, **kwargs)


class MatItemPLRAniso(MatItem, PlasticLevkovitchReuschAniso):
    """MatItem of material type L{material.PlasticLevkovitchReuschAniso}.

    Will be dynamically created by L{dynamicMatItem}.
    """
    _matType = PlasticLevkovitchReuschAniso

    def __init__(self, pLR, **kwargs):
        itemInit(self, pLR, **kwargs)


class MatItemSPMC(MatItem, SinglePlaneMohrCoulomb):
    """MatItem of material type L{material.SinglePlaneMohrCoulomb}.

    Will be dynamically created by L{dynamicMatItem}.
    """
    _matType = SinglePlaneMohrCoulomb

    def __init__(self, sPMC, **kwargs):
        itemInit(self, sPMC, **kwargs)


def dynamicMatItem(pLR, **kwargs):
    """This is used to create a MatItems dynamically by passing the
    MaterialObject and its (additional or required) keyword arguments
    """
    validMatTypesPLRAniso =  ['PlasticLevkovitchReuschAniso',
                              'MatItemPLRAniso']
    validMatTypesPLRIso =  ['PlasticLevkovitchReusch',
                            'MatItemPLR']
    validMatTypesSPMC = ['SinglePlaneMohrCoulomb','MatItemSPMC']

    if pLR.__class__.__name__ in validMatTypesPLRAniso:
        return MatItemPLRAniso(pLR, **kwargs)

    elif pLR.__class__.__name__ in validMatTypesSPMC:
        return MatItemSPMC(pLR, **kwargs)

    elif pLR.__class__.__name__ in validMatTypesPLRIso:
        return MatItemPLR(pLR, **kwargs)

    else:
        print(pLR.__class__.__name__)
        raise TypeError("Can't create a MatItem-object from an object" +
                        " of type %s" % str(type(pLR)) )

##############################################################################
#}#

##############################################################################
#{Collection classes
class MatItemCollection(OrderedDict):
    """Objects of this class describe one kind of material.

    @ivar name: name as provided on initialization. Note: In contrast the
      property L{fullName} provides this name with a prefix identifying
      the L{MatPropsCollection} where it comes from.
    """
    def __init__(self, name, creationStr='new', parent=None, **kwargs):

        super(MatItemCollection, self).__init__(**kwargs)
        self.name = name
        self._meta = MetaData(creationStr=creationStr)
        self._parent = parent
        self._frozen = False


    def __str__(self):
        """Prints a formated string with name and items of MatItemCollection
        """
        rtn = "MatItemCollection %s\n" % self.fullName
        rtn += "\tholding the following %d matItems:\n\t" % len(self)
        rtn += "\n\t".join(self.keys())
        return rtn


    def toDict(self):
        """Stores a more or less complete representation of the MatpropsCollection
        to a dictionary-structure.
        It will hold (parameter)-dictionarys for rock- and fill-MatItemLists and
        for the fillAssignments. Furthermore some/most of the metaData gets stored.

        Example:
         >>> import json
         >>> myCollectionDict = myCollection.toDict()
         >>> fid = open('myCollection.json', 'w')
         >>> json.dump(myCollection, fid)
         >>> # a few weeks later:
         >>> loadedCollectionDict = json.load(open('myCollection.json','r'))
         >>> myRestoredCollection = MatpropsCollection.fromDict(loadedCollectionDict)

        @note: This method is intended to give a mere "simple python" representation
            of the MatpropsCollection that can easily processed by pickle, json, ...
            It's not ment to access and manipulate the stored data. If you do so,
            keep in mind that changing the returned dictionary-structure might
            break its compatibility with L{fromDict}-function.
        """
        return dict(name=self.name,
                    metaData=self._meta.toDict(),
                    frozen=self._frozen,
                    items=[(key,item.toDict()) for key,item in self.iteritems()])


    @classmethod
    def fromDict(cls, input):
        """Creates a MatpropsCollection-instance from a dictionary structure.


        """
        itemCollection = cls(input['name'])
        itemCollection._meta = MetaData.fromDict(input['metaData'])
        itemCollection._frozen = input['frozen']

        for key, itemDict in input['items']:
            itemDict['parent'] = itemCollection
            itemCollection[key] = MatItem.fromDict(itemDict)

        return itemCollection


    def freeze(self, doIt=True):
        """Freezes a MatItemCollection and all its children so it can't be
        edited using matpropscollection-fuctions anymore.
        Use self.freeze(False) to unfreeze.
        """
        self._frozen = doIt

        for mat in self.values():
            mat.freeze(doIt)
        if doIt:
            self._meta.log('MaterialItemCollection got frozen', level=0)
        else:
            self._meta.log('MaterialItemCollection got unfrozen', level=0)


    @property
    def fullName(self):
        """Returns the complete name of a MatItemCollection e.g. 'M01.fill'
        """
        try:
            return self._parent.name + '.' + self.name
        except Exception:
            return self.name


    def getItem(self, name):
        """Returns the requested MatItem. Is identical to
        MatItemCollection.__getitem__"""
        try:
            return self[name]
        except KeyError:
            raise KeyError(('MaterialItemCollection %s has no material' +
                            ' named %s') % (self.name, name))


    def filterAndSort(self, filterStr = '', order = 'name'):
        """Returns a (possibly filtered and) and sorted list of MatItem-keys

        Usage:
            >>> sortedKeys,countOther = materials.filterAndSort(
            ...                            filterStr='MATVOL',
            ...                            order='UCSi')

        returns: filtered/sorted MatItem-keys, count of filtered/removed keys
        """
        import re

        if hasattr(order, '__iter__'):
            order, prop = order[0], order[1]
        try:
            filt = re.compile(filterStr, re.IGNORECASE)
            rtnStr = [ff for ff in self.keys() if filt.search(ff)]
        except:
            rtnStr = []

        invCnt = len(self.keys()) - len(rtnStr)
        if 'name' in order.lower():
            rtnStr = sorted(rtnStr)

        elif 'property' in order.lower() :
            found, notFound = [], []
            for name in rtnStr:
                try:
                    found.append( (self[name].getProperty(prop), name) )
                except:
                    notFound.append(name)

            if notFound:
                msg("Couldn't apply sorting to all items.")
            found = [name for val, name in sorted(found)]
            rtnStr = found + notFound
        else:
            raise ValueError("Don't understand order-type %s" % str(order))

        if 'reverse' in order.lower():
            rtnStr = rtnStr[::-1]

        return rtnStr, invCnt

    def setItemParameter(self, name, prop, val, **kwargs):
        """Sets a material property for requested MatItem.
        See L{MatItem.setProperty} as well.
        @param name: Name of MatItem
        @param prop: property name
        @param val: new property value
        @kwarg idx: index of requested pst-sample.
        """
        item = self.getItem(name)

        oldLogState = item._meta.logState
        noLog = kwargs.pop('disableLogging', False)
        item._meta.logState = not noLog

        item.setProperty(self, prop, val, **kwargs)

        item._meta.logState = oldLogState


    def getItemParameter(self, name, prop, **kwargs):
        """Returns a material property for requested MatItem.
        See L{MatItem.getProperty} as well.
        @param name: Name of MatItem
        @param prop: property name
        @kwarg idx: index of requested pst-sample.
        """
        material = self.getItem(name)
        return material.getProperty(prop, **kwargs)


    def getItemPropertyDict(self, name, idx=None):
        """ Returns property-dictionaries for requested MatItem.
        See L{MatItem.getPropertyDicts} as well.
        @param name: Name of MatItem
        @param idx: index of requested pst-sample. If None, only general
             properties will be returned
        @returns generalPropsDict, samplePropsDict
        """
        material = self.getItem(name)
        return material.getPropertyDicts(idx=idx)


    def removeItem(self, name, logText='removed'):
        """Removes a MatItem from MatItemCollection.
        @param name: name of matItem that will be removed
        @param logText: descriptive string for log-events
        @note: this function will remove all references of MatItem in
            Assignments- structure of a parent-MatPropsCollection as well
        """
        if self._frozen:
            msg('MaterialItemCollection is frozen. No chages allowed.')
            return False

        if name in self.keys():
            self.pop(name, None)
            self._meta.log('%s %s ' % (logText, name))
            if self._parent is not None:
                self._parent.fillAssignments.removeMat(name, self.name)
            return True

        else:
            msg("Couldn't remove material %s" % name +
                "because it's not in MaterialItemCollection %s" % self.name)
            return False

    def _insertItem(self, name, item, force=False):
        """Support fuction that inserts an item to MatItemColltion but dosen't
            trigger alog-entry. It is is used by multible fuctions.
        """
        if self._frozen:
            msg('MaterialItemCollection is frozen. No changes allowed.')
            return False

        if name in self.keys() and not force:
            msg("Material named %s already exists in " % name +
                "MaterialItemCollection %s. Skipping insertion" % self.name)
            return False
        else:
            item.name = name
            self[name] = item
            return True


    def addItem(self, name, item, **kwargs):
        """Adds a copy of a (new) MatItem to MatItemCollection.
        @param name: name of new item in MatItemCollection
        @param item: item to add. Can be a MatItem or a
            L{bae.material_01.material.DefaultPst} (or derived PST-sample list).
        @kwarg logText: default='added', descriptive string for log-events
        @kwarg force: default=False, if True, an existent matItem with same
            name will be overwritten by the new item
        """
        logText = kwargs.pop('logText','added')
        force   = kwargs.pop('force', False)
        if not isinstance(item, (MatItemPLR, MatItemSPMC, MatItemPLRAniso)):
            kwargs['name'] = name
            item = dynamicMatItem(item, **kwargs)

        valid = self._insertItem(name, item.copy(), force)

        if valid:
            item._meta.log('%s to %s' % (logText, self.fullName))
            self._meta.log('%s %s' % (logText,name))
        return valid


    def copyItem(self, name, newName, logText='copied', force=False):
        """Adds a (deep)copy of the requested matItem to MatItemCollection.
        @param name: name of matItem that will be copied
        @param newName: name of copy in MatItemCollection
        @param logText: descriptive string for log-events
        @param force: if True, an existent matItem with same name as newName
            will be overwriten by the copy
        """
        item = self.getItem(name).copy()

        if item._frozen:
            item._meta.clearExceptRootLevel()

        passed = self._insertItem(newName, item, force=force)
        if passed:
            item._meta.log('%s from %s' % (logText, name), level=0)
            self._meta.log('%s %s to %s' % (logText, name, newName))
        return passed


    def copyItemToOtherCollection(self, name, otherCollection, newName=None,
                                  logText='copied', force=False):
        """Copies the requested matItem to different MatItemCollection.
        @param name: name of matItem that will be copied
        @param otherCollection: MatItemCollection to copy MatItem to
        @param newName: name of copy in the other MatItemCollection
        @param logText: descriptive string for log-events
        @param force: if True, an existent matItem with same Name as newName
            will be overwriten by the copy
        """
        item = self.getItem(name).copy()

        if item._frozen:
            item._meta.clearExceptRootLevel()
            item.freeze(False)

        if newName is None:
            newName = name

        item._meta.log('Copied from %s' % self.fullName, level=0)
        passed = otherCollection._insertItem(newName, item, force=force)

        if passed:
            otherCollection._meta.log('%s %s from collection %s to %s' %
                                      (logText, name, self.fullName,
                                       otherCollection.fullName))
        return passed


    def moveItemToOtherCollection(self, name, otherCollection, force=False):
        """Moves the requested matItem to different MatItemCollection.
        @param name: name of matItem that will be copied
        @param otherCollection: MatItemCollection to move MatItem to
        @param logText: descriptive string for log-events
        @param force: if True, an existent matItem with same Name as newName
            will be overwriten by the copy
        @note: this function will remove all references of MatItem in
            Assignments-structure of a parent-MatPropsCollection as well
        """
        if self._frozen:
            msg('MaterialItemCollection is frozen. No chages allowed.')
            return False

        passed = self.copyItemToOtherCollection(name, otherCollection,
                                                force=force, logText='moved')
        if passed:
            self.removeItem(name)
            self._meta.log('moved %s to %s' %
                           (name, otherCollection.fullName) )

        return passed


    def renameItem(self, oldName, newName, logText = 'renamed', force=False):
        """Renames a MatItem in MatItemCollection.
        @param name: name of matItem that will be renamed
        @param newName: new name of MatItem in MatItemCollection
        @param logText: descriptive string for log-events
        @param force: if True, an existent MatItem with same Name as newName
            will be overwriten by the renamed MatItem
        @note: this function will rename all references of MatItem in
            Assignments-structure of a parent-MatPropsCollection as well
        """
        if self._frozen:
            msg('MaterialItemCollection is frozen. No chages allowed.')
            return False

        if newName in self.keys() and not force:
            msg("Material named %s allready exists in " % newName +
                "MaterialItemCollection %s. Skipping rename" % self.name)
            return False
        else:
            item = self.pop(oldName)
            item.name = newName
            self[newName] = item
            item._meta.log('%s from %s to %s'
                           % (logText, oldName, newName))
            self._meta.log('%s from %s to %s'
                           % (logText, oldName, newName))
            if self._parent is not None:
                a = self._parent.fillAssignments
                a.renameMat(oldName, self.name, newName,
                            logText=logText)
            return True


    def importVusubsVarFill(self, vusubsParamPy, force=False,
                            offset=1000, namePattern='varFill_%d'):
        """Imports a varFill materials from vusubs_param.py (state March 2021).
        The vusubs_param.py has to hold the following keywords:
            - useModule_stateCtrl (should be 'varFill'),
            - useAniso,
            - maxNumberFill,
            - fillPropsStartIdx and
            - PROPS_fill,
        where PROPS_fill has the following structure for M different varFill
        materials:
            isotropic:
                M times [UCS, e, 1/a, N, (pst, G, K, s, mb, dilation,),]
            anisotropic:
                M times [UCS, e, 1/a, anisoS, anisoN, N, (pst, G, K, s,
                mb, dilation,),] .

        @param vusubsParamPy: path of vusubs_param.py
        @param force: if True, existing materials in collection having the same
            name get overwriten.
        @param offset: fillTime offset for varFill. Used for naming.
        @param namePattern: namePattern for newly created fill material
        @returns: self
        """
        ### Needs to be improved ironed in future
        content = {}
        extraLocals = {'Var': lambda x: None, 'FloatExpr': lambda x: None}
        execfile(vusubsParamPy, content, extraLocals)
        ###
        try:
            stateCtrl = extraLocals['useModule_stateCtrl']
            useAniso = extraLocals['useAniso']
            maxNumberFill = extraLocals['maxNumberFill']
            fillPropsStartIdx = extraLocals['fillPropsStartIdx']
            PROPS_fill = extraLocals['PROPS_fill']
        except KeyError as e:
            import sys
            message = e.message
            message += ' Badly formated or outdated vusubs_param.py.'
            raise type(e), type(e)(message, sys.exc_info()[2])

        if not stateCtrl == "varFill":
            msg('WARNING! useModule_stateCtrl is "%s" instead of "varFill"'
                % str(stateCtrl))

        fillParams = OrderedDict()
        if maxNumberFill == 1:
            fillParams[namePattern%offset] = [PROPS_fill,]
        else:
            iLower = 0
            fillPropsStartIdx += [len(PROPS_fill) + 1,]
            for ii, iUpper in enumerate(fillPropsStartIdx[1:], start=1):
                iUpper -= 1 # Fortran to python
                name = namePattern % (ii*offset)
                fillParams[name] = PROPS_fill[iLower:iUpper]
                iLower = iUpper
            if not ii == maxNumberFill:
                msg(("WARNING! Counted %d instead of maxNumberFill=%d "
                     "fill items") % (ii, maxNumberFill))

        if useAniso:
            material = PlasticLevkovitchReuschAniso
            sample = material._sampleType
            fixedArgs = ['UCSi', 'e', 'a', 'anisoS', 'ansioN']
        else:
            material = PlasticLevkovitchReusch
            fixedArgs = ['UCSi', 'e', 'a']

        def parser(p):
            ysArgs = dict((name, p.pop(0)) for name in fixedArgs)
            ysArgs['a'] =1./float(ysArgs['a'])
            nSamples = int(p.pop(0))
            matArgs = []

            for ii in range(nSamples):
                pst = p.pop(0)
                G = p.pop(0)
                K = p.pop(0)
                ysArgs['s'] = p.pop(0)
                ysArgs['mb'] = p.pop(0)
                dilation = p.pop(0)
                elastic = ElasticLinear(G=G, K=K)
                ys = sample._yieldType(**ysArgs)
                matArgs.append([pst, elastic, ys, dilation])

            return material(matArgs)

        for name, params in fillParams.iteritems():
            creationStr = ('initially loaded from %s' % vusubsParamPy)
            fill = parser(params)
            mat = dynamicMatItem(fill, parent = self, name = name,
                                 creationStr=creationStr)
            self._insertItem(name, mat, force=force)

        self._meta.log('imported materials from %s' % vusubsParamPy,
                       names = fillParams.keys())

        return self


    def importAbaqusInput(self, inputFiles, force=False):
        """Depending on self.name (can be rock or fill) the material cards for
        rock or fill will be imported from abaqus.inp files.
        @param inpNames: can be a abaqus.inp or list of files
        @param force: if True, existing materials in collection having the same
            name get overwriten.
        @returns: self
        """
        from bae.abq_model_02 import Model
        from bae.material_01 import AbqVumatData

        if self._frozen:
            msg('MaterialItemCollection is frozen. No chages allowed.')
            return False

        if isinstance(inputFiles, basestring):
            inputFiles = [inputFiles,]

        m = Model()
        matKeys = [] # is used to preserve order from inputDeck
        def getName(inpLine):
            for string in ['*MATERIAL',',','NAME=']:
                inpLine = inpLine.replace(string, '').strip()
            return inpLine

        for name in inputFiles:
            m.read(name)
            try:
                ff = open(name,'r')
            except IOError:
                try: #assuming name is StringIO
                    ff = name
                    ff.seek(0)
                except AttributeError:
                    ff = ff.split('\n') #assuming name is text

            keys = [getName(line)
                    for line in ff if line.startswith('*MATERIAL')]
            allReadyExist = [key for key in keys if key in matKeys]
            if allReadyExist:
                raise ValueError(("The following matNames from file %s"
                                  " were already found in a previous file: %s")
                                  % (name, allReadyExist))
            matKeys.extend(keys)

        inCollection = set(self.keys()).intersection(
            set( m.material.keys()) )

        if not force and inCollection:
            msg("Some materials are allready exist in Collection %s:\n" %
                (self.name) + '\n'.join(inCollection) +
                '\n No Material will be loaded from input file. ' +
                ' To overwrite existing use force=True')
            return

        addPropKeys = ["dampingAlpha", "density", "voidStiffnessRatio",
                       "critDamage", "umatDelete", "umatDepvarNum"]

        for name in matKeys:
            mat = m.material[name]
            abq = AbqVumatData(mat)
            if self.name.lower() == 'rock':
                data = abq.plasticLRRock
            elif self.name.lower() == 'fill':
                data = abq.plasticLRFill
            else:
                raise ValueError('Name of MatItemCollection has to be '+
                                 'either "rock" or "fill".')

            addProps = dict( (p, getattr(abq, p)) for p in addPropKeys )

            inpBaseName =  os.path.split(name)[-1]
            creationStr = ('initially loaded from %s' % inpBaseName)

            mat = dynamicMatItem(data, parent = self, name = name,
                          creationStr=creationStr,
                          **addProps)

            self._insertItem(name, mat, force=force)

        self._meta.log('imported materials from %s' %
                       [os.path.basename(n) for n in inputFiles],
                       names = m.material.keys())

        return self


    def analyzeDuplicates(self):
        """Compares the properties of all stored matItems.
        """
        data = OrderedDict( **self )

        allDuplicates = {}
        materialKeys = data.keys()

        for mat in materialKeys:
            try:
                duplicates = []
                master = data[mat]
                for name, matC in data.items():
                    isEqual, delta = master.compare(matC)
                    if isEqual:
                        duplicates.append(name)
                        data.pop(name)
                if len(duplicates) > 1:
                    allDuplicates[mat + '_like'] = duplicates
            except KeyError:
                pass

        return allDuplicates


    def _excelSummary(self, workBook, matNames=None, pstSamples='all'):
        from plot import (pstSampleNames, _getSampleNamesAndIndices,
                          _setExcelFormats, _setExcelSheetTable,
                          descriptiveCellData)
        from collections import defaultdict

        byMatType = defaultdict(list)
        #byMatType.default_factory = None
        for name in matNames:
            matItem = self[name]
            pstS, _ = _getSampleNamesAndIndices(matItem, pstSamples)
            sampleNames = pstSampleNames(matItem)
            _, descriptiveMatType = characteristics(matItem)
            descriptiveMatType = descriptiveMatType.split()[0]
            descriptiveMatType = descriptiveMatType.replace('Computator','PLR')
            descriptiveMatType = descriptiveMatType.replace('arbitrary','PLR')
            cA, cB = descriptiveCellData(matItem, pstSamples=pstS,
                                         latex=False)
            byMatType[descriptiveMatType].append([name, cA, cB, sampleNames])

        if not byMatType:
            return

        sheetName = '_summary_%s' % self.name
        sheet = workBook.add_worksheet(sheetName)

        formats = _setExcelFormats(workBook)

        def puzzleCellData(name, cA, cB, sNames):
            if cA:
                dataA = [ct[1:] for ct in cA]
                dataA = [list(ll) for ll in izip(*dataA)]
                labelsA = [ct[0] for ct in cA]
            else:
                dataA, labelsA = [[],], []

            if cB:
                dataB = [ct[1:] for ct in cB]
                dataB = [list(ll) for ll in izip(*dataB)]
                labelsB = [ct[0] for ct in cB]
            else:
                dataB, labelsB = [[],], []

            if labelsB:
                labels = ['Name',] + labelsA + ['Level',] + labelsB
            else:
                labels = ['Name',] + labelsA

            nRows = max(len(dataA), len(dataB))

            cellData = []
            nameCol = [name,] + ['' for _ in range(nRows-1)]
            for ii in range(nRows):
                try:
                    dA = dataA[ii]
                except IndexError: # pad with empty cells
                    dA = ['' for _ in range(len(dataA[0]))]

                try:
                    dB = dataB[ii]
                except IndexError: # pad with empty cells
                    dB = ['' for _ in range(len(dataB[0]))]

                if dB:
                    cellData.append([nameCol[ii],] + dA + [sNames[ii],] + dB)
                else:
                    cellData.append([nameCol[ii],] + dA)

            return labels, cellData, nRows

        row = 1
        for mType, heading in [('PLR', 'LR2-Materials'),
                               ('MOHR-COULOMB', 'MOHR-COULOMB-like Materials'),
                               ('elastic', 'Purely Elastic Materials'), ]:
            try:
                cellData = byMatType[mType]
                cellData[0]
            except (KeyError, IndexError):
                continue

            _data = [puzzleCellData(*cs) for cs in cellData]
            cells = [c for _,cc,_ in _data for c in cc]
            labels = _data[0][0]
            sepLines = [True if jj+1==matSamples else False
                        for _,_,matSamples in _data
                        for jj in range(matSamples)]

            row,_ =_setExcelSheetTable(sheet, cells, formats,
                                   row=row, col=0, sepLines = sepLines,
                                   rowLabels=None, columnLabels=labels,
                                   title=heading)
            row += 2


    def toXlsx(self, workBook, matNames=None, **kwargs):
        """Creates an excel file/workBook holding formated properties and
        plots of requested materials in MatItemCollection. Each material gets
        stored on a new sheet.

        @params workBook: excel-fileName or xlsxwriter.Workbook-instance
        @params matNames: list of rockNames to be exported. Use None to
            export all materials stored in MatItemCollection
        @note: see L{MatItem.toXlsx} for more keyword-arguments to adjust
            output for your needs.
        """
        try:
            import xlsxwriter as xls
        except ImportError:
            raise ImportError("xlsxwriter is not installed. Use:\n" +
                              ">>> pip install XlsxWriter\n" +
                              "on regular python or\n" +
                              ">>> conda install -c anaconda xlsxwriter\n"
                              "if you are using anaconda")

        if not isinstance(workBook, xls.Workbook):
            workBook = xls.Workbook(workBook, {'nan_inf_to_errors': True})
            closeIfDone = True
        else:
            closeIfDone = False

        if matNames is None:
            matNames = self.keys()

        notFoundNames = [n for n in matNames if n not in self.keys()]
        if notFoundNames:
            raise KeyError("The following materials can't be found in " +
                           "the matpropsCollecion: %s" %
                           ",".join(notFoundNames))
        try:
            self._excelSummary(workBook, matNames)
        except Exception as e:
            print(e)
            msg("ERROR! Failed to write summary sheet for %s. Will continue."
                % self.name)

        for name in matNames:
            msg('Writing %s to %s' % (name, workBook.filename))
            mat = self[name]
            try:
                mat.toXlsx(workBook, **kwargs)
            except Exception as e:
                msg('ERROR!!! While writing %s to %s, the following Error '
                    'occoured: "%s". Will contiue anyway' %
                    (name, workBook.filename, str(e)))
                raise

        if closeIfDone:
            workBook.close()


    def plotYieldSurfs(self, outputDirName=None, fileNamePrefix=None,
                       interactive=False, matNames=None,
                       propertyTables=False, pstSamples=None, postPlotFun=None,
                       **additionalPlotArgs):
        """
        Plots yield surfaces for a set of materials and saves these plots. For each matItem
        the L{MatItem.plotAniso} function is called.

        Example Some default plots:
            >>> myCollection.plotAniso(fileNamePrefix='PRJ2020999_MATPROPS_FAULTS_M99',
            ...                        matNames=['FAULT_A', 'FAULT_B'])

        Example Change y-limits of left axis:
            >>> def myPostPlotFun(fig, axs, plts, legs, mat):
            >>>     axs[0].set_ylim(50, 100)
            >>> myCollection.plotAniso(fileNamePrefix='PRJ2020999_MATPROPS_M99',
            ...                        postPlotFun=myPostPlotFun)

        @param outputDirName: folder where plots get stored. If left None, the name will be
            choosen according to the name of MatPropsCollection.
        @param fileNamePrefix: fileName will be C{fileNamePrefix + '_' + matName + '.png'}
        @param interactive: if you have a pysical screen (not like on nukus) set this True
            to get the plots displayed as matplotlib-windows
        @param matNames: series of material names stored in collection. If None, all will
            be ploted.
        @param propertyTables: if True, the properties of a matItem will appear in tabulars
            at the righthand side of the plots
        @param pstSamples: specifies the pst-samples to be documented (e.g.
            C{['PEAK','RES']}). If None, all samples will be plotted.
        @param postPlotFun: a function that takes figure, axes, plots and legend handles and
            the materialItem to manipulate the plot before it gets stored. See example above.
        @param additionalPlotArgs: dictionary of additional plot-settings that will be passeed
            to L{MatItem.createYieldSurfPlot}.

        """
        if matNames is None:
            matNames = self.keys()

        notFoundNames = [n for n in matNames if n not in self.keys()]
        if notFoundNames:
            raise KeyError("The following materials can't be found in " +
                           "the matItemCollecion: %s" %
                           ",".join(notFoundNames))

        figs = []
        for name in matNames:
            # create output file name
            if fileNamePrefix:
                if outputDirName is None:
                    try:
                        outputDirName = self._parent._parent.name
                    except:
                        pass
                    if not outputDirName:
                        outputDirName = './'
                saveName = os.path.join(outputDirName,
                                        fileNamePrefix + '_' + name + '.png')
            else:
                saveName = None

            if saveName:
                msg('saving yieldsurface plots of %s to %s'
                    % (name, saveName))
            else:
                msg('plotting yieldsurfaces of %s' % name)

            fig = self[name].createYieldSurfPlot(
                              saveName=saveName,
                              interactive=interactive,
                              pstSamples=pstSamples,
                              propertyTables=propertyTables,
                              postPlotFun=postPlotFun,
                              **additionalPlotArgs)
            figs.append(fig)
        return figs


    def plotAniso(self, outputDirName=None, fileNamePrefix=None,
                       interactive=False, matNames=None, postPlotFun=None,
                       **additionalPlotArgs):
        """
        Plots anisotropic behavior for a set of materials and saves these plots. For
        each matItem the L{MatItem.plotAniso} function and thus finally L{plot.plotAniso}
        is called.

        Example Some default plots:
            >>> myCollection.plotAniso(fileNamePrefix='PRJ2020999_MATPROPS_FAULTS_M99',
            ...                        matNames=['FAULT_A', 'FAULT_B'])

        Example Change y-limits of left axis:
            >>> def myPostPlotFun(fig, axs, plts, legs, mat):
            >>>     axs[0].set_ylim(50, 100)
            >>> myCollection.plotAniso(fileNamePrefix='PRJ2020999_MATPROPS_M99',
            ...                        postPlotFun=myPostPlotFun)

        @param outputDirName: folder where plots get stored. If left None, the name will be
            choosen according to the name of MatPropsCollection.
        @param fileNamePrefix: fileName will be C{fileNamePrefix + '_' + matName + '.png'}
        @param interactive: if you have a pysical screen (not like on nukus) set this True
            to get the plots displayed as matplotlib-windows
        @param matNames: series of material names stored in collection. If None, all will
            be ploted.
        @param postPlotFun: a function that takes figure, axes, plots and legend handles and
            the materialItem to manipulate the plot before it gets stored. See example above.
        @param additionalPlotArgs: dictionary of additional plot-settings that will be passeed
            to L{MatItem.plotAniso}.
        """
        # import matplotlibs plt, pylab and gridspec via plot.py
        from plot import plt, pylab, plotParams
        if matNames is None:
            matNames = self.keys()

        notFoundNames = [n for n in matNames if n not in self.keys()]
        if notFoundNames:
            raise KeyError("The following materials can't be found in " +
                           "the matItemCollecion: %s" %
                           ",".join(notFoundNames))

        figs = []
        for name in matNames:
            mat = self[name]
            if not mat.isAnisotropic:
                msg('skipping anisotropic plots for isotropic material %s'
                    % mat.fullName)
                continue

            # create output file name
            if fileNamePrefix:
                if outputDirName is None:
                    try:
                        outputDirName = self._parent._parent.name
                    except:
                        pass
                    if not outputDirName:
                        outputDirName = './'
                saveName = os.path.join(outputDirName,
                                        fileNamePrefix + '_' + name + '.png')
            else:
                saveName = None

            if saveName:
                msg('saving anisotropy plots of %s to %s'
                    % (name, saveName))
            else:
                msg('plotting yieldsurfaces of %s' % name)

            pylab.rcParams.update(plotParams)
            fig = plt.figure(figsize=(17.34/2.54, 11.00/2.54), dpi=150,)
            ax1 = fig.add_subplot(1,2,1)
            ax2 = fig.add_subplot(1,2,2)
            title = (name +
                     r" -- ($s_{aniso}=%.2f$, $n_{aniso}=%.2f$)"
                     % (mat.anisoS, mat.anisoN))
            fig.suptitle(title)

            handles = mat.plotAniso(axh1=ax1, axh2=ax2,
                          overwriteAxesProps=[True, True],
                          **additionalPlotArgs)

            if postPlotFun is not None:
                _, axHandles, pltHandles, legendHandles = handles
                postPlotFun(fig, axHandles, pltHandles, legendHandles, mat)

            if saveName:
                outDir = os.path.dirname(saveName)
                if not os.path.isdir(outDir):
                    os.makedirs(outDir)
                fig.savefig(saveName, bbox_inches='tight')

            if interactive:
                plt.show()

            figs.append(fig)
        return figs

###############################################################################

class Assignment(object):
    """Describes a single rock-to-fill assignment that will be used if no
    varfill=True is used.
    """
    def __init__(self, rockName, fillName, parent=None, name=None,
                 master='rock', creationStr = 'new'):
        """
        @param rockName: name of rock-materialItem in parent.rock
        @param fillName: name of fill-materialItem in parent.fill
        @param parent: MatPropsCollection
        @param name: name of assignment. This will finally be used as key in
            MatItemAssignments and as material name in exported abaqus input
            deck. If None or '' is passed, the name will be created joining
            rock- and fillName or only rockName if it equals the fillName
        @param master: defines if common properties (density, damping,...)
            is taken from 'rock' or 'fill'
        @param creationStr: first (creation) log-event of this item
        """
        self.rock = rockName
        self.fill = fillName
        self._parent = parent

        if name:
            self.name = name
        else:
            if rockName == fillName:
                self.name = rockName
            else:
                self.name = rockName + '_' + fillName

        self.master = master

        # metadata
        self._frozen = False
        self._meta = MetaData(creationStr = creationStr)


    def toDict(self):
        return dict(name=self.name,
                    rockName=self.rock, fillName=self.fill,
                    frozen=self._frozen, master=self.master,
                    metaData=self._meta.toDict())


    @classmethod
    def fromDict(cls, input):
        newAssignment = cls(input['rockName'], input['fillName'],
                            name=input['name'], master=input['master'])

        newAssignment._frozen = input['frozen']
        newAssignment._meta = MetaData.fromDict(input['metaData'])

        return newAssignment

    def freeze(self, doIt=True):
        """Freezes the Assignment so that they can't be edited using
        matpropscollection-fuctions anymore. Use self.freeze(False) to unfreeze.
        """
        self._frozen = doIt
        if doIt:
            self._meta.log('assignment %s got frozen' % self.name)
        else:
            self._meta.log('assignment %s got unfrozen' % self.name)


    def getAssignedMaterials(self):
        """Returns the assigned rock and fill material from parent
        MatPropsCollection.
        """
        try:
            rock = self._parent._parent.rock[self.rock]
            fill = self._parent._parent.fill[self.fill]

        except AttributeError as e:
            raise Exception("Can't get material-items for assignment %s. " %
                            self.name + 'The following Error occured:\n' +
                            str(e))

        return rock, fill

class MatDomainAssignments(dict):
    isCohesiveDomain = lambda s: s.startswith('COHS_')

    def __init__(self, creationStr, parent=None, **kwargs):
        """Creates a new MatItemAssignments object.
        @param parent: set a MatProsCollection as parent if available
        @param creationStr:
        """
        super(MatItemFillAssignments, self).__init__(**kwargs)

        self._parent = parent
        self._meta = MetaData(creationStr = creationStr)
        self._frozen = False

    @property
    def hasAssignments(self):
        return any([len(v) for v in self.values()])


    def addAssignment(self, matItemFillAssignment, domains):
        self[matItemFillAssignment] = domains


class MatItemFillAssignments(OrderedDict):
    """A (ordered) dictionary that stores and manages single rock-to-fill-
    L{Assignment}s.
    """
    def __init__(self, creationStr, parent=None, **kwargs):
        """Creates a new MatItemAssignments object.
        @param parent: set a MatProsCollection as parent if available
        @param creationStr:
        """
        super(MatItemFillAssignments, self).__init__(**kwargs)

        self._parent = parent
        self._meta = MetaData(creationStr = creationStr)
        self._frozen = False


    def toDict(self):
        return dict(assignments=[(key, assignment.toDict())
                                 for key,assignment in self.iteritems()],
                    frozen=self._frozen,
                    metaData=self._meta.toDict())


    @classmethod
    def fromDict(cls, input):
        newAssignments = cls('')
        for key, assignment in input['assignments']:
            newAssignments[key] = Assignment.fromDict(assignment)
        newAssignments._frozen = input['frozen']
        newAssignments._meta = MetaData.fromDict(input['metaData'])

        return newAssignments


    def freeze(self, doIt=True):
        """Freezes the MatItemAssignments and all its childrens (assignments)
        so that they can't be edited using matpropscollection-fuctions anymore.
        Use self.freeze(False) to unfreeze.
        """
        for ass in self.values():
            ass.freeze(doIt)
        if doIt:
            self._meta.log('assignments got frozen')
        else:
            self._meta.log('assignments got unfrozen')


    def addAssignment(self, rockName, fillName,  name=None, force=False,
                 master='rock', creationStr = 'new'):
        """Adds an AssignmentItem if MatItemAssignments is not froozen yet.
        @param rockName: name of rock-materialItem in parent.rock
        @param fillName: name of fill-materialItem in parent.fill
        @param name: This will finally be used as key in
            MatItemAssignments and as material name in exported abaqus input
            deck. If None or '' is passed, the name will be created joining
            rock- and fillName or only rockName if it equals the fillName
        @param master: defines if common properties (density, damping,...)
            is taken from 'rock' or 'fill' for this item
        @param creationStr: first (creation) log-event of this item
        """

        if self._frozen:
            msg('MatItemAssignments is frozen. No chages allowed.')
            return False

        item = Assignment(rockName, fillName, parent=self, name=name,
                          master=master, creationStr = creationStr)
        if not force and item.name in self.keys():
            msg('Assignment %s allready exists.' % item.name +
                ' Use force=True to overwrite')
            return False

        rock, fill = item.getAssignedMaterials()

        if not rock._matType == fill._matType:
            msg('WARNING! %s assigns %s for rock and %s for fill.\n' %
                (item.name, rock._matType.__name__, fill._matType.__name__) +
                'Can your vUmat handle this???')

        self[item.name] = item
        self._meta.log('added %s' % item.name)
        return True


    def renameAssignment(self, name, newName, logText='renamed', force=True):
        """Renames an assignment that is not frozen yet.

            >>> mpc.fillAssignments.renameAssignment('HOST', 'FANCY')

        @param name: old name of assignment
        @param newName: new name of assignment
        @param force: if True, an existant assignment named newName will be
            overwritten
        @param logText: textPrefix that appears before '$name to $newName' in
            log-event of MatItemAssignments (self._meta.logs)
        """
        if self._frozen:
            msg('MatItemAssignments is frozen. No chages allowed.')
            return False

        if newName in self.keys() and not force:
            msg('Assignment %s allready exists.' % newName +
                ' Use force=True to overwrite')
            return False

        if name in self.keys():
            self[newName] = self.pop(name)
            self._meta.log('%s %s to %s'
                           % (logText, name, newName))
            return True

        else:
            msg('No assignment named %s. Renameing skipped.' % name)
            return False


    def removeAssignment(self, name, logText='removed'):
        """Removes an assignment that is not frozen yet.

        @param name: name of assignment that will be removed
        @param logText: textPrefix that appears before '$name' in
            log-event of MatItemAssignments (self._meta.logs)
        """
        if self._frozen:
            msg('MatItemAssignments is frozen. No chages allowed.')
            return False

        if name in self.keys():
            self.pop(name)
            self._meta.log('%s %s' % (logText, name))
            return True
        else:
            msg('No assignment named %s. Removal skipped.' % name)
            return False


    def removeMat(self, name, matType):
        """Removes all assignments holding the requested material(Name).
        @param name: materialName that will be removed from assignments
        @param matType: the materialName referes to "rock" or "fill"
        @note: this function is usually called from parent-MatPropsCollection
            to remove a material (in MatPropItemCollection and
            MatItemAssignments)
        """
        if self._frozen:
            msg('MatItemAssignments is frozen. No chages allowed.')
            return False

        if matType.lower() == 'rock':
            rmKeys = [n for n,v in self.iteritems() if v.rock == name]
        elif matType.lower() == 'fill':
            rmKeys = [n for n,v in self.iteritems() if v.fill == name]
        else:
            msg('matType has to be either "rock" or "fill". Skipped removal')
            return False

        for rmKey in rmKeys:
            self.removeAssignment(rmKey)

        self._meta.log('removed %s:%s from assignments' % (name, matType))
        return True


    def renameMat(self, name, matType, newName,
                  logText='renamed', force=False):
        """Renames all appearences of the requested material(Name). Note
        that assignment names will not be changed (use L{renameAssignment})
        @param name: old materialName that will be replaced in all assignments
        @param matType: the materialName referes to "rock" or "fill"
        @param newName: new materialName
        @note: this function is usually called from parent-MatPropsCollection
            to rename a material (in MatPropItemCollection and
            MatItemAssignments)
        """
        if self._frozen:
            msg('MatItemAssignments is frozen. No chages allowed.')
            return False

        if matType.lower() == 'rock':
            rnKeys = [n for n,v in self.iteritems() if v.rock == name]
        elif matType.lower() == 'fill':
            rnKeys = [n for n,v in self.iteritems() if v.fill == name]
        else:
            msg('matType has to be either "rock" or "fill". Skipped removal')
            return False

        for rnKey in rnKeys:
            setattr(self[rnKey], matType, newName)
            self[rnKey]._meta.log('%s %s from %s to %s' %
                                  (logText, matType, name, newName))



class MatPropsCollection(object):
    """Datastructure to store and manage material properties.

    Structure:
        MatPropsCollection
        E{|--} rock                          # dict-like MatItemCollection
        E{....|--} (rockName, L{MatItem})
        E{|--} fill                          # dict-like MatItemCollection
        E{....|--} (fillName, L{MatItem})
        E{|--} fillAssignments               # dict-like MatItemFillAssignments
        E{....|--} (fillAssignmentName, L{Assignment})

    Usage:
        >>> from bae.material_01 import (MatPropsCollection,
        ...                              PlasticLevkovitchReuschAniso)

        creating a new matPropsCollection
        >>> matProps = MatPropsCollection('Test00')

        adding fill and rock materials
        >>> defaultProps = {"dampingAlpha" : 0.5,
        ...         "density": 2700,
        ...         "voidStiffnessRatio": 5E-6,
        ...         "critDamage" : 1000.,
        ...         "umatDelete":9,
        ...         "umatDepvarNum":16,
        ... }
        >>> rockData = [#name    UCS   GSI  mr(V14)
        ...     ('matA', 40E6,  60, 600),
        ...     ('matB', 120E6, 60, 400),
        ...     ('matC', 80E6,  70, 600),
        ...    ]
        >>> # create a fill
        >>> fill = PlasticLevkovitchReuschAniso()
        >>> fill.fromComputator(40,60,version='V14.3')
        >>> matProps.fill.addItem('defaultFill', fill)
        >>> # lopping all rock materials
        >>> for name, ucs, gsi, mr in rockData:
        >>>     rock = PlasticLevkovitchReuschAniso()
        >>>     rock.fromComputator(ucs, gsi, mr=mr, version='V14.3')
        >>>     matProps.rock.addItem(name, rock, **defaultProps)
        >>>     matProps.fillAssignments.addAssignment(name, 'defaultFill')

        export to abaqus input
         >>> matProps.exportAbaqusInput('test_matV14.inp')

        export to excel
        Note that as long as a material is a (raw) computator material the GSI
        will be added and documented automatically.
        You can change/add GSI manually
         >>> matProps.rock['matA'].GSI = 40
         >>> matProps.rock.toXlsx('testXlsx.xlsx', plotMaxS3=10E6)

        reading from abaqus input
         >>> matPropsIn = MatPropsCollection('Test01')
         >>> matPropsIn.importAbaqusInput('test_matV14.inp')

    @ivar name: name of this MatPropsCollection, will be used to identify
            if you copy matItems from one Collection to an other
    @ivar rock: L{MatItemCollection} holding material props for rock
    @ivar fill: L{MatItemCollection} holding material props for fill

    """

    def __init__(self, name=None, creationStr=None):
        """
        @param name: name of MatPropsCollection. This will be used to identify
            if you copy matItems from one Collection to an other
        @param creationStr: logs some more detailed information
        @note: Would be nice to have a copy constructor here in future
        """
        if name is None:
            msg("WARNING! You should assign a proper name to a new"
                " MatPropsCollection to avoid confusion.")
            name = "NewMatPropsCollection"

        if creationStr is None:
            creationStr = 'new for %s' % name

        self.name = name
        self._meta = MetaData(creationStr=creationStr)

        self.fill = MatItemCollection(
            'fill', parent=self, creationStr=creationStr)
        self.rock = MatItemCollection(
            'rock', parent=self, creationStr=creationStr)

        self.fillAssignments = MatItemFillAssignments(
            parent=self, creationStr=creationStr)

        self._frozen = False

    @property
    def assignments(self):
        msg("DEPRECATED! Accessing the MatItemFillAssignments"
            " (rock-defaultFill) via assignments attribute is deprecated."
            " Use MatpropsCollection.fillAssignments instead")
        return self.fillAssignments

    def freeze(self, doIt=True):
        """Freezes the MatPropsCollection and all its subStructures so that
        theycan't be edited using matpropscollection-fuctions anymore.
        Use self.freeze(False) to unfreeze.
        """
        self._frozen = doIt
        self.rock.freeze(doIt)
        self.fill.freeze(doIt)
        self.fillAssignments.freeze(doIt)
        if doIt:
            self._meta.log('MaterialPropertyCollection got frozen', level=0)
        else:
            self._meta.log('MaterialItemCollection got unfrozen', level=0)


    def toDict(self):
        return dict(name=self.name,
                    metaData=self._meta.toDict(),
                    frozen=self._frozen,
                    rock=self.rock.toDict(),
                    fill=self.fill.toDict(),
                    fillAssignments=self.fillAssignments.toDict())


    @classmethod
    def fromDict(cls, input):
        newMPC = cls(input['name'])

        newMPC.rock = MatItemCollection.fromDict(input['rock'])
        newMPC.fill = MatItemCollection.fromDict(input['fill'])
        newMPC.fillAssignments = MatItemFillAssignments.fromDict(
                                        input['fillAssignments'])
        newMPC._frozen = input['frozen']
        newMPC._meta = MetaData.fromDict(input['metaData'])

        newMPC.rock._parent = newMPC
        newMPC.fill._parent = newMPC
        newMPC.fillAssignments._parent = newMPC

        return newMPC


    def importAbaqusInput(self, inpNames, ignoreFill=False):
        """Imports material cards from abaqus.inp files.
        @param inpNames: can be a abaqus.inp or list of files
        @param ignoreFill: if True, the (default)fill of each abaqus-material
            will not get stored in MatPropsCollection. The related
            fillAssignments will not be created as well.
        """
        if self._frozen:
            msg('MatPropsCollection is frozen. No chages allowed.')
            return False

        if isinstance(inpNames, basestring):
            inpNames = [inpNames,]

        keysAlreadyInMpc = self.rock.keys()
        self.rock.importAbaqusInput(inpNames)
        if ignoreFill:
            return self

        self.fill.importAbaqusInput(inpNames)

        inpBaseNames = '; '.join([os.path.split(name)[-1]
                                 for name in inpNames])
        newNames = self.rock.keys()[len(keysAlreadyInMpc):]
        for name in newNames:
            self.fillAssignments.addAssignment(
                name, name, creationStr='loaded from %s' % inpBaseNames)

        return self


    def exportAbaqusInput(self, inpName, umatDelete=None, umatDepvarNum=None,
                          checkTypes = True, storeMatNamesList=False):
        """Writes the assigned (!) materials (rock-fill-pairs) in
        MatPropsCollection into an abaqus input.
        Names and order will be according to self.fillAssignments.

        @param name: inpName of abaqus input file
        @param umatDelete: abaqus deletion flag. Note that this will overwrite
             stored values of umatDelete in all MatItems
        @param umatDepvarNum: number of Defvars. Note that this will overwrite
             stored values of umatDelete in all MatItems
        @param checkTypes: will check if exported materials are of same type.
        """
        from bae.abq_model_02 import Model
        if isinstance(inpName, basestring):
            fid = open(inpName, 'w')
        else: # assume inpName to be a stringIO
            fid = inpName

        if storeMatNamesList:
            addStar =lambda s: '**' + s
            tt = '\n'.join(map(addStar, self.matNamesList.split('\n')))
            tt += '\n'
            fid.write(tt)

        def setUmatVal(mat, valName, val):
            if val is None and not getattr(mat, valName, False):
                raise KeyError(('No %s specified. You can set this value as'
                                ' keyWordArguement in exportAbaqusInput-'
                                'function') % valName)
            elif val is None:
                return
            else:
                setattr(mat, valName, val)

        for rock in self.rock.values():
            setUmatVal(rock, 'umatDelete', umatDelete)
            setUmatVal(rock, 'umatDepvarNum', umatDepvarNum)

        for fill in self.rock.values():
            setUmatVal(fill, 'umatDelete', umatDelete)
            setUmatVal(fill, 'umatDepvarNum', umatDepvarNum)

        for mat in self.toAbaqusVUmatDataList(checkTypes):
            m = Model()
            mat.updateAbqModel(m)
            m.write(fid)
        fid.close()

    def toAbaqusVUmatDataList(self, checkTypes = True):
        """Service fuction to create a list of L{bae.material_01.AbqVumatData}
        from MatPropsItems.
        @param checkTypes: set True to enforce same material types
        """
        # assignments = [(name, rockMat, fillMat, masterName), ...]
        assignments = [(n,) + a.getAssignedMaterials() + (a.master,)
                        for n, a in self.fillAssignments.iteritems()]

        # check if rock and fillType is the same
        differentRockFillTypes = [q[0] + ':\t ' + q[1]._matType.__name__
                                  + ' vs. ' + q[2]._matType.__name__
                                  for q in assignments
                                  if not q[1]._matType == q[2]._matType]

        if differentRockFillTypes:
            raise ValueError('The following assignments have different' +
                             'MaterialTypes for Rock and Fill: \n' +
                             '\n'.join(differentRockFillTypes) )

        uniqueMatTypes = set([t[1]._matType for t in assignments])
        if len(uniqueMatTypes) > 1 and checkTypes:
            msg('WARNING! Found different MatTypes to be writen. '
                'Take care that your vumat can handle this!!!')

        vUMatList = [self.itemToAbaqusVUmat(*q) for q in assignments]

        return vUMatList


    @staticmethod
    def itemToAbaqusVUmat(name, rock, fill, masterName = 'fill'):
        """Service fuction that converts a rock-fill pair into a  L{bae.material_01.AbqVumatData}
        instance.
        """
        from bae.material_01 import AbqVumatData

        if masterName.lower() == 'rock':
            master = rock
        elif masterName.lower() == 'fill':
            master = fill
        else:
            raise ValueError('Do not know a masterType called %s' % masterName)

        abq = AbqVumatData()
        abq.name            = name
        abq.plasticLRRock   = rock
        abq.plasticLRFill   = fill
        for var in ['density', 'dampingAlpha', 'umatDelete', 'umatDepvarNum',
                    'voidStiffnessRatio', 'critDamage']:
            try:
                setattr(abq, var, getattr(master, var))
            except AttributeError:
                raise AttributeError("Material %s / %s has no attribute %s"
                                     % (masterName, name, var))
        return abq


    @property
    def matNamesList(self):
        """Returns the 'ready-to-copy' python code for matNamesList as used
        by some pre and post scripts.
        Usage:
            >>> print myMatPropsCollection.matNamesList
        """
        matCodes = '\n'.join(map(lambda s: '    "%s",' %s,
                         self.fillAssignments))
        return 'matNamesList=[\n%s]\n' % matCodes

    def plotYieldSurfs(self, rockNames = None, fillNames = None, **kwargs):
        # prepare default output
        _assignedRock = [assign.rock for assign in self.fillAssignments.values()]
        _assignedFill = [assign.fill for assign in self.fillAssignments.values()]
        # remove duplicates but preserve order
        assignedRock, assignedFill = [],[]

        for aR in _assignedRock:
            if aR not in assignedRock:
                assignedRock.append(aR)

        for aF in _assignedFill:
            if aF not in assignedFill:
                assignedFill.append(aF)

        if rockNames is None:
            rockNames = assignedRock
        elif isinstance(rockNames, basestring) and rockNames.lower() == 'all':
            rockNames = self.rock.keys()

        if fillNames is None:
            fillNames = assignedFill
        elif isinstance(fillNames, basestring) and fillNames.lower() == 'all':
            rockNames = self.fill.keys()

        self.rock.plotYieldSurfs(matNames=rockNames, **kwargs)
        self.fill.plotYieldSurfs(matNames=fillNames, **kwargs)


    def plotAniso(self, rockNames = None, fillNames = None, **kwargs):
        # prepare default output
        _assignedRock = [assign.rock for assign in self.fillAssignments.values()]
        _assignedFill = [assign.fill for assign in self.fillAssignments.values()]
        # remove duplicates but preserve order
        assignedRock, assignedFill = [],[]

        for aR in _assignedRock:
            if aR not in assignedRock:
                assignedRock.append(aR)

        for aF in _assignedFill:
            if aF not in assignedFill:
                assignedFill.append(aF)

        if rockNames is None:
            rockNames = assignedRock
        elif isinstance(rockNames, basestring) and rockNames.lower() == 'all':
            rockNames = self.rock.keys()

        if fillNames is None:
            fillNames = assignedFill
        elif isinstance(fillNames, basestring) and fillNames.lower() == 'all':
            rockNames = self.fill.keys()

        self.rock.plotAniso(matNames=rockNames, **kwargs)
        self.fill.plotAniso(matNames=fillNames, **kwargs)


    def toXlsx(self, workBook, rockNames = 'all', fillNames = 'all', **kwargs):
        """Creates an excel file/workBook holding formated properties and
        plots of requested rock and fill materials. Each material gets
        stored on a new sheet.

        @params workBook: excel-fileName or xlsxwriter.Workbook-instance
        @params rockNames: list of rockNames to be exported. Use None to
            export all assigned rock materials and 'all' to export all
            items stored in self.rock
        @params fillNames: list of fillNames to be exported. Use None to
            export all assigned fill materials and 'all' to export all
            items stored in self.fill
        @note: see L{MatItem.toXlsx} for more keyword-arguments to adjust
            output for your needs.
        """
        try:
            import xlsxwriter as xls
        except ImportError:
            raise ImportError("xlsxwriter is not installed. Use:\n" +
                              ">>> pip install XlsxWriter\n" +
                              "on regular python or\n" +
                              ">>> conda install -c anaconda xlsxwriter\n"
                              "if you are using anaconda")

        if not isinstance(workBook, xls.Workbook):
            workBook = xls.Workbook(workBook, {'nan_inf_to_errors': True})
            closeIfDone = True
        else:
            closeIfDone = False

        # prepare default output
        _assignedRock = [assign.rock for assign in self.fillAssignments.values()]
        _assignedFill = [assign.fill for assign in self.fillAssignments.values()]
        # remove duplicates but preserve order
        assignedRock, assignedFill = [],[]

        for aR in _assignedRock:
            if aR not in assignedRock:
                assignedRock.append(aR)

        for aF in _assignedFill:
            if aF not in assignedFill:
                assignedFill.append(aF)

        if rockNames is None:
            rockNames = assignedRock
        elif isinstance(rockNames, basestring) and rockNames.lower() == 'all':
            rockNames = self.rock.keys()

        if fillNames is None:
            fillNames = assignedFill
        elif isinstance(fillNames, basestring) and fillNames.lower() == 'all':
            fillNames = self.fill.keys()

        self.rock.toXlsx(workBook, matNames=rockNames, **kwargs)
        self.fill.toXlsx(workBook, matNames=fillNames, **kwargs)

        if closeIfDone:
            workBook.close()


    @staticmethod
    def getDataFromCsv(csvFile, varNames, colIdxs, headerRows=1,
                       delimiter=',', converterDict={}):
        """Reads abitrary data from a csvFile.

        Usage:
            >>> converters = {'rho': lambda r: 2700 if not r else float(r),}
            >>> for mat in matCollection.getDataFromExcel('MATPROPS.csv',
            ...                      ['matKey', 'rho', 'UCSi', 'GSI'],
            ...                      [2,3,4,13],
            ...                      converterDict=converters):
            ...     name = mat['matKey']
            >>>     rock = PlasticLevkovitchReuschAniso()
            >>>     rock.fromComputator(mat.UCSi*1E6, mat.GSI,
            ...                         version='V14.3')
            >>>     matCollection.rock.addItem(name, rock, density=mat.rho,
            ...                                **defaultProps)

        @param csvFile: path of excel file
        @param varNames: list of (unique) variable names, will appear as keys
            in returned dictionary
        @param colIdxs: list of column indices to read. You can set indexes of
            varNames that are not present in csv-file to None. The values than
            are assumed to be None as well and might get converted later on.
            This might be usefull to create unified interfaces that always
            require a varName/value pair.
        @param converterDict: a dictionary holding converter functions that
            get called on the fly. The default conversion (if no converter
            function is given for a varname) tries to get an integer and than a
            float. If both fail, the cell string is returned.
        @param headerRows: number of headerlines.
        @returns:
        """
        if not len(varNames) == len(colIdxs):
            raise ValueError("List of varNames and colIdxs must be of same"
                             " length.")

        def defaultConverter(dat):
            try:
                return int(dat)
            except ValueError:
                pass
            try:
                return float(dat)
            except ValueError:
                pass
            return dat

        conv = dict((name, converterDict.get(name, defaultConverter))
                    for name in varNames)
        data = []
        with open(csvFile) as fp:
            for cnt, line in enumerate(fp):
                if cnt < headerRows:
                    continue
                cells = map(str.strip, line.split(delimiter))
                dat = []
                for name, colIdx in zip(varNames, colIdxs):
                    if colIdx is None:
                        val = None
                    else:
                        val = cells[colIdx]
                    dat.append( (name, conv[name](val)) )

                data.append(Container(**dict(dat)))

        return data


    @staticmethod
    def getDataFromExcel(xlsxFile, varNames, colIdxs, sheetName='M01',
                         isNewMatColIdx=0, headerRows=None,
                         lastRow=None, converterDict={}):
        """Reads abitrary data from (matprops) excel spread sheets

        Usage:
            >>> for mat in matCollection.getDataFromExcel('MATPROPS.xlsm',
            ...                      ['matKey', 'density', 'UCSi', 'GSI'],
            ...                      [2,3,4,13],
            ...                      sheetName='M02'):
            ...     name = mat.matKey
            >>>     rock = PlasticLevkovitchReuschAniso()
            >>>     rock.fromComputator(mat.UCSi*1E6, mat.GSI, version='V14.3')
            >>>     matCollection.rock.addItem(name, rock,
            ...                                density=mat.density,
            ...                                **defaultProps)

        @param xlsxFile: path of excel file
        @param varNames: list of (unique) variable names, will appear as keys
            in returned dictionary
        @param colIdxs: list of column indices to read
        @param isNewMatColIdx: column that indicates the line to read. Set None
            to get every column.
        @param sheetName: name of sheet (e.g. M04)
        @param headerRows: number of headerlines. If None an good old excel-
            matpros-sheet is assumed, whos last headerRow starts with a
            cell holding the string 'no.'
        @param lastRow: index of last row to be evaluated. Set None to get all
            rows holding (any) data.
        @returns: a list of L{Containers<bae.misc_01.Container>}
            for each row where isNewMatColIdx has a notNone value.
            The Container has varNames as keys/attributes and the related cell
            values as values.
        @note: If the excel file is opened in excel the reader will fail.
        @note: If one materials spans multiple rows you are interested in,
            you'll have to reassign these rows to the material in a second
            step.
        @note: Row and column indexing starts with 0.
        """
        import zipfile
        if not len(varNames) == len(colIdxs):
            raise ValueError("List of varNames and colIdxs must be of same"
                             " length.")

        try:
            import xlrd
        except ImportError:
            raise ImportError("Can't read %s because 'xlrd' (excel reder)"
                              " is not installed." % xlsxFile)

        if not os.path.isfile(xlsxFile):
            raise IOError('No File named %s' % xlsxFile)

        try:
            workbook = xlrd.open_workbook(xlsxFile, on_demand=True)
        except zipfile.BadZipfile:
            raise IOError("Can not read %s." % xlsxFile +
                          " Mostlikely it is currently opened in Excel")

        worksheet = workbook.sheet_by_name(sheetName)

        def defaultConverter(dat):
            return dat

        conv = dict((name, converterDict.get(name, defaultConverter))
                    for name in varNames)

        if lastRow:
            nRows = min(lastRow+1, worksheet.nrows)
        else:
            nRows = worksheet.nrows

        currentRow = 0
        if headerRows is not None:
            currentRow = headerRows
        else:
            while currentRow < nRows:
                row = worksheet.row(currentRow)
                currentRow += 1
                if row[0].value.lower().startswith('no.'):
                    break

        data = list()
        while currentRow < nRows:
            row = worksheet.row(currentRow)
            if isNewMatColIdx is None:
                isNewMatLine = True
            else:
                isNewMatLine = row[isNewMatColIdx].value

            if isNewMatLine:
                dat = []
                for name, colIdx in zip(varNames, colIdxs):
                    if colIdx is None:
                        val = None
                    else:
                        val = row[colIdx].value
                        try:
                            if val.strip() == '':
                                val = None
                        except:
                            pass

                    dat.append( (name, conv[name](val)) ) # apply converter
                data.append(Container(**dict(dat)))
            currentRow += 1

        return data


    @staticmethod
    def getUCSAndGSIFromExcel(xlsxFile, sheetName='M01',
                          matCodeCol=None, GSICol=None, UCSiCol=None,
                          headerRows=None):
        """Reads UCSi and GSI values from Matprops-excel sheet. See
        L{getDataFromExcel} for more flexibility.

        Usage:
            >>> matProps = MatPropsCollection('M02')
            >>> excelData =  matProps.getUCSAndGSIfromExcel('matProps.xlsx',
            ...                                             sheetName='M01')
            >>> for matName, data in excelData.iteritems():
            >>>     rock = PlasticLevkovitchReuschAniso()
            >>>     rock.fromComputator(data.UCSi, data.GSI)
            >>>     matCollection.rock.addItem(matName, rock, **defaultProps)


        @param xlsxFile: path of excel file
        @param sheetName: name of sheet (e.g. M04)
        @param matCodeCol: column to get the matcode from (columnA=0)
        @param GSICol: column to get the GSI from
        @param UCSiCol: column to get the UCSi from
        @param headerRows: number of headerlines
        @note: If the excel file is opened in excel the reader will fail
        """
        import zipfile
        try:
            import xlrd
        except ImportError:
            raise ImportError("Can't read %s because 'xlrd' (excel reder)"
                              " is not installed." % xlsxFile)

        if not os.path.isfile(xlsxFile):
            raise IOError('No File named %s' % xlsxFile)

        try:
            workbook = xlrd.open_workbook(xlsxFile, on_demand=True)
        except zipfile.BadZipfile:
            raise IOError("Can not read %s." % xlsxFile +
                          " Mostlikely it is currently opened in Excel")

        worksheet = workbook.sheet_by_name(sheetName)

        nRows = worksheet.nrows
        currentRow = 0

        if any([c is None for c in [matCodeCol, GSICol, UCSiCol]]):
            msg('At least one of the required columnIndexes is None.'
                ' Will try to detect coloums.')
            #identify and parse header
            while currentRow < nRows:
                row = worksheet.row(currentRow)
                if row[0].value.lower().startswith('no.'):
                    for ii,cell in enumerate(row):
                        val = cell.value.lower()
                        if 'code' in val and 'mat' in val:
                            matCodeCol = ii
                        if 'gsi' in val:
                            GSICol = ii
                        if 'ucs' in val and UCSiCol is None:
                            UCSiCol = ii
                    currentRow += 1
                    break
                currentRow += 1

        if any([c is None for c in [matCodeCol, GSICol, UCSiCol]]):
            print('matCodeCol: ', matCodeCol)
            print('GSICol: ', GSICol)
            print('UCSiCol: ', UCSiCol)
            raise ValueError("Can't identify one of the columnIndexes above."
                             " Try to specify manually")

        if headerRows is not None:
            currentRow = headerRows
        data = OrderedDict()
        while currentRow < nRows:
            row = worksheet.row(currentRow)
            matCode = row[matCodeCol].value
            if matCode:
                GSI = row[GSICol].value
                UCSi = row[UCSiCol].value
                data[str(matCode)] = Container(GSI=GSI, UCSi=UCSi*1E6)
            currentRow += 1

        return data

#}
###############################################################################


if __name__ == '__main__':
    print('No Syntax Errors')
