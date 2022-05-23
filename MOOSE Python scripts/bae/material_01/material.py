# -*- coding: utf-8 -*-
"""classes for material properties.
"""

try:
    import numpy as np
except ImportError:
    np = None

from itertools import izip

from bae.log_01 import msg
from yieldsurfaces import (ExtendedMenetreyWillam,
                           AnisotropicExtendedMenetreyWillam,
                           SinglePlaneMohrCoulomb,
                           QuasiMohrCoulomb,
                           pQThetaToEigs,)
from models import (computatorV12, computatorV13, computatorV14,
                    restoreComputatorInputs)


#{Helper classes and fuctions #################################################

def isInTol(valRef, valCmp, tol=0.001):
    '''Check if two values are equal within given relative tolerance.
    '''
    return abs(valRef-valCmp) <= abs(tol*valRef)


def getValsInTol(valList, testVal, tol=0.001, tryExact=True):
    '''Return matched-list of values out of valList which are close to
    testVal. If the matched-list has more then one items and tryExact
    is set True the tolerance is reduced recusively to force an exact
    match.
    '''
    matched = [val for val in valList
               if isInTol(val, testVal, tol=tol)]

    if tryExact and (len(matched) > 1) and (tol > 1E-12):
        matched = getValsInTol(valList, testVal,
                               tol=0.1*tol, tryExact=True)
    return matched


def deltaDict(sRef, sCmp, ignoreNotFound=False, tol=0.001):
    """Compares dictionary values and returns those that are not
    equal within the specified tolerance.

    @param ignoreNotFound: If a key is not available in either sRef or sCmp
        then the default behaviour is to add an entry to the resulting
        dictionary with the missing item replaced by a message string like
        'no attribute named xyz'. If this argument is True then such an item
        is not generated and the item missing in one of the dictionaries is
        simply ignored.
    @param tol: Tolerance of numeric values relative to values in the first
        argument sRef.

    @returns: A dictionary containing only items that differ in the given input
        dictionaries. The values in this dict are tuples: first the
        corresponding value of the first input sRef and secondly the value from
        the second input sCmp.
    """
    delta = {}
    for name in set(sRef).union(sCmp):
        try:
            if isinstance(sRef[name], basestring):
                if not sRef[name]==sCmp[name]:
                    delta[name] = (sRef[name], sCmp[name])
            else:
                if not isInTol(sRef[name], sCmp[name], tol=tol):
                    delta[name] = (sRef[name], sCmp[name])
        except KeyError:
            if not ignoreNotFound:
                try:
                    delta[name] = (sRef[name], 'no attribute named %s' % name)
                except KeyError:
                    delta[name] = ('no attribute named %s' % name, sCmp[name])

    return delta


class VusubsVarFillMat(list):
    """List of L{material_01.matpropscollection.MatItems<MatItems>} and/or
    L{DefaultPst<pst-sample>} list which will print as a formated
    string of pst-samples that can directly be used in vusubsParam.py.

    Example for usage in vusubs_param.py:
        >>> from bae.material_01.defaults import defaultFillAniso
        >>> from bae.material_01.material import VusubsVarFillMat
        >>> _varFills = VusubsVarFillMat()
        >>> _varFills += defaultFillAniso['FILL_ELAST_20MPA']
        >>> _varFills += defaultFillAniso['FILL_PLAST_100xx']
        >>> maxNumberFill = _varFills.maxNumberFill
        >>> fillPropsStartIdx = _varFills.fillPropsStartIdx
        >>> PROPS_fill = _varFills.PROPS_fill
        >>> print(_varFills) #for debugging or idividual changes afterwards

    @ivar fullOutput: if True, the complete vusubs_param.py section (including
        maxNumberFill and fillPropsStartIdx) will be returned using print/str
    @ivar indent: default indention
    """
    indent = 4*' '
    fullOutput = True
    fillTimeOffset = 1000


    @property
    def maxNumberFill(self):
        """Number of stored items.
        Example for usage in vusubs_param.py:
            >>> maxNumberFill = myVusubsVarFillMat.maxNumberFill
        """
        self.testConsistency()
        return len(self)


    @property
    def fillPropsStartIdx(self):
        """Indexes of new materials in joined varFill tuple.
        Example for usage in vusubs_param.py:
            >>> fillPropsStartIdx = myVusubsVarFillMat.fillPropsStartIdx
        """
        self.testConsistency()
        if not self:
            return []

        if self.isItemAniso(self[0]):
            firstTupSize = 6
        else:
            firstTupSize = 4
        pstTupSize = 6

        idxs = [1,]
        for mat in self:
            idxs.append(idxs[-1] + firstTupSize + len(mat)*pstTupSize)
        idxs.pop() #remove last
        return idxs


    @property
    def PROPS_fill(self):
        """Tuple of material properties of all stored varFill items.
        Example for usage in vusubs_param.py:
            >>> PROPS_fill = myVusubsVarFillMat.PROPS_fill
        """
        return self.toAllTuple()


    def __str__(self):
        """Overwrites list.__str__ function. print() or str() will now return
        a formated string that can directly copy-pasted/used in
        vusubs_param.py. If self.fullOutput is True the complete varfill
        properties section will be returend.
        """
        if not self.fullOutput:
            return self.toAllTupleString()

        output = 'maxNumberFill = %d\n' % self.maxNumberFill
        output += 'fillPropsStartIdx = %s\n' % str(self.fillPropsStartIdx)
        output += 'PROPS_fill = [\n%s\n]' % self.toAllTupleString()
        return output


    def __iadd__(self, other):
        """Redicted to self.__add__."""
        return self.__add__(other)


    def __add__(self, other):
        """Overwrites the list.__add__ function. You can add other
        VusubsVarFillMat instances, single MatItems or single PST-sample
        lists.
        """
        if hasattr(other, 'UCSi'):
            self.append(other)
        elif isinstance(other, self.__class__):
            self.extend(other)
        else:
            raise ValueError('Can only concatenate VusubsVarFillMat instances')

        return self


    def __rmul__(self, N):
        """Redicted to self.__mul__."""
        return self.__mul__(N)


    def __imul__(self,N):
        """Redicted to self.__mul__."""
        return self.__mul__(N)


    def __mul__(self,N):
        return self.__class__([item for i in range(N) for item in self])


    @staticmethod
    def isItemAniso(item):
        """Tests if an item is anisotropic (returns True) or isotropic (False).
        """
        try:
            item.anisoS, item.anisoN
            return True
        except AttributeError:
            return False


    def testConsistency(self):
        """Checks if all stored items/materials do have the same output type.
        Currently, there are only two output types: isotropic and anisotropic.
        @raises: ValueError if items are not consistent
        """
        # test if only aniso or isotropic items are stored
        anisoItems = [self.isItemAniso(item) for item in self]
        if len(self)>1 and not (all(anisoItems) or any(anisoItems)):
            print 'anisotropic items:', anisoItems
            raise ValueError("Some matItems seem to be isotropic and some"
                             " anisotropic.")

    def toAllTupleString(self):
        """Formats all stored materials to one string."""
        self.testConsistency()
        return '\n'.join([self.toSingleTupleString(item) for item in self])


    @classmethod
    def toSingleTupleString(cls, item, name=None):
        """Formats a single matItem into tupleString:
            - Isotropic:   UCS, e, 1/a, N,(PST, G, K, s, mb, d,)
            - Anisotropic: UCS, e, 1/a, s_a, n_a, N,(PST, G, K, s, mb, d,)
        """
        if name is None:
            try:
                name = item.name
            except AttributeError:
                name = ''

        general = '%.4E, %.4f, %.6f, ' % (item.UCSi,item[0].yieldSurf.e,
                                          1/float(item[0].yieldSurf.a))
        if cls.isItemAniso(item):
            general += '%.6f, %.6f, ' % (item[0].yieldSurf.anisoS,
                                         item[0].yieldSurf.anisoN)

        general += '%d, ' % len(item)
        general += '# %s' % name
        general = 2*cls.indent + general
        def sampleStr(sample):
            el, ys = sample.elastic, sample.yieldSurf
            ss = (6*'%.4E, ') % (sample.pst, el.G, el.K, ys.s, ys.mb,
                                 sample.dilation)
            return 3*cls.indent + ss.rstrip()
        return '\n'.join([general,] + map(sampleStr, item))


    def toAllTuple(self):
        """Returns the joined tuple of all stored matItems"""
        self.testConsistency()
        out = ()
        for ii, item in enumerate(self, start=1):
            try:
                name = item.name
            except AttributeError:
                name = 'varFill_%d' % (ii*self.fillTimeOffset)
            out += self.toSingleTuple(item, name)
        return out


    @classmethod
    def toSingleTuple(cls, item):
        """Returns a single matItem tuple"""
        out = item.UCSi, item[0].yieldSurf.e, 1/float(item[0].yieldSurf.a)
        if cls.isItemAniso(item):
            out += item[0].yieldSurf.anisoS, item[0].yieldSurf.anisoN
        for sample in item:
            el, ys = sample.elastic, sample.yieldSurf
            out += sample.pst, el.G, el.K, ys.s, ys.mb, sample.dilation
        return out
#}#End Helper classes and functions
###############################################################################

###############################################################################
#{Elasticity
class ElasticLinear(object):
    """Class to store the elastic properties: simple isotropic linear
    elasticity.

    @ivar E: Youngs modulus
    @ivar nu: Poissons ratio
    @ivar K: bulk modulus
    @ivar G: shear modulus

    @Note: All units SI!
        setting E or nu will update K and G and vice versa
    """

    def __init__(self, *args, **kwargs):
        """Create a ElasticLinear object:
            - as an empty object
            - derived from E and nu as kwargs --> sets K and G automatically
            - derived from K and G as kwargs  --> sets E and nu automatically
            - as a copy of another ElasticLinear object passed

        Usage:
         >>> elastic = ElasticLinear()
         >>> elastic.nu = .2
         >>> elastic.E  = 10E9
         >>> print('K=%.1f' % elastic.K)
         >>> print('G=%.1f' % elastic.G)
         >>> # or via kwargs E and nu
         >>> elastic = ElasticLinear(E=1E9, nu=.2)
         >>> # or via kwargs K and G
         >>> elastic = ElasticLinear(G=4166666666.666, K=5555555555.555)
         >>> # or from copy constructor
         >>> elastic2 = ElasticLinear(elastic)

        Make sure nu != 0.5 (actually nu < 0.5)

        @kwarg E: Youngs modulus
        @kwarg nu: Poissons ratio
        @kwarg K: bulk modulus
        @kwarg G: shear modulus

        @Note: All units SI!
        """
        self.__E = None
        self.__nu = None
        self.__G = None
        self.__K = None

        #Copy constructor
        if not args:
            pass
        elif args and isinstance(args[0], ElasticLinear):
            self.setEnu(args[0].E, args[0].nu)
            return
        else:
            raise ValueError('Only ElasticLinear-objects can be passed as'
                             ' positional arguments (copy constructor).')

        if kwargs:
            if "E" in kwargs and "nu" in kwargs:
                self.setEnu(kwargs["E"], kwargs["nu"])
            elif "G" in kwargs and "K" in kwargs:
                self.setGK(kwargs["G"], kwargs["K"])
            else:
                raise ValueError("Either E and nu or G and K have to be specified"
                                 " for an ElasticLinear-object.")

    def copy(self):
        """Returns a copy of ElasticLinear-object"""
        return self.__class__(self) # using copy constructor

    def toDict(self):
        """Returns a dictionary holding E and nu
        """
        return dict(E=self.E, nu=self.nu)

    @classmethod
    def fromDict(cls, input):
        """Creates an L{ElasticLinear}-instance from dictionary holding
        elastic constants E AND nu or K AND G.

        @note: this is equivalent to calling >>> ElasticLinear(**input)

        """
        return cls(**input)

    @property
    def E(self):
        """Youngs modulus in Pa"""
        return self.__E

    @E.setter
    def E(self, value):
        self.__E = float(value)
        if not self.__nu is None:
            self.setEnu(self.__E, self.__nu)

    @property
    def nu(self):
        """Poissons ratio"""
        return self.__nu

    @nu.setter
    def nu(self, value):
        self.__nu = float(value)
        if not self.__E is None:
            self.setEnu(self.__E, self.__nu)

    @property
    def G(self):
        """Shear modulus in Pa"""
        return self.__G

    @G.setter
    def G(self, value):
        self.__G = float(value)
        if not self.__K is None:
            self.setGK(self.__G, self.__K)

    @property
    def K(self):
        """Bulk modulus in Pa"""
        return self.__K

    @K.setter
    def K(self, value):
        self.__K = float(value)
        if not self.__G is None:
            self.setGK(self.__G, self.__K)

    def setVal(self, **kwargs):
        """DEPRECATED! Set values directly instead.
        Wrapps self.setEnu and self.setGK to a keyword-based method.
        """
        if ( (('E' in kwargs) or ('nu' in kwargs))
             and not
             (('G' in kwargs) or ('K' in kwargs)) ):
            self.E  = kwargs.get('E', self.E)
            self.nu = kwargs.get('nu', self.nu)


        elif ( (('G' in kwargs) or ('K' in kwargs))
               and not
               (('E' in kwargs) or ('nu' in kwargs)) ):
            self.G = kwargs.get('G', self.G)
            self.K = kwargs.get('K', self.K)

        else:
            raise AttributeError("Unknown keywords or ambiguous selection"
                                 " of E&nu / K&G : %s" % str(kwargs.keys()))

    def setEnu(self, E, nu):
        """Sets elastic constants E and nu and updates alternative constants
        K and G.
        @param E: Youngs Modulus
        @param nu: poisson ratio
        """
        self.__E = float(E)
        self.__nu = float(nu)
        self.__G = E/(2.0*(1.0+nu))
        self.__K = E/(3.0*(1.0-2.0*nu))

    def setGK(self, G, K):
        """Sets elastic constants G and K and updates alternative constants.
        @param G: shear modulus
        @param K: bulk modulus
        """
        self.__G = float(G)
        self.__K = float(K)
        self.__E = 9.0*K*G / (3.0*K + G)
        self.__nu = (1.5*K - G) / (3.0*K + G)

    def fromHoekBrown(self, UCSi, GSI, MR, D=0.0, nu=0.2):
        """Sets Youngs modulus from Hoek-Brown formulas.

        @param MR: modulus ratio
        @param D: disturbance factor

        @Note: All units SI!
        @note: from "Empirical estimation of rock mass modulus"
        by Hoek, Diederichs, Int Journal of Rock Mech and Mining Science, 2005
        """
        from math import exp

        # E_intact rock = UCSi*MR
        self.E = UCSi * MR *(0.02+(1.0-0.5*D)/(1.0+exp((60.0+15.0*D-GSI)/11.0)))
        self.nu = float(nu)

###############################################################################
#}# End Elasticity


###############################################################################
#{ Plastic-Strain-Samples
class SampleEMW2010(object):
    """PST-Sample for extended Menetrey-Willam-material

    This is one sample of the complete elasto-plastic behaviour of a
    particular isotropic material, i.e. elastic, plastic and dilation
    properties for one particular plastic strain (pst) value.

    The class holds index information in some hidden (@Tobias: why?) class
    variables describing the position of this data in the VUMAT-version
    introduced in 2010. I.e. those variables describe the structure of the
    *USER MAT card data.

    @ivar pst: plastic strain (pst) value
    @ivar elastic: an L{ElasticLinear}-object
    @ivar yieldSurf: a L{ExtendedMenetreyWillam}-object
    @ivar dilation: the d-parameter of type float
    """
    # Set Class-attributes
    _yieldType = ExtendedMenetreyWillam

    ### varName-to-index-lookUp for umatData

    # umatStructure of dialect vumat2010:
    # [fixSampleRock] + [fixSampleFill] + [fixMat] +
    # [nTupleRock] + nTupleRock*[varSampleRock]
    # [nTupleFill] + nTupleFill*[varSampleFill]
    _vumatDialect = 'vumat2010'

    # i) material specific data (index absolut), like critDamage,...
    _fixedMatIndex = dict(voidStiffnessRatio=6, critDamage=7)

    # ii) material specific data that will appear in each pst sample (index
    # absolut)
    _fixedSampleIndex = dict(UCSi=0, e=1, a=2)

    # iii) pst-tupleData (index relative tupleIndex)
    _varSampleIndex = dict(pst=0, G=1, K=2, s=3, mb=4, dilation=5)

    # iv) convert the following variables after reading
    _conversionRead = dict(a=lambda x: 1./float(x) if x else x, )

    # v) reconvert before writing, inversed _conversionRead
    _conversionWrite= dict(a=lambda x: 1./float(x) if x else x, )

    def __init__(self, pst=None, elastic=None, yieldSurf=None, dilation=None):
        """
        @param pst: plastic strain (pst) value
        @param elastic: an L{ElasticLinear}-object
        @param yieldSurf: a L{ExtendedMenetreyWillam}-object
        @param dilation: the d-parameter of type float
        @note: has no copy constructor yet. Use L{SampleEMW2010.copy}
        """

        if pst is not None:
            self.pst       = float(pst)
        if elastic is not None:
            self.elastic   = elastic
        if yieldSurf is not None:
            self.yieldSurf = yieldSurf
        if dilation is not None:
            self.dilation  = float(dilation)

    def copy(self):
        """Returns a (deep)copy of Plastic-Strain-Samples-object"""
        return self.__class__(pst=self.pst,
                               elastic=self.elastic.copy(),
                               yieldSurf=self.yieldSurf.copy(),
                               dilation=self.dilation)

    def toDict(self):
        """Stores pst, elastic-parameters, yield-parameters and dilation
        to a dictionary (see L{ElasticLinear.toDict} and L{ExtendedMenetreyWillam.toDict})
        as well.
        """
        try:
            elastic = self.elastic.toDict()
        except AttributeError:
            elastic = self.elastic

        try:
            yieldSurf = self.yieldSurf.toDict()
        except AttributeError:
            yieldSurf = self.yieldSurf

        return dict(pst=self.pst, elastic=elastic, yieldSurf=yieldSurf,
                    dilation=self.dilation)

    @classmethod
    def fromDict(cls, input):
        """Creates a Plastic-Strain-Sample-instance from a dictionary holding
        items for
            'pst': float
            'elastic': dict holding elastic properties (see L{ElasticLinear.fromDict})
            'yieldSurf': dict holding required yield-parameters for _yieldType
                (see e.g. L{ExtendedMenetreyWillam.fromDict})
            'dilation' : float
        """
        pst = input.get('pst')
        elasticInp = input.get('elastic', None)
        yieldSurfInp = input.get('yieldSurf', None)
        dilation = input.get('dilation', None)

        if yieldSurfInp is not None:
            yieldSurf = cls._yieldType.fromDict(yieldSurfInp)
        else:
            yieldSurf = None

        if elasticInp is not None:
            elastic = ElasticLinear.fromDict(elasticInp)
        else:
            yieldSurf = None

        return cls(pst=pst, elastic=elastic, yieldSurf=yieldSurf,
                   dilation=dilation)


    def fromHoekBrown(self, UCSi, GSI, mi, D=0.0, a=None, e=0.6, MR=None,
                      nu=0.2):
        """
        Updates elastic properties (only if a value is given for MR) and the
        yield surface according to Hoek-Brown formulas.

        Usage:
         >>> plasticPeak = SampleAEMW2017(
         ...     pst=0.0, dilation=0.2).fromHoekBrown(UCSi, GSI, mi, MR=MR)

        @param a: if not given will be computed by the Hoek-Brown formula
          based on GSI, otherwise usually: a=0.5

        @returns: self
        @Note: All units SI!
        """

        if MR is not None:
            if not hasattr(self, "elastic"):
                self.elastic = ElasticLinear()
            self.elastic.fromHoekBrown(UCSi, GSI, MR, D, nu)

        if not hasattr(self, "yieldSurf"):
            self.yieldSurf = self._yieldType()

        self.yieldSurf.e = e
        self.yieldSurf.fromHoekBrown(UCSi, GSI, mi, D=D, a=a)

        return self


    def compare(self, comp, ignoreNotFound=False, tol=0.001):
        """Compares a PST-sample with another.
        @param comp: PST-sample to compare self with
        @param ignoreNotFound: ignore if a property is not available in one of
            the compare PST-samples
        @param tol: relative tolerance between compared values
        @returns: dictionary for all distinct properties
        """
        delta = {}

        # delta yieldSurface
        deltaYield = deltaDict(self.yieldSurf.yieldParams(),
                               comp.yieldSurf.yieldParams(),
                               ignoreNotFound=ignoreNotFound,
                               tol=tol)

        # delta elastic
        deltaElast = deltaDict({'E':self.elastic.E, 'nu':self.elastic.nu},
                               {'E':comp.elastic.E, 'nu':comp.elastic.nu},
                               ignoreNotFound=ignoreNotFound,
                               tol=tol)

        # delta other (dilation, pst and maybe more)
        allKeys = (self._varSampleIndex.keys() +
                   self._fixedSampleIndex.keys())

        allKeys = (key for key in allKeys if key not in ['G','K'])

        remKeys = ['E', 'nu'] + self.yieldSurf._yieldParameters

        otherKeys = [k for k in allKeys if k not in remKeys]

        deltaOther = {}
        for key in otherKeys:
            try:
                deltaOther.update(deltaDict({key: getattr(self,key)},
                                            {key: getattr(self,key)},
                                            ignoreNotFound=ignoreNotFound,
                                            tol=tol))
            except AttributeError:
                if not ignoreNotFound:
                    deltaOther[key] = (getattr(self, key),
                                       'no attribute named %s' % key)


        delta.update(deltaYield)
        delta.update(deltaElast)
        delta.update(deltaOther)

        return delta


class SampleAEMW2017(SampleEMW2010):
    """PST-Sample for anisotropic extended Menetrey-Willam-material

    This is one sample of the complete elasto-plastic behaviour of a
    particular anisotropic material, i.e. elastic, plastic and dilation
    properties for one particular plastic strain (pst) value.

    The class holds index information in some hidden (@Tobias: why?) class
    variables describing the position of this data in the VUMAT-version
    for anisotropic behaviour introduced in 2017. I.e. those variables describe
    the structure of the corresponding *USER MAT card data.

    @ivar pst: plastic strain (pst) value
    @ivar elastic: an L{ElasticLinear}-object
    @ivar yieldSurf: a  AnisotropicExtendedMenetreyWillam-object
    @ivar dilation: the d-parameter of type float
    """
    _yieldType = AnisotropicExtendedMenetreyWillam

    ### varName-to-index-lookUp for umatData
    # i) material specific data (index absolut), like critDamage,...
    _fixedMatIndex = dict(voidStiffnessRatio=10, critDamage=11)

    # ii) material specific data that will appear in each pst sample (index
    # absolut)
    _fixedSampleIndex = dict(UCSi=0, e=1, a=2, anisoN=3, anisoS=4)

    # _varSampleIndex unchanged
    # _conversionRead unchanged
    # _conversionWrite unchanged


class SampleSPMC2010(SampleEMW2010):
    """PST-Sample for SinglePlaneMohrCoulomb

    The class holds index information in some hidden (@Tobias: why?) class
    variables describing the position of this data in the VUMAT-version
    for single plane Mohr Coulomb behaviour (possibly) introduced in 2010
    (or later). I.e. those variables describe the structure of the
    corresponding *USER MAT card data.

    @ivar pst: plastic strain (pst) value
    @ivar elastic: an L{ElasticLinear}-object
    @ivar yieldSurf: a SinglePlaneMohrCoulomb-object
    @ivar dilation: the d-parameter of type float
    """

    # _fixedMatIndex unchanged
    # _fixedSampleIndex unchanged
    # _varSampleIndex unchanged
    _varSampleIndex = dict(pst=0, G=1, K=2, cohesion=3,
                           friction=4, dilation=5)
    _yieldType = SinglePlaneMohrCoulomb

    _conversionRead = dict()

    _conversionWrite= dict()

    UCSi, e, a = -1., 1., 1.

    def fromHoekBrown(self, *args, **kwargs):
        raise NotImplementedError('SinglePlaneMohrCoulomb-Material has no' +
                                  ' method to be dirived from HoekBrown')

#} #End Plastic-Strain-Samples
###############################################################################

###############################################################################
#{ Errors
class UmatLookUpError(Exception):
    """Throw this on Errors durring parsing of UmatData"""
    pass
#} #End Errors

###############################################################################
#{ PST-dependent Plasticity
class DefaultPst(list):
    """Base-Class for plastic-strain dependent plasticity

    This describes the elasto-plastic behaviour of a particular material:
    The interpolation of several yield surfaces, elastic properties and
    dilation.

    It's basically a list of PST-sample-objects (of type L{SampleEMW2010}
    or similar) in ascending order of pst values (which is not yet enforced
    technically).
    """
    def __init__(self, *args, **kwargs):
        """
        Accepts an iterable of (pst, elast, yield, dilation)-tuples, where:
         - pst is the plastic strain value for this sample (lower boundary),
         - elast is an L{ElasticLinear}-object,
         - yield is a L{yieldsurfaces.ExtendedMenetreyWillam}-object and
         - dilation is the d-parameter of type float.

        Also accepts a PlasticLevkovitchReusch object. In this case it acts as
        a (deep) copy constructor.

        You can pass additional material parameters (like GSI) as
        keyWordArguments. They will be stored as attributes and can be
        accessed later on. But be carefull to not overwrite any functions or
        required parameters.

        Usage:
         >>> from bae.material_01.yieldsurfaces import ExtendedMenetreyWillam
         >>> from bae.material_01 import PlasticLevkovitchReusch, ElasticLinear
         >>> UCSi = 50E6
         >>> dilation = .1
         >>> elast = ElasticLinear(E = 10E9, nu =.2)
         >>> ySurfPeak = ExtendedMenetreyWillam(s=6.E-4, mb=.9, UCSi=UCSi)
         >>> ySurfRes  = ExtendedMenetreyWillam(s=3.E-4, mb=.8, UCSi=UCSi)
         >>> pLR = PlasticLevkovitchReusch([[0, elast, ySurfPeak, dilation],
         ...                                [0.03, elast, ySurfRes, dilation]])
         >>> # create a copy
         >>> pLR2 = PlasticLevkovitchReusch(pLR)

        """
        self.GSI = None    # store an unset GSI for convenience
        #self.__UCSi = None # see self.UCSi-property/setter

        for key, value in kwargs.iteritems():
            setattr(self, key, value)

        if not args:
            return list.__init__(self)

        if len(args) > 1 or not hasattr(args[0], '__iter__'):
            raise ValueError('A PST-List-object excepts up to one positional'
                             ' arguments. That could be a PST-List-object,'
                             ' a sequence of pst samples, or a sequence'
                             ' of (pst, elast, yieldsurf, dilation)-tuples.')
        arg = args[0]

        yieldType = self._sampleType._yieldType
        # args is PST-SampleList --> copy constructor
        if hasattr(arg, '_sampleType'):
            list.__init__(self, [self._sampleType(sample.pst,
                                                  sample.elastic.copy(),
                                                  yieldType(sample.yieldSurf),
                                                  sample.dilation)
                                 for sample in args[0]])
            self.GSI = getattr(args[0], 'GSI', None)
            return

        list.__init__(self)
        # args is list of pst-samples or quadtuples
        for x in arg:
            if hasattr(x, '_yieldType'): #args is list of samples
                self.append(self._sampleType(x.pst, x.elastic.copy(),
                                             yieldType(x.yieldSurf),
                                             x.dilation)) #append pst-sample copy

            elif hasattr(x, "__iter__"): #args is list quadtuples
                self.append(self._sampleType(*x))

            else:
                raise ValueError("Can't convert %s into a PST-sample" % str(x))


    def copy(self):
        """Creates a deep(copy) from PST-Samples-List (same as copy
        constructor)"""
        return self.__class__(self) # use copy-constructor

    def toDict(self):
        """Stores samples and GSI to a dictionary. The samples item
        will be a list of parameter-dicts for each pst-sample (see L{SampleEMW2010.toDict}).
        """
        GSI = getattr(self, 'GSI', None)
        return dict(GSI=GSI, samples=[s.toDict() for s in self])

    @classmethod
    def fromDict(cls, input):
        """Creates a PST-Samples-List instance from a dictionary. The dictionary
        has to have a 'samples'-item that is a list of parameter-dicts for each
        pst-sample (see L{toDict}).
        """
        GSI = input.get('GSI', None)
        samples = input['samples']
        return cls([cls._sampleType.fromDict(s) for s in samples], GSI=GSI)

    def __getattribute__(self, propName):
        """This __getattribute__ function allows to get yieldsurface properties
        that has to be the same for all samples directly by calling
            >>> myPstList.<propName>
        where <propName> could be e.g. UCSi, a, e, anisoN, ... .
        It first checks if the requested property is in _fixedSampleIndex of
        self._sampleType. If so, the property value is taken from the
        yieldSurface of first item in self. Otherwise, the 'traditional'
        __getattribute__ is called.
        """
        # to avoid a recursion error we'll need to call self._sampleType,
        # which calls __getattribute__ itself, in old fashioned way.

        _sampleType = super(DefaultPst, self).__getattribute__('_sampleType')

        if propName in _sampleType._fixedSampleIndex.keys():
            #property has to be the same in all sample yieldsurfaces
            valsInSamples = [getattr(s.yieldSurf, propName) for s in self]
            if any([s is None for s in valsInSamples]):
                raise ValueError('one or more propertys %s is None' %
                                 propName)
            try:
                if any([v-valsInSamples[0] for v in valsInSamples]):
                    raise ValueError("%s-values differ over samples (%s)"
                                     % (propName, valsInSamples))
                return valsInSamples[0]
            except:
                msg('variing props for %s:',propName)
                msg(str(valsInSamples))
                raise

        return super(DefaultPst, self).__getattribute__(propName)


    def __setattr__(self, propName, value):
        """This __setattr__ function allows to set yieldsurface properties
        that has to be the same for all samples directly by using
            >>> myPstList.<propName> = value
        where <propName> could be e.g. UCSi, a, e, anisoN, ... . The
        new value will get set for all yieldsurfaces of items in self.
        """
        if propName in self._sampleType._fixedSampleIndex.keys():
            for sample in self:
                setattr(sample.yieldSurf, propName, value)
        else:
            super(DefaultPst, self).__setattr__(propName, value)


    def _getProp(self, prop, idx=None):
        """Returns a stored property.

        @param prop: name of stored property (attribute)

        @param idx: index of plastic strain sample to take the
            property from. If idx=None, the property will be
            assuemd to be a static/not-pst-dependend one

        @raises ValueError: if prop is not found
        """
        try:
            return getattr(self, prop)
        except:
            pass

        if idx is not None:
            try:
                sample = self[idx]
            except IndexError:
                raise IndexError("Material has only %d pst-samples. " %
                                 len(self) + "You asked for index %d" % idx)
            try:
                return getattr(sample, prop)
            except:
                pass

            try:
                return getattr(sample.elastic, prop)
            except:
                pass

            try:
                return getattr(sample.yieldSurf, prop)
            except:
                pass

            raise ValueError("Can't find Property %s in pst-sample of type %s."
                             % (prop, str( type(sample) )) )

        else:
            raise ValueError("Property %s is not a attribute of %s itself. " %
                             (prop, self.__class__.__name__) +
                             "Maybe you need to specify an pst-sample-index.")


    def getPstDependentParameter(self, prop, pst):
        """Returns the linear interpolated parameter or complete yieldsurface
        for a specific pst value or a series of pst values.

        @param prop: name of stored property (attribute) or 'yieldSurf'. Note
            that, do to performance issues, requesting a yieldSurf-object is
            not supported for a series of pst-values
        @param pst: (series of) plastic strain value(s). Will be clipped below
            zero.
        """
        if hasattr(pst, '__iter__'):
            if 'yield' in prop.lower():
                raise ValueError('Requesting a "yieldSurf"-object for a series'
                                 ' of pst values is not supported.')
            res = self._getPstDependentParameterVector(prop, pst)
            return res

        else:
            return self._getPstDependentParameterSingle(prop, pst)


    def _getPstDependentParameterVector(self, prop, pst):
        """Vectorized, numpy-based part of L{getPstDependentParameter}
        used for iterable pst-values.
        """

        if not self:
            return np.full_like(pst, np.nan)
        pst = np.asarray(pst)
        pstSamples = [sample.pst for sample in self]
        valueSamples = [self._getProp(prop, ii) for ii in range(len(self))]

        # np.interp returns 'last valid' sample value outside of ranges by
        # default
        return np.interp(pst, pstSamples, valueSamples)


    def _getPstDependentParameterSingle(self, prop, pst):
        """Straight forward part  L{getPstDependentParameter} used for single
        pst-values.
        """
        pst = max(pst, 0)
        pstII = [ii for ii,sample in enumerate(self) if sample.pst <= pst]

        if not pstII:
            raise ValueError("Can't find a pst value (equal)below %f in"
                             " pst-sample list" % pst)

        iiL = pstII[-1] #lower pst sample index
        ## if iiL is last sample or pst matches the sample(edge) perfectly
        ## return the sample-parameter
        if (iiL == len(self) - 1) or (self[iiL].pst == pst):
            if 'yield' in prop.lower():
                return self[iiL].yieldSurf
            else:
                return self._getProp(prop, iiL)

        iiU = iiL + 1 #upper pst sample index

        frac = (pst - self[iiL].pst) / (self[iiU].pst - self[iiL].pst)

        if 'yield' in prop.lower():
            props = self[0]._yieldType._yieldParameters
        else:
            props = [prop,]
        interPols = []
        for pp in props:
            valueL = self._getProp(pp, iiL)
            valueU = self._getProp(pp, iiU)
            if valueL == valueU: # no interpolation required
                interPols.append(valueL)
            else:
                interPols.append( (1-frac) * valueL + frac * valueU )

        if 'yield' in prop.lower():
            return self[0]._yieldType(**dict(zip(props, interPols)))
        else:
            return interPols[0]


    def getDamageVariable(self, pst):
        """Returns the damage value which is defined as the following unified
        plastic strain (PST) dependent ratio:
            damVar = (UCS(PST) - UCS_Res) / (UCS_Peak - UCS_Res) ,
        where the UCS is given as:
            UCS = UCSi*s**a .
        Thus, this value ranges from 1 (current UCS has Peak-level) to zero
        (current UCS has Residual-level). If the residal UCS equals the peak
        UCS, e.g. just one pst-sample is given, the ratio is set to 1.

        @param pst: plastic strain. Can be either a single value or a series
            of length N.
        @returns: damage variable as a scalar if pst was scalar or as
            a numpy array of the same size of pst
        """
        try:
            ysPeak = self[0].yieldSurf
            ysRes  = self[-1].yieldSurf
        except IndexError:
            raise IndexError('No PST-samples stored yet.')

        # getting UCS for peak and resudial level
        UCSPeak = ysPeak.UCSi * ysPeak.s**ysPeak.a
        UCSRes  = ysRes.UCSi * ysRes.s**ysRes.a

        # getting UCS for current level
        UCSi = self.getPstDependentParameter('UCSi', pst)
        s = self.getPstDependentParameter('s', pst)
        a = self.getPstDependentParameter('a', pst)

        # everything to numpy for unified handeling
        UCSi, s, a = np.atleast_1d(UCSi), np.atleast_1d(s), np.atleast_1d(a)
        UCS = UCSi * s**a

        UCSspan = UCSPeak - UCSRes
        if not UCSspan: # Peak equals Res
            damVar = np.ones_like(UCS)
        else:
            damVar = (UCS - UCSRes) / UCSspan

        try:
            pst[0]
        except TypeError:
            # scalar in > scalar out
            damVar = damVar[0]

        return damVar


    def _validateUmatData(self, umatData):
        """ Prelimary check if choosen umatDialect and Index-to-Value-LookUp
        fits the dataStructure in umatData. Only the number of values and the
        number of specified tuples is checked --> its not 100% watertight.
        """
        sample = self._sampleType()
        dialect      = sample._vumatDialect
        nFixedMat    = len(sample._fixedMatIndex)
        nFixedSample = len(sample._fixedSampleIndex)
        nPerSample   = len(sample._varSampleIndex)

        if dialect == 'vumat2010':
            # i) checkSampleNumbers
            nVarData = len(umatData) - (2*nFixedSample + nFixedMat + 2)
            nTup, rest = divmod(nVarData, nPerSample)

            if (nTup < 2):
                raise UmatLookUpError(
                    "Found only %d pst-tuples. At least one for Fill and one"
                    " for Rock is required." % nTup)

            if not rest == 0:
                raise UmatLookUpError('Number of pst dependend values'
                                      ' does not fit exactly N*nPSTTuple')

            # ii) check indexPositions
            # number of pst-tuples for rock
            iRockCnt = 2*nFixedSample + nFixedMat
            rockCnt = umatData[iRockCnt]
            if ( (abs(rockCnt - int(rockCnt)) > 1E-12) or
                 (rockCnt > nTup) or (int(rockCnt) == 0) ):
                raise UmatLookUpError(
                    'Value %f (index %d) is not a valid value for '
                    'number of Rock-PST-tuples.'
                    % (rockCnt,iRockCnt))

            # number of pst-tuples for fill
            iFillCnt = iRockCnt + int(rockCnt)*nPerSample + 1
            fillCnt = umatData[iFillCnt]
            if ( (abs(fillCnt - int(fillCnt)) > 1E-12) or
                 (rockCnt > nTup) or (int(rockCnt) == 0) ):
                raise UmatLookUpError(
                    'Value %f (index %d) is not a valid value for '
                    'number of Fill-PST-tuples.'
                    % (fillCnt, iFillCnt))

            return True

        else:
            raise ValueError('Do not understand vumatDialect %s' % dialect)


    def _evalUmatData(self, matState, umatData):
        """Interpreting the umatDataList using the value-to-index-LookUp
        specified in self.__init__
        """
        if not len(self) == 0:
            raise LookupError('Sample List is not empty!')

        self._validateUmatData(umatData)

        sample = self._sampleType()
        _fixedMatIndex      = sample._fixedMatIndex
        _fixedSampleIndex   = sample._fixedSampleIndex
        _varSampleIndex     = sample._varSampleIndex
        _conversionRead     = sample._conversionRead

        nFixedIdx  = 2*len(_fixedSampleIndex) + len(_fixedMatIndex)
        nPerSample = len(_varSampleIndex)

        # getting startIndex for fixed and pst-tupleData
        if matState == 'rock':
            idxA = 0
            idxB = nFixedIdx

        elif matState == 'fill':
            idxA = len(_fixedSampleIndex)
            idxB = (nFixedIdx + int(umatData[nFixedIdx])*nPerSample + 1)

        #number of pst-tuples
        nSamples = int(umatData[idxB])
        #print 'nSamples:', nSamples

        elastPropKeys = ['G', 'K']
        plastPropKeys = [k for k in _varSampleIndex
                         if k not in elastPropKeys]

        def evalUmat(varKey, idx, offset = 0):
            # get associated umatData from index, convert if needed
            #get convert function - if not existent fun just returns arg
            fun = _conversionRead.get(key, lambda x: x)
            return fun(umatData[idx + offset])

        # looping all pst-samples and append a sample to self
        for kk in range(nSamples):
            sample = self._sampleType()
            elastDict, plastDict, toTestDict = {}, {}, {}

            # elastic behavior
            for key in elastPropKeys:
                idx = idxB + 1 + kk*nPerSample + _varSampleIndex[key]
                elastDict[key] = evalUmat(key, idx)

            # plastic behavior
            for key in plastPropKeys:
                idx = idxB + kk*nPerSample + _varSampleIndex[key] + 1
                plastDict[key] = evalUmat(key, idx)

            for key, rIdx in _fixedSampleIndex.iteritems():
                plastDict[key] = evalUmat(key, rIdx, offset = idxA)

            toTestDict.update(plastDict)
            toTestDict.update(elastDict)

            # assembling sample
            # done via setattr, maybe better use self._sampleType.__init__
            for key, val in plastDict.items():
                # add all plastic parameters as attributes that are not
                # directly related to the yieldsurface
                if key not in sample._yieldType._yieldParameters:
                    setattr(sample, key, val)
                    plastDict.pop(key, None)  # remove parameter

            setattr(sample, 'yieldSurf', sample._yieldType(**plastDict))
            setattr(sample, 'elastic', ElasticLinear(**elastDict))

            plastDict.update(elastDict)

            #Test if all specified _ceckUmatVal-functions are fullfilled
            for key, testFun in self._checkUmatVal.iteritems():
                if not testFun(toTestDict[key]):
                    raise UmatLookUpError(
                        "Property %s-%s (%s) doesn't fullfill the specified"
                        " _checkUmatVal-function"
                        % (matState, key, toTestDict[key]))

            self.append(sample)

        matSpecificData = dict( (key, umatData[idx])
                                for key, idx in _fixedMatIndex.iteritems())

        return matSpecificData


    def compare(self, otherSamples, tol=0.001, ignoreNotFound=False,
                verbose=True, firstMismatch=False, tryExact=True):
        """Compares self to a second object of the same type given as argument
        otherSample and lists the differences in a dictionary. All parameter
        pairs (pst, elastic, yieldsurf and dilation) are evaluated using a
        relative tolerance.

         - In a first step, the plastic strain levels are compared
         - In second step all matching plastic strain levels from step one are
           compared. This step is skipped if no pst-pairs were found.
         - In a third step all samples in the same position/index are
           compared. This step is skipped if all pst-values can be paired and
           the result would be the same as in step 2.

        The resulting dictionary holds the differences (as an OrderedDict) for
        each of the comparison steps. If both pLR-samples are identical within
        the relative tolerance then no deltaItems are stored in the resulting
        dicts.

        A deltaItem has a descriptive key and a tuple as value with the
        first tuple item holding the entity of self and the second tuple item
        holding the (differing) item of the compared pLR-SampleSequence.

        Hello Tobias: This section needs a better explanation:
        For the purpose of better readability comments (deltaItems with empty
        tuple as value) can be stored in dicts.

        Example:
         >>> pLRa = PlasticLevkovitchReusch()
         >>> pLRa.fromComputator(65E6, 50)
         >>> pLRb = PlasticLevkovitchReusch()
         >>> pLRb.fromComputator(100E6, 50)
         >>> # apply some changes to pLRb
         >>> # old=pLRb[1].pst*(1+2E-14) #to test close psts
         >>> # setattr(pLRb[2],'pst',old)
         >>> pLRb.pop() #remove last sample
         >>> # print differences
         >>> isEqual, deltas = pLRa.compare(pLRb, verbose=True)
         >>> for key,subs in deltas.iteritems():
         >>>     print '+++ %s +++'%key
         >>>     for sKey,sub in subs.iteritems():
         >>>         if not sub:
         >>>             print '> %s'%sKey
         >>>         else:
         >>>             print '%s: \tpLRa=%s \tpLRb=%s'%(sKey,str(sub[0]),
         ...                                                   str(sub[1]))

        @param otherSamples: PlasticLevkovitchReusch to be compared

        @param tol: relative tolerance for compared values, default = 0.001

        @param verbose: if True, extra comments (only keys, empty values)
            get stored in Dictionary

        @param ignoreNotFound: ignores if an attribute is not found in other
            sample (default = False)

        @returns: bool (true if equal), orderedDict (holding the differences)
        """
        from bae.future_01 import OrderedDict

        rDict = OrderedDict([('pst samples', OrderedDict()),
                             ('comparison by pst', OrderedDict()),
                             ('comparison by position', OrderedDict()),
                            ])

        ### I) compare PST-Levels
        pstRefs = OrderedDict([(s.pst, s) for s in self])
        pstCmps = OrderedDict([(s.pst, s) for s in otherSamples])
        # check variing length of samples
        if not len(pstRefs) == len(pstCmps):
            rDict['pst samples']['number of LR-samples'] = ( len(pstRefs),
                                                             len(pstCmps) )
            if firstMismatch:
                return False, rDict

        # if some pst-Values are very close together, the comparison might be
        # ambiguous. This will be captured in the following ugly lines

        # get ambiguous pstValues in reference
        pstCmpInRef = []
        for pstCmp in pstCmps:
            # try to find identical match
            matched = getValsInTol(pstRefs, pstCmp,
                                   tol=tol, tryExact=tryExact,)
            if len(matched) > 1:
                # found more then one match --> skip comparision by pst
                comment = ("In reference some pst-values are very close "+
                           "together. Comparison by pst will fail")
                rDict['pst samples'][comment] = ()
                rDict['pst samples']["close pst-values in reference"] \
                    = (pstRefs.keys(), pstCmp)

                pstCmpInRef = []
                break

            elif len(matched) == 1:
                # found match
                pstCmpInRef.append(matched[0])

            else:
                pass

        # get ambiguous pstValues in compared sample
        pstRefInCmp = []
        for pstRef in pstRefs:
            #try to find identical match
            f = getValsInTol(pstCmps, pstRef, tol=tol, tryExact=tryExact)
            if len(f) > 1:
                # found more then one match --> skip comparision by pst
                comment = ("In compared LR-samples some pst-values are very "+
                           "close together. Comparison by pst will fail")
                rDict['pst samples'][comment] = ()
                rDict['pst samples']["close pst-values in compared"] \
                    = (pstRefs.keys(), pstCmp)

                pstRefInCmp = []
                break

            elif len(f) == 1:
                # found match
                pstRefInCmp.append(f[0])

            else:
                pass

        # pst-values which can't be found
        pstRefNotInCmp = [pst for pst in pstCmps if pst not in pstRefInCmp]
        pstCmpNotInRef = [pst for pst in pstRefs if pst not in pstCmpInRef]

        if pstRefNotInCmp or pstCmpNotInRef:
            rDict['pst samples']['unmatched pst'] = (pstCmpNotInRef,
                                                     pstRefNotInCmp)
            if firstMismatch:
                return False, rDict

        ### II) found equal pst-values --> compare these pairs
        if pstCmpInRef and pstRefInCmp:
            for pstRef, pstCmp in zip(pstCmpInRef, pstRefInCmp):

                delta = pstRefs[pstRef].compare(pstCmps[pstCmp],
                                                tol=tol,
                                                ignoreNotFound=ignoreNotFound)

                if delta:
                    #found differences
                    comment = 'different values for matched pst=%.3f' % pstRef
                    rDict['comparison by pst'][comment] = ()

                    dd = dict(('pst=%.3f %s' % (pstRef, key),val)
                              for key, val in delta.iteritems())
                    rDict['comparison by pst'].update(dd)
                    if firstMismatch:
                        return False, rDict

            skipCompareByPos = False
            if len(pstRefs) == len(pstCmpInRef):
                # all LR-samples could be compared by pst value
                # --> no need to compare them by postition
                if rDict['comparison by pst']:
                    comment = ('All samples found by pst. ' +
                               'Skipped comparing by position.')
                    rDict['comparison by position'][comment] = ()
                    skipCompareByPos = True

        else:
            comment = "Can't match any pst. Comparison was skipped."
            rDict['comparison by pst'][comment] = ()
            skipCompareByPos = True

        ### III) compare samples by position
        if not skipCompareByPos:
            for ii, (pstRef, pstCmp) in enumerate(zip(pstRefs, pstCmps)):

                delta = pstRefs[pstRef].compare(pstCmps[pstCmp],
                                                tol=tol,
                                                ignoreNotFound=ignoreNotFound)

                if delta:
                    #found differences
                    comment = 'different values sample %d' % ii
                    rDict['comparison by position'][comment] = ()
                    dd = dict(('sample=%d %s' % (ii, key),val)
                              for key, val in delta.iteritems())
                    rDict['comparison by position'].update(dd)

        #clean up if no comments requested
        def remEmptyDict(dd):
            try:
                for key,val in dd.iteritems():
                    if val == ():
                        dd.pop(key,None)
                    else:
                        remEmptyDict(val)
            except AttributeError:
                pass

        if not verbose:
            remEmptyDict(rDict)

        isEqual = not any( map(bool,rDict.values()) )

        return isEqual, rDict


    def toVusubsVarFillMat(self):
        """Converts the PST-sample into a L{VusubsVarFillMat} instance.
        Example:
            >>> from bae.material_01 import PlasticLevkovitchReusch
            >>> myVarFillA = PlasticLevkovitchReusch().fromComputator(50E6,45)
            >>> myVarFillB = PlasticLevkovitchReusch().fromComputator(30E6,45)
            >>> print(    myVarFillA.toVusubsVarFillMat()
            ...       + 5*myVarFillB.toVusubsVarFillMat() )
        """
        return VusubsVarFillMat([self,])


class PlasticLevkovitchReusch(DefaultPst):
    """Describes the elasto-plastic behaviour of a particular material.

    The interpolation of several yield surfaces, elastic properties and
    dilation.

    It's basically a list of PST-sample-objects (of type L{SampleEMW2010}
    or similar) in ascending order of pst values (which is not yet enforced
    technically).

    It can be created
        - as an empty Object
        - as a computator material (L{fromComputator})
        - from a list of pst-samples (L{SampleEMW2010}, L{SampleAEMW2017} or
          L{SampleSPMC2010})
        - from a list of (pst, L{ElasticLinear}, yieldSurface, dilation)-tuples
          (see examples in L{DefaultPst.__init__})
    """
    ### sampleType
    _sampleType = SampleEMW2010

    _checkUmatVal = {'UCSi': lambda ucs: ucs>0, 'a': lambda a: 0.0<a<=1.}


    def fromComputator(self, UCSi, GSI, version='V14.3', **kwargs):
        """Creates a PlasticLevkovitchReusch sampleSequence using the
        Computator definitions.
        If sample-fixed values (a, e, (anisoN, ...) ) where not explicitly
        specified via kwargs, but do already exist, these will be passed to
        new material. Otherwise they'll get their default values specified in
        the associated yieldSurface-Class.

        Usage:
            creating a standard V14.3 computator material
             >>> pLR = PlasticLevkovitchReusch().fromComputator(65E6,50,
             ...                                                version='V14.3')

            setup/derive a user defined computator
             >>> from bae.material_01.models import computatorV14
             >>> # your own computator
             >>> def userComputator(UCSi, GSI, reduceE=1.):
             >>>     props = computatorV14(UCSi, GSI)
             >>>     if UCSi < 40E6 and GSI > 70:
             >>>         for pstSample in props.values():
             >>>             pstSample['E'] *= reduceE
             >>>     return props
             >>>
             >>> # now create a PLR-material
             >>> userPLR = PlasticLevkovitchReusch()
             >>> # note that all userfunction-specific arguments (like reduceE)
             >>> # will be passed via extraArgs-dictionay
             >>> userPLR.fromComputator(30E6, 70, version=userComputator,
             ...                        extraArgs={'reduceE':.75})

        @param UCSi: uniaxial compressive strength of intact rock (in MPa
        or Pa, see UCSisInMPa)

        @param GSI:  Geological strength index (in %)

        @param version: computator version (V12, V13, V14.3) string or userdefined
            computator function. If later, the user has to take care that
            the function returns same structure as the computor-functions defined
            in L{models} (ordered dict, each item a pst-sample, each sample holds
            a dict with all values required to define a plr-material). If you need
            to pass extra (keyword)arguments to the userdefined function you have
            to use the extraArgs-kwarg (see below).

        @kwarg a: Menetrey-Willam exponent. For default value see
            ExtendedMenetreyWillam

        @kwarg e: Menetrey-Willam excentricity. For default value see
            ExtendedMenetreyWillam

        @kwarg anisoS: (PlasticLevkovitchReuschAniso only) anisotropy constant

        @kwarg anisoN: (PlasticLevkovitchReuschAniso only) anisotropy constant

        @kwargs extraArgs: dictionary with abitrary fields that will be passed
            to the computator-function if a computator function is passed as
            version-argument

        @note: all units SI!!!

        """
        if isinstance(version, basestring):
            rawVersion = version
            version = version.strip().lower().replace('.','')
            if version == 'v12':
                pLRs = computatorV12(UCSi, GSI, returns='fullDict')
            elif version == 'v13':
                pLRs = computatorV13(UCSi, GSI, returns='fullDict')
            elif 'v14' in version:
                if version == 'v144':
                    subVersion = 4
                elif version == 'v143':
                    subVersion = 3 # default
                elif version == 'v142':
                    subVersion = 2 #
                elif version == 'v141':
                    subVersion = 1 #
                elif version == 'v140':
                    subVersion = 0
                else:
                    raise ValueError('Do not know subversion %s for computator V14'
                                    % version)
                cArgs = {'mi':kwargs.pop('mi', None),
                        'subVersion':subVersion,}
                if 'mr' in kwargs:
                    cArgs['mr'] = kwargs.pop('mr')
                cArgs['returns'] = 'fullDict'
                pLRs = computatorV14(UCSi, GSI, **cArgs)
            else:
                raise ValueError(
                    "Don't know computatorVersion %s (yet)." % rawVersion)

        elif hasattr(version, '__call__'):
            extraArgs = kwargs.get('extraArgs', {})
            pLRs = version(UCSi, GSI, **extraArgs)

        else:
            raise ValueError('Only strings or callables are allowed as'
                             ' version. You have passed: %s')


        self.GSI = GSI

        for ii in range(len(self))[::-1]:
            del self[ii]

        #pLRs comes as OrderedDict
        for pLR in pLRs.values():
            # assembling sample
            # done via setattr, maybe better use self._sampleType.__init__
            sample = self._sampleType()

            #set elastic behavior
            setattr(sample,'elastic', ElasticLinear(E=pLR['E'], nu=pLR['nu']))
            pLR.pop('E',None)
            pLR.pop('nu',None)

            pLR.update(dict((k,v) for k,v in kwargs.iteritems()
                            if k in sample._yieldType._yieldParameters))
            #set plastic behavior
            for key,val in pLR.copy().iteritems():
                # add all plastic parameters as attributes that are not
                # directly related to the yieldsurface
                if key not in sample._yieldType._yieldParameters:
                    setattr(sample,key,val)
                    pLR.pop(key, None)  # remove parameter

            setattr(sample, 'yieldSurf', sample._yieldType(**pLR))

            self.append(sample)

        for key, val in kwargs.iteritems():
            if hasattr(self, key):
                setattr(self, key, val)

        return self


    def isFromComputator(self, version=None, tol=0.001, firstMismatch=False):
        """Checks if the pLR-samples match a/any computatorVersion.

        Usage:
         >>> pLR = PlasticLevkovitchReusch().fromComputator(65E6,50,
         ...           version='V12')
         >>> print pLR.isFromComputator()
         >>> pLR.pop() #remove last sample
         >>> print pLR.isFromComputator()

        For computator versions V12 and V13 the GSI will be determined using
        the peak_mb. So the returned deltaDict referes to a computator material
        using this GSI.

        @param version: specifies the computator version(s) to be tested. If
            None, all implemented computator versions will be tested.

        @param tol: allowed relative valueTolerance

        @return: matched computatorVersion or None if none matched
            deltaDictionary which will only have empty subDicts if pLR matches
            a computator version. If more than one computator version were
            checked, the deltaDictionary refers to the last one checked
        """
        if len(self) == 0:
            msg("No pst-samples defined yet")
            return None, {}

        refValues = dict(UCSi = self[0].yieldSurf.UCSi,
                         mb   = self[0].yieldSurf.mb,
                         s    = self[0].yieldSurf.s,
                         E    = self[0].elastic.E,)

        if version is None:
            versions = ['V12','V13',] + ['V14.%d' % i for i in range(5)]
        elif isinstance(version,list):
            versions = version
        else:
            versions = [version,]

        found, deltas = False, {}
        for version in versions:
            # get GSI (V14: and optional input arguments) for computator
            inputs = restoreComputatorInputs(version, **refValues)
            if inputs['GSI'] is None:
                # can not restore computor input values
                # this should only happen for V14.1
                continue

            #by design we can reuse the inputs-dictionary to pass any other
            #recalculated kwargs for computators

            GSItest = inputs.pop('GSI') #required as positional argument
            UCSi    = refValues['UCSi'] #remains unchanged in material-object

            try: # extend with anisotropic parameters if existent
                inputs['anisoS'] = self[0].yieldSurf.anisoS
                inputs['anisoN'] = self[0].yieldSurf.anisoN
            except AttributeError:
                pass # is not of anisotropic type

            #create a generic LR-sequence using recalculated GSI
            computLR = type(self)(self)
            computLR.fromComputator(UCSi, GSItest, version=version, **inputs)

            # compare self with generic computator material
            found, deltas = self.compare(computLR, tol=0.001,
                                         ignoreNotFound=True,
                                         firstMismatch=firstMismatch)

            if found:
                # found matching computator version
                self.GSI = GSItest # store gsi
                break

        if found:
            foundVersion = version
        else:
            foundVersion = None

        return foundVersion, deltas


    def convert(self, newType=None, addGeneral={}, addYieldSurf={}):
        """DEPRECATED. Use copyConstructor of newType instead.
        Converts the elastoPlastic behavior to another (similar) one.
        If the new material type needs more parameters to be specified, these
        can be passed via addGeneral-dictionary (will create attributes to
        newMaterial itself) and addYieldSurf-dictionary (will add attributes to
        each new sample).

        Example:
            >>> pLR = PlasticLevkovitchReusch()
            >>> pLR.fromComputator(65E6, 50)
            >>> pLRAniso = pLR.convert(newType = PlasticLevkovitchReuschAniso,
            ...                        addYieldSurf = dict(anisoS=.1,
            ...                                            anisoN=.5,) )

        @param newType: new elastoPlastic material type

        @param addGeneral: dict holding attributes to be set to newMaterial

        @param addYieldSurf: dict holding attributes to be set to each sample

        @note: setting an attribute to each specific sample is not implemented
            yet.
        """
        if newType is None:
            newType = PlasticLevkovitchReuschAniso()
        else:
            newType = newType()

        for sample in self:
            # loop all existing samples and set properties to new
            props = addGeneral

            #copy general behavior
            for key,val in vars(sample).iteritems():
                if not key[0] == '_':
                    props[key] = val

            #set elastic behavior
            props['elastic'] = sample.elastic.copy()

            #plastic behavior
            yieldSurf = sample.yieldSurf.yieldParams()
            yieldSurf.update(addYieldSurf)
            props['yieldSurf'] = newType._sampleType._yieldType(**yieldSurf)

            #append new sample
            newSample = newType._sampleType(**props)
            newType.append(newSample)

        return newType


class PlasticLevkovitchReuschAniso(PlasticLevkovitchReusch):
    """This describes the elasto-plastic behaviour of a particular anisotropic
    material: The interpolation of several yield surfaces, elastic properties
    and dilation.

    It's basically a list of SampleAEMW2017-objects in
    ascending order of pst values (which is not yet enforced technically).

    It can be created
        - as an empty Object
        - as a computator material (L{fromComputator})
        - from a list of pst-samples (L{SampleEMW2010}, L{SampleAEMW2017} or
          L{SampleSPMC2010})
        - from a list of (pst, L{ElasticLinear}, yieldSurface, dilation)-tuples
          (see example in L{DefaultPst.__init__})

    """
    ### sampleType
    _sampleType = SampleAEMW2017

    def __init__(self, *args, **kwargs):
        super(PlasticLevkovitchReuschAniso, self).__init__(*args,**kwargs)
        if args and hasattr(args[0], 'anisoS'):
             self.anisoS = args[0].anisoS
             self.anisoN = args[0].anisoN
        else:
            self.anisoS = kwargs.pop('anisoS', 1)
            self.anisoN = kwargs.pop('anisoN', 1)


    def convert(self, newType=PlasticLevkovitchReusch, **kwargs):
        """DEPRECATED. Use copy constructor from different type instead.
        """
        return PlasticLevkovitchReusch.convert(self, newType=newType, **kwargs)


class SinglePlaneMohrCoulomb(DefaultPst):
    """This describes the elasto-plastic behaviour of a particular material:
    The interpolation of several yield surfaces, elastic properties and
    dilation.

    It's basically a list of PST-sample-objects in ascending
    order of pst values (which is not yet enforced technically).
    """

    ### sampleType
    _sampleType = SampleSPMC2010

    # the UCSi-value in Umat-Data will be negative to identify SPMC-material
    _checkUmatVal = {'UCSi': lambda ucs: ucs < 0}


    def __init__(self, *args, **kwargs):

        DefaultPst.__init__(self, *args, **kwargs)
        self.UCSi = -1.
        del self.GSI



class FactorOfSafetyIso(PlasticLevkovitchReusch):
    """This class holds some approaches to estimate (a local) estimate of a
    Factor Of Safety (FOS) for isotropic materials LR2-materials. The plastic
    strain dependency of a yield surface is evaluated before. So you derive a
    reasonable FOS even for yielded rockmass.
    """

    def _solveYieldSurfQLarge(self, pst, p, theta):
        """Vectorized method to get vonMises stress for a large set of
            pst-p-theta tripels. q on yieldsurfaces is solved iteratively
            using newton slover.
        """
        from scipy.optimize import newton
        ys = self._sampleType._yieldType

        # get strainDependent yield parameters
        params = dict((key, self.getPstDependentParameter(key,pst))
                      for key in ys._yieldParameters)

        # set start value. 2*p works fine but other values might be more
        # appropriate here
        qInter = np.abs(2*p)

        # define function to find root for
        rawYieldFun = ys._yieldSurfaceFun
        def vecYieldFun(q):
            return rawYieldFun(p,q,theta,**params)

        usePrime=True # toggle usage of first deviation here
        if usePrime:
            UCSi = params['UCSi']
            a    = params['a']
            mb   = params['mb']
            e    = params['e']
            ecc = ys.eccentricity(theta, e)
            def vecDYieldFun(q):
                # d(ys)/dq
                return (1./(UCSi*a))*(q/UCSi)**(1/a - 1) + mb*ecc/(3*UCSi)
            # slove root
            q = newton(vecYieldFun, qInter,fprime=vecDYieldFun)
        else:
            # solve root without first deviation
            q = newton(vecYieldFun, qInter)
        return q

    def fosMises(self, S1, S2, S3, pst, pstMax=None, miningConvention=False,
                 infValue=np.inf):
        """This FactorOfSafety uses the vonMises stress to compare the actual
        stress state (described by pressure p, vonMisesStress q and LodeAngle
        theta) with the related stress state of the yieldSurface at a given
        plastic strain.
        Here the actual vonMises stress and the vonMises stress on yield
        surfaces will be compared by:
            fos := q_yield / q
        where p_yield=p and theta_yield=theta is set.

        @param S1, S2, S3: priciple stresses in Pa. Can be either
            a single value or a series of length N.
        @param pst: plastic strain. Can be either a single value or a series
            of length N. Set to zero or N*[0,] to get fos for intact rock.
        @param pstMax: sets fos-values to 0 where pst exceeds pstMax
        @param miningConvention: if True, S1, S2 and S3 are given in mining
            convention (tension is negative).
        @param infValue: by definition, fos becomes infinite when S1 is zero.
            Set infValue to get a definied/fixed value.
        @return: a numpy array of size N. If S1, S3 and pst were given as
            single values, the numpy array will be of size (1,)

        @note: for this method q at yieldSurface needs to be solved
            iteratively, thus it might be slow(er) for many points.
        """
        from bae.material_01.yieldsurfaces import stressToPQTheta
        S1, S2, S3 = S1.copy(), S2.copy(), S3.copy()
        if miningConvention:
            # swap sign and order
            S1, S2, S3 = -S3, -S2, -S1

        p,q,theta = stressToPQTheta(S1,S2,S3)

        if hasattr(pst, '__iter__'):
            qYieldSurf = self._solveYieldSurfQLarge(pst, p, theta)
        else:
            # one yieldsuface for all points. No vectorization required.
            ys = self.getPstDependentParameter('yieldSurf', pst)
            qYieldSurf = ys.ptOnYieldSurf(p=p,theta=theta)

        qYieldSurf[np.isnan(qYieldSurf)] = 0
        if pstMax is not None:
            qYieldSurf[pst>pstMax] = 0

        fos = np.full_like(qYieldSurf, 0)

        iiInf = ~np.isclose(qYieldSurf,0) &  np.isclose(q,0)
        iiVal = ~np.isclose(qYieldSurf,0) & ~np.isclose(q,0)
        fos[iiInf] = infValue
        fos[iiVal] = qYieldSurf[iiVal] / q[iiVal]

        return fos

    def _solveYieldSurfPLarge(self, pst, q, theta):
        """Vectorized method to get pressure for a large set of
            pst-q-theta tripels. p can be solved directly.
        @returns p, tensile strength
        """
        ys = self._sampleType._yieldType

        # get strainDependent yield parameters
        params = dict((key, self.getPstDependentParameter(key,pst))
                      for key in ys._yieldParameters)

        UCSi = params['UCSi']
        a    = 1/params['a']
        mb   = params['mb']
        s    = params['s']
        e    = params['e']

        ii = q >= 0
        p = np.full_like(q, np.nan)
        ecc = ys.eccentricity(theta[ii], e)
        p[ii] = UCSi*((q[ii]/UCSi)**a - s)/mb + q[ii]*ecc/3.
        tensileStrength = -UCSi*s/mb

        return p, tensileStrength


    def fosHydrostatic(self, S1, S2, S3, pst, pstMax=None,
                       miningConvention=False, infValue=np.inf):
        """This FactorOfSafety uses the hydrostatic pressure to compare the
        actual stress state (described by pressure p, vonMisesStress q and
        LodeAngle theta) with the related stress stat of the yieldSurface at
        a given plastic strain.
        Here the actual pressure and the pressure on yield surfaces will be
        compared by:
            fos := (p_yield - tensileStress) / (p - tensileStress)
        where p_yield=p and theta_yield=theta is set and the tensile stress
        describes the lowes valid pressure of yieldSurface (appex).

        @param S1, S2, S3: priciple stresses in Pa. Can be either
            a single value or a series of length N.
        @param pst: plastic strain. Can be either a single value or a series
            of length N. Set to zero or N*[0,] to get fos for intact rock.
        @param pstMax: sets fos-values to 0 where pst exceeds pstMax
        @param miningConvention: if True, S1, S2 and S3 are given in mining
            convention (tension is negative).
        @param infValue: by definition, fos becomes infinite when S1 is zero.
            Set infValue to get a definied/fixed value.
        @return: a numpy array of size N. If S1, S3 and pst were given as
            single values, the numpy array will be of size (1,)
        """
        from bae.material_01.yieldsurfaces import stressToPQTheta

        S1, S2, S3 = S1.copy(), S2.copy(), S3.copy()
        if miningConvention:
            # swap sign and order
            S1, S2, S3 = -S3, -S2, -S1

        p,q,theta = stressToPQTheta(S1,S2,S3,miningConvention=False)

        if hasattr(pst, '__iter__'):
            pYieldSurf,tensStrength = self._solveYieldSurfPLarge(pst,q,theta)
        else:
            ys = self.getPstDependentParameter('yieldSurf', pst)
            pYieldSurf = ys.ptOnYieldSurf(q=q,theta=theta)
            tensStrength = ys.getIsoTensileStrength()

        p -= tensStrength
        pYieldSurf -= tensStrength

        fos = np.full_like(pYieldSurf, 0)

        iiInf =  np.isclose(pYieldSurf,0) & ~np.isclose(p,0)
        iiVal = ~np.isclose(pYieldSurf,0) & ~np.isclose(p,0)
        fos[iiInf] = infValue
        fos[iiVal] = p[iiVal] / pYieldSurf[iiVal]

        if pstMax is not None:
            fos[pst>pstMax] = 0

        return fos


    def fosHoekBrown(self, S1, S3, pst, pstMax=None, miningConvention=True,
                  infValue=np.inf, measureFromBisection=True):
        """Uses Hoek-Brown failure curve:
            sMajorYield = sMinor + UCSi*(mb*sMinor/UCSi + s)**a
        to estimate the Factor Of Safety:
            fos = sMajorYield / S1 .
        The pst-dependency of s and mb will be evaluated first.

        Note that Hoek-Brown yield surface is identical to MeneteryWillam in
        compression meridian (phi=pi/6) and 'below' elsewhere. Thus this FOS
        is a conservative estimate for the True FOS.

        @param S1: major priciple stress in Pa. Can be either a single value
            or a series of length N.
        @param S3: minor priciple stress in Pa. Can be either a single value
            or a series of length N.
        @param pst: plastic strain. Can be either a single value or a series
            of length N. Set to zero or N*[0,] to get fos for intact rock.
        @param pstMax: sets fos-values to 0 where pst exceeds pstMax
        @param miningConvention: if True, S1 and S3 are given in mining
            convention (tension is negative).
        @param infValue: by definition, fos becomes infinite when S1 is zero.
            Set infValue to get a definied/fixed value.
        @return: a numpy array of size N. If S1, S3 and pst were given as
            single values, the numpy array will be of size (1,)
        """
        S1,S3 = S1.copy(), S3.copy()
        if not miningConvention:
            # swap sign and order
            S1, S3 = -S3.copy(), -S1.copy()

        # pst-constant yieldSurf properties
        a = self[0].yieldSurf.a
        UCSi = self[0].yieldSurf.UCSi

        # eventually pst-variable yieldSurf properties
        s  = self.getPstDependentParameter('s', pst)

        mb = self.getPstDependentParameter('mb', pst)

        sMajorYield = self[0].yieldSurf._hoekBrownYieldSurf(S3, UCSi, s, mb, a)

        # sMajorYield == nan --> p is below tensile strength
        sMajorYield[np.isnan(sMajorYield)] = 0
        if pstMax is not None:
            sMajorYield[pst>pstMax] = 0

        fos = np.full_like(sMajorYield, 0)

        if measureFromBisection:
            sMajorYield -= S3
            S1 -= S3

        iiInf = ~np.isclose(sMajorYield,0) &  np.isclose(S1,0)
        iiVal = ~np.isclose(sMajorYield,0) & ~np.isclose(S1,0)
        fos[iiInf] = infValue
        fos[iiVal] = sMajorYield[iiVal] / S1[iiVal]

        return fos

    def fosMohrCoulomb(self, S1, S3, pst, pstMax=None, infValue=np.inf,
                       miningConvention=True, c=None, phi=None,
                       constrainCohesionToUCS=False, fittingRange=[0,None],
                       measureFromBisection=True):
        """Uses the MohrCoulomb criterion:
            sMajorYield = c + sMinor*tan(45 +.5*phi)**2
        to estimate the Factor Of Safety. Here, the cohesion c and friction
        angle phi only represent a MeneteryWillem yield surface only for
        a = 1. So, if not given explicitly, these values will be derived using
        a least-squares fit from the actual (eMW-)yieldSurface in compression
        meridian (theta=60deg). The plastic strain dependency of yield surface
        will be evaluated before.

        @param S1: major priciple stress in Pa. Can be either a single value
            or a series of length N.
        @param S3: minor priciple stress in Pa. Can be either a single value
            or a series of length N.
        @param pst: plastic strain. Can be either a single value or a series
            of length N. Set to zero to get fos for intact rock.
        @param pstMax: sets fos-values to 0 where pst exceeds pstMax
        @param infValue: by definition, fos becomes infinite when S3 is zero.
            Set infValue to get a definied/fixed value.
        @param miningConvention: if True, S1 and S3 are given in mining
            convention (tension is negative).
        @param c: cohesion in Pa. If None and constrainCohesionToUCS is False,
            the cohesion from best linear fit is used. If None and
            constrainCohesionToUCS is True, the UCS=(s**a)*UCSi is used for
            cohesion.
        @param phi: friction angle in deg. If None the friction angle from
            best linear fit is used.
        @param constrainCohesionToUCS: if True fitting will be forced to match
            UCS=(s**a)*UCSi for S3=0.
        @param fittingRange: range of p in Pa where the linear fit of
            MeneteryWillem yield surface will be applied
        @return: a numpy array of size N. If S1, S3 and pst were given as
            single values, the numpy array will be of size (1,)
        @note: Because of the fitting procedures and the unsufficent
            vectorization of this function it might become quite inefficient
            for large sets with varying pst values (e.g. from odb-sampling).

        """
        S1, S3 = S1.copy(), S3.copy()
        if not miningConvention:
            # swap sign and order
            S1, S3 = -S3, -S1

        pst = np.asarray(pst)
        if c is None or phi is None:
            # we need to evaluate the pst-dependency
            # pst-constant yieldSurf properties
            a = self[0].yieldSurf.a
            UCSi = self[0].yieldSurf.UCSi
            e = self[0].yieldSurf.e
            # eventually pst-variable yieldSurf properties
            s  = self.getPstDependentParameter('s', pst)
            mb = self.getPstDependentParameter('mb', pst)

            s, mb = np.atleast_1d(s), np.atleast_1d(mb)

            if c is None and constrainCohesionToUCS:
                c = UCSi*(s**a)

            if c is None or phi is None:
                def getMCparam(ls, lmb, lc):
                    eMW = ExtendedMenetreyWillam(UCSi=UCSi, s=ls,
                                                 mb=lmb, e=e, a=a)
                    cc,pp = QuasiMohrCoulomb.\
                            cohsFrictionfromExtendedMenetreyWillamFit(
                                eMW, fittingRange=fittingRange,
                                constrainedCohesion=lc,
                                nSamples=64, returnFit=False)
                    return cc,pp

                try:
                    tmp = [getMCparam(*x) for x in izip(s, mb, c)]
                except TypeError: #cohesion is single value
                    tmp = [getMCparam(*x) for x in izip(s, mb, len(s)*[c,])]

                if c is None:
                    c = np.asarray([aa for aa,_ in tmp])
                if phi is None:
                    phi = np.asarray([bb for _,bb in tmp])

        # tensile strength:
        ## HoekBrown for S1==S3 >> S3 = -s*UCSi/mb
        ## where s=1, UCSi=c and mb=2sin(phi)/(1-sin(phi))
        sinPhi = np.sin(np.deg2rad(phi))
        mb = (2*sinPhi)/(1-sinPhi)
        tensStrength = -c/mb
        tensFilter = (S3<=tensStrength)

        #sMajorYield = c + S3*(np.tan(np.deg2rad(45 +.5*phi)))**2
        sMajorYield = (1+mb)*S3 + c
        sMajorYield[tensFilter] = 0
        if pstMax is not None:
            sMajorYield[pst>pstMax] = 0

        fos = np.full_like(sMajorYield, 0)

        if measureFromBisection:
            sMajorYield -= S3
            S1 -= S3

        iiInf = ~np.isclose(sMajorYield,0) &  np.isclose(S1,0)
        iiVal = ~np.isclose(sMajorYield,0) & ~np.isclose(S1,0)
        fos[iiInf] = infValue
        fos[iiVal] = sMajorYield[iiVal] / S1[iiVal]
        return fos


class FactorOfSafetyAniso(FactorOfSafetyIso, PlasticLevkovitchReuschAniso):
    # Subclassed from PlasticLevkovitchReuschAniso as second so all
    # specifically anisotropic functions/attributes get updated/overwritten
    """Same as L{FactorOfSafetyIso} but derived from
    L{PlasticLevkovitchReuschAniso}.
    """
    pass
#} #End pst-dependend plasticity
###############################################################################

###############################################################################
#{IO to abaqus-models
class AbqVumatData(object):
    """Holds all data that is contained in a user material for our usual
    Abaqus material.

     - 2 PlasticLevkovitchReusch objects, one for intact rock and one for fill
     - a void ratio
     - critical damage
     - damping - alpha
     - density
     - umatDelete and umatDepvarNum parameter
     - an optional name

    @ivar name: name as in the Abaqus input file
    @ivar vumatType: "pLR" (default,L{PlasticLevkovitchReusch}),
        "pLRAniso" (L{PlasticLevkovitchReuschAniso}),
        "sPMC" (L{SinglePlaneMohrCoulomb})...

    @ivar density: value from the *DENSITY option
    @ivar dampingAlpha: argument of the *DAMPING, ALPHA=... option
    @ivar umatDelete: argument of the *DEPVAR,DELETE=... option: index
        (Fortran style, one-based) of the SDV indicating element deletion
    @ivar umatDepvarNum: value from the data line of the *DEPVAR option:
        number of SDV variables to be used.
    @ivar voidStiffnessRatio: PROPS(7) ELASTIC DEGRADATION RATIO
    @ivar critDamage: PROPS(8) CRITICAL DAMAGE

    @ivar plasticLRRock: PST-Sample object for rock of type
        L{PlasticLevkovitchReusch}, L{PlasticLevkovitchReuschAniso} or
        L{SinglePlaneMohrCoulomb}
    @ivar plasticLRFill: PST-Sample object for fill of type
        L{PlasticLevkovitchReusch} et al.
    """

    _implementedTypes = {'plr'      : PlasticLevkovitchReusch,
                         'plraniso' : PlasticLevkovitchReuschAniso,
                         'spmc'     : SinglePlaneMohrCoulomb,}

    def __init__(self, *args, **kwargs):
        """
        The positional argument might be a
        L{bae.abq_model_02.internal.Material<bae.abq_model_02.container.Material>}-object

        @kwarg inputFile: filename or other suitable argument to
             abq_model_02.Model.read(), requires keywordarg materialName

        @kwarg materialName: material name to read from inputFile, requires
            keywordarg inputFile

        @kwarg vumatType: "pLR" (default,L{PlasticLevkovitchReusch}),
            "pLRAniso" (L{PlasticLevkovitchReuschAniso}),...
        """
        if ('inputFile' in kwargs) and ('materialName' in kwargs):
            from bae.abq_model_02 import Model
            model = Model().read(kwargs['inputFile'])
#            print '>>>', model.material.keys()
#            print '>>>', kwargs['materialName']
#            print '>>>', model.material[kwargs['materialName']]
            mat = model.material[kwargs['materialName']]
            kwargs.pop('inputFile',None)
            kwargs.pop('materialName',None)
            self.__init__(mat, **kwargs)
            return

        elif ('inputFile' in kwargs) or ('materialName' in kwargs):
            raise ValueError('You have to specifiy materialName AND inputFile')


        if args:
            abqMat = args[0]
            # name
            try:
                self.name = abqMat.name
            except AttributeError:
                msg("WARNING: material name missing.")
                pass

            vumatType = kwargs.get("vumatType", None)

            # vumatType from kwargs
            if vumatType is None:
                # autoDetect vumatType from _implementedTypes
                errorTxt = '___UmatLookUpErrors' + 31*'_' + '\n'
                for vumatType in self._implementedTypes:
                    kwargs["vumatType"] = vumatType
                    try:
                        self = self.__init__(abqMat, **kwargs)
                        return
                    except UmatLookUpError as e:
                        errorTxt += ('*Parsing %s-Material failed:\n\t%s\n'
                                     % (vumatType, e))
                        pass
                # tried all implemented vumatTypes --> Error
                errorTxt += 50*'_'
                print errorTxt
                raise UmatLookUpError("Can't match the UmatData to any" +
                                      " of the implemented vumatTypes (%s).\n"
                                      " See parsingErrors above ErrorTraceback."
                                      % str(self._implementedTypes.values())
                                      )

            elif vumatType.lower() not in self._implementedTypes.keys():
                print vumatType.lower(), self._implementedTypes.keys()
                raise ValueError(("Do not know a vumatType %s " +
                                  "(allready implemented: %s)") %
                                 (kwargs["vumatType"],
                                  str(self._implementedTypes)) )

            else:
                self.vumatType = vumatType.lower()


            # density: value from the *DENSITY option
            try:
                self.density = abqMat["Density"]
            except KeyError:
                msg("WARNING: no density property.")

            # dampingAlpha: argument of the *DAMPING, ALPHA=... option
            try:
                self.dampingAlpha = abqMat["DampingAlpha"]
            except KeyError:
                msg("WARNING: no dampingAlpha property.")

            # umatDelete: argument of the *DEPVAR,DELETE=... option
            try:
                self.umatDelete = abqMat["UmatDelete"]
            except KeyError:
                msg("WARNING: no umatDelete property.")

            # umatDepvarNum: value from the data line of the *DEPVAR option
            try:
                self.umatDepvarNum = abqMat["UmatDepvarNum"]
            except KeyError:
                msg("WARNING: no umatDepvarNum property.")

            # values from the *USER MATERIAL option
            try:
                umatData = abqMat["Umat"]
            except KeyError:
                raise ValueError("No Umat values.")


            rock = self._implementedTypes[self.vumatType]()
            fill = self._implementedTypes[self.vumatType]()

            # umatData will be set to rock/fill inplace - the returned value of
            # evalUmatData holds materialinherent properties
            otherProps = rock._evalUmatData('rock',umatData)
            otherProps = fill._evalUmatData('fill',umatData)

            self.plasticLRRock = rock
            self.plasticLRFill = fill

            for key,val in otherProps.iteritems():
                setattr(self,key,val)


    def compare(self, otherMat, valTol = 0.01,):
        """Compares two abaqus-material containers.
        """
        from collections import OrderedDict

        if not isinstance(otherMat, AbqVumatData):
            raise ValueError('Can only compare data of Type %s' %
                             str(type(AbqVumatData)) )

        diffDict = OrderedDict()

        if self.vumatType == otherMat.vumatType:
            diffDict['type'] = {}
        else:
            diffDict['type'] = {'vumatType': (self.vumatType,
                                              otherMat.vumatType)}

        selfMatVars = dict((key, val) for key,val in vars(self).iteritems()
                           if not hasattr(val,'__iter__') )

        otherMatVars = dict((key, val) for key,val in vars(self).iteritems()
                            if not hasattr(val,'__iter__') )

        selfMatVars.pop('vumatType',None)
        otherMatVars.pop('vumatType',None)

        diffDict['general'] = deltaDict(selfMatVars, otherMatVars, tol=valTol)

        rockEqual, diffDict['rockProperties'] = self.plasticLRRock.compare(
            otherMat.plasticLRRock, tol=valTol)

        fillEqual, diffDict['fillProperties'] = self.plasticLRFill.compare(
            otherMat.plasticLRFill, tol=valTol)

        if not ( diffDict['type'] or any(diffDict['general'].values()) or
                 (not rockEqual) or (not fillEqual) ):
            isEqual = True
        else:
            isEqual = False

        return isEqual, diffDict


    def updateAbqModel(self, model):
        """
        Adds this data to the abq_model_02.Model object model.

        Usage:
         >>> model = Model()
         >>> myAbqVumatData.updateAbqModel(model)
         >>> model.write("myMaterial.inp")

        ... would create an Abaqus input file "myMaterial.inp" containing the
        data of the AbqVumatData-object myAbqVumatData.
        """
        def rDict(indexDict):
            # Return reversal Index-Value-Lookups of indexDict for umatData
            dd = {}
            for k,t in indexDict.iteritems():
                if isinstance(t,tuple):
                    # t = (index, convertFunRead, convertFunWrite)
                    ii = t[0]
                    k  = tuple([k] + list(t[1:]))
                else:
                    # t = index
                    ii = t
                dd[ii] = k
            return dd

        # reversal index-to-key-LookUp
        sampleType = self.plasticLRRock._sampleType
        fixedMatIndex    = rDict(sampleType._fixedMatIndex)
        fixedSampleIndex = rDict(sampleType._fixedSampleIndex)
        varSampleIndex   = rDict(sampleType._varSampleIndex)

        # shortcut convertWrite-functions-dict
        convertWrite = self.plasticLRRock._sampleType._conversionWrite

        def setUmatValue(uMat, holder, ii, key):
            """Sets the value of uMat[ii] to FUN(holder.key), where FUN is
            taken from convertWrite if existent.
            The holderObject itself and (if existent) its yieldSurf and elastic
            attribute will be evaluated to 'find' the proper keyValue.
            """
            #get convert function - if not existent fun just returns arg
            fun = convertWrite.get(key, lambda x: x)

            # evaluate holder itself
            try:
                uMat[ii] = fun( getattr(holder, key) )
                return uMat
            except AttributeError:
                pass

            # evaluate yieldSurf attribute
            try:
                uMat[ii] = fun( getattr(holder.yieldSurf, key) )
                return uMat
            except AttributeError:
                pass

            # eval elastic attribute
            try:
                uMat[ii] = fun( getattr(holder.elastic, key) )
                return uMat
            except AttributeError:
                pass

            # proper value not found yet --> Error
            if uMat[ii] is None:
                raise AttributeError("Can't find attribute %s in %s" %
                                         (key,str(holder)) )

        #i) create list of fixed material values (ucs,e,a,dilation,...)
        umatData = (2*len(fixedSampleIndex) + len(fixedMatIndex)) * [None]

        for iiRock,t in fixedSampleIndex.iteritems():
            iiFill = iiRock + len(fixedSampleIndex)
            umatData = setUmatValue(umatData, self.plasticLRRock[0], iiRock, t)
            umatData = setUmatValue(umatData, self.plasticLRFill[0], iiFill, t)

        for ii,t in fixedMatIndex.iteritems():
            umatData = setUmatValue(umatData, self, ii, t)

        #ii) create list for rockSamples
        rSamples = ( len(self.plasticLRRock)*len(varSampleIndex) + 1 )*[None]
        rSamples[0] = len(self.plasticLRRock)

        #loop over rocks pst-samples
        for kk, sample in enumerate(self.plasticLRRock):
            for ii,t in varSampleIndex.iteritems():
                jj = ii + kk*len(varSampleIndex) + 1  # index
                # +1 for sampleCount
                rSamples = setUmatValue(rSamples,sample,jj,t)

        #iii) create list for fillSamples
        fSamples = ( len(self.plasticLRFill)*len(varSampleIndex) + 1 )*[None]
        fSamples[0] = len(self.plasticLRFill)

        #loop over rocks pst-samples
        for kk, sample in enumerate(self.plasticLRFill):
            for ii,t in varSampleIndex.iteritems():
                jj = ii + kk*len(varSampleIndex) + 1  # index
                # +1 for sampleCount
                fSamples = setUmatValue(fSamples,sample,jj,t)

        umatData = umatData + rSamples + fSamples

        model.material.updateItem(self.name, {
                "DampingAlpha": self.dampingAlpha,
                "Density": self.density,
                'UmatDelete': self.umatDelete,
                'UmatDepvarNum': self.umatDepvarNum,
                'Umat': umatData,
                })
#} #End IO to abaqus-models



if __name__ == '__main__':
    print('No syntax Errors')

