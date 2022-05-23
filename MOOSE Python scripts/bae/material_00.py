"""classes for material properties.

IMPORTANT: This module is in ALPHA state. If you want to use it in production
scripts please make a local copy in your working directory. Because this
implementation is not stable yet. Incompatible changes are likely to happen.

Examples
========

Just reads an input file, interprets it an writes it back
 >>> from bae.abq_model_02 import Model
 >>> from material_00 import AbqVumatData
 >>> from bae.log_01 import log
 >>>
 >>> log.setDebugLevel(maxDebugLevel=9999)
 >>>
 >>> m = Model().read("../PERSE2013_material_M22.inp")
 >>>
 >>> material = dict(
 >>>     ( name, AbqVumatData(mat) )
 >>>     for name, mat in m.material.iteritems() )
 >>>
 >>> print sorted(material)
 >>>
 >>> m2 = Model()
 >>> for mat in material.itervalues():
 >>>     mat.updateAbqModel(m2)
 >>>
 >>> m2.write("copy_for_check.inp")

Second Example: Create an input file from material properties:
 >>> from bae.abq_model_02 import Model
 >>> from material_00 import AbqVumatData, PlasticLevkowitchReusch, \\
 >>>     PlasticLevkowitchReuschSample, ElasticLinear
 >>> from bae.log_01 import log
 >>>
 >>> log.setDebugLevel(maxDebugLevel=9999)
 >>>
 >>> # Felsic, UCS=80MPa, GSI=60,
 >>> # ... pe trans, res, cave: 0.005, 0.02, 0.22
 >>> # ... dilation=(0.2,0.2,0.072), s=1.73E-03, mb=1.34, E=14.93GPa, nu=0.3
 >>> # => mi=mb/exp((GSI-100)/28)
 >>> propsHB = dict( UCSi=80E6, GSI=60, mi=5.6)
 >>> elastic = ElasticLinear(E=14.93E9, nu=0.3)
 >>>
 >>> elaPlaPeak = PlasticLevkowitchReuschSample(
 >>>     pst=0.0, dilation=0.2).fromHoekBrown(**propsHB)
 >>> elaPlaPeak.elastic = elastic
 >>> elaPlaTrans = PlasticLevkowitchReuschSample(
 >>>     pst=0.005, dilation=0.2).fromHoekBrown(**propsHB)
 >>> elaPlaTrans.elastic = elastic
 >>>
 >>> propsHB["GSI"] = 35
 >>> elaPlaRes = PlasticLevkowitchReuschSample(
 >>>     pst=0.02, dilation=0.072).fromHoekBrown(**propsHB)
 >>> elaPlaRes.elastic = elastic
 >>>
 >>> propsHB["GSI"] = 25
 >>> elaPlaCave = PlasticLevkowitchReuschSample(
 >>>     pst=0.22, dilation=0.0).fromHoekBrown(**propsHB)
 >>> elaPlaCave.elastic = elastic
 >>>
 >>> felsic = PlasticLevkowitchReusch((
 >>>         elaPlaPeak, elaPlaTrans, elaPlaRes, elaPlaCave ))
 >>>
 >>> # Dolerite, UCS=80MPa, GSI=65,
 >>> # ... pe trans, res, cave: 0.005, 0.02, 0.22
 >>> # ... dilation=(0.2,0.2,0.072), s=1.73E-03, mb=1.45, E=14.93GPa, nu=0.3
 >>> # => mi=mb/exp((GSI-100)/28)
 >>> propsHB = dict( UCSi=80E6, GSI=65, mi=5.06 )
 >>> elastic = ElasticLinear(E=14.93E9, nu=0.3)
 >>>
 >>> elaPlaPeak = PlasticLevkowitchReuschSample(
 >>>     pst=0.0, dilation=0.2).fromHoekBrown(**propsHB)
 >>> elaPlaPeak.elastic = elastic
 >>> elaPlaTrans = PlasticLevkowitchReuschSample(
 >>>     pst=0.005, dilation=0.2).fromHoekBrown(**propsHB)
 >>> elaPlaTrans.elastic = elastic
 >>>
 >>> propsHB["GSI"] = 35
 >>> elaPlaRes = PlasticLevkowitchReuschSample(
 >>>     pst=0.02, dilation=0.072).fromHoekBrown(**propsHB)
 >>> elaPlaRes.elastic = elastic
 >>>
 >>> propsHB["GSI"] = 25
 >>> elaPlaCave = PlasticLevkowitchReuschSample(
 >>>     pst=0.22, dilation=0.0).fromHoekBrown(**propsHB)
 >>> elaPlaCave.elastic = elastic
 >>>
 >>> dolerite = PlasticLevkowitchReusch((
 >>>         elaPlaPeak, elaPlaTrans, elaPlaRes, elaPlaCave ))
 >>>
 >>> # Black_Shale, UCS=65MPa, GSI=40,
 >>> # ... pe trans, res, cave: 0.005, 0.03, 0.23
 >>> # ... dilation=(0.1625,0.163,0.0585),
 >>> # ... s=1.16E-03, mb=0.77, E=12.13GPa, nu=0.3
 >>> # => mi=mb/exp((GSI-100)/28)
 >>> propsHB = dict( UCSi=65E6, GSI=40, mi=6.56 )
 >>> elastic = ElasticLinear(E=12.13E9, nu=0.3)
 >>>
 >>> elaPlaPeak = PlasticLevkowitchReuschSample(
 >>>     pst=0.0, dilation=0.1625).fromHoekBrown(**propsHB)
 >>> elaPlaPeak.elastic = elastic
 >>> elaPlaTrans = PlasticLevkowitchReuschSample(
 >>>     pst=0.005, dilation=0.163).fromHoekBrown(**propsHB)
 >>> elaPlaTrans.elastic = elastic
 >>>
 >>> propsHB["GSI"] = 35
 >>> elaPlaRes = PlasticLevkowitchReuschSample(
 >>>     pst=0.03, dilation=0.0585).fromHoekBrown(**propsHB)
 >>> elaPlaRes.elastic = elastic
 >>>
 >>> propsHB["GSI"] = 25
 >>> elaPlaCave = PlasticLevkowitchReuschSample(
 >>>     pst=0.23, dilation=0.0).fromHoekBrown(**propsHB)
 >>> elaPlaCave.elastic = elastic
 >>>
 >>> blackshale = PlasticLevkowitchReusch((
 >>>         elaPlaPeak, elaPlaTrans, elaPlaRes, elaPlaCave ))
 >>>
 >>> # Faults, UCS=45MPa, GSI=40,
 >>> # ... pe trans, res: 0.005, 0.03
 >>> # ... dilation=(0.1125,0.113,0.0405),
 >>> # ... s=6.74E-04, mb=0.63, E=8.4GPa, nu=0.27
 >>> # => mi=mb/exp((GSI-100)/28); MR might be unreasonable but yields right E
 >>> propsHB = dict( UCSi=45E6, GSI=40, mi=5.37, MR=1169, nu=0.27 )
 >>>
 >>> elaPlaPeak = PlasticLevkowitchReuschSample(
 >>>     pst=0.0, dilation=0.2).fromHoekBrown(**propsHB)
 >>> elaPlaTrans = PlasticLevkowitchReuschSample(
 >>>     pst=0.005, dilation=0.2).fromHoekBrown(**propsHB)
 >>>
 >>> propsHB["GSI"] = 25
 >>> elaPlaRes = PlasticLevkowitchReuschSample(
 >>>     pst=0.03, dilation=0.072).fromHoekBrown(**propsHB)
 >>>
 >>> faults = PlasticLevkowitchReusch((
 >>>         elaPlaPeak, elaPlaTrans, elaPlaRes))
 >>>
 >>> # Fill ... find proper HB param
 >>> # ... then create python material objects
 >>> fill = PlasticLevkowitchReusch((
 >>>         elaPlaPeak, elaPlaTrans, elaPlaRes))
 >>>
 >>> model = Model()
 >>> for name, rock in [
 >>>         ("Felsic", felsic),
 >>>         ("Dolerite", dolerite),
 >>>         ("Black_Shale", blackshale),
 >>>         ]:
 >>>     mat = AbqVumatData()
 >>>     mat.name = name
 >>>     mat.vumatType = "6tuple"
 >>>     mat.density = 2700.0
 >>>     mat.dampingAlpha = 0.5
 >>>     mat.umatDelete = 9
 >>>     mat.umatDepvarNum = 16
 >>>     mat.voidStiffnessRatio = 1e-5
 >>>     mat.critDamage = 1000.0
 >>>     mat.plasticLRRock = rock
 >>>     mat.plasticLRFill = fill
 >>>     mat.updateAbqModel(model)
 >>> model.write("Check_material_M01.inp")

Third example: Create Material from Computator and use linear fill 
 >>> from material_00 import PlasticLevkowitchReusch as PLR
 >>> from bae.abq_model_02 import Model
 >>> from material_00 import AbqVumatData
 >>> 
 >>> #e.g. from csv:
 >>> matList = [#name,ucsi/MPa,gsi
 >>>            ('RockA',60E6,60),
 >>>            ('RockB',65E6,55),
 >>>            ('FaultC',50E6,30),
 >>>           ] 
 >>> 
 >>> def elasticFill(E=.2E7,nu=.2):
 >>>     #create dummy plastic and force to be linear 
 >>>     fill = PLR()
 >>>     dummyUCS,dummyGSI = 1E6,50
 >>>     fill.fromComputator(dummyUCS,dummyGSI)
 >>>     for sample in fill:
 >>>         sample.yieldSurf.UCSi = 10E7
 >>>         sample.elastic.E = E
 >>>         sample.elastic.nu = nu
 >>>     return fill
 >>> 
 >>> model = Model()
 >>> 
 >>> for name, UCSi, GSI in matList:
 >>>     mat = AbqVumatData()
 >>>     mat.name = name
 >>>     mat.vumatType = "6tuple"
 >>>     mat.density = 2700.0
 >>>     mat.dampingAlpha = 0.5
 >>>     mat.umatDelete = 9
 >>>     mat.umatDepvarNum = 16
 >>>     mat.voidStiffnessRatio = 1e-5
 >>>     mat.critDamage = 1000.0
 >>>     rock = PLR().fromComputator(UCSi,GSI,version='V12')
 >>>     mat.plasticLRRock = rock
 >>>     mat.plasticLRFill = elasticFill()
 >>>     mat.updateAbqModel(model)
 >>> model.write("Check_material_M01.inp")

@Note: All units SI!
"""

__version__ = "0.02"

__bugs__ = r"""
 - none (of course)
 - needs completion:
   . AbqVumatData.__init__() with filename and material argument
     reading from file (really needed? Or just a small howto in the docs?)
   . AbqVumatData.__init__() with contents arguments: name, type, density...
   . reverse-Hoek-Brown: derive GSI, mi, MR from s, mb, E
     can be done manually with Rocklab:
       1. UCS as is
       2. choose GSI such that s is right
       3. choose mi such that mb is right
       4. choose MR such that E (Erm) is right
"""

_version_history_ = r"""
0.01 GP: new
"""

from bae.abq_model_02.container import Material as AbqMaterial
from bae.log_01 import msg

from copy import deepcopy
from itertools import imap

class ElasticLinear(object):
    """

    @ivar E: Youngs modulus
    @ivar nu: Poissons ratio
    @ivar K: bulk modulus
    @ivar G: shear modulus

    @Note: All units SI!
    """

    def __init__(self, *args, **kwargs):
        """
        make sure nu != 0.5 (actually nu < 0.5)

        @kwarg E: Youngs modulus
        @kwarg nu: Poissons ratio
        @kwarg K: bulk modulus
        @kwarg G: shear modulus

        @Note: All units SI!
        """
        if "E" in kwargs and "nu" in kwargs:
            self.setEnu(kwargs["E"], kwargs["nu"])
        elif "G" in kwargs and "K" in kwargs:
            self.setGK(kwargs["G"], kwargs["K"])
        else:
            raise ValueError("Either E and nu or G and K have to be specified"
                             " for an ElasticLinear-object.")

    def setEnu(self, E, nu):
        self.E = float(E)
        self.nu = float(nu)
        self.G = E/(2.0*(1.0+nu))
        self.K = E/(3.0*(1.0-2.0*nu))

    def setGK(self, G, K):
        self.G = float(G)
        self.K = float(K)
        self.E = 9.0*K*G / (3.0*K + G)
        self.nu = (1.5*K - G) / (3.0*K + G)

    def fromHoekBrown(self, UCSi, GSI, MR, D=0.0, nu=0.2):
        """
        set Youngs modulus from Hoek-Brown formulas

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
        self.setEnu(self.E, self.nu)


class YieldSurfMenetreyWillam(object):
    """This decribes one yield surface.
    """
    def __init__(self, *args, **kwargs):
        """
        @kwarg UCSi: UCS for intact rock, possible synonym: UCS
        @kwarg e: excentricity, default 0.6
        @kwarg a: Hoek-Brown-exponent, default 0.5
        @kwarg s:
        @kwarg mb:
        @Note: All units SI!
        """
        if "UCSi" in kwargs:
            self.UCSi = float(kwargs["UCSi"])
        elif "UCS" in kwargs:
            self.UCSi = float(kwargs["UCS"])
        if "e" in kwargs:
            self.e = float(kwargs["e"])
        else:
            self.e = 0.6
        if "a" in kwargs:
            self.a = float(kwargs["a"])
        else:
            self.a = 0.5
        if "s" in kwargs:
            self.s = float(kwargs["s"])
        if "mb" in kwargs:
            self.mb = float(kwargs["mb"])

    def fromHoekBrown(self, UCSi, GSI, mi, D=0.0, a=None, e=0.6):
        """
        Usage:
         >>> yield = YieldSurfMenetreyWillam().fromHoekBrown(UCSi, GSI, mi)

        @param a: if not given will be computed by the Hoek-Brown formula
          based on GSI, otherwise usually: a=0.5

        @returns: self
        @Note: All units SI!
        @note: from "HOEK-BROWN FAILURE CRITERION - 2002 EDITION"
        """

        from math import exp

        self.UCSi = float(UCSi)
        self.e = float(e)
        if a is None:
            self.a = 0.5 + (exp(-GSI/15.0)-exp(-20.0/3.0)) / 6.0
        else:
            self.a = float(a)
        self.s = exp((GSI-100.0)/(9.0-3.0*D))
        self.mb = mi * exp((GSI-100.0)/(28.0-14.0*D))
        return self


class PlasticLevkowitchReuschSample(object):
    """This is one sample of the complete elasto-plastic behaviour of a
    particular material, i.e. elastic, plastic and dilation properties for one
    particular plastic strain (pst) value.

    @ivar pst: plastic strain (pst) value
    @ivar elastic: an ElasticLinear-object
    @ivar yieldSurf: a YieldSurfMenetreyWillam-object
    @ivar dilation: the d-parameter of type float
    """
    def __init__(self, pst=None, elastic=None, yieldSurf=None, dilation=None):
        """
        @param pst: plastic strain (pst) value
        @param elastic: an ElasticLinear-object
        @param yieldSurf: a YieldSurfMenetreyWillam-object
        @param dilation: the d-parameter of type float
        """
        if pst is not None:
            self.pst       = float(pst)
        if elastic is not None:
            self.elastic   = elastic
        if yieldSurf is not None:
            self.yieldSurf = yieldSurf
        if dilation is not None:
            self.dilation  = float(dilation)

    def fromHoekBrown(self, UCSi, GSI, mi, D=0.0, a=None, e=0.6, MR=None,
                      nu=0.2):
        """
        Update elastic properties (only if a value is given for MR) and the
        yield surface according to Hoek-Brown formulas.

        Usage:
         >>> plasticPeak = PlasticLevkowitchReuschSample(
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
            self.yieldSurf = YieldSurfMenetreyWillam()
        self.yieldSurf.fromHoekBrown(UCSi, GSI, mi, D, a, e)
        return self


class PlasticLevkowitchReusch(list):
    """This describes the elasto-plastic behaviour of a particular material:
    The interpolation of several yield surfaces, elastic properties and
    dilation.

    It's basically a list of PlasticLevkowitchReuschSample-objects in ascending
    order of pst values (which is not yet enforced technically).
    """

    def __init__(self, *args, **kwargs):
        """
        Accepts an iterable of (pst, elast, yield, dilation)-tuples.
        pst is the plastic strain value for this sample, elast is an
        ElasticLinear-object, yield is a YieldSurfMenetreyWillam-object, and
        dilation is the d-parameter of type float.

        Also accepts a PlasticLevkowitchReusch object, in which case it acts as
        a (deep) copy constructor.
        """
        if len(args)==1 and isinstance(args[0], self.__class__):
            list.__init__(self, imap(deepcopy, args[0]))
        elif len(args)==1 and hasattr(args[0], "__iter__"):
            list.__init__(self)
            for x in args[0]:
                if isinstance(x, PlasticLevkowitchReuschSample):
                    self.append(deepcopy(x))
                elif hasattr(x, "__iter__"):
                    self.append(PlasticLevkowitchReuschSample(*x))

    def fromComputator(self,UCSi,GSI,version='V13',ucsiIsInMPa=False,
                       a=0.5, e=0.6):
        """Creates a PlasticLevkowitchReusch sampleSequence using the
        Computator definitions.

        Usage:
            >>> pLR = PlasticLevkowitchReusch().fromComputator(65,50,
            ...       ucsiIsInMPa=True,version='V13')

        @param UCSi: uniaxial compressive strength of intact rock (in MPa
        or Pa, see UCSisInMPa)

        @param GSI:  Geological strength index (in %)

        @param version: ComputatorVersion

        @param ucsiIsInMPa: if true UCSi is assumed to be in MPa

        @param a: Menetrey-Willam exponent. default = 0.5

        @param e: Menetrey-Willam excentricity. default = 0.6

        """
        if ucsiIsInMPa:
            UCSi = UCSi*1E6
        if version == 'V12':
            pLRs = computatorV12(UCSi,GSI,returns='fullDict',
                                 ucsIsInMPa=False)
        elif version == 'V13':
            pLRs = computatorV13(UCSi,GSI,returns='fullDict',
                                 ucsIsInMPa=False)
        else:
            raise ValueError(
                "Don't know computatorVersion %s (yet)." % version)

        #pLRs comes as OrderedDict
        for pLR in pLRs.values():
            elastic = ElasticLinear(E=pLR['E'],nu=pLR['nu'])
            plastic = YieldSurfMenetreyWillam(e=e,a=a,UCSi=UCSi,**pLR)
            sample = PlasticLevkowitchReuschSample(pst=pLR['pst'],
                                                   elastic=elastic,
                                                   yieldSurf=plastic,
                                                   dilation=pLR['d'])
            self.append(sample)
        return self

    def isFromComputator(self,version=None,a=.5,e=.6):
        """Checks if the pLR-samples match a/any computatorVersion

        Usage:
            >>> pLR = PlasticLevkowitchReusch().fromComputator(65,50,
            ...       ucsiIsInMPa=True,version='V12')
            >>> print pLR.isFromComputator()
            >>> pLR.pop()
            >>> print pLR.isFromComputator()

        For computator versions V12 and V13 the GSI will be determined using
        the peak_mb. So the returned deltaDict referes to a computator material
        using this GSI.

        @param version: specifies the computator version(s) to be tested. If
        None, all implemented computator versions will be tested.

        @param a: Menetrey-Willam exponent. default = 0.5

        @param e: Menetrey-Willam excentricity. default = 0.6

        @return: matched computatorVersion or None if none matched
            deltaDictionary which will only have empty subDicts if pLR matches
            a computator version. If more than one computator version were
            checked, the deltaDictionary refers to the last one checked
        """

        from math import exp

        ucs = list(set([s.yieldSurf.UCSi for s in self]))
        if len(ucs) > 1:
            msg('Warning. Strange Material: LR-Samples have varying UCSi vals')
        if len(ucs) == 0:
            msg('pLR has no Samples or missing UCSi definition.')
            return None,{}

        ucs = ucs[0]

        if version is None:
            versions = ['V12','V13']
        elif isinstance(version,list):
            versions = version
        else:
            versions = [version]

        for ver in versions:
            if ver == 'V12' or ver == 'V13':
                # for convenicence use mbPeak (see Ver13 and Ver12)
                # use linar formula from computatorV13/computatorV12
                # mbPeak=GSI/100.*exp(0.01*UCSi)
                sPeak = self[0].yieldSurf.mb
                GSItest = 100.*sPeak*exp(-0.01*ucs/1E6)
            else:
                raise NotImplementedError(
                    'The computator version %s is not implemented yet!' % ver)

            #create a generic LR-sequence using first UCSi and GSItest
            computLR = PlasticLevkowitchReusch()
            computLR.fromComputator(ucs,GSItest,version=ver,a=a,e=e)
            deltas = self.compare(computLR)
            found = not any([d for d in deltas.values()])
            if found:
                # found matching computator version
                break

        if found:
            foundVersion = ver
        else:
            msg('The GSI-Value for testing was derived using S_PEAK')
            foundVersion = None
            if version is None:
                msg("Can't match any computatorVersions %s to LR-Samples."
                    " Return deltaDict for version %s" % (versions, ver))

        return foundVersion,deltas

    def compare(self,otherSamples,valTol=0.001,comments=True):
        """Compares the current LR-Sample with an other one and lists the
        differences in a dictionary. All parameter pairs (pst, elastic,
        yieldsurf and dilation) are evaluated using a relative tolerance.

           - In a first step, the plastic strain levels are compared
           - In second step all matched plastic strain levels (step one) are
           compared. This step is skipped if no pst-pairs were found.
           - In a third step all samples in the same position/index are
           compared. This step is skipped if all pst-values can be paired and
           the result would be the same as in step 2.

        The returned dictionary holds the differences (as an OrderedDict) for
        each off the comparison steps. If Both pLR-samples are identical within
        the relative tolerance no deltaItems are stored in the returned dicts.

        A deltaItem has a descriptive key and a tuple as value where the
        first tupleItem holds the entity of self and the second the (differing)
        one of the compared pLR-SampleSequence.

        For the purpose of better readability comments (deltaItems with empty
        tuple as value) can be stored in dicts.

        Example:
            >>> pLRa = PlasticLevkowitchReusch()
            >>> pLRa.fromComputator(65,50,ucsiIsInMPa=True)
            >>> pLRb = PlasticLevkowitchReusch()
            >>> pLRb.fromComputator(100,50,ucsiIsInMPa=True)
            >>> # old=pLRb[1].pst*(1+2E-14) #to test close psts
            >>> # setattr(pLRb[2],'pst',old)
            >>> pLRb.pop()
            >>> for key,subs in pLRa.compare(pLRb,comments=True).iteritems():
            >>>     print '+++ %s +++'%key
            >>>     for sKey,sub in subs.iteritems():
            >>>         if not sub:
            >>>             print '> %s'%sKey
            >>>         else:
            >>>             print '%s: \tpLRa=%s \tpLRb=%s'%(sKey,str(sub[0]),
            ...                                                   str(sub[1]))
        """
        from bae.future_01 import OrderedDict
        rDict = OrderedDict()

        def isInTol(valRef,valCmp,tol=valTol):
            '''Check if two values are equal within tolerance'''
            return abs(valRef-valCmp) <= tol*valRef

        def getValsInTol(valList,testVal,tol=valTol,tryExact=True):
            '''Return matched-list of values out of valList which are close to
            testVal. If the matched-list has more then one items and  tryExact
            is set True the tolerance is reduce recusively to force an exact
            match.
            '''
            matched = [val for val in valList if isInTol(val,testVal,tol=tol)]
            if tryExact and len(matched)>1 and tol>1E-12:
                matched = getValsInTol(valList,testVal,
                                       tol=0.1*tol,tryExact=True)
            return matched

        def deltaSample(sRef,sCmp,tol=valTol):
            '''Compares dictionary values and returns, those which are not
            equal within tolerance
            '''
            delta = dict((name,(sRef[name],sCmp[name]))
                         for name in sRef
                         if not isInTol(sRef[name],sCmp[name],tol=tol))
            return delta

        rDict = OrderedDict([('pst samples',OrderedDict()),
                             ('comparison by pst',OrderedDict()),
                             ('comparison by position',OrderedDict()),
                            ])

        pstRefs = OrderedDict([(s.pst,vars(s)) for s in self])
        pstCmps = OrderedDict([(s.pst,vars(s)) for s in otherSamples])
        # check variing length of samples
        if not len(pstRefs) == len(pstCmps):
            rDict['pst samples']['number of LR-samples'] = (len(pstRefs),
                                                            len(pstCmps))

        # get pstValues of compared LR-samples which can be found in reference
        pstCmpInRef = []
        for pstCmp in pstCmps:
            # try to find identical match
            f = getValsInTol(pstRefs,pstCmp)
            if len(f)>1:
                # found more then one match --> skip comparision by pst
                comment = ("In reference some pst-values are very close "+
                           "together. Comparison by pst will fail")
                rDict['pst samples'][comment] = ()
                rDict['pst samples']["close pst-values in reference"] \
                    = (pstRefs.keys(),pstCmp)

                pstCmpInRef = []
                break
            elif len(f)==1:
                # found match
                pstCmpInRef.append(f[0])
            else:
                pass

        # get pstValues of compared LR-samples which can be found in reference
        pstRefInCmp = []
        for pstRef in pstRefs:
            #try to find identical match
            f = getValsInTol(pstCmps,pstRef)
            if len(f)>1:
                # found more then one match --> skip comparision by pst
                comment = ("In compared LR-samples some pst-values are very "+
                           "close together. Comparison by pst will fail")
                rDict['pst samples'][comment] = ()
                rDict['pst samples']["close pst-values in compared"] \
                    = (pstRefs.keys(),pstCmp)

                pstRefInCmp = []
                break
            elif len(f)==1:
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

        # found equal pst-values --> compare these pairs
        if pstCmpInRef and pstRefInCmp:
            for pstRef,pstCmp in zip(pstCmpInRef,pstRefInCmp):
                # delta yieldSurf
                deltaDict = deltaSample(vars(pstRefs[pstRef]['yieldSurf']),
                                        vars(pstCmps[pstCmp]['yieldSurf']))
                # delta other (dilation)
                deltaOther = deltaSample({'d':pstRefs[pstRef]['dilation'],},
                                         {'d':pstCmps[pstCmp]['dilation']})
                # delta elastic
                deltaElast = deltaSample(vars(pstRefs[pstRef]['elastic']),
                                         vars(pstCmps[pstCmp]['elastic']))
                deltaDict.update(deltaOther)
                deltaDict.update(deltaElast)

                if deltaDict:
                    #found differences
                    comment = 'different values for matched pst=%.3f' % pstRef
                    rDict['comparison by pst'][comment] = ()

                    dd = dict(('pst=%.3f %s' % (pstRef,key),val)
                              for key,val in deltaDict.iteritems())
                    rDict['comparison by pst'].update(dd)

            if len(pstRefs) == len(pstCmpInRef):
                # all LR-samples could be compared by pst value
                # --> no need to compare them by postition
                if rDict['comparison by pst']:
                    comment = ('All samples found by pst. ' +
                               'Skipped comparing by position.')
                    rDict['comparison by position'][comment] = ()
                return rDict
        else:
            comment = "Can't match any pst. Comparison was skipped."
            rDict['comparison by pst'][comment] = ()

        for ii,(pstRef,pstCmp) in enumerate(zip(pstRefs,pstCmps)):
            # delta yieldSurf
            deltaDict = deltaSample(vars(pstRefs[pstRef]['yieldSurf']),
                                    vars(pstCmps[pstCmp]['yieldSurf']))
            # delta other (dilation)
            deltaOther = deltaSample({'d':pstRefs[pstRef]['dilation'],
                                      'pst':pstRefs[pstRef]['pst']},
                                     {'d':pstCmps[pstCmp]['dilation'],
                                      'pst':pstCmps[pstCmp]['pst']},)
            # delta elastic
            deltaElast = deltaSample(vars(pstRefs[pstRef]['elastic']),
                                     vars(pstCmps[pstCmp]['elastic']))
            deltaDict.update(deltaOther)
            deltaDict.update(deltaElast)

            if deltaDict:
                #found differences
                comment = 'different values sample %d' % ii
                rDict['comparison by position'][comment] = ()
                dd = dict(('sample=%d %s' % (ii,key),val)
                          for key,val in deltaDict.iteritems())
                rDict['comparison by position'].update(dd)

        if not comments:
            for dd in rDict.values():
                for key,val in dd.iteritems():
                    #remove all items with empty values
                    if val == ():
                        dd.pop(key,None)

        return rDict


class AbqVumatData(object):
    """Holds all data that is contained in a user material for our usual
    Abaqus material.

     - 2 PlasticLevkowitchReusch objects, one for intact rock and one for fill
     - a void ratio
     - critical damage
     - damping - alpha
     - density
     - umatDelete and umatDepvarNum parameter
     - an optional name

    @ivar name: name as in the Abaqus input file
    @ivar vumatType: "6tuple" (default), ...
    @ivar density: value from the *DENSITY option
    @ivar dampingAlpha: argument of the *DAMPING, ALPHA=... option
    @ivar umatDelete: argument of the *DEPVAR,DELETE=... option
    @ivar umatDepvarNum: value from the data line of the *DEPVAR option
    @ivar voidStiffnessRatio: PROPS(7) ELASTIC DEGRADATION RATIO
    @ivar critDamage: PROPS(8) CRITICAL DAMAGE

    @ivar plasticLRRock: PlasticLevkowitchReusch object for rock
    @ivar plasticLRFill: PlasticLevkowitchReusch object for fill
    """

    def __init__(self, *args, **kwargs):
        """
        The positional argument might be a
        L{bae.abq_model_02.internal.Material<bae.abq_model_internal_02.Material>}-object

        @kwarg inputFile: filename or other suitable argument to
             abq_model_02.Model.read() ... not implemented yet
        @kwarg materialName: Name of the material to read from inputFile
              ... not implemented yet
        @kwarg vumatType: "6tuple" (default), ...

        """

        if args and isinstance(args[0], AbqMaterial):
            abqMat = args[0]

            # name
            try:
                self.name = abqMat.name
            except AttributeError:
                msg("WARNING: material name missing.")
                pass

            # name
            try:
                self.vumatType = kwargs["vumatType"].lower()
            except AttributeError:
                raise ValueError(
                    "WARNING: vumatType argument must be a string."
                    " Got: %s" % kwargs["vumatType"])
            except KeyError:
                self.vumatType = "6tuple"

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


            if self.vumatType == "6tuple":
                if len(umatData)<9:
                    raise ValueError(
                        "This set of umat parameters does not have"
                        " enough values (only %d). Will be ignored."
                        % len(umatData))

                # PROPS(7) ELASTIC DEGRADATION RATIO
                self.voidStiffnessRatio = umatData[6]
                # PROPS(8) CRITICAL DAMAGE
                self.critDamage = umatData[7]

                # rock properties
                # PROPS(9) NUMBER (N) OF MATERIAL 6-TUPLES (ROCK)
                inSample = 8
                nSample = int(round(umatData[inSample]))
                if len(umatData)<inSample+1+6*nSample:
                    raise ValueError(
                        "There are only %d umat parameters for %d sample"
                        " points for rock, which is not enough."
                        % (len(umatData), nSample))

                self.plasticLRRock = PlasticLevkowitchReusch(
                    (umatData[inSample+1+i*6],        # pst value
                     ElasticLinear(
                         G=umatData[inSample+2+i*6],
                         K=umatData[inSample+3+i*6],  ),
                     YieldSurfMenetreyWillam(
                         UCSi=umatData[0],
                         e=umatData[1],
                         a=1.0/umatData[2],
                         s=umatData[inSample+4+i*6],
                         mb=umatData[inSample+5+i*6],  ),
                     umatData[inSample+6+i*6] )        # d for dilation
                    for i in range(nSample)  )

                # fill properties
                # PROPS(11+6*N) NUMBER (N) OF MATERIAL 6-TUPLES (FILL)
                inSample = inSample+1+6*nSample
                nSample = int(round(umatData[inSample]))
                if len(umatData)<inSample+1+6*nSample:
                    raise ValueError(
                        "There are only %d umat parameters for %d sample"
                        " points for fill, which is not enough."
                        % (len(umatData), nSample))

                self.plasticLRFill = PlasticLevkowitchReusch(
                    (umatData[inSample+1+i*6],        # pst value
                     ElasticLinear(
                         G=umatData[inSample+2+i*6],
                         K=umatData[inSample+3+i*6],  ),
                     YieldSurfMenetreyWillam(
                         UCSi=umatData[3],
                         e=umatData[4],
                         a=1.0/umatData[5],
                         s=umatData[inSample+4+i*6],
                         mb=umatData[inSample+5+i*6],  ),
                     umatData[inSample+6+i*6] )        # d for dilation
                    for i in range(nSample)  )

                if len(umatData) != inSample+1+6*nSample:
                    msg("WARNING: There are %d umat parameter whereas only %d"
                        " have been expected. Ignoring the rest."
                        % (len(umatData), inSample+1+6*nSample))

            else:  # not self.vumatType == "6tuple":
                raise ValueError("vumatType %s not supported yet."
                                 % self.vumatType)

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

        if self.vumatType == "6tuple":
            umatData = [
                # rock properties
                self.plasticLRRock[0].yieldSurf.UCSi,
                self.plasticLRRock[0].yieldSurf.e,
                1.0/self.plasticLRRock[0].yieldSurf.a,
                # fill props
                self.plasticLRFill[0].yieldSurf.UCSi,
                self.plasticLRFill[0].yieldSurf.e,
                1.0/self.plasticLRFill[0].yieldSurf.a,
                # PROPS(7) ELASTIC DEGRADATION RATIO
                self.voidStiffnessRatio,
                # PROPS(8) CRITICAL DAMAGE
                self.critDamage,
                ]
            # number of 6-tuples for rock
            umatData.append(len(self.plasticLRRock))
            umatData.extend(
                val for mat in self.plasticLRRock
                for val in (  # extend as flat list, last loop iterates fastest
                    mat.pst, mat.elastic.G, mat.elastic.K,
                    mat.yieldSurf.s, mat.yieldSurf.mb, mat.dilation ))
            # number of 6-tuples for fill
            umatData.append(len(self.plasticLRFill))
            umatData.extend(
                val for mat in self.plasticLRFill
                for val in (  # extend as flat list, last loop iterates fastest
                    mat.pst, mat.elastic.G, mat.elastic.K,
                    mat.yieldSurf.s, mat.yieldSurf.mb, mat.dilation ))

            model.material.updateItem(self.name, {
                    "DampingAlpha": self.dampingAlpha,
                    "Density": self.density,
                    'UmatDelete': self.umatDelete,
                    'UmatDepvarNum': self.umatDepvarNum,
                    'Umat': umatData,
                    })

            # check: UCSi, e, a must be the same for all

        else:  # not self.vumatType == "6tuple":
            raise ValueError("vumatType %s not supported yet."
                             % self.vumatType)


def computatorV12(UCSi,GSI,returns='fullDict',ucsIsInMPa=False):
    '''Computator V12 -- as found in matprops.xlsm

    Usage:
    >>> computator(68.,50.,ucsIsInMPa=True)


    @param UCSi: uniaxial compressive strength of intact rock (in MPa or Pa,
       see ucsIsInMPa)

    @param GSI: Geological strength index (in %)

    @param returns: specifies the output:
        - 'fullDict' : orderedDict with the keys 'PEAK', 'TRANS' and 'RES'
        - other : dictionary with props as keys and [peak,trans,res] as values

    @param ucsIsInMPa: if true UCSi is assumed to be in MPa

    '''
    from math import exp

    if returns is 'dependencies':
        a = [0,1,2]
        return dict(UCSi=dict(E=a,nu=a,s=a,mb=a,d=a),
                    GSI=dict(pst=[1,2],mb=[0,1]))

    if not ucsIsInMPa:
        mpa = 1E6
        UCSi = UCSi/mpa
    else:
        mpa = 1.

    #computator formulas taken from matprops.xlsm
    comp = {}
    comp['pst'] = [  0.,
                    (-0.0004*GSI+0.0489)*0.59429056563727,
                    (-0.0004*GSI+0.0489)*1.60458452722063 ]
    comp['E']   = [ 0.18664*UCSi,
                    0.18664*UCSi,
                    0.18664*UCSi ]
    comp['nu']  = [ 0.0006*UCSi+0.2,
                    0.0006*UCSi+0.2,
                    0.0006*UCSi+0.2]
    comp['s']   = [ 0.0002*exp(0.027*UCSi),
                    0.0002*exp(0.027*UCSi),
                    0.00002*UCSi ]
    comp['mb']  = [ GSI/100.*exp(0.01*UCSi),
                    GSI/100.*exp(0.01*UCSi),
                    0.0116*UCSi ]
    comp['d']   = [ 0.0025*UCSi,
                    0.0025*UCSi,
                    0.0009*UCSi ]

    if returns is 'fullDict':
        from bae.future_01 import OrderedDict as oDict
        return oDict((re,dict((name,vals[idx])
                              for name,vals in comp.iteritems() ))
                     for idx,re in enumerate(['PEAK','TRANS','RES']))

    else:
        return comp


def computatorV13(UCSi,GSI,returns='fullDict',ucsIsInMPa=False):
    '''Computator V13 -- as found in matprops.xlsm

    Usage:
        >>> computator(68.,50.,UCSisInMPa=True)

    @param UCSi: uniaxial compressive strength of intact rock (in MPa or Pa,
       see UCSisInMPa)

    @param GSI:  Geological strength index (in %)

    @param returns: specifies the output:
        - 'fullDict' : orderedDict with the keys 'PEAK', 'TRANS' and 'RES'
        - 'dependencies' : dict for UCSi and GSI dependend properties
        - other : dictionary with props as keys and [peak,trans,res] as values

    @param ucsIsInMPa: if true UCSi is assumed to be in MPa
    '''
    from math import exp

    if returns is 'dependencies':
        a = [0,1,2]
        return dict(UCSi = dict(pst=[1,2],E=a,nu=a,s=a,mb=a,d=a),
                    GSI  = dict(s=[0,1],mb=[0,1]))

    if not ucsIsInMPa:
        mpa = 1E6
        UCSi = UCSi/mpa
    else:
        mpa = 1.
    comp = {}

    #computator formulas taken from matprops.xlsm
    comp['pst']  = [ 0.,
                     0.04*exp(-0.006*UCSi),
                     0.07*exp(-0.006*UCSi) ]
    comp['E']    = [ 0.209*UCSi,
                     0.209*UCSi,
                     0.209*UCSi ]
    comp['nu']   = [ 0.0006*UCSi+0.2,
                     0.0006*UCSi+0.2,
                     0.0006*UCSi+0.2 ]
    comp['s']    = [ GSI/350000.*exp(0.027*UCSi),
                     GSI/350000.*exp(0.027*UCSi),
                     4e-08*(UCSi**2.25) ]
    comp['mb']   = [ GSI/100.*exp(0.01*UCSi),
                     GSI/100.*exp(0.01*UCSi),
                     0.0116*UCSi ]
    comp['d']    = [ 0.0025*UCSi,
                     0.0025*UCSi,
                     0.0009*UCSi ]

    if returns is 'fullDict':
        from bae.future_01 import OrderedDict as oDict
        return oDict((re,dict((name,vals[idx])
                              for name,vals in comp.iteritems() ))
                     for idx,re in enumerate(['PEAK','TRANS','RES']))
    else:
        return comp
