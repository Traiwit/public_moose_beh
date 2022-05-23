# -*- coding: utf-8 -*-
"""
Heuristic models to derive material properties
"""

from collections import OrderedDict
from math import exp, log

def hoekBrown( UCSi, GSI, mi, D=0.0, Ei=None, a=None):
    """ HOEK-BROWN FAILURE CRITERION

    Usage:
     >>> hbData = hoekBrown(UCSi, GSI, mi)

    @param UCSi: uniaxial compressive strength of intact rock
    @param GSI: Geological strength index (in %)
    @param D: see "HOEK-BROWN FAILURE CRITERION - 2002 EDITION"
    @param mi: see "HOEK-BROWN FAILURE CRITERION - 2002 EDITION"
    @param a: if not given will be computed by the Hoek-Brown formula
      based on GSI, otherwise usually: a=0.5
    @param Ei: intact rock deformation modulus, if None an other HB-estimation
        is used
    @returns: parameterDictionary
    @Note: All units SI!
    @note: from "HOEK-BROWN FAILURE CRITERION - 2002 EDITION"
    """
    hBdict = {}

    if a is None:
        hBdict['a'] = 0.5 + ( exp(-GSI/15.0) - exp(-20.0/3.0) ) / 6.0
    else:
        hBdict['a'] = float(a)

    if Ei is None:
        hBdict['E'] = 1E5 * (1 - D/2.)/(1 + exp((75+25*D-GSI)/11.))
    else:
        hBdict['E'] = Ei * (0.02 + (1 - D/2.)/(1 + exp((60+15*D-GSI)/11.)))

    hBdict['UCSi'] = float(UCSi)
    hBdict['s']    = exp( (GSI - 100.0) / (9.0 - 3.0*D) )
    hBdict['mb']   = mi * exp( (GSI - 100.0) / (28.0 - 14.0*D) )
    return hBdict



def computatorV12(UCSi, GSI, returns='fullDict'):
    '''Computator V12 -- as found in matprops.xlsm

    Usage:
    >>> params = computatorV12(68E6, 50., )

    @param UCSi: uniaxial compressive strength of intact rock

    @param GSI: Geological strength index (in %)

    @param returns: specifies the output:
        - 'fullDict' : orderedDict with the keys 'PEAK', 'TRANS' and 'RES'
        - 'dependencies' : dict for UCSi and GSI dependend properties/samples
        depends on which input value
        - other : dictionary with props as keys and [peak,trans,res] as values

    @note: All units SI!!!
    '''
    if returns == 'dependencies':
        a = [0,1,2]
        return dict( UCSi = dict(E=a, nu=a, s=a, mb=a, dilation=a),
                     GSI  = dict(pst=[1,2], mb=[0,1]) )

    comp = {}

    #computator formulas taken from matprops.xlsm (there ucs in MPa, E in GPa)
    UCSi /= 1E6
    comp['UCSi']= 3*[UCSi*1E6]
    comp['pst'] = [  0.,
                    (-0.0004*GSI+0.0489)*0.59429056563727,
                    (-0.0004*GSI+0.0489)*1.60458452722063 ]
    comp['E']   = [ 0.18664*UCSi*1E9,
                    0.18664*UCSi*1E9,
                    0.18664*UCSi*1E9 ]
    comp['nu']  = [ 0.0006*UCSi+0.2,
                    0.0006*UCSi+0.2,
                    0.0006*UCSi+0.2]
    comp['s']   = [ 0.0002*exp(0.027*UCSi),
                    0.0002*exp(0.027*UCSi),
                    0.00002*UCSi ]
    comp['mb']  = [ GSI/100.*exp(0.01*UCSi),
                    GSI/100.*exp(0.01*UCSi),
                    0.0116*UCSi ]
    comp['dilation']   = [ 0.0025*UCSi,
                           0.0025*UCSi,
                           0.0009*UCSi ]

    if returns == 'fullDict':
        return OrderedDict((re, dict((name, vals[idx])
                              for name,vals in comp.iteritems() ))
                     for idx,re in enumerate(['PEAK', 'TRANS', 'RES']))

    else:
        return comp



def computatorV13(UCSi, GSI, returns='fullDict'):
    '''Computator V13 -- as found in matprops.xlsm

    Usage:
        >>> params = computatorV13(68., 50.)

    @param UCSi: uniaxial compressive strength of intact rock

    @param GSI:  Geological strength index (in %)

    @param returns: specifies the output:
        - 'fullDict' : orderedDict with the keys 'PEAK', 'TRANS' and 'RES'
        - 'dependencies' : dict for UCSi and GSI dependend properties/samples
        - other : dictionary with props as keys and [peak,trans,res] as values

    @note: All units SI!!!
    '''
    if returns == 'dependencies':
        a = [0,1,2]
        return dict( UCSi = dict(pst=[1,2], E=a, nu=a, s=a, mb=a, dilation=a),
                     GSI  = dict(s=[0,1], mb=[0,1]) )

    comp = {}

    #computator formulas taken from matprops.xlsm (there ucs in MPa, E in GPa)
    UCSi /= 1E6

    comp['UCSi'] = 3*[UCSi*1E6]
    comp['pst']  = [ 0.,
                     0.04*exp(-0.006*UCSi),
                     0.07*exp(-0.006*UCSi) ]
    comp['E']    = [ 0.209*UCSi*1E9,
                     0.209*UCSi*1E9,
                     0.209*UCSi*1E9 ]
    comp['nu']   = [ 0.0006*UCSi+0.2,
                     0.0006*UCSi+0.2,
                     0.0006*UCSi+0.2 ]
    comp['s']    = [ GSI/350000.*exp(0.027*UCSi),
                     GSI/350000.*exp(0.027*UCSi),
                     4e-08*(UCSi**2.25) ]
    comp['mb']   = [ GSI/100.*exp(0.01*UCSi),
                     GSI/100.*exp(0.01*UCSi),
                     0.0116*UCSi ]
    comp['dilation']  = [ 0.0025*UCSi,
                          0.0025*UCSi,
                          0.0009*UCSi ]
    # check for the residual>trans problem
    comp['mb'][2] = min(comp['mb'][1],comp['mb'][2])
    comp['s'][2] = min(comp['s'][1],comp['s'][2])


    if returns == 'fullDict':
        return OrderedDict((re, dict((name, vals[idx])
                              for name, vals in comp.iteritems() ))
                     for idx,re in enumerate(['PEAK', 'TRANS', 'RES']))
    else:
        return comp


def effectiveUCSiV143(UCSi, UCSLim=60, kinkParam=0.05):
    """Intruduced with V14.3 to limit the pst for trans and residual level.

    Example:
     >>> from matplotlib import pyplot as plt
     >>> ucss = [10*ii for ii in range(1,21)]
     >>> ucsEffs = [_effectiveUCSiV143(u) for u in ucss]
     >>> plt.plot(ucss, ucsEffs)
     >>> plt.show()

    """
    effUCSi = UCSi + log(1+exp(-kinkParam*(UCSi-UCSLim)))/float(kinkParam)
    return effUCSi

def computatorV14(UCSi, GSI, mi=None, mr=400., returns='fullDict',
                  subVersion=3):
    '''Computator V14 -- as found in matPropsAniso_VladsComputor_2020Apr M04

    Usage:
        >>> params = computatorV14(68., 50., subVersion=3)

    @param UCSi: uniaxial compressive strength of intact rock

    @param GSI: Geological strength index (in %)

    @param mi: see "HOEK-BROWN FAILURE CRITERION - 2002 EDITION", if None
        mi will be set to 5 + 0.06* UCSi

    @param returns: specifies the output:
        - 'fullDict' : orderedDict with the keys 'PEAK', 'TRANS' and 'RES'
        - other : dictionary with props as keys and [peak,trans,res] as values

    @param mr: scale YOUNG's modulus on PEAK level. If None, mr will be set to
        400 which equals the default input

    @param subversion: default is 3
        0: initial V14 version. Youngs modulus is linear in mr and UCSi
        1: Geros adjustment: Youngs modulus is nonlinear in mr, UCSi and GSI
            Warning: This implementation had a bug, that yields in no
            elastic degradation.
        2: fixed elastic degradation bug from V14.1
        3: limited trans and residual psts by usage of an effective UCSi
            see L{effectiveUCSiV143}
        4: changes nu so that BulkModuls remains constant during softening

    @note: All units SI!!!
    '''

    # constants
    GSI_res = 27.0
    dp1 = 0.07
    if subVersion >= 3:
        dp2=0.08
    else:
        dp2 = 0.10
    elDegRatio = 0.25

    comp = {}

    if mr is None:
        mr = 400.

    # computator formulas assume UCSi in MPa
    UCSi /= 1E6
    comp['UCSi'] = 3*[UCSi*1E6,]

    # approximation for mi from UCSi, based on BE's model calibrations
    if mi is None:
        mi = 5 + 0.06* UCSi

    # effective GSI:
    # peak: from BE's experience our models need a lower GSI
    # trans: heuristical formular mildly depending on UCSi:
    #       low UCSi -> GSI ~= GSI_res...30
    #       high UCSi -> GSI -> 50
    # res: just broken, it's gravel by now
    gsi = [.75*GSI+5,
           GSI_res+(50.0-GSI_res)*( (UCSi**1.5) / (UCSi**1.5 + 200))**2.2,
           GSI_res,]

    # UCS degradation parameter for rock mass
    # modified Hoek Brown formula: exp((GSI-100)/9) --> exp((GSI-100)/7.5)
    comp['s'] = [exp((gsi[0] - 100)/7.5),
                 exp((gsi[1] - 100)/7.5),
                 1E-5,]  # effectively zero cohesion at residual level

    # UCS rock mass
    ucsRm = [UCSi * comp['s'][0]**.5,
             UCSi * comp['s'][1]**.5,
             UCSi * comp['s'][2]**.5,]

    # damage level for trans and res level...
    # A pst increase of dp1 converts intact rock (GSI 100) to trans state.
    # The actual rock starts with a lower GSI and therefore needs
    # proportionally less additional pst to reach trans state.
    # The aforesaid holds for UCS=100. pst_trans is altered
    # antiproportionally with UCS by means of the trailing (100/UCS) factor..
    # A further pst increase of dp2 converts rock of UCS 100 from trans to
    # res state. Again this threshold is altered antiproportionally by UCS.
    if subVersion >= 3:
        UCSEff = effectiveUCSiV143(UCSi)
    else:
        UCSEff = UCSi

    comp['pst'] = [0.,
                   dp1* ((gsi[0] - gsi[1])/(100 - gsi[1])) * (100./UCSEff),]
    comp['pst'].append(dp2* 100./UCSEff + comp['pst'][-1])

    # confinement dependency parameter for rock mass
    # modified Hoek Brown formula: exp((GSI-100)/28) --> exp((GSI-100)/25)
    comp['mb'] = [mi * exp((gsi[0] - 100)/25.),
                  mi * exp((gsi[1] - 100)/25.),
                  mi * exp((gsi[2] - 100)/25.),]

    comp['dilation'] = [comp['mb'][0]/6.,
                        comp['mb'][1]/6.,
                        0.]

    # Hoek-Brown-formula for alpha
    comp['a'] = 3*[0.5 + (exp(-gsi[0]/15.) - exp(-20/3.))/6.,]

    # elastic properties
    # modified Hoek-Brown formula: exp((60-GSI)/11 --> exp((60-GSI)/50
    # Youngs modulus degradation is attached to UCS degradation.
    # elDegRatio=1 means elastic degradation is proportional to UCS
    # degradation. elDegRatio=0 means no elastic degradation.
    if subVersion == 0 or subVersion == '0':
        Epeak = UCSi * mr * (0.02 + 1/(1 + exp((60-gsi[0])/50.))) * 1E6
        ERatioT = ((1-elDegRatio)*ucsRm[0] + elDegRatio*ucsRm[1]) / ucsRm[0]
        ERatioR = ((1-elDegRatio)*ucsRm[0] + elDegRatio*ucsRm[2]) / ucsRm[0]
        comp['E'] = [Epeak,
                     ERatioT * Epeak,
                     ERatioR * Epeak,]
    elif subVersion in [1,2,3,4,'1','2','3', '4']:
        comp['E'] = []
        for ii in range(3):
            # geros excelEstimation
            # E =50*(1-(1-MR*80/50/1000)^(UCS/80))*(0.46+(1.01-0.46)/(1+EXP((60-GSI)/11)))
            if subVersion in [1,'1']:
                cGSI = gsi[0] # that was a bug in V14.1
            else:
                cGSI = gsi[ii]
            fA = 1-mr*80./50000.
            if fA < 0:
                # Note that this heuristic scheme only alows mr up to 625
                raise ValueError("mr must be smaller than 625. Yours is %.3f." % mr)
            fB = 0.46+(1.01-0.46)/(1+exp((60-cGSI)/11.))
            E = 50*(1-fA**(UCSi/80.))*fB
            E *= 1E9
            comp['E'].append(E)
    else:
        raise ValueError('Do not know subversion %s for computator V14')

    if subVersion in [4, '4']:
        refNu = 0.25
        refK = comp['E'][0]/(3.0*(1.0 - 2*refNu))
        comp['nu'] = [0.5*(1 - E/(3*refK)) for E in comp['E']]
    else:
        comp['nu'] = 3*[0.25,]

    # now handles invalid gsi-'orders'
    if gsi[0] <= gsi[2]:  # GSI_peak lower than GSI_res
        samples = [0,]  # take only peak value
    elif gsi[0] <= gsi[1]:  # GSI_peak is between trans and residual
        comp['pst'][2] = dp2*((gsi[0] - gsi[2])/(gsi[1] - gsi[2]))*(100./UCSi)
        samples = [0,2]  # drop transitional
    else:  # peak, trans res as usual
        samples = [0,1,2]

    # dilation of last sample has to be 0
    comp['dilation'][-1] = 0


    for key, val in comp.iteritems():
        comp[key] = [val[ii] for ii in samples]

    if returns == 'fullDict':
        if len(samples) == 1:
            names = ['PEAK',]
        elif len(samples) == 2:
            names = ['PEAK', 'RES']
        else:
            names = ['PEAK', 'TRANS', 'RES']

        return OrderedDict((re, dict((name, vals[idx])
                              for name, vals in comp.iteritems() ))
                           for idx,re in enumerate(names))
    else:
        return comp


def restoreComputatorInputs(version, **kwargs):
    """Recalculates input values that are needed to create computator materials
    but not get (neccessarily) stored in a L{material}-object.
    The values are recalculated from input values of the first (peak) sample!
    Required input-kwargs --> returned values:
        - 'V12'/'V13': mb, UCSi --> GSI
        - 'V14': s, E, UCSi (mb*) --> GSI, mr (mi*); *if mb is passed mi will
            be recalculated and returned if it dosent match BE's default
            approximation (see computatorV14)

    Example:
        >>> data = computatorV14(80E6, 60, mr=600, mi=2)
        >>> dd = data['PEAK']
        >>> print restoreComputatorInputs('v14.1', **dd)

    @param version: specifies computator version ('V12', 'V13', 'V14.0','V14.1')

    @note: all units in SI

    """
    from math import log

    rawVersion = version
    version = version.lower().strip().replace('.','')

    if version in ['v12', 'v13']:
        try:
            mbPeak = kwargs['mb']
            ucs    = kwargs['UCSi']
        except KeyError:
            raise KeyError('For computator V12/V13 UCSi and mb(peak) are'
                           ' required (kwarg-)input arguments')

        GSI = 100.*mbPeak*exp(-0.01*ucs/1E6)
        return {'GSI':GSI}

    elif 'v14' in version:
        try:
            sPeak = kwargs['s']
            ucs   = kwargs['UCSi']
            EPeak = kwargs['E']
        except KeyError:
            raise KeyError('For computator V14 UCSi, s(peak) and E(peak)'
                           'are required (kwarg-)input arguments')
        if ucs <= 0:
            return {'GSI':None, 'mr':None}
        GSI0 = log(sPeak)*7.5 + 100
        GSI = (GSI0-5)/0.75
        if version == 'v140':
            mr = (EPeak/ucs)/ (0.02 + 1/(1 + exp((60-GSI0)/50.)))
        elif version in ['v141', 'v142', 'v143','v144']:
            fB = 0.46+(1.01-0.46)/(1+exp((60-GSI0)/11.))
            fA0 = (1-(EPeak*1E-9) / (50.*fB))
            if fA0 < 0:
                # was not a computator V14.1-Peak sample
                return {'GSI':None, 'mr':None}

            fA = fA0**(80E6/float(ucs))
            mr = (1-fA)*(50000/80.)
        else:
            raise ValueError("Unknown (sub)Version of computator V14: %s"
                             % str(rawVersion))

        result = {'GSI':GSI, 'mr': mr}

        if 'mb' in kwargs:
            mbPeak = kwargs['mb']
            mi = mbPeak / exp((GSI0 - 100)/25.)
            if not mi == (5 + 0.06*(ucs/1E6)):
                # mi-input was not None
                result['mi'] = mi

        return result
    else:
        raise ValueError("Don't know a computatorVersion %s" % str(version))

if __name__ == '__main__':
    print('No syntax errors')
    # dd=computatorV14(80E6,55, subVersion=4)['PEAK']
    # print restoreComputatorInputs('v14.4', **dd)