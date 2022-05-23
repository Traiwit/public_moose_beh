# -*- coding: utf-8 -*-
"""
Functions and Classes to describe YieldSurfaces

Examples:
========
 Test if a (Stress)point yields or not
    >>> from bae.material_01.yieldsurfaces import ExtendedMenetreyWillam
    >>> from bae.material_01.yieldsurfaces import stressToPQTheta
    >>> ySurf = ExtendedMenetreyWillam(s=0.005, mb=.8, UCSi = 50E6,)
    >>> p,q,theta = stressToPQTheta(myStressTensors)
    >>> yieldsOrNot = ySurf.yieldSurfaceFun(p,q,theta)

 Plot a yieldsurface on compression meridian
    >>> S1,_,S3 = ys.sampleMeridianEigs(S3Max=10E6, samples=32)
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> ax.plot(S3/1E6, S1/1E6)
    >>> ax.set_xlabel(r'$S_{minor}$ / MPa')
    >>> ax.set_ylabel(r'$S_{major}$ / MPa')
    >>> fig.show()

 Fit a MOHR-COULOMB material onto EMW in given range
    >>> ysEMW = ExtendedMenetreyWillam(UCSi=40E6, s= 3.365E-4, mb=0.597)
    >>> mcFit = QuasiMohrCoulomb.cohsFrictionfromExtendedMenetreyWillamFit
    >>> c,phi,fit = mcFit(ysEMW, fittingRange=(None,10E6),
    ...                   constrainedCohesion=50E3,)
    >>> ysMC = QuasiMohrCoulomb(c=c,phi=phi)

"""
import sys

from math import pi
import warnings

try:
    import numpy as np
    # raise ImportError # for testing
    np.seterr(all='raise')
except ImportError:
    warnings.warn('\n\n' + 80*'*' + '\n'
                  "This module requires NUMPY which seems to be not " +
                  "installed yet. Many functions may be not avilable.\n" +
                  "You may see AttributeErrors on NoneType-objects, if " +
                  "numpy-dependend functions are called.\n" + 80*'*' +'\n\n')
    np = None

from bae.log_01 import msg



# eval if epydoc is running
# --> suppress decorators, to get reasonable docStrings
if 'epydoc' in sys.modules:
    epyDocBuild = True
else:
    epyDocBuild = False



from models import hoekBrown, computatorV12, computatorV13, computatorV14

### FUCTIONS AND  DECORATORS TO ORGANIZE IN AND OUTPUT ########################
#{Helper functions and decorators
def _streesDictToTensor(d):
    """Maps a dictionary holding stress components or lists of components
    into a numpy array
    """
    tMapFull = lambda s: np.array([[s[0],s[3],s[4]],
                                   [s[3],s[1],s[5]],
                                   [s[4],s[5],s[2]]])

    # check if the dictionanry holds a sequence of stressTensors
    # i) all vaules come as iterables (tuples, lists, ...) --> map
    if all([hasattr(val,'__iter__') for val in d.values()]):
        mapping = True
        if not all([len(d.values()[0]) == len(val) for val in d.values()]):
            raise ValueError('Creation of stress sequence failed. Lengths of '
                             'component sequences differ!')

    # ii) single value for each stressItem --> no mapping required
    elif not any([hasattr(val,'__iter__') for val in d.values()]):
        mapping = False

    # unknown structure
    else:
        raise ValueError('All or not any stress components must ' +
                         'come as iterables')

    # map conversion to all items
    if mapping:
        # i) full stressTensor given
        try:
            return np.array(map(tMapFull,zip(d['S_11'], d['S_22'], d['S_33'],
                                             d['S_12'], d['S_13'], d['S_23'])))
        except KeyError:
            pass

        # ii) eigenValues of stressTensor given
        try:
            sortDiag = lambda s: np.diag(sorted(s)[::-1])
            return np.array(map(sortDiag, zip(d['S1'], d['S2'], d['S3'])))

        # unknown inputType
        except KeyError:
            raise KeyError('StressDict must contain keys ' +
                           'S_11, S_22, S_33, S_12, S_13, S_23 or '+
                           'S1, S2, S3')

    # apply conversation to single item given
    else:
        try:
            return tMapFull([d['S_11'], d['S_22'], d['S_33'],
                             d['S_12'], d['S_13'], d['S_23']])
        except KeyError:
            pass

        try:
            return np.diag(sorted([d['S1'], d['S2'], d['S3']])[::-1])

        except KeyError:
            raise KeyError('StressDict must contain keys ' +
                           'S_11,S_22,S_33,S_12,S_13,S_23 or '+
                           'S1,S2,S3')

#} #End Helper functions and decorators


def stressToPQTheta(*args, **kwargs):
    """Returns pressure, MISES-stress and LODE-angle S{theta} of StressTensor
    S. Here the LODE-Angle is expressed via M{cos(3S{theta})}. Thus S{theta}
    ranges between M{0 S{<=} S{theta} S{<=} S{pi}/3} (if M{S1S{>=}S2S{>=}S3}),
    where the

     - triaxial compression (TXC) mode is at S{pi}/3
     - shear (SHR) mode is at S{pi}/6 and
     - triaxial extension (TXE) is at 0.

    @param S: sequences of N
        - tensors (N X 3 X 3)
        - eigenStresses in one array (N x 3)
        - eigenStresses in 3 array 3*(N)

    @param fullLodeAngle: return 'full' LODE angle range 0-2S{pi} for unsorted
        EigenValues (default = False).

    @return: p, q, theta (in rad)

    """
    args = map(np.asarray, args)
    if len(args) == 1:
        S = args[0]
        if len(S.shape) == 3 and S.shape[1:] == (3,3):
            # full stressTensor as input
            # --> get sorted s1, s2, s3 via eigs
            S = args[0]
            # valid is needed to exclude nan from eigenValue solver
            valid = (~(np.isnan(S))).all(axis=2).all(axis=1)
            E = np.sort(np.linalg.eigvals(S[valid]), axis=1)[:,::-1]
            s1, s2, s3 = 3*[np.full(S.shape[0], np.nan)]
            s1[valid], s2[valid], s3[valid] = E.T
        elif len(S.shape) == 2 and S.shape[1] == 3:
            # eigenVals as one matrix
            s1, s2, s3 = S.T
        else:
            raise ValueError('Only series of nx3x3 tensors or nx3 array'
                             ' [s1,s2,s3] allowed here.')
    elif len(args) == 3:
        # eigenValues as three distinct arrays
        s1, s2, s3 = map(np.asarray, args)
    else:
        raise ValueError('Only series of nx3x3 tensors or nx3 array'
                             ' [s1,s2,s3] allowed here.')

    fullLodeAngle = kwargs.pop('fullLodeAngle', False)

    p = -(s1 + s2 + s3) / 3.

    try:
        q = np.sqrt(.5* ( (s1 - s2)**2
                        + (s2 - s3)**2
                        + (s3 - s1)**2))
    except:
        print('S1,S2,S4', s1, s2, s3)
        raise

    ii = ~np.isclose(q, 0) # calc theta only where mises exists

    # option a: arcCos-argument via invariants
    J2 = (q**2) / 3.
    J3 = (s1 + p) * (s2 + p) * (s3 + p)
    cosThreeTheta = ((27/4.)**0.5) * J3[ii] / J2[ii]**1.5

    # option b: explicitly
    # r3 = ((-2.*s1+s2+s3)**3.+(s1-2.*s2+s3)**3.+(s1+s2-2.*s3)**3.)/6.
    # cosThreeTheta = r3[ii] / q[ii]**3.

    cosThreeTheta = np.clip(cosThreeTheta, -1, 1) # clip to [-1,1]
                                                  #because of rounding errors
    t = np.zeros_like(q)

    t[ii] = np.arccos(cosThreeTheta)/3.

    if not fullLodeAngle:
        theta = t
    else:
        # to overcome sorting issues we should implement a theta correction
        # according to the following lookup
        ## thetaRages vs. majorStressRanges
        #    t==0        >> S1==S2
        #    0<t<60      >> S1>S2>S3
        #    t==60       >> S2==S3
        #    60<t<120    >> S1>S3>S2
        #    t==120      >> S1==S3
        #    120<t<180   >> S3>S1>S2
        #    t=180       >> S1==S2
        #    180<t<240   >> S3>S2>S1
        #    t=240       >> S2==S3
        #    240<t<300   >> S2>S3>S1
        #    t=300       >> S1==S3
        #    300<t<360   >> S2>S1>S3
        pp = 2*np.pi/3
        tt = np.array([0, pp, pp, 2*pp, 2*pp, 3*pp])
        ii = np.array([(0,1,2), #0<t<60      >> S1>S2>S3
                       (1,0,2), #60<t<120    >> S2>S1>S3
                       (1,2,0), #120<t<180   >> S2>S3>S1
                       (2,1,0), #180<t<240   >> S3>S2>S1
                       (2,0,1), #240<t<300   >> S3>S1>S2
                       (0,2,1), #300<t<360   >> S1>S3>S2
                      ])#
        #### And now, some index magic:
        iiOrder = np.array([0, 5, 1, 2, 4, 3]) # indexOrder when ii gets sorted
                                               # used to 'trace' the index
                                               # after unique was used
        # get sorting of eigenvalues (highest first)
        sortI = np.argsort([s1,s2,s3], axis=0)[::-1].T
        uni, inv = np.unique(np.concatenate((ii, sortI)),
                           axis=0,return_inverse=True)
        jj = iiOrder[inv[6:]] # get each index jj of sortI in ii
        ttOff = tt[jj]
        theta = t*(-1)**(jj) + ttOff # toggle symmetric thetas for uneven jj

    return p, q, theta


def pQThetaToEigs(p, q, theta, miningNotation=False):
    """Transforms pressure, misesStress and LODE angle into Eigenstress space.
    Here, the LODE-angle is definied via cos(3*S{theta}) (see L{stressToPQTheta}).
    This transformation is valid for a full 2S{pi} range (unsorted eigen-
    stresses) and for 0-2S{pi}/3 range.

    @param p: hydrostatic pressure (M{-1/3 * Tr(B{S})}) in Pa
    @param q: vonMISES Stress in Pa
    @param theta: LODE angle defined via M{cos(3*S{theta})} in rad
    @param miningNotation: if True, order and sign get swapped
    @returns: (S1, S2, S3) in Pa. Will be sorted (highest first) if theta
        ranges between 0 and 2S{pi}/3.
    """
    if hasattr(q,'__iter__'):
        q = np.asarray(q)[np.newaxis]
        if q.size == 1:
            q = q.flatten()[0]

    if hasattr(p, '__iter__'):
        p = np.asarray(p)[np.newaxis]
        if p.size == 1:
            p = p.flatten()[0]


    twothrdPi = 2./3.*np.pi
    ct = np.array([np.cos(theta),
                   np.cos(theta - twothrdPi),
                   np.cos(theta + twothrdPi),])

    twothrdQ = 2./3.*q
    try:
        s1,s2,s3 = -p + twothrdQ * ct
    except ValueError as e:
        print('probably a dimension mismatch?:')
        print('p.shape:',p.shape)
        print('ct.shape', ct.shape)
        print('q.shape', twothrdQ.shape())
        raise e

    if miningNotation:
        s1,s2,s3 = -s3, -s2, -s1

    return s1,s2,s3


def metToVoigt(M):
    """Shapes sym. matrix (upper triangle) to voigt-vector.
    Component order is according to ABAQUS notation:
        V = [v11, v22, v33, v12, v23, v31] .
    """
    M = np.asarray(M)
    return M[[0,1,2,0,1,0], [0,1,2,1,2,2]]


def voigtToMet(V):
    """Shapes voigt-vector of shape (6xN) to sym. matrix.
    Component order is according to ABAQUS notation:
        V = [v11, v22, v33, v12, v23, v31] .
    """
    return np.array([ [V[0],V[3],V[5]],
                      [V[3],V[1],V[4]],
                      [V[5],V[4],V[2]] ])


def rotVoigtTransversal(V, alpha, beta):
    """Transversal-isotropic (back)rotation of a
    voigt-vector/tensor as used for anisotropy.

    @param V: tensorData in voigtNotation of shape (6xN) where
        order of compoments is according to ABAQUS notation:
        V = [v11, v22, v33, v12, v23, v31] .
    @param alpha: first angle describing the normal vector given with:
        M{S{alpha} = atan2(ny, nx)}
    @param beta: second angle describing the normal vector given with:
        M{S{beta} = acos(nz)}

    @note: code is taken from
        vusubs lrm_c.BasicRoutinesJoinedMaterial.BRotTeCo6
    """
    sA, cA = np.sin(alpha), np.cos(alpha)
    sAq, cAq = sA**2, cA**2
    sB, cB = np.sin(beta), np.cos(beta)
    sBq, cBq = sB**2, cB**2
    vR = 6*[None,]

    #sR[0] = Power(Cb,2)*Power(Sa,2)*s[1] + Power(Sb,2)*s[2] + Cb*(2*Ca*Cb*Sa*s[3]
    #        - 2*Sa*Sb*s[4] - 2*Ca*Sb*s[5] + Power(Ca,2)*Cb*s[0]);
    vR[0] = (cBq*sAq*V[1] + sBq*V[2] + 2*cA*cBq*sA*V[3]
             - 2*cB*sA*sB*V[4] - 2*cA*cB*sB*V[5] + cAq*cBq*V[0])

    #sR[1] = Power(Ca,2)*s[1] + Sa*(-2*Ca*s[3] + Sa*s[0]);
    vR[1] = cAq*V[1] + -2*cA*sA*V[3] + sAq*V[0]

    #sR[2] = Power(Sa,2)*Power(Sb,2)*s[1] + Power(Cb,2)*s[2] + Sb*(2*Ca*Sa*Sb*s[3]
    #        + 2*Cb*Sa*s[4] + 2*Ca*Cb*s[5] + Power(Ca,2)*Sb*s[0]);
    vR[2] = (sAq*sBq*V[1] + cBq*V[2] + 2*cA*sA*sBq*V[3]
             + 2*cB*sA*sB*V[4] + 2*cA*cB*sB*V[5] + cAq*sBq*V[0])

    #sR[3] = Ca*Cb*Sa*s[1] + Cb*(Power(Ca,2) - Power(Sa,2))*s[3]
    #       - Ca*Sb*s[4] + Sa*Sb*s[5] - Ca*Cb*Sa*s[0];
    vR[3] = (cA*cB*sA*V[1] + cB*(cAq - sAq)*V[3]
             - cA*sB*V[4] + sA*sB*V[5] - cA*cB*sA*V[0])

    #sR[4] = Ca*Sa*Sb*s[1] + (Power(Ca,2) - Power(Sa,2))*Sb*s[3]
    #        + Ca*Cb*s[4] - Cb*Sa*s[5] - Ca*Sa*Sb*s[0];
    vR[4] = (cA*sA*sB*V[1] + sB*(cAq - sAq)*V[3]
             + cA*cB*V[4] - cB*sA*V[5] - cA*sA*sB*V[0])

    #sR[5] = Cb*Power(Sa,2)*Sb*s[1] - Cb*Sb*s[2] + 2*Ca*Cb*Sa*Sb*s[3]
    #        + Power(Cb,2)*Sa*s[4] - Sa*Power(Sb,2)*s[4] + Ca*Power(Cb,2)*s[5]
    #        - Ca*Power(Sb,2)*s[5] + Power(Ca,2)*Cb*Sb*s[0];
    vR[5] = (cB*sAq*sB*V[1] - cB*sB*V[2] + 2*cA*cB*sA*sB*V[3]
             + cBq*sA*V[4] - sA*sBq*V[4] + cA*(cBq - sBq)*V[5]
             + cAq*cB*sB*V[0])

    return vR

#} #End TensorMath
###############################################################################


#{ Some rock mechanics functions
def nTauFromBalmersMethod(Smajor, Sminor):
    """Calculates normal and shearStress as displayed in RockLab.
    Found in HOEK, I{A brief history of the development of the Hoek-Brown failure
    criterion}, 2006 (see pdf in RockLab Help).

    From my (Tobias) point of view we should ban this strange 2d approach and
    and use p-q-plot for S{theta}=30deg instead.

    @param Smajor: major principle stress (in mining notation)
    @param Sminor: minor principle stress (in mining notation)
    @note: This method may fail for tesile range of S3, where the gradient can
        be negative. Also, at least three distinct Sminor are required.
    """
    from scipy.interpolate import UnivariateSpline
    spline = UnivariateSpline(Sminor,Smajor)
    d1d3 = np.array([spline.derivatives(n)[1] for n in Sminor])

    n = Sminor + (Smajor-Sminor)/(1+d1d3)
    t = (n - Sminor) * np.sqrt(d1d3)
    return n, t


def sMajorSMinorFromNTau(sigma, tau):
    """Calculates major and minor principle stress from normal and shear
    stress. This is the reversed function of L{nTauFromBalmersMethod}.
    """
    from scipy.interpolate import UnivariateSpline
    spline = UnivariateSpline(sigma, tau)
    derivative = np.array([spline.derivatives(n)[1] for n in sigma])

    d = tau / np.cos(np.arctan(derivative))
    sigmaM = sigma + derivative*tau
    sMajor = sigmaM + d
    sMinor = sigmaM - d
    return sMajor, sMinor


#Barton-Bandis Criterion
def bartonBandis(n, phir, jrc, jcs):
    """
    The Barton-Bandis failure criterion is an empirical relationship widely
    used to model the shear strength of rock discontinuities (e.g. joints).
    It is very useful for fitting a strength model to field or laboratory
    shear test data of discontinuities (see here U{https://www.rocscience.com/help/rocdata/rocdata/Barton-Bandis_Criterion.htm}).

    @param n: normal stress (see L{nTauFromBalmersMethod}) (same unit as jcs)

    @param jrc: joint roughness coefficient

    @param jcs: joint wall compressive strength (same unit as n)

    @param phir: in degree

    @returns: yielding shear stress
    """
    return n*np.tan(np.pi/180.*(phir + jrc*np.log10(jcs/n)))

#} # End RockMechanics

# Default-Setting for Menetry-Willam-Parameters and their bounds
defaultParamMenetreyWillam = dict(
                                a = .5,
                                e = .6,
                                s = 0.005,
                                mb= .8,
                                UCSi = 50E6,
                                anisoS = 1.,
                                anisoN = 1.,
                                )
try:
    inf = np.inf
except AttributeError:
    inf = 1E12

paramBoundsMenetreyWillam = dict(
                            UCSi = (1E-6, 1E12),
                            e = (0.5, 1.),
                            a = (1E-6, 1.),
                            s = (1E-6, 1.),
                            mb = (1E-12, 1000.),
                            anisoS = (1., inf),
                            anisoN = (1., inf),
                            )


#{YieldSurface-Classes
class ExtendedMenetreyWillam(object):
    """Describes an isotropic extendedMenetreyWillam yield surface.
    It provides some usefull mathematical operations for this kind of yield
    surface.
    """
    _type = 'MenetreyWillam'
    _name = 'MenetreyWillam'

    # complete parameter-keys for current yield-surface._type
    _yieldParameters = ['UCSi', 'e', 'a', 's', 'mb']

    # (free) changeable parameters
    _freeYieldParameters = _yieldParameters

    # posible input data formats for fit choose from here:
    #   ['full tensor','full eigs','Sminor and Smajor','Smajor']
    _fitInputRequirements = ['full tensor', 'full eigs']

    def __init__(self, *args, **kwargs):
        """Creates an ExtendedMenetreyWillam yield surface object.
        The parameters can be passed via:

         1. keyWordArguments. Required parameters that are not specified
            in kwargs will be set to None (except e and a).

         2. a yieldSurface object as first positional argument and additional
            keyWordArguments. If a parameter is given in both, the value from
            kwargs will be set. If a parameter is missing either in ysObject
            and kwargs it will be set to None (except e and a).

        The second method can act like a copy constructor or a converter method
        for closely related yieldSurface objects.

        Example:
         >>> from bae.material_01 import ExtendedMenetreyWillam
         >>> from bae.material_01 import AnisotropicExtendedMenetreyWillam
         >>> myParams = dict(UCSi=60E6, mb=.7, s=.02)
         >>> ysIso = ExtendedMenetreyWillam(**myParams)
         >>> ysIso2 = ExtendedMenetreyWillam(ysIso, mb=.8)
         >>> ysAniso = AnisotropicExtendedMenetreyWillam(ysIso,
         ...                                             anisoS=1, anisoN=1.1)

        @kwarg UCSi: UCS for intact rock
        @kwarg e: excentricity, default 0.6
        @kwarg a: Hoek-Brown-exponent, default 0.5
        @kwarg s: degradation parameter of UCS
        @kwarg mb: frictional strength parameter
        """
        params = dict()
        if len(args) == 1:
            for key in self._freeYieldParameters:
                try:
                    params[key] = getattr(args[0], key)
                except AttributeError:
                    pass
        elif len(args) > 1:
            raise ValueError("Only one positional argument (should be a yield"
                         " surface object) allowed, but %d given." % len(args))

        params.update(kwargs)

        try:
            params['e']
        except KeyError:
            params['e'] = .6

        try:
            params['a']
        except KeyError:
            params['a'] = .5

        for key in self._yieldParameters:
            setattr(self, key, params.get(key, None))

        # sets bounds for parameterFit
        self._postInit()

    def toDict(self):
        return dict((attr, getattr(self, attr)) for attr in self._yieldParameters)

    @classmethod
    def fromDict(cls, input):
        return cls(**input)

    def _postInit(self):
        """Include this in the last line of __init__ if a new Yield object is
        created.
        Updates:
            - initital boundSettings
        """
        self.bounds = dict((k, val)
                            for k,val in paramBoundsMenetreyWillam.iteritems()
                            if k in self._freeYieldParameters)

## creation methods
    def fromHoekBrown(self, UCSi, GSI, mi, **kwargs):
        """
        Creates a YieldSurface from Hoek-Browns formula found in "HOEK-BROWN
        FAILURE CRITERION - 2002 EDITION"

        Usage:
            >>> mw = ExtendedMenetreyWillam().fromHoekBrown(50E6, 55, 0.9)

        @param UCSi: uniaxial compressive strength of intact rock (SI!)
        @param GSI: geological strength index (in %)
        @param mi:
        @returns: self
        """
        self.__init__( **hoekBrown( UCSi, GSI, mi, **kwargs ) )

        return self


    def fromComputator(self, UCSi, GSI,
                       sample='PEAK', version='V14.3', **kwargs):
        """
        Creates a YieldSurface from a computator pst-sample. Actualy it just
        provides mb and s. See
        L{bae.material_01.material.PlasticLevkovitchReusch} to create a
        'complete' material from computator.

        Usage:
            >>> mw = ExtendedMenetreyWillam().fromComputator(50E6, 55,
            ...                                 sample=1, a=.52)

        @param UCSi: uniaxial compressive strength of intact rock (SI!)
        @param GSI: geological strength index (in %)
        @param sample: specifies PST-sample, can be a string ('PEAK', 'TRANS',
            'RES') or index
        @param version: computatorVersion (implemented: 12, 13, 14)
        @returns: self
        @note: kwargs for a and e will be passed to self also. But if kwargs
            contains mb or s the computator values will be overwritten by those
        """
        rawVersion = version
        version = version.strip().lower().replace('.','')
        version
        if '13' in version:
            cParam = computatorV13(UCSi, GSI)
        elif '12' in version:
            cParam = computatorV12(UCSi, GSI)
        elif 'v14' in version.lower():
            if version == 'v143':
                subVersion = 3 # default
            elif version == 'v142':
                subVersion = 2
            elif version == 'v141':
                subVersion = 1
            elif version == 'v140':
                subVersion = 0
            else:
                raise ValueError('Do not know subversion %s for computator V14'
                                 % version)
            cArgs = {'mi':kwargs.pop('mi', None),
                     'subVersion':subVersion,}
            cParam = computatorV14(UCSi, GSI, **cArgs)
        else:
            raise ValueError('Do not know a Computator %s' % version)

        if type(sample) == str:
            try:
                param = cParam[sample.upper()]
            except KeyError:
                raise KeyError('Computator %s does not have a pst-sample '+
                               'called %s' % (version, sample.upper()) )
        elif type(sample) == int:
            try:
                param = cParam[cParam.keys()[sample]]
            except IndexError:
                raise IndexError('Computator %s has only %d pst-samples' %
                                 (version, len(cParam)) )

        #append kwargs like e or a to paramDict
        param.update(kwargs)

        self.__init__(**param)

        return self


    def copy(self):
        params = dict((param, getattr(self,param))
                      for param in self._freeYieldParameters)
        return self.__class__(**params)


## description of yieldsurface
    def yieldParams(self):
        """Returns dictionary with parameters of current yieldSurface"""
        return dict((k, getattr(self, k)) for k in self._yieldParameters)


    def freeYieldParams(self):
        """Returns dictionary with free parameters of current yieldSurface"""
        return dict((k, getattr(self, k)) for k in self._freeYieldParameters)


    @staticmethod
    def _hoekBrownYieldSurf(sMinor, UCSi, s, mb, a):
        """HOEK-BROWN yield criterion."""
        # to handle sqrt of negative values (--> stress state is 'beyond'
        # isotropic tensile stregth) we'll use a complex expression here
        # and will only resume the real part
        UCSi += 0. # force UCSi to be a float
        sMajor = sMinor + UCSi*(mb*sMinor/UCSi + s + 0j)**a
        sMajor[np.abs(sMajor.imag) > 0] = np.nan

        return sMajor.real

    def hoekBrownYieldSurf(self, sMinor):
        """Equivalent HOEK-BROWN yield criterion for comparision.
        Stress states beyond isotropic tensile strength will have nan as value.

        @note: Mind that sMinor and sMajor are in mining convention (tension<0)
        @note: see L{_hoekBrownYieldSurf} to have full parameter control
        """
        return self._hoekBrownYieldSurf(sMinor, self.UCSi, self.s, self.mb,
                                        self.a)


    @staticmethod
    def eccentricity(theta, e):
        """convex elliptic function (Klisinski 1985). See _testEccentricity for
        usage and plot
        """
        from math import pi
        #symmetric for each modulo 2*pi/3
        theta = np.mod(theta, 2*pi/3.)
        if hasattr(theta, '__iter__'):
            if not isinstance(theta, np.ndarray):
                theta = np.array(theta)
            theta[theta>pi/3.] = 2*pi/3. - theta[theta>pi/3.]
        elif theta >= pi/3.:
            theta = theta - 2*pi/3.

        cosT  = np.cos(theta) #use np.cos to allow vectorization
        cosSq = cosT**2
        brA = 1.-e**2
        brB = 2.*e-1.
        convEllip = ( 4.*brA*cosSq + brB**2)/ \
                    ( 2*brA*cosT + brB*np.sqrt(4*brA*cosSq +  5*e**2-4*e) )
        return convEllip


    @classmethod
    def _yieldSurfaceFun(cls, p, q, theta, UCSi=None, mb=None, s=None,
                        e=.6, a=.5) :
        """MENETRY-WILLAM-YieldSurface, expressed in pressure, MISESstress
        and LODE-angle. Parameters need to be passed as keyword arguments
        """
        UCSi += 0. # force UCSi to be a float
        a += 0.    # force a to be a float
        extc = cls.eccentricity(theta, e)
        q[np.isclose(q, 0, atol=1E-6)] = 0.
        q[q<0] = 0
        f = (q / UCSi)**(1/a) +  mb * (extc*q  / (3*UCSi) - p / UCSi) - s
        return f


    def yieldSurfaceFun(self, p, q, theta):
        """MENETRY-WILLAM-YieldSurface, expressed in pressure, MISES-stress
        and LODE-angle.

        @note: see L{_yieldSurfaceFun} to have full parameter control
        """
        return self._yieldSurfaceFun(p, q, theta, **self.yieldParams() )


    def _dYielddP(self, p):
        """Deviation of YieldSurfaceFuction in p"""
        return -np.full_like(p, self.mb) / float(self.UCSi)

    def _dYielddQ(self, q, theta):
        """Deviation of YieldSurfaceFunction in MISES-stress q"""
        ecc = self.eccentricity(theta, self.e)
        UCSi, a = float(self.UCSi), float(self.a)
        return (1./(UCSi*a)) * (q/UCSi)**(1/a - 1) + self.mb * ecc / (3*UCSi)


    def _dYielddTheta(self, theta):
        """deviation of YieldSurfaceFunction in theta - Not Implemented Yet!"""
        raise NotImplementedError('Do you really need this? Try yourself!')


## derive properties of yieldsurface
    def getIsoTensileStrength(self, eps=1E-10):
        """Returns the Tesile Strength on isotropic load """
        return -self.UCSi * self.s / self.mb + eps


    def getUCS(self):
        """Returns the UniaxialCompressiveStrength:
            UCS = S1_yield where S2=0 and S3=0
        """
        #return self.ptOnYieldSurfEigen(S2=[0], S3=[0]).max()
        return self.UCSi*(self.s)**self.a


    def getUniaxialShearStrength(self):
        """Returns the UniaxialShearStrength
            USS = S1_yield where S2=S1 and S3=0
        @note: This method is untested yet.
        """
        try:
            from scipy.optimize import newton
        except ImportError:
            raise ImportError('This function requires SciPy')

        def func(S):
            lc = stressToPQTheta(0,S,S)
            return self.yieldSurfaceFun(*lc)

        val = newton(func, 0)

        p,q,t = stressToPQTheta(0, val, val)

        return val


## getting points on yieldsurface
    def ptOnYieldSurfEigen(self, **kwargs):
        """ Returns the sequence of missing eigenStressComponents to a given
        set of eigenStresses. As the YieldSurface usually has two valid values
        for each given EigenStressPair, the returned values will have the shape
        of Srequested.shape = (len(Sgiven),2), where the largest is set first.

        Usage:
            >>> yS = ExtendedMenetreyWillam(UCSi=55E6,mb=...)
            >>> S3 = yS.ptOnYieldSurfEigen(S1=[s1_a,s1_b,...],
            ...                             S2=[s2_a,s2_b,...])

        see also _testMenetreyWillam
        """
        try:
            from scipy.optimize import fsolve, newton
        except ImportError:
            raise ImportError('This function requires SciPy')

        valKeys = ['S1', 'S2', 'S3']
        foundVals = ','.join([valKey for valKey in valKeys
                              if valKey in kwargs])

        if foundVals == 'S1,S2':
            SA, SB =  np.atleast_1d(kwargs['S1']), np.atleast_1d(kwargs['S2'])
            idxs = [0,1,2]
        elif foundVals == 'S1,S3':
            SA, SB =  np.atleast_1d(kwargs['S1']), np.atleast_1d(kwargs['S3'])
            idxs = [0,2,1]
        elif foundVals == 'S2,S3':
            SA, SB =  np.atleast_1d(kwargs['S2']), np.atleast_1d(kwargs['S3'])
            idxs = [2,0,1]
        # flatten inputArrays
        inShape = SA.shape
        SA = SA.ravel()
        SB = SB.ravel()

        # set or estimate startValue
        if 'startVal' not in kwargs:
            # these values work fine while evaluating the compression meridian
            s0 = np.c_[SA, SB].min(axis=1) + self.getIsoTensileStrength()

            s0 = np.full_like(SA, 10000*np.c_[SA, SB].max())
            s1 = np.full_like(SA, -10000*np.c_[SA, SB].max())
        else:
            s0 = kwargs['startVal']
            s1=s0

        SS = np.zeros((3, len(SA)))
        SS[0] = SA
        SS[1] = SB
        def func(SC):
            SS[2] = SC
            lc = stressToPQTheta(SS[idxs].T, fullLodeAngle=True)
            return self.yieldSurfaceFun(*lc)


        SCa = newton(func, s0, tol=1E-5, maxiter=50)
        SCb = newton(func, s1, tol=1E-5, maxiter=50)
        SCa,SCb = np.sort(np.c_[SCa, SCb],axis=1).T
        # SCb = newton(func, -SCa, tol=1E-5, maxiter=50)

        #print np.vstack(stressToPQTheta(SS[idxs].T, fullLodeAngle=True)).T
        return SCa.reshape(inShape),SCb.reshape(inShape)


    def fastQApproximator(self, pRange, includeTheta=False, nSamplePts = 64,
                            distribSamplePts = 'lin'):
        """Returns a fast estimatorFunction f for mises-stress values of known
        p and theta values: q = f(p, theta), where p and theta need to be of
        same shape or one of both is a single value which is to be assumed to
        be constant.
        Its based on the analytical solution of yieldSurface for given q and
        theta. These values are sampled (nSamplePts, distribSamplePts) and
        interpolated. If the theta-dependency is excluded (default), an
        excentricity of e=1 is assumed.

        Example:
            >>> eMW = ExtendedMenetreyWillam().fromComputator(50E6,50)
            >>> qEstimatior = eMW.fastQApproximator([-10E6, 100E6])
            >>> p, theta = linspace(-10E6, 100E6, 100000), 3.141/5.
            >>> qEst = qEstimator(p, theta)

        @param pRange: range [min(xi),max(xi)] or complete p-array

        @param nSamplePts: number of test-samples to setup estimator spline,
            default=16, at least 3

        @param distribSamplePts: distribution of test-samples 'log'(default)
            or 'lin'

        @returns: an estimatorFunction(p, theta) for q which is valid within
            qRange. For p < isotropicTensileStrength it will return nan

        @note: sloving q at the yieldsurface via L{ptOnYieldSurf} is very fast.
            So your code only benifits from this method for a very large number
            of requested points.
        """
        from scipy.interpolate import interp1d, griddata
        from scipy.optimize import newton
        if nSamplePts < 3:
            raise ValueError('SplineEstimation needs at least 3 samplePoints')

        pRange = np.asarray(pRange)
        tStrength = self.getIsoTensileStrength()
        idxs = pRange >= tStrength
        pRange = pRange[idxs]

        #create samplingPoints for spline
        pLower = max(tStrength, pRange.min())
        pUpper = max(pRange.max(), pLower + abs(tStrength)/100.)
        comprAngle = np.full(2, np.deg2rad(60))
        def minFun(q):
            q[q<0] = 0
            pEst = self.ptOnYieldSurf(q=q, theta=comprAngle)
            res = np.abs(pEst - [pLower, pUpper])
            return res
        lower, upper = newton(minFun, np.abs([pLower, pUpper]),
                             tol=.1, rtol=1E-3)

        upper *= 1.01 # overestimate to avoid
        lower = max(lower, abs(tStrength)*1E-6) #has to be > 0
        if distribSamplePts == 'lin' :
            qIp = np.linspace(lower, upper, nSamplePts)
        elif distribSamplePts == 'log':
            # use finer resolution for smaler rho (log-distributed spacing)
            qIp = np.geomspace(lower, upper, nSamplePts)
        else:
            raise ValueError("Don't know distribSamplePts %s" %
                                str(distribSamplePts))
        if includeTheta:
            # logspaced theta distribution get instable
            # thetas = np.arange(-16+1, 16)
            # thetas = np.sign(thetas)*(2**np.abs(thetas))
            # thetas = (np.pi/3.)*(thetas/float(thetas.max()) + 1)
            thetas = np.linspace(0, 2*np.pi/3., 32+1)
            q, t = np.meshgrid(qIp, thetas)
            shape = q.shape
            pIp = self.ptOnYieldSurf(q=q.ravel(),
                                     theta=t.ravel())
            pIp = pIp.reshape(shape)

            def qEstimator(p, theta):
                p, theta = np.asarray(p), np.asarray(theta)

                theta = np.mod(theta, 2*np.pi/3.)

                if p.shape == theta.shape:
                    shape = p.shape
                elif theta.shape == (1,):
                    shape = p.shape
                    theta = np.full_like(p, theta)
                elif  p.shape == (1,):
                    shape = theta.shape
                    p = np.full_like(theta, p)
                else:
                    print('p.shape: ', p.shape)
                    print('theta.shape: ', theta.shape)
                    raise ValueError('p and theta must be (a) of same ' +
                                     'shape or (b) one of both has to be '+
                                     'a single value')
                p = p.ravel()
                theta = theta.ravel()
                idx = (p >= tStrength)
                if p.max() > pRange.max():
                    msg('Warning!!! Range of xi-values extends' +
                        ' inital sample range. Overanging samples' +
                        ' will be set to nan' )
                result = np.nan*np.empty(shape).ravel()
                # result[idx] = rho2dRbf(xi[idx], theta[idx]) #Rbf Fails
                result[idx]= griddata(
                                    np.array([pIp.ravel(), t.ravel()]).T,
                                    q.ravel(),
                                    (p[idx], theta[idx]),
                                    method='cubic', rescale=True,
                                    fill_value=0.)
                try:
                    ax.scatter(p, theta, result, c='k', s=10)
                except NameError:
                    pass

                return result.reshape(shape)

        else:
            pIp = self.ptOnYieldSurf(q=qIp,
                                      theta=np.zeros(qIp.shape))

            qSpline = interp1d(pIp, qIp, kind='cubic',
                                 fill_value="extrapolate")

            def qEstimator(p, theta):
                p, theta = np.asarray(p), np.asarray(theta)
                if not theta.shape:
                    theta = theta * np.ones(q.shape)
                elif not p.shape:
                    p = p * np.ones(theta.shape)
                idx = (p >= self.getIsoTensileStrength())
                if p.max() > pRange.max():
                    msg('Warning!!! Range of q-values extends' +
                        ' inital sample range. Overanging samples' +
                        ' will be set to nan' )
                result = np.zeros(p.shape) * np.nan
                result[idx] = qSpline(p[idx])
                return result

        return qEstimator


    def ptOnYieldSurf(self, **kwargs):
        """ Returns the sequence of missing a value to a given pair of p, q,
        theta (pressure, MISES, LODE angle). If theta is requested
        it will be returned modulo(theta,pi/3). All angles are then given by
        symmetry: allThetas = [theta*n*pi/3 for n in range(6)]

        Usage:
            >>> yS = MenetreyWillam(UCSi=55E6,mb=...)
            >>> q = yS.ptOnYieldSurf(p=myPs, theta=myThetas)

        see also _testMenetreyWillam

        @kwarg p: list or np.array of pressure
        @kwarg q: list or np.array of MISES-stress
        @kwarg theta: list or np.array of LODE-angle
        @kwarg startVal: list or single value as starting point of zeroSolving
        @return: np.array of unset p, q or theta
        """
        try:
            from scipy.optimize import fsolve, newton
        except ImportError:
            raise ImportError('This function requires SciPy')

        valKeys = ['p', 'q', 'theta']
        foundVals = ','.join([valKey for valKey in valKeys
                              if valKey in kwargs])

        p, q, theta = None, None, None

        if foundVals == 'p,q':
            p  = np.atleast_1d(kwargs['p'])
            q = np.atleast_1d(kwargs['q'])

        elif foundVals == 'q,theta':
            q = np.atleast_1d(kwargs['q'])
            theta = np.atleast_1d(kwargs['theta'])

        elif foundVals == 'p,theta':
            p  = np.atleast_1d(kwargs['p'])
            theta = np.atleast_1d(kwargs['theta'])

        else:
            raise ValueError('You have to have to pass (p,q), (p,theta) ' +
                             'or (p,theta) as keyWordArguments. '+
                             'You passed: %s'%','.join([n for n in kwargs]))

        if foundVals == 'p,q':
            startVal = kwargs.get('startVal', 0)
            if not hasattr(startVal, '__iter__'):
                startVal = startVal*np.ones(p.shape[0])

            #limit rhos to HB *[eccent0.5,eccent1] here!!!

            fZero = lambda theta: self.yieldSurfaceFun(p, q, theta)
            val = fsolve(fZero, startVal)
            val = np.mod(val, np.pi/3.)

        elif foundVals == 'q,theta':
            # can be solved directly
            ii = q >= 0
            val = np.full_like(q, np.nan)
            ecc = self.eccentricity(theta[ii], self.e)
            qq = q[ii]
            UCSi = self.UCSi + 0. # force float
            a = 1/(self.a + 0.)   # force float

            val[ii] = UCSi*((qq/UCSi)**a - self.s)/self.mb + qq*ecc/3.

        elif foundVals == 'p,theta':
            #solve only for valid points (xi>=isotropic tensile strength)
            idx = (p >= self.getIsoTensileStrength())

            startVal = kwargs.get('startVal', 2*np.abs(p))

            #probably passing HoekBrownValues as initial guess fails here?
            if not hasattr(startVal,'__iter__'):
                startVal = startVal*np.ones(p.shape[0])

            #create empty/nan array
            val = np.nan*np.ones(p.shape)

            #wrapper for yieldSurfFunction --> fZero = f(q)
            def fZero(q):
                negZero = np.logical_and(q < 0,  np.isclose(q, 0))
                q[negZero] = 0
                return self.yieldSurfaceFun(p[idx], q, theta[idx])

            def fPrime(q):
                negZero = np.logical_and(q < 0,  np.isclose(q, 0))
                q[negZero] = 0
                return self._dYielddQ(q, theta[idx])

            result = newton(fZero, startVal[idx],
                               tol=1E-5,
                               fprime = fPrime,
                               maxiter=50)
            result[result<0] = 0
            val[idx] = result

        return val


    def fitToStressData(self, S, parameters=None, stressWeighted=True,
                       dataType=None, miningNotation=False):
        """
        Fits yieldSurfaceParameters to given yielding stress values using a
        least square fit. The sequence of tressTensor values S, where 'material
        starts to yield' can be arbitrary and needs to be transformed to the
        three Haigh-Westergaard-Coordinates-vectores by the function
        self._transformFitInputData. By default the following StressDataTypes
        flagged by the keywordArg dataType are implemented here:

            - None: (numpy)array of N stressTensors with the shape Nx3x3
            - "full tensor": dictionary with keys S_11,S_22,S_33,S_12,S_13,S_23
            each holding a sorted sequence of stressValues
            - "full eigen": with keys S1,S2,S3 holding sequence of eigen-Values
            - "Smaj,Smin": dictionary with keys S1 and S3, assuming S2=S3

        @param S: sequence of stressTensor values

        @kwarg dataType: string indicating the type of stressfield (see
        description), defaut = None --> S is Array of shape (Nx3x3)

        @kwarg parameter: specifies explicitly the parameters to get fitted,
        default = None --> fit all yieldParameters

        @note: When passing dataType="Smaj,Smin", be aware of correct order
        S1>S3, miningNotation will be forced True
        """
        try:
            from scipy.optimize import least_squares
        except ImportError:
            raise ImportError('This function reqires sciPy')

        if dataType == "Smaj,Smin":
            S1, S2, S3 = S[0], S[1], S[1]
            miningNotation = True
        elif isinstance(S, tuple) and len(S) == 3:
            S1,S2,S3 = S
        else:
            raise NotImplementedError('only input for the three major' +
                                      ' stresses implemented yet')

        if miningNotation:
            S1,S2,S3 = -S3,-S2,-S1

        p,q,theta = stressToPQTheta(S1,S2,S3)
        #determine parameters to fit and store in list to preserve order
        if parameters is None: #use all
            pKeys = self.freeYieldParams().keys()
            pKeys.remove('UCSi')

        else: #use specified
            if not all([para in self._freeYieldParameters
                        for para in parameters]):
                raise KeyError('Some specified parameters are unknown or ' +
                               'not changeable for Class %s'%self._name)
            pKeys = parameters

        #dict holding all parameters (will be updated )
        mwParam = self.yieldParams()

        #initial squared error sum as scaling factor to unifiy convergence crit
        # use _yieldSurfaceFun here to pass arbitrary parameters
        scale = 1#np.nansum((self._yieldSurfaceFun(p,q,theta,**mwParam)**2))
        def minFun(x):
            """ErrorFunction to get minimized"""
            cDict = dict(zip(pKeys,x))  #zip(pKeys,x) will preserve paramOrder
            mwParam.update(cDict)

            # error value for each stress sample
            # use _yieldSurfaceFun here to pass arbitrary parameters
            errors = self._yieldSurfaceFun(p,q,theta,**mwParam)
            errors /= scale

            errors[np.isnan(errors)] = 10E6
            return errors

        #StartValues
        p0 = [self.freeYieldParams()[k] for k in pKeys]
        #Bounds
        bounds = tuple(np.array([list(self.bounds[k]) for k in pKeys]).T)

        #run least squares
        p = least_squares(minFun, p0,
                          # ftol=1E-12,
                          # xtol=1E-12,
                          # loss='soft_l1',
                          bounds=bounds
                          )

        return dict(zip(pKeys, p.x)), p.cost


    def sampleMeridianEigs(self, theta=pi/3., samples=128,
                           distribution='lin', miningNotation=True,
                           **kwargs):
        """Returns S1, S2 and S3 of yieldsurface for a given LODE-angle theta.
        This is very usefull for plotting.

        Example: plot compression meridian (theta=pi/3)
            >>> import matplotlib.pyplot as plt
            >>> S1,_,S3 = ys.sampleMeridianEigs(theta=np.pi/3,
            ...                                 miningNotation=True,
            ...                                 S3Min=15E6)
            >>> plt.plot(S3, S1)
        @param theta: specifies the LODE angle of meridian (e.g. compression
            meridian at pi/3)
        @param samples: number of samples taken along hydrostatic axis (p)
        @param distribution: defines how sampling points get distributet along
            p-axis. 'lin' sets an equidistant, 'log' a logarithmic and 'alpha'
            a root function like spacing. The later gives a kind of natural
            distribution for MeneteryWillam/HoekBrown yieldsurfaces because it
            has a higher density in low confinment regime where the curvature
            is high if alpha is not 1
        @param miningNotation: returns eigenstresses in mining notation
        @kwarg SiMin: with i in (1,2,3) lower bound of eigenstresses to output.
            If more than one component is given, the one with the lowest
            related p-value is choosen. If none is specified, p will be set
            to isotropic tensile strength.
        @kwarg SiMax: with i in (1,2,3) upper bound of eigenstresses to output.
            If more than one component is given, the one with the highest
            related p-value is choosen. If none is specified, p will be set
            to -15-times of the isotropic tensile strength.
        @returns: S1,S2,S3 each of size samples and according to parameter
            miningNotation
        """
        # preparing sampling of p
        if 'lin' in distribution.lower():
            dist = np.linspace(0,1,samples)
        elif 'alpha' in distribution.lower():
            dist = np.arange(samples)**(1/self.a)
        elif 'log' in distribution.lower():
            dist = np.logspace(0, 1, samples)
        else:
            raise ValueError('')
        dist -= dist.min()
        dist /= dist.max()

        # getting the min/max values for p from SiMin and SiMax
        defaultMin = self.getIsoTensileStrength() # lower default p bound
        defaultMax = -15*defaultMin # uper default p bound

        sMins, sMaxs = [], []
        for ii in range(3): # try to find Smins and Smaxs in kwargs
            try:
                kk = dict(miningNotation=miningNotation) #to expand as kwargs
                kk['S%d' % (ii+1)] = [float(kwargs['S%dMin' % (ii+1)]),]
                sMins.append(kk)
            except (KeyError, TypeError):
                pass # SiiMin was None or not in kwargs
            try:
                kk = dict(miningNotation=miningNotation) #to expand as kwargs
                kk['S%d' % (ii+1)] = [float(kwargs['S%dMax' % (ii+1)]),]
                sMaxs.append(kk)
            except (KeyError, TypeError):
                pass # SiiMax was None or not in kwargs

        # decide which Smin, Smax to use
        pLRange, pURange = [], []
        if sMins:
            pLRange.extend([getPAndQFromThetaAndS(self, theta,**jj)[0][0]
                           for jj in sMins])
        else:
            pLRange.append(defaultMin) # no SiMin found --> set default

        if sMaxs:
            pURange.extend([getPAndQFromThetaAndS(self, theta,**jj)[0][0]
                               for jj in sMaxs])
        else:
                pURange.append(defaultMax) # no SiMax found --> set default

        pL = sorted(pLRange)[0]
        pU = sorted(pURange)[-1]
        if pL == pU:
            raise ValueError('Lower and Upper bounds of p are identical.'
                             ' Check input for SiMin and SiMax (i in {1,2,3})')

        p = dist*(pU-pL) + pL # p-samples

        theta = np.full_like(p, theta) # expand theta to same shape as p
        q = self.ptOnYieldSurf(p=p, theta=theta) # get vonMises stress

        # back to eigenStress
        S1,S2,S3 = pQThetaToEigs(p, q, theta, miningNotation=miningNotation)

        return S1,S2,S3


    def _transformFitInputData(self, S, dataType):
        """Transformes a sequence of stressTensorValues to an array of
        pressure, vonMises stress and Lode angle.
        Specify this for special purpose.
        """
        if dataType is None:
            #S is array of shape Nx3X3
            return stressToPQTheta(S)

        elif dataType == 'full tensor' or dataType == 'full eigs':
            return stressToPQTheta(_streesDictToTensor(S))

        elif dataType == 'Smaj,Smin':
            # assumes stresses to be in miningNotation at compression meridian
            Sn = np.empty((S.shape[0],3))
            Sn[:,0] = S[:,0]
            Sn[:,1] = S[:,1]
            Sn[:,2] = S[:,1]
            return stressToPQTheta(_streesDictToTensor(Sn),
                                   miningNotation=True)

        else:
            raise ValueError('FitData InputType not supported for Class %s.'%\
                             self._name + ' Specifiy _transformFitInputData()'+
                             ' if you need somthing very special or pass ' +
                             ' the StressTensor as (native) Nx3x3-Array.')



class QuasiMohrCoulomb(ExtendedMenetreyWillam):
    """
    Mimics MohrCoulomb yield surface via Menetry-Willam. See also
    U{BE-wordpress<http://192.168.31.25/wordpress/archives/2028>}.
    """
    _name = 'QuasiMohrCoulomb'

    #free parameterKeys
    _freeYieldParameters = ['UCSi', 'mb']

    # posible input data formats for fit choose from here:
    #   ['full tensor','full eigs','Sminor and Smajor','Smajor']
    _fitInputRequirements = ['full tensor', 'full eigs', 'Smaj,Smin',]

    def __init__(self, *args, **kwargs):
        """
        Creates a QuasiMohrCoulomb yield surface object.
        The parameters can be passed via:

         1. keyWordArguments. Required parameters that are not specified
            in kwargs will be set to None (except e).

         2. a yieldSurface object as first positional argument and additional
            keyWordArguments. If a parameter is given in both, the value from
            kwargs will be set. If a parameter is missing either in ysObject
            and kwargs it will be set to None (except e).

        The second method can act like a copy constructor or a converter method
        for closely related yieldSurface objects.

        @kwarg c: Mohr-Coulomb cohesion (p, tau-space)
        @kwarg phi: friction angle in deg (p, tau-space)
        @kwarg e: deviatoric eccentricity, range = [.5 - 1.]
        @note: UCSi and mb will be set internally and can not be changed directly.
            If you need to do so, convert the QuasiMohrCoulomb into an
            L{ExtendedMenetreyWillam} before:
              >>> ys = ExtendedMenetreyWillam(QuasiMohrCoulomb(c=c,phi=phi))
              >>> ys.mb *= .5
        @note: You may want to change s to alter the UCS-rockmass value. Take care
            that it does not become zero. Choose a (very) small number instead.
        """
        params = dict()
        if len(args) == 1:
            #params will probably get updated via kwargs!
            for key in ['phi', 'c', 'e']:
                try:
                    params[key] = getattr(args[0], key)
                except AttributeError:
                    pass

        elif len(args) > 1:
            raise ValueError("Only one postional argument (a yield surface)"
                             " allowed.")

        params.update(kwargs)

        self.phi = params.get('phi', None)
        self.c   = params.get('c', None)
        self.e = params.get('e', .6)

        # fixed parameters for QuasiMohrCoulomb
        self.s  = 1.
        self.a  = 1.

        #run _postInit to update other data
        self._postInit()


    @property
    def UCSi(self):
        """UniaxialCompressiveStrength of intact rock. Will be derived from
        stored friction angle and cohesion (see L{frictionAngleCohesionToUCSmb}).
        """
        if not (self.phi is None or self.c is None):
            return self.frictionAngleCohesionToUCSmb(self.phi, self.c)[0]

    @UCSi.setter
    def UCSi(self, value):
        raise AttributeError("UCSi can't be set directly. Change cohesion"
                             " instead or convert QuasiMohrCoulomb into an"
                             " ExtendedMenetreyWillam")


    @property
    def mb(self):
        """Frictional strength parameter. Will be derived from stored friction
         angle and cohesion (see L{frictionAngleCohesionToUCSmb}).
        """
        if not (self.phi is None or self.c is None):
            return self.frictionAngleCohesionToUCSmb(self.phi, self.c)[1]
    @mb.setter
    def mb(self, value):
        raise AttributeError("mb can't be set directly. Change friction angle"
                             " instead or convert QuasiMohrCoulomb into an"
                             " ExtendedMenetreyWillam")

    @staticmethod
    def frictionAngleCohesionToUCSmb(phi, c, minUCSi=1.):
        """Converts the friction angle and cohesion to UCSi and mb of an
        equivalent (for theta=pi/3) extendedMenetreyWillam yield surface having
        alpha = 1 and s = 0. For derivation see U{BE-wordpress<http://192.168.31.25/wordpress/archives/2028>}.

        @param phi: MohrCoulomb friction angle in degree

        @param c: MohrCoulomb Cohesion in Pa

        @returns: UCSi, mb
        """
        from math import sin,cos, pi

        phi = phi*pi/180.
        UCSi = 2*c*cos(phi)/(1.-sin(phi))
        mb = 2*sin(phi)/(1.-sin(phi))
        UCSi = max(UCSi, minUCSi)
        return UCSi, mb

    @staticmethod
    def cohsFrictionfromExtendedMenetreyWillamFit(extMW,
                                    fittingRange=(None,None),
                                    constrainedCohesion=None,
                                    nSamples=256, returnFit=True):
        """Returns the best fitting MohrCoulomb cohesion/friction angle pair
        for a given extended MenetreyWillam yield surface with respect the
        specified pressure range.

        Complete workflow to fit a given ExtendedMenetreyWillam-YieldSurface:
            >>> from bae.material_01 import (ExtendedMenetreyWillam,
            ...                              QuasiMohrCoulomb)
            >>> extMW = ExtendedMenetreyWillam.fromHoekBrown(50E6,60,0.9)
            >>> c,Phi = QuasiMohrCoulomb.cohsFrictionfromExtendedMenetreyWillamFit(
            ...                      extMW, fittingRange=(None,40E6),
            ...                      returnFit=False)
            >>> qMW = QuasiMohrCoulomb(phi=phi,c=c)


        @param extMW: ExtendedMenetreyWillam yieldSurface object or dictionary
            holding all required parameters for a extMW

        @param fittingRange: range of pressure where the fitting will be
            applied. If the lower (first) bound is None, the isotropic tensile
            strength of extMW will be used. If the upper (second) bound is None
            -50 times of isotropic tensile strength of extMW will be used.
            Note that only p-values above the isotropic tensile strength can
            be fitted. So if your full p-range is below the result will be
            useless.

        @param constrainedCohesion: forces the fit to have the given conhesion.
            The returned value for cohesion will equal the constrainedCohesion.

        @param nSamples: number of samples used over p-range.

        @param returnFit: if True a dictionary of
            n : array of fitted normal stress (dim = nSamples-1)
            tauMW : array of shearStress of extMW (dim = nSamples-1)
            tauMC : array of fittedStress of MohrCoulomb model (nSamples-1)

        @returns: cohesion,frictionAngle,(fitDictionary if returnFit=True)
        @note: all units SI, friction angle in deg
        """
        if isinstance(extMW, dict):
            extMW = ExtendedMenetreyWillam(**dict)

        if not isinstance(extMW, ExtendedMenetreyWillam):
            raise ValueError('Expected ExtendedMenetreyWillam object or' +
                             ' dictinoary. Got %s'% type(extMW))

        lower, upper = fittingRange

        if lower is None:
            lower = extMW.getIsoTensileStrength()

        if upper is None:
            upper = -50 * extMW.getIsoTensileStrength()

        p = np.linspace(lower, upper, nSamples)
        theta = np.ones_like(p)*np.pi/3
        q = extMW.ptOnYieldSurf(p=p, theta=theta)

        valid = ~np.isnan(q)
        p = p[valid]
        q = q[valid]
        theta = np.ones_like(p)*np.pi/3.

        S1,S2,S3 = pQThetaToEigs(p, q, theta, miningNotation=True)

        n,tau = nTauFromBalmersMethod(S1, S3)

        if constrainedCohesion is None:
            polyFac = np.polyfit(n,tau,1)
            tanPhi, c = polyFac
        else:
            c = constrainedCohesion
            tau -= c # offset fixed cohesion to solve tau = tanPhi*n + c
            A = np.vstack([n, np.zeros_like(n)]).T
            tanPhi = np.linalg.lstsq(A, tau, rcond=None)[0][0]

        phi = np.arctan(tanPhi) * 180./np.pi

        if returnFit:
            tauFit = tanPhi * n + c
            return c, phi, dict(n=n, tauMW=tau, tauMC=tauFit)
        else:
            return c, phi



class AnisotropicExtendedMenetreyWillam(ExtendedMenetreyWillam):
    """Describes an anisotropic extendedMenetreyWillam yield surface."""

    _type = 'MenetreyWillamAniso'
    _name = 'MenetreyWillamAniso'

    #set free parameterKeys
    _yieldParameters = ['UCSi', 'a', 's', 'e', 'mb', 'anisoS', 'anisoN']

    _freeYieldParameters = _yieldParameters

    # posible input data formats for fit choose from here:
    #   ['full tensor','full eigs','Sminor and Smajor','Smajor']
    _fitInputRequirements = ['full tensor',]

    def __init__(self, *args, **kwargs):
        """
        Creates an anisotropic ExtendedMenetreyWillam yield surface object.
        See L{ExtendedMenetreyWillam} for more details.

        @kwarg UCSi: uniaxial compressive strength of intact rock
        @kwarg mb: frictional strength parameter
        @kwarg s: degradation parameter
        @kwarg a: Hoek-Brown exponent, default=0.5
        @kwarg e: eccentricity, default = 0.6
        @kwarg anisoS: first anisotropic constant, default = 1.
        @kwarg anisoN: second anisotropic constant, default = 1.
        """
        #set related MenetreyWillam Parameters
        if args:
            anisoS = getattr(args[0], 'anisoS', 1.)
            anisoN = getattr(args[0], 'anisoN', 1.)
        else:
            anisoS, anisoN = 1.,1.

        if not 'anisoS' in kwargs:
            kwargs['anisoS'] = anisoS

        if not 'anisoN' in kwargs:
            kwargs['anisoN'] = anisoN

        super(AnisotropicExtendedMenetreyWillam, self).__init__(*args, **kwargs)


    @classmethod
    def _yieldSurfaceFun(cls, p, q, theta, UCSi=None, mb=None, s=None,
                        e=.6, a=.5, anisoN=None, anisoS=None):
        """MENETRY-WILLAM-YieldSurface, expressed in pressure, MISESstress
        and LODE-angle. Parameters need to be passed as keyword arguments
        @note: this is exactly the same as L{ExtendedMenetreyWillam}.
            However,anisoN and anisoS are in parameter list just for convenience

        """
        UCSi += 0. # force UCSi to be a float
        a += 0.    # force a to be a float
        extc = cls.eccentricity(theta, e)
        try:
            f = (q / UCSi)**(1/a) +  mb * (extc*q  / (3*UCSi) - p / UCSi) - s
        except FloatingPointError:
            print('MisesStress q:', q)
            print('UCSi:', UCSi)
            raise FloatingPointError("MisesStress or UCSi is negative."
                                     " Can't calculate root.")
        return f


    @staticmethod
    def morphStressTensor(S, anisoN, anisoS, alpha=None, beta=None,
                          returnPQTheta=True, debug=False):
        """
        Deforms a stress tensor S (expressed in a arbitrary basis B) according
        to the transversal-anisotropic parameters n and s.

        @param S: Stress tensor, can be:
            - a dictionary holding S_11, S_22, S_33, S_12, S_23, S13
            - a full tensor(field), where the last two dims are 3x3
            - a voigt-vectorlike tensor(field), where the last dim is 6
        @param anisoN: anisotropic parameter
        @param anisoS: anisotropic parameter
        @param alpha: first anisotropic angle (acos(n*ez)). Can be None (no
            rotation is applied, beta has to be None as well), scalar or same
            shape as a tensorcomponent.
        @param beta: second anisotropic angle.  Can be None, scalar or same
            shape as a tensorcomponent.
        @param returnPQTheta: if True (default) the mechanical invariants
            will be returned instead of the eigenValues
        @note: The (transversal-)isotropic plane is z-normal.
        """
        # first: some conversion to allow different input types of stress
        # tensor
        if isinstance(S, dict):
            try:
                SS = np.zeros((len(S["S_11"]), 3, 3))
            except TypeError: # S["S_11"] is scalar
                SS = np.zeros(3,3)
            SS[...,0,0] = S["S_11"]
            SS[...,1,1] = S["S_22"]
            SS[...,2,2] = S["S_33"]
            SS[...,1,0] = S["S_12"]
            SS[...,2,1] = S["S_23"]
            SS[...,2,0] = S["S_13"]
        else:
            S = np.asarray(S)
            if len(S.shape) >= 2 and S.shape[-2:] == (3,3):
                SS = S
            elif S.shape[-1] == 6:
                S = np.zeros(S.shape[:-1] + (3,3,))
                SS[...,0,0] = SS[...,0]
                SS[...,1,1] = SS[...,1]
                SS[...,2,2] = SS[...,2]
                SS[...,1,0] = SS[...,3]
                SS[...,2,1] = SS[...,4]
                SS[...,2,0] = SS[...,5]
            else:
                raise ValueError("Shape of S has to be (...,3,3) or (...,6)")

        # second apply rotation if alpha and beta is provided
        if alpha is not None:
            alpha, beta = np.asarray(alpha), np.asarray(beta)
            if alpha.shape and SS.shape == (3,3):
                SS = np.tile(SS[np.newaxis,:,:], alpha.shape + (1,1))

            # Stress to voigt vector
            V = [SS[...,0,0], SS[...,1,1], SS[...,2,2],
                 SS[...,1,0], SS[...,2,1], SS[...,2,0]]
            # rotate
            v11, v22, v33, v12, v23, v31 =  rotVoigtTransversal(V, alpha, beta)
            # back to lower triangle
            SS[...,0,0], SS[...,1,1], SS[...,2,2] = v11, v22, v33
            SS[...,1,0], SS[...,2,1], SS[...,2,0] = v12, v23, v31


        SS[...,0,0] *= anisoN
        SS[...,1,1] *= anisoN
        SS[...,1,0] *= anisoN
        SS[...,2,1] *= anisoS
        SS[...,2,0] *= anisoS

        # eigvalsh calculates eigenVals from Lower triangle (UPLO='L')
        # eigenValues are already sorted when using eigvalsh
        E = np.linalg.eigvalsh(SS, UPLO='L')
        S3, S2, S1 = E[...,0], E[...,1], E[...,2]

        if not returnPQTheta:
            return S1, S2, S3
        else:
            return stressToPQTheta(S1, S2, S3, fullLodeAngle=False)


    def fitToStressData(self, S, dataType=None, parameters=None,
                        stressWeighted=True, miningNotation=False,
                        isotropicOnly=False):
        """
        Fits parameters of anisotripic yield Surface to yielding stress values
        using a least square fit. Compared to the 'regular' procedure the
        stress values will be transformed using the anisotropic transformation
        tensor. While this tensor depends on the parameters anisoN and anisoS,
        this transformation becomes part of the goalFunction.
        """
        try:
            from scipy.optimize import least_squares
        except ImportError:
            raise ImportError('This function reqires sciPy')

        if isotropicOnly:
            ys = ExtendedMenetreyWillam(self)
            return ys.fitToStressData(S, dataType=dataType,
                                      parameters=parameters,
                                      stressWeighted=stressWeighted,
                                      miningNotation=miningNotation)


        #S to np.array with shape Nx3x3
        S = _streesDictToTensor(S)
        #Haigh-Westergaard-Coordinates, needed for scale and weights
        lc = stressToPQTheta(S)

        ###Same Code as in raw ExtendedMenetreyWillam.fitToStressData>>>

        #determine parameters to fit and store in list to preserve order
        if parameters is None:  # use all
            pKeys = self._freeYieldParameters()
            pKeys.remove('UCSi')

        else:  # use specified
            if not all([p in self._yieldParameters for p in parameters]):
                raise KeyError('Some specified parameters are unknown or ' +
                               'not changeable for Class %s' % self._name)
            pKeys = parameters

        #dict holding all parameters (will be updated )
        mwParam = self.yieldParams()

        if stressWeighted:
            #weight absolut errors by 'StressNorm'
            weights = np.linalg.norm(lc[:,:-1], axis=1)
            weights = weights/weights.max()

        #initial squared error sum as scaling factor to unifiy convergence crit
        # use _yieldSurfaceFun here to pass arbitrary parameters
        scale = (self._yieldSurfaceFun(lc[:,0], lc[:,1], lc[:,2],
                                      **mwParam)**2).sum()

        ###<<< END Same Code as in raw ExtendedMenetreyWillam.fitToStressData

        morphS = lambda ns: self.morphStressTensor(S, ns[0], ns[1])

        if any([k in pKeys for k in ['anisoN', 'anisoS']]):
            #apply morphing of StressValues in each itStep
            invokeMorph = True
        else:
            #apply morphing of StressValues only once
            invokeMorph = False
            lc = morphS(self.anisoN, self.anisoS)


        def minFun(x):
            """ErrorFunction to get minimized """
            cDict = dict(zip(pKeys, x))  #zip(pKeys,x) will preserve paramOrder
            mwParam.update(cDict)

            if invokeMorph:
                lc = morphS(mwParam['anisoN'], mwParam['anisoS'])

            mwParam.pop('anisoN', None)
            mwParam.pop('anisoS', None)

            # error value for each stress sample
            # use _yieldSurfaceFun here to pass arbitrary parameters
            errors = self._yieldSurfaceFun(lc[:,0], lc[:,1], lc[:,2],
                                           **mwParam)
            errors = errors/scale

            if stressWeighted:
                return errors/weights
            else:
                return errors

        #StartValues
        params0 = [self.freeYieldParams()[k] for k in pKeys]
        #Bounds
        bounds = tuple(np.array([list(self.bounds[k]) for k in pKeys]).T)

        #run least squares
        params = least_squares(minFun,params0,
                                  ftol=1E-12,
                                  xtol=1E-12,
                                  loss='soft_l1',
                                  bounds=bounds
                                  )

        return dict(zip(pKeys, params.x)), params.cost


    def getOrientedStrength(self, beta, alpha=None, S11=0., S22=None,
                            startValNewton=None):
        """Calculates the strength (S33) of an transversal anisotropic material-
        yieldsurface under a given confinement S11 and S22 with respect to
        a given orientation.
        @param beta: angle angle describing the normal vector given with:
            M{S{beta} = acos(nz)} in degree
        @param alpha: angle between x-axis and xy-projected normal vector
        @param S11: confinement S11 in Pa
        @param S22: confinement S22 in Pa, if None S22=S11
        @param startValNewton: S33 value to start newton iteration with. If
            None (default), it will be set to -100*UCS.
        @returns: strength S3 or np.nan if newton algorithm didn't converge
        @note: all values in mechanical (not mining) notation
        """
        try:
            from scipy.optimize import newton
        except ImportError:
            raise ImportError('This fuction requires Scipy.')

#        alpha = np.deg2rad(alpha)
#        beta = np.zeros_like(alpha)
#
        beta = np.deg2rad(beta)
        if alpha is None:
            alpha = np.zeros_like(beta)
        else:
            alpha = np.deg2rad(beta)

        if S22 is None:
            S22 = S11
        S11, S22 = np.asarray(S11), np.asarray(S22)
        if S11.shape:
            SS = np.zeros(S11.shape + (3,3,))
        elif alpha.shape:
            SS = np.zeros(alpha.shape + (3,3,))
        else:
            SS = np.zeros((3,3))

        SS[...,[0,1,], [0,1,]] = np.array([S11, S22]).T


        def evalFun(S33, debug=False):
            S = SS.copy()
            # set confinments:
            S[...,2,2] = S33

            # apply anisotropic stresstensor deformation
            p, q, theta = self.morphStressTensor(S,
                                            self.anisoN, self.anisoS,
                                            alpha, beta,
                                            returnPQTheta=True,
                                            debug=False)
            if debug:
                print('alpha', np.rad2deg(alpha))
                print('p',p)
                print('q',q)
                print('theta',np.rad2deg(theta))
            ysVal = self.yieldSurfaceFun(p,q,theta)
            return ysVal

        try:
            if S11.shape:
                start = np.full(S11.shape, -100*self.getUCS())
            elif alpha.shape:
                start = np.full(alpha.shape, -100*self.getUCS())
            else:
                start = -100*self.getUCS()
            strength = newton(evalFun, start)

        except RuntimeError:
            strength = np.nan
        return strength


class AnisotropicQuasiMohrCoulomb(AnisotropicExtendedMenetreyWillam,
                                  QuasiMohrCoulomb):
    """An anisotropic ExtendedMenetreyWillam that mimics the MohrCoulomb
    behavoir.

    Example:
        >>> from bae.material_01 import PlasticLevkovitchReuschAniso
        >>> from bae.material_01.yieldsurfaces import AnisotropicQuasiMohrCoulomb
        >>> from bae.material_01.material import ElasticLinear
        >>> faultProps = {'FAULTA' : [#[plasticStrain, cohesion, friction]
        ...                         [0, 1000, 36.],
        ...                         [.05, 100, 32.]
        ...                         ],
        ...              'FAULTB' : [[0, 1000, 36.],],}
        >>> defaultElast = ElasticLinear(E=10E9, nu=.25)
        >>> defaultDilation = 0.1
        >>> for faultName, props in faultProps.iteritems():
        >>>     fault = PlasticLevkovitchReuschAniso(
        >>>                 [(eps, defaultElast,
        ...                   AnisotropicQuasiMohrCoulomb(c=c,phi=phi),
        ...                   defaultDilation)
        ...                   for eps,c,phi in props])
        >>> # see material_01.matpropscollection for ABAQUS.inp-export
    """
    _name = 'AnisotropicQuasiMohrCoulomb'

    #free parameterKeys
    _freeYieldParameters = ['UCSi', 'mb', 'anisoS', 'anisoN']

    # posible input data formats for fit choose from here:
    #   ['full tensor','full eigs','Sminor and Smajor','Smajor']
    _fitInputRequirements = ['full tensor',]


    def __init__(self, *args, **kwargs):
        """Creates an AnisotropicQuasiMohrCoulomb yield surface object.
        See L{QuasiMohrCoulomb} as well.

        @kwarg c: Mohr-Coulomb cohesion (p, tau-space)
        @kwarg phi: friction angle in deg (p, tau-space)
        @kwarg e: deviatoric eccentricity, range = [.5 - 1.], default is 0.6
        @kwarg anisoS: first anisotropic constant, default = 1.
        @kwarg anisoN: second anisotropic constant, default = 1.
        """
        params = dict()
        if len(args) == 1:
            #params will probably get updated via kwargs!
            for key in ['phi', 'c', 'e', 'anisoS', 'anisoN']:
                try:
                    params[key] = getattr(args[0], key)
                except AttributeError:
                    pass

        elif len(args) > 1:
            raise ValueError("Only one postional argument (a yield surface)"
                             " allowed.")

        params.update(kwargs)

        self.phi = params.get('phi', None)
        self.c   = params.get('c', None)
        self.e   = params.get('e', .6)
        self.anisoS = params.get('anisoS', 1.)
        self.anisoN = params.get('anisoN', 1.)

        # fixed parameters for QuasiMohrCoulomb
        self.s  = 1.
        self.a  = 1.

        #run _postInit to update other data
        self._postInit()



class SinglePlaneMohrCoulomb(object):
    """This Class represents a SinglePlaneMohrCoulomb-YieldSurface. By now it
    does not support any fancy math-opperations/fittings.
    """
    _type = 'SinglePlaneMohrCoulomb'
    _name = 'SinglePlaneMohrCoulomb'

    # complete parameter-keys for current yield-surface._type
    _yieldParameters = ['friction', 'cohesion']

    # (free) changeable parameters
    _freeYieldParameters = _yieldParameters
    # posible input data formats for fit choose from here:
    #   ['full tensor','full eigs','Sminor and Smajor','Smajor']
    _fitInputRequirements = []


    def __init__(self, *args, **kwargs):
        """Creates a SinglePlaneMohrCoulomb object. Inputs could be a single
        SinglePlaneMohrCoulomb object (copy constructor) or kwargs 'friction'
        and 'cohesion'.

        @kwarg friction: friction-coeffitient for Mohr-Coulomb
        @kwarg cohesion: cohesion-coeffitient for Mohr-Coulomb
        """
        if len(args) == 1 and not kwargs:
            for key in self._freeYieldParameters:
                setattr(self, key, getattr(args[0], key, None))

        elif not args:
            for key in self._freeYieldParameters:
                setattr(self, key, kwargs.get(key, None))

        else:
            raise ValueError("Only one positional argument (should be a yield"
                         " surface object) allowed, but %d given." % len(args))

    def copy(self):
        params = dict((param, getattr(self,param))
                      for param in self._freeYieldParameters)
        return self.__class__(**params)


    def yieldParams(self):
        """Returns dictionary with parameters of current yieldSurface"""
        return dict((k, getattr(self,k)) for k in self._yieldParameters)


    def freeYieldParams(self):
        """Returns dictionary with free parameters of current yieldSurface"""
        return dict((k, getattr(self,k)) for k in self._freeYieldParameters)


    @staticmethod
    def _frictionAngleToFriction(f):
        """Calculates related mb from frictionAngle (given in deg).
        See also U{BE-wordpress<http://192.168.31.25/wordpress/archives/2028>}
        """
        return QuasiMohrCoulomb._frictionAngleToMb(f)

#} #End YieldSurface-Classes

###############################################################################


def getPAndQFromThetaAndS(ys, theta=np.deg2rad(60), **kwargs):
    """Returns pressure and vonMISES stress at a given eigenstress value and
    lode angle. This can be usefull if you need to pre-estimate the p-range
    up to a corresponding eigen value. e.g. you want to plot the yieldsurface
    at the compression meridian (theta=60deg) up to a given Sminor value:

    Example:
     >>> ys = ExtendedMenetreyWillam().fromComputator(50E6, 55, sample=1)
     >>> SminorMax = 100E6
     >>> theta = np.deg2rad(60)
     >>> pMax,qMax = getPAndQFromThetaAndS(ys, theta,
     ...                                   S1=SminorMax, miningNotation=True)
     >>> pRange = np.linspace(ys.getIsoTensileStrength(), pMax)
     >>> theta = np.full_like(pRange, theta)
     >>> q = ys.ptOnYieldSurf(p=pRange, theta=theta)
     >>> S1,S2,S3 = pQThetaToEigs(pRange,q,theta,miningNotation=True)

     @param ys: yieldSurface object
     @param theta: LODE angle(s) to be evaluated
     @kwarg S1,S2,S3: eigenStress(es) to get q and p for
     @kwarg miningNotation: S1,S2 or S3 given in mining notation
     @kwarg others: other keyWordArguments will be passed to scipys newton
         solver.
    """
    from scipy.optimize import newton
    miningNotation = kwargs.pop('miningNotation', False)

    if 'S1' in kwargs:
        S = np.asarray(kwargs.pop('S1'))
        idx = 0
    elif 'S2' in kwargs:
        S = np.asarray(kwargs.pop('S2'))
        idx = 1
    elif 'S3' in kwargs:
        S = np.asarray(kwargs.pop('S3'))
        idx=2
    else:
        raise ValueError("You have to pass EigenStressValues as kwarg"
                         " 'S1', 'S2' or 'S3'")

    shape = S.shape
    S = S.ravel()
    if isinstance(theta, (int,long,float)):
        theta = np.full(shape, theta)

    def minFun(q):
        if (q < 0).any():
            print q[q<0]
        try:
            p = ys.ptOnYieldSurf(q=q, theta=theta)
        except:
            print q, theta
            raise
        rr = pQThetaToEigs(p,q,theta,miningNotation=miningNotation)
        return rr[idx] - S

    # Mostly abs(S) is a good value except for stresstates nearby tensile
    # strength. We will try different startvalues here. Therefore we will catch
    # Runtime and FloatingPointErrors (which usualy only raise Warnings).
    with np.errstate(all='raise'): # to force FloatingPointErrors
        for k in range(10):
            try:
                startVal = (10**k)*np.abs(S) + 1
                q = newton(minFun, startVal, **kwargs)
                e = None
                break
            except RuntimeError as e:
                pass
            except FloatingPointError as e:
                pass

    if e is not None:
        print('startvalue: ', startVal)
        for name, val in ys.yieldParams().iteritems():
            print(name, val)
        raise e
    p = ys.ptOnYieldSurf(q=q, theta=theta)

    return p.reshape(shape), q.reshape(shape)


def bestComputatorFit(S1, S3, S2=None, version='V14.3',
                      relError=True, miningNotation=True):
    """Returns the UCSi and GSI values that will result in
    PEAK-yieldsurface of a computator material that fits
    the given (sampled) target yield surface best.

    Example:
     >>> version = 'V14.3'
     >>> ysMC = QuasiMohrCoulomb(c=100E3, phi=35.)
     >>> S1,_,S3 = ysMC.sampleMeridianEigs(S3Max=11E6,
     ...                                   distribution='lin',
     ...                                   samples=512)
     >>> UCSi, GSI = bestComputatorFit(S1, S3, version=version)
     >>> ysFit = ExtendedMenetreyWillam().fromComputator(UCSi,GSI,
     ...                               version=version)

    @param S1: sampled highest Eigenstresses on target yieldsurface
    @param S3: sampled lowest Eigenstresses on target yieldsurface
    @param S2: sampled intermed Eigenstresses on target yieldsurface
        or None. If None, S2 is assumed to equal S3 (compression meridian)
    @param version: computator version
    @param relError: if True, the overall relative deviation gets minimized
    @param miningNotaion: if True, copression has a positive sign
    @note: the error is derived from the p-deviation
    """
    from scipy.optimize import least_squares

    if S2 is None:
        S2 = S3

    if miningNotation:
        S1,S2,S3 = -S3,-S2,-S1

    p, q, theta = stressToPQTheta(S1, S2, S3)
    def minFun(arg):
        UCSi, GSI = arg
        ys = ExtendedMenetreyWillam().fromComputator(
                    UCSi, GSI, sample='PEAK', version=version)

        pTest = ys.ptOnYieldSurf(q=q, theta=theta)
        valid = ~np.isnan(pTest)
        if relError:
            error = np.abs((p[valid] - pTest[valid]) / p[valid]).sum()
        else:
            error = np.abs((p[valid] - pTest[valid])).sum()

        return error
    p = least_squares(minFun, (50E6, 50.),
                      bounds=([1000, 5], [400E6, 100]),)
    return p.x



#{some manual test-functions

def _testBestComputatorFit():
    import matplotlib.pyplot as plt

    version = 'V14.3'
    ysMC = QuasiMohrCoulomb(c=100E3, phi=35.)
    S1,_,S3 = ysMC.sampleMeridianEigs(S3Max=11E6,
                                    distribution='lin',
                                    samples=512)

    UCSi, GSI = bestComputatorFit(S1, S3, version=version)
    ysFit = ExtendedMenetreyWillam().fromComputator(UCSi,GSI,
                                    version=version)
    S1fit,_,S3fit = ysFit.sampleMeridianEigs(S3Max=11E6,
                                    distribution='alpha',
                                    samples=128)

    plt.close("all")
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(S3, S1, label='ref')
    ax1.plot(S3fit, S1fit,
             label='fit (%.1fMPa, %.1f)' % (UCSi/1.E6, GSI))
    ax1.legend()
    ax1.grid(True)
    plt.show()


def _test_3d():
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from bae.plot_01.mpl_core import make_get_proj

    miningNotation = False

    ys = ExtendedMenetreyWillam().fromComputator(50E6, 55, sample=1, a=.55)

    tensStr = ys.getIsoTensileStrength()
    pVs = np.linspace(tensStr, 100*abs(tensStr), 32)
    tVs = np.linspace(0, 2*np.pi, 128+1)
    ps,ts = np.meshgrid(pVs, tVs, indexing='ij')
    shape = ps.shape

    #find rhos for zeros of yieldSurf
    qs = ys.ptOnYieldSurf(p=ps.ravel(),
                          theta=ts.ravel()).reshape(shape)


    S1,S2,S3 = pQThetaToEigs(ps,qs,ts, miningNotation=miningNotation)

#
    tts = np.mod(ts, 2*np.pi) * 180/np.pi
    ii = np.logical_and( np.logical_and( S1>=S2, S2>=S3), qs>0)
    label = '%.1f <= theta <= %.1f' % (tts[ii].min(), tts[ii].max())

    plt.close("all")
    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')
    thetaScatter=ax1.scatter(S1[ii].ravel(),S2[ii].ravel(), S3[ii].ravel(),
                c=tts[ii].ravel(), label=label)
    ax1.plot_surface(S1, S2, S3, color='g')
    cb = fig.colorbar(thetaScatter, ax=ax1)
    cb.ax.set_ylabel('theta')
    if False:
        qEstFun = ys.fastQApproximator(pVs, includeTheta=True)
        qEst = qEstFun(ps.ravel(),ts.ravel()).reshape(shape)
        ax1.plot_surface(*pQThetaToEigs(ps,qEst,ts,miningNotation=miningNotation),
                         color='r')
#    fig.show()
#    return
    axLim = qs.max()
    ax1.plot([0,axLim],[0,0],[0,0], 'k')
    ax1.plot([0,0],[0,axLim],[0,0], 'k')
    ax1.plot([0,0],[0,0],[0,axLim], 'k')
    ax1.set_xlabel('$S_1$')
    ax1.set_ylabel('$S_2$')
    ax1.set_zlabel('$S_3$')

    ax1.get_Proj = make_get_proj(ax1, 1, 1, 1)
    ax1.set_aspect(1.0)
    ax1.legend()

    if True:
        # to see whats not working properly

        SS1a,SS1b = ys.ptOnYieldSurfEigen(S2=S2, S3=S3)
        SS2a,SS1b = ys.ptOnYieldSurfEigen(S1=S1, S3=S3)
        SS3a,SS1b = ys.ptOnYieldSurfEigen(S1=S1, S2=S2)

        ax1.plot_wireframe(SS1a, S2, S3, color='b')
        ax1.plot_wireframe(SS1b, S2, S3, color='b')
#        ax1.plot_wireframe(S1, SS2a, S3, color='brown')
#        ax1.plot_wireframe(S1, S2, SS3a, color='yellow')

    if False:
        qVs = np.linspace(0, qs.max(), 64)
        qqs,ts = np.meshgrid(qVs, tVs, indexing='ij')
        pps = ys.ptOnYieldSurf(q=qqs.ravel(), theta=ts.ravel())
        pps = pps.reshape(qqs.shape)
        ax1.plot_surface(*pQThetaToEigs(pps,qqs,ts,
                                        miningNotation=miningNotation),
                         color='pink')

    if False:
        pVs = np.linspace(0, ps.max(), 64)
        pps,ts = np.meshgrid(pVs, tVs, indexing='ij')
        qqs = ys.ptOnYieldSurf(p=pps.ravel(), theta=ts.ravel())
        qqs = qqs.reshape(pps.shape)
        ax1.plot_surface(*pQThetaToEigs(pps,qqs,ts,
                                        miningNotation=miningNotation),
                         color='cyan')
    fig.show()


def _test_2d():
    import matplotlib.pyplot as plt
    plt.close("all")
    fig = plt.figure()


    ax1 = fig.add_subplot(141)
    ax2 = fig.add_subplot(142)
    ax3 = fig.add_subplot(143)
    ax4 = fig.add_subplot(144)

    #ys = ExtendedMenetreyWillam().fromComputator(50E6, 55, sample=1, a=.55)

    ys = ExtendedMenetreyWillam(UCSi=40E6, e=.6, a=.5, s= 3.365E-4, mb=0.597)

    tensStr = ys.getIsoTensileStrength()
    p = np.linspace(tensStr, 50E6, 256)

    # 1st triaxial compression (HB-test)
    # theta = pi/3 (60deg) >> S2==S3
    thetaC = np.ones_like(p)*np.pi/3

    qC = ys.ptOnYieldSurf(p=p, theta=thetaC)

    S1,S2,S3 = pQThetaToEigs(p,qC,thetaC,miningNotation=True)

    sc = 1E6

#    ax1.set_xlim(right=S3.max()/sc)
    ax1.set_title(r'Triaxial compression ($\theta=60\deg$)')

    ax1.set_xlabel('$S_2=S_3$ / MPa')
    ax1.set_ylabel('$S_1$ / MPa')
    # ax1.set_aspect('equal')
    ax1.plot(S3/sc, S1/sc, label='LRX')
    ax1.plot(S3/sc, ys.hoekBrownYieldSurf(S3)/sc,
             label='Hoek-Brown', ls='--')
    ax1.legend()

    # 2nd - shear
    thetaS = np.full_like(p, np.pi/6)
    qS = ys.ptOnYieldSurf(p=p, theta=thetaS)

    ax2.set_title(r'"Shear" ($\theta=30\deg$)')
    ax2.set_xlabel(r'$p$ / MPa')
    ax2.set_ylabel(r'$q$ / MPa')
    # ax2.set_aspect('equal')
    ax2.plot(p/sc, qS/sc, label='LRX')
    # For theta n*60 and (n+1)*60 Hoek-Brown and Menetrey-Willam are equal
    # and between n*60 and (n+1)*60 Hoek-Browns yieldsurface is planar/linear.
    # Thus we can obtain the HB-vals for theta=30deg from the hight of
    # triangle between (Compression-Meridian) - (p-axes) - (extension-meridian)
    # iow: just build the mean vector(length) between both points:
    # meanvec = .5*(qz + cos(pi/3)qc)e1 + .5*sin(pi/3)qc*e2
    # with cos(pi/3)=.5 and sin(pi)=sqrt(3)/2
    # |mean| = sqrt( .25*(qz+.5qc)**2 + (3/16)*qc**2 )
    qZ = ys.ptOnYieldSurf(p=p, theta=np.zeros_like(thetaC))

    qShb = np.sqrt(.25*(qZ+.5*qC)**2 + (3/16.)*qC**2)

    ax2.plot(p/sc, qShb/sc, label='Hoek-Brown', ls='--')

    # testing MC-fit
    c,phi,fit = QuasiMohrCoulomb.cohsFrictionfromExtendedMenetreyWillamFit(
                                ys,
                                fittingRange=(None,10E6),
                                #constrainedCohesion=50E3,
                                )
    print('cohesion %.3f kPa:' % (c/1000.))
    print('friction angle %.2f deg' % phi)
    n = fit['n']
    tauMC = fit['tauMC']
    tauMW = fit['tauMW']

    ax3.set_title(r'MC approximation $\tau/n$')
    ax3.set_xlabel(r'$\sigma_n$ / MPa')
    ax3.set_ylabel(r'$\tau{}$ / MPa')

    ax3.plot(n/sc, tauMW/sc, label='LRX')
    ax3.plot(n/sc, tauMC/sc, label='fitMC')
    # ax3.plot(pMC/sc, qYsMC/sc, label='derivedQMC')
    ax3.legend()

    ysMC  = QuasiMohrCoulomb(c=c,phi=phi)

    thetaMC = np.ones_like(p) * np.pi/3.

    qMC = ysMC.ptOnYieldSurf(p=p, theta=thetaMC)
    S1MC,_,S3MC = pQThetaToEigs(p,qMC,thetaMC,miningNotation=True)

    ax4.set_title(r'MC approximation $S_1/S_3$')
    ax4.set_xlabel(r'$S_3$ / MPa')
    ax4.set_ylabel(r'$S_1$ / MPa')
    ax4.plot(S3/sc, S1/sc, label='LRX')
    ax4.plot(S3MC/sc,S1MC/sc, label='fitMC', ls='--')
    ax4.legend()

    fig.show()
    color = 'tab:red'
    fig, axs = plt.subplots(1,2)
    ax1 = axs[1]
    ax1.set_xlabel(r'upper bound fitting range ($S_3$) / MPa')
    ax1.set_ylabel('phi / deg', color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color2 = 'tab:blue'
    ax2.set_ylabel('cohesion / kPa', color=color2)
    ax2.tick_params(axis='y', labelcolor=color2)

    ax3 = axs[0]
    #ax3.axis('equal')
    ax3.set_title(r'MC approximation $\tau/n$')
    ax3.set_xlabel(r'$\sigma_n$ / MPa')
    ax3.set_ylabel(r'$\tau{}$ / MPa')
    ax3.grid(True)

    sc = 1E6
    cs,phis = [],[]
    fRanges = range(1,51)[::-1]
    for ii,S3Max in enumerate(fRanges):
        S3Max *= sc
        c,phi,fit = QuasiMohrCoulomb.cohsFrictionfromExtendedMenetreyWillamFit(
                                ys,
                                fittingRange=(None,S3Max),
                                #constrainedCohesion=50E3,
                                )
        cs.append(c/1000.)
        phis.append(phi)
        n = fit['n']
        tauMC = fit['tauMC']
        tauMW = fit['tauMW']
        if ii == 0:
            ax3.plot(n/sc, tauMW/sc, label='LRX')
        if not ii % 10:
            ax3.plot(n/sc, tauMC/sc, label='[0,%dMPa]' % (S3Max/1E6), ls='--')
    ax3.legend()
    ax3.set_ylim([0,15])


    ax1.plot(fRanges, phis, color=color)
    ax2.plot(fRanges, cs, color=color2)

    fig.show()

def _test_plotSimple():
    ys = ExtendedMenetreyWillam().fromComputator(50E6, 55, sample=1, a=.55)
    import matplotlib.pyplot as plt
    plt.close("all")
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    S1,_,S3 = ys.sampleMeridianEigs(S3Max=10E6, S1Max=5E7,
                                    distribution='log',
                                    samples=32)
    ax1.scatter(S3, S1, s=5)
    fig.show()

def testHoekBrown():
    pass
#} #End some manual test-functions
###############################################################################

def testStrangeProps():
    ys = AnisotropicExtendedMenetreyWillam(
                # UCSi = 1.282e+08,
                # a = 5.025e-01,
                # s = 2.521e-02,
                # e = 6.000e-01,
                # mb = 2.776e+00,
                # anisoS = 1.000e+00,
                # anisoN = 1.000e+00,
                #
                UCSi = 1.282000e+08,
                a = 5.024595e-01,
                s = 2.521359e-02,
                e = 6.000000e-01,
                mb = 2.775891e+00,
                anisoS = 1.000000e+00,
                anisoN = 1.000000e+00,
        )
    #ys.ptOnYieldSurf(q=[50000001.,], theta=[1.04719755,])

    S1,_,S3 = ys.sampleMeridianEigs()

if __name__ == '__main__':
    print("No syntax errors")
    # _test_3d()
    # _test_2d()
    # _test_plotSimple()
    # _testBestComputatorFit()
    testStrangeProps()
