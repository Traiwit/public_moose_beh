# -*- coding: utf-8 -*-

"""
Package to manage/analyse materials.

It is structured into the following modules:
    - models: heuristic models for material properties (e.g. computator)
    - yieldsurfaces: structures and mathematical functions for yield surfaces
    - material: structures for 'complete' materials (elastic+plastic behavior
      for rock and fill, ...)
    - matpropcollection: provides a more user friendly layer for material
      modules (--> api for graphical beRocker-tool)

Using material_01.yieldsurfaces < handling yield surfaces >:
============================================================
    >>> from bae.material_01.yieldsurfaces import ExtendedMenetreyWillam
    >>> from bae.material_01.yieldsurfaces import hWCoorToEigs
    >>> from mpl_toolkits.mplot3d import Axes3D
    >>> import matplotlib.pyplot as plt
    >>> from matplotlib import cm
    >>> import numpy as np
    >>>
    >>>
    >>> ## create yieldsurface object
    >>> ySurf = ExtendedMenetreyWillam(s=0.005, mb=.8, UCSi = 50E6,)
    >>> ## get some 'special' values
    >>> print 'uniaxialShearStrengs: %.1f'%ySurf.getUniaxialShearStrength()
    >>> print 'UCS: %.1f'%ySurf.getUCS()
    >>> print 'isotropicTensileStrength %.1f'%ySurf.getIsoTensileStrength()
    >>>
    >>> ## zerofinding for points on yieldsurface (here: given p and theta)
    >>> # create a structured grid for p and theta
    >>> pVs = np.linspace(ySurf.getIsoTensileStrength(), .1*ySurf.UCSi, 36)
    >>> tVs = np.linspace(0, 2*np.pi, 36+1)
    >>> ps,ts = np.meshgrid(pVs, tVs)
    >>> shape = ps.shape
    >>>
    >>> # find rhos for zeros of yieldSurf
    >>> qs = ySurf.ptOnYieldSurf(xi=ps.ravel(), theta=ts.ravel())
    >>>
    >>> # retransform xi,rho,theta to eigStresses
    >>> s = hWCoorToEigs(ps.ravel(), qs.ravel(), ts.ravel() ).T
    >>> s = s/1E6
    >>> # remesh eigStresses
    >>> S1 = s[:,0].reshape(shape)
    >>> S2 = s[:,1].reshape(shape)
    >>> S3 = s[:,2].reshape(shape)
    >>>
    >>> # plot yieldsurface
    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111, projection='3d')
    >>> # color is rho
    >>> R  = qs.reshape(shape)
    >>> col= cm.YlGnBu_r((R-R.min())/(R.max()-R.min()))
    >>> ax.plot_surface(S1, S2, S3, facecolors=col)


Usage material_01 < handling single materials >:
================================================

L{AbqVumatData} is a container to hold all material props for one material
defined in a *MATERIAL card in the Abaqus input file: Elastic props,
density, damping alpha, data of the *DEPVAR and *USER MAT cards.

The data can be collected from an L{bae.abq_model_02.Model}-object through
L{AbqVumatData}.__init__() and can also be stored back into an
L{bae.abq_model_02.Model}-object through its
L{updateAbqModel<AbqVumatData.updateAbqModel>} method.


Reading and interpreting data from Abaqus input file:
-----------------------------------------------------

Note that it is strongly recommended to use L{MatPropsCollection} for
every-days abaqus-input creation (see examples below).

The raw data can be read from Abaqus input file through
L{bae.abq_model_02.Model} and will be stored in its
L{material<bae.abq_model_02.Model.material>} dictionary.
Then the raw data can be interpreted by L{AbqVumatData}.__init__():
 >>> from bae.abq_model_02 import Model
 >>> from bae.material_01 import AbqVumatData
 >>>
 >>> model = Model().read("myInputFile.inp")
 >>> rawMat = model.material["HOST"]
 >>> print rawMat["Density"], rawMat["Umat"]
 ... # density and raw values of the *USER MAT card
 (2700.0, [58000000.0, 0.6, 2.0, ...])
 >>> mat = AbqVumatData(rawMat)  # interpreting raw data
 >>> rockPeakUcs = mat.plasticLRRock[0].yieldSurf.UCSi
 >>> rockPeakBulkModulus = mat.plasticLRRock[0].elastic.K


Syntesizing material props and storing in Abaqus input file:
------------------------------------------------------------

Note that it is strongly recommended to use L{MatPropsCollection} for
every-days abaqus-input creation (see examples below).

 >>> from bae.material_01 import AbqVumatData
 >>> from bae.material_01 import PlasticLevkovitchReusch as PLR
 >>> from bae.material_01 import PlasticLevkovitchReuschAniso as anisoPLR
 >>> from bae.material_01.yieldsurfaces import ExtendedMenetreyWillam
 >>> from bae.material_01 import ElasticLinear
 >>> from bae.abq_model_02 import Model
 >>>
 >>> ## create single material from scratch
 >>> ucsi = 60E6
 >>> elasticPeak = ElasticLinear(E=14.93E9, nu=0.3)
 >>> elasticRes = ElasticLinear(E=4.93E9, nu=0.3)
 >>> emwYSPeak = ExtendedMenetreyWillam(mb=1, s=.1, UCSi=ucsi)
 >>> emwYSRes = ExtendedMenetreyWillam(mb=.5, s=.01, UCSi=ucsi)
 >>> pLR0 = PLR([#(pst, elast, yield, dilation),
 ...            (0., elasticPeak, emwYSPeak, .2),
 ...            (0.025, elasticRes, emwYSRes, .1),
 ...            ])
 >>>
 >>> ## create a single material using Computator
 >>> pLR1 = PLR().fromComputator(50E6, 55)
 >>>
 >>> ## comparing materials
 >>> isEqual, deltas = pLR0.compare(pLR1)
 >>>
 >>> ## convert to anisoPLR
 >>> pLRAniso = pLR1.convert(newType = anisoPLR,
 ...                        addYieldSurf = dict(anisoS=.1,
 ...                                            anisoN=.5,) )
 >>>
 >>> ## storing to model
 >>> model = Model()
 >>> mat = AbqVumatData()
 >>> mat.name = 'TEST'
 >>> mat.density = 2700.
 >>> mat.dampingAlpha = 0.5
 >>> mat.umatDelete = 9
 >>> mat.umatDepvarNum = 16
 >>> mat.voidStiffnessRatio = 1e-5
 >>> mat.critDamage = 1000.0
 >>> mat.plasticLRRock = pLR1 #Rock
 >>> mat.plasticLRFill = pLR0 #Fill
 >>> mat.updateAbqModel(model)
 >>> model.write("Check_material_M01.inp")
 >>>
 >>> ## reading from abaqus input file (for testing...)
 >>> mat1 = AbqVumatData(inputFile="Check_material_M01.inp",
 ...                     materialName="TEST")


Usage matpropcollection < handling multiple materials >:
========================================================
    >>> from bae.material_01 import (MatPropsCollection,
    ...                              PlasticLevkovitchReuschAniso)

    importing defaults for general properties (e.g. voidStiffness, ...) and
    a MatPropsCollection with some default fill properties

    >>> from bae.material_01.defaults import generalProps, defaultFillAniso

    creating a new matPropsCollection

    >>> matProps = MatPropsCollection('M01')
    >>> rockData = [#name    UCS   GSI  mr(V14)
    ...     ('matA', 40E6,  60, 600),
    ...     ('matB', 120E6, 60, 400),
    ...     ('matC', 80E6,  70, 600),
    ...    ]
    >>> # adding a default fill from 'defaultFillAniso'-MatPropsCollection
    >>> fill = defaultFillAniso['FILL_PLAST_100xx']
    >>> matProps.fill.addItem('defaultFill', fill, **generalProps)
    >>> # lopping all rock materials
    >>> for name, ucs, gsi, mr in rockData:
    >>>     rock = PlasticLevkovitchReuschAniso()
    >>>     rock.fromComputator(ucs,gsi,mr=mr,version='V14')
    >>>     matProps.rock.addItem(name, rock, **generalProps)
    >>>     # assign current rock material and defaultFill
    >>>     matProps.assignments.addAssignment(name, 'defaultFill')

    adding MohrCoulomb materials

    >>> from material_01 import AnisotropicQuasiMohrCoulomb, ElasticLinear
    >>> mcData = [#name    cohesion   frictionAngle
    ...     ('mcA', 40E4,  30),
    ...     ('mcB', 120E4, 32),
    ...    ]
    >>> elastic = ElasticLinear(E=10E9, nu=0.2)
    >>> dilation = 0.125
    >>> for name, cohesion, friction in mcData:
    >>>     ys = AnisotropicQuasiMohrCoulomb(c=cohesion, phi=friction)
    >>>     mc = PlasticLevkovitchReuschAniso([(0, elastic, ys, dilation)])
    >>>     matProps.rock.addItem(name, mc, **generalProps)
    >>>     matProps.assignments.addAssignment(name, 'defaultFill')

    export to abaqus input
     >>> matProps.exportAbaqusInput('test_matV14.inp',
     ...                            umatDelete=9, umatDepvarNum=16)

    export to excel
     >>> matProps.rock.toXlsx('testXlsx.xlsx')

    reading from abaqus input
     >>> matProps = MatPropsCollection('TestCollection')
     >>> ## reading from abaqus.inp
     >>> matProps.importAbaqusInput('testFile.inp')

    getting and setting parameters
     >>> ## access materials
     >>> rocks = matProps.rock
     >>> fills = matProps.fill
     >>> assignments = matProps.assignments
     >>> hostRock = rocks['MATV_HOST']
     >>> hostFill = fills['MATV_HOST']
     >>> # or
     >>> hostRock,hostfill = assignments['MATV_HOST'].getAssignedMaterials()

     >>> ## access parameters
     >>> #see all parameters
     >>> print hostRock.getPropertyDicts()
     >>> #see only (pst-dependent) parameters of 2nd sample
     >>> print hostRock.getPropertyDicts(idx=1)
     >>> # or
     >>> print rocks.getPropertyDicts('MATV_HOST')
     >>> print rocks.getPropertyDicts('MATV_HOST', idx=1)

     >>> ## edit parameters
     >>> hostRock.setProperty('mb', .1, idx=1)   # mb < specific for each sample
     >>> hostRock.setProperty('UCSi', 50E6)    # UCSi < general material property
     >>> try:
     >>>     hostRock.setProperty('e', .5, idx=1)
     >>>     #e is general but idx isn't None --> Error
     >>> except Exception as e:
     >>>     print str(e)

    editing MatItemCollections (rock and fill) and Assignments
     >>> # renaming matItems
     >>> # note: this will change all appearences of rocks 'MATV_HOST' in
     >>> # assignments as well
     >>> rock.rename('MATV_HOST', 'MATV_HORST')
     >>> assignments.renameAssignment('MATV_HOST', 'MATV_UWE')
     >>> rock.copyItem('MATV_HORST', 'MATV_HELGA')
     >>> rock.removeItem('MATV_HORST') # will also remove Assignment 'MATV_UWE'
"""

__version__ = "1.5"

_version_history_ = """
1.0 TR new
1.1 TR added matpropcollection ("user-layer" for material_01, still
    experimental state),
    added SinglePlaneMohrCoulomb, some minor fixes, reworked Inherence of
    PST(-list) classes
1.2 TR removed version numbers from submodules
1.3 TR MatPropsCollection now holds orderedDicts, fixed excel-docs
1.4 TR added copy function and copy constructors for yieldsurfaces and PST-
    samples, added defaults, large docs update, added FactorOfSafetyIso,
    added more testing in test_bae_package.test_material_01
1.5 TR added rotVoigtTransversal
"""

__todo__ = """
"""

from material import (
    ElasticLinear,
    PlasticLevkovitchReusch,
    PlasticLevkovitchReuschAniso,
    SinglePlaneMohrCoulomb,
    AbqVumatData,
    UmatLookUpError,
    )

from yieldsurfaces import (
    ExtendedMenetreyWillam,
    AnisotropicExtendedMenetreyWillam,
    QuasiMohrCoulomb,
    AnisotropicQuasiMohrCoulomb,
    pQThetaToEigs,
    stressToPQTheta,
    getPAndQFromThetaAndS,
    )

from matpropcollection import MatPropsCollection


# from models import (
#     computatorV12,
#     computatorV13,
#     hoekBrown,
#     )
