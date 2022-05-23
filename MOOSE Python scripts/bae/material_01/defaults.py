# -*- coding: utf-8 -*-
"""default pst-sample-lists for rock and fill

Usage with MatPropsCollection:
    >>> from material_01.defaults import defaultFillAniso, generalProps
    >>> from material_01 import MatPropsCollection
    >>> mPC = MatPropsCollection('M01')
    >>> # add a defaultFill material to your MatPropsCollection
    >>> mPC.fill['myFill'] = defaultFillAniso['FILL_PLAST_050xx']
    >>> # now adding rock materials to your MatPropsCollection
    >>> for name, UCSi, GSI in [('A',50,45), ('B', 70, 50)]:
    >>>     rockMat = PlasticLevkovitchReuschAniso()
    >>>     rockMat.fromComputator(UCS, GSI, version='V14.2')
    >>>     mPC.rock.addItem(name, rockMat, **generalProps)
    >>>     # assign rock material to your fill
    >>>     mPC.assignments(name, 'myFill')

"""
from material import (ElasticLinear,
                       PlasticLevkovitchReusch,
                       PlasticLevkovitchReuschAniso,)

from matpropcollection import MatItemCollection

from yieldsurfaces import (ExtendedMenetreyWillam,
                            QuasiMohrCoulomb,)

#: default general properties
generalProps = {
        "dampingAlpha" : 0.5,
        "density": 2700,
        "voidStiffnessRatio": 5E-6,
        "critDamage" : 1000.,
        "umatDelete":None,
        "umatDepvarNum":None,
}

#: default Rock materials
#:     no members yet
rock = {}


#: default Dump and Cave materials
#:     'DUMP_PLAST_MC_0100_32' - Mohr-Coulomb 100kPa cohesion and 32deg
dump = {
'DUMP_PLAST_MC_0100_32' : {
        'description' : ('One-pst-sample with quasi-Mohr-Coulomb material' +
                         ' with 100kPa cohesion and 32deg friction angle'),
        'type': 'QuasiMohrCoulomb',
        'samples': [{'pst' : 0.0,
                     'cohesion' : 100E3,
                     'phi' : 32.0,
                     'e' : 0.6,
                     'E': 0.2E9,
                     'nu': 0.2,
                     'dilation' : .25},
                    ],
    },
}



#: default Fill materials
#:     FILL_PLAST_050xx - UCSi=0.5MPa, tensile strength ~ UCSi/10
#:     FILL_PLAST_100xx - UCSi=1.0MPa, tensile strength ~ UCSi/10
#:     FILL_PLAST_250xx - UCSi=2.5MPa, tensile strength ~ UCSi/10
#:     FILL_ELAST_20MPA - standard elastic fill for support load evaluation
fill = {
    'FILL_PLAST_050xx' : {
        'description' : ('UCSi is 0.5MPa and tensile strength is approx'
                         ' 1/10th. No elastic softening'),
        'type': 'ExtendedMenetreyWillam',
        'samples': [{'pst' : 0.0,
                     'UCSi' : 0.5E6,
                     'a' : 0.2,
                     'e' : 0.6,
                     's' : 1.0,
                     'mb' : 8.0,
                     'E': 0.2E9,
                     'nu': 0.2,
                     'dilation' : .25},
                    {'pst' : 1.5E-2,
                     'UCSi' : 0.5E6,
                     'a' : 0.2,
                     'e' : 0.6,
                     's' : 0.001,
                     'mb' : 2.0,
                     'E': 0.2E9,
                     'nu': 0.2,
                     'dilation' : .25},
                    {'pst' : 3.5E-2,
                     'UCSi' : 0.5E6,
                     'a' : 0.2,
                     'e' : 0.6,
                     's' : 0.00001,
                     'mb' : 0.025,
                     'E': 0.2E9,
                     'nu': 0.2,
                     'dilation' : .25},
                    ],
    },
    'FILL_PLAST_100xx' : {
        'description' : ('UCSi is 1.0MPa and tensile strength is approx'
                         ' 1/10th. No elastic softening. Standard backfill'
                         ' for stopes'),
        'type': 'ExtendedMenetreyWillam',
        'samples': [{'pst' : 0.0,
                     'UCSi' : 1.0E6,
                     'a' : 0.2,
                     'e' : 0.6,
                     's' : 1.0,
                     'mb' : 8.0,
                     'E': 0.2E9,
                     'nu': 0.2,
                     'dilation' : .25},
                    {'pst' : 1.5E-2,
                     'UCSi' : 1.0E6,
                     'a' : 0.2,
                     'e' : 0.6,
                     's' : 0.001,
                     'mb' : 2.0,
                     'E': 0.2E9,
                     'nu': 0.2,
                     'dilation' : .25},
                    {'pst' : 3.5E-2,
                     'UCSi' : 1.0E6,
                     'a' : 0.2,
                     'e' : 0.6,
                     's' : 0.00001,
                     'mb' : 0.025,
                     'E': 0.2E9,
                     'nu': 0.2,
                     'dilation' : .25},
                    ],
    },
    'FILL_PLAST_250xx' : {
        'description' : ('UCSi is 2.5MPa and tensile strength is approx'
                         ' 1/10th. No elastic softening.'),
        'type': 'ExtendedMenetreyWillam',
        'samples': [{'pst' : 0.0,
                     'UCSi' : 2.5E6,
                     'a' : 0.2,
                     'e' : 0.6,
                     's' : 1.0,
                     'mb' : 8.0,
                     'E': 0.2E9,
                     'nu': 0.2,
                     'dilation' : .25},
                    {'pst' : 1.5E-2,
                     'UCSi' : 2.5E6,
                     'a' : 0.2,
                     'e' : 0.6,
                     's' : 0.001,
                     'mb' : 2.0,
                     'E': 0.2E9,
                     'nu': 0.2,
                     'dilation' : .25},
                    {'pst' : 3.5E-2,
                     'UCSi' : 2.5E6,
                     'a' : 0.2,
                     'e' : 0.6,
                     's' : 0.00001,
                     'mb' : 0.025,
                     'E': 0.2E9,
                     'nu': 0.2,
                     'dilation' : .25},
                    ],
    },
    'FILL_ELAST_20MPA' : {
        'description' : ('STANDARD elastic backfill for drives to show'
                         ' correct support load magnitudes.'),
        'type': 'ExtendedMenetreyWillam',
        'samples': [{'pst' : 0.0,
                     'UCSi' : 999.9E6,
                     'a' : 0.2,
                     'e' : 0.6,
                     's' : 1.0,
                     'mb' : 8.0,
                     'E': 0.02E9,
                     'nu': 0.2,
                     'dilation' : 1},
                    ],
    },
}


def _parseSample(sampleDict):
    """Parses a sampleDict from matDict['samples'] above into a quad-tuple
    (pst, elastic, yieldSurface, dialation) to set up one pst sample of
    PlasticLevkovitchReusch-Object.
    If a sampleDict has 'cohesion'-key it is assumed to be a MohrCoulomb
    sample which will get transformed in its ExtendedMeneteryWillam
    representation first.
    """
    sampleDict = dict(sampleDict)
    if 'cohesion' in sampleDict:
        cohs = sampleDict.pop('cohesion')
        phi = sampleDict.pop('phi')

        UCSi, mb = QuasiMohrCoulomb.frictionAngleCohesionToUCSmb(phi, cohs)
        sampleDict['s'] = 1.
        sampleDict['mb'] = mb
        sampleDict['UCSi'] = UCSi
        sampleDict['a'] = 1.

    E, nu = sampleDict.pop('E'), sampleDict.pop('nu')
    d = sampleDict.pop('dilation')
    pst = sampleDict.pop('pst')
    elastic = ElasticLinear(E=E, nu=nu)
    ys = ExtendedMenetreyWillam(**sampleDict)
    return (pst, elastic, ys, d)


def _addMaterial(collection, name, matData, toAniso=False):
    """Adds a material from dictionarys above to a collection.
    """
    pLR = PlasticLevkovitchReusch([_parseSample(sample)
                                   for sample in matData['samples']])
    if toAniso:
        pLR = pLR.convert(newType = PlasticLevkovitchReuschAniso,
                          addYieldSurf = dict(anisoS=1., anisoN=1.,) )

    collection.addItem(name, pLR, logText=matData['description'])


## create MatItemCollections for rock
#: MatItemCollection of isotropic default-rock-materials (empty yet)
defaultRockIso = MatItemCollection(name='default_rock_iso')
#: MatItemCollection of anisotropic default-rock-materials (empty yet)
defaultRockAniso = MatItemCollection(name='default_rock_aniso')
for _key, _data in rock.iteritems():
    _addMaterial(defaultRockIso, _key, _data, toAniso=False)
    _addMaterial(defaultRockAniso, _key, _data, toAniso=True)


## create MatItemCollections for Dump/cave
#: MatItemCollection of isotropic default-dump/cave-materials
defaultDumpIso = MatItemCollection(name='default_dump_iso')
#: MatItemCollection of anisotropic default-dump/cave-materials
defaultDumpAniso = MatItemCollection(name='default_dump_aniso')
for _key, _data in dump.iteritems():
    _addMaterial(defaultDumpIso, _key, _data, toAniso=False)
    _addMaterial(defaultDumpAniso, _key, _data, toAniso=True)


## create MatItemCollections for fill
#: MatItemCollection of isotropic default-fill-materials
defaultFillIso = MatItemCollection(name='default_fill_iso')
#: MatItemCollection of anisotropic default-fill-materials
defaultFillAniso = MatItemCollection(name='default_fill_aniso')
for _key, _data in fill.iteritems():
    _addMaterial(defaultFillIso, _key, _data, toAniso=False)
    _addMaterial(defaultFillAniso, _key, _data, toAniso=True)



