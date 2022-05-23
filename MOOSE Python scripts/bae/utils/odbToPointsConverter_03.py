"""Converter classes for odbToVtkGridBox and the other
L{odbToPointsData<bae.utils.odbToPointsData_03.odbToPointsData>} classes.

Usage:
======

For points on a regular grid:
 >>> from bae.utils.odbToPointsData_03 import odbToVtkGridBox
 >>> odbToVtkGridBox(...)

Or for unstructured points
 >>> from bae.utils.odbToPointsData_03 import odbToCsvPoints
 >>> odbToCsvPoints(...)
"""

__version__ = "3.01"

_version_history_ = """
3.01 GP new based on odbToPointsData_02 v2.35
  - split off from odbToPointsData
  - removed Cavesim stuff
"""

todo = """
- finish Converter_CombinedOutputOnSurface
- update Converter_StandardAndFOS
- update Converter_PlasticStrainTensor
- add converter function for eigenvalue decomposition (there is one for
  stress, is it desired for the general case?)
"""

import os
from math import log10
from itertools import izip
import zipfile
import time
import traceback

from bae.misc_01 import CheckFileFinished
from bae.field_01 import createFieldClass, createFieldObject
from bae.vtk_02 import FieldsOnStructPointGrid
from bae.log_01 import msg

from bae.material_01 import MatPropsCollection
from bae.material_01.material import FactorOfSafetyIso, FactorOfSafetyAniso

try:
    import numpy as np
except ImportError:
    np = None


class Converter(object):
    pass


class Converter_PlasticStrainTensor(Converter):
    """
    DEPRECATION WARNING:

    This class needs to be re-worked to get the material data from
    the post-data database. This interface is bound for a major change!

    Usage:
     >>> import os
     >>> from bae.abq_model_02 import Model
     >>> from bae.utils.odbToPointsData_03 import odbToVtkGridBox
     >>> from bae.utils.odbToPointsConverter_03 import \\
     >>>     Converter_PlasticStrainTensor as Converter
     >>>
     >>> #-----
     >>> boxName = "box2_mine_h06"
     >>> gridData = ...
     >>> ...
     >>>
     >>> #-----
     >>>
     >>> matNameToElProps = Converter.getMatNameToElPropsFromModel(
     >>>     Model().read("MyProject_material_M01.inp"),
     >>>     iPropsG=10, iPropsK=11)
     >>>
     >>> stateSDV = "SDV8"
     >>> def maskFieldFunc(data, mesh):
     >>>     return [(state<1.5) for state in data[stateSDV]]
     >>>
     >>> converter = Converter(
     >>>     matNameToElProps=matNameToElProps,
     >>>     stateSDV=stateSDV, maskFieldFunc=maskFieldFunc)
     >>>
     >>> odbToVtkGridBox(
     >>>     odbPath, frameList=frameList,
     >>>     gridName=boxName, gridData=gridData,
     >>>     outputDirName=outputDirName, fileNamePrefix=fileNamePrefix,
     >>>     projectPrefix=..., meshVersion=...,
     >>>     **converter.fieldArguments)
     >>>
     >>> # remove the matList.csv file
     >>> os.remove(os.path.join(outputDirName, fileNamePrefix+"_matList.csv"))
    """

    def __init__(
            self,
            matNameToElProps,
            fieldName="EP",
            stateSDV="SDV8",
            maskFieldFunc=(
                lambda data, mesh: [(state<1.5) for state in data["SDV8"]]),
            valueNotDefined=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]):
        """
        @param matNameToElProps: a dict {matName : (G, K)} giving shear and
           bulk modulus. See L{getMatNameToElPropsFromModel}.
        @param fieldName: Name of the resulting field
        @param stateSDV: Name of the state field in the odb.
        @param fieldName: Name of the result field
        @param maskFieldFunc: A function generating a mask field: True where
           the plastic strain shall be calculated and False where we want the
           constant valueNotDefined as output
        @param valueNotDefined: Take this value for points that are masked out
           by maskFieldFunc and for which we don't have results in the odb
           (for example above the surface).
        """
        self.fieldName=fieldName
        self.stateSDV = stateSDV
        self.maskFieldFunc = maskFieldFunc
        self.valueNotDefined = valueNotDefined
        self.fieldPos = None

        # dummy class with arbitrary position
        # will be corrected from first call to self.compose1()
        self.FieldClass = createFieldClass(
            self.fieldName, "structPt", "tensor")

        # prepare material props
        self.matNamesList = sorted(matNameToElProps)
        self.matNbToElProps = dict(
            (i+1, matNameToElProps[matName])
            for i, matName in enumerate(self.matNamesList))

    @staticmethod
    def getMatNameToElPropsFromModel(model, iPropsG=10, iPropsK=11):
        """Return the matNameToElProps dictionary required for the
        Converter_PlasticStrainTensor.__init__ function.

        @param model: a L{bae.abq_model_02.Model} object containing material
            data. E.g.
             >>> model = Model().read("Cadia2016033_material_M21.inp")

        @param iPropsG: zero based index of the shear modulus in the data lines
            of the Abaqus *USER MATERIAL key card
        @param iPropsK: zero based index of the bulk modulus
        """

        matNameToElProps = dict()  # {matName: (G, K)}
        for matName, mat in model.material.iteritems():
            try:
                G = mat['Umat'][iPropsG]
                K = mat['Umat'][iPropsK]
            except KeyError:
                raise KeyError(
                    "Material %s is no user material, other types not"
                    " implemented yet. (Note: this might be very easy, look at"
                    " the comments in odbToPointsData_03."
                    "Converter_PlasticStrainTensor.)" % matName)
                # For materials using the Abaqus key card *ELASTIC query
                # mat["EMod"] and mat["Nue"]. Then use these formulas:
                # G = 0.5 * EMod / (1+Nue)
                # K = EMod / (3 - 6*Nue)
            except IndexError:
                raise IndexError(
                    "Wrong values for parameters iPropsG or iPropsK. Material"
                    " %s has only %d properties." % len(mat['Umat']))
            # store G and K in matNbToElProps
            matNameToElProps[matName] = (G, K)
        msg("Got elastic properties for %d materials. G in [%g, %g], K in"
            " [%g, %g]" % (len(matNameToElProps),
                           min(x[0] for x in matNameToElProps.itervalues()),
                           max(x[0] for x in matNameToElProps.itervalues()),
                           min(x[1] for x in matNameToElProps.itervalues()),
                           max(x[1] for x in matNameToElProps.itervalues()),
                           ))
        return matNameToElProps

    @staticmethod
    def calc_E_plast(E, S, G, K):
        """return e_plast = E - E_elast,
        E_elast = compliance matrix * Stress"""
        a = (3.0*K+G)/(9.0*K*G)  # = 1/3G + 1/9K
        b = (G-1.5*K)/(9.0*K*G)  # = 1/9K - 1/6G
        c = 0.5/G
        return [
            E[0] - a*S[0] + b*S[1] + b*S[2],
            E[1] - b*S[0] + a*S[1] + b*S[2],
            E[2] - b*S[0] + b*S[1] + a*S[2],
            E[3] - c*S[3],
            E[4] - c*S[4],
            E[5] - c*S[5]]

    def composer1(self, cls, dat, mesh):
        "computes plastic strain tensor"

        # get some structural info from result class on very first invocation
        # and prepare class for intermediate storage and componentsDict:
        if self.fieldPos is None:
            self.fieldPos = cls.position
            self.FieldClass = createFieldClass(
                self.fieldName, self.fieldPos, "tensor")
            # self.componentsDict = {"EP_11":0, "EP_22":1, ...}
            self.componentsDict = dict(
                (x,i) for i, x in enumerate(
                    self.FieldClass.getComponentNamesList()))

        # create maskfield and matPropsField
        maskField = self.maskFieldFunc(dat, mesh)
        matPropsField = [self.matNbToElProps.get(matNb, None)
                         for matNb in dat["Material"]]

        # initialize store
        self.store = self.FieldClass(
            (mask and matGK and self.calc_E_plast(E, S, matGK[0], matGK[1])
             or self.valueNotDefined)
            for E, S, mask, matGK in izip(
                dat["E"], dat["S"], maskField, matPropsField))

        # return correct component
        i = self.componentsDict[cls.fieldName]
        return cls(vec[i] for vec in self.store)

    def composer2(self, cls, dat, mesh):
        "return previously computed strain tensor"
        # return correct component
        i = self.componentsDict[cls.fieldName]
        return cls(vec[i] for vec in self.store)

    @property
    def fieldArguments(self):
        # a list of ("EP_11", "scalar", self.composerX) - tuples
        outputFields=[
            (x, "scalar", ((i==0) and self.composer1 or self.composer2))
            for i, x in enumerate(self.FieldClass.getComponentNamesList())]

        return dict(
            fieldNames=["E", "S", self.stateSDV],
            relativeFieldNames=["E", "S"],
            matRegionsOption=True,  # creates a field Material
            matNamesList=self.matNamesList,
            outputFields=outputFields,
            )


class Converter_StressEigenDecomp(Converter):
    """
    Usage:
     >>> import os
     >>> from bae.abq_model_02 import Model
     >>> from bae.utils.odbToPointsData_03 import odbToVtkGridBox
     >>> from bae.utils.odbToPointsConverter_03 import \\
     >>>     Converter_StressEigenDecomp as Converter
     >>>
     >>> #-----
     >>> boxName = "box2_mine_h06"
     >>> gridData = ...
     >>> ...
     >>>
     >>> #-----
     >>>
     >>> converter = Converter(
     >>>     outputFields=[
     >>>         "S1", "Bearing1", "Plunge1", "S3", "Bearing3", "Plunge3"])
     >>>
     >>> odbToVtkGridBox(
     >>>     odbPath, frameList=frameList,
     >>>     gridName=boxName, gridData=gridData,
     >>>     outputDirName=outputDirName, fileNamePrefix=fileNamePrefix,
     >>>     projectPrefix=..., meshVersion=...,
     >>>     **converter.fieldArguments)
     >>>
     >>> # remove the matList.csv file
     >>> os.remove(os.path.join(outputDirName, fileNamePrefix+"_matList.csv"))


    Note:

    This requires numpy >= 1.8.
    See /mnt/boab2/it/install/PYTHON/python-numpy-upgrade4_SLES11SP3.txt
    how to install the required version on nuku SLES11-SP3 servers.
    """

    def __init__(self, outputFields, stressFieldName="S",
                 miningConvention=True):
        """
        @param stressFieldName: optionally override input data odb-field name

        @param outputFields: A list of the desired output components. Which can
           be any combination of the eigenvalues "S1", "S2", "S3", normalized
           eigenvector components "v1x", "v1y", "v1z", "v2x", "v2y", ..., "v3z",
           "Bearing1", "Plunge1", "Bearing2", "Plunge2", "Bearing3", "Plunge3",
           as well as any of the cartesian components "S_11", "S_22", "S_33",
           "S_12", "S13", "S_23". If stressFieldName is given different from
           the default "S" then the real component names of the odb field must
           be used instead of "S_11", ...

        @param miningConvention: if True, the sign and order of eigenValues
           will be inverted
        """
        self.stressFieldName = stressFieldName
        self.outputFields = outputFields
        self.miningConvention = miningConvention

    def composer1(self, cls, data, **_):
        "computes eigenvalues and vectors"

        import numpy as np
        if map(int, np.__version__.split(".")) <= [1,8]:
            raise NotImplementedError(
                "This composer requires numpy v1.8 or greater. Found %s."
                % np.__version__)

        self.store = dict()

        msg("Storing cartesian components.")
        stressField = data[self.stressFieldName]
        Svec = np.array(stressField)
        compStrToIndex = dict()
        for i, name in enumerate(stressField.getComponentNamesList()):
            self.store[name] = Svec[:,i]
            compStrToIndex[name.rsplit("_",1)[-1]] = i

        msg("Preparing stress tensors.")
        N = Svec.shape[0]
        Smat = np.zeros((N, 3,3), float)
        for j, k, compStr in [
                (0,0, "11"), (1,1, "22"), (2,2, "33"),
                (0,1, "12"), (1,0, "12"), (0,2, "13"),
                (2,0, "13"), (1,2, "23"), (2,1, "23")]:
            iVec = compStrToIndex[compStr]
            Smat[:, j, k] = Svec[:, iVec]

        msg("Computing eigenvalues and -vectors.")
        # With using eigh (for symetric matrices) the EigenVectors v will
        # be normalized and the order will be ascending
        w, v = np.linalg.eigh(Smat)

        if self.miningConvention:
            w *= -1     # invert sign for eigenValues but preserve order
        else:
            w = w[:,::-1]   # preserve sign but reversed order
            v = v[:,:,::-1]

        msg("Storing eigenvalues and vectors.")
        self.store["S1"]  = w[:,0]
        self.store["v1x"] = v[:,0,0]
        self.store["v1y"] = v[:,1,0]
        self.store["v1z"] = v[:,2,0]
        self.store["S2"]  = w[:,1]
        self.store["v2x"] = v[:,0,1]
        self.store["v2y"] = v[:,1,1]
        self.store["v2z"] = v[:,2,1]
        self.store["S3"]  = w[:,2]
        self.store["v3x"] = v[:,0,2]
        self.store["v3y"] = v[:,1,2]
        self.store["v3z"] = v[:,2,2]

        # calculate Plunge and Bearing if requested
        needsPlg = any(['plunge' in fld.lower() for fld in self.outputFields])
        needsBrg = any(['bearing' in fld.lower() for fld in self.outputFields])
        if needsPlg or needsBrg:
            msg("Calculating and storing bearing and plunge.")
            for i in range(3):

                # Bearing will be in the x-y plane...
                bearing = np.arctan2(v[:,1,i], v[:,0,i])

                # Plunge is the angle between the line and the x-y plane
                plunge = np.arcsin(-v[:,2,i])

                # Convert back to azimuths in degrees..
                plunge, bearing = np.degrees(plunge), np.degrees(bearing)
                bearing = 90 - bearing
                bearing[bearing < 0] += 360

                # If the plunge angle is upwards, get the opposite end of the
                # line
                upwards = plunge < 0
                plunge[upwards] *= -1
                bearing[upwards] -= 180
                bearing[upwards & (bearing < 0)] += 360

                # storing results
                self.store["Bearing%d" % (i+1)] = bearing
                self.store["Plunge%d" % (i+1)] = plunge

        # return correct component, convert to field_01-cls
        return cls(self.store[cls.fieldName])

    def composer2(self, cls, **_):
        "return previously computed result component, convert to field_01-cls"

        return cls(self.store[cls.fieldName])

    @property
    def fieldArguments(self):
        """This property yields the inputFieldsPerFrame and outputFields
        arguments suitable for L{odbToVtkGridBox} et al. if you want only the
        fields generated by this composer.
        See the outputFields argument of the
        L{constructor<Converter_StressEigenDecomp.__init__>}.
        """
        # A list of ("S1", "scalar", self.composerX) - tuples
        # The first item --the field name-- is taken from self.outputFields.
        outputFields=[
            (x, "scalar", ((i==0) and self.composer1 or self.composer2))
            for i, x in enumerate(self.outputFields)]

        return dict(
            inputFieldsPerFrame=[self.stressFieldName],
            outputFields=outputFields,
            )


class Converter_StandardVtk(Converter):
    """Converter to create standard state and logP values from Abaqus output.

    For cave coupling converter see L{Converter_CombinedVtk_FS4}.

    Usage:
     >>> import os.path
     >>> from bae.utils.odbToPointsData_03 import odbToVtkGridBox
     >>> from bae.utils.odbToPointsConverter_03 import \\
     >>>     Converter_StandardVtk
     >>>
     >>> # from outputVarNames.py
     >>> varNameDamage = "SDV1"
     >>> varNameState = "SDV8"
     >>>
     >>> converter = Converter_StandardVtk(
     >>>     varNameU="U_%d",
     >>>     varNameDamage=varNameDamage, varNameState=varNameState)
     >>>
     >>> odbToVtkGridBox(
     >>>     odbPath, frameList=frameList,
     >>>     gridName=boxName, gridData=gridData,
     >>>     outputDirName=outputDirName, fileNamePrefix=fileNamePrefix,
     >>>     projectPrefix=..., meshVersion=...,
     >>>     inputFieldsPerFrame=[
     >>>         converter.varNameU%1,converter.varNameU%2,converter.varNameU%3,
     >>>         "S_MIN_PRINCIPAL",
     >>>         varNameDamage,
     >>>         (varNameState, "const"),],
     >>>     outputFields=[
     >>>         converter.getU_1, converter.getU_2, converter.getU_3,
     >>>         ("logP", "scalar", converter.getLogP),
     >>>         ("S1", "scalar", converter.getS1),
     >>>         ("state", "scalar", converter.getState),
     >>>     ])
    """

    # The functions L{odbToVtkGridBox} and the like accept a parameter
    # getExtraFields. Thereby a method to read extra data e.g. from Cavesim
    # can be supplied.
    # This converter doesn't need extra fields. Passing this None-valued
    # attribute of the converter to odbToVtkGridBox makes it skip this
    # functionality.
    # This attribute is provided for convenience. So the call to
    # odbToVtkGridBox can have the same set of arguments whether
    # Converter_StandardVtk or Converter_CombinedVtk_FS4 is used.
    getExtraFields = None

    def __init__(self, varNameDamage="SDV1", varNameState="SDV8",
                 varNameU="U_%d"):
        """
        @param varNameDamage: variable name for plastic strain, e.g. "SDV1"
        @param varNameState: variable name for state, e.g. "SDV8"
        @param varNameU: "U_%d" or "UR_%d" or something like that
        """
        self.varNameState = varNameState
        self.varNameDamage = varNameDamage

        self.varNameU = varNameU
        self.getU_1 = varNameU % 1
        self.getU_2 = varNameU % 2
        self.getU_3 = varNameU % 3

    def getLogP(self, cls, data, **_):
        "provides the logP field"
        return cls((
            # 2.0 if cave (state==5) ... even if damage<=0
            2.0 if 4.5<state<5.5 else

            # 0.0 if degrading (state==2) or void (state==3) or deleted ...
            #        ... (state==0) or air (state==6) or negative damage
            0.0 if 1.5<state<3.5 or state<=0.5 or state>=5.5 or damage<=0.0 else

            # logP-formula if rock (state==1) or fill (state==4)
            log10(damage*1000+1))

            # provide the series of values...
            for damage, state in izip(
                    *map(data.get, [self.varNameDamage, self.varNameState])))

    def getS1(self, cls, data, **_):
        return cls(-x for x in data["S_MIN_PRINCIPAL"])

    def getS3(self, cls, data, **_):
        return cls(-x for x in data["S_MAX_PRINCIPAL"])

    def getState(self, cls, data, **_):
        """Convert state into:
         0 : air (deleted elements, outside the mesh)
         1 : void (Abq-state 2, 3, 6)
         2 : rock (Abq-state 1)
         3 : fill (Abq-state 4)
         4 : CAVE (Abq-state 5)
        """
        conv = {
            # note: results must be floats, otherwise Voxler creates two
            # input objects
            0: 0.0,
            1: 2.0,
            2: 1.0,
            3: 1.0,
            4: 3.0,
            5: 4.0,
            6: 1.0}
        return cls(conv.get(int(round(x)),-1.0) for x in data[self.varNameState])


class Converter_CombinedVtk_FS4(Converter_StandardVtk):
    """Converter to create "combined Abaqus/FS4 vtks".

    For Abaqus-only result converter see L{Converter_StandardVtk}.

    Usage:
     >>> import os.path
     >>> from bae.utils.odbToPointsData_03 import odbToVtkGridBox
     >>> from bae.utils.odbToPointsConverter_03 import \\
     >>>     Converter_CombinedVtk_FS4
     >>>
     >>> # from outputVarNames.py
     >>> varNameDamage = "SDV1"
     >>> varNameState = "SDV8"
     >>>
     >>> converter = Converter_CombinedVtk_FS4(
     >>>     csOutputTemplate="../../2_RUN/R02b_RST02/FS4Out/FS4out_%03d.vtk",
     >>>     stepFrameToCsStep={1: [-1,]*8 + range(1,73)},
     >>>     varNameU="U_%d",
     >>>     varNameDamage=varNameDamage, varNameState=varNameState)
     >>>
     >>> odbToVtkGridBox(
     >>>     odbPath, frameList=frameList,
     >>>     gridName=boxName, gridData=gridData,
     >>>     outputDirName=outputDirName, fileNamePrefix=fileNamePrefix,
     >>>     projectPrefix=..., meshVersion=...,
     >>>     fieldNames=[
     >>>         converter.varNameU%1,converter.varNameU%2,converter.varNameU%3,
     >>>         "S_MIN_PRINCIPAL", varNameDamage,
     >>>         (varNameState, "const")],
     >>>     getExtraFields=converter.getExtraFields,
     >>>     outputFields=[
     >>>         ("U_1", "scalar", converter.getU_1),
     >>>         ("U_2", "scalar", converter.getU_2),
     >>>         ("U_3", "scalar", converter.getU_3),
     >>>         ("logP", "scalar", converter.getLogP),
     >>>         ("S1", "scalar", converter.getS1),
     >>>         ("state", "scalar", converter.getState),
     >>>     ])

    If you need special values it's suggested to derive a subclass that
    overrides the corresponding getXXX-function. For example:
     >>> from bae.utils.odbToPointsData_03 import odbToVtkGridBox
     >>> from bae.utils.odbToPointsConverter_03 import \\
     >>>     Converter_CombinedVtk_FS4
     >>>
     >>> class MyConverter(Converter_CombinedVtk_FS4):
     >>>     def getState(self, cls, dat, mesh):
     >>>         normalState = super(MyConverter, self).getState(
     >>>             self, cls, dat, mesh)
     >>>         if not self.csData:
     >>>             return normalState
     >>>
     >>>         return cls(1.0 if ucs<0 else ns
     >>>                    for ns, ucs in izip(normalState, dat["UCS_MAG"]))
     >>>
     >>> converter = MyConverter(
     >>>     ...arguments as for Converter_CombinedVtk_FS4...)
     >>>
     >>> odbToVtkGridBox(
     >>>     odbPath, frameList=frameList, ...
     >>>     ...arguments as usual...)
    """

    def __init__(self, stepFrameToCsStep,
                 csOutputTemplate="./FlowsimOutput/FS4out_step_%03d.vtk",
                 varNameU="U_%d", varNameDamage="SDV1", varNameState="SDV8",
                 varNameUFS4="UPC_%s"):
        """
        @param csOutputTemplate: File path template including a "%03d"
            placeholder to be replaced by the Cavesim/FS4 step number. These
            vtk files must contain the fields identified by the varNameUFS4
            argument.
        @param stepFrameToCsStep: provides the Cavesim/FS4 step for a given odb
            step and frame number: stepFrameToCsStep[odbStep][odbFrame] must
            yield the right Cavesim step number --an integer >0. For frames
            without Cavesim step --i.e. not coupled steps, e.g. before
            cave coupling starts-- this should yield -1.
        @param varNameU: "U_%d" or "UR_%d" or something like that
        @param varNameDamage: variable name for plastic strain, e.g. "SDV1"
        @param varNameState: variable name for state, e.g. "SDV8"
        @param varNameUFS4: variable name for the FS4 displacements to be
            combined with the Abaqus data. Needs a "%s" component placeholder
            that will be replaced by 1, 2, 3, "MAG".
        """
        self.stepFrameToCsStep = stepFrameToCsStep
        self.csOutputTemplate = csOutputTemplate
        self.varNameU = varNameU
        self.varNameState = varNameState
        self.varNameDamage = varNameDamage

        # flag indicating if Cavesim results are there
        # Defaults to True: if getExtraFields is not being called then an error
        # occurs
        self.csData = True
        self.varNameUFS4 = varNameUFS4

    @staticmethod
    def readVTKfromFileOrZipArch(vtkPath, pauseTime=5):
        """Try to find the cavesim output. Either as plain vtk file vtkPath
        or in a zip-archive. In the latter case the vtk file is temporarily
        extracted from the zip-archive.

         - Read the vtk file from vtkPath if present.
         - Check if there is a zip file with the same name, i.e. for
           FS4out_step_064.vtk try FS4out_step_064.vtk.zip and
           FS4out_step_064.zip.
         - Assume that the folder directly holding the vtk file has been zipped,
           i.e. for vtkPath=".../CavesimOutput/VELOCITY_VTK/CAVESIM_VEL_10.vtk"
           try: in archive ".../CavesimOutput/VELOCITY_VTK.zip" members
           "CAVESIM_VEL_10.vtk" or "VELOCITY_VTK/CAVESIM_VEL_10.vtk"
         - Assume that a folder-structure has been zipped and the zip-file has
           the name of the folder plus extension ".zip", i.e.
           for vtkPath=".../OUTPUT.zip/VELOCITY_VTK/CAVESIM_VEL_10.vtk"
           try: in archive ".../OUTPUT.zip" members
           "VELOCITY_VTK/CAVESIM_VEL_10.vtk" or
           "OUTPUT/VELOCITY_VTK/CAVESIM_VEL_10.vtk"
         - Assume that a folder-structure has been zipped and the zip-file has
           the name of the folder plus a prefix like "NCA_Perse_" plus
           extension ".zip", i.e. for
           vtkPath=".../NCA_Perse_OUTPUT.zip/VELOCITY_VTK/CAVESIM_VEL_10.vtk"
           try: in archive ".../NCA_Perse_OUTPUT.zip" members
           "VELOCITY_VTK/CAVESIM_VEL_10.vtk" or
           "OUTPUT/VELOCITY_VTK/CAVESIM_VEL_10.vtk"

        @param vtkPath: vtk-file name
        @param pauseTime: seconds to wait after zip file has changed before
           trying to read again
        @returns: vtk object
        """
        failedPaths = list()

        # simply try the path as vtk file
        if os.path.isfile(vtkPath):
            vtkData = FieldsOnStructPointGrid().fromVtk(vtkPath)
            return vtkData
        else:
            failedPaths.append(vtkPath)

        # Failed to locate the file directly, now trying some heuristics.
        # This while clause will either be left with a break after defining
        # the variables archivePath and extractPathCandidates or raise an error
        while 1:

            # if there is a zip file with the given file name
            # i.e. for FS4out_step_064.vtk try FS4out_step_064.vtk.zip
            archivePath = vtkPath + ".zip"
            if os.path.isfile(archivePath):
                extractPathCandidates = [os.path.basename(vtkPath)]
                break
            else:
                failedPaths.append(archivePath)
                del archivePath

            # if there is a zip file with the given file name ex .vtk
            # i.e. for FS4out_step_064.vtk try FS4out_step_064.zip
            archivePath = os.path.splitext(vtkPath)[0] + ".zip"
            if os.path.isfile(archivePath):
                extractPathCandidates = [os.path.basename(vtkPath)]
                break
            else:
                failedPaths.append(archivePath)
                del archivePath

            # try to detect zipped folders...
            # firstly prepare list of folders in vtkPath
            fileName = os.path.basename(vtkPath)

            dirNames = []
            path = vtkPath
            while 1:
                path, folder = os.path.split(path)
                if folder != "":
                    dirNames.append(folder)
                else:
                    if path != "":
                        dirNames.append(path)
                    break
            dirNames.reverse()
            del path, folder

            # if zip file in path, i.e. ...
            # for vtkPath=".../OUTPUT.zip/VELOCITY_VTK/CAVESIM_VEL_10.vtk"
            # try archive ".../OUTPUT.zip"
            for depth in reversed(range(len(dirNames))):
                dirName = dirNames[depth]
                base, ext = os.path.splitext(dirName)
                if ext.lower()==".zip" and not dirName.startswith("."):
                    archivePath = os.path.join(*dirNames[:depth+1])
                    remainingPath = os.path.join(*dirNames[depth+1:])
                    break

            if ext.lower()==".zip":
                if os.path.isfile(archivePath):
                    extractPathCandidates = [
                        # try VELOCITY_VTK/CAVESIM_VEL_10.vtk
                        remainingPath,
                        # try OUTPUT/VELOCITY_VTK/CAVESIM_VEL_10.vtk
                        os.path.join(base, remainingPath),
                        # for archivePath=".../NCA_Perse_OUTPUT.zip"
                        # try OUTPUT/VELOCITY_VTK/CAVESIM_VEL_10.vtk
                        os.path.join(base.rsplit("_", 1)[-1], remainingPath),
                        ]
                    break
                else:
                    failedPaths.append(archivePath)

            # try folder name as zip-file path:
            # i.e. for vtkPath=".../CsOutput/VELOCITY_VTK/CAVESIM_VEL_10.vtk"
            # try archive ".../CsOutput/VELOCITY_VTK.zip"
            # or for vtkPath=".../OUTPUT/VELOCITY_VTK/CAVESIM_VEL_10.vtk"
            # try archive ".../OUTPUT/VELOCITY_VTK.zip" then ".../OUTPUT.zip"
            for depth in reversed(range(len(dirNames))):
                if dirNames[depth].startswith("."):
                    break
                archivePath = os.path.join(*dirNames[:depth+1]) + ".zip"
                if os.path.isfile(archivePath):
                    remainingPath = os.path.join(*dirNames[depth+1:])
                    extractPathCandidates = [
                        # try CAVESIM_VEL_10.vtk
                        os.path.join(remainingPath, fileName),
                        # try VELOCITY_VTK/CAVESIM_VEL_10.vtk
                        os.path.join(dirNames[depth], remainingPath, fileName),]
                    break
                else:
                    failedPaths.append(archivePath)

            # failed to locate vtk file or zip archive
            raise IOError(
                "Can't find vtkPath or corresponding zip-archive.\nTried: %s"
                % (', '.join(failedPaths)))

        # process the archive we just found
        archive = zipfile.ZipFile(archivePath, 'r')
        stored = archive.namelist()

        extractPath = None
        for extractPath in extractPathCandidates:
            if extractPath in stored:
                break
            else:
                extractPath = None
        if not extractPath:
            raise IOError(
                "Can't identify the vtk %s in zip-archive %s. Tried: %s"
                % (os.path.basename(vtkPath), archivePath,
                   ", ".join(extractPathCandidates)))

        msg("Reading %s from zip archive %s" % (extractPath, archivePath),
            debugLevel=1)

        tmpVtkfile = None
        for cnt in range(20):  # try 20 times to extract from zip file

            checker = CheckFileFinished(archivePath)
            try:
                tmpVtkfile = archive.extract(extractPath)
            except KeyError as exc:
                raise KeyError("Path %s vanished from zip archive %s.\n%s"
                               % (extractPath, archivePath, exc.args[0]))
            except zipfile.BadZipfile:
                if checker.finished():
                    raise
                msg("WARNING: The zip archive %s seems to have changed"
                    " while trying to extract %s. Trying again after"
                    " waiting for %s sec."
                    % (archivePath, extractPath, pauseTime))
                time.sleep(pauseTime)
            else:
                break

        # finished extracting
        if tmpVtkfile is None:
            raise IOError(
                "Can't find the vtk %s in zip-archive %s after %d attempts"
                % (os.path.basename(vtkPath), archivePath), cnt)

        # getting the data and cleaning up
        vtkData = FieldsOnStructPointGrid().fromVtk(tmpVtkfile)
        os.remove(tmpVtkfile)

        return vtkData

    def getExtraFields(self, ctrl, odbStep, odbFrame, topo, **_):
        """Method to read extra data from FS4 output vtk files.
        This method must be supplied to the parameter getExtraFields of
        functions L{odbToVtkGridBox} and the like.
        """

        try:
            csStep = self.stepFrameToCsStep[odbStep][odbFrame]
        except (KeyError, IndexError):
            msg("ERROR: odb steps provided in stepFrameToCsStep: %s"
                "\nRequested step %s, frame %d."
                % (sorted(self.stepFrameToCsStep), odbStep, odbFrame))
            raise
        if csStep<0:
            # no Cavesim/FS4 output for current step/frame
            self.csData = False
            return

        # read the Cavesim/FS4 vtks
        fname = self.csOutputTemplate % csStep
        msg("Reading %s" % fname)
        vtkFlow = self.readVTKfromFileOrZipArch(fname)

        # Cavesim data available
        self.csData = True

        # interpolate fields and "yield"
        # density (0: air, -1: intact rock, >0 : number of particles in a cell)
        # for fieldName in ["du1", "du2", "du3", "dumag", "density"]:
        for x in [1, 2, 3, "MAG"]:
            csFieldName = self.varNameUFS4 % x
            fieldName = "UCS_%s" % x
            msg("Interpolating field %s from FS4-field %s"
                % (fieldName, csFieldName))

            # convert displacement components: x>9998 => x=0
            try:
                fld = vtkFlow.data[csFieldName]
            except KeyError:
                # ignore component if not available
                continue

            # interpolate to output points
            fld = createFieldObject(
                fieldName, "point", "scalar", initArgs=[
                    fld.interpolateToPoints(vtkFlow.mesh, topo),])
            ctrl.storeField(fld)

        # end
        return

    def getU_1(self, cls, data, **_):
        "provides the combined U_1 field"
        if not self.csData:
            return cls(data[self.varNameU % 1])

        return cls((
            # if not in coupling-box: take Abq-values
            # if 0.5 < state_Abq < 4.5: take Abq-values
            u_Abq if (u_CS is None) or (0.5<state_Abq<4.5)

            # if 4.5<abq_state (abq-cave,-air) or 0.5<abq_state (abq-deleted):
            # take FS4-values
            else u_CS if (state_Abq>=4.5 or state_Abq<=0.5)

            # else: return 0.0
            else 0.0)

            # provide the series of values...
            for u_CS, u_Abq, state_Abq in izip(
                    data["UCS_1"],
                    data[self.varNameU%1], data[self.varNameState]))

    def getU_2(self, cls, data, **_):
        "provides the combined U_2 field"
        if not self.csData:
            return cls(data[self.varNameU % 2])

        return cls((
            # if not in coupling-box: take Abq-values
            # if 0.5 < state_Abq < 4.5: take Abq-values
            u_Abq if (u_CS is None) or (0.5<state_Abq<4.5)

            # if 4.5<abq_state (abq-cave,-air) or 0.5<abq_state (abq-deleted):
            # take FS4-values
            else u_CS if (state_Abq>=4.5 or state_Abq<=0.5)

            # else: return 0.0
            else 0.0)

            # provide the series of values...
            for u_CS, u_Abq, state_Abq in izip(
                    data["UCS_2"],
                    data[self.varNameU%2], data[self.varNameState]))

    def getU_3(self, cls, data, **_):
        "provides the combined U_3 field"
        if not self.csData:
            return cls(data[self.varNameU % 3])

        return cls((
            # if not in coupling-box: take Abq-values
            # if 0.5 < state_Abq < 4.5: take Abq-values
            u_Abq if (u_CS is None) or (0.5<state_Abq<4.5)

            # if 4.5<abq_state (abq-cave,-air) or 0.5<abq_state (abq-deleted):
            # take FS4-values
            else u_CS if (state_Abq>=4.5 or state_Abq<=0.5)

            # else: return 0.0
            else 0.0)

            # provide the series of values...
            for u_CS, u_Abq, state_Abq in izip(
                    data["UCS_3"],
                    data[self.varNameU%3], data[self.varNameState]))

    def getState(self, cls, data, **_):
        """
        Creates those values::
          0 : AIR - In Abaqus: outside model or deleted element;
                    in Cavesim / FS4: no particles
          1 : void
          2 : rock
          3 : fill
          4 : CAVE
        """
        abqOnlyStateFld = Converter_StandardVtk.getState(self, cls, data, **_)
        if not self.csData:
            return abqOnlyStateFld

        return cls((
            # Correct the state==0 if there is movement in the FS4 data
            4.0 if (state_Abq==0.0) and (u_CS>0.0)

            # ... else take state from Abaqus
            else state_Abq)

            # provide the series of values...
            for u_CS, state_Abq in izip(
                    data["UCS_MAG"], abqOnlyStateFld))


class Converter_CombinedVtk_FS4_multiBox(Converter_CombinedVtk_FS4):
    """Converter to create "combined Abaqus/FS4 vtks" from multiple coupling
    boxes.

    Works exactly like L{Converter_CombinedVtk_FS4} (the base class) except for
    different arguments to the constructor.
    """
    def __init__(
            self, caveCouplingVtkNb, stepFrameToCsStep, csOutputTemplate,
            varNameDamage, varNameState, varNameU="U_%d",
            varNameUFS4="UPC_%s"):
        """
        @param caveCouplingVtkNb: number of coupling boxes
        @param csOutputTemplate: A list of file path templates including a
            "%03d" placeholder to be replaced by the Cavesim/FS4 step number.
            One for each coupling box. These vtk files must contain the fields
            identified by the varNameUFS4 argument.
        @param stepFrameToCsStep: A list of relations: Abaqus odb-frame to FS4
            mining step number (==Abaqus Transfer Number --ATN). One item
            for each coupling box. Each provides the Cavesim/FS4 step for a
            given odb step and frame number:
            stepFrameToCsStep[iBox][odbStep][odbFrame] must yield the right
            Cavesim step number --an integer >0 for the coupling box with
            index iBox. For frames without Cavesim step --i.e. not coupled
            steps, e.g. before cave coupling starts-- this should yield -1.
        @param varNameU: "U_%d" or "UR_%d" or something like that
        @param varNameDamage: variable name for plastic strain, e.g. "SDV1"
        @param varNameState: variable name for state, e.g. "SDV8"
        @param varNameUFS4: variable name for the FS4 displacements to be
            combined with the Abaqus data. Needs a "%s" component placeholder
            that will be replaced by 1, 2, 3, "MAG".
        """

        # check special arguments
        if (not isinstance(stepFrameToCsStep, (list, tuple))
            or len(stepFrameToCsStep) != caveCouplingVtkNb):
            raise ValueError(
                "ERROR: Argument stepFrameToCsStep must be a list with %d"
                " dictionaries, corresponding to the given caveCouplingVtkNb."
                % (caveCouplingVtkNb, caveCouplingVtkNb))
        if (not isinstance(csOutputTemplate, (list, tuple))
            or len(csOutputTemplate) != caveCouplingVtkNb):
            raise ValueError(
                "ERROR: Argument csOutputTemplate must be a list with %d"
                " dictionaries, corresponding to the given caveCouplingVtkNb."
                % (caveCouplingVtkNb, caveCouplingVtkNb))

        self.caveCouplingVtkNb = caveCouplingVtkNb
        self.stepFrameToCsStep = stepFrameToCsStep
        self.csOutputTemplate = csOutputTemplate
        self.varNameU = varNameU
        self.varNameState = varNameState
        self.varNameDamage = varNameDamage
        self.varNameUFS4 = varNameUFS4

        # flag indicating if Cavesim results are there
        # Defaults to True: if getExtraFields is not being called then an error
        # occurs
        self.csData = True

    def getExtraFields(self, ctrl, odbStep, odbFrame, topo, **_):
        """Method to read extra data from FS4 output vtk files.
        This method must be supplied to the parameter getExtraFields of
        functions L{odbToVtkGridBox} and the like.
        """

        vtkDataList = list()  # list of vtk data potentially for multiple boxes
        self.csData = False
        fileNames = list()  # for diagnostic output
        for iBox in range(self.caveCouplingVtkNb):

            try:
                csStep = self.stepFrameToCsStep[iBox][odbStep][odbFrame]
            except (KeyError, IndexError):
                msg("ERROR: odb steps provided in stepFrameToCsStep[%d]: %s"
                    "\nRequested step %s, frame %d."
                    % (iBox, sorted(self.stepFrameToCsStep), odbStep, odbFrame))
                raise

            if csStep < 0:
                # Cavesim/FS4 output not available for current step/frame
                # for this box
                continue

            self.csData = True

            # read the Cavesim/FS4 vtks
            fname = self.csOutputTemplate[iBox] % csStep
            msg("Reading %s" % fname)
            vtkDataList.append(self.readVTKfromFileOrZipArch(fname))
            fileNames.append(fname)

        # stop here if no coupling box has data for this mining step
        if not self.csData:
            return

        # interpolate fields and "yield"
        # density (0: air, -1: intact rock, >0 : number of particles in a cell)
        # for fieldName in ["du1", "du2", "du3", "dumag", "density"]:
        for x in [1, 2, 3, "MAG"]:
            csFieldName = self.varNameUFS4 % x
            fieldName = "UCS_%s" % x
            msg("Interpolating field %s from FS4-field %s"
                % (fieldName, csFieldName))

            # iterate over coupling boxes
            for iBox, vtkData in enumerate(vtkDataList):

                if csFieldName not in vtkData.data:
                    msg("WARNING: field <%s> not available in %s. (Maybe we"
                        " don't need it?)" % (csFieldName, fileNames[0]))
                    # ignore component if not available
                    continue

                newFld = vtkData.data[csFieldName]

                # interpolate to output points
                if iBox==0:
                    fld = newFld.interpolateToPoints(vtkData.mesh, topo)
                else:
                    # combine: overwrite if possible
                    fld = ((x if y is None else y) for x, y in izip(
                        fld, newFld.interpolateToPoints(vtkData.mesh, topo)))

            # finalize: store in correct field object
            fld = createFieldObject(
                fieldName, "point", "scalar", initArgs=[fld,])
            ctrl.storeField(fld)
            return


class Converter_CombinedOutputOnSurface(Converter_CombinedVtk_FS4):
    """Variant of L{Converter_CombinedVtk_FS4} suitable for generating output
    at the top surface, like surface subsidence.

    Points sitting in a cell that is air (CellParticleState==0) will be
    supplied with the value of the first non-air cell below.

    All other formulas and processes are the same as for
    L{Converter_CombinedVtk_FS4}.

    Usage:
     >>> import os.path
     >>> from bae.utils.odbToPointsData_03 import \\
     >>>     odbVariableFaceCentroidsToCsvPoints, \\
     >>> from bae.utils.odbToPointsConverter_03 import \\
     >>>     Converter_CombinedOutputOnSurface
     >>>
     >>> # from outputVarNames.py
     >>> varNameDamage = "SDV1"
     >>> varNameState = "SDV8"
     >>>
     >>> converter = Converter_CombinedOutputOnSurface(
     >>>     stepFrameToCsStep={"Step-3": [-1,]*8 + range(1,73)},,
     >>>     csOutputDir="../../2_RUN/R02b_RST02/CavesimOutput",
     >>>     varNameU="U_%d",
     >>>     varNameDamage=varNameDamage, varNameState=varNameState)
     >>>
     >>> odbVariableFaceCentroidsToCsvPoints(
     >>>     odbPath, frameList=frameList,
     >>>     initialSurface="top", seqSetsNameExcav="seqElsetsPit",
     >>>     outputDirName=outputDirName, fileNamePrefix=fileNamePrefix,
     >>>     inputFieldsPerFrame=["U_1", "U_2", "U_3", varNameDamage,
     >>>                               (varNameState, "const")],
     >>>     getExtraFields=converter.getExtraFields,
     >>>     outputFields=[
     >>>         ("U_1", "scalar", converter.getU_1),
     >>>         ("U_2", "scalar", converter.getU_2),
     >>>         ("U_3", "scalar", converter.getU_3),
     >>>         ("logP", "scalar", converter.getLogP),
     >>>     ])
    """

    @staticmethod
    def interpolateOnSurface(oldField, filterField, oldMesh, points):
        """Interpolate values from Cavesim (oldData on oldMesh) onto given
        points in a piecewise constant manner. If the corresponding grid point
        (in Cavesim output) doesn't have valid data then take that of the next
        point vertically below that has proper data. This is a generator
        function suitable as initializer for a list/field-object

        @param oldField: Field(pos="structPt") object holding the data to be
           interpolated
        @param filterField: Field is True in cells that contains valid data.
           Otherwise try the cell below...
        @param oldMesh: MeshStructuredPoints object for oldField (and
           filterField)
        @param points: iterable containing the point coordinates where
           interpolated values are seeked for.
        """

        stride = oldMesh.strides[2]  # abbreviation

        # list of point coordinates/L{MeshUnstructuredPoints}
        ptIdsWeights = oldMesh.getPtIdsWeights(points)
        for ptIds, weights in ptIdsWeights:
            if ptIds is None:
                yield None
            else:

                # index of the largest element coordinate / closest corner
                # Note: in pure python this method is faster than:
                # max(enumerate(coords), key=lambda x:x[1])[0]
                # max(ids, key=coords.__getitem__)
                # max((c,i) for i,c in enumerate(coords))[1]
                # max((c,i) for c,i in izip(coords, ids))[1]
                idx = weights.index(max(weights))

                # if not filterField go down one cell
                while not filterField[idx] and idx>=stride:
                    idx -= stride
                yield oldField[idx]

    def getExtraFields(self, ctrl, odbStep, odbFrame, topo, **_):
        """Method to read extra data from Cavesim output vtk files.
        This method must be supplied to the parameter getExtraFields of
        functions L{odbToVtkGridBox} and the like.
        """
        # part of the code has been out-sourced to a common hidden method
        vtkVel, vtkFlow, dispFields = (
            ##### CONSTRUCTION SITE AHEAD ####
            # check out old odbToPointsData_02.Converter_CombinedVtk_Old._getExtraFields_Prepare !!!!
            Converter_CombinedVtk._getExtraFields_Prepare(
                self, odbStep, odbFrame, topo))
        if vtkVel is None:
            return
        for fld in dispFields:
            #### .... really all those fields???
            ctrl.storeField(fld)

        # prepare filter field: True if not air
        msg("Preparing filter field.")
        fld = vtkFlow.data["CellParticleState"]  # air: 0
        airFilter = type(fld)(x>0.5 for x in fld)

        # there is only 0,1,2 values, air: <1, cave: >=1
        fieldName = "CellParticleState"
        msg("Interpolating field %s" % fieldName)
        fld = vtkFlow.data[fieldName]
        fld = createFieldObject(fieldName, "point", "scalar", initArgs=[
            self.interpolateOnSurface(fld, airFilter, vtkFlow.mesh, topo),])
        ctrl.storeField(fld)
        return


class Converter_StandardAndFOS(Converter_StandardVtk):
    """Variant of L{Converter_StandardVtk} suitable for generating output
    for local factor of safety and damage variable (unified pst dependent
    strength ratio).

    The different types of FOS will be calculated with respect to local
    material type, damage (pst), (fill-)state and anisotropic orientation.

    Special settings/arguments for FOS calculation can be set via
    param-dictionarys (ivars, see)

    Usage:
        >>> from bae.utils.odbToPointsData_03 import odbToVtkGridBox
        >>> from bae.utils.odbToPointsConverter_03 import \\
        >>>     Converter_StandardAndFOS
        >>> from bae.misc_01 import Config
        >>> config = Config("configBE.py")
        >>> converter = Converter_StandardAndFOS(
        ...                 isotropic=False,
        ...                 varNameCohsNormals=config.varNameCohsNormals,
        ...                 varNameDamage=config.vumatVerDamage,
        ...                 varNameState=config.vumatVerSTATE,
        ...                 varNameTimeFill=config.varNameTimeFill
        ...                 )
        >>> # something special for FOSMohrCoulomb here:
        >>> #     - use a fixed frictionAngle=35deg
        >>> #     - force cohesion to pass UCS
        >>> #     - compare against intact rock (>> pst = 0)
        >>> converter.paramMC['phi'] = 35.0
        >>> converter.paramMC['constrainCohesionToUCS'] = True
        >>> converter.paramMC['pst'] = 0.0
        >>> outputFields = [
        ...     ("PST", "scalar", config.vumatVerDamage),
        ...     ("STATE", "scalar", converter.getState),
        ...     ("DAMAGEVAR","scalar", converter.damageVariable),
        ...     ("MATERIAL", "scalar", "Material"),
        ...     ("FOS_MC_IntactCohs_phi35", "scalar",converter.fosMohrCoulomb),
        ...     ("FOS_MISES", "scalar", converter.fosMises),
        ...     ("FOS_P", "scalar", converter.fosHydrostatic),
        ...     ]
        >>> boxName = 'myFOSBox'
        >>> fileNamePrefix = "_".join((config.projectPrefix, config.runVersion,
        ...                            config.seqVersion, boxName))
        >>> odbToVtkGridBox(
        ...     odbPath=config.odbPath,
        ...     frameList=config.myFrameList,
        ...     inputFieldsPerFrame=converter.inputFields,
        ...     outputFields=outputFields,
        ...     matRegionsOption=matRegionsOption,
        ...     getExtraFields=converter.getExtraFields,
        ...     gridName=boxName,
        ...     gridData=config.gridDataDict[boxName],
        ...     outputDirName='./VTK',
        ...     fileNamePrefix=fileNamePrefix,
        ...     projectPrefix=config.projectPrefix,
        ...     meshVersion=config.meshVersion,
        ...     abaqusVersion="abq6132",
        ...     outputFormat="binary",
        ...     interpolationDataDir='./',
        ...    )
    """

    #: Parameters for MohrCoulomb approximation of FOS. See
    #: L{bae.material_01.material.FactorOfSafetyIso.fosMohrCoulomb<fosMohrCoulomb>}
    paramMC = dict(
        pst=None,  # if not None: sets a fixed value for pst-evaluation
        c=None,    # if not None: sets a fixed cohesion
        phi=None,  # if not None: sets a fixed friction angle

        # if True and c is None: forces fit to match UCS
        constrainCohesionToUCS=False,
        fittingRange=[0,None],        # restricts fitting range of p
        measureFromBisection=True,    # measures from S1=S3 instead of S1=0
        )

    #: Parameters for HoekBrown approximation of FOS. See
    #: L{bae.material_01.material.FactorOfSafetyIso.fosHoekBrown<fosHoekBrown>}
    paramHB = dict(
        pst=None,  # if not None: sets a fixed value for pst-evaluation
        measureFromBisection=True,   # measures from S1=S3 instead of S1=0
        )

    #: Parameters for Mises approximation of FOS.
    paramMises = dict(
        pst=None,  # if not None: sets a fixed value for pst-evaluation
        )

    #: Parameters for Hydrostatic pressure approximation of FOS.
    paramHydrostatic = dict(
        pst=None,  # if not None: sets a fixed value for pst-evaluation
        )

    #: infinite fos values will set to this value
    infValue = 10

    #: UCSi threshold value; FOS skipped for material assumed to be elastic
    UCSiElastic = 500.E6

    #: offset value for (var)FillTime
    varFillOffset = 1000

    def __init__(self, matInputFiles, **kwargs):
        """
        @kwarg varNameDamage: variable name for plastic strain, e.g. "SDV1"
        @kwarg varNameState: variable name for state, e.g. "SDV8"
        @kwarg varNameU: "U_%d" or "UR_%d" or something like that
        @kwarg fixedFillValue: default is None. If not None, filled cells will
            get this fos-value. No fos-calculation to fill/varFill will be
            applied.
        @kwarg initSDVsVTKPathDummy: default is None. For not-restart jobs the
            fill time and thus the varFill(offset) usually gets stored in the
            first (equilibrium) step and initSDVsVTKPathDummy is not required.
            For restarts you'll have to run odbToVtkGridBox_getFillTime.py first
            and set initSDVsVTKPathDummy accordingly.
        @kwarg varNameTimeFill: variable name for fillTime. Not required if
            initSDVsVTKPathDummy is not None or fixedFillValue is set.
        @kwarg vusubsParamPyPath: default is None. If not specified, the
            vusubs_param.py is assumed to be at odbPath/source/vusubs_param.py.
        @kwarg isotropic: default is True
        @kwarg varNameCohsNormals: variable name for (first) anistropic
            angle(s). Can be None/unset if isotropic is True or
            initSDVsVTKPathDummy is set.
        """
        ## init some attributes
        self.rockLookup = None
        self.fillLookup = None
        self.matField = None
        self.fillField = None
        self.anisoAlpha = None
        self.anisoBeta = None
        self.odbReader = None

        self._currentStepFrame = (None, None)
        self._lastStepFrame = (None, None)
        self._lastS1, self._lastS2, self._lastS3 = None, None, None
        self._lastPst, self._lastState = None, None
        self._lastMesh = None

        ## parse kwargs
        # we will pop the variables from kwargs, then pass the rest to
        # StandardConverter.__init__ later
        self.isotropic = kwargs.pop('isotropic', True)
        self.fixedFillValue = kwargs.pop('fixedFillValue', None)
        self.varNameTimeFill = kwargs.pop('varNameTimeFill', None)
        self.varNameCohsNormals = kwargs.pop('varNameCohsNormals', None)
        self.vusubsParamPyPath = kwargs.pop('vusubsParamPyPath', None)
        self.initSDVsVTKPathDummy = kwargs.pop('initSDVsVTKPathDummy', None)

        ## check inputs
        if self.initSDVsVTKPathDummy:
            raise NotImplementedError(
                'Getting FillTime and anisotropic angles from'
                ' existing vtk is not implemented/tested yet')

        if not self.isotropic and self.varNameCohsNormals is None:
            raise ValueError('varNameCohsNormals (e.g. "SDV15") required for'
                             ' anisotropic fos-calculation')

        if self.fixedFillValue is None and self.varNameTimeFill is None:
            raise ValueError(
                'Set a fixedFillValue to skip fos calculation for fill'
                ' or define the varNameTimeFill (e.g. "SDV13")')

        super(Converter_StandardAndFOS, self).__init__(**kwargs)


    @property
    def inputFields(self):
        """List of input fields required to calculate FOS-data.
        Isotropic materials (or those with anisoS, anisoN = 1, 1) require only
        Eigenstresses. Anisotropic materials require the complete stress
        tensor.

        Note that stress values should be 'interpolated' by using the piecwise
        constant method (takes nearest integration point) to avoid interpolation
        artefacts like negative FOS values.
        """
        if self.isotropic:
            fieldNames = [
                self.varNameDamage, self.varNameState,
                "S_MIN_PRINCIPAL", "S_MID_PRINCIPAL", "S_MAX_PRINCIPAL"]
        else:
            fieldNames = [
                self.varNameDamage, self.varNameState,
                "S_11", "S_22", "S_33", "S_12", "S_23", "S_13"]

        fieldNames = [(f, 'const') if f.startswith('S_') else f
                      for f in fieldNames]
        return fieldNames


    @property
    def standardFOSOutputFields(self):
        return [("REL_STRENGTHFRACTURE", "scalar", self.getDamageVariable),
                ("FOS_vonMISES", "scalar", self.getFosMises),
                ("FOS_HYDROSTATIC", "scalar", self.fosHydrostatic),]


    @property
    def fosType(self):
        if self.isotropic:
            return FactorOfSafetyIso
        else:
            return FactorOfSafetyAniso


    def getExtraFields(self, odbStep, odbFrame, points, odb, **dummy):
        """Enables access to the odbReader which is required to set up the
        material lookup (see L{createLookups}).
        The odbReader is stored in class attribute 'odbReader'.
        """
        self.odbReader = odb
        self._currentStepFrame = (odbStep, odbFrame)
        return
        yield


    def createLookups(self, matNumberList, matField, matProps):
        """Creates the
            - matCode to material lookup dict (self.rockLookup),
            - matCode field on grid (self.matField)
            - *varFillCode to fillMaterial lookup dict (self.fillLookup)
            - *varFill field on grid (self.fillField)
            - *anisotropy-normal-field (angles, self.anisoNormalsField)
        Data tagged with * will only be created if required.

        Usually, this function is just called once during the first
        call of L{getFosInput}.

        The rock and fill materials will be stored as
        L{PST-samples<bae.material_01.material.FactorOfSafetyIso>}
        in lookup dicts.

        @param matNumberList: list of (mat number, mat name)-tuples.
           The composer function gets it as ctrl.matNumberList.
        @param matField: L{Field<bae.field_01.Field>} object holding the
           material numbers / mat-codes. The composer function gets it as
           data["Material"].
        @param matProps: The composer function gets it as postData.material.matProps.
        """
        # init lookup-dicts holding matKey: material
        self.rockLookup, self.fillLookup = dict(), dict()

        # init matCode-field
        self.matField = np.round(matField).astype(int)

        # reading matprops from abaqus input
        mpc = MatPropsCollection.fromDict(matProps)
        rock = mpc.rock
        fill = mpc.fill  # 'default'-fill if varFill-offset is 0
        for rockName in rock:
            if rockName.upper() == rockName:
                continue
            rock.renameItem(rockName, rockName.upper())
            fill.renameItem(rockName, rockName.upper())

        # create rockMaterial lookup dict
        for ii, name in enumerate(matNumberList, start=1):
            msg('    adding %s to rockLookup' % name)
            fos = self.fosType(rock[name])
            setattr(fos, 'name', name)
            self.rockLookup[ii] = fos

            #### Tobias needs to explain this .....
            # assign defaultFills to fillLookup
            fos = self.fosType(fill[name])
            setattr(fos, 'name', name)
            self.fillLookup[ii] = fos

        # create (var)FillMaterial lookup dict
        if self.fixedFillValue is not None:
            # fill lookups not required
            return

        if self.vusubsParamPyPath is None:
            # assume path is odbPath + source/vusubs_param.py
            vusubs = os.path.join(os.path.dirname(self.odbReader.odbPath),
                                  'source', 'vusubs_param.py')
        else:
            vusubs = self.vusubsParamPyPath

        #note that varfill stored in mpc has the same order as found in
        #vusubs_param.py
        for ii, (name, varFillMat) in enumerate(mpc.fill.iteritems(), start=1):
            msg('    adding varFill %s to fillLookup' % name)
            fos = self.fosType(varFillMat)
            setattr(fos, 'name', name)
            self.fillLookup[ii*self.varFillOffset] = fos

        # derive varFillNumber from TimeFill-SDV
        if self.initSDVsVTKPathDummy:
            # getting TimeFill-SDV from existing VTK (required for restart)
            from bae.field_02 import FieldsCollection
            vtkPath = self.initSDVsVTKPathDummy % (self.gridName)
            flds = FieldsCollection.fromVtk(vtkPath)
            fillTime = flds['TimeFill']
        else:
            # getting TimeFill-SDV from equilibrium step
            msg('getting TimeFill from odb')
            step, frame = 1, 0
            fillTime = self.odbReader.getFieldFromOdb(
                step, frame, self.varNameTimeFill, interpolationType='const')
            fillTime = np.array(fillTime)

        self.fillField = self.matField.copy()   #first: default fill everywhere
        mask = (fillTime >= self.varFillOffset) #second: overwrite with varfill
        varFills = (fillTime[mask] // self.varFillOffset)*self.varFillOffset
        self.fillField[mask] = varFills.astype(int)

        if self.isotropic:
            return

        # get normal angles for anisotropy
        if self.initSDVsVTKPathDummy:
            # getting TimeFill-SDV from existing VTK (required for restart)
            from bae.field_02 import FieldsCollection
            vtkPath = self.initSDVsVTKPathDummy % (self.gridName)
            flds = FieldsCollection.fromVtk(vtkPath)
            self.anisoAlpha = flds['ANISO_ALPHA']
            self.anisoBeta = flds['ANISO_BETA']
        else:
            # getting TimeFill-SDV from equilibrium step
            step, frame = 1, 0
            alphaSDV = self.varNameCohsNormals
            try:
                betaSDV = 'SDV%d' % (int(alphaSDV.replace('SDV','')) + 1)
            except Exception as e:
                e.message += (' Expeceted something like "SDV$$" for'
                              ' varNameCohsNormals, but is %s'
                              % alphaSDV)
                msg('ERROR: %s' % e.message)
                raise e

            msg('getting anisotropic directions/angles from odb')
            self.anisoAlpha = np.array(self.odbReader.getFieldFromOdb(
                step, frame, alphaSDV, interpolationType='const'))
            self.anisoBeta = np.array(self.odbReader.getFieldFromOdb(
                step, frame, betaSDV, interpolationType='const'))


    def getFosInput(self, data, topo, ctrl):
        """Returns numpy arrays of principle stresses, pst, state and matcode.
        Make sure all fields are specified as fieldNames for
        odbToPointsData_03.
        """
        if self._lastMesh is not topo:
            # mesh has changed --> new lookup required
            self.createLookups(
                ctrl.matNumberList,
                data["Material"],
                postData.material.matProps)

        if (self._lastMesh is topo
                and self._lastStepFrame == self._currentStepFrame):
            # to avoid overhead if getFosInput gets called multiple times
            # in same step.
            return (
                self._lastS1.copy(), self._lastS2.copy(), self._lastS3.copy(),
                self._lastPst.copy(), self._lastState.copy())

        pst = np.array(data[self.varNameDamage]).astype(float)
        state = np.round(data[self.varNameState]).astype(int)

        if self.isotropic:
            S3 = np.array(data["S_MIN_PRINCIPAL"]).astype(float)
            S2 = np.array(data["S_MID_PRINCIPAL"]).astype(float)
            S1 = np.array(data["S_MAX_PRINCIPAL"]).astype(float)
        else:
            from bae.material_01.yieldsurfaces import rotVoigtTransversal
            # setup full stress tensor and apply the anisotropic
            # distorsion for rock material

            msg("ATTENTION! Anisotropy is only applied to rock material")
            S = np.empty((pst.shape[0], 3, 3))

            isRock = (state==1)

            # voxels that are not in rock material
            # only assing lower triangle of S as required for eigenvalsh
            S[:,0,0] = np.array(data["S_11"])
            S[:,1,1] = np.array(data["S_22"])
            S[:,2,2] = np.array(data["S_33"])
            S[:,1,0] = np.array(data["S_12"])
            S[:,2,1] = np.array(data["S_23"])
            S[:,2,0] = np.array(data["S_13"])

            # loop rock materials and apply anisotropic distorsion of S
            for matIdx in np.unique(self.matField[isRock]):
                material = self.rockLookup[matIdx]
                anisoS, anisoN = material.anisoS, material.anisoN
                ii = isRock & (self.matField == matIdx)
                if not ii.any():
                    continue

                msg(('Transforming StressTensor for %d sampling points'
                     ' with material %s') %(ii.sum(), material.name))
                try:
                    # get normal-vector angles
                    a,b = self.anisoAlpha[ii], self.anisoBeta[ii]

                    # apply rotation to stressTensor
                    V = [S[ii,0,0], S[ii,1,1], S[ii,2,2],
                         S[ii,1,0], S[ii,2,1], S[ii,2,0]]

                    v11, v22, v33, v12, v23, v31 = rotVoigtTransversal(V, a, b)

                    # scale rotated stess with anisotropic parameters
                    # only assing lower triangle of S as required for eigenvalsh
                    S[ii,0,0] = anisoN * v11
                    S[ii,1,1] = anisoN * v22
                    S[ii,2,2] = 1 * v33
                    S[ii,1,0] = anisoN * v12
                    S[ii,2,1] = anisoS * v23
                    S[ii,2,0] = anisoS * v31

                except Exception as e:
                    msg('ERROR: %s, %s' % (type(e), e.message))
                    raise

            try:
                # eigvalsh calculates eigenVals from Lower triangle (UPLO='L')
                # eigenValues already sorted sorted
                msg("getting eigen values of anisotropic deformed stress"
                    " tensor at %d sampling points" % len(isRock))
                E = np.linalg.eigvalsh(S, UPLO='L')
                S3, S2, S1 = E.T  # E.T. home phone

            except Exception as e:
                msg('ERROR during eigenvaluecalculation: %s, %s'
                    % (type(e), e.message))
                raise

            self._lastStepFrame = self._currentStepFrame
            self._lastS1, self._lastS2, self._lastS3 = S1, S2, S3
            self._lastPst, self._lastState = pst, state
            self._lastMesh = topo

        return S1, S2, S3, pst, state

    def _fosMain(self, data, topo, ctrl, pst, fosFun, miningNotation):
        """FOS-'workflow' is the same for all implemented variants of
        FOS. Just the actual calculation, done by fosFun, differs.

        @param data: odb-fields as provided for each converter function
        @param pst: plastic strain. Can be a fixed value or a (odb)field
        @param fosFun: wrapped fos-function that takes
            material, S1, S2, S3, pst
        @param miningNotation: forces reversed sorted negative S1, S2, S3
            -values to be passed to fosFun
        """

        S1, S2, S3, pStrain, state = self.getFosInput(data, topo, ctrl)

        if pst is None:
            pst = pStrain  # evaluate pst dependent yieldSurfaces

        if miningNotation:
            S1, S2, S3 = -S3, -S2, -S1

        fos = np.zeros_like(S1)
        isRock = (state==1)
        isFill = (state==4)
        fos[isRock | isFill] = np.nan  # make sure that we got everything

        ## fos for rock material
        for matIdx in np.unique(self.matField[isRock]):
            material = self.rockLookup[matIdx]
            ii = isRock & (self.matField == matIdx)

            if not ii.any():
                continue

            if material.UCSi < self.UCSiElastic:
                msg('    - calculating fos for %d positions in rock %s'
                    % (ii.sum(), material.name))
                try:                # to get pstField
                    _pst = pst[ii]  # --> pst is a field
                except TypeError:   # fails for a scalar float value
                    _pst = pst      # --> pst is a fixed value

                matFos = fosFun(material, S1[ii], S2[ii], S3[ii], _pst)
            else:  # is elastic rock --> set fos = infValue
                msg('    - rock %s is elastic (UCSi>=%.1f)'
                    % (material.name, self.UCSiElastic))
                matFos = self.infValue

            fos[ii] = matFos

        ## fos for fill material
        if self.fixedFillValue is not None:
            # sets a fixed value to sampling points in fill
            fos[isFill] = self.fixedFillValue
            return fos

        for fillIdx in np.unique(self.fillField[isFill]):
            material = self.fillLookup[fillIdx]
            ii = isFill & (self.fillField == fillIdx)
            if not ii.any():
                continue

            if material.UCSi < self.UCSiElastic:
                msg('    - calculating fos for %d positions in fill %s'
                    % (ii.sum(), material.name))
                try:                # to get pstField
                    _pst = pst[ii]  # --> pst is a field
                except TypeError:   # fails for a scalar float value
                    _pst = pst      # --> pst is a fixed value

                matFos = fosFun(material, S1[ii], S2[ii], S3[ii], _pst)
            else:  # is elastic fill --> set fos = infValue
                msg('    - fill %s is elastic (UCSi>=%.1f)'
                    % (material.name, self.UCSiElastic))
                matFos = self.infValue

            fos[ii] = matFos
        return fos


    def fosMises(self, cls, data, topo, ctrl, **_):
        """FOS that compares the vonMises stress of local stress values
        and their related vonMises stress on the (current) yieldSurface.

        This fos focuses on the influence of shear deformation.

        See L{fosMises<bae.material_01.material.FactorOfSafetyIso.fosMises>}
        for more information. Define further parameters via instance variable
        L{paramMises}.

        @note: This method is time consuming because the vonMises stress
            on yield surface has to be solved iteratively. If you have many
            points and many materials this may become quite painfull.
        """
        msg('## calculating fosMises')
        ## get parameters
        pst = self.paramMises['pst']

        def fosFun(mat, s1, s2, s3, p):
            return mat.fosMises(s1, s2, s3, p, infValue=self.infValue)
        try:
            fos = self._fosMain(data, topo, ctrl, pst, fosFun, False)
        except Exception as e:  # write error and traceback to log
            msg('Error: %s' % type(e))
            msg(traceback.format_exc())
            msg(e.message)
            raise
        return cls(fos)


    def fosHydrostatic(self, cls, data, topo, ctrl, **_):
        """FOS that compares the hydrostatic pressure of local stress values
        and their related pressure on the (current) yieldSurface.

        This fos focuses on the influence of (loose of) confinement.

        See L{fosHydrostatic<bae.material_01.material.FactorOfSafetyIso.fosHydrostatic>}
        for more information. Define further parameters via instance variable
        L{paramHydrostatic}.
        """
        msg('## calculating fosHydrostatic')
        ## get parameters
        pst = self.paramHydrostatic['pst']

        def fosFun(mat, s1, s2, s3, p):
            return mat.fosHydrostatic(s1, s2, s3, p, infValue=self.infValue)

        try:
            fos = self._fosMain(data, topo, ctrl, pst, fosFun, False)
        except Exception as e:  # write error and traceback to log
            msg('Error: %s' % type(e))
            msg(traceback.format_exc())
            msg(e.message)
            raise
        return cls(fos)


    def fosHoekBrown(self, cls, data, topo, ctrl, **_):
        """FOS from HoekBrown approximation. This will call
        L{fosHoekBrown<bae.material_01.FactorOfSafetyIso.fosHoekBrown>}. Define
        further parameters via instance variable L{paramHB}.
        """
        msg('## calculating fosHoekBrown')
        ## get parameters
        pst = self.paramHB['pst']
        measureFromBisection = self.paramHB['measureFromBisection']

        def fosFun(mat, s1, s2, s3, p):
            return mat.fosHoekBrown(s1, s3, p,
                                    measureFromBisection=measureFromBisection,
                                    infValue=self.infValue)

        try:
            fos = self._fosMain(data, topo, ctrl, pst, fosFun, True)
        except Exception as e:  # write error and traceback to log
            msg('Error: %s' % type(e))
            msg(traceback.format_exc())
            msg(e.message)
            raise
        return cls(fos)


    def fosMohrCoulomb(self, cls, data, topo, ctrl, **_):
        """FOS from MohrCoulomb approximation. This will call
        L{fosMohrCoulomb<bae.material_01.FactorOfSafetyIso.fosMohrCoulomb>}.
        Define further parameters via instance variable L{paramMC}.

        Example for outputFields used in odbToVtk and others
         >>> converter = Converter_StandardAndFOS('matprops.inp')
         >>> converter.paramMC['phi'] = 35.0
         >>> converter.paramMC['constrainCohesionToUCS'] = True
         >>> outputFields = [
         ...    ("S3", "scalar", converter.getS3),
         ...    ("FOS_MC35", "scalar", converter.fosMohrCoulomb),
         ...    ]

        @note: If you don't specify phi and c (or constrain c to UCS) this
            method might be very slow because when evaluated on many points
            and materials. In a first step the true yieldSurface has to be
            evaluated and in a second step this will be fitted. This is fast
            for a few points but may be horrible for 25Mio samples.
        """
        msg('## calculating fosMohrCoulomb')

        ## get parameters
        pst = self.paramMC['pst']
        c   = self.paramMC['c']
        phi = self.paramMC['phi']
        constrainCohesionToUCS = self.paramMC['constrainCohesionToUCS']
        fittingRange = self.paramMC['fittingRange']
        measureFromBisection = self.paramMC['measureFromBisection']

        def fosFun(mat, s1, s2, s3, p):
            return mat.fosMohrCoulomb(
                s1, s3, p,
                pstMax=None, infValue=self.infValue,
                miningConvention=True, c=c, phi=phi,
                constrainCohesionToUCS=constrainCohesionToUCS,
                fittingRange=fittingRange,
                measureFromBisection=measureFromBisection)

        try:
            fos = self._fosMain(data, topo, ctrl, pst, fosFun, True)
        except Exception as e:  # write error and traceback to log
            msg('Error: %s' % type(e))
            msg(traceback.format_exc())
            msg(e.message)
            raise
        return cls(fos)


    def damageVariable(self, cls, data, topo, ctrl, **_):
        """Creates a field of pst-dependent damage variable of rock mass which
        is defined as::
           damVar = (UCS(PST) - UCS_Res) / (UCS_Peak - UCS_Res)

        which ranges from 1 (intact strength) to 0 (residual strength). See
        L{bae.material_01.DefaultPst.getDamageVariable} for more details. Fill
        will be flagged with value -1.

        @note: poorly tested yet
        """
        if self._lastMesh is not topo:
            # mesh has changed --> new lookup required
            self.createLookups(
                ctrl.matNumberList,
                data["Material"],
                postData.material.matProps)
            self._lastMesh = topo  # store for next call

        msg('## calculating damage variable')

        # getting required fields
        pst = np.array(data[self.varNameDamage]).astype(float)
        state = np.round(data[self.varNameState]).astype(int)

        # initialize output field
        damVar = np.full_like(pst, -1)
        isRock = (state==1)
        isFill = (state==4)
        damVar[isRock | isFill] = np.nan  # make sure that we got everything

        for matIdx in np.unique(self.matField[isRock]):
            material = self.rockLookup[matIdx]
            ii = isRock & (self.matField == matIdx)
            if not ii.any():
                continue
            damVar[ii] = material.getDamageVariable(pst[ii])

        if self.fixedFillValue is not None:
            damVar[isFill] = -.5
            return cls(damVar)

        for fillIdx in np.unique(self.fillField[isFill]):
            material = self.fillLookup[fillIdx]
            ii = isFill & (self.fillField == fillIdx)
            if not ii.any():
                continue
            damVar[ii] = material.getDamageVariable(pst[ii])

        return cls(damVar)
