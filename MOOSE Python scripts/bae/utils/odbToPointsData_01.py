"""Read interpolated field data from the specified odb and export in different
file formats. Exports regular grid data to legacy vtk files.


Usage:
======

There are two levels of abstraction i.e. two API levels:

Level 1: You can call the all-in-one function odbToVtkGridBox:
 >>> from bae.utils.odbToPointsData_01 import odbToVtkGridBox
 >>> odbToVtkGridBox(...)

Level 2: Or you can call sub tasks individually:
 >>> from bae.utils.odbToPointsData_01 import \
 ...     prepareRemapInterpol, compileRemap, getRemapGridParamString, runRemap
 >>> prepareRemapInterpol(...)
 >>> compileRemap(...)
 >>> remapGridParam = getRemapGridParamString()
 >>> runRemap(..., remapGridParam, ...)
"""

__version__ = "1.14"

_version_history_ = """
1.0 GP: made it a module in bae/utils, derived from:
    3_post/pointsDataFromOdb/odbToVTKGridBox.py version 2.0
1.1 GP: split off getRemapGridParamString()
1.2 GP incompatible change: renamed module
1.3 GP renamed odbToVTKGridBox_wMaxFlds to odbToVTKGridBox. It can now run
  without fieldNamesMaximize options. Added relative fields.
1.4 GP incompatibel change: renamed getMeshCallback to getMeshCallBack
       added matRegionsOption
1.5 GP compileRemap(): don't copy abaqus_v6.env,
       added odbElCentroidsToCsvPoints (from odbCohsvToCsvPoints.py v1.02)
       pipe output from subprocesses to log channel(s)
1.6 AF,GP: made the material name filter optional (and switched off by default)
1.7 GP added rotated grid option
1.8 GP added skipExisting option: skip re-creating existing result files
1.9 GP incompatible change: renamed odbToVTKGridBox to odbToVtkGridBox
1.10 GP changed odbElCentroidsToCsvPoints: boundingBox can be None, filterElset
       can be list of elset names now, added getStructDataValues-function
1.11 GP added frameIdToName argument to odbElCentroidsToCsvPoints
1.12 GP incompatible change: renamed frameIdToName argument to stepFrameToName
        to avoid name clash with vtk_02 functions
1.13 GP in odbElCentroidsToCsvPoints use unique filter elset name
1.14 GP added odbCohsvCentrLocalCoordsToCsv
"""

# Historic note: What follows is an example command-line remap call:
# cd remap; abq6132 remap ../../../2_run/R14/Cadia2015_R14_G01_Q02_S02_D05_M07.odb UEQUIB_1,UEQUIB_2,UEQUIB_3,SDV4,SDV12,S_MIN_PRINCIPAL,SDV15 binary box2_mine_h06.nodalweights box2_mine_h06.gaussptsweights 183,158,251,15089.0,21370.0,4400.0,6.0,6.0,6.0 

import os, sys, shutil, csv, subprocess
from itertools import izip, chain, repeat, count
from bae.mesh_01 import getMesh, InterpolMeshToPoints
from bae.field_01 import createFieldClass, createFieldObject
from bae.vtk_02 import FieldsOnStructPointGrid, FramesFieldsOnUnstructuredPoints
from bae.misc_01 import Container
from bae.log_01 import msg, log


#---------------------------------------

def prepareRemapInterpol(
        gridData, getMeshCallBack, gmCallBackArgs={}, pickleName=None):
    """Make sure the interpolation data files required by remap exist.
    If the files already exist then don't do anything. Otherwise first a
    .interpol file will be read or created. From this the .nodalweights and
    .gaussptsweights files are generated.

    @param gridData: A dict with keys firstPoint, lastPoint and spacing.
      All three values being lists of three floats. This will be passed on as
      keyword-arguments to the L{bae.mesh_01.MeshStructuredPoints} constructor.
    @param getMeshCallBack: Function that returns the mesh (linear tets
      sufficient). Will only be called if required, i.e. if the .interpol file
      is not present.
    @param gmCallBackArgs: arguments dictionary to be passed on to
      getMeshCallBack.
    @param pickleName: Identifier for the interpolation data file names.
      Extension .interpol, .nodalweights and .gaussptsweights will be
      appended automatically to get the actual filenames. May include directory
      names. The path must be relative, i.e. pickleName must not start with a
      slash.
    """
    # prepare grid figures
    # grid = MeshStructuredPoints(**gridData)
    grid = getMesh(**gridData)

    # make sure interpolation data exists
    if not(os.path.isfile("%s.nodalweights" % pickleName)
           and os.path.isfile("%s.gaussptsweights" % pickleName)):
        msg("remap interpolation data files missing."
            " Reading or preparing %s.interpol file." % pickleName)
        interpolation = InterpolMeshToPoints(
            grid=grid,
            getMeshCallBack=getMeshCallBack,
            gmCallBackArgs=gmCallBackArgs,
            pickleName=pickleName
            )

        msg("Converting interpolation data to remap requirements.")
        interpolation.writeRemapParam(pickleName)
    return


def getMeshCallBack_fromOdb(odbPath):
    """Function to get the unstructured FE mesh from the odb.
    Serves as default getMeshCallBack argument in odbToVtkGridBox().
    """
    from bae.odb_02 import OdbReader
    odb = OdbReader(odbPath)
    mesh = odb.getAbqModel(recognizedBlocks="NODE,ELEMENT")
    return mesh


def getMatFieldFromOdb(odbPath, gridData, outputPathMatList,
                       outputPathMatField,
                       pickleName, vtkOutputFormat="binary",
                       nameFilter=None):
    """Read a material number field from the odb.

    @param nameFilter: This is for material name conversion. For example
       if there are region variants that shall be treated as one. I.e. if
       SEDIMENT, SEDIMENT_CAVE and SEDIMENT_UCUT shall all be treated as
       SEDIMENT.

       If this is a function it will be called as a conversion function
       taking the original material name from the odb as only
       argument and returning the material name it should be treated as.

       Example: To take everything up to the first "_" or "-" use:
        >>> matNumberField = getMatFieldFromOdb(
        >>>     ...
        >>>     nameFilter=lambda matname:matname.split("_")[0].split("-")[0])

       That's equivalent to:
        >>> def convert(matname):
        >>>     return matname.split("_")[0].split("-")[0]
        >>> matNumberField = getMatFieldFromOdb(
        >>>     ...
        >>>     nameFilter=convert)

    @returns: material number field
    """

    from bae.odb_02 import OdbReaderInterpolateToGrid
    msg("Opening the odb %s" % odbPath)
    odb = OdbReaderInterpolateToGrid(odbPath)
    grid = odb.initOutputGrid(interpolPickleName=pickleName, **gridData)
    assert len(grid)==len(odb.interpolation.elemCoords)

    (matNamesList, matNumberField) = odb.getMatFieldFromOdb(
        nameFilter=nameFilter)

    # write material number to names list
    msg("Writing material list to the csv file %s"
        % outputPathMatList)
    odb.writeMatNamesToCsv(outputPathMatList, matNamesList)

    # store material number field for subsequent runs
    output = FieldsOnStructPointGrid(mesh=grid)
    # converting to float
    matNumberField = type(matNumberField)(
        float(x) for x in matNumberField)
    output.updateFields(matNumberField)
    msg("Writing material data field %s" % outputPathMatField)
    output.toVtk(outputPathMatField,
                 description="Abq-Grid-Data matfield",
                 outputFormat=vtkOutputFormat)
    return matNumberField


def callAbqPython(abaqusExecutable, pgm, **kwargs):
    """Run some python commands as subprocess with Abaqus Python.

    @param abaqusExecutable: "abaqus" or "abq6132" or such.

    @param pgm: The python code to be executed in the subprocess. Parameter
    templates in the form of "{arg}" will be replaced by a corresponding
    keyword argument with key "arg".

    Note: As this string is being processed by str.format() you have to
    "excape" curly braces by doubling them. I.e. the python expression
    "mydict = {}" must be given as "mydict = {{}}" in this string!

    @param kwargs: following keyword arguments are required as parameter
    in the pgm-string.

    @Note: You should have only one occurence of each parameter template
    in the pgm-string because it is replaced by "pickle.loads(sys.argv[xxx])".
    You may want to decode the pickled argument only once, maybe store it in
    a local variable for multiple uses.
    """
    import cPickle as pickle

    # add two imports and the log initialization at the beginning
    pgm = "import sys\nimport cPickle as pickle\n" + pgm

    # define order of parameters and put replacements into the pgm-string
    parameterNames = sorted(kwargs)
    toReplace = dict()
    for i, name in enumerate(parameterNames):
        toReplace[name] = "pickle.loads(sys.argv[%d])" % (i+1)

    rc = subprocess.call(
        [abaqusExecutable, "python", "-c", pgm.format(**toReplace)]
        + [pickle.dumps(kwargs[name]) for name in parameterNames],
        stdout=log.getLogFileForSubProcess(), stderr=subprocess.STDOUT)
    return rc


def compileRemap(displacementFieldName="U",
                 sdvFieldNames=["SDV4","SDV12","SDV15"],
                 abaqusExecutable="abaqus"):
    """prepare the remap binary. Copy sources, modify hard wired parameters
    and compile.

    @param displacementFieldName: string like "U" or "URF20" or "UR" or "UEQUIB"
    @param sdvFieldNames: list of three strings
    @param abaqusExecutable: like "abq6132"
    """

    if len(sdvFieldNames)!=3:
        raise ValueError(
            "sdvFieldNames must contain 3 valid field names! We have: %s,"
            " len=%d." % (sdvFieldNames, len(sdvFieldNames)))

    # copying the remap source to the current working directory
    srcname = os.path.join(os.path.dirname(__file__), "remap")
    msg("Copying remap source from %s to the current working directory."
        % srcname)
    if not os.path.isdir("remap"):
        os.mkdir("remap")
    for fname in ["remap.cpp", "helper_functions.cpp",
                  "visit_writer_mod.c", "visit_writer.h", "headers.h",
                  # "abaqus_v6.env", "hardCodedFieldNames.cpp",
                  ]:
        shutil.copy(os.path.join(srcname, fname), "remap")

    # compiling/building remap binary
    msg("Compiling/building remap binary.")
    os.chdir("remap")
    output = open("hardCodedFieldNames.cpp", "w")
    output.write(
        'const std::string varName_U = "%s";\n' % displacementFieldName)
    for i in range(3):
        output.write('const std::string varName_U_%1d = "%s_%1d";\n'
                     % (i+1, displacementFieldName, i+1))
    for i, fn in enumerate(sdvFieldNames):
        output.write('const std::string varName_SDV%1d = "%s";\n' % (i+1, fn))
    output.close()
    del output

    subprocess.check_call(
        [abaqusExecutable, "make", "j=remap"],
        stdout=log.getLogFileForSubProcess(), stderr=subprocess.STDOUT)
    os.chdir("..")

    return


def getRemapGridParamString(gridData):
    """Generate the remap parameter that describes the output grid.

    @param gridData: A dict with keys firstPoint, lastPoint and spacing.
      All three values being lists of three floats. This will be passed on as
      keyword-arguments to the L{bae.mesh_01.MeshStructuredPoints} constructor.

    @returns: the string in the remap parameter that describes the output grid.
    This can be used verbatim in the remap.parameter file.
    """

    # prepare grid figures
    # grid = MeshStructuredPoints(**gridData)
    grid = getMesh(**gridData)

    # grid data string: nx,ny,nz,x,y,z,sx,sy,sz
    gridNumbers = list()
    gridNumbers.extend(grid.gridPtNb)
    try:
        origin = grid.originBackRotated
    except AttributeError:
        origin = grid.origin
    gridNumbers.extend(origin)
    gridNumbers.extend(grid.spacing)
    remapGridParam = ",".join(str(x) for x in gridNumbers)
    del gridNumbers
    return remapGridParam


def runRemap(
        odbPath, fieldNames, frameList,
        pickleName, remapGridParam,
        outputDirName, fileNamePrefix,
        abaqusExecutable="abaqus",
        remapParameterFileName="remap.parameter",
        outputFormat="binary",
        skipExisting=False):

    """launch the remap binary. Prepare the necessary parameter file, check
    output directories and call the binary.

    Vtk file names will be fileNamePrefix+"_F2xxxd.vtk" with xxx being the
    frame number and the 2 standing for step-2.

    @param odbPath: complete path to the odb including .odb extension. Can be
      relative or absolute.
    @param fieldNames: example: ["UR_1","UR_2","UR_3","SDV4","S_MIN_PRINCIPAL"]
    @param frameList: example: [ (2, 24), (2, 27) ]
    @param pickleName: Identifier for the interpolation data file names.
      Extension .interpol, .nodalweights and .gaussptsweights will be
      appended automatically to get the actual filenames. May include directory
      names. The path must be relative, i.e. pickleName must not start with a
      slash.
    @param remapGridParam: The string as recieved from L{prepareRemapInterpol}
      (grid data string: nx,ny,nz,x,y,z,sx,sy,sz)
    @param outputDirName: Vtk files will end up in this folder. Can be
      relative or absolute path.
    @param fileNamePrefix: example: "Perse2015_R01_boxD3_h5"
    @param abaqusExecutable: like "abq6132"
    @param outputFormat: "binary" or "ascii"
    @param skipExisting: If False then overwrite possible result files that
      already exist unconditionally. If True then don't re-create result files
      that already exist.
    """

    # file name templates from outputDirName and fileNamePrefix
    outputFileNameTemplate = os.path.join(
        outputDirName, fileNamePrefix+"_F%d%03d.vtk")

    # frames
    if any(x!=2 for x,y in frameList):
        raise ValueError("ERROR: This version can only handle odb step 2!")

    # look for existing files and correct frameList
    if skipExisting:
        newFrameList = [
            (stepNb, frameNb)
            for stepNb, frameNb in frameList
            if not os.path.isfile(outputFileNameTemplate % (stepNb, frameNb))]
        removed = len(frameList) - len(newFrameList)
        if removed:
            msg("WARNING: %d/%d files ordered to be created already exist in"
                " folder %s. Old files will be preserved and have been removed"
                " from the frame list."
                % (removed, len(frameList), outputDirName))
            frameList = newFrameList
        if not frameList:
            msg("Skipping remap alltogether because all files are already"
                " there.")
            return

    # correct outputFileNameTemplate for remap being run in subfolder
    if not os.path.isabs(outputDirName):
        # add .. based on ./remap - subdirectory
        outputFileNameTemplate = os.path.join("..", outputFileNameTemplate)

    # correct odb path for remap being run in subfolder
    if not os.path.isabs(odbPath):
        # add .. based on ./remap - subdirectory
        odbPath = os.path.join("..", odbPath)

    # prepare combined frame, outputfile list
    remapFrameList = [
        ("%d,"+outputFileNameTemplate) % (frameNb, stepNb, frameNb)
        for stepNb, frameNb in frameList]

    cmd = [abaqusExecutable, "remap", remapParameterFileName]

    # create remap.parameter
    f = open(os.path.join("remap", remapParameterFileName), "w")
    f.write("\n".join([
        odbPath,
        ",".join(fieldNames),
        outputFormat,
        "../%s.nodalweights" % pickleName,
        "../%s.gaussptsweights" % pickleName,
        remapGridParam,
        ] + remapFrameList) )
    f.write("\n")
    f.close()
    del f
    msg("wrote parameter file %s" % remapParameterFileName)

    msg("Starting remap. Command:\n%s" % " ".join(cmd))
    os.chdir("remap")
    subprocess.check_call(
        cmd, stdout=log.getLogFileForSubProcess(), stderr=subprocess.STDOUT)
    os.chdir("..")
    msg("Remap finished.")
    return


def odbToVtkGridBox_simple(
        odbPath, outputDirName, fileNamePrefix,
        projectPrefix, meshVersion,
        boxName, gridData,
        getMeshCallBack=None, gmCallBackArgs={},
        fieldNames=None, frameList=None,
        displacementFieldName=None, sdvFieldNames=None,
        abaqusExecutable="abaqus",
        remapParameterFileName="remap.parameter",
        interpolationDataDir="interpolationdata",
        outputFormat="binary"
        ):
    """Simple version for documentation purpose. Use odbToVtkGridBox() instead.

    Read interpolated field data from the specified odb at points on the
    specified regular grid and write the data to legacy vtk files.

    @param odbPath: complete path to the odb including .odb extension. Can be
      relative or absolute.
    @param outputDirName: Vtk files will end up in this folder. Can be
      a relative or an absolute path.
    @param fileNamePrefix: For result files and preliminary data files.
      Example: "Perse2015_R01_boxD3_h5"
    @param projectPrefix: Needed for the interpolation data filename.
    @param meshVersion: Needed for the interpolation data filename.
    @param boxName: Used to identify the regular grid, in particular the
      corresponding interpolation data.

    @param gridData: A dict with keys firstPoint, lastPoint and spacing.
      All three values being lists of three floats. This will be passed on as
      keyword-arguments to the L{bae.mesh_01.MeshStructuredPoints} constructor.

    @param getMeshCallBack: Optional function that returns the mesh (linear tets
      sufficient). Will only be called if required, i.e. if the .interpol file
      is not present. If not specified (or None) then
      L{getMeshCallBack_fromOdb} will be used and gmCallBackArgs will be
      adapted accoringly.

      Example: load the FE mesh from Abaqus input files:
       >>> def getMeshCallBack_fromInp(filename):
       >>>     from bae.abq_model_02 import Model
       >>>     mesh = Model().read(filename, recognizedBlocks="NODE,ELEMENT")
       >>>     return mesh
       >>>
       >>> odbToVtkGridBoxSimple(
       >>>     ...
       >>>     getMeshCallBack=getMeshCallBack_fromInp,
       >>>     gmCallBackArgs={"filename": config.tetLPath},
       >>>     ...

    @param gmCallBackArgs: arguments dictionary to be passed on to
      getMeshCallBack. E.g. {"odbPath": "../../2_run/R01/MyProj2015_R01.odb"}

    @param fieldNames: List of field names to be exported to the output vtk.
      Must include all fields in argument fieldNamesMaximize.
      Example: ["UR_1","UR_2","UR_3","SDV4","S_MIN_PRINCIPAL"]
    @param frameList: List of output frames as (step number, frame number)
      tuples. Example: [ (2, 24), (2, 27) ]

    @param displacementFieldName: string like "U" or "URF20" or "UR" or "UEQUIB"
    @param sdvFieldNames: list of three strings, e.g. ["SDV4","SDV12","SDV15"]

    @param abaqusExecutable: ... to compile and run remap. Example: "abq6132"
    @param remapParameterFileName: Calling the remap executable requires an
      intermediate parameter file. For debugging purposes it's not deleted
      automatically. The file will be in the remap subfolder. This is the file
      name.
    @param interpolationDataDir: .interpol, .nodalweights and .gaussptsweights
      files will be in that subfolder
    @param outputFormat: "binary" or "ascii"
    """

    # process default arguments getMeshCallBack and (possibly) gmCallBackArgs
    if getMeshCallBack is None:
        getMeshCallBack = getMeshCallBack_fromOdb
        gmCallBackArgs = {"odbPath": odbPath}

    # check argument frameList (multi-step not implemented yet.)
    if frameList[0][0] != frameList[-1][0]:
        raise NotImplementedError(
            "ERROR!!! Multi-step results not working yet. Got different"
            " odb steps in the frameList: %s"
            % sorted(set(x for x,y in frameList)))

    # check output directory
    if outputDirName and not os.path.isdir(outputDirName):
        os.makedirs(outputDirName)

    # pickle data name and directory to store it in
    pickleName = "_".join((projectPrefix, meshVersion, boxName))
    if interpolationDataDir:
        if not os.path.isdir(interpolationDataDir):
            os.makedirs(interpolationDataDir)
        pickleName = os.path.join(interpolationDataDir, pickleName)

    # prepare interpolation data
    prepareRemapInterpol(
        gridData=gridData, getMeshCallBack=getMeshCallBack,
        gmCallBackArgs=gmCallBackArgs, pickleName=pickleName)

    # compile remap
    # Note: Compilation is done unconditionally. Saving the compile time could
    # be achieved but we'd have to check that neither the abaqus version nor
    # the hard-wired parameters have changed. I think it's not worth it.
    msg("Compiling remap executable...")
    compileRemap(
        displacementFieldName, sdvFieldNames, abaqusExecutable)

    # run Remap
    msg("Running remap binary...")
    remapGridParam = getRemapGridParamString(gridData)
    runRemap(
        odbPath, fieldNames, frameList,
        pickleName, remapGridParam,
        outputDirName, fileNamePrefix,
        abaqusExecutable,
        remapParameterFileName,
        outputFormat=outputFormat,
        skipExisting=True)


def odbToVtkGridBox(
        odbPath, outputDirName, fileNamePrefix,
        projectPrefix, meshVersion,
        boxName, gridData,
        getMeshCallBack=None, gmCallBackArgs={},
        fieldNames=None, fieldNamesMaximize=None,
        relativeFieldNames=None, relativeToFrame=(2,0),
        matRegionsOption=False, matRegionsNameFilter=None,
        frameList=None,
        displacementFieldName=None, sdvFieldNames=None,
        abaqusExecutable="abaqus",
        remapParameterFileName="remap.parameter",
        interpolationDataDir="interpolationdata",
        scratchDataDir="prelimVTK",
        outputFormat="binary"
        ):
    """Read interpolated field data from the specified odb at points on the
    specified regular grid and write the data to legacy vtk files.

    Option fieldNamesMaximize: Maximize selected fields between output frames.
    Data will be read from the specified odb from the first frame specified in
    frameList to its last including all intermediate frames not just those
    contained in frameList. For frames not specified in frameList read only
    fields in fieldNamesMaximize. For frames listed in frameList -except for
    the first- generate output vtks. Take the maximum value since the previous
    output frame -or the first frame in frameList- for fields listed in
    fieldNamesMaximize.

    Option relativeFieldNames: For certain fields (e.g. displacement) output
    values relative to a reference frame.

    Option matRegionsOption: Add material region number field to every output
    vtk file.

    @param odbPath: Complete path to the odb including .odb extension. Can be
      relative or absolute.
    @param outputDirName: Vtk files will end up in this folder. Can be
      a relative or an absolute path.
    @param fileNamePrefix: For result files and preliminary data files.
      Example: "Perse2015_R01_boxD3_h5"
    @param projectPrefix: Needed for the interpolation data filename.
    @param meshVersion: Needed for the interpolation data filename.
    @param boxName: Used to identify the regular grid, in particular the
      corresponding interpolation data.

    @param gridData: A dict with keys firstPoint, lastPoint and spacing and
      optionally rotString.
      All three values being lists of three floats. This will be passed on as
      keyword-arguments to the L{bae.mesh_01.MeshStructuredPoints} constructor.

    @param getMeshCallBack: Optional function that returns the mesh (linear
      tets sufficient). Will only be called if required, i.e. if the .interpol
      file is not present. If not specified (or None) then
      L{getMeshCallBack_fromOdb} will be used and gmCallBackArgs will be
      adapted accoringly.

      Example: load the FE mesh from Abaqus input files:
       >>> def getMeshCallBack_fromInp(filename):
       >>>     from bae.abq_model_02 import Model
       >>>     mesh = Model().read(filename, recognizedBlocks="NODE,ELEMENT")
       >>>     return mesh
       >>>
       >>> odbToVtkGridBox(
       >>>     ...
       >>>     getMeshCallBack=getMeshCallBack_fromInp,
       >>>     gmCallBackArgs={"filename": config.tetLPath},
       >>>     ...

    @param gmCallBackArgs: arguments dictionary to be passed on to
      getMeshCallBack. E.g. {"odbPath": "../../2_run/R01/MyProj2015_R01.odb"}

    @param fieldNames: List of field names to be exported to the output vtk.
      Must include all fields in argument fieldNamesMaximize.
      Example: ["UR_1","UR_2","UR_3","SDV4","S_MIN_PRINCIPAL","SDV12"]
    @param fieldNamesMaximize: Optional list of field names to take the maximum
      value of intermediate frames of. For maxRER. Example: ["SDV12"]
    @param relativeFieldNames: Optional. From fields listed in this argument
      the value of the time frame relativeToFrame is being subtracted before
      export. E.g. for relative displacements (relative to the equilibrium
      step) let relativeFieldNames=["U_1", "U_2", "U_3"] and
      relativeToFrame=(2.0).
    @param relativeToFrame: Reference time frame for the relativeFieldNames
      option.
    @param matRegionsOption: Write MATERIAL REGIONS data as well?
    @param matRegionsNameFilter: This is for material name conversion. For
      example if there are region variants that shall be treated as one. I.e.
      if SEDIMENT, SEDIMENT_CAVE and SEDIMENT_UCUT shall all be treated as
      SEDIMENT.

      If this is a function it will be called as a conversion function
      taking the original material name from the odb as only
      argument and returning the material name it should be treated as.

      Example: To take everything up to the first "_" or "-" use as
      matRegionsNameFilter:
       >>> lambda matname:matname.split("_")[0].split("-")[0]

      That's equivalent to:
       >>> def convert(matname):
       >>>     return matname.split("_")[0].split("-")[0]
       >>> odbToVtkGridBox(..., matRegionsNameFilter=convert, ...)

    @param frameList: List of output frames as (step number, frame number)
      tuples. Example: [ (2, 24), (2, 27) ]

      In case fieldNamesMaximize is given this list is treated slightly
      different: The first frame in this list is the first frame to
      consider at all. No output is generated for this first frame! For
      subsequent frames output will be generated.
      Example: [ (2, 20), (2, 24), (2, 27) ] will generate output for frame
      24 and frame 27. Not for frame 20! The results for frame 24 will include
      maximum RER from frame 20 to frame 24 (both included). For frame 27 we'll
      have maximum RER from frame 25 to frame 27 (both included). Provided
      that RER (i.e. SDV12 for example) is listed in fieldNamesMaximize.

    @param displacementFieldName: string like "U" or "URF20" or "UR" or "UEQUIB"
    @param sdvFieldNames: list of three strings, e.g. ["SDV4","SDV12","SDV15"]

    @param abaqusExecutable: ... to compile and run remap. Example: "abq6132"
    @param remapParameterFileName: Calling the remap executable requires an
      intermediate parameter file. For debugging purposes it's not deleted
      automatically. The file will be in the remap subfolder. This is the file
      name.
    @param interpolationDataDir: .interpol, .nodalweights and .gaussptsweights
      files will be in that subfolder
    @param scratchDataDir: intermediate / preliminary  results files (e.g. RER
      for intermediate frames) will be in that (sub-)folder. Not needed after
      the run, except for debugging purposes.
    @param outputFormat: "binary" or "ascii"
    """

    # process default arguments getMeshCallBack and (possibly) gmCallBackArgs
    if getMeshCallBack is None:
        getMeshCallBack = getMeshCallBack_fromOdb
        gmCallBackArgs = {"odbPath": odbPath}

    # check arguments fieldNamesMaximize and fieldNames
    if fieldNamesMaximize:
        missingFields = [fn for fn in fieldNamesMaximize
                         if fn not in fieldNames]
        if missingFields:
            raise ValueError(
                "ERROR: %d fields in argument fieldNamesMaximize are not"
                " contained in argument fieldNames. That does not make sense."
                % (len(missingFields)))
        del missingFields

    # check argument relativeFieldNames
    if relativeFieldNames:
        missingFields = [fn for fn in relativeFieldNames
                         if fn not in fieldNames]
        if missingFields:
            raise ValueError(
                "ERROR: %d fields in argument relativeFieldNames are not"
                " contained in argument fieldNames. That does not make sense."
                % (len(missingFields)))
        del missingFields

    # check argument frameList (multi-step not implemented yet.)
    # because remap can't handle it (yet)
    if frameList[0][0] != frameList[-1][0]:
        raise NotImplementedError(
            "ERROR!!! Multi-step results not working yet. Got different"
            " odb steps in the frameList: %s"
            % sorted(set(x for x,y in frameList)))

    # check output directories
    if outputDirName and not os.path.isdir(outputDirName):
        os.makedirs(outputDirName)
    if ((fieldNamesMaximize or relativeFieldNames or matRegionsOption)
        and scratchDataDir and not os.path.isdir(scratchDataDir)):
        os.makedirs(scratchDataDir)

    # write csv file describing the resulting vtk files
    outputFileNameVtkDescription = os.path.join(
        outputDirName, fileNamePrefix+"_vtk_description.csv")
    output = csv.writer(open(outputFileNameVtkDescription, "wb"))
    output.writerow(["index", "component"])
    for i, fieldName in enumerate(fieldNames):
        output.writerow([i+1, fieldName])
    if matRegionsOption:
        output.writerow([len(fieldNames)+1,"matcode"])
    del output
    msg("Wrote vtk file description file %s" % outputFileNameVtkDescription)

    # pickle data name and directory to store it in
    pickleName = "_".join((projectPrefix, meshVersion, boxName))
    if interpolationDataDir:
        if not os.path.isdir(interpolationDataDir):
            os.makedirs(interpolationDataDir)
        pickleName = os.path.join(interpolationDataDir, pickleName)

    # prepare interpolation data
    try:
        prepareRemapInterpol(
            gridData=gridData, getMeshCallBack=getMeshCallBack,
            gmCallBackArgs=gmCallBackArgs, pickleName=pickleName)
    except ImportError, exc:
        # Make sure that the reason is: We are not running Abaqus python.
        if "abaqus" not in exc.args[0]:
            raise
        # in case of failure get the mesh from the odb
        msg("Launching a separate abaqus python process to get"
            " the mesh from the odb.")
        rc = callAbqPython(
            abaqusExecutable, """
from bae.utils.odbToPointsData_01 import \
    prepareRemapInterpol, getMeshCallBack_fromOdb
prepareRemapInterpol(
    gridData={gridData},
    getMeshCallBack=getMeshCallBack_fromOdb,
    gmCallBackArgs={{"odbPath": {odbPath}}},
    pickleName={pickleName})
""",
            gridData=gridData,
            odbPath=odbPath,
            pickleName=pickleName)
        if rc!=0:
            msg("ERROR running prepareRemapInterpol with %s python."
                % abaqusExecutable)
            sys.exit(rc)


    remapGridParam = getRemapGridParamString(gridData)

    # compile remap
    # Note: Compilation is done unconditionally. Saving the compile time could
    # be achieved but we'd have to check that neither the abaqus version nor
    # the hard-wired parameters have changed. I think it's not worth it.
    msg("Compiling remap executable...")
    compileRemap(
        displacementFieldName, sdvFieldNames, abaqusExecutable)

    # run Remap for frameList and fieldNames
    msg("Running remap binary for output frames and all fields ...")
    if fieldNamesMaximize or relativeFieldNames or matRegionsOption:
        runRemap(
            odbPath, fieldNames, frameList,
            pickleName, remapGridParam,
            scratchDataDir, fileNamePrefix+"_prelim",
            abaqusExecutable,
            remapParameterFileName,
            outputFormat=outputFormat,
            skipExisting=True)
    else:
        runRemap(
            odbPath, fieldNames, frameList,
            pickleName, remapGridParam,
            outputDirName, fileNamePrefix,
            abaqusExecutable,
            remapParameterFileName,
            outputFormat=outputFormat,
            skipExisting=True)

    # run Remap for extra data: fields to maximize (like RER)
    if fieldNamesMaximize:
        # extract first entry from frameList as frameStartMaximize
        # (removing this first entry from frameList)
        frameStartMaximize = frameList[0]
        frameList = frameList[1:]

        # create frame list for fields in fieldNamesMaximize
        if frameList[0][0] != frameList[-1][0]:
            msg("WARNING!!! Multi-step results not working yet.")
        stepNb = frameList[0][0]
        frameListMaximize = [
            (stepNb, frNb)
            for frNb in range(frameStartMaximize[1], frameList[-1][1]+1)
            if frNb not in set(y for x,y in frameList)]

        # combined frame list frameListAll
        frameListAll = sorted(frameListMaximize + frameList)

        # prepare frameList for fast lookup in postprocessing main loop
        frameList = set(frameList)

        # run Remap for frameListMaximize and fieldNamesMaximize
        msg("Running remap binary for intermediate frames and fields to"
            " maximize: %s ..." % (", ".join(fieldNamesMaximize)))
        runRemap(
            odbPath, fieldNamesMaximize, frameListMaximize,
            pickleName, remapGridParam,
            scratchDataDir, fileNamePrefix+"_prelim",
            abaqusExecutable,
            remapParameterFileName,
            outputFormat=outputFormat,
            skipExisting=True)
    else:
        frameListAll = frameList

    # run Remap for extra data: reference data for relative fields (like UR)
    if relativeFieldNames:
        msg("Running remap binary for reference frame and fields"
            " %s ..." % (", ".join(relativeFieldNames)))
        runRemap(
            odbPath, relativeFieldNames, [relativeToFrame],
            pickleName, remapGridParam,
            scratchDataDir, fileNamePrefix+"_refData",
            abaqusExecutable,
            remapParameterFileName,
            outputFormat=outputFormat,
            skipExisting=True)

        # read the data into memory
        refDataFileName = os.path.join(
            scratchDataDir,
            "%s_refData_F%d%03d.vtk"
            % (fileNamePrefix, relativeToFrame[0], relativeToFrame[1]))
        msg("Reading reference displacements from %s." % refDataFileName)
        refData = FieldsOnStructPointGrid().fromVtk(refDataFileName).data
    else:
        refData = dict()

    # mat regions preparation
    if matRegionsOption:
        outputPathMatList = os.path.join(
            outputDirName, fileNamePrefix+"_matList.csv")
        outputPathMatField = os.path.join(
            scratchDataDir,
            "%s_matfield.vtk" % fileNamePrefix )

        # Material number field: How do we get it?
        # First try to read a previously created vtk file
        # Second try to read it from the odb. This will fail if we're not
        # running abaqus python.
        # So third: read it from the odb in a separate abaqus python process
        if os.access(outputPathMatField, os.R_OK):
            #-- first option: get it from a vtk file
            msg("Reading material number field from %s"
                % outputPathMatField)
            vtk = FieldsOnStructPointGrid().fromVtk(outputPathMatField)
            matNumberField = vtk.data["Material"]

        else:
            # get the material data from the odb
            msg("Reading material properties from the odb.")
            try:
                #-- second option: get it from the odb
                matNumberField = getMatFieldFromOdb(
                    odbPath, gridData, outputPathMatList,
                    outputPathMatField,
                    pickleName=pickleName, vtkOutputFormat=outputFormat,
                    nameFilter=matRegionsNameFilter)
            except ImportError, exc:
                # Make sure that the reason is:
                # We are not running Abaqus python.
                if "abaqus" not in exc.args[0]:
                    raise
                #-- third option: start a separate process to get it
                msg("Launching a separate abaqus python process to get"
                    " the material data from the odb.")
                rc = callAbqPython(
                    abaqusExecutable, """
from bae.utils.odbToPointsData_01 import getMatFieldFromOdb
getMatFieldFromOdb(
    {odbPath}, {gridData}, {outputPathMatList},
    {outputPathMatField}, {pickleName}, {vtkOutputFormat},
    nameFilter={matRegionsNameFilter})
""",
                    odbPath=odbPath,
                    gridData=gridData,
                    outputPathMatList=outputPathMatList,
                    outputPathMatField=outputPathMatField,
                    pickleName=pickleName,
                    vtkOutputFormat=outputFormat,
                    matRegionsNameFilter=matRegionsNameFilter)
                if rc!=0:
                    msg("ERROR running getMatFieldFromOdb with %s python."
                        % abaqusExecutable)
                    sys.exit(rc)
                msg("Reading material number field from %s"
                    % outputPathMatField)
                vtk = FieldsOnStructPointGrid().fromVtk(outputPathMatField)
                matNumberField = vtk.data["Material"]
                del vtk

    # postprocess/combine preliminary data, loop over all frames
    if fieldNamesMaximize or relativeFieldNames or matRegionsOption:
        msg("Now combining the preliminary results to final ones ...")
        fieldsMaximize = dict()
        for stepNb, frameNb in frameListAll:

            frameName = "step %d, frame %d" % (stepNb, frameNb)
            prelimFileName = os.path.join(
                scratchDataDir,
                "%s_prelim_F%d%03d.vtk" % (fileNamePrefix, stepNb, frameNb))

            msg("Processing %s. Reading vtk file %s"
                % (frameName, prelimFileName))
            vtk = FieldsOnStructPointGrid().fromVtk(prelimFileName)

            if fieldNamesMaximize:
                # calculate maximimum up to now
                for fieldName in fieldNamesMaximize:
                    try:
                        fldMax = fieldsMaximize[fieldName]
                    except KeyError:
                        # no previous max value. Take first as current max.
                        fieldsMaximize[fieldName] = vtk.data[fieldName]
                    else:
                        msg("Updating max value for %s..." % fieldName)
                        fieldsMaximize[fieldName] = type(fldMax)(
                            max(a,b)
                            for a, b in izip(fldMax, vtk.data[fieldName]))

            if fieldNamesMaximize and (stepNb, frameNb) not in frameList:
                # skip output on this frame if it's not requested in frameList
                continue
            else:
                outputFileName = os.path.join(
                    outputDirName,
                    "%s_F%d%03d.vtk" % (fileNamePrefix, stepNb, frameNb))

                for fieldName in fieldNames:

                    if fieldNamesMaximize and (fieldName in fieldsMaximize):
                        # take maximum value if available
                        field = fieldsMaximize[fieldName]
                    else:
                        field = vtk.data[fieldName]

                    # take relative value if reference data is available
                    try:
                        refField = refData[fieldName]
                    except KeyError:
                        pass
                    else:
                        field = type(field)(
                            (a-b) for a, b in izip(field, refField))

                    # remove and add new to force order according to fieldNames
                    del vtk.data[fieldName]
                    vtk.data[fieldName] = field

                if matRegionsOption:
                    # add material number field
                    vtk.updateFields(matNumberField)

                msg("Writing result vtk %s" % outputFileName)
                vtk.toVtk(outputFileName,
                          description="Abq-Grid-Data %s" % frameName,
                          outputFormat=outputFormat)

                # reset all maximum values
                fieldsMaximize.clear()


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


def odbElCentroidsToCsvPoints(
        odbPath, outputDirName, fileNamePrefix, addFrameNames=False,
        boundingBox=None, filterElsetName='COHSV',
        fieldNames=None, fieldNamesMaximize=None,
        frameList=None,
        outputType="CsvPerFieldRowPtColFr",
        stepFrameToName=None,
        structDataFieldNames=None,
        getStructDataValues=None,
        ):
    """Read field data from the specified odb at element centroids
    and write the data to csv files.

    Usage:
     >>> odbPath = config.odbPath
     >>> outputDirName = "CSV/Cadia2016033_R01_COHSV"
     >>> fileNamePrefix = "Cadia2016033_R01_COHSV"
     >>> fieldNames = ["PST", "RER"]
     >>> frameList = [(2, i) for i in range(10, 20)]
     >>>
     >>> def getStructDataValues(elem, centroid, elsetName):
     >>>     if centroid[2]>200:
     >>>         return [0,]  # ignore all above z=200
     >>>     return [int(elsetName=="COHSV_THRUST_FAULT"),]
     >>>
     >>> odbElCentroidsToCsvPoints(
     >>>     odbPath, outputDirName, fileNamePrefix, addFrameNames=False,
     >>>     boundingBox=None, filterElsetName="COHSV",
     >>>     fieldNames=fieldNames, frameList=frameList,
     >>>     outputType="CsvPerFieldRowPtColFr",
     >>>     structDataFieldNames=["ThrustFlt",],
     >>>     getStructDataValues=getStructDataValues,
     >>>     )

    Option fieldNamesMaximize: Maximize selected fields between output frames.
    Data will be read from the specified odb from the first frame specified in
    frameList to its last including all intermediate frames not just those
    contained in frameList. For frames not specified in frameList read only
    fields in fieldNamesMaximize. For frames listed in frameList -except for
    the first- generate output vtks. Take the maximum value since the previous
    output frame -or the first frame in frameList- for fields listed in
    fieldNamesMaximize.

    @param odbPath: Complete path to the odb including .odb extension. Can be
      relative or absolute.
    @param outputDirName: Vtk files will end up in this folder. Can be
      a relative or an absolute path.
    @param fileNamePrefix: For result files and preliminary data files.
      Example: "Perse2015_R01_boxD3_h5"
    @param addFrameNames: Option to add frame names to output file names (in
      case of output types having one csv file per frame) or to header lines
      (in all-frames-in-one-csv's).
    @param boundingBox: L{BoundingBox<bae.misc_01.BoundingBox>}-object (or
      [[xmin,ymin,zmin],[xmax,ymax,zmax]]). Consider only elements whose
      centroid lies within this box. If None then all elements will be
      considered.
    @param filterElsetName: Name of the elset to be exported. Defaults to
      "COHSV". If it's a list then consider all elements from all elsets in
      this list.
    @param fieldNames: List of field names to be exported to the output vtk.
      Must include all fields in argument fieldNamesMaximize.
      Example: ["UR_1","UR_2","UR_3","SDV4","S_MIN_PRINCIPAL","SDV12"]

      Note that "S_MIN_PRINCIPAL" is renamed to "S" in the output and
      "SDV_XXXXX" is renamed to "XXXXX".

    @param fieldNamesMaximize: Optional list of field names to take the maximum
      value of intermediate frames of. For maxRER. Example: ["SDV12"]

    @param frameList: List of output frames as (step number, frame number)
      tuples. Example: [ (2, 24), (2, 27) ]

      In case fieldNamesMaximize is given this list is treated slightly
      different: The first frame in this list is the first frame to
      consider at all. No output is generated for this first frame! For
      subsequent frames output will be generated.
      Example: [ (2, 20), (2, 24), (2, 27) ] will generate output for frame
      24 and frame 27. Not for frame 20! The results for frame 24 will include
      maximum RER from frame 20 to frame 24 (both included). For frame 27 we'll
      have maximum RER from frame 25 to frame 27 (both included). Provided
      that RER (i.e. SDV12 for example) is listed in fieldNamesMaximize.

    @param outputType:
      See L{bae.vtk_02.FramesFieldsOnUnstructuredPoints.initOutput}
       - "OneCsvRowPtColFldFr": one file for the whole data.
         Columns are: x,y,z, constant fields if given..., then
         field 1 frame 1, field 1 frame 2, field 1 frame 3, ...,
         field 2 frame 1, field 2 frame 2, field 2 frame 3, ...
         ... more fields ...
         There is exactly one row per point.
       - "CsvPerFrameRowPtColFld": one file per frame.
         Columns are: x,y,z, then constant fields, then remaining fields
         There is one row per data point and one file per odb frame.
         This format somehow corresponds to our usual vtk-gridbox data.
       - "CsvPerFieldRowPtColFr": one file per field.
         Columns are: x,y,z, then one column per frame
         There is exactly one row per point and one file per output field.
         This is the traditional style of pointsDataCohsv.py. It's
         identical to the "OneCsvRowPtColFldFr"-type in case of only one
         field being asked for.
       - "OneCsvRowFrColPtField": one file for the whole data.
         Columns are: frame, x,y,z, field1, field2,... x,y,z, field1,
         field2,... I.e. one row per frame and different points in columns.
         For each point we have the points coordinates first and then the
         field data.

    @param stepFrameToName: Defines how to build frame names in the csv file
      names or in the csv header line from a (stepNb, frameNb)-tuple.
      Might be a function (stepNb, frameNb)-tuple --> frameName or a template
      string. Default if not given, i.e. None: "F%d%03d".

      For restart analyses use something like this (We only have step 2 and 3,
      the latter sometimes referred to as 1):
       >>> def stepFrameToName(stepNb, frameNb):
       >>>     if stepNb==2:
       >>>         return "F2%03d" % frameNb
       >>>     else:
       >>>         return "F3%03d" % frameNb
       >>>
       >>> odbToVtkGridBox( ..., stepFrameToName=stepFrameToName, ...)

    @param structDataFieldNames: a list of field names of time constant fields
      that will be added to the output. The value for each element is being
      determined by the function given by the argument getStructDataValues.
      None means: no extra fields.
    @param getStructDataValues: A function that takes three arguments for
      each element: element number, centroid and elset name and returns
      a list of values. This list must have the same length as the argument
      structDataFieldNames. Each value of this list defines the value for
      the particular element of one of the time constant fields being added
      to the output.
      None means: no extra fields.
      WE NEED AN EXAMPLE!
    """

    try:
        from abaqusConstants import CENTROID
    except ImportError:
        raise ImportError("Could not import abaqusConstants.CENTROID."
                          " Please run with <abaqus python>.")

    from bae.odb_02 import OdbReader, getOdbFieldClass, odbAccess
    msg("Opening odb %s" % odbPath)
    odb = OdbReader(odbPath)


    #-- create filterElset from elset name and bounding box
    if boundingBox is None:
        extraText = ""
    else:
        extraText = (
            "\n... and filtering elements in bounding box %s" % boundingBox)
    msg("Calculating element centroids%s." % extraText)

    instance = odb.getOdbPartInstance()
    a = odb.odb.rootAssembly

    if isinstance(filterElsetName, basestring):
        elsetNameList = [filterElsetName,]
    else:
        elsetNameList = filterElsetName

    elemIdsInBox = []
    elemCentroids = []
    elsetNameElemNbList = []
    for elsetName in elsetNameList:

        # loop over elements in elsetName...
        # ...calculate centroid and and find those in bounding box
        elementArray = instance.elementSets[elsetName].elements
        cnt = 0
        for elem in elementArray:
            connectivity = elem.connectivity
            # calculate centroid of element
            nodeSet = a.NodeSetFromNodeLabels(
                name='nodeSetForElem_%s' % elem.label,
                nodeLabels=(('PART-1-1',connectivity),))
            xC=[0,0,0]
            for node in nodeSet.nodes[0]:
                xC[0] += node.coordinates[0]
                xC[1] += node.coordinates[1]
                xC[2] += node.coordinates[2]
            xC[0] = xC[0]/len(connectivity)
            xC[1] = xC[1]/len(connectivity)
            xC[2] = xC[2]/len(connectivity)

            if boundingBox is None or (
                    boundingBox[0][0]<xC[0]<boundingBox[1][0]
                    and boundingBox[0][1]<xC[1]<boundingBox[1][1]
                    and boundingBox[0][2]<xC[2]<boundingBox[1][2]):
                elemIdsInBox.append(elem.label)
                elemCentroids.append(xC)
                cnt += 1

        elsetNameElemNbList.append( (elsetName, cnt) )
        msg("Found %d elements in the specified elset %s and %d of them in the"
            " bounding box."
            % (len(elementArray), elsetName, cnt))

    # create filter elset:
    # elems in box and according to parameter filterElsetName
    # filterElsetOdbName = 'ELEMS_IN_BOX'
    # filterElset = instance.ElementSetFromElementLabels(
    #     name=filterElsetOdbName,elementLabels=elemIdsInBox)
    for i in count(1):
        filterElsetOdbName = 'ELEMS_IN_BOX_%02d' % i
        try:
            filterElset = instance.ElementSetFromElementLabels(
                name=filterElsetOdbName,elementLabels=elemIdsInBox)
        except odbAccess.OdbError, exc:
            if exc.args[0].startswith("Duplicate set or surface name"):
                msg("WARNING: Elset %s already exists. Trying different name."
                    % filterElsetOdbName)
                continue
            else:
                raise
        else:
            break
    msg("Created filter elset %s with %d elements all together."
        % (filterElsetOdbName, len(elemIdsInBox)))


    #-- preparing output (directory, file name template and container-object)

    # convert stepFrameToName function:
    # The frameIdToName argument to vtk_02.FramesFieldsOnUnstructuredPoints.
    # ....initOutput --if its a function-- takes a single frameId argument:
    # a (stepNb, frameNb) - tuple.
    # On the other hand all bae.utils/odbToPointsData - methods
    # expect for their stepFrameToName argument a function taking two arguments:
    # stepNb, frameNb
    if callable(stepFrameToName):
        # string don't need conversion, but functions do
        stepFrameToNameTwoArgs = stepFrameToName
        def stepFrameToName(frameId):
            return stepFrameToNameTwoArgs(*frameId)

    # check output directory
    if outputDirName and not os.path.isdir(outputDirName):
        os.makedirs(outputDirName)

    # generate the output file name template
    if outputType.startswith("OneCsv"):
        # all in one file
        # outputFileNameTemplate is the file name (not a template)
        outputFileNameTemplate = os.path.join(
            outputDirName, fileNamePrefix+".csv")
    else:
        # one file per frame
        outputFileNameTemplate = os.path.join(
            outputDirName, fileNamePrefix+"_%s.csv")
        if addFrameNames and stepFrameToName is None:
            def stepFrameToName((stepNb, frameNb),
                                frameNames=odb.postData["frameNames"]):
                # have the frameNames object evaluated at function definition
                # time and stored in a default argument
                return "F%d%03d_%s" % (
                    stepNb, frameNb, frameNames[stepNb][frameNb])
    if not(stepFrameToName):
        stepFrameToName = "F%d%03d"

    output = FramesFieldsOnUnstructuredPoints(points=elemCentroids)
    output.initOutput(
        outputFileNameTemplate, outputType, frameIdToName=stepFrameToName)


    #-- some more preparation

    # prepare structDataFields
    if structDataFieldNames and getStructDataValues:

        msg("Preparing structDataFields %s" % ", ".join(structDataFieldNames))
        elsetNamesIter = chain(*[repeat(x, n) for x,n in elsetNameElemNbList])
        structDataFields = [
            # arguments: element number, centroid and elset name
            getStructDataValues(el, centroid, elsetName)
            for el, centroid, elsetName in izip(
                elemIdsInBox, elemCentroids, elsetNamesIter)]
        structDataFields = zip(*structDataFields)

        output.updateFields(*(
            createFieldObject(name, "point", "scalar", initArgs=[data,])
            for name, data in izip(structDataFieldNames, structDataFields)))
        output.storeConstantData()

    # preparing/modifying frame list in case of fields to maximize
    if fieldNamesMaximize:
        # extract first entry from frameList as frameStartMaximize
        # (removing this first entry from frameList)
        frameStartMaximize = frameList[0]
        frameList = frameList[1:]

        # create frame list for fields in fieldNamesMaximize
        if frameList[0][0] != frameList[-1][0]:
            msg("WARNING!!! Multi-step results not working yet.")
        stepNb = frameList[0][0]
        frameListAll = [
            (stepNb, frNb)
            for frNb in range(frameStartMaximize[1], frameList[-1][1]+1) ]

    else:
        frameListAll = frameList
    lastFrame = frameList[-1]
    frameList = set(frameList)

    # initialize fieldPersistentData, store field classes
    fieldPersistentData = dict()
    lastFrame = odb.getFrame2(*lastFrame)
    for fieldName in fieldNames:
        c = Container()
        fieldPersistentData[fieldName] = c

        # class to read data from the odb
        c.OdbFieldClass = getOdbFieldClass(fieldName, lastFrame)

        # possible name changes for output
        # remove SDV_ from fieldName from now on
        if fieldName.startswith("SDV_"):
            fieldName = fieldName[4:]
        if fieldName == "S_MIN_PRINCIPAL":
            fieldName = "S"

        # class to store result data
        c.FieldClass = createFieldClass(
            fieldName, "point", c.OdbFieldClass.dataType)


    # function to get the data from the odb
    def getOdbField(odbFrame, OdbFieldClass, FieldClass,
                    filterElset, elemIdsInBox):
        """
        @param FieldClass: Field class for the resulting field
        """
        # get result Subset - CENTROID and COHSV_IN_BOX elset
        odbField = OdbFieldClass().odbFieldFromFrame(odbFrame)
        fieldSubset = OdbFieldClass()
        fieldSubset.position="element"
        fieldSubset.readFromOdbField(
            odbField.getSubset(position=CENTROID)
            .getSubset(region=filterElset))

        # sort according to elemIdsInBox, store in list
        field = FieldClass(
            fieldSubset[elemNb] for elemNb in elemIdsInBox)
        return field

    #-- main loop over all frames
    for stepNb, frameNb, odbFrame in odb.getFrameIterator3(frameListAll):

        # main frames for output: read all variables
        if (stepNb, frameNb) in frameList:
            msg("Now processing field output for frame F%d%03d"
                % (stepNb, frameNb))

            # loop over field variables and extract values to .csv
            for fieldName in fieldNames:

                persist = fieldPersistentData[fieldName]
                field = getOdbField(
                    odbFrame, persist.OdbFieldClass, persist.FieldClass,
                    filterElset, elemIdsInBox)

                try:
                    oldMax = persist.maxData
                except AttributeError:
                    pass
                else:
                    msg("Computing max value for %s..." % fieldName)
                    field = type(field)(
                        max(a,b) for a, b in izip(field, oldMax))
                    del persist.maxData

                msg("adding data for %s, %s, %s"
                    % (stepNb, frameNb, field.fieldName))
                output.updateFields(field)

            output.storeFrame((stepNb, frameNb))

        # intermediate frames for RER calculation only
        else:
            msg("Now processing frame F%d%03d, max field(s) only..."
                % (stepNb, frameNb))

            for fieldName in fieldNamesMaximize:

                persist = fieldPersistentData[fieldName]
                field = getOdbField(
                    odbFrame, persist.OdbFieldClass, persist.FieldClass,
                    filterElset, elemIdsInBox)

                try:
                    oldMax = persist.maxData
                except AttributeError:
                    persist.maxData = field
                else:
                    msg("Computing max value for %s..." % fieldName)
                    persist.maxData = type(field)(
                        max(a,b) for a, b in izip(field, oldMax))

    # end of the framelist-loop
    msg("Got data from %d output frames from the odb." % len(frameList))

    output.closeCsvOutput()
    msg("Wrote data to csv file(s).")
    odb.close()


#------------------------------------------------------------------------------


def odbCohsvCentrLocalCoordsToCsv(
        odbPath, outputDirName, outputFileName,
        boundingBox=None, filterElsetName='COHSV', fieldName="S", frame=(2,0),
        ):
    """Read local coordinate system vectors and element centroids from the
    specified odb and write the data to csv files.

    Usage:
     >>> odbPath = config.odbPath
     >>> outputDirName = "CSV/Cadia2016033_R01_COHSV"
     >>> outputFileName = "Cadia2016033_G01_COHSV_localCoords.csv"
     >>> fieldName = "E"
     >>> frame = (2, 0)
     >>>
     >>> odbElCentroidsToCsvPoints(
     >>>     odbPath, outputDirName, outputFileName,
     >>>     boundingBox=None, filterElsetName="COHSV",
     >>>     fieldName=fieldName, frame=frame,
     >>>     )

    @param odbPath: Complete path to the odb including .odb extension. Can be
      relative or absolute.
    @param outputDirName: Vtk files will end up in this folder. Can be
      a relative or an absolute path.
    @param outputFileName:
    @param boundingBox: L{BoundingBox<bae.misc_01.BoundingBox>}-object (or
      [[xmin,ymin,zmin],[xmax,ymax,zmax]]). Consider only elements whose
      centroid lies within this box. If None then all elements will be
      considered.
    @param filterElsetName: Name of the elset to be exported. Defaults to
      "COHSV". If it's a list then consider all elements from all elsets in
      this list.
    @param fieldName: Field name from which to get the local coordinate system.
      Must be a tensor field stored in a local coordinate system like "E" or
      "S".
    @param frame: (step number, frame number) tuple from which to get the local
      coordinate system. This frame must contain the specified field.
      Example: (2, 0)
    """

    try:
        from abaqusConstants import CENTROID
    except ImportError:
        raise ImportError("Could not import abaqusConstants.CENTROID."
                          " Please run with <abaqus python>.")

    from bae.odb_02 import OdbReader, getOdbFieldClass, odbAccess
    msg("Opening odb %s" % odbPath)
    odb = OdbReader(odbPath)


    #-- create filterElset from elset name and bounding box
    if boundingBox is None:
        extraText = ""
    else:
        extraText = (
            "\n... and filtering elements in bounding box %s" % boundingBox)
    msg("Calculating element centroids%s." % extraText)

    instance = odb.getOdbPartInstance()
    a = odb.odb.rootAssembly

    if isinstance(filterElsetName, basestring):
        elsetNameList = [filterElsetName,]
    else:
        elsetNameList = filterElsetName

    elemIdsInBox = []
    elemCentroids = []
    elsetNameElemNbList = []
    for elsetName in elsetNameList:

        # loop over elements in elsetName...
        # ...calculate centroid and and find those in bounding box
        elementArray = instance.elementSets[elsetName].elements
        cnt = 0
        for elem in elementArray:
            connectivity = elem.connectivity
            # calculate centroid of element
            nodeSet = a.NodeSetFromNodeLabels(
                name='nodeSetForElem_%s' % elem.label,
                nodeLabels=(('PART-1-1',connectivity),))
            xC=[0,0,0]
            for node in nodeSet.nodes[0]:
                xC[0] += node.coordinates[0]
                xC[1] += node.coordinates[1]
                xC[2] += node.coordinates[2]
            xC[0] = xC[0]/len(connectivity)
            xC[1] = xC[1]/len(connectivity)
            xC[2] = xC[2]/len(connectivity)

            if boundingBox is None or (
                    boundingBox[0][0]<xC[0]<boundingBox[1][0]
                    and boundingBox[0][1]<xC[1]<boundingBox[1][1]
                    and boundingBox[0][2]<xC[2]<boundingBox[1][2]):
                elemIdsInBox.append(elem.label)
                elemCentroids.append(xC)
                cnt += 1

        elsetNameElemNbList.append( (elsetName, cnt) )
        msg("Found %d elements in the specified elset %s and %d of them in the"
            " bounding box."
            % (len(elementArray), elsetName, cnt))

    # create filter elset:
    # elems in box and according to parameter filterElsetName
    # filterElsetOdbName = 'ELEMS_IN_BOX'
    # filterElset = instance.ElementSetFromElementLabels(
    #     name=filterElsetOdbName,elementLabels=elemIdsInBox)
    for i in count(1):
        filterElsetOdbName = 'ELEMS_IN_BOX_%02d' % i
        try:
            filterElset = instance.ElementSetFromElementLabels(
                name=filterElsetOdbName,elementLabels=elemIdsInBox)
        except odbAccess.OdbError, exc:
            if exc.args[0].startswith("Duplicate set or surface name"):
                msg("WARNING: Elset %s already exists. Trying different name."
                    % filterElsetOdbName)
                continue
            else:
                raise
        else:
            break
    msg("Created filter elset %s with %d elements all together."
        % (filterElsetOdbName, len(elemIdsInBox)))


    #-- preparing output (directory, file name template and container-object)

    # check output directory
    if outputDirName and not os.path.isdir(outputDirName):
        os.makedirs(outputDirName)

    #-- get the coordinate system from the odb
    msg("Reading local coordinate systems.")
    odbFrame = odb.getFrame2(*frame)
    odbField = odbFrame.fieldOutputs[fieldName].getSubset(region=filterElset)
    fieldLocCoordSys = dict(
        (odbFldVal.elementLabel, odbFldVal.localCoordSystem)
        for odbFldVal in odbField.values)

    # write output
    outputPath = os.path.join(outputDirName, outputFileName)
    msg("Now writing results to %s" % outputPath)
    output = csv.writer(open(outputPath, "wb"))
    output.writerow(["x","y","z",
                     "t1_x","t1_y","t1_z",
                     "t2_x","t2_y","t2_z",
                     "n_x","n_y","n_z"])
    for elemNb, centroid in izip(elemIdsInBox, elemCentroids):
        row = list(centroid)
        for x in fieldLocCoordSys[elemNb]:
            row.extend(x)
        output.writerow(row)
    del output
    msg("Wrote results to %s" % outputPath)
    odb.close()
