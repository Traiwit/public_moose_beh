#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""This module contains functions to import and export field data to and from
various file formats. They are supplied to the main
L{FieldsCollection<bae.field_02.FieldsCollection>} class
through the Mixin-class L{FieldsCollectionReaderWriter} defined in this module.

The actual algorithms are spread over format-specific sub-modules.
"""

# Note: version is stored in bae/field_02/__init__.py

import numpy as np
import os
import glob
from collections import OrderedDict

from .vtk import read_buffer as readVtk, write as writeVtk
from .iv import IvWriter
from ..topo.structuredgrid import StructuredGrid    ## check: still needed?
from ..topo.unstructuredmesh.meshbase import \
    HomoMesh, HeteroMesh, \
    MeshNodes, MeshElems, MeshGaussPts
from bae.log_01 import msg


__todo__ = """
- vtk interface: fc = FieldsCollection(...), FC = FieldsCollection
  * fc = FC.fromVtk(...) can take...
    . a single vtk file path
    . a filename-glob for a time series
    . an open stream like a StringIO-object...
  * FC.fromVtk() calls specific subroutines that might call each other again:
    . fromVtkFiles() reads many files into a time series
    . fromVtkStream() (or -Buffer?) reads one file which is already open and
      can also be a stream like a StringIO-object
  * FC.toVtk() similar...
    . can take single path or an open stream if it has only one time point,
      otherwise --if it's got more than one time point-- it must stop
    . the path can contain a wildcard that is to be filled with an item of
      fc.times. Which item to take can be specified by an additional argument
      which defaults to "timeName".


- pvf / csv read write


        ... for many files, time points
        ===============================
         >>> from field_02 import FieldsCollection
         >>> from field_02.readerwriter import globTimeSeries
         >>> flds = FieldsCollection()
         >>> for i, time, fname in enumerate(globTimeSeries("mydata_*.vtk")):
         >>>     flds[:,i] = FieldsCollection.fromVtk(fname)  # or flds[:,:,i] ?
         >>>     flds .... nee, wir muessen vorher wissen, wie viele Dateien wir haben...
         >>>     
         >>>     

         >>> 
         >>> 
         >>> 

- include something like...
  bae.utils.odbToPointsData_02.Converter_CombinedVtk.readVTKfromFileOrZipArch()
  and possibly also merge with globTimeSeries.
  ... Don't know if this module is the right place...
"""

vtkModuleInstallInstructions = """
Install instructions for vtk module
===================================

linux/pip: "pip install --update pip
            python -m pip install --upgrade --user vtk"
Note: the pip-USERPATH on nukus should be set to /usr/local/pip for all
users. So you may run the previous lines:
 - as su (sudo ...)
 - and change owner (should be tobias) and rights recursivelly:
     sudo chown -R tobias /usr/local/pip
     sudo chmod -r g+rx /usr/local/pip

conda (not suggested for our Linux servers): 
    "conda install -c anaconda vtk"
or for conda python 2.7
    "conda install -c anaconda vtk=8.1"

(you may need to run anaconda-console/prompt as admin)

windows/raw-python: Not sure what to do here. Fred tried for one full day and
finally used anaconda. Any further suggestions/experiences welcome.
just a hint: https://vtk.org/download/#latestcand
"""


#{ Errors
class ReadError(Exception):
    r"""Can be raised by FieldsOnStructPointGrid.fromVtk() if the vtk file
    to be read contains unexpected data or structure.
    """
    pass


class WriteError(Exception):
    r"""Can be raised by FieldsOnStructPointGrid.toVtk().
    """
    pass
#} end of Errors


def globTimeSeries(path, filePattern="*"):
    """Service function to parse path expressions that can refer to one or
    many files. Expands path name pattern like 'mydata_*.vtk' and finds files
    in a specified folder.

    Examples
    ========

    file name pattern
    -----------------
     >>> fnames, ids = globTimeSeries("mydata_*.vtk")
     >>> print fnames
     ["mydata_Y2010.vtk", "mydata_Y2013.vtk", "mydata_Y2017.vtk"]
     >>> print ids
     ["Y2010", "Y2013", "Y2017"]

    directory name
    --------------
     >>> fnames, ids = globTimeSeries("./data", filePattern="*.vtk")
     >>> print fnames
     ["mydata_Y2010.vtk", "mydata_Y2013.vtk", "mydata_Y2017.vtk"]
     >>> print ids
     ["mydata_Y2010", "mydata_Y2013", "mydata_Y2017"]

    single file
    -----------
     >>> fnames, ids = globTimeSeries("./data/mydata_Y2017.vtk")
     >>> print fnames
     ["mydata_Y2017.vtk"]
     >>> print ids
     []

    @param path: Can be a single filename, a directory name or a file name
    pattern containing an asterisk '*'.

    @param filePattern: Optional, if path refers to a directory then
    filePattern acts as a filter. E.g. filePattern="*.vtk" would restrict
    the result to vtk files.

    @returns: a tuple: list of expanded paths and list of identifier strings.
    The identifier strings contain -for each file- the wildcard part of the
    filename, i.e. the part that is replaced by the asterisk ('*') and can be
    thought as identifier for the different files. If path is a directory then
    this refers to the wildcard part of filePattern.

    If path doesn't contain a wildcard then the second item will be an empty
    list.
    """

    if isinstance(path, basestring):
        # single stringinput --> parse inputFileNames
        if '*' in path:
            # is a pattern --> get all files matching pattern
            fileNames = glob.glob(path)
            pos1 = path.index("*")
            pos2 = pos1 - len(path) + 1
            ids = [fn[pos1:pos2] for fn in fileNames]

        elif os.path.isdir(path):
            # is a directory --> get all vtk-files in directory
            fileNames = glob.glob(os.path.join(path, filePattern))
            pos1 = filePattern.index("*")
            pos2 = pos1 - len(filePattern) + 1
            ids = [os.basename(fn)[pos1:pos2] for fn in fileNames]
        else:
            # is a single filename
            fileNames = [path,]
            ids = []
    else:
        fileNames = path
        ids = []

    # check if files exist
    notFoundFiles = [f for f in fileNames if not os.path.isfile(f)]
    if notFoundFiles:
        raise ReadError(
            '%d file(s) not found: %s'
            % (len(notFoundFiles), ', '.join(notFoundFiles)))
    return fileNames, ids


class FieldsCollectionReaderWriter(object):
    """Mixin for FieldsCollection supplying methods to import and export
    field data to various file formats.
    """

    ######################################################
    #{ file paths
    def getPathListForExport(self, outputPath, timesField="timeName"):
        """Service function for L{toVtk} and others: Determines the list of
        actual file names (full path) for export e.g. as vtk file(s).

        @param outputPath: Can be a file name for a single file or a
          file name pattern expression (containing an asterisk '*' or a '%s'.
          If the placeholder exists it will be replaced by the contents of
          the times field identified by the timesField argument.

        @param timesField: key in self.times supplying the variable part in the
          outputPath.
        """
        # nb of files to be written and file names -> pathList
        nbTimes = self.getShape()[1]
        if "*" in outputPath:
            try:
                dateList = self.times[timesField]
            except KeyError:
                dateList = ["%03d"%i for i in range(nbTimes)]
            pathList = [outputPath.replace("*", x) for x in dateList]
        elif nbTimes==1:
            pathList = [outputPath]
        else:
            # create
            prefix, ext = os.path.splitext(outputPath)
            try:
                dateList = self.times[timesField]
            except KeyError:
                dateList = ["%03d" % i for i in range(nbTimes)]
            pathList = ["%s_%s%s" % (prefix, x, ext) for x in dateList]
        return pathList
    #} end of file paths

    ######################################################
    #{ read and write VTK
    @classmethod
    def fromVtk(cls, inputPath):
        r"""Creates Fields-Container from a single vtk file or a time series
        of vtk files.

        Usage:
        ======
         >>> from field_02 import FieldsCollection
         >>> flds = FieldsCollection.fromVtk("mydata.vtk")

        @param inputPath: Can be a single filename, a directory name or a file
           name pattern containing an asterisk '*'.

        @returns: object of class
           L{FieldsCollection<bae.field_02.FieldsCollection>} (or derived
           class). If inputPath is a file name pattern or a directory then
           self.times["timeName"] will contain time identifiers: In case of a
           directory this will be the file name excluding the extension ".vtk".
           In case of a file name pattern it will contain the actual substring
           replacing the "*" asterisk in the file name.
        """

        fileList, timeName = globTimeSeries(inputPath, filePattern="*.vtk")
        nbTimes = len(fileList)
        descriptions = list()  # collect descriptions in list first

        for ii, fileName in enumerate(fileList):

            #--- read single vtk-file
            inputFile = open(fileName, "rb")
            vtkData, mesh = readVtk(inputFile)

            #--- process the topo
            if vtkData.point_data and vtkData.cell_data_raw:
                raise ReadError(
                    "Input file %s contains point data and cell data at"
                    " the same time. This requires a function (yet to be"
                    " implemented) that returns two or three fieldsCollection"
                    " items."% (fileName, inputPath))

            if ii==0:
                # only process the topo for the first vtk file in the series
                if isinstance(mesh, (HomoMesh, HeteroMesh)):
                    # mesh is an unstructured mesh
                    if vtkData.point_data:
                        topo = MeshNodes(mesh)
                    elif vtkData.cell_data_raw:
                        topo = MeshElems(mesh)
                    else:
                        raise NotImplementedError(
                            "Currently only cell and point data are already"
                            " implemented for unstructured grids.")
                elif isinstance(mesh, StructuredGrid):
                    topo = mesh
            else:
                # compare data from current file to stored topo
                if (vtkData.point_data
                    and not isinstance(topo, (MeshNodes, StructuredGrid))):
                    raise ReadError(
                        "Input file %s in the file series %s has a point data"
                        " as opposed to earlier file(s)."
                        % (fileName, inputPath))

                if vtkData.cell_data_raw and not isinstance(topo, MeshElems):
                    raise ReadError(
                        "Input file %s in the file series %s has a point data"
                        " as opposed to earlier file(s)."
                        % (fileName, inputPath))

                if topo.mesh != mesh:
                    raise ReadError(
                        "Input file %s in the file series %s has a different"
                        " mesh" % (fileName, inputPath))

            #--- data field from point or cell data?
            if vtkData.point_data:
                newFlds = vtkData.point_data
            elif vtkData.cell_data_raw:
                newFlds = vtkData.cell_data_raw

            #--- in first iteration:
            # initialize FieldsCollection object as result
            if ii==0:
                times = dict()
                if timeName:
                    times["timeName"] = np.array(timeName)
                result = cls(topo=topo, times=times)
                for name, newFld in newFlds.iteritems():
                    shape = (len(topo), nbTimes) + newFld.shape[1:]
                    result[name] = np.empty(shape=shape, dtype=newFld.dtype)

            #--- store the field data (copy)
            for name, newFld in newFlds.iteritems():
                result[name][:,ii,...] = newFld

            descriptions.append(vtkData.title)

        # Only now store descriptions. So the ndarray will hold all strings with
        # minimal memory.
        result.times["description"] = np.array(descriptions)

        return result


    @classmethod
    def fromVtkBuffer(cls, inputFile):
        r"""Creates Fields-Container from an open vtk file object.

        Usage:
        ======
         >>> from field_02 import FieldsCollection
         >>> flds = FieldsCollection.fromVtkBuffer(open("mydata.vtk", "rb"))

        @param inputFile: An open file-like object.

        @returns: object of class
           L{FieldsCollection<bae.field_02.FieldsCollection>} (or derived
           class).
        """

        vtkData, topo = readVtk(inputFile)
        times = dict(description=np.array([vtkData.title,]))

        
    @classmethod
    def fromVtk_OLD(cls, inputPath):
        r"""Old version, DEPRECATED! For reference only.

        Creates Fields-Container from a single vtk file.
        """

        ###GP: is it possible to read from input streams like zipfile-objects?
        ###GP: check out vtkDataSetReader.SetInputArray. takes  vtkCharArray
        ###GP: which might take a string as initializer.

        # postponed import of vtk module because it takes ~1sec.
        try:
            import vtk
            from vtk.util.numpy_support import vtk_to_numpy
            # from vtk.numpy_interface import dataset_adapter as dsa
        except ImportError:
            raise ImportError("Could not import vtk module.\n\n%s"
                              % vtkModuleInstallInstructions)

        fileList, timeName = globTimeSeries(inputPath, filePattern="*.vtk")
        nbTimes = len(fileList)
        fields = OrderedDict()
        descriptions = list()

        for ii, fileName in enumerate(fileList):

            # read single vtk-file
            reader = vtk.vtkDataSetReader()
            reader.SetFileName(fileName)

            reader.ReadAllScalarsOn()
            reader.ReadAllVectorsOn()
            reader.ReadAllTensorsOn()

            # actually read the file
            reader.Update()
            output = reader.GetOutput()
            
            ###GP: here we need to check and differentiate what kind of topo we
            ###GP: actually have...

            ptData = output.GetPointData()
            cellData = output.GetCellData()
            fldData = output.GetFieldData()

            # check what's in the file
            if cellData.GetNumberOfArrays()>0:
                msg("WARNING! Found CellData. This will be ignored!")
            if fldData.GetNumberOfArrays()>0:
                msg("WARNING! Found FieldData. This will be ignored.")
            if ptData.GetNumberOfArrays()<=0:
                msg("WARNING! No PointData found in %s." % fileName)

            descriptions.append(reader.GetHeader())

            # create StructuredGrid-object
            gridParam = dict(
                gridPtNb=output.GetDimensions(),
                spacing=output.GetSpacing(),
                origin=output.GetOrigin(),
                )
            grid = StructuredGrid(**gridParam)

            if ii == 0:
                # store grid from first vtk file
                refGrid = grid
            else:
                # compare in subsequent iterations
                if not grid == refGrid:
                    raise ValueError(
                        "In time series of vtk files nb %d (%s) has different"
                        "grid parameters." % (ii, fileName))

            # collect pointdata fields
            for fieldCnt in range(ptData.GetNumberOfArrays()):
                fieldName = ptData.GetArrayName(fieldCnt)
                vtkArray = ptData.GetArray(fieldCnt)
                arr = vtk_to_numpy(vtkArray)  # or GetArray(fldName)
                if vtkArray.GetNumberOfComponents()==9:
                    arr = arr.reshape((-1, 3, 3))
                if ii == 0:
                    # init emtpy field in first iteration
                    newArr = np.zeros((nbTimes,) + arr.shape)
                    newArr = np.moveaxis(newArr, 0, 1)
                    fields[fieldName] = newArr

                # assign data to array
                fields[fieldName][:,ii] = arr

        times = dict()
        times["description"] = np.array(descriptions)
        if timeName:
            times["timeName"] = np.array(timeName)

        # first create new fieldscollection then append fields in correct order
        result = cls(topo=grid, times=times)
        for fieldName, arr in fields.iteritems():
            result[fieldName] = arr
        return result

    def toVtk(self, outputPath,
              timesField="timeName",
              description=None,
              descrField="description",
              outputFormat='binary', writeSinglePrecision=True):
        r"""Write stored data to vtk file of the given name.

        @param outputPath: Can be a file name for a single vtk file or a
          file name pattern expression (containing an asterisk '*' or a '%s'.
          If the placeholder exists it will be replaced by the contents of
          the times field identified by the timesField argument.

        @param timesField: key in self.times supplying the variable part in the
          outputPath.

        @param description: String for the second line of the vtk data file
          (title field). If not specified then a description times-field is
          tried (see descrField argument). If that doesn't apply then the
          description defaults to "Field data".
        @param descrField: field name of a field in self.times that holds
          strings to be used as description/vtk-header for the particular
          time.

        @param outputFormat: One of "binary", "ascii", "ascii_1row".
          "ascii_1row" means scalar fields are written as one row with space as
          delimiters instead of newline characters. (This is the pyvtk default
          ascii implementation but Glenn for example needs "ascii" with newline
          as delimiter.)
        @param writeSinglePrecision: If True then all binary data will be
          exported as single precision
        """

        pathList = self.getPathListForExport(outputPath, timesField)

        # loop over times / individual vtk files
        for iTime, fileName in enumerate(pathList):

                # set header / description
                if description is not None:
                    currentDescr = description
                elif descrField in self.times:
                    currentDescr = self.times[descrField][iTime]
                else:
                    currentDescr = "Field data"

                if writeSinglePrecision:
                    # convert to single precision
                    fields = dict((k, v.astype(dtype=np.float32))
                                  for k, v in self.fields.iteritems())
                else:
                    fields = self.fields

                # write a vtk file
                writeVtk(fileName, self.topo, fields, iTime,
                         description=currentDescr)
        return


    def toVtk_OLD(self, outputPath,
              timesField="timeName",
              description=None,
              descrField="description",
              outputFormat='binary', writeSinglePrecision=True):
        r"""Old version, DEPRECATED! For reference only.

        Write stored data to vtk file of the given name.

        @param outputPath: Can be a file name for a single vtk file or a
          file name pattern expression (containing an asterisk '*' or a '%s'.
          If the placeholder exists it will be replaced by the contents of
          the times field identified by the timesField argument.

        @param timesField: key in self.times supplying the variable part in the
          outputPath. times-field that contains

        @param description: String for the second line of the vtk data file
          (title field). If not specified then a description times-field is
          tried (see descrField argument). If that doesn't apply then the
          description defaults to "Field data".
        @param descrField: field name of a field in self.times that holds
          strings to be used as description/vtk-header for the particular
          time.

        @param outputFormat: One of "binary", "ascii", "ascii_1row".
          "ascii_1row" means scalar fields are written as one row with space as
          delimiters instead of newline characters. (This is the pyvtk default
          ascii implementation but Glenn for example needs "ascii" with newline
          as delimiter.)
        @param writeSinglePrecision: If True then all binary data will be
          exported as single precision
        """

        # postponed import of vtk module because it takes ~1sec.
        try:
            import vtk
            from vtk.util.numpy_support import numpy_to_vtk
        except ImportError:
            raise ImportError("Could not import vtk module.\n\n%s"
                              % vtkModuleInstallInstructions)

        # nb of files to be written and file names -> pathList
        nbTimes = self.getShape()[1]
        if "*" in outputPath:
            try:
                pathList = [
                    outputPath.replace("*", x) for x in self.times[timesField]]
            except KeyError:
                pathList = [
                    outputPath.replace("*", "%03d"%i)
                    for i in range(nbTimes)]
        elif nbTimes==1:
            pathList = [outputPath]
        else:
            # create
            prefix = outputPath.rsplit(".vtk", 1)[0]
            try:
                pathList = [
                    "%s_%s.vtk" % (prefix, x) for x in self.times[timesField]]
            except KeyError:
                pathList = [
                    "%s_%03d.vtk" % (prefix, i)
                    for i in range(nbTimes)]

        # loop over times / individual vtk files
        for iTime in range(nbTimes):

            # depending on type of topo: init vtk and store topo
            if isinstance(self.topo, StructuredGrid):
                out = vtk.vtkStructuredPoints()
                out.SetDimensions(*self.topo.gridPtNb)
                out.SetSpacing(*self.topo.spacing)
                out.SetOrigin(*self.topo.origin)

                # vtk data entity to store fields at
                # note: if we ever want to store 'cell' or 'field' data then
                # take out.GetCellData() resp. out.GetFieldData()
                data = out.GetPointData()

            else:
                raise NotImplementedError(
                    "Writing data on topo of type %s to vtk not implemented,"
                    " yet. Currently only structured grid supported."
                    % type(self.topo))


            for fieldName in self:
                arr = self.fields[fieldName]

                # set the type of field
                dtype = arr.dtype
                if dtype.kind in 'f':
                    vtkType = vtk.VTK_FLOAT

                elif dtype.kind in 'ui':
                    vtkType = vtk.VTK_INT

                else:
                    raise ValueError('dtype %s is not supported' % str(dtype))

                # depending on data dimensions (scalar, vector, tensor)...
                # - reshape data
                # - choose storage method
                if len(arr.shape)==2:
                    # scalar

                    # initialize vtkArray, transfer data and set field name
                    valsVtk = numpy_to_vtk(
                        # num_array=arr[:,iTime],
                        num_array=np.ascontiguousarray(arr[:,iTime]),
                        deep=False,  # the numpy array persists until writing
                        array_type=vtkType)
                    valsVtk.SetName(fieldName)

                    data.AddArray(valsVtk)

                elif len(arr.shape)==3:
                    # vector
                    if not arr.shape[2] == 3:
                        raise ValueError(
                            "Expects last axis of array to have length 3."
                            " Yours has %d. This does not work."
                            % arr.shape[2])

                    # initialize vtkArray, transfer data and set field name
                    valsVtk = numpy_to_vtk(
                        # num_array=arr[:,iTime],
                        num_array=np.ascontiguousarray(arr[:,iTime]),
                        deep=False,
                        array_type=vtkType)
                    valsVtk.SetName(fieldName)

                    data.SetVectors(valsVtk)

                elif len(arr.shape)==4:
                    # tensor
                    if arr.shape[2:4] != (3,3):
                        raise ValueError(
                            "Expects last two axes of array to have length 3.",
                            " Your shape: %s. This does not work."
                            % str(arr.shape))

                    # initialize vtkArray, transfer data and set field name
                    valsVtk = numpy_to_vtk(
                        # To avoid -by all means- an error concering
                        # discontiguous data we follow the following
                        # --supposedly less efficient-- procedure:
                        # 1. reshape so that we have the components and space
                        #    dim in correct order
                        # 2. make it contiguous in memory because vtk needs it
                        # 3. store in vtk structure with "deep=True" because
                        #    the array created "on the fly" might not stay
                        #    alive long enough
                        num_array=np.ascontiguousarray(
                            arr[:,iTime].reshape(-1,9)),
                        deep=True,
                        # supposedly more efficient method below, might be less
                        # secure:
                        # num_array=arr[:,iTime].reshape(-1,9),
                        # deep=False,
                        array_type=vtkType)
                    valsVtk.SetName(fieldName)

                    data.SetTensors(valsVtk)

                # for future extensions: consider these methods...
                #vectors = vtkDoubleArray()
                #vectors.SetNumberOfComponents(3)
                #for vec in map(tuple(vecField)):
                #    vectors.insertNextTuple(*vec)
                #out.GetPointData().SetVectors(vectors)
                #tensors = vtkDoubleArray()
                #tensors.SetNumberOfTuples(2);
                #tensors.SetNumberOfComponents(9);
                #tensors.InsertTuple9(*...)

            # init new vtk file
            fileName = pathList[iTime]
            writer = vtk.vtkDataSetWriter()
            writer.SetFileName(fileName)

            # set header / description
            if description is not None:
                writer.SetHeader(description)
            elif descrField in self.times:
                writer.SetHeader(self.times[descrField][iTime])
            else:
                writer.SetHeader("Field data")

            if outputFormat=='binary':
                writer.SetFileTypeToBinary()

            writer.SetInputData(out)
            err = writer.Write()
            if not err:
                raise WriteError("Writing %s failed." % fileName)

            # modify VTK version to 2.0 (for Voxler)
            ff = open(fileName, "r+b")
            ff.seek(23)
            ff.write("2.0")
            ff.close()

    #} end of read and write VTK

    ####################################################
    #{ write OpenInventor (.iv) files
    def exportColourCodesAsIv(
            self, outputPath, colours, categoryField=None,
            timesField="timeName", outputFormat="ascii"):
        r"""Create iv files with surfaces coloured based on an element based
        field of integer colour indexes and a corresponding colours array.
        self.topo must be a L{MeshElems} object based on a L{TriMesh}.

        Example:
         >>> # five colours ...
         >>> colours = np.array([
         >>>     # rgb, alpha  (alpha: 0=opaque)
         >>>     [0.0, 0.0, 1.0, 0.0],  # blue
         >>>     [0.0, 1.0, 0.9, 0.0],  # light blue
         >>>     [0.1, 1.0, 0.0, 0.0],  # green
         >>>     [0.9, 1.0, 0.0, 0.0],  # yellow
         >>>     [1.0, 0.0, 0.0, 0.0],  # red
         >>> ])
         >>> # ... for five categories separated by following four values
         >>> bounds = np.array([0.5, 1, 2, 5])
         >>> inpField = FieldsCollection.fromVtk(inputPath)
         >>> category = np.searchsorted(bounds, inpField.array[:,0])
         >>> inpField["category"] = category[:,np.newaxis]
         >>> inpField.exportColourCodesAsIv(
         >>>     outputPath, colourCodes, categoryField="category")

        @param outputPath: Can be a file name for a single iv file or a
          file name pattern expression (containing an asterisk '*' or a '%s').
          If the placeholder exists it will be replaced by the contents of
          the times field identified by the timesField argument.

        @param categoryField: field name (key in self.fields) for a field of
          integers to be used as keys in the colours array to assign
          colours to the triangle faces.

        @param colours: array Nx3 or Nx4 for N categories (i.e. N > max
          value of the category field. For each category have a RGB-tuple,
          three floats between 0 and 1. Or an RGBA - tuple, four floats between
          0 and 1, the last for transparency (0 = opaque).

        @param timesField: key in self.times supplying the variable part in the
          outputPath.

        @param outputFormat: One of "binary", "ascii".
           "binary" not implemented yet. Default is subject to change from
           the current "ascii" to "binary" as soon as it's implemented.
        """

        elType = self.topo.mesh.elType  # abbreviation
        if (self.topo.position != "element" and elType.shape=="TRI"):
            raise ValueError(
                "exportColourCodesAsIv only works for element based values on"
                " a triangle mesh. We have self.topo.position='%s' and element"
                " shape='%s'" % (self.topo.position, elType.shape))

        if outputFormat != "ascii":
            raise NotImplementedError("Output format '%s' not implemented yet."
                                      % outputFormat)

        # select and check the category field
        if categoryField:
            catFldAllTimes = self[categoryField]
        else:
            catFldAllTimes = self.array

        if len(catFldAllTimes.shape)!=2:
            raise ValueError(
                "The category field '%s' must be an integer scalar. Its shape"
                " is expected to have length two. Instead the shape is: %s"
                % (categoryField, catFldAllTimes.shape))

        if np.min(catFldAllTimes)<0 or np.max(catFldAllTimes)>=len(colours):
            raise ValueError(
                "We have colour codes between %s and %s in the category field"
                " '%s' and the colour codes array is of length %d. This"
                " doesn't fit well."
                % (np.min(catFldAllTimes), np.max(catFldAllTimes),
                   categoryField, len(colours)))

        # prepare colour and transparency
        if colours.shape[1]==3:
            # if no transparency given then append a zero transparency value
            ccnew = np.zeros((len(colours.shape[0]), 4), colours.dtype)
            ccnew[:,:3] = colours
            colours = ccnew
            del ccnew
            
        pathList = self.getPathListForExport(outputPath, timesField)

        # loop over times / individual vtk files
        for iTime, fileName in enumerate(pathList):

            catFld = catFldAllTimes[:,iTime]
            writer = IvWriter(open(fileName, "w"))

            for catId in np.unique(catFld):
                elset = (catFld==catId)
                subTopo = self.topo.getSubTopo(elset)
                colour = colours[catId,:3]
                transparency = colours[catId,3]
                writer.writeTriMeshToIv(
                    subTopo.mesh, colour=colour, transparency=transparency)

            # close the iv file
            del writer
            
        return
    
    #} end of write OpenInventor (.iv) files
