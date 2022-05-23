"""Module to read and write vtk data sets and csv point data files.

For data on structured grids (e.g. vtk files for Voxler):
 - L{FieldsOnStructPointGrid}

For data on unstructured points clouds (csv point data files):
 - L{FieldsOnUnstructuredPoints}: read from a csv file (data for one time frame)
 - L{FramesFieldsOnUnstructuredPoints}: write to csv file/s

For data on FEM-meshes:
 - L{FieldsOnUnstructMesh}

@todo: build testmodule
 - especially: FramesFieldsOnUnstructuredPoints.closeCsvOutput() seems to mix
   up data columns. They should be in the order given by update but this does
   not seem to work with the script odbToCSVPoints.py...

@todo: Read vector and tensor fields with
FieldsOnUnstructuredPoints.fromCsvFile: rsplit fieldnames by "_" and group by
first part. If we have a second part for all required components then create
a vector or tensor. Make this behaviour optional.
"""

__version__ = "1.19"

_version_history_ = r"""
1.0 GP: first version
1.1 GP: added FramesFieldsOnUnstructuredPoints
1.2 GP: added writeSinglePrecision option to VtkData.toFile, default=True
1.3 GP: fixed bug in FramesFieldsOnUnstructuredPoints.closeCsvOutput()
        added field output of type str for FramesFieldsOnUnstructuredPoints
           (csv files)
1.4 GP: added outputType "OneCsvRowFrColPtField" to
        FramesFieldsOnUnstructuredPoints
1.5 GP added outputType "OneCsvPerField" to FramesFieldsOnUnstructuredPoints
1.6 GP: cleaned up API: new initOutput() method in
        FramesFieldsOnUnstructuredPoints
1.7 GP fixed: bae.vtk_02.FieldsOnStructPointGrid.fromVtk now reads point data
        attributes of type FIELD
1.8 GP added FieldsOnUnstructMesh
1.9 GP added outputType "CsvPerFieldRowFrColPtName"
1.10 GP added FieldsOnStructPointGrid.copyToOtherGrid()
1.11 GP added FieldsOnUnstructuredPoints
1.12 GP added outputType "OneCsvRowFrColField"
1.13 GP added FieldsOnStructPointGrid.toVtkAsCellData()
1.14 TR added in/output-functionality of meshgrid-numpy-arrays to
    FieldsOnStructPointGrid
1.15 TR added gridPointCloudData to FieldsOnStructPointGrid
1.16 GP changed (incompatibly) renamed multipleDataFun argument of
    FieldsOnStructPointGrid.gridPointCloudData()
1.17 fixed 'ij'-ordering for gridToMeshGrid and fieldToMeshGrid, unique
    fallback from bae.future_01 for gridPointCloudData
1.18 GP added FieldsOnStructPointGrid.description, content read from vtk file
1.19 GP added optionally use vtk module for reading
"""

import csv
from itertools import izip

import pyvtk
from bae.future_01 import OrderedDict
from bae.field_01 import Field
from bae.mesh_01 import MeshStructuredPoints, MeshStructuredPointsRot, \
    MeshUnstructuredPoints, Mesh

from bae.log_01 import msg, MsgTicker


class ReadError(Exception):
    r"""Can be thrown by FieldsOnStructPointGrid.fromVtk() if the vtk file
    to be read contains unexpected data or structure.
    """
    pass


class FieldsOnTopology(object):
    r"""Base class for container for one or many fields defined on the same
    topology (i.e. mesh).

    (Equivalent to one vtk file)

    @ivar notDefinedValue: A dictionary {field name: default value} stating
      default values for positions where there is no corresponding value in
      the odb. Or for points outside the FEM-mesh...
      For example:
       >>> ... notDefinedValue = { "SDV4": -1E20, "U": None }

      It's up to you whether the type of those values conforms to the other
      ordinary values from the odb. I.e. "U" would generally be something like
      [0.1, 0.0, 0.0], whereas None is not a list.

      This dictionary is empty by default.

      For values not explicitly listed in this dictionary the class defaults
      from L{FieldsOnTopology.defaultNotDefinedValue} apply.
    """

    defaultNotDefinedValue = {
        "scalar": 0.0,
        "vector": [0.0, 0.0, 0.0],
        "tensor": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        }

    def __init__(self, notDefinedValue=None):
        if notDefinedValue is None:
            self.notDefinedValue = dict()
        else:
            self.notDefinedValue = notDefinedValue

    def updateFields(self, *args):
        r"""Add field data for one or more fields. Replace values in each of
        the fields according to the instance variable
        L{FieldsOnStructPointGrid.notDefinedValue}.

        @param args: An arbitrary number of L{Field<bae.field_01.Field>}
          objects. The fieldName attribute of each item serves as key
          in self.data. Therefore the fieldname must be unique in order
          not to overwrite existing fields.

        @returns: self. So you can do the following:
          >>> vtk = FieldsOnStructPointGrid(mesh=mygrid).updateFields(mydata)

        @note: If a field with the same name already exists the old data will
          be overwritten without warning. The field order is not changed, i.e.
          an output file will contain the new data at the position of the
          initial insertion of data with the same fieldName attribute.
        """

        for newfield in args:

            fieldName = newfield.fieldName
            DataType = newfield.__class__

            # init replacing undefined values
            try:
                thisNotDefinedValue = self.notDefinedValue[fieldName]
            except KeyError:
                try:
                    thisNotDefinedValue = self.defaultNotDefinedValue[
                        newfield.dataType]
                except KeyError:
                    thisNotDefinedValue = None

            # copy values to new Field object replacing undefined values
            if thisNotDefinedValue is None:
                newfield = DataType(newfield)
            else:
                def replaceNotDefined(val):
                    if val is None:
                        return thisNotDefinedValue
                    else:
                        return val
                if issubclass(DataType, list):
                    newfield = DataType(replaceNotDefined(val)
                                        for val in newfield)
                elif issubclass(DataType, dict):
                    newfield = DataType(
                        (key, replaceNotDefined(val))
                        for key, val in newfield.iteritems())
                else:
                    raise NotImplementedError(
                        "Data type <%s> with position %s not implemented."
                        % (DataType, DataType.position))

            # a fieldName instance attribute may have overridden the class
            # attribute fieldName in the original field object passed in
            # we want the stored field to have the same name as the original
            # (possibly differing from the classes fieldName attribute)
            if newfield.fieldName != fieldName:
                newfield.fieldName = fieldName

            # actually store data
            self.data[fieldName] = newfield

        return self


class FieldsOnStructPointGrid(FieldsOnTopology):
    r"""
    A class to store field data (for one or more scalar, vector or tensor
    fields) on a structured points grid. With methods to read and write this
    data to/from vtk files.

    Stores only point data. Only for one time point. (equivalent to one vtk
    file)

    Example how to create a vtk file:
     >>> from bae.mesh_01 import MeshStructuredPoints
     >>> from bae.vtk_02 import FieldsOnStructPointGrid as Vtk
     >>> from bae.field_01 import Field
     >>> from math import sin
     >>> from bae.vecmath_01 import length, norm, vector_scale
     >>>
     >>> # create grid and field data (lists)
     >>> grid = MeshStructuredPoints(origin=[0,0,0], spacing=[0.05,0.05,0.05],
     >>>                             gridPtNb=[100,50,20])
     >>> FieldT = Field.classFromPosType(
     >>>     fieldName="T", position="structPt", dataType="scalar")
     >>> FieldU = Field.classFromPosType(
     >>>     fieldName="U", position="structPt", dataType="vector")
     >>> fldT = FieldT( ((z+1) * sin(x**2+y**2))
     >>>                for x,y,z in grid.getPointsIter() )
     >>> fldU = FieldU( grid.getPointsIter() )
     >>> for i, u in enumerate(fldU):
     >>>     if length(u)>2.0:
     >>>         fldU[i] = vector_scale(norm(u), 2.0)
     >>>
     >>> # export to vtk
     >>> vtk = Vtk(mesh=grid)
     >>> vtk.updateFields(fldT, fldU)
     >>> vtk.toVtk("mydata.vtk")

    And then read this very same vtk file:
     >>> from bae.vtk_02 import FieldsOnStructPointGrid as Vtk
     >>> vtk = Vtk().fromVtk("mydata.vtk")

    Now we have:
     - vtk.data["U"][100] is a list of the three components of the vector U at
       the 101st grid point.
     - vtk.data["T"].fieldName == "T"
     - len(list(vtk.mesh.getPointsIter())) == len(vtk.data["T"])

    @ivar mesh: A L{mesh_01.MeshStructuredPoints} object representing the
      output grid.
    @ivar data: An OrderedDict containing actual data fields.
    @ivar description: the vtk file header as read by L{fromVtk}()

    """

    def __init__(self, mesh=None, notDefinedValue=None):
        r"""

        @param mesh: A L{MeshStructuredPoints<bae.mesh_01.MeshStructuredPoints>}
          object representing the output grid
        @param notDefinedValue: Initializer for the corresponding instance
          variable. For a description see
          L{FieldsOnTopology.notDefinedValue}.

        @Note: Please supply the mesh argument as keyword arguments
          only since a positional argument is reserved for a possible
          implementation of a copy constructor.
        """
        FieldsOnTopology.__init__(self, notDefinedValue)
        self.mesh = mesh
        self.data = OrderedDict()

    def toVtk(self, outputFileName, description=None,
              outputFormat='binary', writeSinglePrecision=True):
        r"""Write stored data to vtk file of the given name.

        @param description: String for the second line of the vtk data file
          (title field). Defaults to self.description or (if that doesn't
          exist) "Field data".

        @param outputFormat: One of "binary", "ascii", "ascii_1row".
          "ascii_1row" means scalar fields are written as one row with space as
          delimiters instead of newline characters. (This is the pyvtk default
          ascii implementation but Glenn for example needs "ascii" with newline
          as delimiter.)
        @param writeSinglePrecision: If True then all binary data will be
          exported as single precision
        """

        if description is None:
            if hasattr(self, "description"):
                description = self.description
            else:
                description = "Field data"
        
        if isinstance(self.mesh, MeshStructuredPointsRot):
            origin = self.mesh.originBackRotated
        else:
            origin = self.mesh.origin

        vtkGrid = pyvtk.StructuredPoints(
            self.mesh.gridPtNb, origin, self.mesh.spacing )

        vtkPtData = pyvtk.PointData()
        for fieldName, field in self.data.iteritems():

            try:
                field = field.tolist()
            except AttributeError:
                pass

            if field.dataType == 'scalar':
                vtkPtData.append(
                    pyvtk.Scalars(field,name=fieldName,lookup_table="default"))
            elif field.dataType == 'vector':
                vtkPtData.append(
                    pyvtk.Vectors(field,name=fieldName))
            elif field.dataType == 'tensor':
                field = [ [ [x[0], x[3], x[4]],
                            [x[3], x[1], x[5]],
                            [x[4], x[5], x[2]] ]
                          for x in field ]
                vtkPtData.append(
                    pyvtk.Tensors(field,name=fieldName))

        v = pyvtk.VtkData(vtkGrid, description, vtkPtData)
        v.tofile(outputFileName, format=outputFormat,
                 writeSinglePrecision=writeSinglePrecision)

    def fromVtk(self, inputFileName):
        r"""Read field data on structured points grid from vtk file.
        Stores the vtk file header in self.description

        Usage:
         >>> from bae.vtk_02 import FieldsOnStructPointGrid as Vtk
         >>> vtk = Vtk().fromVtk("mydata.vtk")

        @returns: self. So you can do the following:
          >>> vtk = FieldsOnStructPointGrid().fromVtk("mydata.vtk")

        @raises ReadError: If the vtk file to be read contains unexpected data
          (e.g. cell data) or structure (i.e. dataset is not
          pyvtk.StructuredPoints).
        """

        try:
            from vtk import vtkStructuredPointsReader
            vtk = True
        except ImportError:
            vtk = None

        if vtk:
            # read the vtk file and do some checks
            reader = vtkStructuredPointsReader()
            reader.SetFileName(inputFileName)
            reader.ReadAllTensorsOn()
            reader.ReadAllVectorsOn()
            reader.ReadAllScalarsOn()
            reader.Update()
            output = reader.GetOutput()
            self.description = reader.GetHeader()
            if output.GetCellData().GetNumberOfArrays()>0:
                msg("WARNING: Found cell data in vtk file %s. Ignoring this."
                    % (inputFileName))
            if output.GetFieldData().GetNumberOfArrays()>0:
                msg("WARNING: Found field data in vtk file %s. Ignoring this."
                    % (inputFileName))

            # get the mesh/grid data
            gridParam = dict(
                gridPtNb=list(output.GetDimensions()),
                spacing=list(output.GetSpacing()),
                origin=list(output.GetOrigin()),
                )
            self.mesh = MeshStructuredPoints(**gridParam)

            # get the field i.e point data
            ptData = output.GetPointData()
            self.data = OrderedDict()
            for fieldCnt in range(ptData.GetNumberOfArrays()):
                fieldName = ptData.GetArrayName(fieldCnt)
                fieldData = ptData.GetArray(fieldCnt)  # or GetArray(fieldName)
                nComp = fieldData.GetNumberOfComponents()
                try:
                    # Nov 2019 on nukuD:
                    # fieldData is of type vtkCommonCorePython.vtkIntArray.
                    # This has a method GetNumberOfValues which delivers the
                    # right number. There is also a method GetSize but it
                    # yields a wrong number.
                    nVal = fieldData.GetNumberOfValues()
                except AttributeError:
                    # Nov 2019 on cassowary (Debian 10), vtk 6.3.0:
                    # fieldData is of type vtkobject. There is no member
                    # GetNumberOfValues, but Getsieze delivers the right number.
                    nVal = fieldData.GetSize()
                if nComp==1:
                    FieldType = Field.classFromPosType(
                        fieldName=fieldName, position="structPt",
                        dataType="scalar")
                    field = FieldType(
                        fieldData.GetValue(i) for i in xrange(nVal))
                    self.data[fieldName] = field  # store the field

                elif nComp==3:
                    FieldType = Field.classFromPosType(
                        fieldName=fieldName, position="structPt",
                        dataType="vector")
                    datIter = (fieldData.GetValue(i) for i in xrange(nVal))
                    field = FieldType( list(x) for x in izip(*[datIter]*nComp) )
                    self.data[fieldName] = field  # store the field

                elif nComp==9:
                    FieldType = Field.classFromPosType(
                        fieldName=fieldName, position="structPt",
                        dataType="tensor")
                    datIter = (fieldData.GetValue(i) for i in xrange(nVal))
                    fullTensor = zip(*[datIter]*nComp)
                    # check symmetry
                    isSymm = all(
                        # x[0][1]==x[1][0] ^ x[0][2]==x[2][0] ^ x[1][2]==x[2][1]
                        x[1]==x[3] and x[2]==x[6] and x[5]==x[7]
                        for x in fullTensor )
                    if not isSymm:
                        msg("WARNING: Tensor data for field %s in file %s"
                            " contains non symmetric values. Only the upper"
                            " triangle part of the matrix is being considered."
                            % (fieldName, inputFileName))
                    field = FieldType( (
                        # [x[0][0], x[1][1], x[2][2],
                        #  x[0][1], x[0][2], x[1][2]]
                        [x[0], x[4], x[8],
                         x[1], x[2], x[5]]
                        for x in fullTensor ) )
                    self.data[fieldName] = field  # store the field

            return self
        else:  # if not vtk
            # read the vtk file and do some checks
            v = pyvtk.VtkData( inputFileName )
            self.description = v.header
            if not(isinstance(v.structure, pyvtk.StructuredPoints)):
                raise ReadError(
                    "Improper dataset in file %s. Expected"
                    " pyvtk.StructuredPoints but found %s."
                    % (inputFileName, v.structure.__class__.__name__))
            if len(v.cell_data.data)>0:
                msg("WARNING: Found cell data in vtk file %s. Ignoring this."
                    % (inputFileName))

            # get the mesh/grid data
            self.mesh = MeshStructuredPoints(
                gridPtNb=list(v.structure.dimensions),
                origin=list(v.structure.origin),
                spacing=list(v.structure.spacing))

            # get the field i.e point data
            self.data = OrderedDict()
            for fieldData in v.point_data.data:
                fieldName = fieldData.name

                if isinstance(fieldData, pyvtk.Scalars):
                    FieldType = Field.classFromPosType(
                        fieldName=fieldName, position="structPt",
                        dataType="scalar")
                    field = FieldType( fieldData.scalars )
                    self.data[fieldName] = field  # store the field

                elif isinstance(fieldData, pyvtk.Vectors):
                    FieldType = Field.classFromPosType(
                        fieldName=fieldName, position="structPt",
                        dataType="vector")
                    field = FieldType( list(x) for x in fieldData.vectors )
                    self.data[fieldName] = field  # store the field

                elif isinstance(fieldData, pyvtk.Tensors):
                    FieldType = Field.classFromPosType(
                        fieldName=fieldName, position="structPt",
                        dataType="tensor")
                    # check symmetry
                    isSymm = all(
                        x[0][1]==x[1][0] and x[0][2]==x[2][0]
                        and x[1][2]==x[2][1]
                        for x in fieldData.tensors )
                    if not isSymm:
                        msg("WARNING: Tensor data for field %s in file %s"
                            " contains non symmetric values. Only the upper"
                            " triangle part of the matrix is being considered."
                            % (fieldName, inputFileName))
                    field = FieldType( (
                        [x[0][0], x[1][1], x[2][2],
                         x[0][1], x[0][2], x[1][2]]
                        for x in fieldData.tensors ) )
                    self.data[fieldName] = field  # store the field

                elif isinstance(fieldData, pyvtk.Field):
                    for fieldName in sorted(fieldData.data):
                        data = fieldData.data[fieldName]
                        if not(len(data)):
                            continue
                        try:
                            ndim = len(data[0])
                        except TypeError:
                            # if nb of components == 1 then (our modified) pyvtk
                            # delivers flat array
                            FieldType = Field.classFromPosType(
                                fieldName=fieldName, position="structPt",
                                dataType="scalar")
                            field = FieldType(data)
                        else:
                            if ndim==1:
                                FieldType = Field.classFromPosType(
                                    fieldName=fieldName, position="structPt",
                                    dataType="scalar")
                                field = FieldType( x[0] for x in data )
                            else:
                                msg("WARNING: The point data attribute %s"
                                    " contains the array %s with %d components."
                                    " It will be stored as vector field."
                                    % (fieldData.name, fieldName, ndim) )
                                FieldType = Field.classFromPosType(
                                    fieldName=fieldName, position="structPt",
                                    dataType="vector")
                                field = FieldType( list(x) for x in data )
                        self.data[fieldName] = field  # store the field

                else:
                    msg("WARNING: Point data attribute %s of type %s will be"
                        " ignored. This type is not recognized yet in"
                        " bae.vtk_02.FieldsOnStructPointGrid.fromVtk()."
                        % (fieldName, str(fieldData.__class__)))
                    continue

            return self

    def copyToOtherGrid(self, newGrid, defaultValueDict={}):
        """
        @param newGrid: a L{MeshStructuredPoints} object.

        @param defaultValueDict: dict { fieldName : default value } defines
           values to be taken for new grid points that have no corresponding
           point in the old grid. I.e. for grid points outside the old grid.

        @returns: new L{FieldsOnStructPointGrid} object with the same data on
           the new grid.

        @Note: Currently only works if the spacing is identical for old and new
           grids.
        @Note: Use L{MeshStructuredPoints.getPointAlignedGrid}() to make the
           new grid align to old grid points.
        """

        # check: spacing is the same
        if newGrid.spacing != self.mesh.spacing:
            raise ValueError(
                "ERROR in FieldsOnStructPointGrid.copyToOtherGrid():"
                " spacings differ! Old %s, new %s"
                % (self.mesh.spacing, newGrid.spacing))

        # get the closest grid point in new grid for the origin of self
        oldToNew = newGrid.getClosestGridPointIdx(self.mesh.origin)
        msg("Offset for the old grid in the new one: %s" % oldToNew,
            debugLevel=1)

        # each point in old (old) grid needs its index in the new grid
        msg("Calculating new grid indexes for grid points of old grid.")
        oldRanges = list()
        newRanges = list()
        for dim in range(3):
            # first assume old is completely in new
            old = [0, self.mesh.gridPtNb[dim]]
            new = [oldToNew[dim], oldToNew[dim]+self.mesh.gridPtNb[dim]]

            # now check that new boundaries don't exceed the grid
            if new[0]<0:
                old[0] -= new[0]
                new[0] = 0
            if new[1]>newGrid.gridPtNb[dim]:
                old[1] -= new[1]-newGrid.gridPtNb[dim]
                new[1] = newGrid.gridPtNb[dim]

            # store
            oldRanges.append(old)
            newRanges.append(new)

        oldStrides = self.mesh.strides  # abbreviation
        newStrides = newGrid.strides  # abbreviation
        oldToNewGridIds = [
            (oi*oldStrides[0] + oj*oldStrides[1] + ok*oldStrides[2],
             ni*newStrides[0] + nj*newStrides[1] + nk*newStrides[2])
            for ok, nk in izip(xrange(oldRanges[2][0], oldRanges[2][1]),
                               xrange(newRanges[2][0], newRanges[2][1]))
            for oj, nj in izip(xrange(oldRanges[1][0], oldRanges[1][1]),
                               xrange(newRanges[1][0], newRanges[1][1]))
            for oi, ni in izip(xrange(oldRanges[0][0], oldRanges[0][1]),
                               xrange(newRanges[0][0], newRanges[0][1]))]

        # convert data to new new grid
        msg("Converting data to new grid...")
        newVtk = FieldsOnStructPointGrid(mesh=newGrid)

        for oldField in self.data.itervalues():
            msg(" ... Processing field %s" % oldField.fieldName)

            # assign default value to the whole field
            try:
                defaultValue = defaultValueDict[oldField.fieldName]
            except KeyError:
                defaultValue = self.defaultNotDefinedValue[oldField.dataType]
            newField = type(oldField)([defaultValue]*len(newGrid))

            # transfer actual values
            for i, j in oldToNewGridIds:
                newField[j] = oldField[i]

            newVtk.updateFields(newField)

        # fini
        return newVtk

    def toVtkAsCellData(self, outputFileName, description="Field data",
                        outputFormat='binary', writeSinglePrecision=True):
        """Shift the mesh by half a cell size and export all fields as cell
        data instead of point data to a vtk file of the given name.

        @param outputFileName: of the new vtk file to be created.

        @param description: String for the second line of the vtk data file
          (title field).

        @param outputFormat: One of "binary", "ascii", "ascii_1row".
          "ascii_1row" means scalar fields are written as one row with space as
          delimiters instead of newline characters. (This is the pyvtk default
          ascii implementation but Glenn for example needs "ascii" with newline
          as delimiter.)
        @param writeSinglePrecision: If True then all binary data will be
          exported as single precision
        """

        if isinstance(self.mesh, MeshStructuredPointsRot):
            origin = self.mesh.originBackRotated
        else:
            origin = self.mesh.origin
        origin = [x-0.5*h for x,h in zip(origin, self.mesh.spacing)]
        gridPtNb = [i+1 for i in self.mesh.gridPtNb]

        vtkGrid = pyvtk.StructuredPoints(
            gridPtNb, origin, self.mesh.spacing )

        vtkCellData = pyvtk.CellData()
        for fieldName, field in self.data.iteritems():

            if field.dataType == 'scalar':
                vtkCellData.append(
                    pyvtk.Scalars(field,name=fieldName,lookup_table="default"))
            elif field.dataType == 'vector':
                vtkCellData.append(
                    pyvtk.Vectors(field,name=fieldName))
            elif field.dataType == 'tensor':
                field = [ [ [x[0], x[3], x[4]],
                            [x[3], x[1], x[5]],
                            [x[4], x[5], x[2]] ]
                          for x in field ]
                vtkCellData.append(
                    pyvtk.Tensors(field,name=fieldName))

        v = pyvtk.VtkData(vtkGrid, description, vtkCellData)
        v.tofile(outputFileName, format=outputFormat,
                 writeSinglePrecision=writeSinglePrecision)

    def _xyzSubGrid(self, xRange, yRange, zRange, returnType='linIdx'):
        """
        Parse the range arguments to returnType =:
            - 'linIdx' (d) --> (min,max)-indexTuple in linear dim-sequence
            - 'linIdxSlice' --> linIdx as a slice-object
            - 'coordinates' --> linear sub-sequence for x,y,z dimension

        Each range can be specified as:
            - a single value --> the nearest index will be returned
            - a range (min, max) --> all indizes within the range
                (including bounds) will be returned

        If a range is given as None then take all points along that axis. A min
        or max value of None represents infinity. I.e. None acts like a ":" as
        array index.

        Example:
         >>> # x = nearestTo(15.5), y <= 12., z = all
         >>> idx = myVtk._xyzSubGrid( 15.5, (None, 12.), None )

        @note: this function requires numpy
        """
        try:
            import numpy as np
        except ImportError:
            raise ImportError("This function needs numpy")

        from math import floor, ceil

        if isinstance(self.mesh, MeshStructuredPointsRot):
            origin = self.mesh.originBackRotated
        else:
            origin = self.mesh.origin

        returnVals = []

        # all directions to (lower, upper)-tuple
        for iq, q in enumerate(['x', 'y', 'z']):
            qRange = locals()[q + 'Range']

            # deltaQ
            dq = self.mesh.spacing[iq]

            # pysical length of coordinate range
            lq = (self.mesh.gridPtNb[iq] - 1) * dq

            if not hasattr(qRange, '__iter__'):
                qRange = (qRange, qRange)

            qRange = list(qRange)

            ## replace None by max/min-Range
            if qRange[0] is None:
                qRange[0] = origin[iq]

            if qRange[1] is None:
                qRange[1] = origin[iq] + lq

            if qRange[0] > qRange[1]:
                raise ValueError("Requested range for %s is unsorted" % q)


            ## get related (upper,lower)-indextuple
            if qRange[0] == qRange[0]:
                iqRange = (round((qRange[0] - origin[iq]) / float(dq)),
                           round((qRange[1] - origin[iq]) / float(dq)))
            else:
                iqRange = (ceil((qRange[0] - origin[iq]) / float(dq)),
                           floor((qRange[1] - origin[iq]) / float(dq)))

            iqRange = map(int, iqRange)

            ## check index-range and reset to valid
            if iqRange[0] < 0:
                msg("Requested Minimum for %s is below gridBound. " % q +
                    "Minimum will be set to gridBound.")

                if iqRange[0] == iqRange[1]:
                    iqRange = [0, 0]
                else:
                    iqRange[0] = 0

            if iqRange[1] >= self.mesh.gridPtNb[iq]:
                msg("Requested Maximum for %s is below gridBound. " % q +
                    "Maximum will be set to gridBound.")
                ii = self.mesh.gridPtNb[iq] - 1

                if iqRange[0] == iqRange[1]:
                    iqRange = [ii, ii]
                else:
                    iqRange[0] = [ii]

            ### specify return
            if returnType == 'linIdx':
                # will return [(ixmin,ixmax), (iymin,iymax), (izmin,izmax)]
                returnVals.append(iqRange)

            elif returnType == 'linIdxSlices':
                # will return [ixmin:ixmax, iymin:iymax, izmin:izmax]
                returnVals.append(slice(iqRange[0], iqRange[1] + 1, 1))

            elif returnType == 'coordinates':
                # will return [np.array(x[:]), np.array(y[:]), np.array(z[:])]
                coor = np.arange(iqRange[0], iqRange[1] + 1) * dq + origin[iq]

                returnVals.append(coor)

            else:
                raise ValueError("Don't understand returnType %s."
                                 % returnType)

        return returnVals

    def gridToMeshGrid(self, xRange=None, yRange=None, zRange=None,
                       preserveDim=True):
        """
        Returns the meshgrid-numpy arrays for coordinates. I.e. three arrays
        X[i,j,k], Y[i,j,k], Z[i,j,k] with the coordinates of each point of
        self.mesh, points identified through their index tuple (i,j,k).

        See also 
        L{self.mesh.getIdsInBox<bae.mesh_01.MeshStructuredPoints.getIdsInBox>}.

        If you whant to get the whole grid, this is equivalent to:
            >>> m = myVtk.mesh
            >>> o, s, n = m.origin, m.spacing, m.gridPtNb
            >>> x = o[0] + s[0] * np.arange(n[0])
            >>> y = o[1] + s[1] * np.arange(n[1])
            >>> z = o[2] + s[2] * np.arange(n[2])
            >>> X,Y,Z = numpy.meshgrid(x, y, z, indexing='ij')
            >>> #Mind the indexing='ij' kwarg!

        If you just need a subset of the grid, you can specify ranges.
        Each range can be defined as:
            - a single value --> the nearest gridIndex will be used
            - a range (min, max) --> all gridIndizes within the range
                (including bounds) will used

        If a range is given as None then take all points along that axis. A min
        or max value of None represents infinity. I.e. None acts like a ":" as
        array index.

        Example:
         >>> # get x = nearestTo(15.5), y <= 12., z = all
         >>> xyz = myVtk.gridToMeshGrid(xRange = 15.5,
         ...                            yRange = (None, 12.),
         ...                            zRange = None)

        @param xRange: x-range as float (take nearest) or list/tuple
            holding min and max
        @param yRange: y-range
        @param zRange: z-range
        @param preserveDim: if False, the fields will be squeezed
            (default True)
        @returns: meshgridded X, Y, Z
        @note: this function requires numpy
        """
        try:
            import numpy as np
        except ImportError:
            raise ImportError("This function needs numpy")

        #Create a Nx-Ny-Nz-dimed meshgrid. Mind the indexing='ij' kwarg!
        meshGrid = np.meshgrid(*self._xyzSubGrid(xRange, yRange, zRange,
                                                 returnType='coordinates'),
                              indexing='ij')

        if not preserveDim:
            for i in range(3):
                meshGrid[i] = meshGrid[i].squeeze()

        return meshGrid[0], meshGrid[1], meshGrid[2]

    def fieldToMeshGrid(self, fieldName,
                        xRange=None, yRange=None, zRange=None,
                        preserveDim=True):
        """
        Returns the meshgrid-numpy arrays for coordinates.
        If you whant to get the whole grid, this is equivalent to:
         >>> field = myVtk.data[fieldName]
         >>> gridPtNb = myVtk.mesh.gridPtNb
         >>> field = np.asarray(field).reshape(*gridPtNb, order = 'F')
         >>> field = np.swapaxes(field, 1, 0)


        If you just need a subset of the grid, you can specify ranges.
        Each range can be defined as:
            - a single value --> the nearest gridIndex will be used
            - a range (min, max) --> all gridIndizes within the range
                (including bounds) will be used

        If a range is given as None then take all points along that axis. A min
        or max value of None represents infinity. I.e. None acts like a ":" as
        array index.

        Example:
         >>> # get x = nearestTo(15.5), y <= 12., z = all
         >>> field = myVtk.fieldToMeshGrid('myField',
         ...                            xRange = 15.5,
         ...                            yRange = (None, 12.),
         ...                            zRange = None)

        @param fieldName: valid field-key
        @param xRange: x-range as float (take nearest) or list/tuple
            holding min and max
        @param yRange: y-range
        @param zRange: z-range
        @param preserveDim: if false, the fields will be squeezed
            (default True)
        @returns: meshgridded numpy array of fieldValues
        @note: this function requires numpy
        @warning: only scalar-fields implemented yet
        """
        try:
            import numpy as np
        except ImportError:
            raise ImportError("This function needs numpy")

        iSlices = self._xyzSubGrid(xRange, yRange, zRange,
                                   returnType='linIdxSlices')

        gridPtNb = self.mesh.gridPtNb

        # reshape to an X x Y x Z - array
        # numpys-axes-order is C-type. To reshape the flat vtk-fieldData-list
        # to a meshgrid-ordered array have to use order = 'F'
        field = self.data[fieldName]

        if not field.dataType == 'scalar':
            raise NotImplementedError("Implemented for scalar fields only.")

        field = np.asarray(field).reshape(*gridPtNb, order='F')

        field = field[iSlices[0], iSlices[1], iSlices[2]]

        if not preserveDim:
            field = field.squeeze()

        return field

    def meshGridToField(self, newFieldName, meshGrid):
        """
        Stores a meshgrided numpy array of field values in self.

        Example:
         >>> field = myVtk.fieldToMeshGrid('myField')
         >>> grad = np.gradient(field, *myVtk.mesh.spacing)
         >>> myVtk.meshGridToField('myGradient', grad)

        @param newFieldName: name of the field
        @param meshGrid: meshgrided numpy array
        @note: this function requires numpy
        """

        try:
            import numpy as np
        except ImportError:
            raise ImportError("This function needs numpy")

        from bae.field_01 import createFieldObject

        if not tuple(self.mesh.gridPtNb) == tuple(meshGrid.shape):
            raise ValueError("Dimension-mismatch: " +
                             "meshGrid.shape=%s " % str(meshGrid.shape) +
                             "vs.  gridPtNb=%s" % str(self.mesh.gridPtNb) +
                             "Did you forgot to use indexing='ij' kwarg " +
                             "while calling numpy.meshgrid?")

        # numpys-axes-order is C-type. To ravel an meshgrid-ordered array to a
        # vtk-fieldData-list we have to use order = 'F'
        meshGrid = meshGrid.ravel(order='F').tolist()

        out = createFieldObject(newFieldName, "structPt", "scalar",
                                initArgs=[meshGrid,])
        self.updateFields(out)

    def gridPointCloudData(self, newFieldName,
                           xyz, inData, noDataValue=0.,
                           multipleDataFun='mean'):
        """
        Grids unstructured PointCloudData to current grid and stores
        the new field.
        No interpolation is performed here. If a voxel contains
        one or more points of the point cloud it will get a value
        (see multipleDataFun for details) otherwise the voxel value
        will be set to the value of noDataValue (default 0).

        @param newFieldName: name of the new gridded field

        @param xyz: list/numpyArray of x,y,z point data

        @param inData: data at xyz

        @param noDataValue: points not captured by xyz will get this
            value (default = 0.)

        @param multipleDataFun: defines how to handle multiple
            point data in single voxel ('mean' (default), 'min',
            'max')
        """

        try:
            import numpy as np
        except ImportError:
            raise ImportError("This function needs numpy")

        from bae.future_01 import unique

        #prepare input data
        xyz, inData = np.atleast_2d(xyz), np.array(inData)

        #prepare current mesh data
        gridPtNb = self.mesh.gridPtNb
        spacing = np.array(self.mesh.spacing)
        if isinstance(self.mesh, MeshStructuredPointsRot):
            origin = np.array(self.mesh.originBackRotated)
        else:
            origin = np.array(self.mesh.origin)

        # voxelize coordinates --> [(iX,iY,iZ), ...]
        xyzV = np.floor((xyz-origin)/spacing).astype(int)

        # use only values within current grid
        validIdx = (np.all(xyzV >= 0, axis=1) &
                    np.all(xyzV < gridPtNb, axis=1))
        xyzV = xyzV[validIdx]
        inData = inData[validIdx]

        # get unique voxel-indices
        uniqueIdx, idx, inv = unique(
            xyzV, axis=0, return_index=True, return_inverse=True)

        # create 'empty/noDataValue' array
        data = np.ones(gridPtNb) * noDataValue

        # loop valid voxels
        # i do not know how to do this in a vectorized way
        for ii, uIdx in enumerate(uniqueIdx):
            iX,iY,iZ = uIdx

            datIdx = (inv == ii)

            if multipleDataFun == 'mean':
                dat = inData[datIdx].mean()

            elif multipleDataFun == 'max':
                dat = inData[datIdx].max()

            elif multipleDataFun == 'min':
                dat = inData[datIdx].min()

            else:
                raise ValueError(('Do not know a function %s ' +
                                  'to reduce multiple data in ' +
                                  'single voxel.') % multipleDataFun)

            data[iX,iY,iZ] = dat

        self.meshGridToField(newFieldName, data)


class FieldsOnUnstructuredPoints(FieldsOnTopology):
    """A class to store field data (for one or more scalar, vector or tensor
    fields) on unstructured points. With methods to read (and eventually write)
    this data to/from csv files.

    Currently this class only reads from a csv file, see L{fromCsvFile}.

    Usage:
     >>> from bae.vtk_02 import FieldsOnUnstructuredPoints
     >>> 
     >>> fields = FieldsOnUnstructuredPoints.fromCsvFile("mydata.csv")
     >>> fields.mesh       # a MeshUnstructuredPoints object
     >>> fields.data["T"]  # a Field-object

    @ivar mesh: L{MeshUnstructuredPoints<bae.mesh_01.MeshUnstructuredPoints>}
        object
    @ivar data: dictionary {field name: field} with field being a scalar field
        object, see L{bae.field_01}.
    """

    @classmethod
    def fromCsvFile(cls, inp, joinTensors=False):
        """Initialize object with data from csv file

        Reads from a csv file with a single header row that specifies the
        field names. First three columns must be the point coordinates
        following columns are field values according to the header row.

        I.e.::
          x,    y,    z,    T,    U
          10.1, -3.1, 5.1, 293,   0.0
          10.1, -3.1, 5.2, 293.5, 1.3
          10.1, -3.1, 5.3, 293.7, 5.7
          10.1, -3.1, 5.4, 294.5, 2.5
          10.1, -3.1, 5.5, 295.2, 0.2
          10.1, -3.1, 5.6, 297.7, -2.4
        """

        if joinTensors:
            raise NotImplementedError(
                "FieldsOnUnstructuredPoints.fromCsvFile joinTensors=True not"
                " implemented yet.")

        self = cls()

        if isinstance(inp, basestring):
            inp = open(inp, "rb")
        tab = csv.reader(inp)
        self.headLine = tab.next()
        self.fieldNames = self.headLine[3:]
        if len(set(self.fieldNames)) != len(self.fieldNames):
            msg("ERROR: Found ambigous field names.")
            return
        N = len(self.fieldNames)
        points = list()
        data = list()
        ignoredLines = 0
        for row in tab:
            # ignore lines with wrong number of items
            if len(row) != 3+N:
                ignoredLines += 1
                continue
            points.append(map(float, row[0:3]))
            data.extend(map(float, row[3:3+N]))

        # diagnostic output
        msg("Read %d valid data lines." % len(points))
        if ignoredLines:
            msg("WARNING: Ignored %d lines with wrong number of columns."
                % ignoredLines)

        # create fields
        self.data = OrderedDict()
        for i, fieldName in enumerate(self.fieldNames):
            self.data[fieldName] = Field.classFromPosType(
                fieldName, "point", "scalar")(data[i::N])

        self.mesh = MeshUnstructuredPoints(points)
        return self


# Note: this should be split into a base class with two child classes, one for
# the "OneCsvRowPtColFldFr" case and one for the "CsvPerFrameRowPtColFld" case.
# Note: reading from csv tables not implemented yet.

class FramesFieldsOnUnstructuredPoints(FieldsOnTopology):
    r"""A class to store field data (for one or more scalar, vector or tensor
    fields, see module L{bae.field_01}) on unstructured points for multiple
    time frames. With methods to read and write this data to/from csv files.

    After creating this object, you need to call self.initOutput. This
    determines the type of output.

    Usage:
     >>> from bae.vtk_02 import FramesFieldsOnUnstructuredPoints
     >>> from bae.field_01 import createFieldObject
     >>> 
     >>> points = [ [123.3, 32.5, -10.0], [145.3, 32.5, -10.0], ...]
     >>> framesList = [ (2, i) for i in [1, 18] + range(20, 73, 2) ]
     >>> fieldNames = [ "U_1", "U_2", "U_3", "SDV4", "S_MIN_PRINCIPAL" ]
     >>>
     >>> # keyword args please
     >>> output = FramesFieldsOnUnstructuredPoints(points=points)
     >>> output.initOutput("MyResultsOnPoints.csv", "OneCsvRowPtColFldFr")
     >>>
     >>> # constant values first
     >>> zFlagField = createFieldObject(
     ...     fieldName="zFlag", position="point", dataType="scalar",
     ...     initArgs=[ ((z>200) for z in points), ])
     >>> output.updateFields(zFlagField)
     >>> output.storeConstantData()
     >>>
     >>> # fields that change with frame
     >>> for stepNb, frameNb, odbFrame in odb.getFrameIterator3(framesList):
     >>>     for fieldName in fieldNames:
     >>>         field = odb.getFieldFromOdb(
     >>>             fieldName=fieldName, odbFrame=odbFrame)
     >>>         output.updateFields(field)
     >>>     output.storeFrame((stepNb, frameNb))
     >>> output.closeCsvOutput()

    @ivar mesh: a L{MeshUnstructuredPoints} object. (Typically) initialized
         with all points passed in through self.__init__()
    @ivar data: dict containing the fields for the current frame. Being emptied
         by self.storeFrame()
    @ivar frameData: {frame id : {field id : field data } }
         field data might be something like a list of tuples or so.
         Exists only if self.outputType == "OneCsvRowPtColFldFr".
    @Note: Design principles:
       - different output formats (csv files), common API, select output format
         by single parameter
       - big data: only keep in memory what's necessary, write as early as
         possible
    """

    # default: leave not defined values as None, yielding empty fields in the
    # csv output
    defaultNotDefinedValue = dict()

    def __init__(self, points=None, notDefinedValue=None):
        r"""
        @param points: A list of point coordinates or a
          L{MeshUnstructuredPoints<bae.mesh_01.MeshUnstructuredPoints>}
          object representing the output points

          This argument is optional, you must supply the self.mesh object as
          L{MeshUnstructuredPoints}-object at any time before
          the first call to storeFrame().

        @param notDefinedValue: Initializer for the corresponding instance
          variable. For a description see
          L{FieldsOnTopology.notDefinedValue}.

        @Note: Please supply the points argument as keyword argument only
          since a positional argument is reserved for a possible
          implementation of a copy constructor.
        """
        FieldsOnTopology.__init__(self, notDefinedValue)

        if isinstance(points, MeshUnstructuredPoints):
            self.mesh = points
        else:
            try:
                self.mesh = MeshUnstructuredPoints(points)
            except TypeError:
                self.mesh = points

        self.data = OrderedDict()
        self.fieldClasses = None

    def _storeFrameIdToName(self, frameIdToName):
        """used by the initCsvOutput-functions
        store a self.frameIdToName method from arbitrary input

        @param frameIdToName: Defines the relation between the frame id (used
           as key in the self.frameData dictionary) and the frame name used in
           the csv output file.

           Might be a function frameId -> frameName or a template string (like
           the default value).

           Note that the type of the frame ids is not determined by this class.
           Frame ids will be (stepNb, frameNb)-tuples preferably and only then
           the default value for frameIdToName makes much sense.
           But you could as well use other hashable types like frame names
           (i.e. strings) or plain frame numbers. Then you have to specify a
           adequate frameIdToName parameter here.
        """
        if isinstance(frameIdToName, basestring):
            self.frameIdToNameTemplate = frameIdToName
            self.frameIdToName = (
                lambda frameId: (self.frameIdToNameTemplate % frameId))
        else:
            self.frameIdToName = frameIdToName

    def initOutput(
            self, fileNameTemplate, outputType="CsvPerFrameRowPtColFld",
            filterEmpty=False, frameIdToName="F%d%03d"):
        r"""
        Initialize output.

        @param fileNameTemplate: Output file name or template string for
           all output files.

           If outputType starts with "OneCsvRow..." then this argument shall
           contain the complete output file name.

           If outputType="CsvPerFrameRowPtColFld": Template string for
           all csv files. Must contain one %s placeholder for the frame
           name string.

           If outputType="CsvPerFieldRowPtColFr": Template string for all
           csv files. Must contain one %s placeholder for the field name
           string.

        @param outputType:
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
           - "CsvPerFieldRowFrColPtName": one file per field.
             Columns are: frame name then one column for each point having the
             name in the column header. One row per frame.
             There has to be a constant field "ptName" (case matters) for this
             outputType to work. Supply it by means of the L{updateFields} and
             L{storeConstantData} methods.
              >>> from bae.field_01 import createFieldObject
              >>> output = FramesFieldsOnUnstructuredPoints(points=points)
              >>> output.initOutput("result.csv", "CsvPerFieldRowFrColPtName")
              >>> ptNames = createFieldObject("ptName", "point", "str",
              ...               initArgs=[... list of names ...])
              >>> output.updateFields(ptNames)
              >>> output.storeConstantData()
              >>> # ... now store the usual frame field data

             In case you have a list of [x,y,z,name]-tuples you can separate
             them like this:
              >>> for itertools import izip
              >>> points, ptNames = izip(*((x[:3], x[3])
              ...                          for x in pointsWithNames))

             All other constant fields are being treated like a separate frame
             written at the end of each file.
           - "OneCsvRowFrColPtField": one file for the whole data.
             Columns are: frame, x,y,z, field1, field2,... x,y,z, field1,
             field2,... I.e. one row per frame and different points in columns.
             For each point we have the points coordinates first and then the
             field data.
           - "OneCsvRowPtFldColFr": one file for the whole data.
             Columns are: x,y,z, var, frame1, frame2, frame3, ...
             The var-column contains the field name. Order of rows is all points
             field 1 then all points field 2 and so on. Vector and tensor
             components come on succesive rows for the same point. I.e.::
               x,  y,  z,  var,  F2001, F2002, ...
               x1, y1, z1, logP, 0.0,   0.0,  ...
               x2, y2, z2, logP, 0.0,   0.1,  ...
               x3, y3, z3, logP, 0.0,   0.2,  ...
               x1, y1, z1, U_1,  0.0,  -0.1,  ...
               x1, y1, z1, U_2,  0.0,   0.2,  ...
               x1, y1, z1, U_3,  0.0,  -0.3,  ...
               x2, y2, z2, U_1,  0.0,  -0.1,  ...
               x2, y2, z2, U_2,  0.0,   0.2,  ...
               x2, y2, z2, U_3,  0.0,  -0.3,  ...
               ...

           - "OneCsvRowFrColField": one file for all; constant fields (usually:
             ptName), then frames in rows, no point coordinates.
             Columns are: frame name then first field with one column for each
             point. Then second field for all points.
             In the header line we have the field names.
             One row per frame.
             Usually there will be a constant field "ptName" to identify the
             points supplied by means of the L{updateFields} and
             L{storeConstantData} methods. See description of outputType
             "CsvPerFieldRowFrColPtName".

        @param filterEmpty: supress rows (points) without any data. This only
           works if None-values are *not* converted (to zeros for example)
           Only implemented for outputTypes CsvPerFrameRowPtColFld,
           OneCsvRowPtColFldFr and CsvPerFieldRowPtColFr so far.

        @param frameIdToName: Defines the relation between the frame id (used
           as key in the self.frameData dictionary) and the frame name used in
           the csv output file.

           Might be a function frameId -> frameName or a template string (like
           the default value).

           Note that the type of the frame ids is not determined by this class.
           Frame ids will be (stepNb, frameNb)-tuples preferably and only then
           the default value for frameIdToName makes much sense.
           But you could as well use other hashable types like frame names
           (i.e. strings) or plain frame numbers. Then you have to specify a
           adequate frameIdToName parameter here.
        """
        self.outputType = outputType
        self._storeFrameIdToName(frameIdToName)
        self.filterEmpty = filterEmpty

        # date structure to store intermediate data
        # constant fields container needed for all output types
        self.constData = OrderedDict()

        if outputType in ("OneCsvRowPtColFldFr", "OneCsvRowPtFldColFr"):
            # data structure to store intermediate data
            self.frameData = OrderedDict()
            # results file
            self.fileName = fileNameTemplate
            self.output = csv.writer(open(self.fileName, "wb"))

        elif outputType == "CsvPerFrameRowPtColFld":
            # no data structure needed to store intermediate data
            self.headline = None
            # results file name template
            self.fileNameTemplate = fileNameTemplate

        elif outputType.startswith("CsvPerFieldRow"):
            # data structure to store intermediate data
            self.frameData = OrderedDict()
            self.headline = None
            # results file name template
            self.fileNameTemplate = fileNameTemplate

        elif outputType in ("OneCsvRowFrColPtField", "OneCsvRowFrColField"):
            # no data structure needed to store intermediate data
            self.fileName = fileNameTemplate
            ot = outputType
            if outputType=="OneCsvRowFrColPtField (TimeLine)":
                ot = ot + ""
            msg("Results will go to %s in the %s-format" % (self.fileName, ot))
            # results file
            self.output = csv.writer(open(self.fileName, "wb"))
            # make sure the header is only written once
            self.writeHeaderFlag = True
            self.filterEmpty = None  # Ignored, does not work!
        else:
            raise ValueError(
                "Unknown outputType %s specified. (Maybe this bae.vtk_02-module"
                " of version %s is outdated?)" % (outputType, __version__))

    def initCsvOutputAllInOne(
            self, fileName, filterEmpty=False, frameIdToName="F%d%03d"):
        r"""Deprecated alias for initOutput(outputType="OneCsvRowPtColFldFr")
        """
        self.initOutput(fileName, outputType="OneCsvRowPtColFldFr",
                        filterEmpty=filterEmpty, frameIdToName=frameIdToName)

    def initCsvOutputPerFrame(
            self, fileNameTemplate, filterEmpty=False, frameIdToName="F%d%03d"):
        r"""Deprecated alias for initOutput(outputType="CsvPerFrameRowPtColFld")
        """
        self.initOutput(fileNameTemplate, outputType="CsvPerFrameRowPtColFld",
                        filterEmpty=filterEmpty, frameIdToName=frameIdToName)

    @staticmethod
    def _appendVal(dataCols, FieldClass, val):
        """Append val to the list dataCols. If FieldClass is a multi-component
        type then use list.extend(). If val is None then add a number of Nones
        (for empty fields)
        """
        if (FieldClass.dataType=="scalar"
                or FieldClass.dataType=="str"):
            dataCols.append(val)
        elif val is None:
            # vector or tensor value missing for this point
            dataCols.extend((None,)*FieldClass.nbComponents)
        else:
            dataCols.extend(val)

    def storeFrame(self, frameId):
        """Finalize the data for one frame after it has been collected by one
        or more calls to self.updateFields(). In the "CsvPerFrameRowPtColFld"
        case the data will be written to an external csv file. In the
        "OneCsvRowPtColFldFr" it will be stored for later processing.

        @param frameId: frame identifier, usually a (step nb, frame nb)-tuple.
        But can be anything corresponding (compatible) to the frameIdToName
        argument to L{initOutput}.

        @Note: In the current version it is assumed that all frames contain the
        same number and order of fields. This has to be assured by ordering the
        calls to self.updateFields() accordingly.
        """

        # need to store the classes (not only the names) because classes also
        # know their type (vector, scalar) and component names
        if self.fieldClasses is None:
            self.fieldClasses = [
                field.__class__ for field in self.data.itervalues()]
        elif ( len(self.fieldClasses) != len(self.data)
                or any(
                    f1.fieldName != f2.fieldName
                    for f1,f2 in izip(self.fieldClasses,self.data.itervalues())
                    ) ):
            msg("ERROR: Fields in frame %s %s differ. Before we had %s.\n"
                "Only fields of this latest frame might be in the output."
                % (frameId, self.data.keys(),
                   [f1.fieldName for f1 in self.fieldClasses]))
            self.fieldClasses = [
                field.__class__ for field in self.data.itervalues()]

        if (self.outputType in ("OneCsvRowPtColFldFr", "OneCsvRowPtFldColFr")
                or self.outputType.startswith("CsvPerFieldRow")):
            # for those output types the file(s) will be written at the end
            # now we just store self.data in self.frameData
            self.frameData[frameId] = self.data
            self.data = OrderedDict()

        elif self.outputType == "CsvPerFrameRowPtColFld":
            # for this output type a new file is created for this frame
            # then the frame's data is written to this file immediately

            # prepare the head line
            if self.headline is None:
                self.headline = ["x", "y", "z"]
                for field in self.constData.itervalues():
                    FieldClass = field.__class__
                    self.headline.extend(FieldClass.getComponentNamesList())
                for FieldClass in self.fieldClasses:
                    self.headline.extend(FieldClass.getComponentNamesList())

            # open output file
            frameName = self.frameIdToName(frameId)
            fileName = self.fileNameTemplate % frameName
            output = csv.writer(open(fileName, "wb"))

            # write lines to output file
            msg("Writing results for frame '%s' to %s" % (frameName, fileName))
            output.writerow(self.headline)
            msg("Wrote headline: %s" % self.headline, debugLevel=10)

            ticker = MsgTicker("...wrote data for %d points.")
            numFilteredPts = 0
            for ptIndex, pt in enumerate(self.mesh):
                dataCols = []

                for field in self.constData.itervalues():
                    FieldClass = field.__class__
                    self._appendVal(dataCols, FieldClass, field[ptIndex])
                for FieldClass in self.fieldClasses:
                    fieldName = FieldClass.fieldName
                    self._appendVal(dataCols, FieldClass,
                                    self.data[fieldName][ptIndex])
                if self.filterEmpty and all(x is None for x in dataCols):
                    numFilteredPts += 1
                else:
                    output.writerow(pt+dataCols)
                ticker.msg(ptIndex+1)
            del output
            del ticker
            if numFilteredPts>0:
                msg("Filtered out %d points without data." % numFilteredPts)

        elif self.outputType=="OneCsvRowFrColPtField":
            # for this output type the data for this frame is appended
            # to an already open output file. I.e. it's being written now.

            if self.writeHeaderFlag:
                # make sure the header is only written once
                self.writeHeaderFlag = False

                # prepare the head line
                headline=["frame number"]
                for ptIndex,pt in enumerate(self.mesh):
                    if ptIndex>0:
                        # separator: empty column
                        headline.append("")
                    headline.extend(["x","y","z"])
                    for field in self.constData.itervalues():
                        FieldClass = field.__class__
                        headline.extend(FieldClass.getComponentNamesList())
                    for FieldClass in self.fieldClasses:
                        headline.extend(FieldClass.getComponentNamesList())
                self.output.writerow(headline)

            # write actual data line for this frame
            frameName = self.frameIdToName(frameId)
            msg("Writing results for frame '%s' to %s"
                % (frameName, self.fileName))

            numFilteredPts = 0
            dataCols = [frameName]
            for ptIndex, pt in enumerate(self.mesh):
                if ptIndex>0:
                    # separator: empty column
                    dataCols.append("")
                dataCols.extend(pt)
                for field in self.constData.itervalues():
                    FieldClass = field.__class__
                    self._appendVal(dataCols, FieldClass, field[ptIndex])
                for FieldClass in self.fieldClasses:
                    fieldName = FieldClass.fieldName
                    self._appendVal(dataCols, FieldClass,
                                    self.data[fieldName][ptIndex])

            self.output.writerow(dataCols)
            msg("Wrote data for %d points." % len(self.mesh))

        elif self.outputType=="OneCsvRowFrColField":
            # for this output type the data for this frame is appended
            # to an already open output file. I.e. it's being written now.

            if self.writeHeaderFlag:
                # make sure the header is only written once
                self.writeHeaderFlag = False

                # prepare the head line
                headline=["frame number"]
                for FieldClass in self.fieldClasses:
                    headline.extend(FieldClass.getComponentNamesList()
                                    *len(self.mesh))
                self.output.writerow(headline)

                # write constant data
                for constfield in self.constData.itervalues():
                    ConstFieldClass = constfield.__class__
                    try:
                        # nb of components of Const Field
                        nbCompCF = ConstFieldClass.nbComponents
                    except AttributeError:
                        nbCompCF = 1
                    dataCols = [constfield.fieldName]

                    for FieldClass in self.fieldClasses:
                        try:
                            # nb of components of Per Frame Field
                            nbCompPFF = FieldClass.nbComponents
                        except AttributeError:
                            nbCompPFF = 1

                        numEmptyCols = (nbCompPFF - nbCompCF)
                        for val in constfield:
                            if val is None:
                                val = [None]*FieldClass.nbComponents
                            elif numEmptyCols<0:
                                # Const value has more components than the
                                # per-frame-fields to follow in subsequent
                                # rows. Truncate value.
                                val = val[:numEmptyCols]
                            elif (ConstFieldClass.dataType=="scalar"
                                  or ConstFieldClass.dataType=="str"):
                                val = [val] + [None]*numEmptyCols
                            else:
                                val = val + [None]*numEmptyCols
                            dataCols.extend(val)

                    # write this row
                    self.output.writerow(dataCols)

            # write actual data line for this frame
            frameName = self.frameIdToName(frameId)
            msg("Writing results for frame '%s' to %s"
                % (frameName, self.fileName))

            numFilteredPts = 0
            dataCols = [frameName]

            for FieldClass in self.fieldClasses:
                fieldName = FieldClass.fieldName
                for val in self.data[fieldName]:
                    self._appendVal(dataCols, FieldClass, val)

            self.output.writerow(dataCols)
            msg("Wrote data for %d points." % len(self.mesh))

    def storeConstantData(self):
        r"""Store all previously added fields as constant fields, i.e.
        data that does not change with time. This will be written
        once only in the "OneCsvRowPtColFldFr" case. And it will be written
        identically to each individual frame's csv file in the other cases.

        Usage:
         >>> points = [ [123.3, 32.5, -10.0], [145.3, 32.5, -10.0], ...]
         >>> framesList = [ (2, i) for i in [1, 18] + range(20, 73, 2) ]
         >>> fieldNames = [ "U_1", "U_2", "U_3", "SDV4", "S_MIN_PRINCIPAL" ]
         >>>
         >>> # keyword args please
         >>> output = FramesFieldsOnUnstructuredPoints(points=points)
         >>> output.initCsvOutputAllInOne("MyResultsOnPoints.csv")
         >>> (matNamesList, matNumberField) = odb.getMatFieldFromOdb()
         >>> output.updateFields(matNumberField)
         >>> output.storeConstantData()
         >>>
         >>> for stepNb, frameNb, odbFrame in odb.getFrameIterator3(framesList):
         >>>     for fieldName in fieldNames:
         >>>         field = odb.getFieldFromOdb(
         >>>             fieldName=fieldName, odbFrame=odbFrame)
         >>>         output.updateFields(field)
         >>>     output.storeFrame((stepNb, frameNb))
         >>> output.closeCsvOutput()

        @Note: In the "CsvPerFrameRowPtColFld" and the "OneCsvRowFrColPtField"
        case this method must not be called after storeFrame() has been called
        for the first time.
        """
        self.constData.update(self.data)
        self.data = OrderedDict()

    def closeCsvOutput(self):

        if self.outputType == "OneCsvRowPtColFldFr":

            msg("Starting to write data for all frames to the csv file %s."
                % self.fileName)

            # prepare the head line
            headline = ["x", "y", "z"]
            for field in self.constData.itervalues():
                FieldClass = field.__class__
                headline.extend(FieldClass.getComponentNamesList())
            for FieldClass in self.fieldClasses:
                componentNames = FieldClass.getComponentNamesList()
                for frameId in self.frameData:
                    frameName = self.frameIdToName(frameId)
                    headline.extend(("%s-%s" % (cn, frameName))
                                    for cn in componentNames)
            self.output.writerow(headline)

            ticker = MsgTicker("...wrote data for %d points.")
            numFilteredPts=0
            for ptIndex, pt in enumerate(self.mesh):
                dataCols = []

                for field in self.constData.itervalues():
                    FieldClass = field.__class__
                    self._appendVal(dataCols, FieldClass, field[ptIndex])
                for FieldClass in self.fieldClasses:
                    fieldName = FieldClass.fieldName
                    for frameId, data in self.frameData.iteritems():
                        self._appendVal(dataCols, FieldClass,
                                        data[fieldName][ptIndex])
                if self.filterEmpty and all(x is None for x in dataCols):
                    numFilteredPts += 1
                else:
                    self.output.writerow(pt+dataCols)
                ticker.msg(ptIndex+1)
            del self.output
            del ticker
            if numFilteredPts>0:
                msg("Filtered out %d points without data." % numFilteredPts)

        if self.outputType == "OneCsvRowPtFldColFr":

            msg("Starting to write data for all frames to the csv file %s."
                % self.fileName)

            # prepare the head line
            headline = ["x", "y", "z", "var"]
            for frameId in self.frameData:
                frameName = self.frameIdToName(frameId)
                headline.append(frameName)
            self.output.writerow(headline)

            # initialize main output loop
            data = self.frameData.itervalues().next()
            ticker = MsgTicker("...wrote data for field %%d/%d."
                               % (len(data)+len(self.constData)))
            framesList = self.frameData.keys()
            numFrames = len(framesList)
            numFilteredPts = 0

            # write per frame data
            for FieldClass in self.fieldClasses:
                fieldName = FieldClass.fieldName
                componentNames = FieldClass.getComponentNamesList()
                if len(componentNames)>1:
                    for ptIdx, pt in enumerate(self.mesh):
                        for iCmp, cn in enumerate(componentNames):
                            self.output.writerow(pt+[cn]+[
                                self.frameData[frameId][fieldName][ptIdx][iCmp]
                                for frameId in framesList] )
                else:
                    for ptIdx, pt in enumerate(self.mesh):
                        self.output.writerow(pt+[fieldName]+[
                            self.frameData[frameId][fieldName][ptIdx]
                            for frameId in framesList] )
                ticker.tick()

            # write constant fields
            for field in self.constData.itervalues():
                FieldClass = field.__class__
                componentNames = FieldClass.getComponentNamesList()
                if len(componentNames)>1:
                    for iCmp, cn in enumerate(componentNames):
                        self.output.writerows(
                            pt+[cn]+[field[ptIndex][iCmp]]*numFrames
                            for ptIndex, pt in enumerate(self.mesh) )
                else:
                    fieldName = field.fieldName
                    self.output.writerows(
                        pt+[fieldName]+[field[ptIndex]]*numFrames
                        for ptIndex, pt in enumerate(self.mesh) )
                ticker.tick()

            del self.output
            del ticker

        if self.outputType == "CsvPerFieldRowPtColFr":

            # write constant fields
            for fieldName, field in self.constData.iteritems():
                FieldClass = field.__class__
                componentNames = FieldClass.getComponentNamesList()

                headline = ["x", "y", "z"]
                if len(componentNames)>1:
                    headline.extend(componentNames)

                # open output file
                fileName = self.fileNameTemplate % fieldName
                output = csv.writer(open(fileName, "wb"))

                # write lines to output file
                msg("Writing results for field %s to %s"
                    % (fieldName, fileName))
                output.writerow(headline)

                ticker = MsgTicker("...wrote data for %d points.")
                numFilteredPts = 0
                for ptIndex, pt in enumerate(self.mesh):
                    dataCols = []
                    self._appendVal(dataCols, FieldClass, field[ptIndex])
                    if self.filterEmpty and all(x is None for x in dataCols):
                        numFilteredPts += 1
                    else:
                        output.writerow(pt+dataCols)
                    ticker.msg(ptIndex+1)
                del output
                del ticker
                if numFilteredPts>0:
                    msg("Filtered out %d points without data." % numFilteredPts)

            # write framedata files, one for each field
            for FieldClass in self.fieldClasses:
                fieldName = FieldClass.fieldName
                componentNames = FieldClass.getComponentNamesList()

                headline = ["x", "y", "z"]
                if len(componentNames)>1:
                    for frameId in self.frameData:
                        frameName = self.frameIdToName(frameId)
                        headline.extend(("%s-%s" % (cn, frameName))
                                        for cn in componentNames)
                else:
                    headline.extend( self.frameIdToName(frameId)
                                     for frameId in self.frameData )

                # open output file
                fileName = self.fileNameTemplate % fieldName
                output = csv.writer(open(fileName, "wb"))

                # write lines to output file
                msg("Writing results for field %s to %s"
                    % (fieldName, fileName))
                output.writerow(headline)

                ticker = MsgTicker("...wrote data for %d points.")
                numFilteredPts = 0
                for ptIndex, pt in enumerate(self.mesh):
                    dataCols = []

                    for frameId, data in self.frameData.iteritems():
                        self._appendVal(dataCols, FieldClass,
                                        data[fieldName][ptIndex])
                    if self.filterEmpty and all(x is None for x in dataCols):
                        numFilteredPts += 1
                    else:
                        output.writerow(pt+dataCols)
                    ticker.msg(ptIndex+1)
                del output
                del ticker
                if numFilteredPts>0:
                    msg("Filtered out %d points without data." % numFilteredPts)
                msg("Finished writing result file %s" % (fileName))

        if self.outputType == "CsvPerFieldRowFrColPtName":
            # constant fields currently ignored
            # Gero suggest as future extension: Treat them like a separate
            # frame

            try:
                ptNames = self.constData["ptName"]
            except KeyError:
                raise KeyError(
                    "FramesFieldsOnUnstructuredPoints.closeCsvOutput(): There"
                    " has to be a constant field ptName for outputType"
                    " CsvPerFieldRowFrColPtName to work. Supply it by means of"
                    " the storeConstantData-method.")

            # write framedata files, one for each field
            for FieldClass in self.fieldClasses:
                fieldName = FieldClass.fieldName
                componentNames = FieldClass.getComponentNamesList()

                # open output file
                fileName = self.fileNameTemplate % fieldName
                output = csv.writer(open(fileName, "wb"))

                # prepare the head line (same for each file)
                headline=["frame number"]

                if len(componentNames)==1:
                    headline.extend(ptNames)
                else:
                    headline.extend(
                        ("%s-%s" % (ptName, cn))
                        for ptName in ptNames for cn in componentNames)

                # write headline to output file
                msg("Writing results for field %s to %s"
                    % (fieldName, fileName))
                output.writerow(headline)

                # write data lines for frames
                if len(componentNames)==1:
                    for frameId, data in self.frameData.iteritems():
                        frameName = self.frameIdToName(frameId)
                        output.writerow( [frameName] + data[fieldName] )
                else:
                    for frameId, data in self.frameData.iteritems():
                        frameName = self.frameIdToName(frameId)
                        dataCols = [frameName]
                        dataCols.extend(
                            x for val in data[fieldName] for x in val)
                        output.writerow( dataCols )

                # constant values as separate rows at the end
                if len(self.constData)>1:
                    # write empty line as separator
                    output.writerow(["" for x in headline])

                    for constFieldName, cField in self.constData.iteritems():
                        # ignore ptName: already in the header line
                        if constFieldName=="ptName":
                            continue
                        cFieldCompNames = cField.getComponentNamesList()
                        if len(cFieldCompNames)==1:
                            dataCols = [constFieldName] + [
                                x for val in cField
                                for x in ([val]+[""]*(len(componentNames)-1))]
                            output.writerow(dataCols)
                        else:
                            for i, cFieldName in enumerate(cFieldCompNames):
                                dataCols = [cFieldName] + [
                                    x for val in cField
                                    for x in ([val[i]]
                                              +[""]*(len(componentNames)-1))]
                                output.writerow(dataCols)

                # finished this file / field
                del output
                msg("Finished writing result file %s" % (fileName))

        if self.outputType in (
                "OneCsvRowFrColPtField", "OneCsvRowFrColField"):
            del self.output
            msg("Closed results file %s" % self.fileName)

        if self.outputType == "CsvPerFrameRowPtColFld":
            pass  # all done in self.storeFrame(), already


class FieldsOnUnstructMesh(FieldsOnTopology):
    r"""
    A class to store field data (for one or more scalar, vector or tensor
    fields) on an unstructured mesh. With methods to read and write this
    data to/from vtk files.

    Only for one time point. (equivalent to one vtk file)

    Usage:
     >>> from bae.abq_model_02 import Model
     >>> from bae.vtk_02 import FieldsOnUnstructMesh as Vtk
     >>> from bae.field_01 import Field
     >>> from math import sin
     >>> from bae.vecmath_01 import length, norm, vector_scale
     >>>
     >>> # create mesh and field data (lists)
     >>> mesh = Model().read("MyMesh.inp")
     >>> FieldT = Field.classFromPosType(
     ...     fieldName="T", position="node", dataType="scalar")
     >>> FieldU = Field.classFromPosType(
     ...     fieldName="U", position="node", dataType="vector")
     >>> fldT = FieldT( (node, (z+1) * sin(x**2+y**2))
     ...                for node, (x,y,z) in mesh.nodeCoords.iteritems() )
     >>> fldU = FieldU( mesh.nodeCoords )
     >>> for node, u in fldU.iteritems():
     >>>     if length(u)>2.0:
     >>>         fldU[node] = vector_scale(norm(u), 2.0)
     >>>
     >>> # export to vtk
     >>> vtk = Vtk(mesh=mesh)
     >>> vtk.updateFields(fldT, fldU)
     >>> vtk.toVtk("mydata.vtk")

    @ivar mesh: A L{mesh_01.Mesh} object representing the output
      mesh.
    @ivar data: An OrderedDict containing actual data fields.
    """

    shapeToVtkElType = {
        # shape : (pyvtk elem type, node index correction)
        "TET_L": ("tetra", None),
        "TET_Q": ("quadratic_tetra", None),
        "TRI_L": ("triangle", None),
        "TRI_Q": ("quadratic_triangle", None),
        "QUAD_L": ("quad", None),
        "QUAD_Q": ("quadratic_quad", None),
        "HEX_L": ("hexahedron", None),
        "HEX_Q": ("quadratic_hexahedron", None),
        "WEDGE_L": ("wedge", [0,2,1,3,5,4]),
        # "": ("", None),
        }

    def __init__(self, mesh=None, notDefinedValue=None):
        r"""

        @param mesh: A L{MeshStructuredPoints<bae.mesh_01.MeshStructuredPoints>}
          object representing the output mesh
        @param notDefinedValue: Initializer for the corresponding instance
          variable. For a description see
          L{FieldsOnTopology.notDefinedValue}.

        @Note: Please supply the mesh argument as keyword arguments
          only since a positional argument is reserved for a possible
          implementation of a copy constructor.
        """
        FieldsOnTopology.__init__(self, notDefinedValue)
        self.mesh = mesh
        self.data = OrderedDict()

    def toVtk(self, outputFileName, description="Field data",
              outputFormat='binary', writeSinglePrecision=True):
        r"""Write stored data to vtk file of the given name.

        @param description: String for the second line of the vtk data file
          (title field).

        @param outputFormat: One of "binary", "ascii", "ascii_1row".
          "ascii_1row" means scalar fields are written as one row with space as
          delimiters instead of newline characters. (This is the pyvtk default
          ascii implementation but Glenn for example needs "ascii" with newline
          as delimiter.)
        @param writeSinglePrecision: If True then all binary data will be
          exported as single precision
        """

        #--- store the mesh in an pyvtk.UnstructuredGrid - object
        nodeList = sorted(self.mesh.nodeCoords)
        nodeToIndex = dict(
            (node, idx) for idx, node in enumerate(nodeList))

        elemList = list()
        vtkElems = dict()
        for shape, elems in self.mesh.shapeEl.iteritems():
            vtkElType, nodeShift = self.shapeToVtkElType[shape]
            elems = sorted(elems)
            nodes = [[nodeToIndex[node] for node in self.mesh.elNodes[el]]
                     for el in elems]
            if nodeShift:
                nodes = [[nnn[i] for i in nodeShift]
                         for nnn in nodes]
            vtkElems[vtkElType] = nodes

            # bookkeeping
            elemList.extend(elems)

        vtkMesh = pyvtk.UnstructuredGrid(
            [self.mesh.nodeCoords[node] for node in nodeList],
            **vtkElems)

        #--- store the fields
        vtkPtData = pyvtk.PointData()
        vtkCellData = pyvtk.CellData()
        for fieldName, field in self.data.iteritems():

            if field.position == 'element':
                data = vtkCellData
                fieldAsList = [field[el] for el in elemList]
            elif field.position == 'node':
                data = vtkPtData
                fieldAsList = [field[node] for node in nodeList]
            else:
                msg("WARNING: Cannot write field %s at position %s. Not"
                    " implemented yet." % (field.name, field.position))
                continue

            if field.dataType == 'scalar':
                data.append(pyvtk.Scalars(
                    fieldAsList, name=fieldName,lookup_table="default"))
            elif field.dataType == 'vector':
                data.append(pyvtk.Vectors(fieldAsList, name=fieldName))
            elif field.dataType == 'tensor':
                fieldAsList = [ [ [x[0], x[3], x[4]],
                                  [x[3], x[1], x[5]],
                                  [x[4], x[5], x[2]] ]
                                for x in fieldAsList ]
                data.append(pyvtk.Tensors(fieldAsList,name=fieldName))

        #--- write the vtk file
        data = list()
        if vtkPtData.data:
            data.append(vtkPtData)
        if vtkCellData.data:
            data.append(vtkCellData)

        v = pyvtk.VtkData(vtkMesh, description, *data)
        v.tofile(outputFileName, format=outputFormat,
                 writeSinglePrecision=writeSinglePrecision)

    def fromVtk(self, inputFileName):
        r"""Not Implemented yet. Read field and mesh data from vtk file.

        Use FieldsOnStructPointGrid.fromVtk() as template.
        """
        raise NotImplementedError("FieldsOnUnstructMesh.fromVtk not implemented yet.")
        return self
