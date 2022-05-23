#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""This module contains functions to import and export field data to and from
various file formats.

It's not intended to be used directly. Use corresponding methods of the
L{FieldsCollection<bae.fieldData_00.fieldscollection.FieldsCollection>} class:
 - L{FieldsCollection.fromVtks
     <bae.fieldData_00.fieldscollection.FieldsCollection.fromVtks>}
 - L{FieldsCollection.toVtks
     <bae.fieldData_00.fieldscollection.FieldsCollection.toVtks>}

Created on Wed Apr  3 08:52:44 2019

@author: tobias
"""
import numpy as np
import os, glob

from bae.log_01 import msg

try:
    import vtk
    from vtk.util.numpy_support import (vtk_to_numpy,
                                        numpy_to_vtk,)
    from vtk.numpy_interface import dataset_adapter as dsa
except ImportError:
    vtk = None
    dsa = None
    msg('Warning. The vtk-package (not pyvtk!) is not installed.',
        ' See code (bae\fieldData_00\readerWriter.py) for install instructions')

# install instructions
# ====================

# conda: "conda install -c anaconda vtk"
# (you may need to run anaconda-console/promt as admin)

# windows/raw-python: Not sure what to do here. Fred tried one day long and
# finally used anaconda. Any further suggestions/experiences welcome.
# just a hint: https://vtk.org/download/#latestcand

# linux/pip: "pip install --update pip
            # python -m pip install --upgrade --user vtk"
# Note: the pip-USERPATH on nukus should be set to /usr/local/pip for all
# users. So you may run the previous lines:
    # - as su (sudo ...)
    # - and change owner (should be tobias) and rights recursivelly:
        # sudo chown -R tobias /usr/local/pip
        # sudo chmod -r g+rx /usr/local/pip



from bae.log_01 import msg

from fields import (ScalarField, VectorField, TensorField)
from fieldscollection import FieldsCollection
from frames import (Frames,)
from topo import (PointCloud, StructuredGrid,)

#{Errors
class ReadError(Exception):
    r"""Can be raised by FieldsOnStructPointGrid.fromVtk() if the vtk file
    to be read contains unexpected data or structure.
    """
    pass


class WriteError(Exception):
    r"""Can be raised by FieldsOnStructPointGrid.fromVtk() if the vtk file
    to be read contains unexpected data or structure.
    """
    pass
#} #End Errors


#{StructuredGridVTK

def appendNpArrayToVtk(vtkObj, name, npArray, storeAt='PointData',
                       storeAs='scalar', forceFloat=False):
    """Appends numpy-arrays to a vtkObject

    Example:
    ========
     >>> # gridSpecs
     >>> gridPtNb = (10, 11, 12)
     >>> origin = (101,102,103)
     >>> spacing = (.1, .2, .3)
     >>>
     >>> # create numpyArray
     >>> x = spacing[0] * np.arange(gridPtNb[0]) + origin[0]
     >>> y = spacing[0] * np.arange(gridPtNb[1]) + origin[1]
     >>> z = spacing[0] * np.arange(gridPtNb[2]) + origin[2]
     >>>
     >>> xx,yy,zz = np.meshgrid(x,y,z, indexing='ij')
     >>>
     >>> ss = (xx*yy*zz)**1/3. # scalarField of shape gridPtNb
     >>>
     >>> vv = np.array([xx,yy,zz])
     >>> vv = np.moveaxis(vv,0,-1) #vectorField of shape gridPtNb + (3,)
     >>>
     >>> #setUp Vtk
     >>> myVtk = vtk.vtkStructuredPoints()
     >>>
     >>> myVtk.SetDimensions(*gridPtNb)
     >>> myVtk.SetSpacing(*spacing)
     >>> myVtk.SetOrigin(*origin)
     >>>
     >>> # append-numpy array to vtk
     >>> appendNpArrayToVtk(myVtk, 'myScalar', ss,
     ...                   storeAt='pointData',
     ...                   storeAs='scalar')
     >>>
     >>> appendNpArrayToVtk(myVtk, 'myVector', vv,
     ...                   storeAt='pointData',
     ...                   storeAs='vector')


    @param vtkObj: vtkObject to store fields in

    @param name: name of the field to be stored

    @param npArray: numpyArray of field to be stored

    @param storeAt: specifies the field will be stored as 'PointData'(default)
        'CellData' or 'FieldData'

    @param storeAs: specifies wether the field will be stored as a
        'scalar', 'vector' or 'tensor'

    @param forceFloat: store int-fields as floats also (for Voxler)

    @note: only integer- and float-fields allowed
    """

    npArray = np.asarray(npArray)

    #Set the type of field
    dtype = npArray.dtype
    if dtype.kind in 'f' or forceFloat:
        vtkType = vtk.VTK_FLOAT

    elif dtype.kind in 'ui':
        vtkType = vtk.VTK_INT

    else:
        raise ValueError('dtype %s is not supported' % str(dtype))

    #Transform numpyArray to flat or 2-D field. x,y,z-dims will be flattend
    npArray = npArray.ravel(order='F')
    if 'sca' in storeAs.lower():
        npArray = npArray.ravel(order='F')

    elif 'vec' in storeAs.lower():
        if not npArray.shape[-1] == 3:
            msg("Expects last axisDim of array to be 3.",
                " Yours is %d" % npArray.shape[-1],
                " VTKoutput will be crap.")
        npArray = npArray.reshape(-1,3, order='F')

    elif 'ten' in storeAs.lower():
        if not npArray.shape[-1] == 3:
            msg("Expects last axisDim of array to be 3.",
                " Yours is %d" % npArray.shape[-1],
                " VTKoutput will be crap.")
        npArray = npArray.reshape(-1,9, order='F')

    #Transfer to vtkArray
    valsVtk = numpy_to_vtk(
                        num_array=np.ascontiguousarray(npArray),
                        deep=True,
                        array_type=vtkType)

    valsVtk.SetName(name)

    if hasattr(vtkObj, 'VTKObject'):
        #if vtkObj is numpy_interface.dataset_adapter
        vtkObj = vtkObj.VTKObject

    # Choose topo entity to store field at
    if 'point' in storeAt.lower():
        data = vtkObj.GetPointData()
    elif 'cell' in storeAt.lower():
        data = vtkObj.GetCellData()
    elif 'field' in storeAt.lower():
        data = vtkObj.GetFieldData()
    else:
        raise ValueError("Don't understand 'storeAs=%s'"%str(storeAt))

    # Choose the array type to store field in
    if 'sca' in storeAs.lower():
        data.AddArray(valsVtk)
    elif 'vec' in storeAs.lower():
        data.SetVectors(valsVtk)
    elif 'ten' in storeAs.lower():
        data.SetTensors(valsVtk)

    #vectors = vtkDoubleArray()
    #vectors.SetNumberOfComponents(3)
    #for vec in map(tuple(vecField)):
    #    vectors.insertNextTuple(*vec)
    #vtkObj.GetPointData().SetVectors(vectors)
    #tensors = vtkDoubleArray()
    #tensors.SetNumberOfTuples(2);
    #tensors.SetNumberOfComponents(9);
    #tensors.InsertTuple9(*...)

def readVtkStructuredPoints(fileName, fieldNames=None):
    """Reads a StructuredPoints-Vtk using vtk-pythonWrapper.

    >>> from bae.fieldData_00.readerWriter import readVtkStructuredPoints
    >>> gridParameters, fields = readVtkStructuredPoints("mydata.vtk")

    @param fileName: location of vtk-file

    @param fieldNames: list of PointData-fieldnames to read. None will return
        all stored fields

    @returns: gridParameters, fields
        - gridParameters is a dict holding gridPtNb, spacing and origin
        - fields is a dictionary holding (requested) numpyArrays

    @note: stored Cell- and FieldData will be ignored
    """

    reader = vtk.vtkDataSetReader()
    reader.SetFileName(fileName)

    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.ReadAllTensorsOn()

    reader.Update()

    raw = reader.GetOutput()
    data = dsa.WrapDataObject(raw)

    fields = {}

    if data.CellData.keys():
        msg("Found CellData. This will be ignored!")

    if data.FieldData.keys():
        msg("Found FieldData. This will be ignored!")

    if not data.PointData.keys():
        msg("WARNING! No PointData found in %s." % fileName)

    ptData = data.PointData

    if not fieldNames:
        fieldNames = ptData.keys()

    for name, vals in zip(ptData.keys(), ptData.values()):
        if name in fieldNames:
            fields[name] = vtk_to_numpy(vals)
            fieldNames.remove(name)

    if fieldNames:
        msg("Can't find these fields:", str(fieldNames))
        raise KeyError('Some specified fieldNames not found.')

    topoSpecs = dict(
        gridPtNb = data.VTKObject.GetDimensions(),
        spacing  = data.VTKObject.GetSpacing(),
        origin   = data.VTKObject.GetOrigin(),
        )

    return topoSpecs, fields


def fieldsFromStructuredPointsVtk(fileNames, namesToFrames='default',
                                  fieldNames = None):
    r"""Creates FieldsCollection from  a set of structured-grid vtk files.

    Usage:
    ======
     >>> from fields_00.readerWriter import fieldsFromStructuredGridVtk as rVtk

    Single file without frame information:
     >>> fields = rVtk("mydata.vtk", namesToFrames=None)

    Multiple files without frame information -> sets frameIds according to
    position in frameList:
     >>> fields = rVtk(["mydataA.vtk","mydataB.vtk"],
     ...               namesToFrames=None)
     >>> print 'frameIds: 'fields.frames.ids
     frameIds: [[2,0],[2,1]]

    Multiple files with standard frame information in names:
     >>> fields = rVtk(["mydata_F000.vtk","mydata_F001.vtk"])
     >>> # or to get all "mydata_*.vtk"
     >>> fields = rVtk(mydata_*.vtk)


    @param fileNames: can be:
        - a single vtk-fileName
        - a list of vtk-fileNames
        - a fileNamePattern (eg. 'myVtkFiles_*.vtk')

    @param namesToFrames: specifies the 'translation' of fileNames to frameIds
        - None creates a frameId according to position in fileNames-list
        - 'default' parses default frame-named vtks like 'myData_F001.vtk'
        - userdefined function that returns frameId,frameName from fileName.
          e.g.:
           >>> namesToFrams = lambda fName: (doFId(fName), doFName(fName))

    @returns: FieldsCollection object on StructuredGrid topography
    """

    if isinstance(fileNames, str):
        # single stringinput --> parse inputFileNames
        if '*' in fileNames:
            # is a pattern --> get all files matching pattern
            fileNames = glob.glob(fileNames)

        elif os.path.isdir(fileNames):
            # is a directory --> get all vtk-files in directory
            fileNames = glob.glob(fileNames + '*.vtk')
        else:
            # is a single filename
            fileNames = [fileNames,]

    # check if files exist
    notFoundFiles =  [f for f in fileNames if not os.path.isfile(f)]
    if notFoundFiles:
        msg('!!! Warning, can not find: ' + '; '.join(notFoundFiles))
        raise ReadError('%d file(s) not found.' % len(notFoundFiles))


    ## name/specify frames
    if namesToFrames == None:
        # use index of file in list to set default frameidx
        frameNames = ['F2%03d' % d for d in range(len(fileNames))]
        frameIds = [(2,d) for d in range(len(fileNames))]

    elif namesToFrames == 'default':
        #parse default vtk-fileNames
        def parser(name):
            """Parses default frame-listed file names to frameIds:
                fileBase_F2005.vtk --> (2,5)
                fileBase_F005.vtk --> (2,5)
            """
            name, _ = os.path.splitext(name)
            base = name.split('_')[-1]
            base = base.replace('F', '')
            if len(base) == 4:
                return (int(base[0]), int(base[1:]))
            elif len(base) == 3:
                return (2, int(base))
            else:
                raise ValueError("Can't parse %s to frameIdTuple" % base)

        frameIds = map(parser, sorted(fileNames))
        frameNames = ['F%d03%d' % fId for fId in sorted(frameIds)]

    else:
        #map userspecific translation function to fileNamse
        frameIds, frameNames = map(namesToFrames, sorted(fileNames))

    nFrames = len(frameIds)

    ## store fields
    fields = {}
    # loop all files
    for ii, fileName in enumerate(fileNames):
        # read single vtk-file
        gridParam, sFields = readVtkStructuredPoints(fileName)

        # create StructuredGrid-object early to use its inbuild compare
        # function to check consistency of grids
        grid = StructuredGrid(**gridParam)

        # store first grid to compare against all others
        if ii == 0:
            refGrid = grid
        else:
            if not grid == refGrid:
                raise ValueError('The parameters for structured grid differ!')


        # loop all fields stored in vtk
        for name, arr in sFields.iteritems():
            if ii == 0:
                # init emtpy field
                newArr = np.zeros((nFrames,) + arr.shape)
                newArr = np.moveaxis(newArr, 0, 1)
                fields[name] = newArr

            # assign data to array
            fields[name][:,ii] = arr

    ## setup FieldsCollection
    if fields:
        frames = Frames(ids=frameIds, names=frameNames)
        grid = StructuredGrid(**gridParam)
        outFields = FieldsCollection(topo=grid, frames=frames, **fields)

    else:
        raise Exception('Something went wrong here - fields is emtpy.')

    return outFields


def writeVtkStructuredPoints(fileName, gridPtNb, origin, spacing, fields,
                             binary=True):
    """Stores numpyArrays to StructuredPoints-VTK

    @param fileName: name of VTK-file

    @param gridPtNb: gridDimension (nX, nY, nZ)

    @param origin: 'lower-right' point of grid

    @param spacing: equidistant spacing (dX, dY, dZ)

    @param fields: dictionary holding names and values (numpy-arrays) to store
        dictionary values may come as:
            - numpy arrays -> which are assumed to be scalar-fields than
            - (numpyArray,<arrayType>) -> where arrayType is 'scalar',
            'vector' or 'tensor'

        Numpy array can be 'ij-meshgrided' (shape = (nX,nY,nZ,...)) or with
        raveled positional axis (order='F', shape=(nX*nY*nZ,...))

    @param binary: if True (default) vtk-files will be binaries

    """

    out = vtk.vtkStructuredPoints()

    out.SetDimensions(*gridPtNb)
    out.SetSpacing(*spacing)
    out.SetOrigin(*origin)

    for name, vals in sorted(fields.iteritems()):
        if isinstance(vals, ScalarField):
            if type(vals) == ScalarField:
                storeAs = 'scalar'
            elif type(vals) == VectorField:
                storeAs = 'vector'
            elif type(vals) == TensorField:
                storeAs = 'tensor'

            vals = vals.topo.fullValuesGrid(vals,
                            emptyValues=0.,
                            rollPositionAxis=False)

        elif isinstance(vals, tuple):
            vals, storeAs = vals

        else:
            storeAs = 'scalar'

        if not isinstance(vals, np.ndarray):
            raise ValueError('Allowed inputs: Field-object, numpyArray or',
                             ' tuple(numpyArray, fieldTypeString).',
                             ' You passed %s' %  str(type(vals)))

        appendNpArrayToVtk(out, name, vals,
                           storeAt='pointData',
                           storeAs=storeAs)

    writer = vtk.vtkDataSetWriter()
    writer.SetFileName(fileName)
    if binary:
        writer.SetFileTypeToBinary()
    writer.SetInputData(out)
    err = writer.Write()
    if not err:
        raise WriteError("Writing %s failed because of <lame excuse>" %
                         fileName)

    # modify VTK version to 2.0 (for Voxler)
    ff = open(fileName, "r+b")
    ff.seek(23)
    ff.write("2.0")
    ff.close()



def fieldsToStructuredPointsVtk(fields, fileNamePrefix,
                              framesToNames='default',
                              componentsAsScalarFields=True):
    """Stores data from FieldsCollection on a StructuredGrid topography to
    a (frame-)series of StructuredPoint vtk-files.

    @param fields: FieldsCollection on a StructuredGrid topography

    @param fileNamePrefix: frame identifiers will be appended automatically

    @param framesToNames: specifies the frame identifier:
        - None - no frameIndex will be appended (raises IOError if more than
          one frame is stored
        - 'default' - default frameNaming (e.g. fileNamePrefix + '_F2003')
        - userfunction - will be mapped to frameIds

    @param componentsAsScalarFields: if True all components of Vector- and
        TensorFields will stored as single ScalarFields (needed for Voxler)

    @Note: Storing components as scalarfields requires a deepcopy of the
        fields-container. If this causes time or memory issues we'll need to
        rework this.
    """

    from copy import deepcopy

    fileNamePrefix = os.path.splitext(fileNamePrefix)[0]

    grid = fields.topo
    if not hasattr(grid, 'gridPtNb'):
        raise ValueError("FieldsCollection needs to have a valid" +
                         " StructuredGrid-topo-object")

    frameIds = fields.frames.ids
    if framesToNames is None:
        if len(frameIds) > 1:
            raise IOError('You are trying to store more than one frame'
                          ' into the same vtk-file (without a FrameId'
                          ' in fileName). To enforce this stupid behavior' +
                          ' set framesToNames=lambda(s): ""')
        postfixes = ['' for dd in frameIds]
    elif framesToNames == 'default':
        postfixes = ['_F%d%03d' % tuple(dd) for dd in frameIds]
    else:
        try:
            postfixes = map(framesToNames, frameIds)
        except Exception as e:
            print str(e)
            raise Exception('Something went wrong while mapping the',
                            ' framesToNames-function to frameIds.',
                            ' See error-output above.')

    if componentsAsScalarFields:
        # find all Vector- and TensorFields
        fieldsToFlatten = [name for name, field in fields.iteritems()
                           if not type(field) == ScalarField]

        if fieldsToFlatten:
            # To simplify the following code, reshaping will be done inplace.
            # To avoid changing the original fields-object a deepcopy is used.
            fields = deepcopy(fields)

        for name in fieldsToFlatten:
            msg('Reshaping field %s to scalar fields' % name)
            fields.reshape.toScalars(name, inplace=True)

    fileNames = []
    for iFrame, postfix in enumerate(postfixes):
        fieldsFrame = fields[:,:,iFrame]

        fileName = fileNamePrefix + postfix + '.vtk'

        writeVtkStructuredPoints(fileName,
                                 grid.gridPtNb, grid.origin, grid.spacing,
                                 fieldsFrame, binary=True)

        fileNames.append(fileName)

    return fileNames

#} #End StructuredGridVTK

###############################################################################
#{PointCloud
def writePvfCsv(fields, csvFileName, **writerArgs):
    """Save the fields in Point Variable Frame (PVF) format. This is a csv file
    with columns: x, y, z, var, frmId1, frmId2, ...

    This format will also be written by L{odbToCsvPoints<bae.utils
    .odbToPointsData_02.odbToCsvPoints>} and friends with parameter
    csvOutputType set to 'OneCsvRowPtFldColFr'.

    @param fields: L{FieldsCollection} object

    @param csvFileName: filename of csv-file

    @kwarg writerargs: keywordArguments will be passed to writer
      (pandas.DataFrame.to_csv if present or else numpy.savetxt)
      WARNING: Using this parameter compomises the portability of the
      program. Don't use in libraries.
    """

    try:
        import pandas as pd
    except ImportError:
        pd = None

    fields._checkConsistency()

    outFieldVals, outFieldNames = (), ()

    # loop all fields
    for name, field in sorted(fields.iteritems()):
        try:
            # disassemble Vector/TensorField to ScalarField and append
            tmpFields = field.toScalars(name)
        except AttributeError:
            # is already ScalarField
            tmpFields = {name : field}
        for tmpName, tmpField in sorted(tmpFields.iteritems()):
            outFieldVals = outFieldVals + (tmpField.values,)
            outFieldNames = outFieldNames + (tmpName,)

    # shortcut coordinates
    coords = fields.topo.coords

    # in case of structured grids
    if callable(coords):
        coords = coords()

    nPts = coords.shape[0]
    nFields = len(outFieldNames)

    #repeat the fieldNames for each point -->['a', 'a', 'a', ... 'z','z']
    outFieldNames = np.repeat(outFieldNames, nPts)[:,np.newaxis]

    #repeat coordinates for each field --> [pt1, pt2, ..., pt1, pt2]
    coords = np.tile(coords, (nFields,1))

    if hasattr(fields.topo, 'ids'):
        header = ['x', 'y', 'z', 'ptId', 'var'] + fields.frames.names
        # fmt = 3*['%f'] + 2*['%s'] + len(fields.frames.names)*['%f']
        ptIds = np.tile(fields.topo.ids, nFields)[:,np.newaxis]
        firstCols = np.hstack((coords, ptIds, outFieldNames))

    else:
        header = ['x', 'y', 'z', 'var'] + list(fields.frames.names)
        # fmt = 3*['%f'] + 1*['%s'] + len(fields.frames.names)*['%f']
        firstCols = np.hstack((coords, outFieldNames))

    #stacking first columns and data columns
    outData = np.hstack((firstCols, np.vstack(outFieldVals)))

    msg('### storing FieldData to %s' % csvFileName)
    if pd is not None:
        # write csv using pandas
        df = pd.DataFrame(data=outData, columns=header)
        df.to_csv(csvFileName, index=False)
    else:
        # write using numpy.savetxt
        # ToDo:
        # for specifiing output format for each column we need to transform
        # outData to a structured array.
        np.savetxt(csvFileName, outData,
                   # fmt=fmt, #outData as structured array?!
                   fmt='%s',
                   delimiter=',',
                   header=','.join(header))



def readPvfCsv(csvFileName, varNames=[], frameList=[],
               skipColumnList=[],
               pointCoordList=[], pointIdList=[], pointIdxList=[],
               pointsTol = 1E-3,
               **readerArgs):
    """Reads data from a Point Variable Frame (PVF) format csv file with
    columns: x, y, z, var, frmId1, frmId2, ...

    This format will be written by L{readPvfCsv}. Additionally you get it from
    L{odbToCsvPoints<bae.utils.odbToPointsData_02.odbToCsvPoints>} and
    friends with parameter csvOutputType set to 'OneCsvRowPtFldColFr'.

    Example:
    ========
     >>> sLabels = ['S11','S22','S33','S23','S13','S12']
     >>> uF = FieldContainer.fromPvfCsv(csvName,varNames=slabels,
     ...                         frameList=['F000','F014',], delimiter=';')

    @param csvFileName: full path to csv-file

    @param varNames: list of variable names to read (default=[]: read all)

    @param frameList: list of frame names to read (default=[]: read all)

    @param pointCoordList: point filter - get only data for points listed
      in this point filter list: a list of point coordinates [(x,y,z), ...].
      Points must match by the tolerance given by argument pointsTol.

    @param pointsTol: tolerance for the pointCoordList argument

    @param pointIdList: point filter - list of pointIds to be filtered

    @param pointIdxList: pointFilter - list of point indices to be filtered

    @param readerArgs: additional keyword arguments to be passed to the reader
      L{npReadCsv<bae.misc_01.npReadCsv>}.
      WARNING: Using this parameter compomises the portability of the
      program. Don't use in libraries.

    @note: If the filter options result in no points being actually read
      then this will NOT raise an error or warning.

    @returns: self, so you can do myField = FieldsCollection().readPvfCsv(name)
    """
    from bae.misc_01 import npReadCsv, readCsvHeader
    from topo import KDTree

    try:
        import pandas as pd
        hasPandas = True
    except ImportError:
        hasPandas = False
        if any([pointCoordList, pointIdList, pointIdxList,]):
            msg("Warning!!! OnTheFly-Filtering by points is only " +
                "supported for pandas (which you havn't installed). " +
                "You need to remove unwanted Points manually.")


    if bool(pointCoordList) + bool(pointIdList) + bool(pointIdxList) > 1:
        raise ValueError('Ambiguous setting of point-filter. ' +
                         'pointCoordList, pointIdList and ' +
                         'pointIdxList can only be used exclusively.')

    ## setUp rules for (row-)conditional reading
    ## if all(rowConditions): append row to data
    rowConditions = []
    if hasPandas:
        if bool(varNames):  # @ given varNames
            rowConditions.append(
                lambda data: np.isin(data['var'], varNames))

        if bool(pointCoordList):  # @ given coordinates
            def dist(testPt):
                return np.norm(np.asarray(pointCoordList) - testPt, axis=1)

            rowConditions.append(
                lambda data: dist(data[['x','y','z']] < pointsTol))

        if bool(pointIdList):  # @ given point ids
            rowConditions.append(
                lambda data: np.isin(data['ptId'], pointIdList))

        if bool(pointIdxList):  # @ given point indexs
            rowConditions.append(
                lambda data: np.isin(data['ptIdx'], pointIdxList))

    ## cummulate all rowContitions to one function
    if rowConditions:
        rowCondition = lambda data: all([c(data) for c in rowConditions])
    else:
        rowCondition = None

    ## specify columns to read
    delimiter = readerArgs.pop('delimiter', ',')

    header = readCsvHeader(csvFileName, delimiter=delimiter)

    iCols = ['x', 'y', 'z', 'var']

    if 'ptId' in header:
        iCols.append('ptId')

    if 'ptIdx' in header:
        iCols.append('ptIdx')

    if not frameList:
        #read all
        frameList = [h for h in header if h not in iCols]

    readCols = iCols + frameList
    readCols = [c for c in readCols if c not in skipColumnList]

    ## read data
    data = npReadCsv(csvFileName, columns=readCols,
                     rowCondition=rowCondition,
                     delimiter=delimiter, **readerArgs)

    if not data.size:
        raise ValueError("Dataset is empty. Maybe the "
                         "rowContition isn't defined correctly.")

    ## create FrameStructure
    frames = Frames()
    frames.names = frameList
    try:
        frames._defaultIdsFromNames()
    except ValueError as e:
        raise ValueError("Error while parsing FrameNames. The following " +
                         "Error was raised: \n \t%s\n" % str(e) +
                         "Maybe the csv contains extraColumns which do " +
                         "not represent FrameData. Use keywordArg " +
                         "'skipColumnList' to exclude these.")
    nFrs = len(frameList)

    ## create PointsStructure
    coords = np.vstack([data[q] for q in ['x','y','z']]).T
    nFlatPts = coords.shape[0]

    # do not get confused: PointCloud.__init__ will also build a KD-Tree
    # to check the (hopefully) unique coords-points, which we'll need to
    # identify before by using the following lines of code
    tree = KDTree.KDTree(coords)
    groups = KDTree.uniqueTreeIdx(tree, atol=pointsTol)

    idxUnique = [ii[0] for ii in groups]

    coords = coords[idxUnique,:] #to get the points in the same order as they
                           #appear in csv
    nPts = coords.shape[0]

    points = PointCloud()
    points.tol = pointsTol
    points.addVariable('coords', coords)
    try:
        pointIds = data['ptId'][idxUnique]
        points.addVariable('ids', pointIds)
    except:
        pass

    ## create Fields
    outFields = FieldsCollection(topo=points, frames=frames)

    ## create a indexLookup to unique Points
    ## use length of coords to raise IndexError on failed lookup
    lookUp = nPts * np.ones(nFlatPts, dtype=int)
    for ii, look in enumerate(groups):
        lookUp[look] = ii

    ## get uniqueFields
    fieldNames = np.unique(data['var'])

    fieldData = np.vstack([data[q] for q in frameList]).T

    ## storing fields
    for fieldName in fieldNames:
        valueIdx = data['var'] == fieldName
        idx = lookUp[valueIdx]
        vals = np.nan * np.ones((nPts, nFrs))

        vals[idx] = fieldData[valueIdx,:]

        outFields[fieldName] = vals

    return outFields
#} #End PointCloud

def _test_vtkStructuredPoints():
    from topo import StructuredGrid
    gridPtNb = (10, 11, 12)
    origin = (101,102,103)
    spacing = (.1, .2, .3)

    x = spacing[0] * np.arange(gridPtNb[0]) + origin[0]
    y = spacing[0] * np.arange(gridPtNb[1]) + origin[1]
    z = spacing[0] * np.arange(gridPtNb[2]) + origin[2]

    xx,yy,zz = np.meshgrid(x,y,z, indexing='ij')

    ss = (xx*yy*zz)**1/3.
    vv = np.array([xx,yy,zz])
    vv = np.moveaxis(vv,0,-1)

    ss = ss.ravel(order='F')
    vv = vv.reshape(-1, 3, order='F')

    # stack some frames
    sf = np.moveaxis(np.array([ss,ss + 1000]), 0, 1)
    vf = np.moveaxis(np.array([vv, vv + 1000]), 0, 1)


    topo = StructuredGrid(gridPtNb=gridPtNb, origin=origin, spacing=spacing)
    frames = Frames([(2,d) for d in range(2)])

    fields = FieldsCollection(topo=topo, frames=frames, myVector=vf,
                              myScalar=sf)
    #writeVtkStructuredPoints(fileName, gridPtNb, origin, spacing, fields)



    fieldsToStructuredPointsVtk(fields, 'someTestVtks',
                              componentsAsScalarFields=False)

    fff = fieldsFromStructuredPointsVtk('someTestVtks*.vtk')
    return fff


if __name__ == "__main__":
    print 'No syntax errors'
