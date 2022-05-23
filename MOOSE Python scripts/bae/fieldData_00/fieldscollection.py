###GP:
__suggestions__ = """

current concept in fieldscollection.py
======================================

- FieldsCollection derives from dict, ordinary fields are stored in the dict
- timeless "space" fields (point names, mat codes) are stored as attributes
  in topo (variant: in a separate attribute "space") by the addVariable method
- spaceless "time" fields (frame names, ...) are stored as attributes in
  frames (to be renamed to "time" or so)
- timeless and spaceless attributes have one dimension less than ordinary
  fields, i.e. scalarfields are 1D - arrays
- is field name index consistent with ndarray and list indexing? In the sense
  that giving one single index --fld["u"]-- returns a VectorField and giving an
  index-tuple --fld[("U", "V")] returns a FieldsCollection?


single dict concept
===================

- fields / values attribute: dict {fieldName:ndarray}. I prefer: "values"
- timeless attributes (point names, mat codes) have shape (n,1) are stored in
  fields dict as well
- spaceless attributes (frame names, ...) have shape (1,m) are stored in
  fields dict as well
- topo attribute has __getitem__ and supports numpy-style slicing
- name conventions for frameIds, frameNames (and other special timeless or
  spaceless attributes to be used as keys in self.fields
- self.shape gets updated from newly added fields (ndarrays) if
  it's (n,1) or (1,m) or not yet assigned (0,0).
- self.shape is being checked agains newly added fields
- __getitem__ always returns a FieldsCollection.
- FieldsCollection can also serve as single field container (at least to
  acertain degree)
- We could implement arithmetic operations that work if there is only one
  field in self.fields (self.values ?) with shape==self.shape --i.e.
  field.shape==(n,m) in general--. The operation would only apply to this
  field. Problem: in case of timeless field data this field data can't be
  distinguished from timeless data like point names or so.

fancy functions:
fld[] ... fld.__getitem__, __setitem__

flds = FieldsCollection(topo, U(nxm), S(nxm), pointname(nx1), framename(1xm))
f2 = flds["U"]   =>  U(nxm), pointname(nx1), framename(1xn)
... how do we recognize pointname and framename to be taken over automatically?

f2[1:10, :]  =>  works as expected

npArr = f2.values["U"]   +> returns a numpy array


three dict concept
==================

fields            ... np.arrays
time   (spaceless-fields)
space  (timeless-fields)
topo ... separate attribute

fancy functions:
fld[] ... fld.__getitem__, __setitem__

fld.time["framenames"]
fld.space["pointname"]  ... np.arrays
fld.topo.getPoints()

fld.space["centr"] = fld.topo.getPoints()

"""
###GP: end of suggestions

import sys
import numpy as np

from bae.log_01 import msg
from fields import (ScalarField, VectorField, TensorField)
from topo import (PointCloud, StructuredGrid)
from frames import Frames

from helper import (getBaseDataType, )

#check module is loaded by epydoc-build
#--> _epyDocBuild=True to disable decorators
if 'epydoc' in sys.modules:
    _epyDocBuild = True
else:
    _epyDocBuild = False


implementedTopos = [PointCloud, StructuredGrid]


class FieldsCollection(dict):
    """Structure to store and manage time dependent scalar, vector
    and tensor fields which share the topography and time points.

    Structure of FieldsCollection:
     - topo:   topo attribute of stored fields is a reference
     - frames: frames attribute of stored fields is a reference
     - keys:   field names
     - values: field objects of type L{ScalarField}, L{VectorField} or
               L{TensorField}
    """

    def __init__(self, **kwargs):
        """
        Field-objects or/and numpy arrays have to be passed as
        keyword arguments. The keywords 'topo' and 'frames' are exclusively
        for topo- and frames-objects.

        Usage:
        ======
         >>> f = FieldsCollections(
         ...         topo=myTopo, frames=myFrames,
         ...         aField=aScalarFieldObject,
         ...         bField=aNumpyArrayThatFits,
         ...         )

        """
        self.reshape = Reshape(self)

        ## Topo-object from keywordArguments
        topo = kwargs.pop('topo', PointCloud())
        if any([isinstance(topo, t) for t in implementedTopos]):
            self.topo = topo
        else:
            raise ValueError("KeywordArgument 'topo' must be a"
                             " Topo-object object")


        ## Frames-object from keywordArguments
        frames = kwargs.pop('frames', Frames())
        if isinstance(frames, Frames):
            self.frames = frames
        else:
            raise ValueError("KeywordArgument 'frames' has to be "
                             "a Frames object")

        self.pointsTol = kwargs.pop('pointTol', 1E-3)


        ## storing fields
        for fieldName, vals in kwargs.iteritems():
            self[fieldName] = vals


    def _checkConsistency(self, fieldNames=None):
        """
        Checks if all (fieldNames=None) or specified fields are
        consistent with stored topo and frames object
        """

        if len(self) == 0:
            return

        if fieldNames is None:
            fieldNames = self.keys()
        elif isinstance(fieldNames, str):
            fieldNames = [fieldNames, ]

        notFoundFields = set(fieldNames) - set(self.keys())
        if notFoundFields:
            raise KeyError("The following fields can't be found: %s"
                           % ", ".join(notFoundFields))


        for name, field in self.iteritems():
            try:
                field._checkConsistency()
            except Exception as e:
                raise ValueError(
                    "There was an Error in field %s:\n%s" % (name, e))


#{Set and get fields
    def __setitem__(self, name, field):
        """ This overrides the dictionary __setitem__ function, so you can do
        the following:

            >>> f = Fields(topo=myTopo, frames=myFrames, someFieldsInADict)
            >>> # set a new field from numpy array
            >>> f['viaNumpy'] = np.zeros((nPt,nFr))
            >>> # set a new field from existing field
            >>> f['viaField'] = anExistingField

        """
        if name in self.keys():
            raise KeyError("A field named '%s' already exists" % name)

        fieldTypes = {2: ScalarField, 3: VectorField, 4: TensorField}

        if isinstance(field, np.ndarray):
            # input is numpy-ndArray
            try:
                fieldType = fieldTypes[len(field.shape)]
            except KeyError:
                raise ValueError('The field has no valid shape (%s)' %
                                 str(field.shape))

            field = fieldType(name, field, frames=self.frames, topo=self.topo)

        elif any([isinstance(field, t) for t in fieldTypes.values()]):
            # input is Scalar-, Vector- or Tensor-Field
            if not field.topo == self.topo:
                # no additional topo-comparison implenented yet
                msg("WARNING! The Topo-object of field to append may "
                    "differ. It will be overridden by the Topo-object "
                    "of the Fields-class. So be carefull!")

            if not field.frames == self.frames:
                # no additional topo-comparison implenented yet
                msg("WARNING! The Frames-object of field to append may "
                    "differ. It will be overridden by the Frames-object "
                    "of the Fields-class. So be carefull!")

            field.topo = self.topo
            field.frames = self.frames
            field._checkConsistency()

        else:
            raise ValueError("The field has to be a numpyArray or type of %s."
                             % str(fieldTypes.values()))

        # self[name] = field #can't be used since __setitem__ is overridden
        dict.__setitem__(self, name, field)


    def __getitem__(self, item):
        """ This enriches the dict __getitem__ fuction with the functionality
        of the getSubSet function.

        See also L{fields.ScalarField.__getitem__}.

        Usage:
        ======
            >>> field =  myFields['a'] # returns stored Field object
            >>> single = myFields[['a',]] # returns a Fields object
            >>> multible = myFields[['a', 'b']]
            >>> subPtAll = myFields[:, 10:15]
            >>> subPtSome = myFields[['a', 'b'], 10:15]
            >>> subFrAll = myFields[:, :, 1:3]
            >>> subAbr = myFields[['a', 'b'], 10:15, 1:3]

        Note that indexing of components is not allowed on Fields-level.
        Note that field names have to be lists not tuples.
        """
        if isinstance(item, str):
            # string
            return dict.__getitem__(self, item)

        elif (isinstance(item, list)
              and all([isinstance(it,str) for it in item])):
            # list of strings
            return self._getSubSet(fieldNames=item)

        elif isinstance(item, tuple):
            # parse item and use getSubSet
            try:
                fieldItem = item[0]

                # allow :-indexing Frames[:,ptIdx,frIdx]
                if fieldItem == slice(None, None, None):
                    fieldItem = None

            except IndexError:
                fieldItem = None

            try:
                topoItem = item[1]
            except IndexError:
                topoItem = None

            try:
                frameItem = item[2]
            except IndexError:
                frameItem = None

            return self._getSubSet(fieldNames=fieldItem,
                                   pointIdx=topoItem,
                                   frameIdx=frameItem)

        else:
            msg('Failed to parse index(es):\n %s' % str(item))
            raise ValueError("Can't parse the requested index (above). "
                             "Allowed indexMethods for field: "
                             "singleField <- fields[str], "
                             "multFields <- fields[list(str)], "
                             "subSet(general) <- fields[fNames, ptIdx, frIdx]")


    def update(self, *args, **kwargs):
        """This overrides the dictionary update function, so you can do
        the following:

            >>> f = Fields(topo=myTopo, frames=myFrames, someFieldsInADict)
            >>> # set a new field from numpy array
            >>> f.update({#numpy arrays or fields that fit nFrames and nPoints
            ...           'npScalar' : np.zeros((nPt,nFr)),
            ...           'npVector' : np.zeros((nPt,nFr,3)),
            ...           'yetAnotherFieldObject' : aField,
            ...           })

        """
        if len(args) > 1:
            raise TypeError("update expected at most 1 arguments, "
                            "got %d" % len(args))

        other = dict(*args, **kwargs)

        for key, value in other.iteritems():
            self[key] = value

    def _getSubSet(self, fieldNames=None, pointIdx=None, frameIdx=None):
        """
        """
        if fieldNames is None:
            fieldNames = self.keys()

        elif isinstance(fieldNames, str):
            fieldNames = [fieldNames,]

        if pointIdx is None:
            pointIdx = slice(None)

        if frameIdx is None:
            frameIdx = slice(None)

        values = {}
        for name in fieldNames:
            d = self[name]
            values[name] = d[pointIdx, frameIdx].values

        subSet = FieldsCollection(topo=self.topo[pointIdx],
                                  frames=self.frames[frameIdx],
                                  **values)
        return subSet
#} #End Set and get Fields


#{Operations on/using all stored fields
    def filterEmptyPoints(self, filterValue=0, inplace=False,
                          returnIndex=True):
        '''Filters empty pointData (e.g. sampling points above topography).
        A point will be removed if local values of each field* equals the
        filterValue. (* If filterValue is numeric/string, only numeric/string
        fields will be checked.)

        @param filterValue: point haveing these values will be excluded

        @param inplace: if True, point will be removed from all stored fields

        @param returnIndex: if True: pointIndex will be returned, if False:
            a new fields-object will be returned
        '''

        fType = getBaseDataType(filterValue)
        fieldKeys = [key for key, data in self.iteritems()
                     if fType == getBaseDataType(data.values)]

        nPts = self.topo.coords.shape[0]
        mask = np.ones((nPts, len(fieldKeys)), dtype=bool)
        for ii, fieldKey in enumerate(fieldKeys):
            idx = (self[fieldKey].values == filterValue)
            mask[:,ii] = idx.reshape(nPts,-1).all(axis=1)

        mask = ~mask.all(axis=1)

        if inplace:
            self.topo = self.topo[mask]
            for field in self.values():
                field.values = field.values[mask,:]
                field.topo = self.topo

        if returnIndex:
            return mask
        else:
            vals = dict((val.values[ii,:]) for name, val in self.values)
            newField = FieldsCollection(topo=self.topo[ii],
                                        frames=self.frames[ii],
                                        **vals)
            return newField
#} #End Operations on/using all fields

#{Read and Write
    ### PointClouds
    def writePvfCsv(self, csvFileName, **writerArgs):
        """
        Saves the fieldData in following format::
            x, y, z, var, frmId1, frmId2, ...

        to csv file

        @param csvFileName: filename of csv-file

        @kwarg writerargs: keywordArguments will be passed to writer
            (if present pandas.DataFrame.to_csv or elsewise numpy.savetxt)

        @note: just a wrapper for bae.fieldData_00.readerWriter.writePvfCsv
        """
        from readerWriter import writePvfCsv

        writePvfCsv(self, csvFileName, **writerArgs)


    @staticmethod
    def fromPvfCsv(csvFileName, varNames=[], frameList=[],
                   skipColumnList=[],
                   pointCoordList=[], pointIdList=[], pointIdxList=[],
                   pointsTol=1E-3,
                   **readerArgs):
        '''Reads data from  PointVariableFrame-ordered csv

        Example:
        ========
        >>> sLabels = ['S11','S22','S33','S23','S13','S12']
        >>> uF = FieldContainer.fromPvfCsv(csvName,varNames=slabels,
        ...                         frameList=['F000','F014',], delimiter=';')

        @param csvFileName: full path to csv-file

        @param varNames: list of variable names to read (default=[]: read all)

        @param frameList: list of frame names to read (default=[]: read all)

        @param pointCoordList: pointFilter -list of pointCoordinats
            [(x,y,z), ...] to be filtered (within self.pointsTol)

        @param pointIdList: pointFilter - list of pointIds to be filtered

        @param pointIdxList: pointFilter - list of point indices to be filtered

        @param readerArgs: additional reader arguments. See L{npReadCsv}

        @note: if points specified in pointFilters are not in csv-data this
            will NOT raise an Error/Warning

        @returns: Fields-container with topo is of type PointCloud

        @note: just a wrapper for bae.fieldData_00.readerWriter.readPvfCsv
        '''
        from readerWriter import readPvfCsv

        fields = readPvfCsv(csvFileName, varNames=[], frameList=[],
                            skipColumnList=[],
                            pointCoordList=[], pointIdList=[], pointIdxList=[],
                            pointsTol=1E-3,
                            **readerArgs)

        return fields


    ### Rectlinear-VTK-Grids
    @staticmethod
    def fromVtks(vtkNames, namesToFrames=None, fieldNames=None,
                 vtkType='structuredGrid'):
        r"""Creates Fields-Container from  a set of structured-grid vtk files.

        Usage:
        ======
         >>> from fields_00 import Fields
         >>> fields =  Fields.fromVtks("mydata.vtk",
         ...                           namesToFrames=None,
         ...                           vtkType='structuredGrid')

        This is just a wrapper for
        L{fields_00.readerWriter.fieldsFromStructuredPointsVtk}. So see here
        for more examples.

        @param fileNames: can be:
            - a single vtk-fileName
            - a list of vtk-fileNames
            - a fileNamePattern (eg. 'myVtkFiles_*.vtk')

        @param namesToFrames: specifies the 'translation' of fileNames to
          frameIds:
           - None creates a frameId according to position in fileNames-list
           - 'default' parses default frame-named vtks like 'myData_F001.vtk'
           - userdefined function that returns frameId,frameName from fileName.
             e.g.:
              >>> namesToFrams = lambda fName: (doFId(fName), doFName(fName))

        @returns: Fields-container object on StructuredGrid topography
        """
        from readerWriter import fieldsFromStructuredPointsVtk
        if vtkType.lower() == 'structuredgrid':
            fields = fieldsFromStructuredPointsVtk(
                vtkNames,
                namesToFrames=namesToFrames,
                fieldNames=fieldNames)
        else:
            raise ValueError("Can't understand vtkTyp %s" % str(vtkType))

        return fields


    def toVtks(self, filePrefix, framesToNames='default',
               vtkType='structuredGrid', componentsAsScalarFields=True):
        """Stores data from Fields-container on a StructuredGrid topography
        to a (frame-)series of vtk-files.

        This is just a wrapper for
        L{fields_00.readerWriter.fieldsToStructuredPointsVtk}. So see here
        for examples.

        @param fields: Fields-container on a StructuredGrid topography

        @param fileNamePrefix: frame identifiers will be appended
            automatically

        @param framesToNames: specifies the frame identifier:
            - None - frameIndex will be appended (e.g. prefix + '_003')
            - 'default' - default frameNaming (e.g. prefix + '_F2003')
            - userfunction - will be mapped to frameIds

        @param componentsAsScalarFields: if True all components of Vector-
            and TensorFields will stored as single ScalarFields (needed
            for Voxler)
        """
        ###ok: derive vtkType from self.topo!
        ###ok: derive default framesToNames from self.frames

        ###maybe: Can we split a class definition over multiple files in Python?
        ###maybe: If so I'd suggest to not have a separate function but define
        ###maybe: FieldsCollection.toVtks and .fromVtks in a separate file.
        ###maybe: We could define a base class for FieldCollection that has all
        ###maybe: the IO related methods in module readerWriter.
        ###ok: what about the slow vtk import? We only need vtk for I/O.
        ###ok: move all vtk import s into functions. But have a comment at the
        ###ok: top of the module
        from readerWriter import fieldsToStructuredPointsVtk

        if vtkType.lower() == 'structuredgrid':
            vtkFiles = fieldsToStructuredPointsVtk(self, filePrefix,
                         framesToNames=framesToNames,
                         componentsAsScalarFields=componentsAsScalarFields )
        else:
            raise ValueError("Can't understand vtkTyp %s" % str(vtkType))

        return vtkFiles


#} #End Read and write

###maybe: remove class Reshape and convert to reshapeToScalars, ... methods of
###maybe: FieldsCollections?

def _passSelfFields(func):
    """If self.fields is None, assume first argument to be a Fields-object"""
    if _epyDocBuild:
        return func

    def func_wrapper(self, *args, **kwargs):
        if self.fields is not None:
            return func(self.fields, *args, **kwargs)
        else:
            return func(*args, **kwargs)
    return func_wrapper


class Reshape(object):
    '''Class to rearrange and/or modify Scalar-, Vector- and/or TensorFields
    stored in a Field-object.

    This will be an attribute of Field-class: Field.reshape

    Usage:
    ======
     >>> # returns sequence of ScalarFields
     >>> ff = myFields.reshape.toScalars('aVector')
     >>> # stores ScalarFields
     >>> myFields.reshape.toScalars('aVector', store=True)
     >>> # stores new VectorField
     >>> myFields.reshape.toVector(['u1', 'u2', 'u3'], store=True)
    '''

    def __init__(self, *args):
        """Initialization can be done in two ways:
            >>> reA = Reshape()
            >>> scalars = reA.toScalars(myFields, 'myTensor')
            >>> # or set Fields in __init__
            >>> reB = Reshape(myFields)
            >>> scalars = reA.toScalars('myTensor')
        """
        if len(args) > 1:
            raise ValueError('Init of Reshape accepts up to one Fields-object')

        if len(args) == 1:
            if not isinstance(args[0], FieldsCollection):
                raise ValueError('Init of Reshape excepts up to Fields-object')

            self.fields = args[0]
        else:
            self.fields = None

    @_passSelfFields
    def toScalars(fields, fieldName, store=False, inplace=False):
        """Disassembles a Tensor- or VectorField to a series of ScalarFiels
        holding the components.

        Usage:
         >>> components = myField.toScalars('T')

        @param fieldName: key of Field

        @param store: ScalarFields of components will be stored in self
            (default = False)

        @param inplace: removes initial Field from self (default = False)
            if True kwarg store will be overridden as True

        @returns: dictionary holding ScalarFields of components
        """
        field = fields[fieldName]

        if type(field) == ScalarField:
            # msg('%s is already a ScalarField.' % fieldName)
            return {fieldName: fields,}

        rtnFields = field.toScalars()

        if store or inplace:
            existent = set(fields.keys()).intersection(set(rtnFields.keys()))
            if existent:
                raise KeyError('The following fields are already defined. '
                               'Remove them first or assign them manually.\n%s'
                               % ', '.join(existent))

            fields.update(rtnFields)

        if inplace:
            # remove initial dataset from fields
            fields.pop(fieldName)

        return rtnFields

    @_passSelfFields
    def toVector(fields, fieldNames, vectorName='', inplace=False):
        """Reshapes a sequence of ScalarFields to a VectorField

        @param fieldNames: sequence* of fieldNames used to create the new
             VectorField (*to use all VectorField-Functions, it should be 3)

        @param inplace: if True a fields in fieldNames will be removed from self

        @param vectorName: name of the new VectorField in self, if empty
            (default) the VectorField will not be appended to self

        @returns: VectorField
        """
        fields._checkConsistency(fieldNames)

        newDim = fields[fieldNames[0]].values.shape[:2] + (len(fieldNames),)
        newValues = np.empty(newDim)

        try:
            for kk, field in enumerate(fieldNames):
                newValues[:,:,kk] = fields[field].values
        except ValueError as e:
            raise ValueError("The following Error occured. Did you try to "
                             "stack any vector- and/or tensorlike fields?\n%s"
                             % str(e))

        if vectorName:
            # append to fields
            if vectorName in fields.keys():
                raise KeyError(
                    "A field named %s already exists. Remove first or assign"
                    " manually." % vectorName)

            vField = VectorField(vectorName, newValues,
                                 topo=fields.topo,
                                 frames=fields.frames)
            fields[vectorName] = vField

        if inplace:
            for field in fieldNames:
                # remove initial dataset from fields
                fields.pop(field, None)

        return newValues

    @_passSelfFields
    def toTensor(fields, fieldNames, order='A', inplace=False, tensorName=''):
        """Reshapes 3, 6 or 9 ScalarFields to a TensorField

        Example:
         >>> T = myFields.toTensor(['S_1', 'S_2', 'S_3'], order = 'D')

        @param fieldNames: sequence of 3, 6 or 9 fieldNames (see order) used to
            create the new TensorField
        @param order: The order of varLabels determines its position in 3x3-
            tensor ('A'==ABAQUS, default)
             1. diagonal
              - 'D'  [0,1,2]       --> diag([0,1,2])
             2. symetric ordering
              - 'A'  [0,1,2,3,4,5] --> [[0 3 5] [3 1 4] [5 4 2]]
              - 'RS' [0,1,2,3,4,5] --> [[0 1 2] [1 3 4] [2 4 5]]
              - 'CS' [0,1,2,3,4,5] --> [[0 1 3] [1 2 4] [3 4 5]]
              - 'V'  [0,1,2,3,4,5] --> [[0 5 4] [5 1 3] [4 3 2]]
             3. non-symetric ordering:
              - 'R'  [0,1,2,3,4,5,6,7,8] --> [[0 1 2] [3 4 5] [6 7 8]]
              - 'C'  [0,1,2,3,4,5,6,7,8] --> [[0 3 6] [1 4 7] [2 5 8]]

        @param inplace: if True a fields in fieldNames will be removed from self

        @param tensorName: name of the new TensorField in self, if empty
            (default) the TensorField will not be appended to self

        @returns: TensorField
        """
        order = order.upper()
        if len(fieldNames) == 3:
            order = 'D'
        elif len(fieldNames) == 6:
            if order in ['R','RS','C','CS']:
                order = 'RS'
            elif order in ['C','CS']:
                order = 'CS'
            elif not order == 'V':
                e = "%s  no valid type for 'order' of nonsym tensor" % order
                raise ValueError(e)

        elif len(fieldNames) == 9:
            if not (order in ['R','C']):
                e = "%s is no valid type for 'order' of nonsym tensor" % order
                raise ValueError(e)

        else:
            e = ("A tensor needs 3 (diag), 6 (symetrical) or 9 scalar "
                 "fields to be defined. You passed %d." % len(fieldNames))

        #rearrange varLabels with respect to their position in tensor
        if order == 'D':
            srcIdx = [0,1,2]
        elif order == 'A':                  # ABAQUS symmetric
            srcIdx = [0,3,5,3,1,4,5,4,2]
        elif order == 'RS':                 # rowwise symmetric
            srcIdx = [0,1,2,1,3,4,2,4,5]
        elif order == 'CS':                 # columnwise symmetric
            srcIdx = [0,1,3,1,2,4,3,4,5]
        elif order == 'V':                  # Voigt symmetric
            srcIdx = [0,5,4,5,1,3,4,3,2]
        elif order == 'R':                  # rowwise not symmetric
            srcIdx = [0,1,2,3,4,5,6,7,8]
        elif order == 'C':                  # columnwise not symmetric
            srcIdx = [0,3,6,1,4,7,2,5,8]

        fieldNames = np.array(fieldNames)[srcIdx]  # ordered varNames

        fields._checkConsistency(fieldNames)

        nPts, nFrames = fields[fieldNames[0]].values.shape[:2]
        newDim = (nPts, nFrames, 3, 3)

        #stacking data to nPointsXnFramesX3X3 numpy array
        arrayTensor = np.zeros(newDim)

        ii = 0
        try:
            for iR in range(3):
                if order == 'D':
                    arrayTensor[:,:,iR,iR] = fields[fieldNames[ii]].values
                    ii += 1
                else:
                    for iC in range(3):
                        arrayTensor[:,:,iR,iC] = fields[fieldNames[ii]].values
                        ii += 1

        except ValueError as e:
            raise ValueError("The following Error occured. Did you try to "
                             + "stack a vector- or tensorlike field?\n"
                             + str(e))

        if tensorName:
            if tensorName in fields.keys():
                raise KeyError(
                    "A field named %s already exists. Remove first or assign"
                    " manually." % tensorName)

            tField = TensorField(tensorName,arrayTensor,
                                 topo=fields.topo,
                                 frames=fields.frames)
            fields[tensorName] = tField

        if inplace:
            for field in fieldNames:
                # remove initial dataset from fields
                fields.pop(field, None)

        return arrayTensor

    @_passSelfFields
    def anglesToVector(fields, angleNameA, angleNameB,
                       angleDef='bearing+plunge',
                       vectorName='', inplace=False):
        """ Calculates a (unit-)vector from a common 2-angle-definition
        specified by angleDef

        Example:
         >>> vU = myFields('bearing(U)', 'plunge(U)',
         ...               angleDef = 'bearing+plunge' )

        @param angleNameA: field name of a ScalarField holding the 1st angle

        @param angleNameB: field name of a ScalarField holding the 2nd angle

        @param angleDef: type of angle definition of angleA and angleB
            - 'bearing+plunge' (default)
            - 'strike+dip' (needs a check)
            - 'polar+azimut'

        @returns: VectorField holding directions as UnitVector

        @warning: angleDef = 'strike+dip' needs to be checked
        """

        from transform import line,sph2cart

        # getting angular fields
        angleA = fields[angleNameA]
        angleB = fields[angleNameB]

        # retransform to vector, depending on angular definition
        if angleDef.lower() == 'polar+azimut':
            n = np.array([np.cos(angleA) * np.sin(angleB),
                          np.sin(angleA) * np.sin(angleB),
                          np.cos(angleB)])

        elif angleDef.lower() in [ 'strike+dip', 'bearing+plunge']:
            if angleDef.lower() == 'strike+dip':
                # strike to bearing
                angleB = angleB - np.pi/2.

            # bearing+plunge to coords
            x,y,z = sph2cart( line(angleB, angleA) )
            n = np.dstack( (x, y, z) )

        # create a new vectorfield
        vField = VectorField(vectorName, n,
                             topo=fields.topo, frames=fields.frames)

        if vectorName:
            # append to fields
            if vectorName in fields.keys():
                raise KeyError('A field named %s already exists. '
                               % vectorName +
                               'Remove first or assign manually.')
            else:
                fields[vectorName] = vField

        if inplace:
            # remove initial dataset from fields
            fields.pop(angleNameA, None)
            fields.pop(angleNameB, None)

        return vField


if __name__ == '__main__':
    print('No syntax errors')
