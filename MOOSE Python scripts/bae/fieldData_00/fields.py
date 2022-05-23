# -*- coding: utf-8 -*-
"""Module containing classes to store
    - scalar (L{ScalarField}),
    - vector (L{VectorField}) and
    - tensor (L{TensorField}) fields
that are given on a
    - spacial structure / L{topology<bae.fieldData_00.topo>}
at a given
    - time / L{frame<bae.fieldData_00.frames>}.

General structure of a Field-object::
 Field
    |- name
    |- topo     <- L{topo<bae.fieldData_00.topo>}-object of size nPoints
    |- frames   <- L{frames<bae.fieldData_00.frames>}-object of size nFrames
    |- values   <- numpyArray of shape:
                         ScalarField: (nPoints, nFrames)
                         VectorField: (nPoints, nFrames, 3)
                         TensorField: (nPoints, nFrames, 3, 3)

The field objects behave like multidimansional arrays, the first index
represents space and the second represents time. A particular index or index
range simultaneously addresses field values, frames object and topo object.

Indices can be:
    - single integers
        >>> subField = myField[1,8] # second point, 9th frame
    - slices
        >>> subField = myField[:,:3] # all points, first 3 frames
    - index lists: second, 7th, 9th point, first and eleventh time point
        >>> subField = myField[[1,6,8], [0,10]]
    - boolean index arrays
        >>> validPoints = myField[:,0].values > 0 # returns [True, False, ...]
        >>> subField = myField[validPoints, :]

Different fields sharing the same time and spacial structure can be stored
together in a L{FieldsCollection<bae.fieldData_00.fieldscollection>} object.
"""

###GP: Do we need separate classes ScalarField, VectorField? Why?
###GP: - VectorField and TensorField have special methods
###GP:   and Tobias prefers fld.norm() and fld.eigen() (implemented as methods)
###GP:   rather than norm(fld) ...

import sys
import numpy as np
from copy import deepcopy

from bae.log_01 import msg

from .topo import (PointCloud, StructuredGrid)
from .frames import Frames
###ok: from .time import FrameList
from .commonfieldplots import (CommonScalarPlots,
                              CommonVectorPlots,)
from .helper import (getBaseDataType, )


#check module is loaded by epydoc-build
#--> epyDocBuild=True to disable decorators
if 'epydoc' in sys.modules:
    _epyDocBuild = True
else:
    _epyDocBuild = False


implementedTopos = [PointCloud, StructuredGrid]

#{ Helper functions and decorators
def _passValuesToOperator(func):
    '''This wrapper (decorator) is needed to pass the values-attribute
    instead of the Scalar/Vector/TensorField itself to mathematical
    operators.
    '''
    if _epyDocBuild:
        return func

    def func_wrapper(self, other):
        if hasattr(other, 'topo'):
            other = other.values
        return func(self, other)
    return func_wrapper

#} # End Helper functions and decorators

class ScalarField(object):
    """
    Structure to store a time dependent scalar field on a topography.

    The time-dependent scalar field will be stored as a numpy array with
    the shape of (nPoints, nTimePoints) in ScalarField.values.
    The time information can be stored in a FrameList-object ScalarField.times
    having the same length nTimePoints.

    It allows simultaneously indexing of values, frame- and topo-object.
    For example to get all points starting with the third; first two frames:
     >>> myField[2:,:2]

    """

    def __init__(self, *args, **kwargs):
        """
        Name and values can be passed as positional arguments or as keyword
        arguments

        Examples:
         >>> A = ScalarField('A', vals)
         >>> B = ScalarField(name='A', values=vals)
         >>> C = ScalarField('A', values=vals)

        FrameList- and Topo-object has to be passed as keywordArg. If no frames
        structure is passed, an empty one will be created. If no topo structure
        is passed an empty PointCloudObject will be created.

        Examples:
         >>> nFrames, nPoints = 5, 20
         >>> myFrames = Frames(ids=[(2,ii) for ii in range(nFrames)])
         >>> myTopo = PointCloud(np.random.randn(3,nPoints))
         >>> vals = np.random.randn(nPoints,nFrames)
         >>> D = ScalarField('A', vals, frames=myFrames, topo=myTopo)

        Additional keyword arguments:
         - enablePlot: True (default) enables/appends commonScalarPlots-Object
        """

        #name from positional/keywordArguments
        if len(args):
            self.name = args[0]
        else:
            self.name = kwargs.pop('name', None)

        if self.name is None:
            raise ValueError("You need to define a fieldName")

        ## values from positional/keywordArguments
        if len(args) == 2:
            self.values = np.asarray(args[1])
        else:
            self.values = kwargs.pop('values', np.array([]))

        ## Topo-object from keywordArguments
        topo = kwargs.pop('topo', PointCloud())
        if any([isinstance(topo, t) for t in implementedTopos]):
            self.topo = topo
        else:
            raise ValueError("KeywordArgument 'topo' must be a Frames object")


        ## Frames-object from keywordArguments
        frames = kwargs.pop('frames', Frames())
        if isinstance(frames, Frames):
            self.frames = frames
        else:
            raise ValueError("KeywordArgument 'frames' has to be " +
                             "a Frames object")

        ## additional keywordArguments
        if kwargs.get('enablePlot', True):
            self.plot = CommonScalarPlots(self)

        self._checkConsistency()


    def _slicedName(self, componentIndexs):
        '''Helper function to get proper names for component-indexed fields.
        Parses the 'old' name of Field and returns the current components
        of indexed Field.

        See L{__getitem__}-Naming convention for details.

        @componentIndexes: list holding the requested component indices/slices
        '''
        # no name set --> return None
        if self.name is None:
            return [self.name]

        # no indexing on components --> unchanged name
        if not componentIndexs:
            return [self.name]

        inShape = self.values.shape

        if len(inShape) == 4 and len(componentIndexs) == 1:
            componentIndexs.append(slice(None,None,None))

        # quick check to get the dimensions after indexing
        newDims = [len(np.empty(inShape[ii+2])[cpIdx])
                   for ii,cpIdx in enumerate(componentIndexs)]

        # see which dimensions will change through indexing
        dimRemains = [self.values.shape[ii+2] == newDim
                      for ii,newDim in enumerate(newDims)]

        # no change of componennts --> unchanged name
        if all(dimRemains):
            return [self.name]

        # possible structures of initial name:
        # cases:
        # 1,4) baseName           # wasn't indexed yet
        # 2) baseName[0,...]    # vector was indexed in some components
        # 3) baseName[1][0,...] # tensor was indexed in some components of
        #                    # one direction
        # *) baseName[0,...][0,...] # should never happen because
        #                        # TensorDirections have to be assigned
        #                        # direction by direction


        # get the basename (everything before '[')
        # works even if no '[' is in self.name
        name = self.name.split('[',1)[0]

        # get the first old indexSpecifier
        try:
            slice1 = self.name.split('[', 1)[1].split(']')[0]
        except:
            slice1 = None

        # get second old indexSpecifier
        if not '][' in self.name:
            slice2 = None
        else:
            try:
                slice2 = self.name.rsplit('[', 1)[1].rsplit(']')[0]
            except:
                slice2 = None

        newSlice2 = np.array([])
        # case 1 - was unsliced vector
        if len(inShape)==3 and slice1 is None and slice2 is None:
            newSlice1 = np.arange(inShape[2])[componentIndexs[0]]

        # case 2 - was sliced vector
        elif len(inShape)==3 and slice1 is not None and slice2 is None:
            oldSlice1 = np.array(map(int, slice1.split(',')))
            newSlice1 = oldSlice1[componentIndexs[0]]

        # case 3 - was a (vector)sliced tensor
        elif len(inShape)==3:
            newSlice1 = np.array(map(int, slice1.split(',')))
            oldSlice2 = np.array(map(int, slice2.split(',')))
            newSlice2 = oldSlice2[componentIndexs[0]]

        # case 4 - was a tensor field
        # by design a tensorfield can only be indexed in the first component
        # and only with a singular integer (see case *)
        elif len(inShape)==4:
            newSlice1 = np.arange(inShape[2])[componentIndexs[0]]
            newSlice2 = np.arange(inShape[3])[componentIndexs[1]]

        return [name, newSlice1.tolist(), newSlice2.tolist()]



    def __getitem__(self, item):
        """This __getitem__ function allows (sub-) indexing of the field
        including values, topo-object and frames-object simultaneously.

         >>> # myField is a ScalarField on 10 points over 5 frames
         >>> print 'myField.values.shape: ', myField.values.shape
         myField.values.shape: (10,5)
         >>> print 'len(myField.topo): ', len(myField.topo)
         len(myField.topo): 10
         >>> print 'len(myField.frames): ', len(myField.frames)
         len(myField.topo): 5
         >>> # now pick the first 3 points of the last frame
         >>> mySub = myField[:3, -1]
         >>> print 'mySub.values.shape: ', mySub.values.shape
         mySub.values.shape: (3,1)
         >>> print 'len(mySub.topo): ', len(mySub.topo)
         len(mySub.topo): 3
         >>> print 'len(mySub.frames): ', len(mySub.frames)
         len(mySub.frames): 1

        B{Allowed index types} are:
            - integer values: myField[1,2]
            - slices: myField[1:,:-4]
            - index flat lists or np.arrays: myField[[1,3,5], [0,1]]
            - flat boolean arrays with same length of indexed dimension:
            myField[[True,False,....]]

        Please note that multidimensional indexing is slightly different from
        numpy-syntax. In numpy indexing like myField[[1,4],[3,5]] will return
        the two values myField13 and myField45. Where here it would return
        a 2Pointsx2Frames data field.

        B{Indexing of components}:
        The components of Vector- and TensorFields can be indexed as well.
        This __getitem__ function will return the appropriate FieldType and
        fieldNames as well.

        Example:
         >>> type(myVectorField[:,:,1]) # returns a ScalarField
         >>> type(myTensorField[:,:,1,0]) # returns a ScalarField
         >>> type(myTensorField[:,:,1,:]) # returns a VectorField

        B{Naming convention}:
        The name-attribute of indexed fields will be set automatically. This is
        quite tricky. If it fails or you do not like the resulting names please
        set the name-attribute manually after indexing.
         - indexing only points and/or frames will not effect the name of the
         resulting field
         - if just a single component is choosen from Tensor- or VectorFields
         the well known 'underscoreIndexes' format starting with index=1 is
         used:
             >>> print 'tensorFieldName: ', myTensor.name
             tensorFieldName: S
             >>> print 'name of tensorComponent: ', myTensorField[:,:,1,2].name
             name of tensorComponent: S_23
             >>> print 'vectorFieldName: ', myVector.name
             vectorFieldName: U
             >>> print 'name of vectorComponent: ', myVectorField[:,:,0].name
             name of vectorComponent: U_1
         - if subsets of components (like first and third) are requested these
         (zeroStarting) indices will appear in the resulting field name in
         brackets:
             >>> print 'name of subSetTens: ', myTensField[:,:,1].name
             name of subSetTens: S[1][1,2,3]
             >>> print 'name of subSetVec: ', myVectorField[:,:,[0,2]].name
             name of subSetVec: U[0,2]

        If a subset of components gets indexed once more the naming function
        will try to resolve the correct components by parsing the subset name:
         >>> subSetTens = myTensField[:,:,1]
         >>> print 'name of subSetTens: ', subSetTens.name
         name of subSetTens: S[1][0,1,2]
         >>> print 'subSetTens[:,:,0].name: ', subSetTens[:,:,0].name
         subSetTens[0].name: S_21
         >>> print 'subSetTens[:,:,[0,2]].name: ', subSetTens[:,:,[0,2]].name
         subSetTens[0,2].name: S[1][0,2]

        @note: avoid ']' or '[' in when setting initial fieldNames because
        it will confuse the automated name creation of sliced fields

        @warning: if you want to index only points (first dimension), do not use
        an index tuple because its entrys will be interpreted as indexdimension
        each
        """
        #1) prepare index input
        if isinstance(item, (int, long, list, np.ndarray, slice),):
            # only a point-item is given
            item = (item,)
        elif isinstance(item, tuple):
            pass
        else:
            raise IndexError("Can't determine the index-logic from %s"
                             % str(item))

        # 2) indexing values, topo and frames
        # fieldData to be indexed
        vals = self.values
        topo = self.topo
        frames = self.frames

        # initial field type
        if len(vals.shape) == 2:
            inType = 'scalar'
        elif len(vals.shape) == 3:
            inType = 'vector'
        elif len(vals.shape) == 4:
            inType = 'tensor'
        else:
            raise ValueError("Can not identify a field of dim %s"
                             % str(vals.shape))

        componentIndex = [] # to collect componentSlices
        # the field/field.values will be indexed squentially
        for ii, idx in enumerate(item):
            if isinstance(idx, (int,long)):
                # force single index to be iterable to preserve dimension
                idx = np.array([idx])

            # create 'complete' slice object over all indexed dimensions
            cSlice = len(item) * [slice(None,None,None),]
            # replace the sliced dimension with current index
            cSlice[ii] = idx
            vals = vals[tuple(cSlice)]

            if not vals.size:
                # return None if indexed array has no data stored
                msg('WARNING! Indexed field has no values. Returning None')
                return None

            if ii == 0: # spacial index
                topo = topo[idx]
            elif ii == 1: # time index
                frames = frames[idx]
            else: # indexing a component
                componentIndex.append(idx)

        # 3) try to resolve a proper name for indexed field
        try:
            stringData = self._slicedName(componentIndex)
            def slicesToStr(s0=[], s1=[]):
                if not s0 and not s1:
                    return ''

                ## single tensor component
                if len(s0) == 1 and len(s1) == 1:
                    return '_%d%d' % (s0[0] + 1, s1[0] + 1)

                ## single vector component
                elif len(s0) == 1 and len(s1) < 1:
                    return '_%d' % (s0[0] + 1)

                ## sliced tensor/vector
                else:
                    s = ['','']
                    if s0:
                        s[0] = '[%s]' % ','.join(map(str,s0))

                    if s1:
                        s[1] = '[%s]' % ','.join(map(str,s1))
                    return ''.join(s)

            newName = stringData[0] + slicesToStr(*stringData[1:])

        except:
            raise
            msg(("Can't parse the name of the indexed field for %s correctly "
                 % self.name) + " Rename the indexed field manualy.")
            stringData = [self.name + '_sliced']

        # 4) create returned field
        if inType == 'vector' and vals.shape[2] == 1:
            rtn = ScalarField(newName, vals[:,:,0],
                              topo=topo, frames=frames)

        elif inType == 'tensor' and vals.shape[2:] == (1,1):
            rtn = ScalarField(newName, vals[:,:,0,0],
                              topo=topo, frames=frames)

        elif inType == 'tensor' and vals.shape[2] == 1:
            rtn = VectorField(newName, vals[:,:,0,:],
                              topo=topo, frames=frames)

        else:
            rtn = self.__class__(newName, vals,
                                 topo=topo, frames=frames)
        return rtn


    def _checkValuesConsistency(self):
        '''Checks the consistency of values attribute.'''
        dims = self.values.shape
        if not len(dims) == 2:
            raise ValueError("The ScalarField.values has to be a 2D-numpy-"
                             "array of shape Npoints x Nframes")


    def _checkConsistency(self):
        """Simple dimension-check of values, points and frames"""
        if len(self.values) == 0:
            return

        # call self._checkValuesConsistency() here, so only this needs to
        # be updated for inherent VectorField and TensorField classes
        self._checkValuesConsistency()

        dims = self.values.shape

        # check topo
        try:
            self.topo._checkConsistency()

        except Exception as e:
            raise ValueError("There is a problem in the stored Topo-object: "
                             + str(e))

        lenTopo = len(self.topo)

        if lenTopo == 0:
            msg('WARNING! The Topo-Object has no (valid) data assigned')

        elif not dims[0] == lenTopo:
            raise ValueError(("The number of values (%d) does not match the " +
                              "number of points (%d)") % (dims[0], lenTopo))

        # check frames
        try:
            self.frames._checkConsistency()
        except Exception as e:
            print self.frames, self.frames.ids
            raise ValueError("There is a problem in the stored Frames-object: "
                             + str(e))

        lenFrames = len(self.frames)

        if lenFrames == 0:
            msg('WARNING! The Frames-Object has no (valid) data assigned')

        elif not dims[1] == lenFrames:
            raise ValueError(("The number of values (%d) does not match the " +
                              "number of Frames (%d)") % (dims[1], lenFrames))



class VectorField(ScalarField):
    """
    Structure to store time (frame) dependend vector fields on a topography.

    The timedependend (frames) vector field will be stored as a numpy array
    with the shape of (nPoints, nFrames, 3) in VectorField.values.
    """
    def __init__(self, *args, **kwargs):
        ScalarField.__init__(self, *args, **kwargs)

        if kwargs.get('enablePlot', True):
            self.plot = CommonVectorPlots(self)

    def _checkValuesConsistency(self):
        '''Checks the consistency of values attribute.'''
        dims = self.values.shape
        if not len(dims) == 3:
            raise ValueError("The VectorField.values has to be a 3D-numpy-"
                             "array of shape Npoints x Nframes x 3")


#{Common functions for VectorFields
    def toScalars(self, nameFun = None):
        """Returns VectorField-components as ScalarFields

        @param nameFun: namingFunction(self.name,i) -
            default: nameFun = lambda s: '%s_%d' % (s[0], s[1]+1)

        @returns: dictionary with ScalarFields of components
        """

        if nameFun is None:
            nameFun = lambda s: '%s_%d' % (s[0], s[1]+1)

        rtnFields = {}
        for i in range(self.values.shape[2]):
                cName = nameFun((self.name,i))
                rtnFields[cName] = ScalarField(cName, self.values[:,:,i],
                                               topo = self.topo,
                                               frames = self.frames)
        return rtnFields


    def norm(self):
        """Returns the Euklidic norm (length) of vector as ScalarField"""

        length = np.linalg.norm(self.values, axis=2)

        length = ScalarField('norm(%s)' % self.name, length,
                            points = self.topo, frames = self.frames)

        return length


    def vectorToAngles(self, angleDef='bearing+plunge', pointDownwards=False):
        '''Returns the tow angular field values (scalarFields) describing the
        direction of stored vectorField.

        @param angleDef: specifies the angle definition: 'bearing+plunge'
            (default), 'polar+azimut' and 'strike+dip' (needs a check!!!)

        @param pointDownwards: if True, all vectors having an positive
            z-component will be mirrored before

        @note: returned values in radians
        '''
        from numpy_geom import vector2plunge_bearing

        vec = self.values

        if pointDownwards:
            sign = np.sign(vec[:,:,2])[:,:,np.newaxis]
            sign[sign==0] = 1
            vec.values *= -sign

        n = np.linalg.norm(vec, axis=2)
        iValid = n > 1E-12
        vec[iValid,:] = vec[iValid] / n[iValid,np.newaxis]

        if angleDef.lower() == 'polar+azimut':
            phi   = np.arctan2(vec[:,:,1], vec[:,:,0])
            theta = np.arccos(vec[:,:,2])

        elif angleDef.lower() == 'strike+dip':
            # check if this is what we/our MEs excpect to be strick+dip
            phi,theta = vector2plunge_bearing(vec[:,:,0], vec[:,:,1],
                                              vec[:,:,2], asRad=True)
            theta = theta + np.pi/2.

        elif angleDef.lower() == 'bearing+plunge':
            phi,theta = vector2plunge_bearing(vec[:,:,0], vec[:,:,1],
                                              vec[:,:,2], asRad=True)

        # preserve Dim if nPoints and/or nFrames == 1
        phi = phi.reshape(self.values.shape[:-1])
        theta = theta.reshape(self.values.shape[:-1])

        names = lambda s: '%s(%s)' % (s, self.name)
        names = map(names, angleDef.split('+'))

        A = ScalarField(names[0], phi,
                        points = self.topo, frames = self.frames)

        B = ScalarField(names[1], theta,
                        points = self.topo, frames = self.frames)

        return A, B
#} #End Common functions for VectorFields


class TensorField(ScalarField):
    """
    Structure to store time (frame) dependend tensor fields on a topography.

    The timedependend (frames) tensor field will be stored as a numpy array with
    the shape of (nPoints, nFrames, 3, 3) in TensorField.values.
    """
    def __init__(self, *args, **kwargs):
        ScalarField.__init__(self, *args, **kwargs)

        if kwargs.get('enablePlot', True):
            self.plot = None #none defined yet

    def _checkValuesConsistency(self):
        '''Checks the consistency of values attribute.'''
        dims = self.values.shape
        if not len(dims) == 4 or not dims[2:] == (3,3):
            raise ValueError("The TensorField.values has to be a 4D-numpy-"
                             "array of shape Npoints x Nframes x 3 x 3")


#{Common functions for TensorFields
    def toScalars(self, nameFun=None):
        """Returns TensorField-components as ScalarFields

        @param nameFun: namingFunction(self.name,i,j) -
            default: nameFun = lambda s: '%s_%d%d' % (s[0], s[1]+1, s[2]+1)

        @returns: dictionary with ScalarFields of components
        """

        if nameFun is None:
            nameFun = lambda s: '%s_%d%d' % (s[0], s[1]+1, s[2]+1)

        tensorIsSym = self.isSymmetric()
        tensorIsDiag = True

        rtnFields = {}
        #loop rows
        for i in range(3):
            #loop columns
            for j in range(3):
                if tensorIsSym and (i,j) in [(0,1), (1,2), (0,2)]:
                    # skip if tensor is symmetric
                    continue

                cName = nameFun((self.name, i, j))
                values = self.values[:,:,i,j]

                # detect if tensor is diagonal
                if tensorIsSym and not i == j and np.any(values):
                    tensorIsDiag = False

                rtnFields[cName] = ScalarField(cName, values,
                                               topo=self.topo,
                                               frames=self.frames)
        if tensorIsSym and tensorIsDiag:
            # remove nondiagonal fields
            for s in [(0,1), (1,2), (0,2)]:
                rtnFields.pop(nameFun((self.name, s[0], s[1])))

        return rtnFields

    def eigen(self, directionType='n', pointDownwards=False):
        """ Returns the ordered EigenValues and their directions

        @param directionType:
            - 'n' (default) returns the EigenVectors
            - 'bearing+plunge'
            - 'strike+dip' (needs a check)
            - 'polar+azimut'

        @param pointDownwards: (default = False), forces the eigenDirection to
            point downwards (southernHemisphere)

        @returns: sorted 3-Element-list of (eigenval, *direction)-tuples
        """

        from numpy_geom import eigen as getEigen

        t = self.values
        # get sorted eigenValues and related eigenVectors
        eigs, eigsV = getEigen(t)

        rtn = []
        # loop the three eigenComponents
        for ii, cc in zip(range(3), ['I', 'II', 'III']):
            # ScalarField to save eigenValue
            s = ScalarField('%s_%s' % (self.name, cc), eigs[:,:,ii],
                            topo=self.topo, frames=self.frames)

            # VectorField to save/calculate eigenVector
            n = VectorField('%sDir_%s' % (self.name, cc), eigsV[:,:,ii,:],
                            topo=self.topo, frames=self.frames)

            if directionType == 'n':
                # return 'pure' eigenVector
                nn = n.values
                if pointDownwards:
                    sign = np.sign(nn[:,:,2])[:,:,np.newaxis]
                    sign[sign==0] = 1
                    n.values *= -sign
                rtn.append((s, n))

            else:
                # return (direction) angles of eigenVector
                rtn.append((s,) + n.vectorToAngles(
                    angleDef=directionType,pointDownwards=pointDownwards))

        return rtn


    def isSymmetric(self):
        """Checks if TensorField is symmetric"""
        return np.any( self.values - np.transpose(self.values,(0,1,3,2) ) )
#} #End Common functions for TensorFields


def _test_indexing():
    sField = [[ii + 10**jj for jj in range(3)] for ii in range(4)]
    pts = PointCloud(np.arange(4*3).reshape(4,3))
    frs = Frames(ids=[(2,ii) for ii in range(3)])
    s = ScalarField('scalarField', sField, topo=pts, frames=frs)

    print (s[1,-1].values == [[101]]).all()
    print (s[[0,1],[0,1]].values == [[1,10],[2,11]]).all()
    print (s[1:].values == [[2,11,101],[3,12,102],[4,13,103]]).all()
    print (s[:,[True, False, True]].values ==[[1, 100],[2,101],[3,102],[4,103]]
            ).all()
    print s[[2],[]] is None

    vField =  [[[ii + 10**jj + .1*kk
                 for jj in range(3)]
                 for kk in range(3)]
                 for ii in range(4)]

    v = VectorField('vectorField', vField, topo=pts, frames=frs)
    print (v[:,:,0].values == [[1. , 1.1, 1.2],
                                [2. , 2.1, 2.2],
                                [3. , 3.1, 3.2],
                                [4. , 4.1, 4.2]]).all()

    print (v[1,:,[1,2]].values == [[[ 11. , 101. ],
                                    [ 11.1, 101.1],
                                    [ 11.2, 101.2]]]).all()
    # test naming
    print v[:,:,[0,2]].name == 'vectorField[0,2]'
    print v[:,:,[0,2]][:,:,0].name == 'vectorField_1'

    tField = np.random.randn(4,3,3,3)
    t = TensorField('tensorField', tField, topo=pts, frames=frs)
    print t[:,:,:].name == 'tensorField'
    print t[:,:,1].name == 'tensorField[1][0,1,2]'
    print t[:,:,1,[True, False,False]].name == 'tensorField_21'
    print t[:,:,[True, False,False]][:,:,1].name == 'tensorField_12'
    print t[:,:,1][:,:,[0,2]].name == 'tensorField[1][0,2]'
    print t[:,:,1][:,:,[0,2]][:,:,0].name == 'tensorField_21'
    print t[:,:,-1,-1].name == 'tensorField_33'


def _test_fieldsPointCloud():
    #
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from bae.plot_01.mpl_core import southernHemisphere,HarrissonAxes,MapAxes
    from cStringIO import StringIO
    from fieldscollection import FieldsCollection as Fields

    ### generate some RandomData / nothing special...
    # Make the grid
    x, y, z = np.meshgrid( np.arange(-0.8, 1, 0.11),
                           np.arange(-0.8, 1, 0.12),
                           np.arange(-0.8, 1, 0.83))

    x, y, z = x.ravel(), y.ravel(), z.ravel()



    # Make the direction data for the arrows
    u = np.sin(np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z)
    v = -np.cos(np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z)
    w = (np.sqrt(2.0 / 3.0) * np.cos(np.pi * x) * np.cos(np.pi * y) *
         np.sin(np.pi * z))

    coords = np.vstack((x, y, z)).T

    nFrames = 8
    frameIds = [(2,i) for i in range(nFrames)]

    val = np.zeros((u.size, nFrames, 3))

    np.random.seed(seed=123456)

    for ii, random in zip(range(nFrames), np.random.randn(nFrames,3)):
        val[:,ii,0] = (10 + 5*random[0]) * u
        val[:,ii,1] = (1 + .3*random[1]) * v
        val[:,ii,2] = (1 + 1*random[2]) * w

    ### start reading here
    points = PointCloud()
    points.tol = 1E-10
    points.addVariable('coords', coords)
    frames = Frames(ids=frameIds)

    #crate a Fields-container and store a (vector-)field named 'U'
    data = Fields(frames=frames, topo=points, U=val)

    # test writing to csv
    tmpOut = StringIO()
    data.writePvfCsv(tmpOut)

    # test reading from csv
    tmpOut.seek(0)
    data2 = Fields.fromPvfCsv(tmpOut)
    data2.reshape.toVector(['U_1', 'U_2', 'U_3'], vectorName='U', inplace=False)

    # compar I/O
    msg(('IO same Vals: ',
         np.all(np.isclose(data['U'].values, data2['U'].values))))
    msg(('IO same Points: ',
         np.all(np.isclose(data.topo.coords, data2.topo.coords))))


    data.frames.regularSequence(unit = 'monthly', start = (2019, 'Feb'),
                                 unitTimes=[3*i for i in range(nFrames)])

    msg('frame labels:', data.frames.seq.labels )

    if True:
        fig = plt.figure()

        ## plunge and bearing
        axh = fig.add_subplot(3, 2, 1, projection='shemisphere')
        ## also test indexing by FrameId
        data['U'][:,1].plot.plungeBearingDensity(axh = axh)
        data['U'][:,1].plot.plungeBearingScatter(axh = axh)

        ## 3d quiver-plot
        axh = fig.add_subplot(3, 2, 2, projection='3d')
        data['U'][:,:].plot.quiver3D(axh = axh)

        ## plot over frames
        axh = fig.add_subplot(3, 2, 3)
        # also test indexing by integer
        data['U'][5,:,0].plot.plotOverFrames(currentIndex = 4, axh = axh, label='one')
        data['U'][5,:,1].plot.plotOverFrames(currentIndex = 4, axh = axh, label='two')
        data['U'][5,:,2].plot.plotOverFrames(currentIndex = 4, axh = axh, label='three',
           seqTicks = True)
        axh.legend()

        ## plot hist
        axh = fig.add_subplot(3, 2, 4)
        # also test indexing by indexlist
        ptIdx = np.arange(data.topo.coords.shape[0]-5) #.tolist() #list works also
        data['U'][ptIdx,:,0].plot.hist(axh = axh)

        ## plot harrisson(background)
        axh = fig.add_subplot(3, 2, 5, projection='harrisson')
        if True:
           axh._harrissonLabels[-1]= 'BOOM!'
           axh._harrissonBackgroundCmap= plt.get_cmap('Reds')
           axh.updateHarrisson()
        axh.plot(np.arange(7), .25*np.arange(7)*(1 + .2*np.random.rand(7)))

        # plot map(background)
        axh = fig.add_subplot(3, 2, 6, projection='mapaxes')
        axh.setMap('../plot_01/mietkuendigung.png', [-150, 100, -100, 100])

        plt.show()

def _test_fieldsPointStructuredGrid():
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from bae.plot_01.mpl_core import southernHemisphere, HarrissonAxes, MapAxes
    from cStringIO import StringIO

    ### generate some RandomData / nothing special...
    # Make the grid
    xx, yy, zz = np.meshgrid(np.arange(-0.8, 1, 0.11),
                          np.arange(-0.8, 1, 0.12),
                          np.arange(-0.8, 1, 0.83))

    x, y, z = xx.ravel(), yy.ravel(), zz.ravel()


    # Make the direction data for the arrows
    u = np.sin(np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z)
    v = -np.cos(np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z)
    w = (np.sqrt(2.0 / 3.0) * np.cos(np.pi * x) * np.cos(np.pi * y) *
         np.sin(np.pi * z))

    coords = np.vstack((u,v,w)).T

    nFrames = 8
    frameIds = [(2,i) for i in range(nFrames)]

    val = np.zeros((u.size, nFrames, 3))

    np.random.seed(seed=123456)

    for ii, random in zip(range(nFrames), np.random.randn(nFrames,3)):
        val[:,ii,0] = (10 + 5*random[0]) * u
        val[:,ii,1] = (1 + .3*random[1]) * v
        val[:,ii,2] = (1 + 1*random[2]) * w

    ### start reading here
    grid = StructuredGrid(
                gridPtNb=xx.shape,
                origin=(x[0],y[0],z[0]),
                spacing=(.11, .12, .83)
                )
    frames = Frames(frameIds)


if __name__ == "__main__":
    print 'No syntax errors.'
    # _test_fieldsPointCloud()
    # _test_fieldsPointStructuredGrid()
    # _test_indexing()


