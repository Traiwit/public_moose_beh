# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 10:07:15 2018

@author: Tobias
"""
import sys

import numpy as np
from bae.log_01 import msg

#check module is loaded by epydoc-build
#--> epyDocBuild=True to disable decorators
if 'epydoc' in sys.modules:
    epyDocBuild = True
else:
    epyDocBuild = False


class Points(object):
    def __init__(self, *args, **kwargs):

        if len(args) >= 1 and not args[0] is None:
            self.xyz = np.asarray( args[0] )
            self._checkDim(self.xyz)

        else:
            self.xyz = np.zeros( (0,3) )

        if len(args) >= 2 and not args[1] is None:
            self.ids = np.asarray( args[1] )
        else:
            self.ids = np.array([], dtype='S64')

        self._parent = kwargs.pop('parent', None)


    @staticmethod
    def _checkDim(xyz):
        shape = xyz.shape
        if len(shape) == 2 and shape[1] == 3:
            return
        else:
            raise ValueError('PointVector needs to be of dimension Npts x 3' +
                             'Yours is %s.' % str(shape))

    @staticmethod
    def _KDTree(xyz):
        from scipy.spatial import cKDTree

        Points._checkDim(xyz)
        if xyz.shape[0] > 1E6:
            msg('Building KDTree for %d 3D-Points' % xyz.shape[0])

        tree = cKDTree(xyz)

        return tree


    @staticmethod
    def _uniqueTreeIdx(tree, atol = 1E-6, normType = 'sphere'):

        if normType == 'sphere':
            pNorm = 2.
        elif normType == 'box':
            pNorm = np.inf
        else:
            raise ValueError("Don't know a normType %s" % str(normType))

        groups = tree.query_ball_tree(tree, atol, p = pNorm)
        #  group = [[self1, 1st neighbour, 2nd, ...],       #first pt in list
        #           [self2, 1st neighbour, 2nd, ...],...    #second pt
        #           [selflast, 1st neighbour, 2nd, ...]]
        # --> group holds redundant for each neigbour --> use np.unique

        #sort each group
        groups = np.sort(groups, axis=0)

#        groups, rIdx  = np.unique(groups, axis=0, return_index=True)
#        # np.unique returns a sorted list --> use return_index to undo this
#        # sorting
#        groups = groups[rIdx]
        try:
            # 'balanced'-Points, len(group) for group in groups is equal
            # eg. each field has same pointdef
            groups  = np.unique(groups, axis=0)
        except TypeError:
            # 'unbalanced'-Points, len(group) for group in groups differ
            # eg. each fields are definend on different point(sub)sets
            groups = np.unique(groups)

        return groups


    @staticmethod
    def _ensureUnique(arg, atol = 1E-6):
        import scipy
        if isinstance(arg, np.ndarray):
            Points._checkDim(arg)
            tree = Points._KDTree(arg)
        elif isinstance(arg, scipy.spatial.ckdtree.cKDTree):
            tree = arg
        else:
            raise ValueError("Input must be numpy-array or cKDTree (scipy)")

        idx, _ = Points._uniqueTreeIdx(tree, atol = atol)

        if len(idx) == tree.n:
            return True
        else:
            return False


    def subSetByPoints(self, points, atol = 1E-6, normType = 'sphere',
                       strict = True):
        import itertools as IT

        if normType == 'sphere':
            pNorm = 2.
        elif normType == 'box':
            pNorm = np.inf
        else:
            raise ValueError("Don't know a normType %s" % str(normType))

        rstTree = self._KDTree(points)
        xyzTree = self._KDTree(self.xyz)

        groups = rstTree.query_ball_tree(xyzTree, atol, p = pNorm)

        if strict:
            ident = np.array( map(len, groups) )
            if not np.all(ident == 1):
                raise ValueError("Never found %d points and found %d points" %
                                 ((ident == 0), (ident > 1)) +
                                 'multible times')

        indices = np.unique( IT.chain.from_iterable(groups) )

        return indices, groups



    def __getitem__(self, item):
        from copy import deepcopy
        rtn = deepcopy(self)

        rtn.xyz = rtn.xyz[item]
        if rtn.ids.size:
            rtn.ids = rtn.ids[item]

        return rtn


    def subSetByIds(self, pattern, strict=True):
        if self.ids is None:
            raise ValueError("No point ids defined yet.")

        if not len(self.ids) == self.xyz.shape[0]:
            raise IndexError("Id- or point-definition is corrupt.\n" +
                             "Length id: %d\n Length points: %d" %
                             (len(self.ids), self.xyz.shape[0]) )


        if isinstance(pattern, str):
            pattern = [pattern,]
        pattern = np.asarray(pattern)


        checkFun = lambda p: lambda r: all([q in r.split('_')
                                            for q in p.split('_')])
        inMask = np.array([map(checkFun(p), self.ids) for p in pattern]).T

        if strict:
            patternFound = np.any(inMask, axis = 1)
            if not np.all(patternFound):
                raise ValueError("The following patterns couldn't be found" +
                                 "in point ids: \n" +
                                 '\n'.join(pattern[~patternFound]))

        bIndices = inMask.any(axis=1)

        return bIndices


class Frames(object):
    def __init__(self, *args):

        if len(args) >= 1 and not args[0] is None:
            self.ids = args[0]
        else:
            self.ids = []

        if len(args) >= 2 and not args[1] is None:
            self.names = np.asarray( args[1] )
        else:
            self.names = np.array([], dtype='S64')

    def idsToIdx(self, ids):
        if isinstance(ids, tuple) and isinstance(ids[0], int):
            ids = [ids,]

        try:
            index = [self.ids.index(id) for id in ids]
        except ValueError:
            raise ValueError("Can't find frame with id %s." % str(id))

        return index

    def namesToIdx(self, names):
        if isinstance(names, str):
            names = [names,]

        try:
            index = [self.names.index(name) for name in names]
        except ValueError:
            raise ValueError("Can't find frame with id %s." % str(id))

        return index


    def _defaultNamesFromIds(self):
        if not self.ids:
            return

        self.names = np.array(['F%d03%d' % id for id in self.ids])


    def _defaultIdsFromNames(self, defaultStep=2, frameDigits = 3):
        import re
        def parser(s):
            nums = re.findall(r'\d+', s)
            if len(nums) == 1:
                if int(nums[0]) > 10**frameDigits - 1:
                    # Type of F2014
                    return (int(nums[0][0]), int(nums[0][1:]))
                else:
                    # Type of F014
                    return (defaultStep, int(nums))
            elif len(nums) == 2:
                #Type of Step2Frame014
                return (int(nums[0]), int(nums[1]))
            else:
                raise ValueError("Can't parse a valid frameId-tuple from "
                                 "'%s'" % s)

        frameIds = map(parser, self.names)
        self.ids = frameIds

    def __getitem__(self, item):
        from copy import deepcopy
        rtn = deepcopy(self)

        rtn.ids = map(tuple, np.array(self.ids)[item])
        if np.array(rtn.names).size:
            rtn.names = list( np.array(self.names[item]) )
        return rtn


class FieldValues(dict):
    def __init__(self, parent = None):
        pass


    def _checkConsistency(self, fieldNames = None):
        if len(self) == 0:
            return

        if fieldNames is None:
            fieldNames = self.keys()
        elif isinstance(fieldNames, str):
            fieldNames = [fieldNames, ]

        notFoundFields = set(fieldNames) - set(self.keys())
        if notFoundFields:
            raise KeyError("The following fields can't be found:\n" +
                           ", ".join(notFoundFields))

        shapes = [f.shape[0:2] for n, f in self.iteritems()
                    if  n in fieldNames]

        if not all([s == shapes[0] for s in shapes]):
            raise ValueError("Found inconsitent frame and/or point "
                             "definition for field value data.")



    def toVector(self, newName, fieldNames, inplace = False):
        self._checkConsitency(fieldNames)
        newDim = self[fieldNames[0]].shape[:2] + (len(fieldNames),)
        newValues = np.empty(newDim)
        try:
            for kk, field in enumerate(fieldNames):
                newValues[:,:,kk] = field
        except ValueError as e:
            raise ValueError("The following Error occured. Did you try to " +
                             "stack a vector- or tensorlike field?\n" +
                             str(e))

        if inplace:
            for field in fieldNames:
                self.pop(field)

        return newValues


    def toTensor(self, fieldNames, inplace, order= 'R') :

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
            e = ("A tensor needs 3 (diag), 6 (symetrical) or 9 scalar " +
                 "fields to be defined. You passed %d." % len(fieldNames) )

        #rearrange varLabels with respect to their position in tensor
        if order == 'D':
            srcIdx = [0,1,2]
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

        self._checkConsistency(fieldNames)

        nPts, nFrames = self[fieldNames[0]].shape[:2]
        newDim = (nPts, nFrames, 3, 3)

        #stacking data to nPointsXnFramesX3X3 numpy array
        arrayTensor = np.zeros(newDim)

        ii = 0
        try:
            for iR in range(3):
                if order == 'D':
                    arrayTensor[:,:,iR,iR] = self[fieldNames[ii]]
                    ii += 1
                else:
                    for iC in range(3):
                        arrayTensor[:,:,iR,iC] = self[fieldNames[ii]]
                        ii += 1

        except ValueError as e:
            raise ValueError("The following Error occured. Did you try to "
                             + "stack a vector- or tensorlike field?\n"
                             + str(e))

        if inplace:
            for field in fieldNames:
                self.pop(field)

        return arrayTensor



def _checkDataConnection(func, passThrough=True):
    '''Checks if data.points.xyz and data.frames.ids is available
    '''
    if epyDocBuild:
        return func
    else:
        def func_wrapper(self, *args, **kwargs):
            try:
                hasPoints = bool(self.points.xyz)
            except:
                hasPoints = False

            try:
                hasFrames = bool(self.frames.ids)
            except:
                hasFrames = False

            if hasPoints and hasFrames:
                return func(self, *args, **kwargs)

            else:
                if not hasPoints:
                    msg('WARNING! data has no valid points assigned!')
                if not hasFrames:
                    msg('WARNING! data has no valid frames assigned!')

                if passThrough:
                    return func(self, *args, **kwargs)
                else:
                   msg('Function %s skipped' % func.__name__)

        return func_wrapper


import matplotlib.pyplot as plt

class CommonScalarPlots(object):
    def __init__(self, data):
        self.data = data

    @staticmethod
    def _initPlotHandles( **kwargs):
        from mpl_toolkits.mplot3d import Axes3D
        # get existing axes handle
        # if existent: this will remove axh from kwargs
        axh = kwargs.pop('axh', None)
        projection = kwargs.pop('projection', None)
        if axh is None:
            # create new fig and axes
            fig = plt.figure()
            axh = fig.add_subplot(1,1,1, projection = projection)
        else:
            # get fig/canvas from axh
            fig = axh.figure
            if projection == '3d':
                if not hasattr(axh, 'get_zlim'):
                    msg('Warning!!! The axes you provided does not support ' +
                        'a %s projection.' % projection)

                    axh = fig.add_axes(axh.get_position(),
                                        projection = projection)




        return fig, axh

    @staticmethod
    def _uniqueMarkers(n):
        markers = [(2+ii, 1+ii%2, ii/n*90.0) for ii in range(1, n+1)]
        return markers


    @_checkDataConnection
    def hist(self, sortByTime = True, **kwargs):
        fig, axh = self._initPlotHandles(**kwargs)
        handles = []
        if sortByTime:
            for frame in xrange(self.data.shape[1]):
                pltH = axh.hist(self.data[:,frame].ravel(), **kwargs)
                handles.append(pltH)
        else:
            pltH = axh.hist(self.data.ravel(), **kwargs)
            handles.append(pltH)

        return pltH

    def scatter3D(self, sortByTime = True, **kwargs):
        fig, axh = self._initPlotHandles(projection='3d', **kwargs)
        handles = []
        if sortByTime:
            x,y,z = (self.data.points.xyz[:,ii] for ii in range(3))
            markers = self._uniqueMarkers(self.data.shape[1])

            for frame in xrange(self.data.shape[1]):

                try:
                    label = 'frame: %s' % str(self.frames.ids[frame])
                except:
                    label = 'frame: %d' % frame

                pltH = axh.scatter(x, y, z,
                                   c = self.data[:,frame].ravel(),
                                   marker = markers[frame],
                                   label = label,
                                   **kwargs)
                handles.append(pltH)

        else:
            xyz = np.tile(self.data.points.xyz, (self.data.shape[1],1))
            x,y,z = (xyz[:,ii] for ii in range(3))
            pltH = axh.scatter(x, y, z,
                   c = self.data.ravel(),
                   **kwargs)

            handles.append(pltH)

        return pltH





class ScalarArray(np.ndarray):
    def __new__(cls, field, fieldName='', points=None, frames=None):
        obj = np.asarray(field).view(cls)

        obj.points    = points
        obj.frames    = frames
        obj.fieldName = fieldName

        obj.plot    = CommonScalarPlots(obj)

        return obj


    def __getitem__(self, item):
        try:
            if isinstance(item, (slice, int)):
                self._new_pt_index = item
            else:
                self._new_pt_index = item[0]
        except:
            pass

        try:
            if isinstance(item, (slice, int)):
                self._new_fr_index = None
            else:
                self._new_fr_index = item[1]
        except:
            pass

        return super(ScalarArray, self).__getitem__(item)


    def __array_finalize__(self, obj):
        self.fieldName = getattr(obj, 'fieldName', '')
        self.frames    = getattr(obj, 'frames', None)
        self.points    = getattr(obj, 'points', None)

        self.plot      = getattr(obj, 'plot', None)

        try:
            self.points = self.points[obj._new_pt_index]
        except:
            pass


        try:
            if not obj._new_fr_index is None:
                self.frames = self.frames[obj._new_fr_index]
        except:
            pass


class VectorArray(ScalarArray):

    def direction(self, angles = 'plungeBearing', pointDownwards = True):
        rtn = self
        if pointDownward:
            vec = self
#            vec = np.sign(vec[])
        pass

    def length(self):
        norm = np.linalg.norm(self, axis=2)
        norm = ScalarArray(norm, fieldName=self.fieldName,
                           points=self.points, frames = self.frames)
        return norm


class FieldContainer(FieldValues):
    def __init__(self, *args, **kwargs):
        super(FieldValues, self).__init__()

        ## point/pointId input
        points   = kwargs.pop('points', None)
        pointIds = kwargs.pop('pointIds', None)


        if isinstance(points, Points):
            self.points = points
        else:
            self.points = Points(points, pointIds)

        self.pointsTol = kwargs.pop('pointTol', 1E-3)

        ## frame/frameId input
        frames   = kwargs.pop('frames', None)
        frameIds = kwargs.pop('frameIds', None)

        if isinstance(frames, Frames):
            self.frames = frames
        else:
            self.frames = Frames(frames, frameIds)

        ## storing fields
        for fieldName, vals in kwargs.iteritems():
            self[fieldName] = ScalarArray(vals, fieldName, points=self.points)

        self._checkConsitency()


    def _checkConsitency(self):
        pass


    def readPvfCsv(self, csvFileName, varNames = [], frameList = [],
                   skipColumnList = [],
                   pointCoordList = [], pointIdList = [], pointIdxList = [],
                   **readerArgs):
        '''Reads data from  PointVariableFrame-ordered csv

        Example:
        ========
        >>> sLabels = ['S11','S22','S33','S23','S13','S12']
        >>> uF = FieldContainer().readFromCsv(csvName,varNames=slabels,
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

        '''
        from bae.misc_01 import npReadCsv, readCsvHeader
        try:
            import pandas as pd
            hasPandas = True
        except ImportError:
            hasPandas = False
            if any([pointCoordList, pointIdList, pointIdxList,]):
                msg("Warning!!! OnTheFly-Filtering by points is only " +
                    "supported for pandas (which you havn't installed). " +
                    "You need need to remove unwanted Points manually.")


        if bool(pointCoordList) + bool(pointIdList) + bool(pointIdxList) > 1:
            raise ValueError('Ambiguous setting of point-filter. ' +
                             'pointCoordList, pointIdList and ' +
                             'pointIdxList can only be used exclusively.')

        ## setUp rules for (row-)conditional reading
        ## if all(rowConditions): append row to data
        rowConditions = []
        if hasPandas:
            if bool(varNames): # @ given varNames
                rowConditions.append(
                        lambda data: np.isin(data['var'], varNames) )

            if bool(pointCoordList): # @ given coordinates
                dist = lambda testPt: np.norm(
                    np.asarray(pointCoordList) - testPt, axis=1)

                rowConditions.append(
                    lambda data: dist(data[['x','y','z']] < self.pointsTol) )

            if bool(pointIdList): # @ given point ids
                rowConditions.append(
                        lambda data: np.isin(data['ptId'], pointIdList) )

            if bool(pointIdxList): # @ given point indexs
                rowConditions.append(
                        lambda data: np.isin(data['ptIdx'], pointIdxList) )

        ## cummulate all rowContitions to one function
        if rowConditions:
            rowCondition = lambda data: all([c(data) for c in rowConditions])
        else:
            rowCondition = None

        ## specify columns to read
        delimiter = readerArgs.pop('delimiter',',')

        header = readCsvHeader(csvFileName, delimiter=delimiter)

        iCols = ['x','y','z','var']

        if 'ptId' in header:
            iCols.append('pdId')

        if 'ptIdx' in header:
            iCols.append('pdIdx')

        if not frameList:
            #read all
            frameList = [h for h in header if h not in iCols]

        readCols = iCols + frameList
        readCols = [c for c in readCols if c not in skipColumnList]

        ## read data
        data = npReadCsv(csvFileName, columns=readCols,
                         rowCondition = rowCondition,
                         delimiter = delimiter, **readerArgs)

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
        xyz = np.vstack((data[q] for q in ['x','y','z'])).T
        nFlatPts = xyz.shape[0]

        points = Points()
        tree = Points._KDTree(xyz)
        groups = points._uniqueTreeIdx(tree, normType='box',
                                            atol=self.pointsTol)

        idxUnique = [ii[0] for ii in groups]

        xyz = xyz[idxUnique,:]
        nPts = xyz.shape[0]

        try:
            pointIds = data['ptId'][idxUnique]
            points = Points(xyz, pointIds)
        except:
            points = Points(xyz)

        self.points = points
        self.frames = frames

        ## create Fields

        ## create a indexLookup to unique Points
        ## use length of xyz to raise IndexError on failed lookup
        lookUp = nPts * np.ones(nFlatPts, dtype=int)
        for ii, look in enumerate(groups):
            lookUp[look] = ii
#        print max(lookUp), nPts
        ## get uniqueFields
        fields = np.unique(data['var'])

        fieldData = np.vstack((data[q] for q in frameList)).T

        ## storing fields
        for field in fields:
            valueIdx = data['var'] == field
            idx = lookUp[valueIdx]
            vals = np.nan * np.ones((nPts, nFrs))
#            print valueIdx, len(valueIdx), fieldData.shape[0], len(idx), valueIdx.sum()
            vals[idx] = fieldData[valueIdx,:]

            self[field] = ScalarArray(vals, field,
                                    xyz = self.points,
                                    frames = self.frames)

        return self




if __name__ == "__main__":
#    nPts = 12
#    nFrs = 4
#    nFld = 8
#
#    data = {'points' : np.random.randn(nPts,3),
#            'frames': [(2,ii) for ii in range(nFrs)],
#            }
#    fields = dict( ('field%d' % d, np.random.randn(nPts,nFrs,1,1))
#                    for d in range(nFld))
#    data.update(fields)
#    pd = FieldContainer(**data)
    csv = 'ESC2017_R26_2017NOV16.csv'
    data = FieldContainer().readPvfCsv(csv)
    u = data['U_3']
