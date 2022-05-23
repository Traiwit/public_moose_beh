"""Class that holds a collection of fields. Possibly only one field.

fancy functions:
fld[] ... fld.__getitem__, __setitem__, __delitem__ (the last is less fancy...)

Some data arrays:
 >>> U = np.array(...).reshape((n,m))
 >>> S = np.array(...).reshape((n,m))
 >>> pointnames = np.array(["pt_%02d" % i for i in range(n)])  # shape = (n,)
 >>> frames = np.array([(2,i) for i in range(m)])  # shape = (m,2)
 >>> frameNames = np.array(["step_%03d" % (i+1) for i in range(m)])  # sh.=(m,)

Now creating the actual fields collection object. You can pass arrays to the
constructor or assign later:
 >>> flds = FieldsCollection(
 >>>     topo=topo,
 >>>     times={"frames": frames, "frameNames": frameNames},
 >>>     U=U(nxm), S=S(nxm))
 >>> flds.space["pointnames"] = pointnames

Extracting one of the fields:
 >>> fldU = flds["U"]  # just the numpy array
 >>> fldU_FC = flds[["U"]]  # fldU_FC contains all times and space arrays

Adding a new persistent field
 >>> fld.space["coords"] = fld.topo.getPoints()

Adding ordinary fields
 >>> fld["T"] = np.array((...))  # (nxm) array == scalar field all times
 >>> fld["T"][:,0:2] = 0.5*fld["T"][:,0:2]  # (nxm) array == scalar, some times
 >>> fld["U"][:,0] = refDispl  # (nx3) ndarray == vector field, one time

Slicing
 >>> fldU[1:10, :]  # ... points 1 to 9, all times
 >>> fldU[1:10, 0:3]  # ... points 1 to 9, some times

To get the plain numpy array from a Fieldscollection that contains only one
field:
 >>> fldU_FC = flds[["U"]]
 >>> npArr = fldU_FC["U"]  # usual access of a field returns the numpy array
 >>> npArr = fldU.array  # more convenient, don't need the field name again
"""
###GP: Consolitate and shorten docs of module and class and __getitem__ method!

# Note: version is stored in __init__.py

import numpy as np
from collections import OrderedDict

# from .topo.structuredgrid import StructuredGrid
from bae.log_01 import msg


class FieldsCollectionCore(object):
    """Structure to store and manage time dependent scalar, vector
    and tensor fields which share the topography and time points.

    There are "persistent" fields (stored in dict-attributes times and space)
    that are automatically transfered to new FieldCollection-objects by the
    L{__getitem__} method --i.e. on component access.

    @ivar fields: A dict {name: np.ndarray} with the ordinary space-time fields
       in this container.
       The first two dimensions must match self.shape or have size one in
       either --not both. I.e. We have N locations and M time points then
       the following shapes are valid: (N,M), (N,M,3,3), (1,M), (N,1,3)
       The following shapes would not be ok: (1,1), (N,3) --assuming M!=3--
    @ivar times: A dict of ndarrays with persistent space-less fields like
       time values, frame names and the like. Must have shape (M,...).
       They are automatically transfered to results of any slicing operations.
    @ivar space: A dict of ndarrays with persistent time-less fields like
       point coords, point names and the like. Must have shape (N,...)
       They are automatically transfered to results of any slicing operations.
    @ivar topo: The topology for the fields like L{PointCloud},
       L{StructuredGrid}, Mesh??
    """

    def __init__(self, topo=None, times=None, **kwargs):
        """Constructor.
        All arguments have to be passed as keyword arguments. Positions of the
        arguments are subject to change.

        Usage:
        ======
         >>> f = FieldsCollections(
         ...         topo=myTopo, times=myFrames,
         ...         aField=aScalarFieldObject,
         ...         bField=aNumpyArrayThatFits,
         ...         )

        @param topo: a (optional) L{PointCloud}, L{StructuredGrid} or other
           topo-object.
        @param times: can be a dict {field name: ndarray} for initialization of
           self.times. Note that the space dimension of those fields must not
           be there. I.e. if we have N points in space and M time points then
           the ndarrays for the times parameter must have shape (M,), (M,2),
           or so. Not (1,N,...)!

           Or an iterable that can be converted into an ndarray with first
           dimension==N. It will be stored to self.times with default name
           "times".
        @param kwargs: fields for initializing the fields collection.
        """
        # initialize main attributes
        self.topo = topo
        self.space = dict()
        self.times = dict()
        self.fields = OrderedDict()


        # times is a dict
        if isinstance(times, dict):
            self.times.update(times)
        elif times is not None:
            self.times["times"] = np.asarray(times)


        # storing fields
        for fieldName, field in kwargs.iteritems():
            # for something like FieldsCollection(T=np.array(...))
            # or other iterables like lists and the like...
            field = np.asarray(field)

            # special case: 1D array as time-less / space- field
            if len(field.shape)==1:
                shape = self.getShape()
                if shape[0]<=1 or field.shape[0]==shape[0]:
                    field = field[:, np.newaxis]

            # store the data
            if len(field.shape) < 2:
                raise ValueError(
                    "The field '%s' passed to the FieldsCollection-init doesn't"
                    " have both space and time dimensions (but shape: %s)"
                    % (fieldName, field.shape))
            self.fields[fieldName] = field


        # finalize: check consistency
        self.checkConsistency()

    @property
    def array(self):
        """In case this FieldsCollection object contains only one field
        then this is a shortcut to it.

        @raises IndexError: when requested for a FieldsCollection with more
            than one fields
        """
        if len(self.fields) != 1:
            raise IndexError("FieldsCollection.array attribute not valid for"
                             " FieldsCollections of length other than one.")
        return self.fields.itervalues().next()

    def __len__(self):
        """len(fieldscollection) returns the number of fields in the
        fieldscollection-object.
        """
        return len(self.fields)

    def getShape(self):
        """Return the tuple of (number of locations (spacial field size),
        number of time-points).

        This tuple will generally comply with the shape of a single scalar
        field in self.

        How it is determined:
         - space dimension from 1. topo, 2. space, 3. fields  (check all
           unless != 1)
         - time dimension from 1. times, 2. fields (check all unless != 1)
        """
        # block to determine spacial size N (nb. of space points)
        N = 0
        while 1:
            if self.topo:
                N = len(self.topo)
                break
            for fld in self.space.itervalues():
                N = fld.shape[0]
                if N>1:
                    break
            for fld in self.fields.itervalues():
                N = fld.shape[0]
                if N>1:
                    break
            break
        # block to determine number of time points
        M = 0
        while 1:
            for fld in self.times.itervalues():
                M = fld.shape[0]
                if M>1:
                    break
            for fld in self.fields.itervalues():
                M = fld.shape[1]
                if M>1:
                    break
            break
        return (N, M)

    def checkConsistency(self, fieldNames=None):
        """
        Checks if all (fieldNames=None) or specified fields are
        consistent with stored topo, times and space objects

        Checking sizes of all dimensions.

        @param fieldNames: list of field names. If given, only check those
           fields for efficiency. (Not implemented yet, currently this
           parameter is being ignored.)
        """

        # Check space dimensions with topo.
        # This check is somehow redundant because getShape would get the nb of
        # space dimensions from topo in the first place
        shape = self.getShape()
        if self.topo is not None and (
                shape[0] != len(self.topo)  # nb of space dims must match...
                # ... or equal one (if we only have space constant fields)!
                and shape[0] != 1):
            raise ValueError(
                "Mismatching shape: The number of values (N=%d) does not match"
                " the number of points in the topo (%d)."
                % (shape[0], len(self.topo)))

        # check space dimensions of self.fields and self.space
        invalidDims = [
            (n, a.shape[0])
            for n, a in self.fields.iteritems()
            if a.shape[0] != shape[0] and a.shape[0] != 1]
        invalidDims.extend(
            ("%s (persistent time-less)" % n, a.shape[0])
            for n, a in self.space.iteritems()
            if a.shape[0] != shape[0])
        if invalidDims:
            raise ValueError(
                "There is a mismatch in space dimensions of stored fields"
                " (size should be %d):\n%s"
                % (shape[0],
                   "\n".join(('size %s: %d' % nameLen)
                             for nameLen in invalidDims)))

        # check time dimensions
        invalidDims = [
            (n, a.shape[1])
            for n, a in self.fields.iteritems()
            if a.shape[1] != shape[1] and a.shape[1] != 1]
        invalidDims.extend(
            ("%s (persistent space-less)" % n, a.shape[0])
            for n, a in self.times.iteritems()
            if a.shape[0] != shape[1])
        if invalidDims:
            raise ValueError(
                "There is a mismatch in time dimensions of stored fields"
                " (size should be %d):\n%s"
                % (shape[1],
                   "\n".join(('size %s: %d' % nameLen)
                             for nameLen in invalidDims)))

    def __iter__(self):
        return self.fields.__iter__()

#{Set and get fields
    def __setitem__(self, name, field):
        """Mimics the dictionary __setitem__ function, so you can do
        the following:
         >>> f = Fields(topo=myTopo, times=myFrames, someFieldsInADict)
         >>> # set a new field from numpy array
         >>> f['viaNumpy'] = np.zeros((nPt,nFr))
         >>> # set a new field from existing field
         >>> f['viaField'] = anExistingField
        """

        if isinstance(name, basestring):

            shape = self.getShape()
            # field value is a numpy array
            if isinstance(field, np.ndarray):
                # special case: 1D array as time-less / space- field
                if len(field.shape)==1 and field.shape[0]==shape[0]:
                    field = field[:, np.newaxis]
                if (not (
                        field.shape[:2]==shape
                        # or field is time-less / space- field
                        or (field.shape[0]==shape[0] and field.shape[1]==1)
                        # or field is a space-less / times- field
                        or (field.shape[1]==shape[1] and field.shape[0]==1)
                        # adding a time series to a (so far) time-less collection
                        or (field.shape[0]==shape[0] and shape[1]<=1)
                        # adding space values to a (so far) space-less collection
                        or (field.shape[1]==shape[1] and shape[0]<=1))):
                    raise ValueError(
                        "Shape of the field to add %s is %s and does not"
                        " match that of the fields collection %s."
                        % (name, field.shape, shape))

                # add the new field
                self.fields[name] = field

            else:
                raise ValueError(
                    "The field (%s) to be added to fields collection must be a"
                    " numpy ndarray. Got %s." % (name, type(field)))

        else:
            raise ValueError(
                "FieldsCollection.__setitem__ currently only accepts single"
                " strings as keys. Got %s." % type(name))

        # consistency check is redundant here. Checked new field before inserting it
        # self.checkConsistency(name)

    def update(self, *args, **kwargs):
        """Mimics the dictionary update function.

        Add np.ndarrays:
         >>> f = FieldsCollection(topo=myTopo, times=myTimes)
         >>> # add numpy arrays as new fields
         >>> f.update({
         ...     # numpy arrays or fields with matching nb of times and points
         ...     'npScalar' : np.zeros((nPt,nTimes)),
         ...     'npVector' : np.zeros((nPt,nTimes,3)),
         ...     'timeLess' : np.zeros((nPt,1)),
         ...     'spaceLess': np.zeros((1,nTimes)),
         ...     })

        Add all fields of an existing FieldsCollection:
         >>> flds1 = FieldsCollection(U=..., V=...)
         >>> flds2 = FieldsCollection(PST=..., matCode=...)
         >>> flds2.update(flds1)
         >>> print sorted(flds2.fields)
         ... ["matCode", "PST", "U", "V"]
        """
        if len(args) > 1:
            raise TypeError("update expected at most 1 arguments, "
                            "got %d" % len(args))

        # process another FieldsCollection as single argument
        if len(args)==1 and isinstance(args[0], FieldsCollectionCore):
            for key, value in args[0].fields.iteritems():
                self[key] = value

        # process other iterables with key-value pairs
        else:
            other = dict(*args, **kwargs)
            for key, value in other.iteritems():
                self[key] = value

    def __delitem__(self, item):
        """Remove fields by their name(s)

        Usage:
         >>> del flds["U"]  # remove just one field
         >>> del flds[["V", "W"]]  # remove many fields at once
        """

        # item is a string or list of strings
        isStr = isinstance(item, basestring)
        if isStr or isinstance(item, list):
            if isStr:
                item = [item,]
            elif not all(isinstance(it, basestring) for it in item):
                raise IndexError(
                    "FieldsCollection.__delitem__: List of index-values must"
                    " contain only strings but got %s." % item)

            for it in item:
                del self.fields[it]

        else:
            raise IndexError(
                "FieldsCollection.__delitem__: Index-value must be a string or"
                " a list of strings. Instead got %s." % item)

    @staticmethod
    def _sliceTakesAll(index, length):
        """If the specified <index>-slice would select all of the <length>
        items then return None to indicate that subslicing with that index is a
        no-op and is not necessary. Otherwise return <index> untouched.
        """
        if ((index.start is None or index.start==0)
            and (index.stop is None or index.stop>=length)
            and (index.step is None or index.step==1)):
            return None
        else:
            return index

    def __getitem__(self, indexes):
        """Get the FieldsCollection with only the subset according to the
        indexes given as parameter.

        There eare generally two cases:

        If a single string is used as index then the corresponding numpy array
        is retrieved and returned. This allows in-place data modification, i.e.
        the returned object is not a copy.

        In any other case a FieldsCollection object is returned. A list of
        strings and/or numeric or boolean indexes can be passed simultaneously.

        It's not possible to combine a single string index with numeric or other
        indexes.

        The topo attribute as well as the persistent fields in self.times and
        self.space are automatically being indexed simultaneously in the first
        (space) and second (time) dimension. Should there be further dimensions
        in fields from self.space or self.times those will not be considered.
        I.e. vector or tensor fields in self.times and self.space keep all
        their components even if the indexes argument contains component
        indexes.

        Usage to get a numpy array:
        ===========================
         >>> array =  myFields['U']  # ndarray all pts, all times for just 'U'
         >>> subArr =  myFields['U'][:, 10:15]  # 'U'-ndarray all pts, 5 times

        Usage to get a new FieldsCollection object:
        ===========================================
         >>> single = myFields[['U',]]  # all pts, all times for just 'U'
         >>> multiple = myFields[['U', 'S']]  # two fields: 'U' and 'S'
         >>> someTimes = myFields[:, 10:15]  # all positions, 5 times
         >>> somePts = myFields[10:15]  # 5 positions
         >>> subPtSome = myFields[['U', 'S'], 10:15]  # 2 flds, 5 pts
         >>> Uxy = myFields[:, :, 1:3]  # all pts, all times, 2 components
         >>> subAbr = myFields[['U', 'S'], 10:15, 1:3]

        B{Allowed index types} are:
         - a string myFields["U"]  ... get the single numpy array
         - list of strings: myFields[["U","S"]]  ... select fields
           -> new FieldsCollection
         - integer values: myFields[1,2]  ... select (space-) point 1, time 2
         - slices: myFields[1:,:-4]   ... all but first pts, times to 4th last
         - flat index lists or np.arrays: myFields[[1,3,5], [0,1]]
           ... points 1,3,5 and times 0,1
         - flat boolean arrays with same length of indexed dimension:
           myFields[[True,False,....]]  ... useful for conditional expressions
           like myFields[myFields["status"]

        Please note that multidimensional indexing is slightly different from
        numpy-syntax. In numpy indexing like myFields[[1,4],[3,5]] will return
        the two values myField13 and myField45. Here it would return
        a 2 points x 2 times data field (similar to npArr[[[1],[4]], [[3,5]]]).

        B{Indexing of components}:
        The components of vector- and tensor fields can be indexed as well.
        Contrary to the space- and times- indexes, in this case the
        corresponding dimension is collapsed in the result. I.e. myFields[1,:]
        will return a FieldsCollection object for all times and only the second
        point. The space/point dimension has only length one but still exists.
        On the other hand: myFields["U",:,:,0] would return a FieldsCollection
        object with all points and all times and only the first component of
        the vector "U". After U1 = myFields["U",:,:,0] requesting U1[:,:,0]
        would then raise an index error because U1 contains only a (is a)
        scalar field.

        @Note: Field names have to be lists not tuples. That's how they are
        differentiated from a multidimensional index.

        @Note: Field names should be valid Python names. Others might work but
        are not guaranteed to do!

        @raises KeyError: If the field name (first and string-valued index)
        does not exist in self.fields.

        @raises IndexError: If numeric indexes (location, time, tensor) are
        wrong.
        """

        # special case #1: indexes argument is single string
        # -> return numpy array
        if isinstance(indexes, basestring):
            return self.fields[indexes]

        # special case #2: indexes argument contains a single string but there
        # more... --> Error
        if isinstance(indexes, tuple) and isinstance(indexes[0], basestring):
            raise ValueError(
                "Selecting a single field results in a numpy array. Subslicing"
                " space, time or tensor-dimesion at the same time is"
                " prohibited.")

        # convert into tuple if indexes is a singleton
        # (be it a list of fieldNames, space index(es))
        if not isinstance(indexes, tuple):
            indexes = (indexes,)

        # ... 1. field names (optional, only if it contains a string)
        if isinstance(indexes[0], basestring):
            raise ValueError(
                "We should never get here. This case is treated earlier."
                " Program error in FieldsCollection.__getitem__() with"
                "indexes = %s" % repr(indexes))
        elif (hasattr(indexes[0], '__iter__') and
              isinstance(indexes[0][0], basestring)):
            # ... list of field names
            fieldNames = indexes[0]
            if (fieldNames and
                not all([isinstance(n,basestring) for n in fieldNames])):
                raise IndexError("FieldsCollection.__getitem__: All requested"
                                 " fieldNames have to be strings")
            indexes = indexes[1:]
        else:
            fieldNames = list(self.fields)

        if not indexes:
            # shortcut if only fieldNames requested
            fieldArgs = dict((n, self.fields[n]) for n in fieldNames)
            result = type(self)(
                 topo=self.topo, times=self.times, **fieldArgs)
            result.space.update(self.space)
            return result

        shape = self.getShape()

        # ... 2. space index
        try:
            spaceIndex = indexes[0]
        except IndexError:
            spaceIndex = None
        else:

            if isinstance(spaceIndex, (int, long)):
                # convert single index into index list to prevent dim collapse
                spaceIndex = [spaceIndex,]

            elif isinstance(spaceIndex, slice):
                # convert a slice that selects all to None
                spaceIndex = self._sliceTakesAll(spaceIndex, shape[0])

            else:
                spaceIndex = np.squeeze(spaceIndex)
                if len(spaceIndex.shape)>1:
                    raise IndexError(
                        "FieldsCollection.__getitem__ got a space index of"
                        " shape %s with more than one dimension."
                        " If you use boolean indexing make sure to apply the"
                        " condition to a time-constant scalar field."
                        % spaceIndex.shape)

        # ... 3. times index
        try:
            timesIndex = indexes[1]
        except IndexError:
            timesIndex = None
        else:
            if isinstance(timesIndex, (int, long)):
                # convert single index into index list to prevent dim collapse
                timesIndex = [timesIndex,]

            elif isinstance(timesIndex, slice):
                # convert slice that selects all to None
                timesIndex = self._sliceTakesAll(timesIndex, shape[1])

            else:
                timesIndex = np.squeeze(timesIndex)
                if len(timesIndex.shape)>1:
                    raise IndexError(
                        "FieldsCollection.__getitem__ got a times index of"
                        " shape %s with more than one dimension."
                        " If you use boolean indexing make sure to apply the"
                        " condition to a space-constant scalar field."
                        % timesIndex.shape)

        # ... 4. tensor index
        tensorIndex = indexes[2:]

        # ... prepare topo, space and times attributes for new FieldsCollection
        if spaceIndex is None:
            topo = self.topo
            space = self.space
        else:
            if self.topo is None:
                topo = None
            else:
                topo = self.topo.getSubTopo(spaceIndex)
            space = dict(
                (it, fld[spaceIndex])
                for it, fld in self.space.iteritems())
        if timesIndex is None:
            times = self.times
        else:
            times = dict(
                (it, fld[timesIndex])
                for it, fld in self.times.iteritems())

        # ... prepare fields
        fields = dict()
        for it in fieldNames:
            # may raise KeyError or IndexError...
            arr = self.fields[it]

            # determine if and which tensor indexes can be considered
            if len(indexes) > len(arr.shape):
                raise IndexError(
                    "FieldsCollection.__getitem__ got %d indexes (%s, field"
                    " names excluded) but field %s has shape %s."
                    % (len(indexes), indexes, it, arr.shape))

            # rename this field if all tensor indexes are numbers
            if (tensorIndex and len(indexes) == len(arr.shape)
                and all(isinstance(idx, (int,long)) for idx in tensorIndex)):
                it = "%s_%s" % (it, "".join("%d" % (i+1) for i in tensorIndex))

            # subslicing the array
            # But first some numpy-background: if a multidimensional np.array
            # is accessed using index lists like:
            # >>> arr[[1,2],[5,6]]
            # you'll get a flat array [arr_15, ar_26]. But we expect to get
            # the 2nd and 3rd point at 6th and 7th frame. To get this result
            # from numpy you'll need to index the array like:
            # >>> arr[ [[1],[2]], [[5,6]] ]
            # or with a trick
            # >>> arr[np.ix_([1,2],[5,6])]
            # To make the np.ix_ function work properly with slices we'll need
            # to replace slices by an index-list (a) or we'll need to iterate
            # the array dimensions susequently (b)

            # (a) replace slices and None by indexList
            idxs = [spaceIndex, timesIndex,]
            for ii, idx in enumerate(idxs):
                if arr.shape[ii] == 1:
                    idxs[ii] = [0,]

                elif idx is None:
                    idxs[ii] = np.arange(arr.shape[ii])

                elif isinstance(idx, slice):
                    idxs[ii] = np.arange(arr.shape[ii])[idx]

            # reformat all indexes and get subslice of field
            idxs = np.ix_(*tuple(idxs))
            fields[it] = arr[idxs + tensorIndex]

        # ... create new FieldsCollection object

        result = type(self)(
            topo=topo,
            times=times, **fields)
        result.space.update(space)
        return result


    def getArrays(self, *fieldNames):
        """Returns a tuple containing the fields requested by the fieldNames
        argument, i.e. ndarrays.

        Usage:
         >>> U, V = flds.getArrays("U", "V")
         >>> W = U*V
         >>> flds["W"] = W
        """
        if len(fieldNames)>1:
            return (self.fields[x] for x in fieldNames)
        elif len(fieldNames):
            return self.fields[fieldNames[0]]
        else:
            raise ValueError(
                "FieldsCollection.getArrays() requires at least one field name"
                " as argument.")

    ###GP: don't know yet if we really want this...
    # def getGridShapedArray(self, fieldName=None):
    #     """self.topo must be a L{StructuredGrid}. Return the requested scalar
    #     field(s) (only its first time point) reshaped to three (x,y,z) space
    #     dimensions. I.e. result[i,j,k] would be the value of the specified
    #     field at the grid positions with i,j,k being the x,y,z-grid indexes.

    #     Usage:
    #      >>> Tgrid = flds.getGridShapedArray("T")
    #      >>> T_ij_row = Tgrid[i,j,:]

    #     Features:
    #      >>> flds["T"].shape
    #      (24000,1)
    #      >>> flds.topo.gridPtNb
    #      (20,30,40)
    #      >>> Tgrid = flds.getGridShapedArray("T")
    #      >>> Tgrid.shape
    #      (20,30,40)
    #      >>> flds["T"][:,0][:10] == Tgrid[:10,0,0]
    #      True

    #     If the field contains more than one time point then the first is
    #     taken.

    #     For fields having more than one times or non-scalar fields first
    #     get a scalar for a single time-point. E.g. the 10th time point of U_3:
    #      >>> U3grid = flds[["U"],:,9,2].getGridShapedArray()

    #     @param fieldName: (optional) Can be omitted if there is only one field
    #         in self.

    #     @return: an NxMxP array with (N,M,P) := self.topo.gridPtNb
    #     """
    #     if not isinstance(self.topo, StructuredGrid):
    #         raise ValueError(
    #             "FieldsCollection.getGridShapedArray called for fields on"
    #             " other topo than StructuredGrid.")
    #     if fieldName is None:
    #         fld = self.array[:,0]
    #     else:
    #         fld = self[fieldName][:,0]

    #     if len(fld.shape)>1:
    #         raise ValueError(
    #             "FieldsCollection.getGridShapedArray called for non-scalar"
    #             "field. (Suggestion: flds[['U'],:,9,2].getGridShapedArray())")

    #     # The order argument affects how fld is initially flattened. And at the
    #     # same time affects the order of the resulting indexes. In this case
    #     # fld is already flat so order="F" only makes sure that the data first
    #     # iterates over the x-index then y then z.
    #     return fld.reshape(self.topo.gridPtNb, order="F")

#} #End Set and get Fields


#{Operations on/using all stored fields
###GP: needs test:

    def filterEmptyPoints(self, filterValue=0, returnIndex=True):
        '''Filters empty pointData (e.g. sampling points above topography).
        A point will be removed if local values of each numeric field equals the
        filterValue. (* If filterValue is numeric/string, only numeric/string
        fields will be checked.)

        @param filterValue: float number; all points with this value will be
           filtered out

        @param returnIndex: if True: pointIndex will be returned, if False:
            a new FieldsCollection-object will be returned
        '''

        fieldNames = [
            key for key, array in self.fields.iteritems() if all(
                hasattr(array, attr) for attr in
                # python3: '__div__' --> '__truediv__'
                ['__add__', '__sub__', '__mul__', '__div__', '__pow__'])]

        nPts = self.getShape()[0]
        mask = np.ones((nPts, len(fieldNames)), dtype=bool)
        for ii, fieldName in enumerate(fieldNames):
            idx = (self[fieldName] == filterValue)
            mask[:,ii] = idx.reshape(nPts,-1).all(axis=1)

        mask = ~mask.all(axis=1)

        if returnIndex:
            return mask
        else:
            return self[mask]
#} #End Operations on/using all fields


#{ convert fields to/from scalars, vectors, tensors
# rearrange and/or modify scalar-, vector- and/or tensor fields
# stored in a FieldsCollection-object.
    def toScalars(self, fieldName):
        """Disassembles a Tensor- or VectorField to a series of ScalarFields
        holding the components.

        If the named field is already a scalar then returns self[fieldName]
        and ignores the store argument.

        Usage:
         >>> components = myField.toScalars('T')
         >>> # store components inplace of 'T'
         >>> myField.update(components)
         >>> del myField['T']

        @param fieldName: name of field.

        @returns: A new FieldsCollection object holding the scalar fields of
            all components.

        @raises KeyError: If the field name does not exist in self.fields.
        """
        field = self.fields[fieldName]  # may raise KeyError

        if len(field.shape) <= 2:
            # msg('%s is already a ScalarField.' % fieldName)
            result = self[[fieldName,]]
            return result

        tensorShape = field.shape[2:]

        # collect component fields
        fldArgs = dict()
        if len(tensorShape)==1:
            for i in range(tensorShape[0]):
                newName = "%s_%d" % (fieldName, i+1)
                newFld = field[:,:,i]
                fldArgs[newName] = newFld

        else:
            for i in range(tensorShape[0]):
                for j in range(tensorShape[1]):
                    newName = "%s_%d%d" % (fieldName, i+1, j+1)
                    newFld = field[:,:,i,j]
                    fldArgs[newName] = newFld


        # assemble resulting FieldsCollection
        result = type(self)(
            topo=self.topo, times=self.times, **fldArgs)
        result.space.update(self.space)

        return result


    @staticmethod
    def _guessTensorName(fieldNames, tensorName, default):
        """try to guess new name for vector or tensor field from components
        Service for toVector() and toTensor()
        """
        if tensorName is None:
            names = set(name.rsplit("_",1)[0] for name in fieldNames)
            if len(names)==1:
                tensorName = names.pop()
        if tensorName is None:
            tensorName = default
        return tensorName


    def toVector(self, fieldNames, vectorName=None):
        """Reshapes a sequence of scalar fields in self to a vector field

        @param fieldNames: sequence* of fieldNames to create the new
            vector field (*to use all vector-functions, it should be a triple)

        @param vectorName: optional name of the new VectorField in self. If not
            given then some heuristics try to guess, if that doesn't work then
            it defaults to "vector"

        @returns: a FieldsCollection with just the vector field
        """
        self.checkConsistency(fieldNames)

        # combine field data
        newDim = self.fields[fieldNames[0]].shape[:2] + (len(fieldNames),)
        newValues = np.empty(newDim)

        try:
            for kk, field in enumerate(fieldNames):
                newValues[:,:,kk] = self.fields[field]
        except ValueError as e:
            raise ValueError(
                "FieldsCollection.toVector() raised the following exception."
                " Could be caused by not all given component fields being"
                " scalars.\n%s" % e)

        # try to guess name
        vectorName = self._guessTensorName(fieldNames, vectorName, "vector")

        # assemble resulting FieldsCollection
        fldArgs = {vectorName: newValues}
        result = type(self)(
            topo=self.topo, times=self.times, **fldArgs)
        result.space.update(self.space)

        return result

    def toTensor(self, fieldNames, tensorName=None, order="A"):
        """Reshapes 3, 6 or 9 scalar fields to a tensor field

        Example:
         >>> T = myFields.toTensor(['S_1', 'S_2', 'S_3'], "S")

        @param fieldNames: sequence of 3, 6 or 9 fieldNames (see order) used to
            create the new tensor field

        @param tensorName: optional name of the new TensorField in self. If not
            given then some heuristics try to guess, if that doesn't work then
            it defaults to "tensor"

        @param order: Determines the association of the scalar fields to
            the component of the resulting 3x3-tensor ('A'==ABAQUS, default)
             1. diagonal
              - 'D'  [0,1,2]       --> diag([0,1,2])
             2. symmetric ordering
              - 'A'  [0,1,2,3,4,5] --> [[0 3 5] [3 1 4] [5 4 2]]
              - 'RS' [0,1,2,3,4,5] --> [[0 1 2] [1 3 4] [2 4 5]]
              - 'CS' [0,1,2,3,4,5] --> [[0 1 3] [1 2 4] [3 4 5]]
              - 'V'  [0,1,2,3,4,5] --> [[0 5 4] [5 1 3] [4 3 2]]
             3. non-symmetric ordering:
              - 'R'  [0,1,2,3,4,5,6,7,8] --> [[0 1 2] [3 4 5] [6 7 8]]
              - 'C'  [0,1,2,3,4,5,6,7,8] --> [[0 3 6] [1 4 7] [2 5 8]]

        @returns: a FieldsCollection with just the tensor field
        """
        order = order.upper()
        if len(fieldNames) == 3:
            order = 'D'
        elif len(fieldNames) == 6:
            if order in ['R','RS']:
                order = 'RS'
            elif order in ['C','CS']:
                order = 'CS'
            elif order not in ['V', 'A']:
                e = "%s is no valid order for a symmetric tensor" % order
                raise ValueError(e)

        elif len(fieldNames) == 9:
            if not (order in ['R','C']):
                e = "%s is no valid order for a non-symmetric tensor" % order
                raise ValueError(e)

        else:
            e = ("A tensor needs 3 (diag), 6 (symmetrical) or 9 scalar"
                 " fields to be defined. You passed %d." % len(fieldNames))

        # rearrange varLabels with respect to their position in tensor
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

        self.checkConsistency(fieldNames)

        nPts, nTimes = self.fields[fieldNames[0]].shape[:2]
        newDim = (nPts, nTimes, 3, 3)

        #stacking data to nPointsXnTimesX3X3 numpy array
        arrayTensor = np.zeros(newDim)

        ii = 0
        try:
            for iR in range(3):
                if order == 'D':
                    arrayTensor[:,:,iR,iR] = self.fields[fieldNames[ii]]
                    ii += 1
                else:
                    for iC in range(3):
                        arrayTensor[:,:,iR,iC] = self.fields[fieldNames[ii]]
                        ii += 1

        except ValueError as e:
            raise ValueError(
                "FieldsCollection.toTensor() raised the following exception."
                " Could be caused by not all given component fields being"
                " scalars.\n%s" % e)

        tensorName = self._guessTensorName(fieldNames, tensorName, "tensor")

        # assemble resulting FieldsCollection
        fldArgs = {tensorName: arrayTensor}
        result = type(self)(
            topo=self.topo, times=self.times, **fldArgs)
        result.space.update(self.space)

        return result

    def anglesToVector(self, angleNameA, angleNameB, vectorName="vector",
                       angleDef='bearing+plunge'):
        """ Calculates a (unit-)vector from a common 2-angle-definition
        specified by angleDef

        Example:
         >>> vU = myFields('bearing(U)', 'plunge(U)', "direction",
         ...               angleDef='bearing+plunge' )

        @param angleNameA: field name of a ScalarField holding the 1st angle

        @param angleNameB: field name of a ScalarField holding the 2nd angle

        @param vectorName: name of the new VectorField in self, defaults to
            "vector"

        @param angleDef: type of angle definition of angleA and angleB
            - 'bearing+plunge' (default)
            - 'strike+dip' (needs a check)
            - 'polar+azimut'

        @returns: a FieldsCollection with just the vector field

        @warning: angleDef = 'strike+dip' needs to be checked
        """

        from bae.field_02.numpy_geom import line, sph2cart

        angleDef = angleDef.lower()
        # getting angular fields
        angleA = self.fields[angleNameA]
        angleB = self.fields[angleNameB]

        # retransform to vector, depending on angular definition
        if angleDef == 'polar+azimut':
            n = np.array([np.cos(angleA) * np.sin(angleB),
                          np.sin(angleA) * np.sin(angleB),
                          np.cos(angleB)])

        elif angleDef in ['strike+dip', 'bearing+plunge']:
            if angleDef == 'strike+dip':
                # strike to bearing
                angleB = angleB - np.pi/2.

            # bearing+plunge to coords
            x,y,z = sph2cart(line(angleB, angleA))
            n = np.dstack((x, y, z))

        # assemble resulting FieldsCollection
        fldArgs = dict(vectorName=n)
        result = type(self)(
            topo=self.topo, times=self.times, **fldArgs)
        result.space.update(self.space)

        return result
