# -*- coding: utf-8 -*-
"""
Tested using:
 - numpy 1.16.2
 - scipy 1.2.1
"""

# creating new topo objects:
# ===========================
# A topo object has to store different time invariant fields called
# variables. These field values have to come as numpy-convertible
# iterables (list, tuples, np.arrays) or as objects which can be indexed
# similar as np.arrays (e.g. own objects with a decided __getitem__ function.

# For downstream compability reasons a new topo object needs to have a:
#   * __getitem__ function - which returns a sub topo object according to the
#   passed (point-, element-, ... ) index
#   * __eq__ fuction - to compare two topo objects that returns True/False
#   * _checkConsistency fuction - to check the 'length' of variables
# Try to subclass StaticEntities for new topo objects.

###GP:
__suggestions__ = """
Main questions:
- in unstructured mesh: nodeCoords and elNodes: two plain ndarrays labels and
  coords

    acces via label:
    self.coords[self.labels == 3214]
    self.getIndexFromLabel(attr/kind/type = "elem"/"node", iterable of labels):
      using np.intersect1d()

  ok: two plain ndarrays nodeLabels, coords; elemLabel, elNodes


- Gauss-pt topo: 1-D

- how do we establish an element series e.g. matching fields


further questions:
- make topo a package. subpackages: topo.pointcloud, topo.structuredgrid
  topo.mesh, topo.interpolation, possibly: topo.base. Main classes
  imported in topo/__init__.py
- topo classes common interface:
  . getPoints(): (nx3) array of point coordinates
  . __cmp__ or (__eq__ and __ne__): check if two grids are identical
  . __len__: nb of points
  . __str__: readable description for diagnostic purposes
  . getBoundingBox():


unstructuredmesh:
=================

- coords: (nx3) array
- elNodes:
  . have to be separated according to elType, but...
  . we want a sequential index e.g. for elsets
  . maybe: have one flat array and then views as parts of this with
    correct shape (nxm, m depending on elType)
    ... hmm doesn't help with the sequential index...
  . problem of dict {eltype: sub-elNodes}: has no defined order

  . have a list of elNodes arrays -- one for each el-type.
    For a global element index how do we determine the index in the right
    partial elNodes array?
    def getTyNoIndex(self, elemIds):
        "return (type-idx, elNodes-idx) tuple(s) for sequential elemIdx"

        if isinstance(elemIds, int):
            # selecting just one element
            offset = 0
            for typeId, elNodes in enumerate(typeElNodes):
                if elemIds < offset+len(elNodes):
                    return (typeId, elemIds-offset)
                offset += len(elNodes)
            raise IndexError("...")
        elif isinstance(elemIds, slice):
            # selecting a slice
            offset = 0
            for elNodes in typeElNodes:
                # check if slice intersects interval [ offset, offset+len(elNodes) )
                # if so add to resulting index-list
                offset += len(elNodes)
        else:
            # other iterable: list, tuple, ndarray of indexes

"""
###GP: end of suggestions


import numpy as np

from bae.log_01 import msg

from helper import (
                    # _checkConsistencyBefore,
                    _checkConsistencyAfter)

from numpy_geom import KDTree
import copy


############## Topography classes ########################################
#{Topology Classes

### Base class for Topo and Frame-Classes
class StaticEntities(object):
    """Base class for Topo and Frames and other array-like objects that have
    a consistent dimension.

    For setting new fields use self.addVariable(name,field). These fields will
    be indexed on calling the __getitem__-function. The names of these iterable
    fields will be stored in self.variables and the field-values will be
    attributes of StaticEntites.
    For setting entity-specific/unmutable values (like gridPtNb, origin,
    spacing) for StructuredGrid-object) use self.addFixed(name, entity). These
    values will eventually be copied unmodified to the new topo or frames
    object by the __getitem__-function.

    Usage:
    ======
        >>> se = StaticEntities(oldField=range(6))
        >>> se.addVariable('newField', [2**d for d in range(6)])
        >>> se.addFixed('fix', ('a','b',3))
        >>> print se.newField
        [ 1  2  4  8 16 32]
        >>> print len(se)
        6
        >>> for var in se.variables:
        >>>     print var, getattr(se, var)
        'oldField' [ 0 1 2 3 4 5]
        'newField' [ 1  2  4  8 16 32]
        >>> print se.fix
        ('a', 'b', 3)
        >>> print se[2:5].oldField, se[2:5].fix
        [2, 3, 4] ('a', 'b', 3)

    @ivar variables: list of attribute names. Only mutable attributes are
       listed here.
    @ivar fixed: list of immutable attributes.

    """
    ###ok: Will we ever really use the StaticEntities as such? Or only
    ###ok: subclasses? If only subclasses then better have different examples
    ###ok: examples of actual subclasses or refer to them...
    ###maybe: check if hiding some methods makes sense?
    def __init__(self, **kwargs):
        """Entity-fields can be passed as keyword-arguements

            >>> se = StaticEntities(fieldA=range(11), fieldB=range(11))

        """
        self.variables = []
        self.fixed = []

        for key, values in kwargs.iteritems():
            self.addVariable(key, values)


    def __getitem__(self, item):
        """Returns a sub-indexed copy of StaticEntities-object. Index-items
        can be either slices, indexlists or boolean index-arrays.

        For subclassing you may overload the
        L{_parseItem,<StaticEntities._parseItem>}-function. It is called before
        indexing and can be used to interpret the index-item in a suitable way
        or to check it first.

        @note: needs a cleanup since it got tricky handling true numpy-
        arrays and custom numpylike array-like objects
        """

        rtn = self.__class__()
        for fixed in self.fixed:
            value = copy.deepcopy(getattr(self, fixed))
            rtn.addFixed(fixed, value)

        item = self._parseItem(item)

        for eKey in self.variables:
            e = getattr(self, eKey)
            try:
                size = e.size
                shape = e.shape
            except AttributeError:
                #for non-numpy-arrays
                size = len(e)

            ee = e[item]
            if size:
                try:
                    # to preserve dimension of numpy-arrays
                    ee = ee.reshape((1,) + shape[1:])
                except:
                    ###GP: what's this for? can we accept wrong shape?
                    ###GP: or specify kind of exception
                    pass
                rtn.addVariable(eKey, ee)

        return rtn

    ###GP: can we get rid of this as well, like _checkVariable?
    def _parseItem(self, item):
        """This a convenience-function. Currently it does nothing but it can
        be specified for subclassing.
        It is called from L{__getitem__} before indexing the StaticEntities
        object can be specified for transforming or/and checking items.
        """
        ### example if first item is a list of 'frame-'strings
        # if isinstance(item, tuple) and isinstance(item[0], list):
        #    if all[isinstance(it, str) for it in item[0]]:
        #        item = [int(it[-3:] for it in item[0])]

        return item

    def __eq__(self, other):
        """To compare to StaticEntities-objects.

        >>> myPointCloudA == myPointCloudB  # True/False?
        """
        if not self.__class__ == other.__class__:
            # different object classes
            msg('different object classes:\n\t self: %s\n\t other: %s'
                % (self.__class__, other.__class__))
            return False

        if not self.fixed == other.fixed:
            # different set of registred data
            msg('different registered fixed entities:\n\t self:%s\n\t other:%s'
                % (self.fixed, other.fixed))
            return False

        def _checkEq(a,b):
            if not isinstance(a, type(b)):
                return False

            if isinstance(a, np.ndarray):
                if not a.shape == b.shape:
                    return False
                else:
                    return (a == b).all()
            else:
                return a == b

        for e in self.fixed:
            selfVar, otherVar = getattr(self, e), getattr(other, e)

            if not _checkEq(selfVar, otherVar):
                return False

        if not sorted(self.variables) == sorted(other.variables):
            # different set of registred data
            msg('different registered variables: \n\t self:%s \n\t other: %s'
                % (self.variables, other.variables))
            return False

        for e in self.variables:
            selfVar, otherVar = getattr(self, e), getattr(other, e)

            if not _checkEq(selfVar, otherVar):
                return False

        return True

    @_checkConsistencyAfter
    def addVariable(self, name, values):
        """Adds a mutable attribute to self, e.g. a topo field (i.e. point
        list) or a time series.

        @param values: must be a numpy array or convertible to such. Or
        it must have a __getitem__ and a __len__ function. __getitem__ must
        then accept arguments like a numpy-array.

        @Note: To implement checks overload this method in your child class.
           Do the check then call StaticEntities.addVariable(self, name, val).

        @Note: Silently overwrites an already existing attribute of the same
           name!
        """
        ###maybe: Do we need the option that values is something that is not
        ###maybe: convertible to ndarray?
        ###maybe: My suggestion: setattr(self, name, np.asarray(values))
        ###maybe: Con: one possible values argument would be a StaticEntity
        ###maybe: object.
        if isinstance(values, (np.ndarray, list, tuple)):
            values = np.asarray(values)

        elif not all(hasattr(values, aa) for aa in ("__len__", "__getitem__")):
            raise ValueError(
                "Values must be numpy-convertible iterables (list, tuples,"
                " np.arrays) or objects which can be indexed like np.arrays.")

        self.variables.append(name)
        setattr(self, name, values)

    def iterVariables(self):
        """Returns iterator for variableNames and variables.

        Example:
            >>> for name, field in myStaticEntity.iterVariables():
            >>>    print name, field[:5]
        """
        return ((n, getattr(self, n)) for n in self.variables)


    def addFixed(self, name, values):
        """Adds fixed entities to like (origin, spacing, gridPtNb) for
        StructuredGrid as attributes to TopoObject.
        These attributes will get copied to the returned object when
        __getitem__ is called.

        @Note: Silently overwrites an already existing attribute of the same
           name!
        """
        self.fixed.append(name)
        setattr(self, name, values)


    def iterFixed(self):
        """Returns iterator for fixedValueNames and fixedValues.

        Example:
            >>> for name, fixedAttr in myStaticEntity.iterFixed():
            >>>    print name, fixedAttr
        """
        return ((n, getattr(self, n)) for n in self.fixed)


    def _checkConsistency(self, allowEmpty=True):
        """Checks the inner dimensions of all stored/registered entity-fields.

        @note: cause of this function, all stored fields need to have a
            - shape-attribute or
            - __len__ function
        """

        def getDim(e):
            d = getattr(self, e)
            if d is None:
                return 0

            try:
                return d.shape[0]
            except AttributeError:
                return len(d)

        allDims = [getDim(e) for e in self.variables]
        validDims = [allDims[0] == d for d in allDims
                     if not (d == 0 and allowEmpty)]

        if not all(validDims):
            raise ValueError(
                'There is a Dimension missmatch in stored entities:\n%s'
                % "\n".join('Dim %s: %d' % (e,d)
                            for e, d in zip(self.variables, allDims)))

        return True


    def __len__(self):
        """Specifies length of StaticEntities-object

            >>> se = StaticEntities(testField=range(4))
            >>> print len(se)
            4

        """
        if self.variables == []:
            return 0
        else:
            e = getattr(self, self.variables[0])
            try:
                l = e.shape[0]
            except AttributeError:
                l = len(e)
            return l


###############################################################################

### Unstructured Point Cloud
class PointCloud(StaticEntities):
    """Class to manage unstructured points. Points is
    --at least-- specified by their coordinates (coords).
    The coordinates have to be unique within a tolerance tol.
    """

    def __init__(self, coords=None, **kwargs):
        """A PointCloud object can be created as an empty PointCloud (len will
        be zero) or by passing point dependent fields/arrays of the same
        length as keyword arguments. The coordinates (keyword: coords) can also
        be passed as positional argument.

        Example:
        ========
         >>> nPts = 5
         >>> # some points
         >>> x = np.arange(nPts)
         >>> y = np.arange(nPts) + 10
         >>> z = np.arange(nPts) + 100
         >>> coords = np.vstack((x, y, z)).T
         >>> # somePointIds
         >>> someIds = ['set%d' % s for s in range(nPts)]
         >>> points = PointCloud(coords, ids=someIds)
         >>> # index points
         >>> subPoints = points[1:3]
         >>> for id, coor in zip(subPoints.ids, subPoints.coords):
         >>>    print id, coor
        """
        ###ok: - docs on possible keywords: @kwarg coords: blabla
        ###ok: - don't store tolerance in object but have a tol argument
        ###ok:   (with default) for each method that needs it
        ###ok:   (like idxFromPoints)?
        ###GP:   but: addVariable still uses self.tol
        self._KDTree = None
        self._coordsHash = None
        self.tol = .1
        
        if coords is not None:
            kwargs['coords'] = coords
        StaticEntities.__init__(self, **kwargs)

        # for storring KDTree and hash (to check if update req) of coords
        


    ###GP: I think we might want to accept duplicate points here.
    ###GP: this would also remove the reqirement of global self.tol
    def addVariable(self, name, values):
        """Adds a mutable attribute to self, e.g. a topo field (i.e. point
        list) or a time series.

        If name=='coords' then check if coordinates are valid and unique within
        the tolerance (self.tol)

        @param values: must be a numpy array or convertible to such. Or
        it must have a __getitem__ and a __len__ function. __getitem__ must
        then accept arguments like a numpy-array.

        @Note: Silently overwrites an already existing attribute of the same
           name!
        """
        if name == 'coords':
            # 1st check if Coordinates are valid and unique
            values = np.asarray(values)
            kdTree = KDTree.KDTree(values)
            isUnique = KDTree.ensureUnique(kdTree, atol=self.tol)
            if not isUnique:
                raise ValueError('The coordinate sequence you passed is not'
                                 ' unique within the (absolute) tolerance %f.'
                                 % self.tol)
            else:
                # 2nd store KDTree and current hash
                self._KDTree = kdTree
                self._coordsHash = hash(values.tostring())

        StaticEntities.addVariable(self, name, values)


    def _coordsHasChanged(self):
        """Updates hash of coords"""
        return not(self._coordsHash == hash(self.coords.tostring()))


    def KDTree(self):
        """Creates or updates the KDTree of stored coords-data and store its
        current (coords) hash
        """
        if self._coordsHasChanged():
            # update
            self._coordsHash = hash(self.coords.tostring())
            self._KDTree = KDTree.KDTree(self.coords)

        return self._KDTree


    def idxFromPoints(self, points, atol=1E-6, shape='sphere',
                      strict=True):
        """Returns the indices of stored points which are close to
        the given points.

        Example
        =======
         >>> nPts=101
         >>> x = np.random.randn(nPts)
         >>> y = np.random.randn(nPts) + 10
         >>> z = np.random.randn(nPts) - 10
         >>> coords = np.vstack((x, y, z)).T
         >>> p = PointCloud(coords)
         >>> noiseMax = 0.1
         >>> noisyPoints = coords[[5,6]] + noiseMax*np.random.random((2,3))
         >>> ptIdx = p.idxFromPoints(noisyPoints,
         ...                         atol=np.sqrt(3)*noiseMax)

        @param points: set of point coordinates to get from self

        @param atol: absolute tolerance for ball query

        @param shape: search stored points in a 'sphere' or in a 'box' around
            each point in points.

        @param strict: if True (default) a one by one search will be applied.
            If the result is not complete and identically (all points found
            once) an ValueError will be raised

        @returns: indices of found points and list of lists (groups) of related
            stored points for each requested points. Indices and groups are the
            same if a search strict=True dosen't raise an error.

        @raises AttributeError: if self does not have a "coords" attribute.
        """
        import itertools as IT

        #update internal KDtree if needed
        self.KDTree()

        if shape == 'sphere':
            pNorm = 2.
        elif shape == 'box':
            pNorm = np.inf
        else:
            raise ValueError("Don't know a search shape %s" % str(shape))

        rstTree = KDTree.KDTree(points)
        coordsTree = self._KDTree

        groups = rstTree.query_ball_tree(coordsTree, atol, p = pNorm)

        if strict:
            ident = np.array( map(len, groups) )
            if not np.all(ident == 1):
                raise ValueError("Never found %d points and found %d points " %
                                 ((ident == 0).sum(), (ident > 1).sum()) +
                                 "multible times. If you don't a one-to-one " +
                                 "match, set strict=False")

        indices = np.unique(groups)

        if indices.dtype == object:
            # groups have different length --> dtypye == object, axis=0 fails
            indices = np.unique(list(IT.chain.from_iterable(groups)))
            # alternativly but slower:
            # indices = np.unique(np.concatenate(groups))

        return indices, groups


    def boundingBox(self):
        """Returns bounding box of stored coordinates

        @note: requires the attributes 'coords' to exist
        """
        return [self.coords.min(axis=0), self.coords.max(axis=0)]


###############################################################################
### Structured rectlinear Grid
class StructuredGrid(StaticEntities):
    """Represents a rectilinear grid which has regularly uniformly spaced
    points. It is defined by grid dimensions gridPtNb=(nx, ny, nz),
    grid origin=(x0, y0, z0) and point spacing=(dx, dy, dz). Thus all
    coordinates can be derived from the integer index (see below).

    The values of related fields can be stored in a sparse way. Only
    points with (valid) data will be stored and their postion is defined by
    this integer index. In case of a full/not sparse grid the integer index
    identifies the position in the (numpy) array.

    Nomenclature:
     - integer index: plain integer number starting at zero. This is the
       numpy index to access a particular value in a corresponding field
       (if it's not sparse).
     - grid index triple: a tuple of three zero-based indexes identifying
       the position in the grid in x-, y-, z-direction.

    The integer index sequence (order of items in corresponding fields)
    corresponds to iterating first over the x-index then y then z.
    """

    def __init__(self, *args, **kwargs):
        """
        A StructuredGrid-object can be constructed
            - as an empty object
            - by specifying the grid characteristics in the order gridPtNb,
              origin, spacing
            - by specifying a sparse (or complete) index as a last argument
              or as a keywordargument

        Example:
        ========
            >>> gridPtNb = (10,11,12)
            >>> orignin = (0.0, 100.0, 2000.0)
            >>> spacing = (101, 102, 103)
            >>> knownValIdx = [120, 177,100,1001]
            >>> sparseGrid = StructuredGrid(knownValIdx,
            ...                             gridPtNb, origin, spacing)
        """

        StaticEntities.__init__(self)

        if not len(args) in [0, 3, 4]:
            raise TypeError(
                ('__init__ of StructuredGrid no, 3 (gridPtNb, origin,'
                 ' spacing), or 4 (..., idx) positional arguments (%d given)')
                % len(args))

        # gridParameters if present
        for ii, key in enumerate(['gridPtNb', 'origin', 'spacing']):
            if len(args) > 0:
                val = args[ii]
            else:
                val = kwargs.pop(key, None)

            self.addFixed(key, val)

        # gridPoint-idx if present
        if len(args) == 4:
            idx = args[3]
        else:
            idx = kwargs.pop('idx', None)

        if idx is not None:
            self.addVariable('idx', idx)

        # use linear list if gridPtNb is defined but idx is not
        if idx is None and self.gridPtNb is not None:
            linRange = np.arange(self.nFullPts())
            self.addVariable('idx', linRange)

        self._checkConsistency()

    def _checkConsistency(self, allowEmpty=True):
        """Tests proper grid definition
        """
        test = [self.gridPtNb, self.spacing, self.origin]

        err = False
        if all(a is None for a in test):  # not all unset
            pass
        else:
            if any(len(a) != 3 for a in test):  # all 3-Tuples
                err = True

            if not all([aa > 0 for a in test[:2] for aa in a]):
                err = True  # nbs,spaces >0

            if not all([isinstance(b, int) for b in test[0]]):  # nbs are int
                err = True

        if err:
            msg("Some of the grid parameters are not defined properly.\n"
                "gridPtNb: " + str(self.gridPtNb) + " ? (uint,uint,uint)\n"
                "spacing: " + str(self.spacing) + " ? (f+,f+,f+)\n"
                "origin:" + str(self.origin) + " ? (f,f,f)")
            raise ValueError("Grid not properly defined.")

        # check maxIndex of idx
        if hasattr(self, 'idx'):
            if self.idx.max() + 1 > self.nFullPts():
                raise ValueError('Some integer indexes (max=%d) are larger'
                                 ' than the total number of grid points (%d).'
                                 % (self.idx.max(), self.nFullPts()))

        # check length of all stored staticFields
        StaticEntities._checkConsistency(self)


    ###GP: rename to __len__ ?
    def nFullPts(self):
        """Returns the point number of the fully filled grid"""
        if self.gridPtNb is not None:
            return self.gridPtNb[0] * self.gridPtNb[1] * self.gridPtNb[2]
        else:
            return 0


    ###GP: rename flat -> grid index triple
    def idxToFlatIdx(self):
        """Unravels integer indexes to a sequence of grid index triples

        @note: requires the attributes 'gridPtNb', 'origin' and
        'spacing' to exist
        """
        return np.unravel_index(self.idx, self.gridPtNb, order='F')


    def flatIdxToCoordComp(self, idx, component):
        """Returns the coordinate component from grid index triple.

        @param idx: integer index/-es of the grid point(s)
        @param component: 0 for x-coordinate, 1 for y and 2 for z

        @note: requires the attributes 'gridPtNb', 'origin' and
        'spacing' to exist
        """
        return idx * self.spacing[component] + self.origin[component]


    def idxToCoord(self, idx):
        """Returns the coordinates from integer indexes
        """
        coord = np.vstack([self.flatIdxToCoordComp(cIdx, ii)
                   for ii, cIdx in enumerate(self.idxToFlatIdx())]).T
        return coord


    def coordToIdx(self, coords, bounds=('clip', 'wrap')):
        """Returns the integer index from a list of coordinates with. The
        bounds kwarg specifies behavior for coordinates which are not inside
        the grid.

        @param coords: list of coordinates

        @param bounds: specifies how out-of-bounds indices are handled. Can
            specify either one mode or a tuple of modes, one mode per index.
                - raise - raise an error
                - wrap - wrap around
                - clip - clip to the range

            In 'clip' mode, a negative index which would normally wrap will
            clip to 0 instead.
            Default is ('clip', 'wrap'), so outside coordinates will get
            indexed at the 'nearest' bound.

        @returns: flattened index list in fortran order

        @note: requires the attributes 'gridPtNb', 'origin' and
        'spacing' to exist
        """
        coords = np.asarray(coords)
        origin = np.asarray(self.origin)
        gridPtNb =np.asarray(self.gridPtNb)
        spacing = np.asarray(self.spacing)

        # grid index triples
        fIdx = np.round((coords - origin) / spacing).astype(int)

        # ravel to integer indexes
        tIdx = np.ravel_multi_index(fIdx, tuple(gridPtNb),
                                    mode=bounds, order='F')

        return tIdx


    def coords(self):
        """Returns coordinates of stored grid points

        @note: requires the attributes 'idx', 'gridPtNb', 'origin' and
        'spacing' to exist
        """
        return self.idxToCoord(self.idx)

    def flatX(self):
        """Returns (full) list of x-coordinates of dimension
        self.gridPtNb[0]

        @note: requires the attributes 'gridPtNb', 'origin' and
        'spacing' to exist
        """
        return self.flatIdxToCoordComp(np.arange(self.gridPtNb[0]), 0)

    def flatY(self):
        """Returns (full) list of y-coordinates of dimension
        self.gridPtNb[1]

        @note: requires the attributes 'gridPtNb', 'origin' and
        'spacing' to exist
        """
        return self.flatIdxToCoordComp(np.arange(self.gridPtNb[1]), 1)


    def flatZ(self):
        """Returns (full) list of z-coordinates of dimension
        self.gridPtNb[3]

        @note: requires the attributes 'gridPtNb', 'origin' and
        'spacing' to exist
        """
        return self.flatIdxToCoordComp(np.arange(self.gridPtNb[2]), 2)


    def fullXYZGrid(self):
        """Returns the numpy.meshgrid fields of coordinates

        @returns: (X,Y,Z), each of dimension self.gridPtNb

        @note: requires the attributes 'gridPtNb', 'origin' and
        'spacing' to exist
        """
        return np.meshgrid(self.flatX(), self.flatY(), self.flatZ(),
                           indexing='ij')


    def fullValuesGrid(self, values, emptyValues=np.nan, rollPositionAxis=True):
        """Arranges field values (flat) on given topo to full 3D X-Y-Z field
        as created from numpy.meshgrid.

        Some usefull scipy and numpy fuctions assume the last 3 array-dimension
        to represent the coordinates, while our fields framework referes to the
        first dimension. Thus, by default, the positional axis will be rolled
        to last dimensions. Set rollPositionAxis=False to avoid this behaviour.

        Example:
        ========
         >>> print 'type of myField: ', myField
         type of myField: <type 'ScalarField'>
         >>> print 'sparse field values shape: ', myField.values.shape
         sparse field values shape: (100,3)
         >>> print 'grid size: ', myField.topo.gridPtNb
         grid size: (101,102,103)
         >>> npArrA = myField.topo.fullValuesGrid(myField)
         >>> npArrB = myField.topo.fullValuesGrid(myField, rollPositionAxis=False)
         >>> print 'npArrA.shape=',npArrA.shape
         npArrA.shape=(3,101,102,103)
         >>> print 'npArrB.shape=',npArrB.shape
         npArrA.shape=(101,102,103,3)

        @param values: can be a L{Scalar<bae.fieldData_00.fields.ScalarField>}-,
            L{Vector<bae.fieldData_00.fields.VectorField>}- or
            L{Tensor<bae.fieldData_00.fields.TensorField>} or a numpy-array of
            shape (nPts, nFrames, :, :), where nPts has to be len(self).

        @param emptyValues: value to fill undefined points

        @param rollPositionAxis: if true the positional axes will get rolled to
            last array dimension (as expected by many scipy/numpy functions)

        @returns: numpy.meshgrid-like array

        @note: requires the attributes gridPtNb', 'origin' and 'spacing'
        to exist
        """
        if hasattr(values, 'topo') and hasattr(values, 'values'):
            if (values.topo is not None) and not (self == values.topo):
                raise ValueError('The topo-object of the fields-object you'
                                 ' passed differ.')
            else:
                values = values.values


        values = np.asarray(values)
        if not values.shape[0] == len(self):
            raise ValueError(
                "The first dimension of the values array (shape %s) doesn't"
                " match the length of topo-object %d." % (values.shape, len(self)))

        # create empty field
        out = emptyValues * np.empty((self.nFullPts(),) + values.shape[1:])

        # set values to ozt-field at stored pointIdx
        out[self.idx] = values

        if rollPositionAxis:
            #roll pointIndex to last dimension
            out = np.rollaxis(out, -1)
            out = out.reshape(out.shape[:-1] + tuple(self.gridPtNb),
                              order='F')

        else:
            out = out.reshape(tuple(self.gridPtNb) + out.shape[1:],
                              order='F')

        return out


#} #End Topology Classes


###############################################################################

###ok: tests go in bae_utils/test_bae_package/test_fieldData_00.py
### TestFuctions
def _someTests():
    gridA = StructuredGrid()

    gridB = StructuredGrid((10,11,12), (1,2,3), (.1,.2,.3))

    idx = np.arange(10*11*12)

    gridC = StructuredGrid((10,11,12), (1,2,3), (.1,.2,.3), idx)

    idxD = np.arange(9*11*10)
    gridD = StructuredGrid((10,11,12), (1,2,3), (.1,.2,.3), idxD)

    # check comparison
    se = StaticEntities(oldField=range(6))
    se.addFixed('fix', ('a','b',3))
    se.addFixed('bla', np.array([1,2]))
    sb = copy.deepcopy(se)
    print se == sb

    sb.addFixed('blub',(None,))
    print se == sb

    sb = copy.deepcopy(se)
    sb.bla = 'q'
    print se == sb

    sb = copy.deepcopy(se)
    sb.bla = np.array([1,2])
    print se == sb


if __name__ == '__main__':
    print 'No syntax errors.'

