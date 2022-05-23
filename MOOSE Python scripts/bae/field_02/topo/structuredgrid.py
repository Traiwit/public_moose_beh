# -*- coding: utf-8 -*-
"""Topology classes for structured grids.
"""

import numpy as np

from .base import TopoBase
from .pointcloud import PointCloud
from bae.log_01 import msg

###############################################################################
#{Topology Classes

###GP: rotated option missing completely, yet
###GP: check for other functions to be included:
###GP:   - fieldData_00.StructuredGrid
###GP:   - mesh_01.MeshStructuredPoints, ...Rot
class StructuredGrid(TopoBase):
    """Represents a rectilinear grid which has regularly uniformly spaced
    points. It is defined by grid dimensions gridPtNb=(nx, ny, nz),
    grid origin=(x0, y0, z0) and point spacing=(dx, dy, dz). Thus all
    coordinates can be derived from the integer index (see below).

    Usage:
     >>> from bae.field_02 import StructuredGrid
     >>> gridPtNb = (10,11,12)
     >>> origin = (0.0, 100.0, 2000.0)
     >>> spacing = (101, 102, 103)
     >>> grid = StructuredGrid(
     ...     origin=origin, spacing=spacing, gridPtNb=gridPtNb)

    ...or...
     >>> from bae.field_02 import StructuredGrid
     >>> box = [[-10.0, 22.5, 83.0], [50.0, 77.5, 200.0]]
     >>> topo = StructuredGrid.fromBox(box[0], box[1], spacing=5)

    B{Indexes:}
     - serial integer index: plain integer number starting at zero. This is the
       numpy index to access a particular value in a corresponding field
       (if it's not sparse).
     - grid index triple: a tuple of three zero-based indexes identifying
       the position in the grid in x-, y-, z-direction.

    The serial integer index (order of items in corresponding fields)
    corresponds to iterating first over the x-index then y then z. This is
    corresponding to the order in vtk files.

    To determine the serial integer index for a given grid index triple do the
    following. To calculate many indexes at once make the input parameter
    gridIndexTuple a triple of ndarrays for the i, j, k indexes (obviously all
    of the same shape). Then you'll get a ndarray of serial integer indexes of
    the same shape.
     >>> gridIndexTuple = ([0,3], [1,0], [4,2])  # (i-indexes, j-id's, k-id's)
     >>> np.ravel_multi_index(gridIndexTuple, self.gridPtNb, order="F")

    To determine the grid index tuple for serial integer indexes --given in an
    array serialIndex-- do the following. The result will be a triple (x-ids,
    y-ids, z-ids) of ndarrays of the same shape as the input array serialIndex.
    This result can be used to directly access a field stored in "grid shape":
     >>> gridIds = np.unravel_index(serialIndex, self.gridPtNb, order="F")
     >>> print gridIds
     ([4,1,...], [5,3,...], [1,9,...])

    If you prefer a list of index-triples then you want the transpose of this:
     >>> gridIds = np.transpose(
     ...     np.unravel_index(serialIndex, self.gridPtNb, order="F"))
     >>> print gridIds
     array([[4,5,1], [1,3,9], ...])

    B{Indexes and "serial shaped" vs "grid shaped" fields:}
    Fields in a FieldsCollection have a "serial" shape and are accessed with a
    serial index: There is only one array-dimension to identify the location.
    I.e. a scalar has a two-dimensional shape, in this particular case we have
    many point locations and one point in time:
     >>> flds = FieldsCollection(....)
     >>> flds["T"].shape
     (185856, 1)

    You can convert the shape to a "grid shape" so that you'd use a grid index
    to access certain locations:
     >>> T_serial = flds["T"][:,0]
     >>> T_grid = T_serial.reshape(flds.topo.gridPtNb, order="F")
     >>> T_grid.shape
     (64, 24, 121)

    B{Strides:}
    The offsets for the serial index of the next grid point in x,y,z
    direction are commonly called strides. They can be calculated like so:
     >>> strides = np.array(
     >>>     [1, self.gridPtNb[0], self.gridPtNb[0]*self.gridPtNb[1]])

    (Note that ndarray.strides in contrast give the *byte*-offsets along the
    respective axes.)

    For the grid point with the serial index i the neighbour grid point in
    +x (+y,+z) direction would be i+strides[0] (i+strides[1], i+strides[2]).

    The serial index of a point with grid indexes gridIndex=(i,j,k) can then
    be calculated like so:
     >>> seriaIndex = sum( i*s for i, s in zip(gridIndex, strides) )

    @ivar gridPtNb: A triple of integers specifying the number of grid
      points in x, y, z direction
    @ivar origin: First point of the grid (min x, y, z)
    @ivar spacing: A triple of floats: Point spacing in x, y, z direction
    @cvar position: "structPt". For L{FieldsCollection}.topo objects the
      position attribute identifies the position of values ("node", "element",
      "elemIP", "point" or "structPt").
    """

    position = "structPt"

    def __init__(self, *args, **kwargs):
        """
        A StructuredGrid-object can be constructed
            - from another StructuredGrid object (copy constructor)
            - by specifying the grid characteristics origin, spacing, gridPtNb
              as keyword arguments

        Example:
        ========
         >>> gridPtNb = (10,11,12)
         >>> origin = (0.0, 100.0, 2000.0)
         >>> spacing = (101, 102, 103)
         >>> sparseGrid = StructuredGrid(
         ...     origin=origin, spacing=spacing, gridPtNb=gridPtNb)

        @param args: can take one positional argument only: another
          StructuredGrid
        @kwarg origin: First point of the grid (min x, y, z).
        @kwarg spacing: A triple of floats: Point spacing in x, y, z direction
          Or just a single number that will then be used for x, y and z
          directions.
        @kwarg gridPtNb: A triple of integers specifying the number of grid
          points in x, y, z direction.
        """

        if len(args)==3 and not kwargs:

            # For fromXXX() class methods only: take origin, spacing, gridPtNb
            # This code path is not meant for external use!
            # For efficiency we don't want to add checks or conversions here!
            # And we don't make a copy of the args here!
            self.origin = args[0]
            self.spacing = args[1]
            self.gridPtNb = args[2]

        elif len(args)==1 and not kwargs:

            # copy constructor:
            init = args[0]
            if isinstance(init, StructuredGridSparse):
                raise NotImplementedError()
            elif isinstance(init, StructuredGrid):
                # make a deep copy!
                self.origin = np.array(init.origin)
                self.spacing = np.array(init.spacing)
                self.gridPtNb = np.array(init.gridPtNb)

        elif not args and set(kwargs)==set(("origin", "spacing", "gridPtNb")):
            # grid characteristics origin, spacing, gridPtNb

            # make a deep copy!
            self.origin = np.array(kwargs["origin"], dtype=float)
            spacing = kwargs["spacing"]
            if isinstance(spacing, (int, long, float)):
                spacing = [spacing,spacing,spacing]
            self.spacing = np.array(spacing, dtype=float)
            self.gridPtNb = np.array(kwargs["gridPtNb"], dtype=int)

        elif not args and set(kwargs)==set(
                ("firstPoint", "lastPoint", "spacing")):
            # grid characteristics firstPoint (=origin), lastPoint, spacing

            # make a deep copy!
            self.origin = np.array(kwargs["firstPoint"], dtype=float)

            spacing = kwargs["spacing"]
            if isinstance(spacing, (int, long, float)):
                spacing = [spacing,spacing,spacing]
            self.spacing = np.array(spacing, dtype=float)

            gridBoxSize = (np.array(kwargs["lastPoint"], dtype=float)
                           - np.array(kwargs["firstPoint"], dtype=float))
            self.gridPtNb = np.rint(gridBoxSize/self.spacing).astype(int) + 1

        else:
            # invalid arguments
            raise ValueError(
                "Improper arguments to constructor of StructuredGrid:\nargs=%s"
                " kwargs=%s" % (args, kwargs))

        if not all(x>0 for x in self.gridPtNb):
            raise ValueError('GridSizes have to be positive integers:\n'
                             'gridPtNb: %s' % str(self.gridPtNb))

        if not all(x>0 for x in self.spacing):
            raise ValueError('Spacesizes have to be positive:\n'
                             'spacing: %s' % str(self.spacing))

    def __nonzero__(self):
        """Checks if L{StructuredGrid} is valid or not emtpy:
            - spacing, gridPtNb and origin are set and hold three items each
            - all spacings and gridPtNbs are different from 0

            Example:
                >>> if not myStructGrid:
                >>>     raise ValueError('StructuredGrid is not set properly.'
                >>>                      + str(myStructGrid))
        """
        tt = [self.spacing, self.gridPtNb, self.origin]

        if not all([hasattr(t, '__iter__') for t in tt]):
            return False

        if not all([len(m)==3 for m in tt]):
            return False

        if not all(self.spacing) or not all(self.gridPtNb):
            return False

        return True

    def __eq__(self, other):
        """Compares L{StructuredGrid}instance with another object
            Example:
                >>> if not myGridA == myGridB:
                >>>     raise ValueError('StructuredGrids do not have equal'
                ...                      ' topology: %s vs. %s' %
                ...                      (myGridA, myGridB))
        """
        if not isinstance(other, type(self)):
            return False

        ts = [self.spacing, self.gridPtNb, self.origin]
        to = [other.spacing, other.gridPtNb, other.origin]

        def same(a,b):
            """ To handle diffent (non)enumerables """
            # both are iterable or both are not iterable
            if not (hasattr(a, '__iter__') == hasattr(b, '__iter__')):
                return False

            if hasattr(a, '__iter__'):
                if not len(a) == len(b):
                    return False
                if not all(ia==ib for ia, ib in zip(a,b)):
                    return False
                return True
            else:
                try:
                    return bool(a==b)
                except:
                    return False

        return all([same(s,o) for s,o in zip(ts,to)])


    def __len__(self):
        """len(mygrid) returns the number of grid points"""
        return self.gridPtNb[0]*self.gridPtNb[1]*self.gridPtNb[2]

    def __str__(self):
        """return a descriptive string of self.
        """
        return (
            "StructuredGrid %d points: %d, %d, %d points in"
            " x,y,z-direction, origin %s, spacing %s"
            % (len(self), self.gridPtNb[0], self.gridPtNb[1], self.gridPtNb[2],
               self.origin, self.spacing))

    @classmethod
    def fromBox(cls, firstPoint=None, lastPoint=None,
                spacing=None, gridPtNb=None):
        """Create a new StructuredGrid object from two corner points --i.e. a
        box-- and spacing or alternatively gridPtNb.

        Note that if spacing is given then lastPoint will be taken
        approximately and adjusted according to the other parameters.

        Example:
        ========
         >>> grid = StructuredGrid.fromBox(box[0], box[1], spacing=5)
        """
        # make spacing a new ndarray if given
        if isinstance(spacing, (int, long, float)):
            spacing = np.array([spacing,spacing,spacing], dtype=float)
        elif spacing is not None:
            spacing = np.asarray(spacing, dtype=float)  # still make a copy!

        # prepare origin
        origin = np.array(firstPoint, dtype=float)
        gridBoxSize = np.array(lastPoint) - origin

        # prepare gridPtNb or spacing from the other
        if gridPtNb is None and spacing is not None:
            gridPtNb = np.array(
                [int(round(x/dx))+1 for x, dx in zip(gridBoxSize, spacing)],
                dtype=int)
        elif spacing is None and gridPtNb is not None:
            gridPtNb = np.array(gridPtNb, dtype=int)
            spacing = gridBoxSize / gridPtNb
        else:
            raise ValueError(
                "StructuredGrid.fromBox() requires either <spacing> or"
                " <gridPtNb> parameter. We got:\nspacing=%s, gridPtNb=%s"
                % (repr(spacing), repr(gridPtNb)))

        # create the new object
        return cls(origin, spacing, gridPtNb)

    def getSubTopo(self, index):
        """Service method for FieldsCollection.__getitem__() selecting from the
        space dimension.

        Returns a L{StructuredGridSparse} --or if possible a StructuredGrid--
        of only the specified points.

        @param index: A serial integer index --anything suitable to have numpy
           choose from one dimension.
        """
        ###GP: check if structgrid is possible, sparsegrid implementation missing
        if isinstance(index, (int, long)):
            index = [index,]
        elif isinstance(index, slice):
            index = np.arange(len(self), dtype=int)[index]
        return PointCloud(self.getPoints(index))

    def __getitem__(self, index):
        """Return a L{StrucutredGridSparse} --or if possible a StrucutredGrid--
        of only the specified points.

        @param index: A serial integer index --anything suitable to have numpy
           choose from one dimension.

           Docs missing: Apparently can also take other types of indexes....
        """
        ###GP: Docs missing
        if not isinstance(index, tuple):
            index = (index,)

        slices = 3*[slice(None),]
        if all(isinstance(idx, slice) for idx in index):
            for ii, sl in enumerate(index):
                slices[ii] = sl

            newO, newPtNb, newSp = (), (), ()
            for ii,(nPt, sl) in enumerate(zip(self.gridPtNb, slices)):
                startI, stopI, stepI = sl.start, sl.stop, sl.step
                if startI is None:
                    startI = 0
                if stopI is None:
                    stopI = nPt
                if stepI is None:
                    stepI = 1

                if startI < -nPt or startI > nPt:
                    print 'startIndex: ', startI
                    raise IndexError('Index (dim %d) is out of bounds' % nPt)

                if stopI < -nPt or stopI > nPt-1:
                    print 'endIndex: ', stopI
                    raise IndexError('Index (dim %d) is out of bounds' % ii)

                if startI < 0:
                    startI = nPt + startI
                if stopI < 0:
                    stopI = nPt + stopI

                newPtNb += (int((stopI - startI) / stepI),)
                newO    += (self.origin[ii] + startI * self.spacing[ii],)
                newSp   += (stepI * self.spacing[ii],)

            return type(self)(origin=newO, spacing=newSp, gridPtNb=newPtNb)

        elif len(index) == 1 and map(int, np.__version__.split(".")) < [1, 13]:

            # first create a 3D-mask of indexes on current grid
            toCheck = np.zeros(len(self))
            toCheck[index[0]] = True
            toCheck = self.serialToGrid(toCheck)

            # now get indexed ijk-triples
            where = np.argwhere(toCheck)

            # differences in each dir should only be
            #   - the step and/or a
            #   - negative multiple of it
            #   - zero
            diff = np.diff(where, axis=0)
            startsI = where[0]
            stopsI = where.max(axis=0)
            #print stopsI

            uni = np.unique(diff, axis=0)
            if len(uni) <= 3:
                stepsI = uni.max(axis=0)
                slices = ()
                for start,stop,step in zip(startsI, stopsI, stepsI):
                    slices += (slice(start, stop+step, step),)

                return self[slices]

        raise NotImplementedError("Returning a SparseStructuredGrid is not"
                                  " implemented yet.")

    def getGridDataString(self):
        """Return a Python command string that can be used as arguments for the
        constructor to rebuild this object.

         >>> kwargStr = grid.getGridDataString()
         >>> print kwargStr
         dict(
             origin=[10.00,3.00,1.00],
             spacing=[5,5,5],
             gridPtNb=[20,20,2])
         >>> kwargs = eval(kwargStr)
         >>> mesh = StructuredGrid(**kwargs)
        """
        gridDataString = (
            'dict(\n'
            '    origin=[%.2f,%.2f,%.2f],\n'
            '    spacing=[%s,%s,%s],\n'
            '    gridPtNb=[%d,%d,%d])'
            % (self.origin[0], self.origin[1], self.origin[2],
               self.spacing[0], self.spacing[1], self.spacing[2],
               self.gridPtNb[0], self.gridPtNb[1], self.gridPtNb[2]))
        return gridDataString

    def getPoints(self, index=None):
        """return Nx3 array of point coordinates

        @param index: If not given (None) then return the coordinates of *all*
           points. Otherwise it's got to be an array of indexes. The indexes
           might be serial integer indexes (len(index.shape)==1) or grid indexes
           (index.shape==(N, 3))

           A serial integer index reflects the order of items in corresponding
           fields. It corresponds to iterating first over the x-index then y
           then z. This is also the order in vtk files.

           Grid indexes are triples of zero based integers. The first (x-)
           index changes fastest.

           index can be any iterable convertible to a numpy array (or None).
        """
        if index is None:
            # index not given, return all points
            gridCoords = tuple([np.arange(nPt)[:np.newaxis]*spacing + origin
                                for spacing, nPt, origin in
                                zip(self.spacing, self.gridPtNb, self.origin)])

            X,Y,Z = np.meshgrid(*gridCoords, indexing='ij')
            xyz = np.vstack((X.ravel("F"),Y.ravel("F"),Z.ravel("F"))).T
        else:
            index = np.asarray(index)
            if len(index.shape)==2 and index.shape[1]==3:
                # index is grid index
                # convert possible negative indexes
                index = np.where(index<0, index+self.gridPtNb, index)
                if np.any(index<0) or np.any(index>=self.gridPtNb):
                    raise IndexError(
                        "StructuredGrid.getPoints: grid index out of bounds:"
                        " index=%s, self.gridPtNb=%s." % (index, self.gridPtNb))
            else:
                # index is a serial index
                # convert possible negative indexes
                N = len(self)
                index = np.where(index<0, index+N, index)
                if np.any(index<0) or np.any(index>=N):
                    raise IndexError(
                        "StructuredGrid.getPoints: serial index out of bounds:"
                        " index=%s, len(self)=%s." % (index, N))
                # convert index to grid index
                index = np.array(np.unravel_index(index, self.gridPtNb[::-1])
                                 [::-1]).transpose()
            # calculate coordinates
            xyz = self.origin + index*self.spacing

        return xyz

    ###GP: The methods _gridToSerial might be converted to methods of the
    ###GP: fieldsCollection class: e.g. fromGridField, asGridField
    ###GP: otherwise just use numpy.ndarray.reshape
    def reshapeGridToSerial(self, gridField):
        """Reshapes the gridField argument.
        Returns the input reshaped with the first three dimensions flattened.

        After the conversion...
         >>> serialField = self.reshapeGridToSerial(gridField)

        ... for any grid index (i,j,k) and the corresponding serial index si,
        gridField[i,j,k] should be equal to serialField[si].

        @param gridField: a numpy array of shape (N,M,P,...) with (N,M,P)==
        self.gridPtNb.

        @Returns: the input reshaped with the first three dimensions i,j,k
        flattened.

        @Note:
        Should be identical to gridField.reshape(-1, gridField.shape[3:]).
        For scalar fields: gridField.reshape(-1) or gridField.flatten() (always
        yields a copy) or gridField.ravel().
        """

        # check inputShape
        if not gridField.shape[:3] == tuple(self.gridPtNb):
            raise IndexError(
                "StructuredGrid.gridToSerial got a field shaped %s which is"
                " incompatible with topo.gridPtNb %s."
                % (gridField.shape, self.gridPtNb))

        newShape = (-1,) + gridField.shape[3:]
        return gridField.reshape(newShape)

    def reshapeSerialToGrid(self, serialField):
        """Reshapes the serialField argument. The first dimension must match
        the length of self (i.e. the number of grid points). Then this first
        dimension will be reshaped to three dimensions according to
        self.gridPtNb.

        After the conversion...
         >>> gridField = self.reshapeSerialToGrid(serialField)

        ... for any serial index si and the corresponding grid index (i,j,k),
        serialField[si] should be equal to gridField[i,j,k].

        @param serialField: a numpy array of shape (N,...) with N=len(self)

        @Note:
        Should be identical to serialField.reshape(tuple(grid.gridPtNb)) for
        scalar fields and
        serialField.reshape(tuple(grid.gridPtNb)+serialField.shape[1:]) for
        vector/tensor fields.
        """

        # check input shape
        if serialField.shape[0] != len(self):
            raise IndexError(
                "StructuredGrid.serialToGrid got a field shaped %s which is"
                " incompatible with the topo of length %d. The first dimension"
                " of the field must match the length of the topo."
                % (serialField.shape, len(self)))

        newShape = tuple(self.gridPtNb) + serialField.shape[1:]
        return serialField.reshape(newShape)

    def getIdsInBox(self, box):
        """Returns an iterable (currently implemented: a list) of serial
        integer indexes for all points inside the given box. Points on the
        border included.

        @param box: [[min x, min y, min z], [max x, max y, max z]]. If self is
            a L{MeshStructuredPointsRot}-object then the given box coords must
            be backrotated into the rotated coordinate system. This is exactly
            what we have in the use case of L{Mesh.getElemCoordsForGrid()}.
            So don't overload this function in the MeshStructuredPointsRot
            class.

        @Note: WARNING: uses the rotated coordinate system for
            L{MeshStructuredPointsRot}-objects. And this is the intention!

            ###GP: rotation feature needs to be checked: I think the general
            behaviour should be that the box is in global coords. But we need
            this alternative for something important, but I forgot for what...
        """
        box = np.asarray(box)

        # determine grid index boundaries
        ximin = np.ceil((box[0] - self.origin) / self.spacing).astype(int)
        ximax = np.floor((box[1] - self.origin) / self.spacing).astype(int)

        ximin = np.maximum(ximin, 0)
        ximax = np.minimum(ximax+1, self.gridPtNb)

        # sequential indexes for the whole grid
        ijk = np.arange(np.prod(self.gridPtNb)).reshape(self.gridPtNb)

        # apply index boundaries
        ijk = ijk[ximin[0]:ximax[0], ximin[1]:ximax[1], ximin[2]:ximax[2]]

        # flat array in the order of the first index variing fastest
        return ijk.ravel('F')

    def getIdsFromCoords(self, coords, bounds="clip"):
        """Returns the integer serial indexes for a list (an array) of points.
        Each point given by it's 3D coordinates is assigned to the nearest grid
        point.

        More precisely: The grid (self) is interpreted as consisting of cells
        of size self.spacing having the grid points as cell centres. Each of
        the given points is assigned to the cell it lies in, if any.

        @param coords: Nx3 array of point coordinates

        @param bounds: specifies how out-of-bounds points are handled. Can
            specify either one mode (e.g. "raise") or a tuple of modes, one
            mode per index (e.g. ("clip","raise","raise")).
                - raise - raise an error
                - clip - clip to the range
                - wrap - wrap around / use only if you know what you're doing

            In 'clip' mode, points outside the grid yield the index of the
            nearest boundary. A negative index which would normally wrap will
            clip to 0 instead. An index equal to or larger than the
            corresponding grid dimension (gridPtNb) will clip to the largest
            valid index of that axis.

            Default is 'clip' for all axes.

        @returns: list of serial indexes

        @note: To get grid of indexes for given point coords use this:
         >>> fIdx = np.round((coords - self.origin) / self.spacing).astype(int)
         >>> #... and now check for out-of-bounds conditions
        """
        coords = np.asarray(coords)

        # grid index triples
        fIdx = np.round((coords - self.origin) / self.spacing).astype(int)

        # ravel_multi_index expects a tuple of integer arrays, one array for
        # each dimension. Like this: (i-array, j-array, k-array)
        fIdx = fIdx.T
        # ravel to serial integer indexes
        # order="F" makes the first index i change fastest then j then k
        tIdx = np.ravel_multi_index(
            fIdx, self.gridPtNb, mode=bounds, order='F')

        return tIdx

    def getBoundingBox(self):
        """Returns a L{BoundingBox<bae.misc_01.BoundingBox>} with self.origin
        being one corner and the last point of the point grid being the
        opposite.
        """
        return np.array([
            self.origin,
            [x0 + (n-1)*dx for x0, n, dx in zip(
                self.origin, self.gridPtNb, self.spacing)]])

    @staticmethod
    def sharedIndicesFromAlignedGrids(topoA, topoB, returnSerialIndex=True,
                                returnBool=True, enforceAlignment=True):
        """Returns the grid indices and the boolean union of two overlapping
        structured grids.

        Example 1: get mean values of two overlapping fields
            >>> from bae.field_02 import FieldsCollection, StructuredGrid
            >>> fldColA = FieldsCollection().fromVtk('./myVtkBoxA.vtk')
            >>> fldColB = FieldsCollection().fromVtk('./myVtkBoxB.vtk')

            >>> iA,iB,boolTopo = StructuredGrid.sharedIndicesFromAlignedGrids(
            ...                             fldColA.topo, fldColB.topo)
            >>> if not boolTopo:
            >>>     raise ValueError('VTK-boxes do not overlap!')

            >>> fldColBool = FieldsCollection(topo=boolTopo)
            >>> fldColBool['UmeanAB'] = .5*(fldColA['U'][iA] + fldColB['U'][iB])
            >>> fldColBool.toVtk('./myVtkBoxBool_meanU.vtk')

        Example 2: insert masked data of setA into setB
            >>> iA,iB,_ = StructuredGrid.sharedIndicesFromAlignedGrids(
            >>>                             fldColA.topo, fldColB.topo)
            >>> nfld = fldColB['U'].copy()
            >>> # discard/mask all values above 0.2
            >>> nfld[nfld > .2] = np.nan
            >>> # paste field from boxB in copied U-field on boxA
            >>> fldColA['UJoined'] = fldColA['U'].copy()
            >>> fldColA['UJoined'][iA] = nfld[iB]
            >>> # get indices of masked pasted voxels
            >>> invMask = np.isnan(fldColA['UJoined'][0]) & iA
            >>> # (re)set to old values from boxA
            >>> fldColA['UJoined'][invMask] = fldColA['U'][invMask]

        @param topoA: L{StructuredGrid} object
        @param topoB: L{StructuredGrid} object with the same spacing as topoA.
        @param enforceAlignment: if True, topoA and topoB have to be perfectly
            aligned, what means that x,y and z-distace between the grid origins
            has (integer) multiple of gridspacing.
            If False, the origin of topoB gets 'rounded' onto the grid of
            topoA. So here the returned boolean topo is perfectly aligned to
            topoA, not topoB.
        @param returnSerialIndex: If True the serial indices and if False the
            grid-indices will be returned.
        @param return bool: If True, the returned indices will be boolean
            masks
        @returns: shared indices in topoA and topoB, L{StructuredGrid} that
            represents the
            the boolean union of both

        """
        if not np.isclose(topoA.spacing - topoB.spacing, 0).all():
            raise ValueError("Grid spacings are not the same.")

        spacing = np.asarray(topoA.spacing)
        gridPtNbA = np.asarray(topoA.gridPtNb).astype(int)
        gridPtNbB = np.asarray(topoB.gridPtNb).astype(int)
        originA = np.asarray(topoA.origin)
        originB = np.asarray(topoB.origin)
        deltaL = (originA - originB) / spacing

        if enforceAlignment and not np.isclose(deltaL % 1, 0).all():
            raise ValueError('Grid origins are not aligned.')

        deltaL = np.round(deltaL).astype(int)
        deltaU = deltaL + (gridPtNbA - gridPtNbB)
        idxLA, idxUA = np.zeros(3, dtype=int), gridPtNbA
        idxLB, idxUB = np.zeros(3, dtype=int), gridPtNbB
        for ii, (dL, dU) in enumerate(zip(deltaL, deltaU)):
            if dL < 0:
                idxLA[ii] = abs(dL)
            else:
                idxLB[ii] = abs(dL)
            if dU < 0:
                idxUB[ii] = gridPtNbB[ii] - abs(dU)
            else:
                idxUA[ii] = gridPtNbA[ii] - abs(dU)
            if ((idxUA < 0).any() or (idxUB < 0).any() or
                (idxLA > gridPtNbA).any() or (idxLB > gridPtNbB).any()):
                return [],[]

        xyzRangesA = [np.arange(l,u).astype(int) for l,u in zip(idxLA, idxUA)]
        xyzRangesB = [np.arange(l,u).astype(int) for l,u in zip(idxLB, idxUB)]
        xyzSliceA = [slice(l,u) for l,u in zip(idxLA, idxUA)]
        xyzSliceB = [slice(l,u) for l,u in zip(idxLB, idxUB)]

        if xyzRangesA[0].size:
            newOrigin = originA + idxLA*gridPtNbA
            newGridPtNb = [len(rr) for rr in xyzRangesA]
            booleanTopo = StructuredGrid(origin=newOrigin,
                                         spacing=spacing,
                                         gridPtNb=newGridPtNb)
        else:
            booleanTopo = None


        if not returnSerialIndex:
            if returnBool:
                for ii in range(3):
                    aa = np.zeros(topoA.gridPtNb[ii], dtype=bool)
                    bb = np.zeros(topoB.gridPtNb[ii], dtype=bool)
                    aa[xyzRangesA[ii]] = True
                    bb[xyzRangesB[ii]] = True
                    xyzRangesA[ii] = aa
                    xyzRangesB[ii] = bb
            return xyzRangesA, xyzRangesB, booleanTopo

        else:
            llA = np.zeros(topoA.gridPtNb,dtype=bool)
            llA[xyzSliceA[0], xyzSliceA[1], xyzSliceA[2]] = True
            idxA = np.ravel(llA, order='F')

            llB = np.zeros(topoB.gridPtNb,dtype=bool)
            llB[xyzSliceB[0], xyzSliceB[1], xyzSliceB[2]] = True
            idxB = np.ravel(llB, order='F')

            if not returnBool:
                idxA = np.argwhere(idxA).ravel()
                idxB = np.argwhere(idxB).ravel()

        return idxA, idxB, booleanTopo



class StructuredGridSparse(StructuredGrid):
    """Represents a rectilinear grid which has regularly uniformly spaced
    points of which only a subset is actually active.
    The grid is defined by grid dimensions gridPtNb=(nx, ny, nz), grid
    origin=(x0, y0, z0) and point spacing=(dx, dy, dz). Thus all coordinates
    can be derived from the serial integer index (see below).

    Nomenclature:
     - grid index triple: a tuple of three zero-based indexes identifying
       the position in the grid in x-, y-, z-direction.
     - serial integer index of the full grid, in short "full serial index":
       plain integer number starting at zero, max possible value: nx*ny*nz-1.
       Not all index values in this range correspond to an active point.
       Therefore not all indexes in the range are valid.
     - sparse serial integer index, in short "sparse serial index":
       plain integer number starting at zero. This is the
       numpy index to access a particular value in a corresponding field.

    The values of related fields are stored sparsely: Only points with
    (valid) data will be stored and their postion is defined by
    the sparse serial index integer index.

    A StructuredGridSparse-object can be constructed
     - from another StructuredGridSparse or StructuredGrid object (copy
       constructor)
     - by specifying the grid characteristics origin, gridPtNb, spacing
       as keyword arguments ..... and more....
     - ###GP: by specifying a sparse (or complete) index as a last argument
       or as a keywordargument
    """

    def __init__(self, *args, **kwargs):
        """
        A StructuredGrid-object can be constructed
            - as an empty object
            - from another StructuredGrid object (copy constructor)
            - by specifying the grid characteristics origin, gridPtNb, spacing
              as keyword arguments

            - by specifying a sparse (or complete) index as a last argument
              or as a keywordargument

        Example:
        ========
            >>> gridPtNb = (10,11,12)
            >>> origin = (0.0, 100.0, 2000.0)
            >>> spacing = (101, 102, 103)
            >>> knownValIdx = [120, 177,100,1001]
            >>> sparseGrid = StructuredGrid(knownValIdx,
            ...                             gridPtNb, origin, spacing)
        """
        raise NotImplementedError


#@@@@@@ CONSTRUCTION SITE AHEAD @@@@@@
#
#    def _checkConsistency(self, allowEmpty=True):
#        """Tests proper grid definition
#        """
#        test = [self.gridPtNb, self.spacing, self.origin]
#
#        err = False
#        if all(a is None for a in test):  # not all unset
#            pass
#        else:
#            if any(len(a) != 3 for a in test):  # all 3-Tuples
#                err = True
#
#            if not all([aa > 0 for a in test[:2] for aa in a]):
#                err = True  # nbs,spaces >0
#
#            if not all([isinstance(b, int) for b in test[0]]):  # nbs are int
#                err = True
#
#        if err:
#            msg("Some of the grid parameters are not defined properly.\n"
#                "gridPtNb: " + str(self.gridPtNb) + " ? (uint,uint,uint)\n"
#                "spacing: " + str(self.spacing) + " ? (f+,f+,f+)\n"
#                "origin:" + str(self.origin) + " ? (f,f,f)")
#            raise ValueError("Grid not properly defined.")
#
#        # check maxIndex of idx
#        if hasattr(self, 'idx'):
#            if self.idx.max() + 1 > self.nFullPts():
#                raise ValueError('Some integer indexes (max=%d) are larger'
#                                 ' than the total number of grid points (%d).'
#                                 % (self.idx.max(), self.nFullPts()))
#
#        # check length of all stored staticFields
#        StaticEntities._checkConsistency(self)
#
#
#
#    ###GP: rename flat -> grid index triple
#    def idxToFlatIdx(self):
#        """Unravels integer indexes to a sequence of grid index triples
#
#        @note: requires the attributes 'gridPtNb', 'origin' and
#        'spacing' to exist
#        """
#        return np.unravel_index(self.idx, self.gridPtNb, order='F')
#
#
#    def flatIdxToCoordComp(self, idx, component):
#        """Returns the coordinate component from grid index triple.
#
#        @param idx: integer index/-es of the grid point(s)
#        @param component: 0 for x-coordinate, 1 for y and 2 for z
#
#        @note: requires the attributes 'gridPtNb', 'origin' and
#        'spacing' to exist
#        """
#        return idx * self.spacing[component] + self.origin[component]
#
#
#    def idxToCoord(self, idx):
#        """Returns the coordinates from integer indexes
#        """
#        coord = np.vstack([self.flatIdxToCoordComp(cIdx, ii)
#                           for ii, cIdx in enumerate(self.idxToFlatIdx())]).T
#        return coord
#
#
#    def coordToIdx(self, coords, bounds=('clip', 'wrap')):
#        """Returns the integer index from a list of coordinates with. The
#        bounds kwarg specifies behavior for coordinates which are not inside
#        the grid.
#
#        @param coords: list of coordinates
#
#        @param bounds: specifies how out-of-bounds indices are handled. Can
#            specify either one mode or a tuple of modes, one mode per index.
#                - raise - raise an error
#                - wrap - wrap around
#                - clip - clip to the range
#
#            In 'clip' mode, a negative index which would normally wrap will
#            clip to 0 instead.
#            Default is ('clip', 'wrap'), so outside coordinates will get
#            indexed at the 'nearest' bound.
#
#        @returns: flattened index list in fortran order
#
#        @note: requires the attributes 'gridPtNb', 'origin' and
#        'spacing' to exist
#        """
#        coords = np.asarray(coords)
#        origin = np.asarray(self.origin)
#        gridPtNb =np.asarray(self.gridPtNb)
#        spacing = np.asarray(self.spacing)
#
#        # grid index triples
#        fIdx = np.round((coords - origin) / spacing).astype(int)
#
#        # ravel to integer indexes
#        tIdx = np.ravel_multi_index(fIdx, tuple(gridPtNb),
#                                    mode=bounds, order='F')
#
#        return tIdx
#
#
#    def coords(self):
#        """Returns coordinates of stored grid points
#
#        @note: requires the attributes 'idx', 'gridPtNb', 'origin' and
#        'spacing' to exist
#        """
#        return self.idxToCoord(self.idx)
#
#    def flatX(self):
#        """Returns (full) list of x-coordinates of dimension
#        self.gridPtNb[0]
#
#        @note: requires the attributes 'gridPtNb', 'origin' and
#        'spacing' to exist
#        """
#        return self.flatIdxToCoordComp(np.arange(self.gridPtNb[0]), 0)
#
#    def flatY(self):
#        """Returns (full) list of y-coordinates of dimension
#        self.gridPtNb[1]
#
#        @note: requires the attributes 'gridPtNb', 'origin' and
#        'spacing' to exist
#        """
#        return self.flatIdxToCoordComp(np.arange(self.gridPtNb[1]), 1)
#
#
#    def flatZ(self):
#        """Returns (full) list of z-coordinates of dimension
#        self.gridPtNb[3]
#
#        @note: requires the attributes 'gridPtNb', 'origin' and
#        'spacing' to exist
#        """
#        return self.flatIdxToCoordComp(np.arange(self.gridPtNb[2]), 2)
#
#
#    def fullXYZGrid(self):
#        """Returns the numpy.meshgrid fields of coordinates
#
#        @returns: (X,Y,Z), each of dimension self.gridPtNb
#
#        @note: requires the attributes 'gridPtNb', 'origin' and
#        'spacing' to exist
#        """
#        return np.meshgrid(self.flatX(), self.flatY(), self.flatZ(),
#                           indexing='ij')
#
#
#    def fullValuesGrid(self, values, emptyValues=np.nan, rollPositionAxis=True):
#        """Arranges field values (flat) on given topo to full 3D X-Y-Z field
#        as created from numpy.meshgrid.
#
#        Some usefull scipy and numpy fuctions assume the last 3 array-dimension
#        to represent the coordinates, while our fields framework referes to the
#        first dimension. Thus, by default, the positional axis will be rolled
#        to last dimensions. Set rollPositionAxis=False to avoid this behaviour.
#
#        Example:
#        ========
#         >>> print 'type of myField: ', myField
#         type of myField: <type 'ScalarField'>
#         >>> print 'sparse field values shape: ', myField.values.shape
#         sparse field values shape: (100,3)
#         >>> print 'grid size: ', myField.topo.gridPtNb
#         grid size: (101,102,103)
#         >>> npArrA = myField.topo.fullValuesGrid(myField)
#         >>> npArrB = myField.topo.fullValuesGrid(
#         >>>              myField, rollPositionAxis=False)
#         >>> print 'npArrA.shape=',npArrA.shape
#         npArrA.shape=(3,101,102,103)
#         >>> print 'npArrB.shape=',npArrB.shape
#         npArrA.shape=(101,102,103,3)
#
#        @param values: can be a L{Scalar<bae.fieldData_00.fields.ScalarField>}-,
#            L{Vector<bae.fieldData_00.fields.VectorField>}- or
#            L{Tensor<bae.fieldData_00.fields.TensorField>} or a numpy-array of
#            shape (nPts, nFrames, :, :), where nPts has to be len(self).
#
#        @param emptyValues: value to fill undefined points
#
#        @param rollPositionAxis: if true the positional axes will get rolled to
#            last array dimension (as expected by many scipy/numpy functions)
#
#        @returns: numpy.meshgrid-like array
#
#        @note: requires the attributes gridPtNb', 'origin' and 'spacing'
#        to exist
#        """
#        if hasattr(values, 'topo') and hasattr(values, 'values'):
#            if (values.topo is not None) and not (self == values.topo):
#                raise ValueError('The topo-object of the fields-object you'
#                                 ' passed differ.')
#            else:
#                values = values.values
#
#
#        values = np.asarray(values)
#        if not values.shape[0] == len(self):
#            raise ValueError(
#                "The first dimension of the values array (shape %s) doesn't"
#                " match the length of topo-object %d."
#                % (values.shape, len(self)))
#
#        # create empty field
#        out = emptyValues * np.empty((self.nFullPts(),) + values.shape[1:])
#
#        # set values to ozt-field at stored pointIdx
#        out[self.idx] = values
#
#        if rollPositionAxis:
#            #roll pointIndex to last dimension
#            out = np.rollaxis(out, -1)
#            out = out.reshape(out.shape[:-1] + tuple(self.gridPtNb),
#                              order='F')
#
#        else:
#            out = out.reshape(tuple(self.gridPtNb) + out.shape[1:],
#                              order='F')
#
#        return out

#} end of Topology Classes






if __name__ == '__main__':
    print("No Syntax Errors")
    #test_commenIndicesFromAlignedGrids()
