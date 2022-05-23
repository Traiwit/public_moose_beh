#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""This module contains functions to interpolate field data to other topo's.
They are supplied to the main L{FieldsCollection<bae.field_02.FieldsCollection>}
class through the Mixin-class L{FieldsCollectionInterpolation} defined in this
module.


# lowlevel
for targetTopo in targetTopos:
    interpol = mySrcPtCloud.interpolate(targetPtCloud, method='nearst',
                                        cutoff=5.)
    myIpFlds['U'] = interpol(srcFlds['U'])

# highlevel
myIpFlds = myFlds.interpolate(myPtCloud, method='nearst', cutoff=5.)


"""

# Note: version is stored in bae/field_02/__init__.py

todo = """
- for interpolation from unstructured meshes to arbitrary points finding the
  position in the mesh (which element and elemCoords) is expensive. See
  mesh_01.Mesh.getElemCoordsForGrid and mesh_01.InterpolMeshToPoints.

  Suggestion: Again have an InterpolMeshToPoints class that holds the data.
  Provide a FieldsCollectionInterpolation - method (or HomoMesh/HeteroMesh
  -class method) that creates such an interpolation. See mesh_01 and
  suggestions at the top of the module for implementation details. Then
  optionally accept an InterpolMeshToPoints objects as input to
  FieldsCollectionInterpolation.interpolate method instead of otherTopo.
- unittests
"""
import numpy as np
try:
    from scipy.interpolate import (NearestNDInterpolator,
                                   LinearNDInterpolator,)
    from scipy.spatial import Delaunay, qhull
    from scipy.interpolate.interpnd import _ndim_coords_from_arrays
except:
    NearestNDInterpolator = None
    LinearNDInterpolator = None


from itertools import izip, product




#from .topo.unstructuredmesh.meshbase import \
#    HomoMesh, MeshElems
#
#from .topo.structuredgrid import StructuredGrid
#from .topo.pointcloud import PointCloud

from topo.unstructuredmesh.meshbase import \
    HomoMesh, MeshElems

from topo.structuredgrid import StructuredGrid
from topo.pointcloud import PointCloud


from bae.log_01 import msg

class FieldsCollectionInterpolation(object):
    """Mixin for FieldsCollection supplying methods to interpolate field data
    to other topo's.
    """

    def interpolate(self, otherTopo):
        """Interpolate self to the specified other topo.

        @param otherTopo: an object of a subclass of L{TopoBase}.
        """

        # special cases: nodal field to elem centroid (topo-type: MeshElems)
        if (otherTopo.position == "element" and self.topo.position == "node"
            and otherTopo.mesh is self.topo.mesh):

            mesh = self.topo.mesh  # abbreviation
            newFlds = type(self)(topo=MeshElems(mesh), times=self.times)

            if isinstance(mesh, HomoMesh):
                # process all fields
                for fieldName, field in self.fields.iteritems():
                    # average nodal values for each element
                    newFlds[fieldName] = np.average(
                        # take value for each node...
                        np.take(field, mesh.elNodes, axis=0),
                        # and average over
                        axis=1)
                # process self.space
                for fieldName, field in self.space.iteritems():
                    # average nodal values for each element
                    newFlds.space[fieldName] = np.average(
                        # take value for each node...
                        np.take(field, mesh.elNodes, axis=0),
                        # and average over
                        axis=1)
            else:
                raise NotImplementedError(
                    "Interpolation to element centroids of HeteroMeshes not"
                    " implemented yet.")

            # done...
            return newFlds

        else:
            raise NotImplementedError(
                "Normal interpolation not implemented yet.")


#%% fieldscollection interpolator classes
###############################################################################

#%% from topo to PointCloud
class ToPointsInterpolator(object):
    fromPointInterpolationMethods = ['closest', 'linear', 'radial']

    def __init__(self, topo, otherTopo, **kwargs):
        if isinstance(topo, (list, np.ndArray)):
            topo = PointCloud(coords=topo)

        _testTopoType(topo, PointCloud, 'Target topology')

        if isinstance(otherTopo, StructuredGrid):
            return self.fromStructuredGrid(topo, otherTopo, **kwargs)

        if isinstance(otherTopo, (np.ndArray, PointCloud)):
            return self.fromPoints(topo, otherTopo, **kwargs)

        raise NotImplementedError("Can't interpolate to a topology of type %s"
                                  % type(otherTopo))

    @classmethod
    def fromPoints(cls, topo, otherTopo, **kwargs):
        _testTopoType(topo, (PointCloud, np.array, list), 'Target topology')
        _testTopoType(otherTopo, (PointCloud, np.array, list), 'Source topology')

        if isinstance(topo, PointCloud):
            topo = topo.getPoints()
        if isinstance(otherTopo, PointCloud):
            otherTopo = otherTopo.getPoints()

        topo = np.atleast_2d(topo)
        otherTopo = np.atleast_2d(otherTopo)

        if not topo.shape[1] == otherTopo.shape[1]:
            raise ValueError("PointClouds do not have the same dimensions: "
                             "target: %s, source %s"
                             % (topo.shape, otherTopo.shape))


    @classmethod
    def fromStructuredGrid(cls, topo, otherTopo, **kwargs):
        # testing inupts
        _testTopoType(topo, PointCloud, 'Target topology')
        _testTopoType(otherTopo, StructuredGrid, 'Source topology')

        method = kwargs.get('method', 'linear')
        if method not in cls.fromPointInterpolationMethods:
            raise ValueError("Don't know interpolation method %s for"
                             " StructuredGrid-to-StructuredGrid interpolation"
                              % method)


#%% from topo to structured grid
class ToStructuredGridInterpolator(object):
    fromPointInterpolationMethods = ['closest', 'linear', 'radial']

    def __init__(self, topo, otherTopo, **kwargs):
        if not isinstance(topo, StructuredGrid):
            raise ValueError("Target topology has to be of type StructuredGrid"
                             " You passed %s" + type(topo))

        if isinstance(otherTopo, StructuredGrid):
            return self.fromStructuredGrid(topo, otherTopo, **kwargs)

        if isinstance(otherTopo, (np.ndArray, PointCloud)):
            return self.fromPoints(topo, otherTopo, **kwargs)

        raise NotImplementedError("Can not interpolate to a topology of type %s"
                             % type(otherTopo))


    @classmethod
    def fromStructuredGrid(cls, topo, otherTopo, **kwargs):
        _testTopoType(topo, StructuredGrid, 'Target topology')
        _testTopoType(otherTopo, StructuredGrid, 'Source topology')

        method = kwargs.get('method', 'linear')
        if method not in cls.fromPointInterpolationMethods:
            raise ValueError("Don't know interpolation method %s for"
                             " StructuredGrid-to-StructuredGrid interpolation"
                              % method)


#%% Base interpolator classes - mainly derived from scipy.interpolation

#%% PointCloud as source
class LinearNDInterpolator(object):
    """
    Linear interpolation from an abitrary set of points to an other
    abitrary set of points.
    In a first step a Delaunay triangulation of the source points will be
    created.

    See U{scipy.interpolation<https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.LinearNDInterpolator.html>}
    for original implementation.


    """
    def __init__(self, points):
        if len(points) > 1E5:
            msg('processing DELAUNAY triangulation')
        self.tri = Delaunay(points)

        self._weights = None
        self._simplexIdx = None
        self._ptsHash = None

    def __call__(self, needlePoints, values):
        """
        Interpolation at coordinates
        @param needlePoints: ndarray of shape (..., ndim)
            The coordinates to sample the gridded data at
        @param values : array_like, shape (m1, ..., mn, ...)
            The source data/field on the regular grid in n dimensions that will
            be sampled
        @note:
            If the target points (xi) changed during the last call, the indices
            and weights need to be updated. So best practice is to alter
            source-values (and frames) first before you change the target
            points
        """

        values = np.asarray(values)
        try:
            needlePoints = needlePoints.getPoints()
        except AttributeError:
            pass
        needlePoints = np.asarray(needlePoints)
        k,n = needlePoints.shape

        ptsHash = hash(needlePoints.tostring())
        if self._weights is None or not (self._ptsHash==ptsHash):
            msg("finding needle points in triangulation")
            # find simplexes that contain interpolated points
            s = self.tri.find_simplex(needlePoints)
            validIdx = (s >= 0)
            # get transform matrices for each simplex (see explanation bellow)
            m = self.tri.transform[s[validIdx]]

            # for each interpolated point p, mutliply the transform matrix by
            # vector p-r, where r=m[:,n,:] is one of the simplex vertices to which
            # the matrix m is related to (again, see bellow)
            b = np.einsum('ijk,ik->ij', m[:,:n,:n],
                          needlePoints[validIdx]-m[:,n,:])

            # get the weights for the vertices; `b` contains an n-dimensional vector
            # with weights for all but the last vertices of the simplex
            # (note that for n-D grid, each simplex consists of n+1 vertices);
            # the remaining weight for the last vertex can be copmuted from
            # the condition that sum of weights must be equal to 1
            w = np.c_[b, 1-b.sum(axis=1)]

            weights = np.full((k,n+1), np.nan)
            weights[validIdx] = w

            self._weights = weights
            self._simplexIdx = s
            self._ptsHash = ptsHash

        v = self.tri.vertices[self._simplexIdx]
        extraDim =  tuple([1 for _ in range(values.ndim-1)])
        w = self._weights.reshape((-1,n+1,) + extraDim)
        return (w*values[v,...]).sum(axis=1)


class NearestNDInterpolator(NearestNDInterpolator):
    """
    NearestNeighbor-"Interpolation" from an abitrary set of points to an other
    abitrary set of points.

    See U{scipy.interpolation<https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.NearestNDInterpolator.html>}
    for original implementation.



    """
    def __init__(self, points, **kwargs):
        # cached results
        self._xiHash = None
        self._dist, self._idx = None, None

        dummy = np.empty(len(points))
        super(NearestNDInterpolator,self).__init__(points, dummy, **kwargs)

    def __call__(self, needlePoints, values, cutoffDist=None,
                 fill_value=np.nan):
        """
        @param needlePoints : ndarray of shape (..., ndim)
            The coordinates to sample the gridded data at
        @param values : array_like, shape (m1, ..., mn, ...)
            The source data/field on the regular grid in n dimensions that will
            be sampled
        @param: cutoffDist : number, optional
            If provided, only needlePoints with a distance equal-closer than
            cutoffDist will get values assigned. The fill_value will be
            assigned to all points outside cutOffdist.
        @param: fill_value : number, optional
            If provided, the value to use for points outside of cutoffDist.
            Default is nan.
        @note:
            If the target points (xi) changed during the last call, the indices
            and weights need to be updated. So best practice is to alter
            source-values (and frames) first before you change the target
            points
        """
        values = np.asarray(values)
        xi = _ndim_coords_from_arrays(needlePoints, ndim=self.points.shape[1])

        # test if output-topo is new
        xiHash = hash(xi.tostring())
        if not (xiHash == self._xiHash):
            self._xiHash = xiHash
            xi = self._check_call_shape(xi)
            xi = self._scale_x(xi)
            self._dist, self._idx = self.tree.query(xi)

        if cutoffDist is None:
            return values[self._idx]
        else:
            outShape = values.shape
            A = np.fill(outShape, fill_value)
            mask = (self._dist <= cutoffDist)
            A[mask] = values[self._idx][mask]
            return A




class RegularGridInterpolator(object):
    """
    Interpolation from a regular grid to abitrary points
    The data must be defined on a regular grid; the grid spacing however may be
    uneven. Linear and nearest-neighbor interpolation are supported. After
    setting up the interpolator object, the interpolation method (*linear* or
    *nearest*) may be chosen at each evaluation.

    See U{scipy.interpolation<https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html>}
    for original implementation.

    @param points : tuple of ndarray of float, with shapes (m1, ), ..., (mn, )
        The points defining the regular grid in n dimensions or
        L{PointCloud}object. Note that spacing can be irregural here. Thus,
        x,y and z vectors instead of origin, gridPtNb and spacing are required.
    @param method : str, optional
        The method of interpolation to perform. Supported are "linear" and
        "nearest". This parameter will become the default for the object's
        L{__call__} method. Default is "linear".
    @param bounds_error : bool, optional
        If True, when interpolated values are requested outside of the
        domain of the input data, a ValueError is raised.
        If False, then `fill_value` is used.
    @param: fill_value : number, optional
        If provided, the value to use for points outside of the
        interpolation domain. If None, values outside
        the domain are extrapolated.

    @notes:
    Contrary to LinearNDInterpolator and NearestNDInterpolator, this class
    avoids expensive triangulation of the input data by taking advantage of the
    regular grid structure.
    If any of `points` have a dimension of size 1, linear interpolation will
    return an array of nan values. Nearest-neighbor interpolation will work
    as usual in this case.

    Examples
    --------
    Evaluate a simple example function on the points of a 3-D grid:
     >>> from scipy.interpolate import RegularGridInterpolator
     >>> def f(x, y, z):
     ...     return 2 * x**3 + 3 * y**2 - z
     >>> x = np.linspace(1, 4, 11)
     >>> y = np.linspace(4, 7, 22)
     >>> z = np.linspace(7, 9, 33)
     >>> data = f(*np.meshgrid(x, y, z, indexing='ij', sparse=True))
     >>> myInterpolator = RegularGridInterpolator((x, y, z))
     >>> pts = np.array([[2.1, 6.2, 8.3], [3.3, 5.2, 7.1]])
     >>> myInterpolator(pts, data)
     array([ 125.80469388,  146.30069388])
    """

    def __init__(self, points, method="linear", bounds_error=True,
                 fill_value=np.nan):
        if method not in ["linear", "nearest"]:
            raise ValueError("Method '%s' is not defined" % method)
        self.method = method
        self.bounds_error = bounds_error

        if isinstance(points, PointCloud):
            os, sps, nPts = points.origin, points.spacing, points.gridPtNb
            points = tuple([o + sp*np.arange(nPt)
                            for o, sp, nPt in zip(os,sps,nPts)])

        self.gridPtNb = tuple([len(p) for p in points])

        for i, p in enumerate(points):
            if not np.all(np.diff(p) > 0.):
                raise ValueError("The points in dimension %d must be strictly "
                                 "ascending" % i)
            if not np.asarray(p).ndim == 1:
                raise ValueError("The points in dimension %d must be "
                                 "1-dimensional" % i)

        self.grid = tuple([np.asarray(p) for p in points])

        self.fill_value = fill_value

        self._xiHash = None
        self._idx_res = None
        self._weights = None
        self._valuesDim = None
        self._valuesGridShaped = True


    def _testValueInput(self, values):
        if not hasattr(values, 'ndim'):
            # allow reasonable duck-typed values
            values = np.asarray(values)

        if (not (values.shape[:3] == self.gridPtNb)
            and values.shape[0] == np.prod(self.gridPtNb)):
            self._valuesGridShaped = False
            values = values.reshape(self.gridPtNb + values.shape[1:],
                                    order = 'F')
        else:
            self._valuesGridShaped = True

        nPts = len(self.grid)
        if nPts > values.ndim:
            raise ValueError("There are %d point arrays, but values has %d "
                             "dimensions" % (nPts, values.ndim))


        for i,nPt in enumerate(self.gridPtNb):
            if not values.shape[i] == len(nPt):
                raise ValueError("There are %d points and %d values in "
                                 "dimension %d"
                                 % (len(nPt), values.shape[i], i))

        if hasattr(values, 'dtype') and hasattr(values, 'astype'):
            if not np.issubdtype(values.dtype, np.inexact):
                values = values.astype(float)

        fill_value = self.fill_value
        if fill_value is not None:
            fill_value_dtype = np.asarray(fill_value).dtype
            if (hasattr(values, 'dtype') and not
                    np.can_cast(fill_value_dtype, values.dtype,
                                casting='same_kind')):
                raise ValueError("fill_value must be either 'None' or "
                                 "of a type compatible with values")
        return values

    def __call__(self, needlePoints, values, method=None):
        """
        Interpolation at coordinates
        @param needlePoints: ndarray of shape (..., ndim)
            The coordinates to sample the gridded data at
        @param values : array_like, shape (m1, ..., mn, ...)
            The source data/field on the regular grid in n dimensions that will
            be sampled
        @param method : str
            The method of interpolation to perform. Supported are "linear" and
            "nearest".
        @note:
            If the target points (xi) changed during the last call, the indices
            and weights need to be updated. So best practice is to alter
            source-values (and frames) first before you change the target
            points
        """
        xi = needlePoints #abrav.

        values = self._testValueInput(values)

        method = self.method if method is None else method
        if method not in ["linear", "nearest"]:
            raise ValueError("Method '%s' is not defined" % method)

        ndim = len(self.grid)
        xi = _ndim_coords_from_arrays(xi, ndim=ndim)
        if xi.shape[-1] != len(self.grid):
            raise ValueError("The requested sample points xi have dimension "
                             "%d, but this RegularGridInterpolator has "
                             "dimension %d" % (xi.shape[1], ndim))

        xi_shape = xi.shape
        xi = xi.reshape(-1, xi_shape[-1])

        if self.bounds_error:
            for i, p in enumerate(xi.T):
                if not np.logical_and(np.all(self.grid[i][0] <= p),
                                      np.all(p <= self.grid[i][-1])):
                    raise ValueError("One of the requested needlePoints"
                                     " is out of bounds in dimension %d" % i)

        # check if target topology has changed
        xiHash = hash(xi.tostring())
        if not (xiHash == self._xiHash):
            ## has changed -> recalculate indices and weights
            indices, norm_distances, out_of_bounds = self._find_indices(xi.T)
            self._idx = indices, norm_distances, out_of_bounds
            # reset stored data
            self._idx_res, self._weights, self._valuesDim = None, None, None
            self._xiHash = xiHash
        else:
            ## no change -> reuse old indices and weights
            indices, norm_distances, out_of_bounds = self._idx

        if method == "linear":
            result = self._evaluate_linear(indices,
                                           norm_distances,
                                           out_of_bounds,
                                           values)
        elif method == "nearest":
            result = self._evaluate_nearest(indices,
                                            norm_distances,
                                            out_of_bounds,
                                            values)
        if not self.bounds_error and self.fill_value is not None:
            result[out_of_bounds] = self.fill_value

        if self._valuesGridShaped:
            return result.reshape(xi_shape[:-1] + values.shape[ndim:])
        else:
            raise NotImplementedError('Correct result shape not tested yet')
            return result

    def _evaluate_linear(self, indices, norm_distances, out_of_bounds, values):
        # slice for broadcasting over trailing dimensions in self.values
        vslice = (slice(None),) + (None,)*(values.ndim - len(indices))

        # find relevant values
        # each i and i+1 represents a edge
        edges = product(*[[i, i + 1] for i in indices])

        vals = 0.
        if self._weights is None or not self._valuesDim == values.ndim:
            weights = []
            for edge_indices in edges:
                weight = 1.
                for ei, i, yi in zip(edge_indices, indices, norm_distances):
                    weight *= np.where(ei == i, 1 - yi, yi)
                weights.append(weight[vslice])
            self._weights = weights
            self._valuesDim = values.ndim

        vals = 0.
        edges = product(*[[i, i + 1] for i in indices])
        for edge_indices, weight in izip(edges, self._weights):
            vals += np.asarray(values[edge_indices]) * weight
        return vals

    def _evaluate_nearest(self, indices, norm_distances, out_of_bounds, values):
        if self._idx_res is None:
            self._idx_res = [np.where(yi <= .5, i, i + 1)
                       for i, yi in zip(indices, norm_distances)]
        return values[tuple(self._idx_res)]

    def _find_indices(self, xi):
        # find relevant edges between which xi are situated
        indices = []
        # compute distance to lower edge in unity units
        norm_distances = []
        # check for out of bounds xi
        out_of_bounds = np.zeros((xi.shape[1]), dtype=bool)
        # iterate through dimensions
        for x, grid in zip(xi, self.grid):
            i = np.searchsorted(grid, x) - 1
            i[i < 0] = 0
            i[i > grid.size - 2] = grid.size - 2
            indices.append(i)
            norm_distances.append((x - grid[i]) /
                                  (grid[i + 1] - grid[i]))
            if not self.bounds_error:
                out_of_bounds += x < grid[0]
                out_of_bounds += x > grid[-1]
        return indices, norm_distances, out_of_bounds


###############################################################################
def _createSomeTestData(nx, ny, nz, asGrid=True):
    def f(x, y, z):
        return 2 * x**3 + 3 * y**2 - z

    x = np.linspace(0, 100, nx)
    y = np.linspace(0, 100, ny)
    z = np.linspace(0, 100, nz)

    data = f(*np.meshgrid(x, y, z, indexing='ij', sparse=True))
    if asGrid:
        xyz = (x,y,z)
    else:
        data = data.ravel(order='F')
        xyz = np.array([cc.ravel(order='F')
                        for cc in np.meshgrid(x, y, z, indexing='ij')]).T
    return xyz, data

def _testRegularGridInterpolator():
    from scipy.interpolate import RegularGridInterpolator as ScipyIp
    nx, ny, nz = 200, 250, 300
    (x,y,z), data = _createSomeTestData(nx, ny, nz, asGrid=True)
    scipyIpFun = ScipyIp((x, y, z), data)
    ipFun = RegularGridInterpolator((x, y, z))

    nTarget = int(1E6)
    pts = np.random.rand(nTarget, 3) * 100
    msg('%d sourcePoints vs %d targetPoints' % (nx*ny*nz, nTarget))
    msg('pure scipy')
    rscipy = scipyIpFun(pts)
    msg('first own')
    ipFun(pts, data)
    msg('second own')
    rOwn = ipFun(pts, data)
    msg('same result: %s' % (rscipy == rOwn).all())


def testNearestNdInterpolator():
    from scipy.interpolate import NearestNDInterpolator as ScipyIp
    nx, ny, nz = 150, 50, 120
    xyz, data = _createSomeTestData(nx, ny, nz, asGrid=False)

    scipyIpFun = ScipyIp(xyz, data)
    ipFun = NearestNDInterpolator(xyz)

    nTarget = int(1E6)
    pts = np.random.rand(nTarget, 3) * 100
    msg('%d sourcePoints vs %d targetPoints' % (nx*ny*nz, nTarget))
    msg('pure scipy')
    rscipy = scipyIpFun(pts)
    msg('first own')
    ipFun(pts, data)
    msg('second own')
    rOwn = ipFun(pts, data)
    msg('same result: %s' % (rscipy == rOwn).all())

def _testTopoType(topo, targetTypes, descText="Topology"):
    if isinstance(topo, targetTypes):
        return # all good

    #raise Error
    try:
        targetTypes[0]
        isIter = True
    except TypeError:
        isIter = False

    if not isIter or len(targetTypes) == 1:
        try:
            targetTypes = targetTypes[0]
        except:
            pass
        errTxt = (descText + ' has to be of type %s. You passed %s.'
                  % ( targetTypes.__name__, type(topo)))
    else:
        targetNames = [t.__name__ for t in targetTypes]
        tt = ', '.join(targetNames[:-1]) + ' or %s' % targetNames[-1]
        errTxt = (descText + ' has to be of type %s. You passed %s.'
                  % (tt, type(topo)))
    raise ValueError(errTxt)


def _test_LinearNDInterpolator():
    points = [[0,0,0], [1,0,0], [0,1,0], [0,0,1]]
    needlePoints = points + [[0.25, 0.25, 0.25], [1,1,1]]
    values = [1,2,3,4,]
    lIp = LinearNDInterpolator(points)
    print lIp(needlePoints, values)

    nn, mm = int(1E3), int(1E3)
    points = np.random.rand(nn, 3)

    needlePoints = np.random.rand(mm, 3)
    values = np.random.rand(nn,5,3)
    msg('process triangulation')
    lIp = LinearNDInterpolator(points)
    lIp(needlePoints, values)
    lIp(needlePoints, values)
    msg('done')


if __name__ == '__main__':
    print('No syntax errors')
    #testRegularGridInterpolator()
    #testNearestNdInterpolator()
    test_LinearNDInterpolator()