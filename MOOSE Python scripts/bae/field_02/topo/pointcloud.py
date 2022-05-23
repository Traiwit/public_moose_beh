# -*- coding: utf-8 -*-
"""Topology classes for unstructured points, i.e. point clouds.

"""

from itertools import chain
import numpy as np

from bae.misc_01 import BoundingBox

from bae.field_02.topo.base import TopoBase
# from bae.log_01 import msg


###############################################################################
#{Topology Classes

class PointCloud(TopoBase):
    """Class to manage unstructured points. Points are
    --at least-- specified by their coordinates (coords).
    The coordinates have to be unique within a tolerance tol.

    @ivar coords: Nx3 array of point coords. Don't access directly, use
       self.L{getPoints}()
    @cvar position: "point". For L{FieldsCollection}.topo objects the
      position attribute identifies the position of values ("node", "element",
      "elemIP", "point" or "structPt").

    @Note: It's not feasible to modify the coordinates of a PointCloud object.
    """

    position = "point"

    def __init__(self, coords=None):
        """A PointCloud object can be created as an empty PointCloud (len will
        be zero) or by passing the coordinates or another topo object.

        Example:
        ========
         >>> nPts = 5
         >>> # some points
         >>> x = np.arange(nPts)
         >>> y = np.arange(nPts) + 10
         >>> z = np.arange(nPts) + 100
         >>> coords = np.vstack((x, y, z)).T
         >>> points = PointCloud(coords)
         >>> # index points
         >>> subPoints = points[1:3]
         >>> print subPoints.getPoints()

        @param coords: Nx3 array of point coords
        """

        if isinstance(coords, TopoBase):
            self.coords = coords.getPoints()
        else:
            self.coords = coords

        # for storring KDTree and hash (to check if update req) of coords
        self._KDTree = None

    def getSubTopo(self, index):
        """Service method for FieldsCollection.__getitem__() selecting from the
        space dimension.

        Returns a L{PointCloud} of only the specified points.

        @param index: an index or indexes --anything suitable to have numpy
           choose from one dimension.
        """
        return PointCloud(self.coords[index])

    def getPoints(self, index=None):
        """Returns a copy of the point coordinates.

        It's a copy to make sure that the stored coords don't change!

        @param index: If not given (None) then return the coordinates of *all*
           points. Otherwise an index or indexes --anything suitable to have
           numpy choose from one dimension.
        """
        if index is None:
            return np.array(self.coords)
        else:
            return self.coords[index]

    def __len__(self):
        return len(self.coords)

    def __getitem__(self, index):
        """Return a PointCloud of only the specified points.

        @param index: An index suitable to have numpy choose from one
           dimension.
        """
        return PointCloud(self.coords[index])

    def KDTree(self):
        """Returns a KDTree (see U{https://docs.scipy.org/doc/scipy/reference/
        generated/scipy.spatial.cKDTree.html#scipy.spatial.cKDTree}) of the
        point coords. The KDTree is stored internally to save its costly
        creation.
        """
        if self._KDTree is None:
            from scipy.spatial import cKDTree
            self._KDTree = cKDTree(self.coords)

        return self._KDTree

    ###GP: needs a more realistic example in the doc string
    def getIdsFromPoints(
            self, points, atol=1E-6, shape='sphere', strict=True):
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
        from scipy.spatial import cKDTree

        if shape == 'sphere':
            pNorm = 2.
        elif shape == 'box':
            pNorm = np.inf
        else:
            raise ValueError("Don't know search shape '%s'" % str(shape))

        coordsTree = self.KDTree()
        rstTree = cKDTree(points)

        groups = rstTree.query_ball_tree(coordsTree, atol, p=pNorm)

        if strict:
            ident = np.array(map(len, groups))
            if not np.all(ident == 1):
                raise ValueError(
                    "%d points have not been found at all and %d points have"
                    " been found multiple times. If you don't need a"
                    " one-to-one match then set strict=False"
                    % ((ident == 0).sum(), (ident > 1).sum()))

        indices = np.unique(groups)

        if indices.dtype == object:
            # groups have different length --> dtypye == object, axis=0 fails
            indices = np.unique(list(chain.from_iterable(groups)))
            # alternativly but slower:
            # indices = np.unique(np.concatenate(groups))

        return indices, groups

    def getIdsInBox(self, box):
        """Returns a list of point indexes for all points inside the given
        box. Points on the border included (might actually be a matter of
        good luck due to rounding errors).

        @param box: [[min x, min y, min z], [max x, max y, max z]].
        """
        from scipy.spatial import cKDTree

        p1 = np.asarray(box[0])
        p2 = np.asarray(box[1])
        centre = 0.5*(p1+p2)
        radius = 0.5*(p2-p1)
        normCoords = (self.coords - centre) / radius
        tree = cKDTree(normCoords)
        return tree.query_ball_point(centre, 1.0, p=np.inf)

    def getBoundingBox(self):
        """Returns bounding box of stored coordinates

        @note: requires the attributes 'coords' to exist
        """
        return BoundingBox([self.coords.min(axis=0), self.coords.max(axis=0)])

#} end of Topology Classes
