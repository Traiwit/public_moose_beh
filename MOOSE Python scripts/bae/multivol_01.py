"""Classes that describe a multitude of volumes or a volume that is divided
into parts / subvolumes.

Mainly used to effectively assign points or elements of a mesh to those volume
parts.
"""

__version__ = "1.09"
_version_history_ = """
Versions:

1.01 new     - no version documented so far
1.02 changed - incompatible interface change: msgTicker argument to
     intersectionPointIdTuples()
1.03 changed - to bae.log_01 interface for message logging
1.04 GP added MultiVolRegularBoxGrid.pointIsInsideOnePart()
        added MultiVolRegularGridPoints, MultiVolClosestPoint
1.05 GP: added MultiVolUnique base class for speed improvements
     (less quecking required, simpler data structures)
1.06 GP added MultiVolRegularBoxGrid.getPartCentroid()
1.07 GP changed MultiVolSkewBoxGrid based on MultiVolUnique
        added numpy variant of MultiVolClosestPoint
1.08 GP fixed: old versions of scipy.spatial.cKDTree.query() have no n_jobs
        argument
1.09 GP added MultiVolBoundedRegularGrid class
"""

ideas = """
Possible alternative for class MultiVolClosestPoint.
Idea: Instead of just taking the closest point take a number of nearby points
and sum up the weighted distance --e.g. square_dist ** (-1.5)-- grouped by
volume-part id. Then assign the volume-part id with the largest sum.

Here is a sketch of the algorithm:
---snip---

# get data
data = np.genfromtxt(
    blockModelName, delimiter = ',', names=True,
    usecols=["X", "Y", "Z", "Rock_Type"],
    dtype=[float, float, float, int])

points = np.vstack((data["X"], data["Y"], data["Z"])).T
type_ = data["Rock_Type"]
msg("Read %d points and rock types from %s."
    % (len(points), blockModelName))

# init K-D-tree
kdt = KDTree(points)

# evaluate...
matCodeDict = defaultdict(float)
otherPtVals = list()
ticker = MsgTicker("Processed %%d/%d points." % len(otherPoints))
for pt in otherPoints:
    sqdIds = kdt.knnSearch(pt, 10)
    matCodeDict.clear()
    for sqd, idx2 in sqdIds:
        matCodeDict[type_[idx2]] += sqd**(-1.5)
    sortedMatCodes = sorted(matCodeDict.iteritems(), key=lambda x:x[1])
    otherPtVals.append(sortedMatCodes[-1][0])
    ticker.tick()
del ticker
--- snap ---

-------------------------------------------------------------------------------

MultiVolRegularGridPoints shall be incorporated into MultiVolBoundedRegularGrid
... and only the class interface be left for compatibilty purposes.

-------------------------------------------------------------------------------

remove the notFoundItem argument of MultiVol.intersectionElset(). Shall be
replaced by an optional outsideId-argument to __init__() like for
MultiVolBoundedRegularGrid and MultiVolRegularGridPoints.
"""


from math import pi, sin, cos, floor
from itertools import izip

try:
    import numpy as np
    import scipy.spatial
except ImportError:
    np = None

from bae.future_01 import defaultdict
from bae.vecmath_01 import vector_plus, vector_minus, mat_inverse, mat_multvec
from bae.volume_02 import VolumeFromMesh
from bae.abq_model_02 import checkModelModuleVersion, Model
from bae.misc_01 import KDTree, selectattrversion
from bae.log_01 import msg, MsgTicker


class MultiVol(object):
    """
    Base class for all types of MultiVol... classes
    """

    def __init__(self, logfile=None):
        """
        MultiVol Constructor

        @param logfile: deprecated, no function, for compatibility only
        """
        pass

    def intersectionPointList(self, points):
        r"""
        Assign each point in the given list to the volume parts it is in.

        @param points: list of points (each point being a coordinates triple)

        @returns: A dict {volume part identifier: list of indices of
        points in the given list points}.
        The kind of volume part identifier being used depends on the actual
        child class.

        Example with MultiVolRegularBoxGrid:

        >>> vol = MultiVolRegularBoxGrid(xyzStep=[10,10,10], xyzOrigin=[0,0,0])
        >>> points = [[12,3,-5], [1,0.5,45], ...]  # list of points
        >>> values = ["red", "green", "blue", ...] # corresponding values
        >>> volPartToIndex = vol.intersectionPointList(points)
        >>> # get a sorted list of all grid cells with points:
        >>> gridCellsWithPoints = sorted(volPartToIndex)
        >>> for multiIdx in gridCellsWithPoints:
        >>>     print "Grid cell %s contains:" % multiIdx
        >>>     print ", ".join(values[ptIdx]
        ...                     for ptIdx in volPartToIndex[multiIdx])
        """
        volPartToIndex = defaultdict(list)
        for idx, point in enumerate(points):
            for vol_part in self.pointIsInsideParts(point):
                volPartToIndex[vol_part].append(idx)
        volPartToIndex.default_factory = None  # convert to an ordinary dict
        return volPartToIndex

    def intersectionPointIdTuples(
            self, pointIdTuples,
            notFoundItem=None,
            msgTicker=None):
        r"""
        Assign each point with its corresponding identifier to the volume
        parts it is in. Takes a list of (point, point identifier) tuples and
        returns a dictionary {subvol id: set of point ids}.

        Example:
         >>> vol = MultiVolMeshSets("myVols.inp", ["DEV01", "DEV02", "DEV03"])
         >>> meshSetToSeqId = vol.intersectionPointIdTuples(
         ...     pointIdTuples=[([0.1,0.2,0.1], "yes"),  # in DEV01
         ...                    ([0.1,0.1,0.2], "yes"),  # in DEV01
         ...                    ([0.5,0.1,0.7], "no")],  # in none
         ...     notFoundItem="notFound")
         >>> notFoundPts = [tup[0] for tup in meshSetToSeqId["notFound"]]
         >>> del meshSetToSeqId["notFound"]
         >>> print notFoundPts
         [(0.5, 0.1, 0.7)]
         >>> print meshSetToSeqId
         {'DEV01': set(['yes'])}

        Example, sequence nodes to be embedded in sequenced host elsets:
         >>> from bae.multivol_01 import MultiVolMeshSets
         >>> from bae.abq_model_02 import Model
         >>>
         >>> hostMesh = Model().read("...")
         >>> hostElsetNames = hostMesh.regexpNames("elset", "CONCRETE_")
         >>>
         >>> muvol = MultiVolMeshSets(hostMesh, hostElsetNames)
         >>> pointIdTuples = [ (mesh2.nodeCoords[node], node)
         >>>                   for node in nodesToBeSeqd ]
         >>> seqNsets = muvol.intersectionPointIdTuples(
         >>>     pointIdTuples, notFoundItem="NOT_IN_HOST")
         >>>
         >>> if ("NOT_IN_HOST" in seqNsets):
         >>>     msg("WARNING: %d nodes not inside host elements."
         >>>         % len(seqNsets["NOT_IN_HOST"]))
         >>>     seqNsets["NOT_IN_HOST"] = set(
         >>>         y for x,y in seqNsets["NOT_IN_HOST"])
         >>>
         >>> mout = Model()
         >>> mout.nset.update(seqNsets)
         >>> mout.write("sequenced_nodes.inp")

        @param pointIdTuples: A list of (point, identifier)-tuples. point
          being a three-floats-tuple and identifier any hashable value to show
          up in the resulting dictionary

        @param notFoundItem: if not None specifies the key in the resulting
          dictionary under which all pointIdTuples are stored that are not
          found in any of the volume parts. Contrary to the ordinary entries
          in the resulting dictionary the set associated with the notFoundItem
          key does not contain just the ids from the pointIdTuples but the
          whole tuple of point and id. I.e:

          >>> vol = MultiVolMeshSets("myVols.inp", ["DEV01", "DEV02", "DEV03"])
          >>> print vol.intersectionPointIdTuples(
          ...     pointIdTuples=[([0.1,0.2,0.1], "yes"),([0.5,0.1,0.7], "no")],
          ...     notFoundItem="notFound")
          {'DEV01': set(['yes']), 'notFound': set([((0.5, 0.1, 0.7), 'no')])}

        @param msgTicker: DEPRECATED, no function, for compatibilty only

        @returns: a dict {volume part identifier : set of identifiers}
        The kind of volume part identifier being used depends on the actual
        child class.
        """
        msgTicker = MsgTicker(
            "processed %%d/%d points so far." % len(pointIdTuples))
        volPartToId = defaultdict(set)
        for point, ident in pointIdTuples:
            volParts = self.pointIsInsideParts(point)
            for thisVolPart in volParts:
                volPartToId[thisVolPart].add(ident)
            if notFoundItem is not None and len(volParts)==0:
                volPartToId[notFoundItem].add((tuple(point), ident))
            msgTicker.tick()
        del msgTicker
        volPartToId.default_factory = None  # convert to an ordinary dict
        return volPartToId

    def intersectionElset(self, model, elset=None, method='centroid',
                          notFoundItem=None):
        r"""
        Find all elements from the elset in model that are within any part
        of the volume.

        Usage (with a MultiVolMeshSets subclass):

        >>> from bae.multivol_01 import MultiVolMeshSets
        >>> from bae.abq_model_02 import Model
        >>> vol = MultiVolMeshSets("targets.inp", ["DEV01", "DEV02", "DEV03"])
        >>> newsets = vol.intersectionElset("source.inp")
        >>> mout = Model()
        >>> mout.elset.update(newsets)
        >>> mout.write("elsets_dev.inp")

        @param model: an abq_model_02.Model object or a filename (a string)
              or an open file object of an abaqus input file.

        @param elset: states the "searchset". Ignore the rest of model.
              might be anything that model.getUnionSet() accepts as input:
              an elset name, a set of element numbers, a list of elset names
              and element numbers
              elset might also be None, then all elements of model are tested.
              Element numbers in elset that are not defined in the model are
              silently ignored.

        @param method: optional method parameter. May be 'centroid'/'c' or
              'gauss point'/'gp'.

        @param notFoundItem: if not None specifies the elset name for elements
          not found in any of the volume parts. If None (the default), those
          elements are silently ignored.

        @return: a dict {volume part identifier: set of element numbers}
              The type of volume part identifier depends on the actual subclass
              of MultiVol. It is an elset name in the case of MultiVolMeshSets.
              It is a grid index triple in the case of MultiVolRegularBoxGrid.

        @note: bugs/limitations:
          - gauss point method not implemented yet.
        """

        method = method.lower()

        if isinstance(model, (basestring, file)):
            # is model any kind of string?
            m = Model()
            m.read(model)
            model = m
        else:
            checkModelModuleVersion(model, "%s.%s.intersectionElset"
                                    % (__name__, self.__class__.__name__))

        if elset is None:
            # caution, this is fast but dangerous if you modify the code
            # in this case here, elset is a dictionary, not a set
            elset = model.elNodes
        else:
            elset = model.getUnionSet("elset", elset, onlyValid=True)

        # initalize the result dict { volume part identifier: set of elnums }
        if method in ('centroid', 'c'):
            msg("Calculating centroids for %d elements." % len(elset))
            ptElemTup = [
                (centr, el)
                for el, centr in model.getElCentroids(elset).iteritems()]

        elif method in ('gauss point', 'gp'):
            msg("Calculating Gauss points for %d elements." % len(elset))
            ptElemTup = [
                (pt, el)
                for el, gaussPts in model.getElGaussPoints(elset).iteritems()
                for pt in gaussPts]

        else:
            raise Exception('MultiVol.intersectionElset: method <%s> not'
                            ' implemented.' % method)

        msg("Assigning elsets to %d elements." % len(elset))
        try:
            volPartToElset = self.intersectionPointIdTuples(
                pointIdTuples=ptElemTup, notFoundItem=notFoundItem)
        except TypeError:
            # MultiVolUnique.intersectionPointIdTuples has no argument
            # notFoundItem
            volPartToElset = self.intersectionPointIdTuples(
                pointIdTuples=ptElemTup)

        # fini
        if notFoundItem is not None and notFoundItem in volPartToElset:
            notFoundElset = set(id for pt, id in volPartToElset[notFoundItem])
            volPartToElset[notFoundItem] = notFoundElset

        return volPartToElset

    def writeIntersectionElset(self, output,
                               *intersectionElsetArgs,
                               **intersectionElsetKwArgs):
        r"""
        Find all elements from the elset in model that are within any part
        of the volume. Write those elsets to an abaqus input file.

        Usage (with a MultiVolMeshSets subclass):

        >>> from bae.multivol_01 import MultiVolMeshSets
        >>> vol = MultiVolMeshSets("targets.inp", ["DEV01", "DEV02", "DEV03"])
        >>> vol.writeIntersectionElset(
        >>>     model="source.inp", output="elsets_dev.inp")

        This is a convenience function that basically does the following:

        >>> newsets = self.intersectionElset(model, elset, method)
        >>> mout = Model()
        >>> mout.elset.update(newsets)
        >>> mout.write(output)

        All Arguments beside output (the first argument) are passed on to
        self.intersectionElset(). See their describtion there. Note that it is
        the intersectionElset-method of the actual MultiVol-subclass that is
        being called.

        I.e. for the arguments of MultiVolRegularBoxGrid.writeIntersectionElset
        look at MultiVolRegularBoxGrid.intersectionElset(). Plus the first
        that is described here:

        @param output: a file name (a string) or an open file object

        @note: You may (re)define the function self.getElsetNameFromId() to
          specify how the volume part identifiers returned by
          self.intersectionElset() should be converted to elset names.

          E.g.

          >>> from bae.multivol_01 import MultiVolMeshSets
          >>> vol = MultiVolMeshSets("targets.inp", ["DEV01", "DEV02", "DEV03"])
          >>> vol.getElsetNameFromId = lambda n: ("NEW_%s" % n)
          >>> vol.writeIntersectionElset(
          >>>     model="source.inp", output="elsets_dev.inp")

          Or:

          >>> from bae.multivol_01 import MultiVolRegularBoxGrid
          >>> vol = MultiVolRegularBoxGrid([99999, 100, 10], [-59999,None,None])
          >>> def myconvert(xyzTuple):
          >>>     return "BOX_%s%02d" % (chr(xyzTuple[1]+ord("A")), xyzTuple[2])
          >>> vol.getElsetNameFromId = myconvert
          >>> vol.writeIntersectionElset(
          >>>     model="source.inp", output="elsets_dev.inp")

        @note: uses MultiVol.intersectionElset() with all its limitations
        """
        newsets = self.intersectionElset(*intersectionElsetArgs,
                                         **intersectionElsetKwArgs)
        if not newsets:
            msg("WARNING: writeIntersectionElset did not find any element."
                " Elsets file $s not written." % output)
            return

        mout = Model()
        try:
            getElsetNameFromId = getattr(self, "getElsetNameFromId")
        except AttributeError:
            # use volume part names as new elset names
            if not isinstance(newsets.iterkeys().next(), basestring):
                raise ValueError(
                    "ERROR: writeIntersectionElset cannot determine elset name"
                    " from subvolume identifier. Supply a getElsetNameFromId()"
                    " method for the class %s or this object. Identifier"
                    " example: %s"
                    % (type(self), newsets.iterkeys().next()))
            mout.elset.update(newsets)
        else:
            # use self.getElsetNameFromId()
            for key, elset in newsets.iteritems():
                mout.elset[getElsetNameFromId(key)] = elset
        mout.write(output, withSummary=True)


class MultiVolUnique(MultiVol):
    """
    Base class for MultiVol types that consist of mutually disjoint
    parts / subvolumes that fill the whole 3D space. (That can assign
    each and every point in the 3D space to exactly one volume part.)

    Examples are L{MultiVolRegularBoxGrid} and L{MultiVolClosestPoint}

    Subclasses define a method self.pointIsInsideOnePart(point) that
    should be prefered over self.pointIsInsideParts(point)
    """

    def intersectionPointIdTuples(self, pointIdTuples):
        r"""
        Assign each point with its corresponding identifier to the volume
        parts it is in. Takes a list of (point, point identifier) tuples and
        returns a dictionary {subvol id: set of point ids}.

        Example:

        >>> vol = MultiVolRegularBoxGrid([1, 1, 1])
        >>> cellToColour = vol.intersectionPointIdTuples(
        ...     [([0.1,0.2,0.1], "green"),  # in (0,0,0)
        ...      ([2.1,1.1,1.2], "green"),  # in (2,1,1)
        ...      ([1.8,1.1,0.7], "red"),    # in (2,1,1)
        ...      ([1.9,0.9,0.7], "blue"),   # in (2,1,1)
        ...     ] )
        >>> print cellToColour
        {(0,0,0): set(['green']), (2,1,1): set(['green','red','blue'])}

        @param pointIdTuples: A list of (point, identifier)-tuples. point
          being a three-floats-tuple and identifier any hashable value to show
          up in the resulting dictionary

        @returns: a dict {volume part identifier : set of identifiers}
        The kind of volume part identifier being used depends on the actual
        child class.
        """
        msgTicker = MsgTicker(
            "processed %%d/%d points so far." % len(pointIdTuples))
        volPartToId = defaultdict(set)
        for point, ident in pointIdTuples:
            volPart = self.pointIsInsideOnePart(point)
            volPartToId[volPart].add(ident)
            msgTicker.tick()
        del msgTicker
        volPartToId.default_factory = None  # convert to an ordinary dict
        return volPartToId

    def intersectionElset(self, model, elset=None, method='centroid'):
        r"""
        Find all elements from the elset in model that are within any part
        of the volume.

        Usage (with a MultiVolMeshSets subclass):

        >>> from bae.multivol_01 import MultiVolMeshSets
        >>> from bae.abq_model_02 import Model
        >>> vol = MultiVolMeshSets("targets.inp", ["DEV01", "DEV02", "DEV03"])
        >>> newsets = vol.intersectionElset("source.inp")
        >>> mout = Model()
        >>> mout.elset.update(newsets)
        >>> mout.write("elsets_dev.inp")

        @param model: an abq_model_02.Model object or a filename (a string)
              or an open file object of an abaqus input file.

        @param elset: states the "searchset". Ignore the rest of model.
              might be anything that model.getUnionSet() accepts as input:
              an elset name, a set of element numbers, a list of elset names
              and element numbers
              elset might also be None, then all elements of model are tested.
              Element numbers in elset that are not defined in the model are
              silently ignored.

        @param method: optional method parameter. May be 'centroid'/'c' or
              'gauss point'/'gp'.
              currently being ignored

        @return: a dict {volume part identifier: set of element numbers}
              The type of volume part identifier depends on the actual subclass
              of MultiVol. It is a grid index triple in the case of
              MultiVolRegularBoxGrid.
        """
        # call super-class method but hide its notFoundItem parameter!
        return MultiVol.intersectionElset(
            self, model, elset=elset, method=method)


class MultiVolMeshSets(MultiVol):
    """
    Volume parts defined by certain elsets in an abq_model_02.Model

    Usage:

    >>> from bae.multivol_01 import MultiVolMeshSets
    >>> from bae.abq_model_02 import Model
    >>> vol = MultiVolMeshSets("targets.inp", ["DEV01", "DEV02", "DEV03"])
    >>> vol.writeIntersectionElset("source.inp", "elsets_dev.inp")

    @ivar allPartIds: Names of all elsets / volume parts

    @Note: This should be the preferred means over multiple
    volume_02.VolumeFromMesh instances if the volumes / volume parts intersect.
    If they don't intersect, multiple volume_02.VolumeFromMesh instances may be
    preferable.
    """
    def __init__(self, model, elsetNames, logfile=None):
        """
        @param model: an abq_model_02.Model object or a filename of an
              abaqus input file (a string in the latter case).
        @param elsetNames: a list of elsetnames defining the volume parts.
        @param logfile: deprecated, no function, for compatibility only

        @note: If you want all elsets of model, use model.elset.keys() as
              argument elsetNames. None is not recognized as all elsets, in
              contrast e.g. to the elset argument of VolumeFromMesh().
        """
        MultiVol.__init__(self)

        self.allPartIds = elsetNames
        self.volume = VolumeFromMesh(model, elset=elsetNames)

        self.elemElset = defaultdict(list)
        ticker = MsgTicker("Initializing element to elset relation, finished"
                           " %s/%d elsets.\n" % ("%d", len(self.allPartIds)))
        for cnt, elsetName in enumerate(self.allPartIds):
            newElset = model.elset[elsetName]
            for elem in newElset:
                self.elemElset[elem].append(elsetName)
            ticker.msg(cnt+1)
        del ticker
        # convert to tuples so they may not accidently be modified
        for elem in self.elemElset.keys():
            self.elemElset[elem] = tuple(self.elemElset[elem])

    def initializePointSearch(self):
        """This may be called to trigger the point search initialization.
        This is not necessary as it is done automatically on first invocation
        of self.pointIsInsideParts().
        """
        self.volume.initializePointSearch()

    def pointIsInsideParts(self, point):
        """
        returns a list of elset names (volume part identifiers) in which
        the specified points may be found.
        """
        elNum = self.volume.pointIsInside(point, returnVolElem=True)
        try:
            volName = self.elemElset[elNum]
        except KeyError:
            volName = []
        return volName

    def scale(self, scale_factor, scale_origin):
        """
        Scale the volume by a certain factor relative to scale_origin

        scale_factor may be a single number in which case it is applied to all
        three dimensions or a list of three numbers with one factor for each
        dimension.

        Example:

        >>> # scale this volume by 2 in each direction (volume := volume*8)
        >>> v.scale(scale_factor=2, scale_origin=[0,0,0])

        Scaling is done by simply moving the point coordinates for this
        volume. Other volumes or the original mesh from which the volume
        originates are not affected.
        """
        return self.volume.scale(scale_factor, scale_origin)

    def translate(self, move_vector):
        """
        Move the volume in space by move_vector

        Example:

        >>> # move this volume up by 10
        >>> vol.translate(move_vector=[0,0,10])

        Translating is done by simply moving the point coordinates for this
        volume. Other volumes or the original mesh from which the volume
        originates are not affected.
        """
        return self.volume.translate(move_vector)


class MultiVolRegularBoxGrid(MultiVolUnique):
    """
    Volume Parts defined as disjoint boxes with identical shape and size
    aligned to the coordinate axes. They virtually strech out to infinity (in
    all directions)

    The identifier of the volume parts is a tuple of three integer indices.

    This class is basically a grid aligned to the coordinate axes with a
    certain origin and mesh size (xyzStep). It's capable of identifying the
    grid cell a given query point lies in.

    In accordance with the other MultiVol classes a grid cell is also referred
    to as volume part or subvolume. The index triple identifying a grid cell is
    also referred to as subvolume identifier.
    """

    def __init__(self, xyzStep, xyzOrigin=[None, None, None]):
        """
        @param xyzStep: x-y-z tuple of the width of a single volume part
        @param xyzOrigin: origin of the first box; the box (volume part) with
          the identifier (0,0,0) will start at x,y,z = xyzOrigin and end at
          x,y,z = xyzOrigin + xyzStep.

          Any component specified as None defaults to minus half the
          corresponding xyzStep component so you just have to specifiy a
          ridiculously big xyzStep component to not discriminate along this
          axis.

          Defaults to [None, None, None].

        @note: If you don't want to seperate points along a certain axis, you
          would specify a large component in xyzStep for that axis. Then
          self.pointIsInsideParts() only discriminates points according to the
          other coordinates.
        """
        MultiVolUnique.__init__(self)
        self.xyzStep = tuple(float(x) for x in xyzStep)
        self.xyzOrigin = list(
            ((x is not None) and float(x)) or x
            for x in xyzOrigin)
        for i in range(3):
            if self.xyzOrigin[i] is None:
                self.xyzOrigin[i] = -0.5*self.xyzStep[i]

    def pointIsInsideParts(self, point):
        """
        Check in which subvolumes (parts) the specified point is located

        @param point: triple of coordinates

        @returns: A list containing one(!) tuple of three integer indexes
        identifying the grid box the point lies in. Why is it a list and not
        just the single index-triple? To make it compatible with the
        corresponding methods of other MultiVol-classes.
        """
        ixyz = [int(floor(float(point[i]-self.xyzOrigin[i])/self.xyzStep[i]))
                for i in range(3)]
        return [tuple(ixyz), ]

    def pointIsInsideOnePart(self, point):
        """
        Check in which subvolume (part) the specified point is located.
        This is special to Multivol classes that have mutually disjoint
        subvolumes and is equivalent to self.pointIsInsideParts(point)[0]

        @param point: triple of coordinates

        @returns: A tuple of three integer indexes identifying the grid box
        the point lies in.
        """
        return tuple(
            int(floor(float(point[i]-self.xyzOrigin[i])/self.xyzStep[i]))
            for i in range(3))

    def getPartCentroid(self, ijk):
        """Return the centroid of the specified multivol-part a.k.a. cell.

        @param ijk: part/cell index. A tuple of three integers as returned by
           L{pointIsInsideParts} or L{pointIsInsideOnePart}.

        @returns: point coordinates. I.e. a list of three floats for x,y,z.
        """
        if len(ijk)!=3:
            raise ValueError(
                "ijk-argument to getPartCentroid must be a three-tuple."
                " Instead we got: %s" % ijk)
        return [x0 + dx*(ix+0.5)
                for x0, dx, ix in izip(self.xyzOrigin, self.xyzStep, ijk)]

    def intersectionElset(self, model, elset=None, method='centroid',
                          startIndices=[None, None, None]):
        """
        Associate all elements from the elset in model to the volume part
        they are in

        Usage:

        >>> from bae.multivol_01 import MultiVolRegularBoxGrid
        >>> vol = MultiVolRegularBoxGrid([1, 5, 99999])
        >>> newsets = vol.intersectionElset(
        >>>     "source.inp", startIndices=[0, None, None])
        >>> mout = Model()
        >>> mout.elset.update(newsets)
        >>> mout.write("elsets_dev.inp")

        @param model: an abq_model_02.Model object or a filename of an
              abaqus input file (a string in the latter case).

        @param elset: states the "searchset". Ignore the rest of model.
              might be anything that model.getUnionSet() accepts as input:
              an elset name, a set of element numbers, a list of elset names
              and element numbers
              elset might also be None, then all elements of model are tested.
              Element numbers in elset that are not defined in the model are
              silently ignored.

        @param method: optional method parameter
              May be 'centroid'/'c' or 'gauss point'/'gp'

        @param startIndices: Lowest index values for the x,y,z index of the
          boxes. That means, the xyzOrigin is virtually adapted (in steps
          according to xyzSteps) so that the indices match the specified
          condition. A component None leads to indices according to the
          specified xyzOrigin on the corresponding axis.

          (xyzOrigin and xyzSteps are arguments to the constructor)

        @return: a dict {volume part identifier: set of element numbers}
              The volume part identifier depends on the actual subclass of
              MultiVol. It is an elset name in the case of MultiVolMeshSets.

        @note: bugs/limitations:
          - gauss point method not implemented yet.
        """

        method = method.replace(' ', '')

        if isinstance(model, basestring):  # is model any kind of string?
            m = Model()
            m.read(model)
            model = m
        else:
            checkModelModuleVersion(model, "%s.%s.intersectionElset"
                                    % (__name__, self.__class__.__name__))

        if elset is None:
            # caution, this is fast but dangerous if you modify the code
            # in this case here, elset is a dictionary, not a set
            elset = model.elNodes
        else:
            elset = model.getUnionSet("elset", elset, onlyValid=True)

        # initalize check-point elem association centroidElemTup
        if method in ('centroid', 'c'):
            centroidElemTup = list()
            for el in elset:
                centroid = model.getElCentroid(el)
                centroidElemTup.append((centroid, el))

        elif method in ('gauss point', 'gp'):
            raise Exception('MultiVolRegularBoxGrid.intersectionElset: method'
                            ' <%s> not implemented so far.' % method)
            #----- the following might help you if you want to implement this
            # from volume_02
            # create list with tuples (el numer, [points to check])
            el_pts = [(el, [model.getElGaussPoints(el), "no, other points!"])
                      for el in elset]
            # check if all checkpoints are in the volume
            els_inside = set()
            for el, points in el_pts:
                if all([self.pointIsInside(pt) for pt in points]):
                    els_inside.add(el)
            return els_inside

        else:
            raise Exception('MultiVolRegularBoxGrid.intersectionElset: method'
                            ' <%s> not implemented so far.' % method)

        # create the result dict { volume part identifier: set of elnums }
        volPartToElset = self.intersectionPointIdTuples(centroidElemTup)

        # correct the volume part indices to meet the requirement defined
        # by startIndices
        # this is done the straight way: 1. compute the current start indices
        # 2. compute offsets to the supposed values 3. correct the indices
        indexOffsets = [0,0,0]
        if any([idx is not None for idx in startIndices]):
            currentStart = None
            for ixyz in volPartToElset:
                try:
                    for i in range(3):  # x,y,z
                        if ixyz[i]<currentStart[i]:
                            currentStart[i]=ixyz[i]
                except TypeError:
                    currentStart=list(ixyz)

            for i in range(3):
                if startIndices[i] is not None:
                    indexOffsets[i] = startIndices[i]-currentStart[i]

        if any(indexOffsets):
            newVolPartToElset = dict()
            for ixyz, thisElset in volPartToElset.iteritems():
                newkey = tuple(vector_plus(ixyz,indexOffsets))
                newVolPartToElset[newkey] = thisElset
            volPartToElset = newVolPartToElset

        # fini
        return volPartToElset

    def _intersectionElset_startpointmethod(
            self, model, elset=None, method='centroid',
            startIndices=[None, None, None]):
        """
        associate all elements from the elset in model to the volume part
        they are in

        This function is the older version of
        MultiVolRegularBoxGrid.intersectionElset(). It does not work with the
        subclass MultiVolRegularBoxGridRotZ and a nondefault value for
        parameter startIndices. But it may be faster, I don't know. It saves
        the rearranging of all indices to meet the specified startIndices
        condition but on the other hand has to calculate the actual startpoint
        of all elsets to intersect.

        The result is also different: this Version starts its first block at
        the beginning of the actual elset, whereas the other version starts
        at the specified xyzOrigin + N*xyzSteps where N is such that the first
        index equals the specified startIndex.

        You can use this function but remember that it does not work with the
        subclass MultiVolRegularBoxGridRotZ.

        Usage:

        >>> from bae.multivol_01 import MultiVolRegularBoxGrid
        >>> vol = MultiVolRegularBoxGrid([1, 5, 99999])
        >>> newsets = vol.intersectionElset(
        >>>     "source.inp", startIndices=[0, None, None])
        >>> mout = Model()
        >>> mout.elset.update(newsets)
        >>> mout.write("elsets_dev.inp")

        @param model: an abq_model_02.Model object or a filename of an
              abaqus input file (a string in the latter case).

        @param elset: states the "searchset". Ignore the rest of model.
              might be anything that model.getUnionSet() accepts as input:
              an elset name, a set of element numbers, a list of elset names
              and element numbers
              elset might also be None, then all elements of model are tested.
              Element numbers in elset that are not defined in the model are
              silently ignored.

        @param method: optional method parameter
              May be 'centroid'/'c' or 'gauss point'/'gp'

        @param startIndices: Lowest index values for the x,y,z index of the
          boxes. That means, the xyzOrigin is virtually adapted so that the
          indices match the specified condition. A component None leads to
          indices according to the specified xyzOrigin (constructor argument)
          on the corresponding axis.

        @return: a dict {volume part identifier: set of element numbers}
              The volume part identifier depends on the actual subclass of
              MultiVol. It is an elset name in the case of MultiVolMeshSets.

        @note: bugs/limitations:
          - gauss point method not implemented yet.
        """

        method = method.replace(' ', '')

        if isinstance(model, basestring):  # is model any kind of string?
            m = Model()
            m.read(model)
            model = m
        else:
            checkModelModuleVersion(model, "%s.%s.intersectionElset"
                                    % (__name__, self.__class__.__name__))

        if elset is None:
            # caution, this is fast but dangerous if you modify the code
            # in this case here, elset is a dictionary, not a set
            elset = model.elNodes
        else:
            elset = model.getUnionSet("elset", elset, onlyValid=True)

        # initalize check-point elem association centroidElemTup
        if method in ('centroid', 'c'):
            centroidElemTup = list()
            startCorner = None
            for el in elset:
                centroid = model.getCentroid(el)
                centroidElemTup.append((centroid, el))
                try:
                    for i in range(3):  # x,y,z
                        if ((self.xyzStep[i]>0)==(centroid[i]<startCorner[i])):
                            startCorner[i]=centroid[i]
                except TypeError:
                    startCorner=list(centroid)

        elif method in ('gauss point', 'gp'):
            raise Exception('MultiVolRegularBoxGrid.intersectionElset: method'
                            ' <%s> not implemented so far.' % method)
            #----- the following might help you if you want to implement this
            # from volume_02
            # create list with tuples (el numer, [points to check])
            el_pts = [(el, [model.getCentroid(el), "no, other points!"])
                      for el in elset]
            # check if all checkpoints are in the volume
            els_inside = set()
            for el, points in el_pts:
                if all([self.pointIsInside(pt) for pt in points]):
                    els_inside.add(el)
            return els_inside

        else:
            raise Exception('MultiVolRegularBoxGrid.intersectionElset: method'
                            ' <%s> not implemented so far.' % method)

        # save original xyzOrigin and modify according to startIndices
        saveOrigin = self.xyzOrigin
        self.xyzOrigin = list(self.xyzOrigin)
        for i in range(3):
            if startIndices[i] is not None:
                self.xyzOrigin[i] = (
                    startCorner[i] - (1E-3 + startIndices[i])*self.xyzStep[i])

        # create the result dict { volume part identifier: set of elnums }
        volPartToElset = self.intersectionPointIdTuples(centroidElemTup)

        # restore original xyzOrigin
        self.xyzOrigin = saveOrigin

        # fini
        return volPartToElset

    getElsetNameFromIdPattern = "BOX_%02d_%02d_%02d"
    def getElsetNameFromId(self, partId):
        """Create an elset name from the part id.
        Used by self.L{writeIntersectionElset}(), see description there.

        @returns: self.getElsetNameFromIdPattern % partId
        """
        return self.getElsetNameFromIdPattern % partId


class MultiVolRegularBoxGridRotZ(MultiVolRegularBoxGrid):
    """
    Volume Parts defined as disjoint boxes with identical shape and size
    aligned to a coordinate system turned around the z-axis by a specified
    angle. The boxes virtually strech out to infinity (in all directions)

    The identifier of the volume parts is a tuple of three integer indices.

    This class is basically a grid with a certain origin and mesh size (xyzStep)
    aligned to the axes of a coordinate system rotated by a certain angle
    around the z axis. It's capable of identifying the grid cell a given query
    point lies in.
    """

    def __init__(self, xyzStep, xyzOrigin=[None, None, None], rotZ=0.0):
        """
        @param xyzStep: x-y-z tuple of the width of a single volume part

        @param xyzOrigin: origin of the first box; the box (volume part) with
          the identifier (0,0,0) will start at x,y,z = xyzOrigin and end at
          xyz = xyzStep.

          Any component specified as None defaults to minus half the
          corresponding xyzStep component so you just have to specifiy a
          ridiculously big xyzStep component to not disciminate along this axis.

          xyzOrigin also serves as rotation centre, around which the box is
          rotated by rotZ degrees. If the x- or y-component of xyzOrigin is
          None, the corresponding component of the rotation centre is set to
          zero.

          Defaults to [None, None, None].

        @param rotZ: angle in degrees from the (global) x-axis to the axis of
          the first index of the subvolumes (parts)

        @note: If you don't want to seperate points along a certain axis, you
          would specify a large component in xyzStep for that axis. Then
          self.pointIsInsideParts() only discriminates points according to the
          other coordinates.
        """
        self.rotCenter = xyzOrigin[:2]
        for i in range(2):
            if xyzOrigin[i] is None:
                self.rotCenter[i] = 0.0

        MultiVolRegularBoxGrid.__init__(
            self, xyzStep=xyzStep, xyzOrigin=xyzOrigin)
        alpha = rotZ*pi/180
        self.ca=cos(alpha)
        self.sa=sin(alpha)
        # rotated vec: [ca*v[0]+sa*v[1], -sa*v[0]+ca*v[1], v[2]]

    def pointIsInsideParts(self, point):
        """
        check in which subvolumes (parts) the specified point is located

        @param point: 3 tuple of coordinates

        @returns: list of subvolume identifiers (in this case the list has
        always exactly one item). The subvolume identifier is a tuple of three
        integers specifying the the three indexes of the subvolume.
        """
        v = [point[0]-self.rotCenter[0],
             point[1]-self.rotCenter[1],
             point[2]-self.xyzOrigin[2]]
        v = [self.ca*v[0]+self.sa*v[1],
             -self.sa*v[0]+self.ca*v[1],
             float(v[2])]
        ixyz = (int(floor(v[0]/self.xyzStep[0])),
                int(floor(v[1]/self.xyzStep[1])),
                int(floor(v[2]/self.xyzStep[2])))
        return [ixyz, ]

    def pointIsInsideOnePart(self, point):
        """
        Check in which subvolume (part) the specified point is located.
        This is special to Multivol classes that have mutually disjoint
        subvolumes and is equivalent to self.pointIsInsideParts(point)[0]

        @param point: triple of coordinates

        @returns: A tuple of three integer indexes identifying the grid box
        the point lies in.
        """
        v = [point[0]-self.rotCenter[0],
             point[1]-self.rotCenter[1],
             point[2]-self.xyzOrigin[2]]
        v = [self.ca*v[0]+self.sa*v[1],
             -self.sa*v[0]+self.ca*v[1],
             float(v[2])]
        return (int(floor(v[0]/self.xyzStep[0])),
                int(floor(v[1]/self.xyzStep[1])),
                int(floor(v[2]/self.xyzStep[2])))


class MultiVolBoundedRegularGrid(MultiVolUnique):
    """
    Volume Parts defined as regular rectangular grid aligned to the coordinate
    axes. In contrast to L{MultiVolRegularBoxGrid} the cells are confined to a
    certain bounding box.

    The identifier of the volume parts is a tuple of three integer indices.
    Optionally a field of values serving as volume identifiers can be passed.
    In this case volume part identifiers are taken from this field.

    Points outside the bounding box are assigned to "outsideId" as given by the
    corresponding argument to __init__().

    In accordance with the other MultiVol classes a grid cell is also referred
    to as volume part or subvolume. The index triple identifying a grid cell is
    also referred to as subvolume identifier.

    Example:
     >>> grid = MeshStructuredPoints(
     >>>     origin=[1000, 200, -500],
     >>>     spacing=[5, 5, 5],
     >>>     gridPtNb=[12, 24, 13])
     >>> subVolIdField = createFieldObject(
     >>>     "partVol", "structPt", "str",
     >>>     initArgs=[(("box %d %d %d" % (int(i/3), int(j/6), int(k/2)))
     >>>                for i, j, k in grid.getGridIdxIter()),])
     >>> muvol = MultiVolBoundedRegularGrid(
     >>>     pointsGrid=grid, subVolIdField=subVolIdField)
     >>>
     >>> pointIdTuples = [
     >>>     (pt, i+1) for i, pt in enumerate((
     >>>         [1010, 240, -430], ... ))]
     >>>
     >>> volPartToId = muvol.intersectionPointIdTuples(pointIdTuples)
     >>>
     >>> for volPart in sorted(vorPartToId):
     >>>     ptIds = volPartToId[volPart]
     >>>     print "In %s we have points nb %s" % (volPart, sorted(ptIds))

    Note: For efficiency reasons points on the lower/left grid cell boundary
    are considered in this cell whereas points on the upper/right boundary are
    considered out (and in the next cell if there is one).
    """

    def __init__(self, pointsGrid, subVolIdField=None, outsideId=None):
        """
        @param pointsGrid: a L{bae.mesh_01.MeshStructuredPoints} object
          identifying the cell centres of the regular grid. Pass as keyword
          argument, see note below.

        @param subVolIdField: a L{bae.field_01.Field} object containing the
          subvolume identifiers. Actually a list of values is suffient as well.

        @param outsideId: Value to be returned for points outside any cell for
          which an id has been specified. self.pointIsInsideOnePart() would
          return this value whereas self.pointIsInsideParts() would return an
          empty list if the point is not in any of the cells identified by
          the points argument. Defaults to None.

        @Note: IMPORTANT: It's advised to supply all arguments as keyword
          arguments because later versions might make the pointsGrid argument
          optional and support different ways to specify the grid.

        @Note: If you don't want to seperate points along a certain axis, then
          make the xyzStep-component of that exis (this particular spacing
          component of the grid) very large.
        """
        MultiVolUnique.__init__(self)
        self.grid = MultiVolRegularBoxGrid(
            xyzStep=pointsGrid.spacing,
            xyzOrigin=[
                x-0.5*dx
                for x, dx in zip(pointsGrid.origin, pointsGrid.spacing)])
        self.gridPtNb = pointsGrid.gridPtNb
        self.strides = pointsGrid.strides
        self.topCorner = [
            self.grid.xyzOrigin[i]  # this is the lower corner
            + pointsGrid.gridPtNb[i]*pointsGrid.spacing[i]  # plus N cell widths
            for i in range(3)]
        # to get the field (sequential) index from the index tuple:
        # sum( i*s for i, s in zip(index, self.strides) )

        self.subVolIdField = subVolIdField
        self.outsideId = outsideId

    def pointIsInsideParts(self, point):
        """
        Check in which subvolumes (parts) the specified point is located

        @param point: triple of coordinates

        @returns: A list containing one(!) tuple of three integer indexes
        identifying the grid box the point lies in. Or the corresponding
        item from the field of subvolume identifiers if it has been passed
        to __init__(). Or an empty list if the point is outside the bounding
        box of the grid.

        Why is it a list and not just the single index-triple? To make it
        compatible with the corresponding methods of other MultiVol-classes.
        """
        origin = self.grid.xyzOrigin
        spacing = self.grid.xyzStep
        if not all(origin[i]<=point[i]<self.topCorner[i] for i in range(3)):
            return []
        id_ = tuple(int((point[i]-origin[i])/spacing[i])
                    for i in range(3))
        if self.subVolIdField:
            id_ = sum(i*s for i, s in zip(id_, self.strides))
            id_ = self.subVolIdField[id_]
        return [id_, ]

    def pointIsInsideOnePart(self, point):
        """
        Check in which subvolume (part) the specified point is located.
        This is special to L{MultiVolUnique} subclasses that have mutually
        disjoint subvolumes and is roughly equivalent to
        self.pointIsInsideParts(point)[0].

        @param point: triple of coordinates

        @returns: A tuple of three integer indexes identifying the grid box
        the point lies in. Or the corresponding item from the field of
        subvolume identifiers if it has been passed to __init__(). Or
        the outsideId argument to __init__ if the point is outside the
        bounding box of the grid.
        """
        origin = self.grid.xyzOrigin
        spacing = self.grid.xyzStep
        if not all(origin[i]<=point[i]<self.topCorner[i] for i in range(3)):
            return self.outsideId
        id_ = tuple(int((point[i]-origin[i])/spacing[i])
                    for i in range(3))
        if self.subVolIdField:
            id_ = sum(i*s for i, s in zip(id_, self.strides))
            id_ = self.subVolIdField[id_]
        return id_

    def intersectionPointIdTuples(self, pointIdTuples):
        r"""
        Assign each point with its corresponding identifier to the volume
        parts it is in. Takes a list of (point, point identifier) tuples and
        returns a dictionary {subvol id: set of point ids}.

        Example:
         >>> grid = ...
         >>> vol = MultiVolRegularBoxGrid(grid)
         >>> cellToColour = vol.intersectionPointIdTuples(
         ...     [([0.1,0.2,0.1], "green"),  # in (0,0,0)
         ...      ([2.1,1.1,1.2], "green"),  # in (2,1,1)
         ...      ([1.8,1.1,0.7], "red"),    # in (2,1,1)
         ...      ([1.9,0.9,0.7], "blue"),   # in (2,1,1)
         ...     ] )
         >>> print cellToColour
         {(0,0,0): set(['green']), (2,1,1): set(['green','red','blue'])}

        @param pointIdTuples: A list (iterator is not sufficient) of
          (point, identifier)-tuples. point being a three-floats-tuple and
          identifier any hashable value to show up in the resulting dictionary.

        @returns: a dict {volume part identifier : set of identifiers}
          Each volume part identifier might be an index-triple or an item of
          the sub-volume-identifier-field passed in as the subVolIdField
          argument to self.__init__(). Or outsideId.
        """

        # abbreviations
        origin = self.grid.xyzOrigin
        spacing = self.grid.xyzStep

        volPartToId = defaultdict(set)

        # find ousiders
        outsiders = set(
            id_
            for point, id_ in pointIdTuples
            if any(point[i]<origin[i] for i in range(3))
            or any(point[i]>=self.topCorner[i] for i in range(3)))
        volPartToId[self.outsideId] = outsiders

        # find the rest
        if self.subVolIdField:
            for point, ptId in pointIdTuples:
                if ptId not in outsiders:
                    volId = sum(
                        int((point[i]-origin[i])/spacing[i])*self.strides[i]
                        for i in range(3))
                    volPartToId[self.subVolIdField[volId]].add(ptId)
        else:
            for point, ptId in pointIdTuples:
                if ptId not in outsiders:
                    volId = tuple(int((point[i]-origin[i])/spacing[i])
                                  for i in range(3))
                    volPartToId[volId].add(ptId)

        # end
        volPartToId.default_factory = None  # convert to an ordinary dict
        return volPartToId


class MultiVolSkewBoxGrid(MultiVolUnique):
    """
    Volume Parts defined as disjoint boxes with identical shape and size
    aligned to three given base vectors. They virtually stretch out to
    infinity (in all directions).

    The identifier of the volume parts is a tuple of three integer indices.

    This class is basically a skew grid with a certain origin. The grid axes
    and mesh sizes are given by three base vectors. Objects of this class are
    suitable for identifying the grid cell a given query point lies in.
    """

    def __init__(self, baseVectors, xyzOrigin):
        """
        @param baseVectors: a list of three base vectors
        @param xyzOrigin: origin of the first box; the box (volume part) with
          the identifier (0,0,0) will have one corner at x,y,z = xyzOrigin and
          the opposite corner at x,y,z = xyzOrigin + sum_i(baseVectors[i]).
        """
        MultiVol.__init__(self)

        self.mat = mat_inverse(baseVectors)
        self.xyzOrigin = list(xyzOrigin)

    def pointIsInsideParts(self, point):
        """
        check in which subvolumes (parts) the specified point is located

        @param point: 3-tuple of coordinates

        @returns: list of subvolume identifiers (in this case the list has
        always exactly one item). The subvolume identifier is a tuple of three
        integers specifying the three indexes of the subvolume.
        """
        dx = vector_minus(point, self.xyzOrigin)
        ixyz = map(int, mat_multvec(self.mat,dx))
        return [tuple(ixyz), ]

    def pointIsInsideOnePart(self, point):
        """
        Check in which subvolume (part) the specified point is located.
        This is special to Multivol classes that have mutually disjoint
        subvolumes and is equivalent to self.pointIsInsideParts(point)[0]

        @param point: triple of coordinates

        @returns: A tuple of three integer indexes identifying the grid box
        the point lies in.
        """
        dx = vector_minus(point, self.xyzOrigin)
        return tuple(int(x) for x in mat_multvec(self.mat,dx))

    getElsetNameFromIdPattern = "BOX_%02d_%02d_%02d"
    def getElsetNameFromId(self, partId):
        """Create an elset name from the part id.
        Used by self.writeIntersectionElset(), see description there.

        @returns: self.getElsetNameFromIdPattern % partId
        """
        return self.getElsetNameFromIdPattern % partId


class MultiVolRegularGridPoints(MultiVolUnique):
    """
    Volume parts defined by point data. The points are given trough their
    individual coordinates but still are expected to be on a regular grid
    aligned to the coordinate axes. Each point carries a value serving as
    subvolume identifier.

    The points are assumed to be at the centre of a particular grid cell. Each
    subvolume consists of the grid cells with the corresponding subvolume
    identifier.
    """

    def __init__(self, points, ids,
                 xyzStep=[None, None, None], xyzOrigin=None,
                 outsideId=None):
        """
        @param points: list of coordinate tuples, i.e. [[x1,y1,z1], [x2,y2,z2],
          ...]
        @param ids: list of subvolume identifiers.
        @param xyzStep: x-y-z tuple of the width of a single grid cell
          Currently this argument has to be supplied. In later version it might
          be calculated from the points data, if someone finds a clever
          algorithm for that.
        @param xyzOrigin: centre of a reference cell; all cells will have their
          centre at [x0+i*dx, y0+j*dy, z0 + k*dz] with xyzOrigin==[x0,y0,z0],
          xyzStep==[dx,dy,dz] and i,j,k being three integers.

          If xyzOrigin is a list or tuple then any component specified as None
          defaults to minus half the corresponding xyzStep component so you just
          have to specifiy a ridiculously big xyzStep component to not
          discriminate along this axis.

          If xyzOrigin is None then the first point in points is being taken.

          Defaults to None.

        @param outsideId: Value to be returned for points outside any cell for
          which an id has been specified. self.pointIsInsideOnePart() would
          return this value whereas self.pointIsInsideParts() would return an
          empty list if the point is not in any of the cells identified by
          the points argument. Defaults to None.

        @note: If you don't want to seperate points along a certain axis, you
          would specify a large component in xyzStep for that axis. Then
          self.pointIsInsideParts() only discriminates points according to the
          other coordinates.
        """
        if len(points)!=len(ids):
            raise ValueError(
                "points and ids arguments don't appear to have the same"
                " length. points: %d, ids: %d" % (len(points), len(ids)))

        MultiVol.__init__(self)
        try:
            xyzStep = tuple(float(x) for x in xyzStep)
        except (TypeError, ValueError):
            raise ValueError("Currently the xyzStep argument must be provided."
                             " xyzStep=%s" % xyzStep)
        if xyzOrigin is None:
            xyzOrigin = list(xyzStep)

        self.grid = MultiVolRegularBoxGrid(
            xyzStep=xyzStep,
            xyzOrigin=[x-0.5*dx for x, dx in zip(xyzOrigin, xyzStep)])

        # {cell-id (i,j,k) : [list of indexes in points]}
        cellIdPtId = self.grid.intersectionPointList(points)
        if any(len(pts)>1 for pts in cellIdPtId.itervalues()):
            raise ValueError("Found multiple points in the same grid cell.")
        assert not(any(len(pts)<1 for pts in cellIdPtId.itervalues())), (
            "This should just not happen! Well it did apparently. What"
            " now? Look at the code. And search for the error.")
        self.cellIdVolId = dict((k, ids[v[0]])
                                for k,v in cellIdPtId.iteritems())
        del cellIdPtId

        self.outsideId = outsideId

    def pointIsInsideParts(self, point):
        """
        Check in which subvolumes (parts) the specified point is located

        @param point: triple of coordinates

        @returns: A list containing one(!) subvolume identifier corresponding
        to the grid cell the point lies in. If the point does not lie in any
        grid cell with a subvolume identifier assigned to it then an empty
        list is being returned.

        Note: This method returns a list in any case and not just the single
        value in order to be compatible with the corresponding methods of other
        MultiVol-classes. The alternative L{pointIsInsideOnePart} method can be
        used if you want just one subvolume identifier.
        """
        gridIdx = self.grid.pointIsInsideOnePart(point)
        try:
            volId = self.cellIdVolId[gridIdx]
        except KeyError:
            return []
        else:
            return [volId,]

    def pointIsInsideOnePart(self, point):
        """
        Check in which subvolume (part) the specified point is located.
        This is special to Multivol classes that have mutually disjoint
        subvolumes and is equivalent to self.pointIsInsideParts(point)[0]
        (if the point is in any of the subvolumes at all).

        See also L{pointIsInsideParts}.

        @param point: triple of coordinates

        @returns: The subvolume identifier assigned to the grid cell the point
        lies in if any. Or outsideId as spicified by the corresponding argument
        to self.__init__() if not.
        """
        gridIdx = self.grid.pointIsInsideOnePart(point)
        return self.cellIdVolId.get(gridIdx, self.outsideId)


class MultiVolClosestPoint(MultiVolUnique):
    """
    Volume parts defined by point data. Each of those constitutive points
    carries a subvolume identifier.

    An arbitrary query point belongs to the subvolume identified by the
    closest of the constitutive points.
    """

    def __init__(self, points, ids):
        """
        @param points: list of coordinate tuples, i.e. [[x1,y1,z1], [x2,y2,z2],
          ...]
        @param ids: list of subvolume identifiers.
        """
        if len(points)!=len(ids):
            raise ValueError(
                "points and ids arguments don't appear to have the same"
                " length. points: %d, ids: %d" % ( len(points), len(ids)))

        MultiVolUnique.__init__(self)
        if np and hasattr(scipy.spatial, "cKDTree"):
            self.kdt = scipy.spatial.cKDTree(points)
            self.ids = np.array(ids)
        else:
            self.kdt = KDTree(points)
            self.ids = ids

    def _nonp_pointIsInsideOnePart(self, point):
        idx = self.kdt.knnSearch(point, K=1)[0][1]
        return self.ids[idx]
    def _np_pointIsInsideOnePart(self, point):
        try:
            dists, idx = self.kdt.query(point, k=1, eps=1e-3, n_jobs=-1)
        except TypeError:
            # old versions of scipy don't have the n_jobs argument!
            dists, idx = self.kdt.query(point, k=1, eps=1e-3)
        return self.ids[idx]
    pointIsInsideOnePart = selectattrversion(
        np and hasattr(scipy.spatial, "cKDTree"),
        _np_pointIsInsideOnePart,
        _nonp_pointIsInsideOnePart, """
        Check in which subvolume (part) the specified point is located.
        This is special to Multivol classes that have mutually disjoint
        subvolumes and is equivalent to self.pointIsInsideParts(point)[0]

        @param point: triple of coordinates

        @returns: A tuple of three integer indexes identifying the grid box
        the point lies in.
        """)

    def _nonp_pointIsInsideParts(self, point):
        idx = self.kdt.knnSearch(point, K=1)[0][1]
        return [self.ids[idx],]
    def _np_pointIsInsideParts(self, point):
        try:
            dists, idx = self.kdt.query(point, k=1, eps=1e-3, n_jobs=-1)
        except TypeError:
            # old versions of scipy don't have the n_jobs argument!
            dists, idx = self.kdt.query(point, k=1, eps=1e-3)
        return [self.ids[idx],]
    pointIsInsideParts = selectattrversion(
        np and hasattr(scipy.spatial, "cKDTree"),
        _np_pointIsInsideParts,
        _nonp_pointIsInsideParts, """
        Check in which subvolumes (parts) the specified point is located

        @param point: triple of coordinates

        @returns: A list containing one(!) subvolume identifier corresponding
           to the closest constitutive point.

           Why is it a list and not just the single value? To make it compatible
           with the corresponding methods of other MultiVol-classes.
        """)

    if np and hasattr(scipy.spatial, "cKDTree"):
        def intersectionPointIdTuples(self, pointIdTuples):
            r"""
            Assign each point with its corresponding identifier to the volume
            parts it is in. Takes a list of (point, point identifier) tuples and
            returns a dictionary {subvol id: set of point ids}.

            This is a numpy/scipy enhanced variant.
            See L{MultiVolUnique.intersectionPointIdTuples} for description.
            """
            points = np.array([xyz for xyz, id_ in pointIdTuples])
            try:
                dists, idsInSelfIds = self.kdt.query(
                    points, k=1, eps=1e-3, n_jobs=-1)
            except TypeError:
                # old versions of scipy don't have the n_jobs argument!
                dists, idsInSelfIds = self.kdt.query(points, k=1, eps=1e-3)
            volParts = self.ids[idsInSelfIds]

            volPartToId = defaultdict(set)
            for volPart, (pt, ptId) in izip(volParts, pointIdTuples):
                volPartToId[volPart].add(ptId)
            volPartToId.default_factory = None  # convert to an ordinary dict
            return volPartToId
