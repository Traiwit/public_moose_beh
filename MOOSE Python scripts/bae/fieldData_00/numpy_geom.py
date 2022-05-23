# -*- coding: utf-8 -*-
"""
Collection of classes and functions for solving geometric problems
unsing numpy.
"""
import numpy as np
import sys, os

#check module is loaded by epydoc-build
#--> epyDocBuild=True to disable decorators
if 'epydoc' in sys.modules:
    _epyDocBuild = True
else:
    _epyDocBuild = False

### import scipy cKDTree
try:
    from scipy.spatial import cKDTree
except ImportError:
    cKDTree = None

### import rtree
try:
    from rtree import index, Rtree
    import cPickle

    class FastRtree(Rtree):
        """Faster R-TreeVersion - documented in Rtree manuel"""
        def dumps(self, obj):
            return cPickle.dumps(obj,-1)
except:
    FastRtree = None

### import numba
try:
    from numba import njit, jit
except:
    def jit(func):
        # Just Do nothing. Decorated function might be very slow.
        return func

    def njit(func):
        # Just Do nothing. Decorated function might be very slow.
        return func

from bae.log_01 import msg

#### Decorator functions ######################################################
#{decorator functions
#adds the additional keyword-argument 'asRad' to function, which converts
#angle-like inputData from rad to deg and outputData from deg to rad
#respectively
def _inAsRad(func):
    """Adds asRad to kwarg. If asRad=True the input values will be converted
    from rad to degree before they get passed to decorated function.
    """
    #to get correct args/kwargs in epydoc
    if _epyDocBuild:
        return func

    #decorator
    def func_wrapper(*args, **kwargs):
        if 'asRad' in kwargs:
            if kwargs['asRad']:
                args = [arg*(180./np.pi) for arg in args]
                del kwargs['asRad']
        return func(*args, **kwargs)
    return func_wrapper


def _outAsRad(func):
    """Adds asRad to kwarg. If asRad=True the returned values will be converted
    from degree to rad.
    """
    #to get correct args/kwargs in epydoc
    if _epyDocBuild:
        return func

    #decorator
    def func_wrapper(*args, **kwargs):
        transform = False
        if 'asRad' in kwargs:
            if kwargs['asRad']:
                del kwargs['asRad']
                transform = True
        if transform:
            return (out*np.pi/180. for out in func(*args,**kwargs))
        else:
            return func(*args,**kwargs)
    return func_wrapper


def _inAndOutAsRad(func):
    """Adds asRad to kwarg. If asRad=True the input values will be converted
    from rad to degree before they get passed to decorated function. And
    output values will get converted from deg to rad.
    """
    if _epyDocBuild:
        return func

    #decorator
    def func_wrapper(*args, **kwargs):
        transform = False
        if 'asRad' in kwargs:
            if kwargs['asRad']:
                del kwargs['asRad']
                transform = True
                args = [arg*(180/np.pi) for arg in args]
        if transform:
            return (out*np.pi/180. for out in func(*args,**kwargs))
        else:
            return func(*args,**kwargs)
    return func_wrapper

#} #end decorator functions
###############################################################################


###############################################################################
#{ mplstereonet
"""Utilities to convert between strike/dip, etc and points/lines in lat, long
space.

A stereonet in <long,lat> coordinates::
              <0,90>
               ***
            *       *
   <-90,0> *         *<90,0>
           *         *
            *       *
               ***
             <0,-90>

If strike=0, plotting lines, rakes, planes or poles to planes is simple.  For a
plane, it's a line of constant longitude at long=90-dip.  For a line, it's a
point at long=0,lat=90-dip.  For a rake, it's a point at long=90-dip,
lat=90-rake.  These points can then be rotated to the proper strike. (A
rotation matrix around the X-axis is much simpler than the trig otherwise
necessary!)
All of these assume that strikes and dips follow the "right-hand-rule".
In other words, if we're facing in the direction given for the strike, the plane
dips to our right.

Taken from:
U{https://github.com/joferkington/mplstereonet/blob/master/mplstereonet}
and added decorators to use radians
"""

@_inAsRad
def sph2cart(lon, lat):
    """
    Converts a longitude and latitude (or sequence of lons and lats) given in
    _radians_ to cartesian coordinates, 'x', 'y', 'z', where x=0, y=0, z=0 is
    the center of the globe.

    @param lon: array-like, longitude in radians

    @param lat: array-like, latitude in radians

    @return: 'x', 'y', 'z': Arrays of cartesian coordinates
    """
    x = np.cos(lat)*np.cos(lon)
    y = np.cos(lat)*np.sin(lon)
    z = np.sin(lat)
    return x, y, z


@_outAsRad
def cart2sph(x, y, z):
    """
    Converts cartesian coordinates 'x', 'y', 'z' into a longitude and latitude.
    x=0, y=0, z=0 is assumed to correspond to the center of the globe.
    Returns lon and lat in radians.

    @param x: array of x-cartesian coordinates

    @param y: array of y-cartesian coordinates

    @param z: array of z-cartesian coordinates

    @return: (longitude, latitude) both in radians
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    lat = np.arcsin(z/r)
    lon = np.arctan2(y, x)
    return lon, lat


def _rotate(lon, lat, theta, axis='x'):
    """
    Rotate "lon", "lat" coords (in _degrees_) about the X-axis by "theta"
    degrees.  This effectively simulates rotating a physical stereonet.

    @Returns: rotated lon, lat coords in _radians_.
    """
    # Convert input to numpy arrays in radians
    lon, lat = np.atleast_1d(lon, lat)
    lon, lat = map(np.radians, [lon, lat])
    theta = np.radians(theta)

    # Convert to cartesian coords for the rotation
    x, y, z = sph2cart(lon, lat)

    lookup = {'x':_rotate_x, 'y':_rotate_y, 'z':_rotate_z}
    X, Y, Z = lookup[axis](x, y, z, theta)

    # Now convert back to spherical coords (longitude and latitude, ignore R)
    lon, lat = cart2sph(X,Y,Z)
    return lon, lat  # in radians!


def _rotate_x(x, y, z, theta):
    X = x
    Y = y*np.cos(theta) + z*np.sin(theta)
    Z = -y*np.sin(theta) + z*np.cos(theta)
    return X, Y, Z


def _rotate_y(x, y, z, theta):
    X = x*np.cos(theta) + -z*np.sin(theta)
    Y = y
    Z = x*np.sin(theta) + z*np.cos(theta)
    return X, Y, Z


def _rotate_z(x, y, z, theta):
    X = x*np.cos(theta) + -y*np.sin(theta)
    Y = x*np.sin(theta) + y*np.cos(theta)
    Z = z
    return X, Y, Z


@_inAsRad
def antipode(lon, lat):
    """
    Calculates the antipode (opposite point on the globe) of the given point or
    points. Input and output is expected to be in radians.

    @param lon: number or sequence of numbers, Longitude in radians

    @param lat: number or sequence of numbers, Latitude in radians

    @return: lon, lat as arrays
        Sequences (regardless of whether or not the input was a single value or
        a sequence) of longitude and latitude in radians.
    """
    x, y, z = sph2cart(lon, lat)
    return cart2sph(-x, -y, -z)


@_inAndOutAsRad
def plane(strike, dip, segments=100, center=(0, 0)):
    """
    Calculates the longitude and latitude of 'segments' points along the
    stereonet projection of each plane with a given 'strike' and 'dip' in
    degrees.  Returns points for one hemisphere only.

    @param strike: number or sequence of numbers
        The strike of the plane(s) in degrees, with dip direction indicated by
        the azimuth (e.g. 315 vs. 135) specified following the "right hand
        rule".

    @param dip: number or sequence of numbers
        The dip of the plane(s) in degrees.

    @param segments: number or sequence of numbers
        The number of points in the returned 'lon' and 'lat' arrays.  Defaults
        to 100 segments.

    @param center: sequence of two numbers (lon, lat)
        The longitude and latitude of the center of the hemisphere that the
        returned points will be in. Defaults to 0,0 (approriate for a typical
        stereonet).

    @return: lon, lat as arrays
        'num_segments' x 'num_strikes' arrays of longitude and latitude in
        radians.
    """
    lon0, lat0 = center
    strikes, dips = np.atleast_1d(strike, dip)
    lons = np.zeros((segments, strikes.size), dtype=np.float)
    lats = lons.copy()
    for i, (strike, dip) in enumerate(zip(strikes, dips)):
        # We just plot a line of constant longitude and rotate it by the strike.
        dip = 90 - dip
        lon = dip * np.ones(segments)
        lat = np.linspace(-90, 90, segments)
        lon, lat = _rotate(lon, lat, strike)

        if lat0 != 0 or lon0 != 0:
            dist = angular_distance([lon, lat], [lon0, lat0], False)
            mask = dist > (np.pi / 2)
            lon[mask], lat[mask] = antipode(lon[mask], lat[mask])
            change = np.diff(mask.astype(int))
            ind = np.flatnonzero(change) + 1
            lat = np.hstack(np.split(lat, ind)[::-1])
            lon = np.hstack(np.split(lon, ind)[::-1])

        lons[:,i] = lon
        lats[:,i] = lat

    return lons, lats


@_inAndOutAsRad
def pole(strike, dip):
    """
    Calculates the longitude and latitude of the pole(s) to the plane(s)
    specified by 'strike' and 'dip', given in degrees.

    @param strike: number or sequence of numbers
        The strike of the plane(s) in degrees, with dip direction indicated by
        the azimuth (e.g. 315 vs. 135) specified following the "right hand
        rule".

    @param dip: number or sequence of numbers
        The dip of the plane(s) in degrees.

    @return: lon, lat: Arrays of longitude and latitude in radians.
    """
    strike, dip = np.atleast_1d(strike, dip)
    mask = dip > 90
    dip[mask] = 180 - dip[mask]
    strike[mask] += 180
    # Plot the approriate point for a strike of 0 and rotate it
    lon, lat = -dip, 0.0
    lon, lat = _rotate(lon, lat, strike)
    return lon, lat


@_inAndOutAsRad
def rake(strike, dip, rake_angle):
    """
    Calculates the longitude and latitude of the linear feature(s) specified by
    'strike', 'dip', and 'rake_angle'.

    @param strike: number or sequence of numbers
        The strike of the plane(s) in degrees, with dip direction indicated by
        the azimuth (e.g. 315 vs. 135) specified following the "right hand
        rule".

    @param dip: number or sequence of numbers
        The dip of the plane(s) in degrees.

    @param rake_angle: number or sequence of numbers
        The angle of the lineation on the plane measured in degrees downward
        from horizontal. Zero degrees corresponds to the "right- hand"
        direction indicated by the strike, while 180 degrees or a negative
        angle corresponds to the opposite direction.

    @return: lon, lat: Arrays of longitude and latitude in radians.
    """
    strike, dip, rake_angle = np.atleast_1d(strike, dip, rake_angle)
    # Plot the approriate point for a strike of 0 and rotate it
    dip = 90 - dip
    lon = dip

    rake_angle = rake_angle.copy()
    rake_angle[rake_angle < 0] += 180
    lat = 90 - rake_angle

    lon, lat = _rotate(lon, lat, strike)
    return lon, lat


@_inAndOutAsRad
def line(plunge, bearing):
    """
    Calculates the longitude and latitude of the linear feature(s) specified by
    'plunge' and 'bearing'.

    @param plunge: number or sequence of numbers
        The plunge of the line(s) in degrees. The plunge is measured in degrees
        downward from the end of the feature specified by the bearing.

    @param bearing: number or sequence of numbers
        The bearing (azimuth) of the line(s) in degrees.

    @return: lon, lat: Arrays of longitude and latitude in radians.
    """
    plunge, bearing = np.atleast_1d(plunge, bearing)
    # Plot the approriate point for a bearing of 0 and rotate it
    lat = 90 - plunge
    lon = 0
    lon, lat = _rotate(lon, lat, bearing)
    return lon, lat


@_inAndOutAsRad
def cone(plunge, bearing, angle, segments=100):
    """
    Calculates the longitude and latitude of the small circle (i.e. a cone)
    centered at the given *plunge* and *bearing* with an apical angle of
    *angle*, all in degrees.

    @param plunge: number or sequence of numbers
        The plunge of the center of the cone(s) in degrees. The plunge is
        measured in degrees downward from the end of the feature specified by
        the bearing.

    @param bearing: number or sequence of numbers
        The bearing (azimuth) of the center of the cone(s) in degrees.

    @param angle: number or sequence of numbers
        The apical angle (i.e. radius) of the cone(s) in degrees.

    @param segments: int, optional
        The number of vertices in the small circle.

    @return: lon, lat: arrays
        'num_measurements' x 'num_segments' arrays of longitude and latitude in
        radians.
    """
    plunges, bearings, angles = np.atleast_1d(plunge, bearing, angle)
    lons, lats = [], []
    for plunge, bearing, angle in zip(plunges, bearings, angles):
        lat = (90 - angle) * np.ones(segments, dtype=float)
        lon = np.linspace(-180, 180, segments)
        lon, lat = _rotate(lon, lat, -plunge, axis='y')
        lon, lat = _rotate(np.degrees(lon), np.degrees(lat), bearing, axis='x')
        lons.append(lon)
        lats.append(lat)
    return np.vstack(lons), np.vstack(lats)


@_inAndOutAsRad
def plunge_bearing2pole(plunge, bearing):
    """
    Converts the given 'plunge' and 'bearing' in degrees to a strike and dip
    of the plane whose pole would be parallel to the line specified. (i.e. The
    pole to the plane returned would plot at the same point as the specified
    plunge and bearing.)

    @param plunge: number or sequence of numbers
        The plunge of the line(s) in degrees. The plunge is measured in degrees
        downward from the end of the feature specified by the bearing.

    @param bearing: number or sequence of numbers
        The bearing (azimuth) of the line(s) in degrees.

    @return: strike, dip: arrays
        Arrays of strikes and dips in degrees following the right-hand-rule.
    """
    plunge, bearing = np.atleast_1d(plunge, bearing)
    strike = bearing + 90
    dip = 90 - plunge
    strike[strike >= 360] -= 360
    return strike, dip


@_inAndOutAsRad
def pole2plunge_bearing(strike, dip):
    """
    Converts the given *strike* and *dip* in dgrees of a plane(s) to a plunge
    and bearing of its pole.

    @param strike: number or sequence of numbers
        The strike of the plane(s) in degrees, with dip direction indicated by
        the azimuth (e.g. 315 vs. 135) specified following the "right hand
        rule".

    @param dip: number or sequence of numbers
        The dip of the plane(s) in degrees.

    @return: plunge, bearing: arrays
        Arrays of plunges and bearings of the pole to the plane(s) in degrees.
    """
    strike, dip = np.atleast_1d(strike, dip)
    bearing = strike - 90
    plunge = 90 - dip
    bearing[bearing < 0] += 360
    return plunge, bearing


@_inAndOutAsRad
def mean_vector(lons, lats):
    """Returns the resultant vector from a series of longitudes and latitudes.

    @param lons: array-like, A sequence of longitudes (in radians)

    @param lats: array-like, A sequence of latitudes (in radians)

    @return: (mean_vec, r_value) - tuple:

    mean_vec: (lon, lat) - tuple in radians.

    r_value:
    The magnitude of the resultant vector (between 0 and 1) This
    represents the degree of clustering in the data.
    """
    coords = sph2cart(lons, lats)
    coords = np.vstack(coords).T
    mean_vec = coords.mean(axis=0)
    r_value = np.linalg.norm(mean_vec)
    mean_vec = cart2sph(*mean_vec)
    return mean_vec, r_value


@_inAndOutAsRad
def fisher_stats(lons, lats, conf=95):
    """
    Returns the resultant vector from a series of longitudes and latitudes. If
    a confidence is set the function additionally returns the opening angle
    of the confidence small circle (Fisher, 19..) and the dispersion factor
    (kappa).

    @param lons: array-like
        A sequence of longitudes (in radians)

    @param lats: array-like
        A sequence of latitudes (in radians)

    @param conf: confidence value
        The confidence used for the calculation (float). Defaults to None.

    @return: mean vector, a tuple. The point that lies in the center of a set
    of vectors. (Longitude, Latitude) in radians.

    If one vector is passed to the function it returns two None-values. For
    more than one vector the following 3 values are returned as a tuple:

    r_value (float):
    The magnitude of the resultant vector (between 0 and 1) This represents
    the degree of clustering in the data.

    angle (float):
    The opening angle of the small circle that corresponds to confidence
    of the calculated direction.

    kappa (float):
    A measure for the amount of dispersion of a group of layers. For
    one vector the factor is undefined. Approaches infinity for nearly
    parallel vectors and zero for highly dispersed vectors.
    """
    coords = sph2cart(lons, lats)
    coords = np.vstack(coords).T
    mean_vec = coords.mean(axis=0)
    r_value = np.linalg.norm(mean_vec)
    num = coords.shape[0]
    mean_vec = cart2sph(*mean_vec)

    if num > 1:
        p = (100.0 - conf) / 100.0
        vector_sum = coords.sum(axis=0)
        result_vect = np.sqrt(np.sum(np.square(vector_sum)))
        fract1 = (num - result_vect) / result_vect
        fract3 = 1.0 / (num - 1.0)
        angle = np.arccos(1 - fract1 * ((1 / p) ** fract3 - 1))
        angle = np.degrees(angle)
        kappa = (num - 1.0) / (num - result_vect)
        return mean_vec, (r_value, angle, kappa)
    else:
        return None, None


@_inAndOutAsRad
def geographic2pole(lon, lat):
    """
    Converts a longitude and latitude (from a stereonet) into the strike and dip
    of the plane whose pole lies at the given longitude(s) and latitude(s).

    @param lon: array-like
        A sequence of longitudes (or a single longitude) in radians

    @param lat: array-like
        A sequence of latitudes (or a single latitude) in radians

    @return: (strike, dip) -tuple of arrays.
    strike being a sequence of strikes in degrees.
    dip being a sequence of dips in degrees.
    """
    plunge, bearing = geographic2plunge_bearing(lon, lat)
    strike = bearing + 90
    strike[strike >= 360] -= 360
    dip = 90 - plunge
    return strike, dip


@_inAndOutAsRad
def geographic2plunge_bearing(lon, lat):
    """
    Converts longitude and latitude in stereonet coordinates into a
    plunge/bearing.

    @param lon: Longitudes in radians as measured from a lower-hemisphere

    @param lat: Longitudes in radians as measured from a lower-hemisphere

    @return: (plunge, bearing) -tuple of arrays holding the plunge of the
        vector in degrees downward from horizontal. And the bearing of the
        vector in degrees clockwise from north.
    """
    lon, lat = np.atleast_1d(lon, lat)
    x, y, z = sph2cart(lon, lat)

    # Bearing will be in the y-z plane...
    bearing = np.arctan2(z, y)

    # Plunge is the angle between the line and the y-z plane
    r = np.sqrt(x*x + y*y + z*z)
    r[r == 0] = 1e-15
    plunge = np.arcsin(x / r)

    # Convert back to azimuths in degrees..
    plunge, bearing = np.degrees(plunge), np.degrees(bearing)
    bearing = 90 - bearing
    bearing[bearing < 0] += 360

    # If the plunge angle is upwards, get the opposite end of the line
    upwards = plunge < 0
    plunge[upwards] *= -1
    bearing[upwards] -= 180
    bearing[upwards & (bearing < 0)] += 360

    return plunge, bearing


@_inAndOutAsRad
def plane_intersection(strike1, dip1, strike2, dip2):
    """
    Finds the intersection of two planes. Returns a plunge/bearing of the linear
    intersection of the two planes.
    Also accepts sequences of strike1s, dip1s, strike2s, dip2s.

    @param strike1, dip1: numbers or sequences of numbers
        The strike and dip (in degrees, following the right-hand-rule) of the
        first plane(s).

    @param strike2, dip2: numbers or sequences of numbers
        The strike and dip (in degrees, following the right-hand-rule) of the
        second plane(s).

    @return: plunge, bearing: arrays
        The plunge and bearing(s) (in degrees) of the line representing the
        intersection of the two planes.
    """
    norm1 = sph2cart(*pole(strike1, dip1))
    norm2 = sph2cart(*pole(strike2, dip2))
    norm1, norm2 = np.array(norm1), np.array(norm2)
    lon, lat = cart2sph(*np.cross(norm1, norm2, axis=0))
    return geographic2plunge_bearing(lon, lat)


@_inAndOutAsRad
def project_onto_plane(strike, dip, plunge, bearing):
    """
    Projects a linear feature(s) onto the surface of a plane. Returns a rake
    angle(s) along the plane.
    This is also useful for finding the rake angle of a feature that already
    intersects the plane in question.

    @param strike, dip: numbers or sequences of numbers
        The strike and dip (in degrees, following the right-hand-rule) of the
        plane(s).

    @param plunge, bearing: numbers or sequences of numbers
        The plunge and bearing (in degrees) or of the linear feature(s) to be
        projected onto the plane.

    @return: rake: array
        A sequence of rake angles measured downwards from horizontal in
        degrees.  Zero degrees corresponds to the "right- hand" direction
        indicated by the strike, while a negative angle corresponds to the
        opposite direction. Rakes returned by this function will always be
        between -90 and 90 (inclusive).
    """
    # Project the line onto the plane
    norm = sph2cart(*pole(strike, dip))
    feature = sph2cart(*line(plunge, bearing))
    norm, feature = np.array(norm), np.array(feature)
    perp = np.cross(norm, feature, axis=0)
    on_plane = np.cross(perp, norm, axis=0)
    on_plane /= np.sqrt(np.sum(on_plane**2, axis=0))

    # Calculate the angle between the projected feature and horizontal
    # This is just a dot product, but we need to work with multiple measurements
    # at once, so einsum is quicker than apply_along_axis.
    strike_vec = sph2cart(*line(0, strike))
    dot = np.einsum('ij,ij->j', on_plane, strike_vec)
    rake = np.degrees(np.arccos(dot))

    # Convert rakes over 90 to negative rakes...
    rake[rake > 90] -= 180
    rake[rake < -90] += 180
    return rake


@_inAndOutAsRad
def azimuth2rake(strike, dip, azimuth):
    """
    Projects an azimuth of a linear feature onto a plane as a rake angle.

    @param strike: ... of the plane in degrees

    @param dip: ... of the plane in degrees.
        Strike and dip follow the right-hand-rule.

    @param azimuth: numbers
        The azimuth of the linear feature in degrees clockwise from north (i.e.
        a 0-360 azimuth).

    @return: rake: number
        A rake angle in degrees measured downwards from horizontal. Negative
        values correspond to the opposite end of the strike.
    """
    plunge, bearing = plane_intersection(strike, dip, azimuth, 90)
    rake = project_onto_plane(strike, dip, plunge, bearing)
    return rake


@_outAsRad
def coords2stereonet(x, y, z):
    """
    Converts x, y, z in _world_ cartesian coordinates into lower-hemisphere
    stereonet coordinates.

    @param x, y, z: array-likes
        Sequences of world coordinates

    @return: lon, lat: arrays
        Sequences of longitudes and latitudes (in radians)
    """
    x, y, z = np.atleast_1d(x, y, z)
    return cart2sph(-z, x, y)


@_inAsRad
def stereonet2coords(lon, lat):
    """
    Converts a sequence of longitudes and latitudes from a lower-hemisphere
    stereonet into _world_ x,y,z coordinates.

    @param lon, lat: array-likes
        Sequences of longitudes and latitudes (in radians) from a
        lower-hemisphere stereonet

    @return: x, y, z: arrays
        The world x,y,z components of the vectors represented by the lon, lat
        coordinates on the stereonet.
    """
    lon, lat = np.atleast_1d(lon, lat)
    x, y, z = sph2cart(lon, lat)
    return y, z, -x


@_outAsRad
def vector2plunge_bearing(x, y, z):
    """
    Converts a vector or series of vectors given as x, y, z in world
    coordinates into plunge/bearings.

    @param x: number or sequence of numbers
        The x-component(s) of the normal vector

    @param y: number or sequence of numbers
        The y-component(s) of the normal vector

    @param z: number or sequence of numbers
        The z-component(s) of the normal vector

    @return: (plungs, bearing) tuple of arrays
        - plunge: The plunge of the vector in degrees downward from horizontal.
        - bearing: The bearing of the vector in degrees clockwise from north.
    """
    return geographic2plunge_bearing(*coords2stereonet(x,y,z))


@_outAsRad
def vector2pole(x, y, z):
    """
    Converts a vector or series of vectors given as x, y, z in world
    coordinates into the strike/dip of the planes whose normal vectors are
    parallel to the specified vectors.  (In other words, each xi,yi,zi is
    treated as a normal vector and this returns the strike/dip of the
    corresponding plane.)

    @param x: number or sequence of numbers
        The x-component(s) of the normal vector

    @param y: number or sequence of numbers
        The y-component(s) of the normal vector

    @param z: number or sequence of numbers
        The z-component(s) of the normal vector

    @return: (strike, dip) tuple of arrays
        - strike: The strike of the plane, in degrees clockwise from north.  Dip
        direction is indicated by the "right hand rule".
        -dip: The dip of the plane, in degrees downward from horizontal.
    """
    return geographic2pole(*coords2stereonet(x, y, z))


@_inAndOutAsRad
def angular_distance(first, second, bidirectional=True):
    """
    Calculate the angular distance between two linear features or elementwise
    angular distance between two sets of linear features. (Note: a linear
    feature in this context is a point on a stereonet represented
    by a single latitude and longitude.)


    Examples
    ========

    Calculate the angle between two lines specified as a plunge/bearing
        >>> angle = angular_distance(line(30, 270), line(40, 90))
        >>> np.degrees(angle)
        array([ 70.])

    Let's do the same, but change the "bidirectional" argument:
        >>> first, second = line(30, 270), line(40, 90)
        >>> angle = angular_distance(first, second, bidirectional=False)
        >>> np.degrees(angle)
        array([ 110.])

    Calculate the angle between two planes.
        >>> angle = angular_distance(pole(0, 10), pole(180, 10))
        >>> np.degrees(angle)
        array([ 20.])

    @param first: (lon, lat) 2xN array-like or sequence of two numbers
        The longitudes and latitudes of the first measurements in radians.

    @param second: (lon, lat) 2xN array-like or sequence of two numbers
        The longitudes and latitudes of the second measurements in radians.

    @param bidirectional: boolean
        If True, only "inner" angles will be returned. In other words, all
        angles returned by this function will be in the range [0, pi/2]
        (0 to 90 in degrees).  Otherwise, 'first' and 'second'
        will be treated as vectors going from the origin outwards
        instead of bidirectional infinite lines.  Therefore, with
        'bidirectional=False', angles returned by this function
        will be in the range [0, pi] (zero to 180 degrees).

    @return: dist: array
        The elementwise angular distance between each pair of measurements in
        (lon1, lat1) and (lon2, lat2).

    """
    lon1, lat1 = first
    lon2, lat2 = second
    lon1, lat1, lon2, lat2 = np.atleast_1d(lon1, lat1, lon2, lat2)
    coords1 = sph2cart(lon1, lat1)
    coords2 = sph2cart(lon2, lat2)
    # This is just a dot product, but we need to work with multiple
    #  measurements at once, so einsum is quicker than apply_along_axis.
    dot = np.einsum('ij,ij->j', coords1, coords2)
    angle = np.arccos(dot)

    # There are numerical sensitivity issues around 180 and 0 degrees...
    # Sometimes a result will have an absolute value slighly over 1.
    if np.any(np.isnan(angle)):
        rtol = 1e-4
        angle[np.isclose(dot, -1, rtol)] = np.pi
        angle[np.isclose(dot, 1, rtol)] = 0

    if bidirectional:
        mask = angle > np.pi / 2
        angle[mask] = np.pi - angle[mask]

    return angle

#} #end mplstereonet


###############################################################################
### functions to calculate and transform eigenvalues and directions
#{ eigenvalues and eigendirections
def eigen(dataTensor):
    '''Calculate local, framewise principle stresses and their directions.
    They are ordered according to the largest eigenvalue.

    @param dataTensor: numpy array of dimension [nPoints,nFrames,3,3]

    @return: S = numpy array holding ordered (largest first) eigenvalues
              (dim = [nPoints,nFrames,3])
              V = eigenvectors (dim = [nPoints,nFrames,3,3])
    '''
    #calculating eigenvalues and eigenvectors
    #msg("Computing eigenvalues and -vectors")
    S, V = np.linalg.eigh(dataTensor)
    # remember: eig retuns eigenvector where last Index == Index eigenvalue
    #           so _,ev=eigh(mat) --> ev.T.dot(mat).dot(ev) = I
    # swap dims to be consistent -- > last Index == coords:
    V = np.swapaxes(V, -1, -2)

    #reversed order: largest eigval first
    def apply_argsort(a, axis=-1, reverse=False):
        # numpy.argsort is a bit tricky for nD-arrays
        i = list(np.ogrid[[slice(x) for x in a.shape]])
        if reverse:
            a = -a
            i[axis] = a.argsort(axis)
        else:
            i[axis] = a.argsort(axis)
        return tuple(i)

    oIdx = apply_argsort(S, axis=-1, reverse=True)
    S = S[oIdx]
    V = V[oIdx]

    return S, V


def transformEigen(S, V,
                   eigenFmt='polar+azimut',
                   angleFmt='rad'):
    '''
    @param S: numpy-array containing eigenvalues (dim: [nPoints,nFrames,3])

    @param V: numpy-array containing S-related eigenvectors
        (dim: [nPoints,nFrames,3,3])

    @param eigFmt: format of eigenvectors:
        'polar+azimut', 'strike+dip' or 'bearing+plunge'

    @param angleFmt: output angles in degree ('deg') or radians ('rad')

    @return: eigenValues and directions
         [nPoints,nFrames,3, [ magnitude, 1stDirAngle, 2ndDirAngle ] ]

    '''
    # switch outputformat: degree or radian
    if angleFmt == 'deg':
        fac = 180./np.pi
    else:
        fac = 1

    [nPoints, nFrames, _] = S.shape
    Strans = np.zeros((nPoints, nFrames, 3, 3))

    # first transform to spherical coordinates:
    # phi   = rotation of x-axis in direction of/around (+) z-axis
    #       phi = arctan2(ny,nx)
    #           = pi/2 - alpha_strike
    #           = pi - alpha_trend
    #           = pi - alpha_bearing
    # theta = rotation in direction/around rotated y-axis = alpha_dip
    #       theta = arccos(nz)-pi/2

    # force vectors to point downwards pi/2<theta<pi --> rot phi by 180deg
    #idx        = theta<np.pi/2
    #phi[idx]   = phi[idx] + np.pi
    #theta[idx] = np.pi/2 - theta[idx]

    # transform directions to desired form
    if eigenFmt.lower() == 'polar+azimut':
        R = np.linalg.norm(V, ord=2, axis=-1)
        R[R==0] = 1E-12
        phi   = np.arctan2(V[:,:,:,1], V[:,:,:,0])
        theta = np.arccos(V[:,:,:,2]/R)
    elif eigenFmt.lower() == 'strike+dip':
        #phi   = np.mod(.5*np.pi - phi,2*np.pi)*fac  #strike
        #theta = ( theta-0.5*np.pi)*fac              #dip
        phi,theta = vector2plunge_bearing(V[:,:,:,0], V[:,:,:,1],
                                          V[:,:,:,2], asRad=True)
        theta = theta+np.pi/2
        #phi,theta = plunge_bearing2pole(phi,theta,asRad=True)
        # phi=strike, theta=dip
    elif eigenFmt.lower() == 'bearing+plunge':
        #phi   =   phi - np.pi       #bearing
        #theta = (theta-0.5*np.pi)   #plunge=dip
        #idx = theta<0
        #phi[idx] = np.mod(phi[idx]-np.pi,2*np.pi)
        #theta[idx]=-theta[idx]
        phi,theta = vector2plunge_bearing(V[:,:,:,0], V[:,:,:,1],
                                          V[:,:,:,2], asRad=True)
        # phi=plunge,theta=bearing #as returned from vector2plunge_bearing
    else:
        raise Exception('%s is not a valid id for eigenvalue fmt' % eigenFmt)

    #stack data for output [:points,:frames,:eig123,[mag,dir1,dir2]]
    Strans[:,:,:,0] = S
    Strans[:,:,:,1] = theta * fac
    Strans[:,:,:,2] = phi * fac

    return Strans
#} #end eigenvalues and eigendirections

###############################################################################
### Statistical functions for angular fields
#{Statistics for angular fields
def meanAngle(angle, axis=-1, igNan=False, isDeg=False):
    '''Calculates the average angle for a set of angles using complex pointer.

    np.mean([2*np.pi,0]) = pi but meanAngle([2*np.pi,0]) = 0

    @param angle : array of angles

    @param axis  : dimension to apply mean (default = last)

    @param igNan : True --> ignores nan-values

    @param isDeg : True --> input/output-values in degrees
    '''
    if isDeg:
        angle = np.radians(angle)

    if not igNan:
        avrg = np.arctan2(np.mean(np.sin(angle), axis=axis),
                          np.mean(np.cos(angle), axis=axis))
    else:
        avrg = np.arctan2(np.nanmean(np.sin(angle), axis=axis),
                          np.nanmean(np.cos(angle), axis=axis))

    if isDeg:
        avrg = np.degrees(avrg)

    return avrg


def stdAngle(angle, axis=-1, igNan=False, isDeg=False):
    '''Calculates the standard deviation for a set of angles using complex
    pointer. See also meanAngle.

    @param angle : array of angles

    @param axis  : dimension to apply mean (default = last)

    @param igNan : True --> ignores nan-values

    @param isDeg : True --> input/output-values in degrees
    '''
    if isDeg:
        angle = np.radians(angle)

    if not igNan:
        std = np.arctan2(np.std(np.sin(angle), axis=axis),
                         np.std(np.cos(angle), axis=axis))
    else:
        std = np.arctan2(np.nanstd(np.sin(angle), axis=axis),
                         np.nanstd(np.cos(angle), axis=axis))

    if isDeg:
        std = np.degrees(std)

    return std
#} #End Statistics for angular fields

#{Rotations and projections
def rodrigues(n,theta):
    """RotationMatrix from direction (n) and angle (theta) (Rodrigues)"""
    def S(n):
        Sn = np.array([[0,-n[2],n[1]],[n[2],0,-n[0]],[-n[1],n[0],0]])
        return Sn
    if theta > 1e-30:
        Sn = S(n)
        R = np.eye(3) + np.sin(theta)*Sn + (1-np.cos(theta))*np.dot(Sn,Sn)
    else:
        Sr = S(n)
        theta2 = theta**2
        R = np.eye(3) + (1-theta2/6.)*Sr + (.5-theta2/24.)*np.dot(Sr,Sr)
    return R


def projected2D(n, coords):
    """Calculates the orthogonal projection of coords onto a n-normal plane
    through origin. The returned 2-D-coordinates are expressed with respect to
    the (also returned) inplane unit vectors.

    @param n: normal vector of plane (not necessarily normalized)

    @param coords: 3D-set of coordiantes ([(x,y,z), (xx,yy,zz),...])

    @returns: 2D-coordinate-set [(u,v), (uu, vv),...]
        inplane normal vectors in coords-coordinates nu, nv
    """
    n, coords = np.asarray(n), np.atleast_2d(coords)

    norm = np.linalg.norm(n)

    tol = 1E-10
    if norm < tol:
        raise ValueError('Norm of normal vector is to small.')
    n = n/norm

    # shortcut for axes-parallel projection
    if np.abs(n.dot([0,0,1]) - 1) < tol:
        return np.delete(coords, 2, axis=1), np.array([1,0,0]), np.array([0,1,0])

    elif np.abs(n.dot([0,1,0]) - 1) < tol:
        return np.delete(coords, 1, axis=1), np.array([1,0,0]), np.array([0,0,1])

    elif np.abs(n.dot([1,0,0]) - 1) < tol:
        return np.delete(coords, 0, axis=1), np.array([0,1,0]), np.array([0,0,1])

    else:
        pass

    #1st orthoVector:
    m = np.cross(n, [0,0,1])
    m /= np.linalg.norm(m)
    #2nd orthoVector:
    k = np.cross(n,m)

    #projector matrix
    prj = np.vstack((n,m,k))

    coords = np.array(map(lambda v: np.dot(prj,v), coords))

    return np.delete(coords, 0, axis=1), m, k

#} #End Rotaions and projections

###############################################################################
#{Finding algorithms
def createRTree(nodes, treeBaseName='rtreeMesh', force=False):
    """Creates and stores or reads rTree for elements specified by nodes.

    @param nodes: elementNodes of dim [nCornerNodes=4, nElems, 3DCoords=3],
    not needed if treeBaseName.idx and treeBaseName.dat already exist

    @param treeBaseName: base-filename for rtree-data (creates '.dat' and
    '.idx'). If treeBaseName is None, the rtree will not be loaded/saved from/
    to file. It will be created by a quite fast generator function.

    @param force: if True (re)creation of rtree will be forced even if data
    (treeBaseName.idx/.dat) exists

    @note: the elements-ids in rtree will be the indexes (appearence in .inp)
    and not the element labels! INDEXING STARTS WITH 1! (0 is not allowed by
    Rtree-package)

    @note: this requires the Rtree-package and the lib/dll spatialindex(-dev).
    """
    if not FastRtree:
        msg('ERROR. It seems that rtree is not installed properly.')
        raise ImportError('Can not import rtree.')

    # set rTree properties
    p = index.Property()
    p.dimension = 3

    #some tweaks to speed up
    p.leaf_capacity = 96#1024
    #p.fill_factor = .9

    # concanate min/max of nodes (first idx loops elnodes)
    bounds = np.concatenate((nodes.min(axis=0),
                             nodes.max(axis=0)), axis=1)
    bounds = bounds.tolist()

    if treeBaseName is not None:
        # load tree if exists
        if (os.path.exists(treeBaseName + '.dat') or
            os.path.exists(treeBaseName + '.idx')):
            if not force:
                msg('Using stored rTree for %s' % treeBaseName)
                return index.Index(treeBaseName, properties=p)
            else:
                # clean up so rtree will be (re)generated
                try:
                    os.remove(treeBaseName + '.dat')
                except:
                    pass
                try:
                    os.remove(treeBaseName + '.idx')
                except:
                    pass

        nodes = np.asarray(nodes)
        if not nodes.size:
            raise ValueError('If no existing rTree data path is specified one'
                             ' nodes array is required.')

        # create or load rTree
        idx3d = FastRtree(treeBaseName, properties = p)

        # insert elements to rTree
        # warning!! the elementId will be index and not elementlabel!!!
        for cnt,bound in enumerate(bounds):
            if (cnt % 250000) == 0:
                msg('\t insert done for %d elements' % cnt)
            idx3d.insert(cnt+1, bound)

    else:
        ### filefree generator solution - see rtree docs -- performance
        def generator():
            for cnt, bound in enumerate(bounds):
                yield (cnt+1, bound, None)
        idx3d = index.Index(generator(), properties=p)

    return idx3d


@jit
def rTreeQuery(rTree, points):
    """Querys points in rTree

    @param rTree: rTree object

    @param points: sequence of points to search for

    @note: install numba to speed up
    """

    points = np.atleast_2d(points)

    # points need to "extended" to objects of zero-sized bounding box
    points = np.concatenate((points, points), axis=1)

    # keep in mind that the stored index starts with 1
    indexes = [[tid - 1 for tid in rTree.intersection(tuple(pt))]
                for pt in points]

    return indexes


class KDTree(object):
    """ Some usefull pointSearch functions using scipy's cKDTree implementation
    """
    @staticmethod
    def checkDim(coords):
        """ Forces coords to be a 3-D PointField """
        coords = np.asarray(coords)
        shape = coords.shape
        if len(shape) == 2 and shape[1] == 3:
            return
        else:
            raise ValueError('PointVector needs to be of dimension Npts x 3' +
                             'Yours is %s.' % str(shape))


    @classmethod
    def KDTree(cls, coords):
        """ Builds a scipy cKDTree from coords """
        if not cKDTree:
            msg('ERROR. It seems that scipy is not installed properly.')
            raise ImportError('Can not import cKDTree from scipy.spatial')

        cls.checkDim(coords)
        if coords.shape[0] > 1E6:
            msg('Building KDTree for %d 3D-Points' % coords.shape[0])

        tree = cKDTree(coords)

        return tree


    @staticmethod
    def uniqueTreeIdx(tree, atol=1E-6, normType='sphere'):
        """ Returns unique Points in tree which are unique within a tolerance
        """
        if normType == 'sphere':
            pNorm = 2.
        elif normType == 'box':
            pNorm = np.inf
        else:
            raise ValueError("Don't know a normType %s" % str(normType))

        groups = tree.query_ball_tree(tree, 1E-10, p=pNorm)
        #  group = [[self1, 1st neighbour, 2nd, ...],       #first pt in list
        #           [self2, 1st neighbour, 2nd, ...],...    #second pt
        #           [selflast, 1st neighbour, 2nd, ...]]
        # --> group holds redundant for each neigbour --> use np.unique

        #sort each group
        try:
            groups = np.sort(groups, axis=0)
        except Exception:
            import scipy
            print("This may has failed because of an old Scipy-Version. \n"+
                  "tested for scipy 1.1.0; yours: %s" % scipy.__version__)
            raise


        # groups, rIdx  = np.unique(groups, axis=0, return_index=True)
        # # np.unique returns a sorted list --> use return_index to undo this
        # # sorting
        # groups = groups[rIdx]
        try:
            # 'balanced'-Points, len(group) for group in groups is equal
            # eg. each field has same pointdef
            groups = np.unique(groups, axis=0)
        except TypeError:
            # 'unbalanced'-Points, len(group) for group in groups differ
            # eg. each fields are definend on different point(sub)sets
            groups = np.unique(groups)

        return groups


    @classmethod
    def ensureUnique(cls, arg, atol=1E-6):
        """ checks if a PointCloud (or its KDTree) has only unique Poitns
        within a tolerance
        """
        if not cKDTree:
            msg('ERROR. It seems that scipy is not installed properly.')
            raise ImportError('Can not import cKDTree from scipy.spatial')

        if isinstance(arg, np.ndarray):
            cls.checkDim(arg)
            tree = cls.KDTree(arg)
        elif isinstance(arg, cKDTree):
            tree = arg
        else:
            raise ValueError("Input must be numpy-array or cKDTree (scipy)")

        idx = cls.uniqueTreeIdx(tree, atol=atol)

        if len(idx) == tree.n:
            return True
        else:
            return False


def getPointsAlongPath(start, end, radius, coords):
    """Returns the index of points coords which are enveloped by the cylinder with
    a radius between start and end.

    @param start: start point

    @param end: end point

    @param radius: radius of cylinder

    @param coords: 3D-set of coordiantes ([(x,y,z), (xx,yy,zz),...])

    @returns: indices of coords inside cylinder
    """
    if not cKDTree:
        msg('ERROR. It seems that scipy is not installed properly.')
        raise ImportError('Can not import cKDTree from scipy.spatial')
    # ToDo:
    # extend for a sequence of points (recursive calling + handeling caps)

    ## 1st step: clip coords to maximum bounding box of path+radius
    start, end = np.asarray(start), np.asarray(end)

    # preselect by BoundingBox
    bb = np.vstack((start, end))
    ebb = (bb.min(axis=0) - 2*radius, bb.max(axis=0) + 2*radius)

    # indices in boundingBox
    iib = np.logical_and(np.all(np.greater_equal(coords, ebb[0]), axis=1),
                         np.all(np.less_equal(coords, ebb[1]), axis=1))

    iibN = np.argwhere(iib).flatten()

    ## 2nd step: project (prefiltered) points to n-normal plane
    # normal vector (not normalized yet)
    n = end - start

    # project points to plane

    uv, _, _ = projected2D(n, coords[iib])

    uvStart, _, _ = projected2D(n, start)
    ## 3rd step: get all projected points in radius around startPoint
    tree = cKDTree(uv)
    pts = tree.query_ball_point(uvStart, radius)[0]
    idxR = np.zeros(iibN.shape, dtype=bool)
    idxR[pts] = True

    ## 4th step: find points between start and end
    dotN = lambda v: n.dot(v)
    idxU = np.array(map(dotN, coords[iib] - end)) <= 0
    idxL = np.array(map(dotN, coords[iib] - start)) >= 0

    ## 5th step: get all (sub)points in radius and between start and end
    idxC = np.all(np.vstack((idxL, idxU, idxR)), axis=0)

    ## return indices with respect to complete coords
    return iibN[np.argwhere(idxC).flatten()]

#} #End Finding algorithms

if __name__ == '__main__':
    print('No syntax errors.')