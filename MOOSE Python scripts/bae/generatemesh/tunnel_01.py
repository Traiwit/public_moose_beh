"""tunnel_01.py

This module defines a class that describes properties of a tunnel excavation
(not the possible GS installation in it).

Usage
=====
  >>> from bae.generatemesh.tunnel_01 import XSectCircBack
  >>> m = XSectCircBack()
  >>> help(m) # for more information
"""

__version__ = "1.1"

_version_history_ = """\
Versions:
=========

tunnel_01
1.0 : GP: derived from Perilya GS project scripts
1.1 : GP added: normal to the tunnel contours
"""

from math import asin, sqrt, sin, cos, pi
# from bae.future_01 import *
# from bae.misc_01 import MsgTicker, Container

class Tunnel(object):
    """geometry of a tunnel excavation
    """

class StraightTunnel(Tunnel):
    """A straight piece of tunnel

    The global coordinate system defaults to
     - x_g pointing right and being x of the cross sectional coordinate system
     - y_g along the tunnel
     - z_g pointing up, being y of the cross sectional coordinate system.

    @ivar crossSection: object of type TunnelCrossSection
    @ivar length: ... of (this part of) the tunnel
    @ivar orientation: Unit vector specifying the orientation of the tunnel
    main direction in the global coordinate system. Defaults to [0,1,0],
    the tunnel stretching in north-south direction.
    @ivar origin: Start coordinate of the tunnel (centre of the bottom ground)
    """
    def __init__(self, crossSection, length,
                 heading=0.0, origin=[0.0, 0.0, 0.0]):
        """
        @param crossSection: of type TunnelCrossSection
        @param length: length of this piece of tunnel
        @param heading: orientation of the tunnel axis in degrees clockwise
        from the global y-axis (north).
        @param origin: start coordinate of the tunnel (centre of the bottom
        ground)
        """
        self.crossSection = crossSection
        self.length = length
        self.heading = heading
        self.orientation = [sin(heading*pi/180), cos(heading*pi/180), 0.0]
        self.origin = origin

    def ys2xyzContour(self, y, s):
        """Get the global x,y,z coordinates of a point on the wall or back of
        the tunnel identified by its cross sectional s coordinate and the
        position y along the tunnel.

        s=0 the is top centre point of the tunnel cross section. Positive s 
        direction goes clockwise along the tunnel cross section when looking
        along the tunnel in positive y direction.

        y is the distance from self.origin, positive in the direction specified
        by self.orientation.
        """
        x_xs, y_xs = self.crossSection.xyContour(s)
        ds, dc, dummy = self.orientation
        p0 = self.origin
        return [ x_xs*dc + y*ds + p0[0],
                -x_xs*ds + y*dc + p0[1],
                y_xs + p0[2]]

    def ys2normal(self, y, s):
        """Get the normal vector perpendicular to the tunnel contour at a
        of a point on the wall or back of the tunnel identified by its cross
        sectional s coordinate and the position y along the tunnel. The normal
        points outwards.

        s=0 the is top centre point of the tunnel cross section. Positive s 
        direction goes clockwise along the tunnel cross section when looking
        along the tunnel in positive y direction.

        y is the distance from self.origin, positive in the direction specified
        by self.orientation.

        @returns: global x,y,z components of the normal vector.
        """
        x_xs, y_xs = self.crossSection.xyNormal(s)
        ds, dc, dummy = self.orientation
        return [ x_xs*dc,
                -x_xs*ds,
                 y_xs]

    def getOffsetContour(self, offset):
        """Return a new Tunnel-object with a contour that is offset to
        the current one by the given parameter. Positive values means the new
        contour is bigger.
        """
        otherXSect = self.crossSection.getOffsetContour(offset)
        otherTunnel = StraightTunnel(
            otherXSect, length=self.length, heading=self.heading)
        return otherTunnel

    def sFromBottom(self, sInv):
        """Return the s coordinate corresponding to an inverted s-coordinate
        sInv counting from the bottom of the drive (y_xs=0) upwards along the
        contour.

        @Note: Only sensible for sInv>=0. Refer to the documentation of the
        sFromBottom() method of the corresponding cross section
        self.crossSection.sFromBottom()
        """
        return self.crossSection.sFromBottom(sInv)

##########################################################################

class TunnelCrossSection(object):
    """
    The tunnel cross sectional coordinate system x_xs to the right, y_xs up,
    with the origin at the centre of the bottom::

         ^ y_xs
         |
         |

       ,--.
      /    \
      |    |
      |    |       x_xs
      +----+  --->

    Additionally there is the curved coordinate s defined along the cross
    section contour. The origin is in the middle of the back and it counts
    positive to the right (initially parallel to x).
    """


class XSectCircBack(TunnelCrossSection):
    """Tunnel contour with vertical walls and a circular section at the back
    (top).

    The coordinates used are described in the parent class TunnelCrossSection.

    The current implementation ignores the bottom, which is assumed to always
    be at y=0.

    @ivar width: width of the tunnel cross section
    @ivar height: heigth of the tunnel cross section
    @ivar backCircCy: centre of the circular back, y coord
    @ivar sMaxBack: (positive) s coordinate where the circular section of
    the back ends and the vertical walls start.
    """

    def __init__(self, width, **kwargs):
        """Constructor

        You will have to specify the width and either of the following:
         - height and backRadius
         - height and wallHeight
         - wallHeight and backRadius

        Only the first option is implemented already.

        @param width: width of the tunnel cross section
        @keyword height: (optional) heigth of the tunnel cross section
        @keyword backRadius: (optional) specify the radius of a circular top
            (back) section.
        """

        self.width = float(width)

        if (("height" in kwargs)
            and ("backRadius" in kwargs)):
            self.backRadius = float(kwargs["backRadius"])
            self.height = float(kwargs["height"])
            self.backCircCy = self.height - self.backRadius
        else:
            raise ValueError(
                "XSectCircBack did not get sufficient arguments.")

        self.sMaxBack = self.backRadius*asin(0.5*self.width/self.backRadius)
        self.yMaxWall = self.backCircCy + sqrt(self.backRadius**2
                                               - (0.5*self.width)**2)

    def xyContour(self, s):
        """Return contour of the tunnel. s=0 is the top centre point. Positive
        s values lead to positive x values.

        @Returns: a list [x,y]
        """
        if abs(s)<self.sMaxBack:
            # circular back section
            alpha = float(s)/self.backRadius
            x = self.backRadius*sin(alpha)
            y = self.backCircCy + self.backRadius*cos(alpha)
        else:
            # vertical section
            if s<0:
                x = -0.5*self.width
            else:
                x = 0.5*self.width
            y = self.yMaxWall-(abs(s)-self.sMaxBack)
        return [x,y]

    def xyNormal(self, s):
        """Return normal vector in the tunnel x_xs-y_xs plane. s=0 is the top
        centre point. Positive s values lead to positive x values.

        @Returns: a list [x,y] stating the normal vector components
        """
        if abs(s)<self.sMaxBack:
            # circular back section
            alpha = float(s)/self.backRadius
            return [sin(alpha), cos(alpha)]
        else:
            # vertical section
            if s<0:
                x = -1.0
            else:
                x = 1.0
            return [x, 0.0]

    def getOffsetContour(self, offset):
        """Return a new XSectCircBack-object with a contour that is offset to
        the current one by the given parameter. Positive values means the new
        contour is bigger.
        """
        other = XSectCircBack(
            width = self.width + 2.0*offset,
            height = self.height + offset,
            backRadius = self.backRadius + offset
            )
        return other

    def sFromBottom(self, sInv):
        """Return the s coordinate corresponding to an inverted s-coordinate
        sInv counting from the (s>0) bottom corner of the drive (y_xs=0)
        upwards along the contour (i.e. along the side wall, then further
        across the back).

        @Note: Only sensible for sInv>=0. Not sensible for positions on the
        floor of the tunnel.
        """
        return self.yMaxWall + self.sMaxBack - sInv


class XSectCircCorner(TunnelCrossSection):
    """Tunnel contour with vertical walls. The back is horizontal in the middle
    with quarter circular corners.

    The coordinates used are described in the parent class TunnelCrossSection.

    The current implementation ignores the bottom, which is assumed to always
    be at y=0.

    @ivar width: width of the tunnel cross section
    @ivar height: heigth of the tunnel cross section
    @ivar sMaxHoriz: (positive) s coordinate where the horizontal part of the
    back ends
    @ivar sMaxBack: (positive) s coordinate where the circular sections of
    the back ends and the vertical walls start.
    @ivar yMaxWall: height of the vertical wall sections and y coordinate of
    the centre of the corner circles
    """

    def __init__(self, width, height, cornerRadius):
        """Constructor
        @param width: width of the tunnel cross section
        @param height: heigth of the tunnel cross section
        @param cornerRadius: radius of the quarter circular corners
        """

        if 2.0*cornerRadius>width:
            raise ValueError(
                "XSectCircCorner(): width<2*cornerRadius")

        self.width = float(width)
        self.height = float(height)
        self.cornerRadius = float(cornerRadius)

        self.sMaxHoriz = 0.5*self.width - self.cornerRadius
        self.sMaxBack = 0.5*self.cornerRadius*pi + self.sMaxHoriz
        self.yMaxWall = self.height - self.cornerRadius

    def xyContour(self, s):
        """Return contour or the tunnel. s=0 the is top centre point. Positive
        s values lead to positive x values.

        @Returns: a list [x,y]
        """
        if abs(s)<=self.sMaxHoriz:
            # horizontal back section
            x = s
            y = self.height
        elif abs(s)<self.sMaxBack:
            # circular corner section
            if s>=0:
                sHor = self.sMaxHoriz
            else:
                sHor = -self.sMaxHoriz
            alpha = float(s-sHor)/self.cornerRadius
            x = self.cornerRadius*sin(alpha) + sHor
            y = self.yMaxWall + self.cornerRadius*cos(alpha)
        else:
            # vertical section
            if s<0:
                x = -0.5*self.width
            else:
                x = 0.5*self.width
            y = self.yMaxWall+self.sMaxBack-abs(s)
        return [x,y]

    def xyNormal(self, s):
        """Return normal vector in the tunnel x_xs-y_xs plane. s=0 is the top
        centre point. Positive s values lead to positive x values.

        @Returns: a list [x,y] stating the normal vector components
        """
        if abs(s)<=self.sMaxHoriz:
            # horizontal back section
            return [0.0, 1.0]
        elif abs(s)<self.sMaxBack:
            # circular corner section
            if s>=0:
                sHor = self.sMaxHoriz
            else:
                sHor = -self.sMaxHoriz
            alpha = float(s-sHor)/self.cornerRadius
            return [sin(alpha), cos(alpha)]
        else:
            # vertical section
            if s<0:
                x = -1.0
            else:
                x = 1.0
            return [x, 0.0]

    def getOffsetContour(self, offset):
        """Return a new XSectCircCorner-object with a contour that is
        offset to the current one by the given parameter. A positive value
        means the new contour is bigger.
        """
        other = XSectCircCorner(
            width = self.width + 2.0*offset,
            height = self.height + offset,
            cornerRadius = self.cornerRadius + offset
            )
        return other

    def sFromBottom(self, sInv):
        """Return the s coordinate corresponding to an inverted s-coordinate
        sInv counting from the (s>0) bottom corner of the drive (y_xs=0)
        upwards along the contour (i.e. along the side wall, then further
        across the back).

        @Note: Only sensible for sInv>=0. Not sensible for positions on the
        floor of the tunnel.
        """
        return self.yMaxWall + self.sMaxBack - sInv
