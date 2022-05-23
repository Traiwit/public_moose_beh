"""create rhino objects

"""

import itertools

class RhinoScript(file):
    """Creates a rhino script file. The class is derived from file.

    Basic usage:

    >>> from bae.rhino_01 import RhinoScript
    >>> rh = RhinoScript("myscript.rhino.txt")
    >>> rh.useWorldCoords()
    >>> rh.cmd("Point 0.0,10.0,5.5") # discouraged, ignores useWorldCoords()
    >>> # or (yielding the same result)
    >>> rh.cmd("Point", [0.0,10.0,5.5])
    ...
    """

    fmtfloat = "%.15g"

    def __init__(self, fname, fast=True, clean=False):
        file.__init__(self, fname, 'w')

        self.finalcmds = list()
        if fast:
            self.cmd("SetRedrawOff NoEcho")
            self.finalcmds.append("SetRedrawOn Echo")
        if clean:
            self.cleanAll()

        self.flagWorldCoords = False

    def close(self):
        if self.closed:
            return
        self.cmd("_SelNone")
        self.finalcmds.reverse()
        map(self.cmd, self.finalcmds)
        file.close(self)

    def __del__(self):
        self.close()

    def useWorldCoords(self, onFlag=True):
        """All coordinates specified as parameter to subsequent commands shall
        be treated as world coordinates."""
        self.flagWorldCoords = onFlag

    def cmd(self, *args):
        """Place a rhino command in the script.

        Writes the argument(s) and adds a newline character at the end.
        Strings and integers are written as such, floats are formated according
        to self.fmtfloat. Lists, tuples, iterables of floats are interpreted as
        point coordinates (number of coordinates is not checked) and seperated
        by commas.

        Examples:

        >>> rh = RhinoScript("myscript.txt")
        >>> pts = [[0,0,0], [0,0,1], [1,1,0]]
        >>> rh.cmd("point", pts[0])
        >>> rh.cmd("line", pts[1], pts[2])
        """
        first = True
        for arg in args:

            # seperate parts by spaces
            if first:
                first = False
            else:
                self.write(" ")

            # write content
            if type(arg) is str:
                # string argument is passed to the rhino script unchanged
                self.write(arg)
            elif type(arg) is int:
                # int argument is converted to string
                self.write("%d" % arg)
            elif type(arg) is float:
                # float argument is converted to string
                self.write(self.fmtfloat % arg)
            else:
                try:
                    # list of floats is treated as point coordinates
                    res = ",".join([
                            self.fmtfloat % x for x in arg])
                except TypeError:
                    raise ValueError(
                        "Unrecognized type %s of argument to RhinoScript.cmd()"
                        " or items of this argument not suitable for conversion"
                        " to float." % (str(type(arg))))
                if self.flagWorldCoords:
                    res = "w"+res
                self.write(res)

        # add newline
        self.write("\n")

    def cleanAll(self):
        """Deletes all objects and layers (also hidden and invisible)
        Sets the cplane to the world origin (now world and cplane coordinates
        are the same).
        """
        self.cmd("_-Layer On * Unlock * _Enter")
        self.cmd("_Show")
        self.cmd("_SelAll")
        self.cmd("_Delete")
        self.cmd('_-Layer New "Layer 01" Current "Layer 01" Delete * _Enter')
        self.cmd("_-CPlane _World _Top")


    ##----
    ## Higher level commands


    #----
    # Layer manipulation
    def createNewLayer(self, layerName, color):
        """Create new colored layer (if it does not exist already) and change
        current layer to this.

        >>> from bae.rhino_01 import RhinoScript, layerColors
        >>> rh = RhinoScript("rhino.txt")
        >>> for ob in someobjects:
        >>>    rh.createNewLayer(ob.layerName, layerColors.next())
        >>>    ... create object

        @param color: must be a rgb - tuple: three integers between 0 and 255
        """
        self.cmd("_-Layer new %s current %s color %s"
                 % ((layerName,)*3))
        self.cmd(",".join(map(str,color))+" _Enter")

    def changeLayer(self, layerName):
        """Change current layer
        """
        self.cmd("_-Layer Current %s _Enter" % layerName)



    #----
    # Create geometry
    def addSrfPt(self,ptQuadrupel):
        """Create a nurbs surface of three or four points.
        The points should be ordered in clockwise succession around the border.
        """
        if len(ptQuadrupel)<3 or len(ptQuadrupel)>4:
            raise ValueError(
                "RhinoScript.writeSrfPt() needs three or four points as"
                " argument. %d specified" % len(ptQuadrupel))
        self.cmd("SrfPt")
        for pt in ptQuadrupel:
            self.cmd(pt)
        if len(ptQuadrupel)<4:
            self.cmd(ptQuadrupel[0])

    def writeSrfPt(self,ptQuadrupel):
        """Deprecated alias for RhinoScript.addSrfPt()
        for compatibility reasons.
        """
        self.addSrfPt(ptQuadrupel)

    def addBrick(self,points):
        """Creates a polysurface of six nurbs surfaces
        points are ordered according to abaqus C3D8 brick elements:
        This is a literal block::
             7 ---- 6
            /|     /|
           / |    / |
          4 ---- 5  |
          |  |   |  |
          |  3 --+- 2
          | /    | /
          |/     |/
          0 ---- 1
        
        Not implemented yet.
        """
        raise Exception("Not implemented yet.")


    def addBox(self,box):
        """Create a box aligned to the coordinate axes
        box is a list of two points:
        [[xmin, ymin, zmin], [xmax, ymax, zmax]]
        """
        p1 = box[0]
        p2 = list(box[1])
        p2[2] = p1[2]
        h = box[1][2]-box[0][2]
        self.cmd("Box", p1, p2, h)


    # Rhino - command to create sphere
    # sphere_cmd = '_-Sphere %g,%g,%g %g\n'
    sphere_cmd_first = ('_-MeshSphere VerticalFaces=6 AroundFaces=8'
                        ' %g,%g,%g %g\n')
    sphere_cmd = '_-MeshSphere %g,%g,%g %g\n'
    colorcmd = ('_-SelLast _Enter _-Properties Object Color Object %g,%g,%g'
                ' _Enter _Enter\n')

    def addSpheres(self, points, sizes, colors=None, colormap=None):
        """Creates mesh spheres from a list of points

        @param points: List of point coordinates each item a tuple of three
        floats
        @param sizes: May be a float or a list of the same length as points specifying
        the sizes of the spheres
        @param colors: May be list of values specifying the value to look up in
           colormap to get the rgb-colour for the sphere
        @param colormap: Colormap object from bae.colormap_01.py
        """
        if type(sizes)==float or type(sizes)==int:
            sizes = itertools.repeat(sizes)
        for i, pt, size in itertools.izip(itertools.count(0), points, sizes):

            coord_size = tuple(pt)+(size,)
            if i==0:
                self.cmd(self.sphere_cmd_first % coord_size)
            else:
                self.cmd(self.sphere_cmd % coord_size)

            if colors != None:
                self.cmd(self.colorcmd % colormap.get_rgb_color(colors[i]))

    def addPoints(self, points, colors=None, colormap=None):
        """Creates point objects from a list of point coordinates

        @param points: List of point coordinates each item a tuple of three
        floats

        @param colors: May be list of values specifying the value to look up in
        colormap to get the rgb-colour for the sphere

        @param colormap: Is a Colormap object from bae.colormap_01.py
        """
        if colors != None:
            for pt, col in itertools.izip(points, colors):
                self.cmd("point", pt)
                self.cmd(self.colorcmd % colormap.get_rgb_color(col))
        else:
            for pt in points:
                self.cmd("point", pt)


    ##----
    ## possibly deprecated commands 

    def writeSelAllDelete(self):
        """Deprecated function. Not needed, don't use, will eventuelly
        disapear. Only left for compatibility reasons.
        Use the clean option of the RhinoScript constructor instead."""
        self.cmd("SelAll")
        self.cmd("Delete")

    def writePoint(self,pt):
        """Deprecated function. Not needed, don't use, will eventuelly
        disapear. Only left for compatibility reasons.
        Use self.cmd("Point", pt) instead."""
        self.cmd("Point", pt)

    def writePointWorld(self,pt):
        """deprecated.
        use self.useWorldCoords(); self.cmd("Point", pt) instead"""
        self.cmd("Point w%.15g,%.15g,%.15g" % tuple(pt))

    def writeLineWorld(self,ptPair):
        """deprecated."""
        self.cmd("Line")
        for pt in ptPair:
            self.cmd("w%.15g,%.15g,%.15g"%tuple(pt))

    def writeBoxWorld(self,box):
        """Deprecated function. Not needed, don't use, will eventuelly
        disapear. Only left for compatibility reasons.
        Use self.useWorldCoords(); self.addBox(...) instead."""
        x1,y1,z=box[0][0],box[0][1],box[0][2]
        x2,y2=box[1][0],box[1][1]
        h=box[1][2]-box[0][2]
        self.cmd("Box")
        self.cmd("w%.15g,%.15g,%.15g"%(x1,y1,z))
        self.cmd("w%.15g,%.15g,%.15g"%(x2,y2,z))
        self.cmd("%.15g"%(h))

    def writeSrfPtWorld(self,ptQuadrupel):
        """deprecated.
        use self.useWorldCoords(); self.addSrfPt(...) instead"""
        self.cmd("SrfPt")
        for pt in ptQuadrupel:
            self.cmd("w%.15g,%.15g,%.15g"%tuple(pt))
        if len(ptQuadrupel)<4:
            self.cmd("_Enter")

    def writeDelete(self):
        """Deprecated function. Not needed, don't use, will eventuelly
        disapear. Only left for compatibility reasons.
        Use self.cmd("Delete") instead."""
        self.cmd("Delete")

    def writeSelLast(self):
        """Deprecated function. Not needed, don't use, will eventuelly
        disapear. Only left for compatibility reasons.
        Use self.cmd("SelLast") instead."""
        self.cmd("SelLast")

    def writeExtrudeSrf(self,height):
        """Deprecated function. Don't use, will eventuelly
        disapear. Only left for compatibility reasons."""
        self.cmd("ExtrudeSrf",
                 "d",
                 "d",        # direction
                 "0.,0.,0.", # base point for direction
                 "0.,0.,1.", # second point for direction
                 height)


class ColorsIter(itertools.cycle):
    """Loop infinitely over a set of rgb colors.

    usage of the predefined layerColors:

    >>> from bae.rhino_01 import RhinoScript, layerColors
    >>> rh = RhinoScript("rhino.txt")
    >>> for ob in someobjects:
    >>>    rh.createNewLayer(ob.layerName, layerColors.next())
    >>>    ... create object
    ... or:  (even better: use izip from module itertools)
    >>> for ob, col in zip(someobjects, layerColors):
    >>>    rh.createNewLayer(ob.layerName, col)
    >>>    ... create object

    define another colour sequence

    >>> colorsRGB = ColorsIter([(255,0,0),(0,255,0),(0,0,255)])
    """
    pass

layerColors = ColorsIter([
        (127,  0,0),   # brown
        (255,  0,0),   # red
        (255,127,0),   # orange
        (  0,255,0),   # green
        (  0,127,0),   # dark green
        (  0,255,255), # cyan
        (  0,  0,255), # blue
        (127,  0,255), # purple
        (255,  0,255), # magenta
        (255,191,191), # pink
        (255,255,255), # white
        ])
