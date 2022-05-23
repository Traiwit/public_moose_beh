"""volume_01 - module: DEPRECATED, use L{bae.volume_02}

some functionality concerning a distinct volume in space

# initialization
from bae.volume_01 import Volume
v = Volume()
elset = [1234, 4567, ...]  # elset contains elemnt numbers
v.create_from_tet_elem(node_coords, el_nodes, elset)

# move this volume up by 10
v.translate(move_vector=[0,0,10])

# scale this volume by 2 in each direction (volume := volume*8)
v.scale(scale_factor=2, scale_origin=[0,0,0])

# test if nodes are in that volume
v.initialize_point_search()
for point in [[0,0,0], [1.0, 0.4, 1.6], [-0.1, 2.0, 0.3]]:
   print point, "is inside:", v.point_is_inside(point)

# bounding box of that volume
bb = v.getboundingbox()
print "lower left front corner:", bb[0]
print "upper right back corner:", bb[1]

# 'findAtMesh':
# return an elset from model with all its elements inside the volume
# model is an abq_model instance,
# elset may be an elset name or a list of element numbers specifying the subset
# of model to search elements,
# optional method parameter
v.intersection_elset(model, elset)

"""


"""
Versions:

1.01 new   - pt_is_inside() and initialize_point_search() derived
             from findAtMesh 1.17
             scale, copy, getboundingbox methods
1.02 added - translate method, more doc, added smooth mode for point finding
1.03 added - intersection_elset,
             don't need to specify elset for create_from_tet_elem
             fixed: all from future
             added: VolumeBase class, VolumeExtrudeOutline,
               Volume.get_exterior_surface() method
             changed: introduced tolerance in point_is_inside methods
             

known bugs/problems:
- intersection_elset lacks the gauss point method

todo:
* volume_02

  Volume = base class not to be used
  gets the intersection_elset method

  child class VolumeFromTetMesh takes abq_model.Model object with optional elset
  as arguments to the constructor

  child classes sphere, box (xyz-aligned)

* rotate method

* point_is_inside ... faster!

"""

from bae.vecmath_01 import *
from bae.future_01 import *
from bae.abq_model_01 import QuietLog, Model
from bae.surface_01 import TriangleSurface
from sys import stdout
import copy


class VolumeBase(object):
    def intersection_elset(self, model, elset = None, method = 'centroid'):
        """find all elements from the elset in model that are within the volume

        model is an abq_model_01.Model - object with el_node, node_coord and
              elset member attributes (elset member only needed if argument
              elset is a string)
              and an calc_centroid(element_id)-method
        elset is a string (name of the elset in model) or a list of element
              numbers in model or None. If None all elements of model are
              tested, otherwise only those in the set.
        method may be 'centroid'/'c' or 'gauss point'/'gp'

        bugs/limitations:
        - elements are assumed to be tet elements with only the corner nodes
          (first four) considered
        """

        method = method.replace(' ', '')

        # hot fix for abq_model_02
        if model.__class__.__module__.split('.')[-1]=='abq_model_02':
            model.el_node = model.elNodes
            model.node_coord = model.nodeCoords

        if elset==None:
            # caution, this is fast but dangerous if you modify the code
            elset = model.el_node
        elif type(elset) == str:
            if not hasattr(model, 'elset'):
                raise Exception('Volume.intersection_elset: model has no elset'
                                ' member attribute. Cannot intersect elset %s.'
                                % elset)
            elset = model.elset[elset]
        else:
            try:
                elset = iter(elset)
            except TypeError:
                raise TypeError('Volume.intersection_elset: elset is neither'
                                ' None nor a string. It cannot be converted to'
                                ' an iterable either. Cannot intersect this'
                                ' elset.')

        if not hasattr(model, 'el_node'):
            raise AttributeError('Volume.intersection_elset: model has'
                                 ' no el_node member attribute.')
        if not hasattr(model, 'node_coord'):
            raise AttributeError('Volume.intersection_elset: model has'
                                 ' no node_coord member attribute.')

        # create list with tuples (el numer, [points to check])
        el_pts = list()
        if method in ('centroid', 'c'):
            for el in elset:
                try:
                    nodes = model.el_node[el]
                except KeyError:
                    raise KeyError('Volume.intersection_elset: element %d'
                                   ' is not in the model.' % el)
                try:
                    x=y=z=0
                    for node in nodes[:4]:
                        coords=model.node_coord[node]
                        x=x+coords[0]/4.0
                        y=y+coords[1]/4.0
                        z=z+coords[2]/4.0
                except KeyError:
                    raise KeyError('Volume.intersection_elset: node %d'
                                   ' is not in the model.' % node)
                el_pts.append((el, [[x,y,z], ]))
        else:
            raise Exception('Volume.intersection_elset: method <%s> not'
                            ' implemented so far.' % method)

        # check if all checkpoints in the volume
        els_inside = list()
        for el, points in el_pts:
            if all([self.point_is_inside(pt) for pt in points]):
                els_inside.append(el)

        return els_inside

###################################################################
###################################################################
###################################################################
###################################################################

class Volume(VolumeBase):
    """volume class

    You may supply an open logfile for diagnostic output or None to suppress
    it to the constructor.
    """

    cellSizeFactor=0.5      # cell Size for array relative to 'avgSize'
    ident_El=None


    def __init__(self, logfile = stdout):
        "if logfile == None, no diagnostic output"
        if logfile == None:
            self.logfile = QuietLog()
        else:
            self.logfile = logfile
        
        # internal variable, use getboundigbox method! (initialize it)
        self.boundingBox = None

        # variables used for finding points in the volume already set?
        self.point_search_initialized = False
        return

    def copy(self):
        """make a (shallow) copy of this volume

        only copies references to coords, el_nodes, elset
        """
        new_vol = Volume()
        new_vol.coords = self.coords
        new_vol.el_nodes = self.el_nodes
        new_vol.elset = self.elset
        new_vol.logfile = self.logfile

        new_vol.boundingBox = None
        new_vol.point_search_initialized = False

        return new_vol

    def create_from_tet_elem(self, node_coords, el_nodes, elset=None):
        """creates a volume from some tet elements
        node_coords is a dict {node_id: node-coord-tuple}
        el_nodes is a dict {element number: node_id-tuple}
        elset is a list of element numbers

        no deep copy of the arguments is performed, if the objects change
        outside this object, it will affect this object.

        node_coords and el_nodes may contain more than the nodes and elements
        needed for the actual set elset.
        
        Elements in elset that are not found in el_nodes are silently ignored.
        """
        self.coords = node_coords
        self.el_nodes = el_nodes
        if elset==None:
            self.elset = el_nodes.keys()
        else:
            self.elset = set(elset).intersection(el_nodes)

        if len(self.elset)==0:
            self.logfile.write('WARNING: Volume initialized with empty element'
                               ' set.\n')
        return


    ## ------------------------------------------------------------------------
    ## some info

    def getboundingbox(self):
        """Bounding Box of the whole volume.
        Returns a list of two coordinate tuples (rather lists)
        The first coord tuple states the min value, the last the max.
        I.e. self.getboundingbox()[i][j]:
        i \in {0,1} for min/max value, j \in {0,1,2} coordinate index
        """
        if self.boundingBox:
            return self.boundingBox

        # nodes set (contains unique nodes)
        nodes = set()
        for element in self.elset:
            nodes.update(self.el_nodes[element])

        # initialize bounding box with first node
        self.boundingBox=[[0,0,0],[0,0,0]]

        # update bounding box
        first = True
        for node in nodes:
            if first:
                for i in range(3): # x,y,z
                    self.boundingBox[0][i]=self.coords[node][i]
                    self.boundingBox[1][i]=self.coords[node][i]
                first = False
                continue
            coords=self.coords[node]
            for i in range(3): # x,y,z
                if coords[i]<self.boundingBox[0][i]: # min
                    self.boundingBox[0][i]=coords[i]
                if coords[i]>self.boundingBox[1][i]: # max
                    self.boundingBox[1][i]=coords[i]

        return self.boundingBox




    ## ------------------------------------------------------------------------
    ## get the exterior surface of this volume as TriangleSurface object

    def get_exterior_surface(self):
        """return a TriangleSurface object with all outer triangles of this
        volume.
        
        outer triangles are those with exactly one tet connected to.
        
        IMPORTANT NOTE: Node coordinates are passed by reference, if the volume
        is moved or scaled, the surface changes accordingly."""
        
        faceToTet = defaultdict(list)
        for thisTet in self.elset:
            nodesOfThisTet = self.el_nodes[thisTet]
            for face in self.tetFacesIter(nodesOfThisTet):
                faceToTet[face].append(thisTet)

        model = Model()
        last_elnum = 0
        triList = list()
        for face, tets in faceToTet.iteritems():
            if len(tets)!=1: continue
            last_elnum += 1
            model.update_elem(last_elnum, 'S2', list(face))

        surf = TriangleSurface('EXT_SURF', model.el_node.keys(),
                               model.el_node, self.coords, self.logfile)
        return surf



    ## ------------------------------------------------------------------------
    ## methods for testing if a point lies inside the volume or outside

    def initialize_point_search(self, smooth_findradius=0.0):
        """call this function once before point_is_inside
        (it is done automatically)

        smoothmode: Additionally consider points to be inside the volume that
        are in fact outside the volume but within the smooth_findradius to the
        closest centre point of any of the elements of the volume.
        """

        # list of centroid-coord-tuples
        self.tetCentroids=dict()

        if len(self.elset)==0:
            return

        # list of [[xmin,ymin,zmin],[xmax,ymax,zmax]] for each element
        elementBox=dict()

        self.boundingBox = None
        avgSize=0
        for element in self.elset:
            nodes = self.el_nodes[element]
            xyz=[[],[],[]]
            x=y=z=0
            for node in nodes:
                coords=self.coords[node]
                for i in range(3): # for min max
                    xyz[i].append(coords[i])
                x=x+coords[0]/4.0 # for centroid
                y=y+coords[1]/4.0
                z=z+coords[2]/4.0
            self.tetCentroids[element]=[x,y,z]

            elementBox[element]=[[min(xyz[0]),min(xyz[1]),min(xyz[2])],
                                 [max(xyz[0]),max(xyz[1]),max(xyz[2])]]
            
            # update bounding box
            if self.boundingBox:
                for i in range(3): # x,y,z
                    if elementBox[element][0][i]<self.boundingBox[0][i]: # min
                        self.boundingBox[0][i]=elementBox[element][0][i]
                    if elementBox[element][1][i]>self.boundingBox[1][i]: # max
                        self.boundingBox[1][i]=elementBox[element][1][i]
            else:
                self.boundingBox=[[0,0,0],[0,0,0]]
                for i in range(3): # x,y,z
                    self.boundingBox[0][i]=elementBox[element][0][i]
                    self.boundingBox[1][i]=elementBox[element][1][i]

            avgSize+=dist(elementBox[element][1],elementBox[element][0])

        avgSize /= len(self.elset)

        # in method point_is_inside() tolerance is compared to
        # dot(cross(...)...) of vectors of avgSize magnitude
        self.tolerance = (1E-4 * avgSize)**3

        self.logfile.write(
            'bounding box of the volume: [(%.1f,%.1f,%.1f),(%.1f,%.1f,%.1f)]\n'
            % (self.boundingBox[0][0],self.boundingBox[0][1],self.boundingBox[0][2],
               self.boundingBox[1][0],self.boundingBox[1][1],self.boundingBox[1][2]))
        self.logfile.write('avg. element size: %g\n' % avgSize)


        ## -----------------------------------------------------------------
        ## map elements to cells

        self.cellOrigin=self.boundingBox[0]
        cellDimensions=[0,0,0]
        self.cellSize=avgSize*self.cellSizeFactor
        self.cellNumber=[0,0,0]

        # x,y,z number of cells
        self.logfile.write('cell array dimensions:')
        for i in range(3):
            cellDimensions[i]=self.boundingBox[1][i]-self.boundingBox[0][i]
            self.logfile.write(" "+str(cellDimensions[i]))

        self.logfile.write('\ncell numbers:')
        for i in range(3):
            self.cellNumber[i]=int(cellDimensions[i]/self.cellSize)+1
            self.logfile.write(" "+str(self.cellNumber[i]))
        
        self.logfile.write(' cells: %d\n' % (
                self.cellNumber[0]*self.cellNumber[1]*self.cellNumber[2]))

        # create empty array
        for i in range(self.cellNumber[0]):
            if i==0: self.cells=[]
            for j in range(self.cellNumber[1]):
                if j==0: LJ=[]
                for k in range(self.cellNumber[2]):
                    if k==0: LK=[]
                    LK.append([])
                LJ.append(LK)
            self.cells.append(LJ)

        # map elements 
        for element in self.elset:
            # get index ranges
            cellIndex=[[0,0,0],[0,0,0]] # low,high
            for l in range(2):
                for i in range(3):
                    cellIndex[l][i]=int((elementBox[element][l][i]
                                         -self.cellOrigin[i])/self.cellSize)
            # append element ID to cells that overlap the element's bounding box
            for i in range(cellIndex[0][0],cellIndex[1][0]+1):
                for j in range(cellIndex[0][1],cellIndex[1][1]+1):
                    for k in range(cellIndex[0][2],cellIndex[1][2]+1):
                        self.cells[i][j][k].append(element)

        # for smooth mode store radius
        self.smoothmode = smooth_findradius>0
        self.smooth_findradius = float(smooth_findradius)

        # finished initializing
        self.point_search_initialized = True
        return

    def point_is_inside(self, point):
        """test if point is inside the volume
        point is tuple of coordinates
        """

        if len(self.elset)==0:
            return False

        if not self.point_search_initialized:
            self.initialize_point_search()

        # cell check (includes bounding box test from cell index)
        index=[0,0,0]
        for i in range(3):
            index[i]=int((point[i]-self.cellOrigin[i])/self.cellSize)
            if index[i]<0 or index[i]>=self.cellNumber[i]:
                return False

        # if cell check passed ... go on
        elems_in_cell=self.cells[index[0]][index[1]][index[2]]
        centroid_to_point_distance=[]
        # sort by distance
        for element in elems_in_cell:
            l=dist(point,self.tetCentroids[element])
            centroid_to_point_distance.append([l,element])
        centroid_to_point_distance.sort()

        # full check if point inside any tet of elems_in_cell
        inside=0

        if self.smoothmode:
            if (centroid_to_point_distance and
                (centroid_to_point_distance[0][0]<self.smooth_findradius)):
                inside = 1

        for [l,element] in centroid_to_point_distance:
            nodes=self.el_nodes[element]
            # vectors
            v12=vector(self.coords[nodes[0]],
                       self.coords[nodes[1]])
            v13=vector(self.coords[nodes[0]],
                       self.coords[nodes[2]])
            v14=vector(self.coords[nodes[0]],
                       self.coords[nodes[3]])
            v23=vector(self.coords[nodes[1]],
                       self.coords[nodes[2]])
            v24=vector(self.coords[nodes[1]],
                       self.coords[nodes[3]])
            # normals (element face inward)
            n1=cross(v12,v13)
            n2=cross(v14,v12)
            n3=cross(v24,v23)
            n4=cross(v13,v14)
            # test vectors (corner-centroid)
            c1=vector(self.coords[nodes[0]],point)
            c2=vector(self.coords[nodes[1]],point)
            # test inside
            if ((dot(n1,c1)>-self.tolerance) and
                (dot(n2,c1)>-self.tolerance) and
                (dot(n3,c2)>-self.tolerance) and
                (dot(n4,c1)>-self.tolerance)):
                inside=1
                Volume.ident_El=element
                break
        return inside

    ## ------------------------------------------------------------------------
    ## ------------------------------------------------------------------------
    ## transforming the volume

    def scale(self, scale_factor, scale_origin):
        """Scaling is done by simply moving the point coordinates for this
        volume. Other volumes or the original mesh from which the volume
        originates are not affected.

        This creates a new node-dictionary, which is refered to from now on.
        It contains only those nodes connected to any element in self.elset.

        After scaling point search has to be initialized again,
        self.point_search_initialized is set to False.
        """

        all_nodes = set()
        for element in self.elset:
            all_nodes.update(self.el_nodes[element])
        
        new_coords = dict()
        for node_id in all_nodes:
            orig2old = vector(scale_origin, self.coords[node_id])
            orig2new = vector_scale(orig2old, scale_factor)
            new_coords[node_id] = vector_plus(scale_origin, orig2new)
        
        self.coords = new_coords

        # after scaling point search has to be initialized again
        self.point_search_initialized = False


    def translate(self, move_vector):
        """Translating is done by simply moving the point coordinates for this
        volume. Other volumes or the original mesh from which the volume
        originates are not affected.

        This creates a new node-dictionary, which is refered to from now on.
        It contains only those nodes connected to any element in self.elset.

        The point search initialization is modified and doesn't have to be done
        again.
        """

        all_nodes = set()
        for element in self.elset:
            all_nodes.update(self.el_nodes[element])
        
        new_coords = dict()
        for node_id in all_nodes:
            new_coords[node_id] = vector_plus(self.coords[node_id],move_vector)
        
        self.coords = new_coords

        # adapt point search properties
        if self.point_search_initialized and len(self.elset)>0:
            
            # tet Centroids
            for centroid in self.tetCentroids.itervalues():
                vector_modif_add(centroid, move_vector)

            # BoundingBox
            vector_modif_add(self.boundingBox[0], move_vector)
            vector_modif_add(self.boundingBox[1], move_vector)

            # search cells: nothing to modify here
            # self.cellOrigin is just an alias for self.boundingBox[0]


    ## ------------------------------------------------------------------------
    ## some service functions: iterator function

    @staticmethod
    def tetFacesIter(nodesOfThisTet):
        """Returns an iterator of the faces of the given tet element
        nodesOfThisTet gives the node ids, it is converted to a tuple

        use as in
        myTet = mdb.el_node[5731]  # list of node ids, e.g. [3, 43, 12, 36]
        for faces in tetFacesIter(myTet):
           ... faces will be frozenset((3, 43, 12)), next time
               frozenset((43, 3, 36)) and so on.

        An extra comment on the index ordering used in the for loop below: If
        the points order would stay that of the used tuple, each faces normal
        points inwardly.
        """
        nodesOfThisTet = tuple(nodesOfThisTet)
        for i1,i2,i3 in ((0,1,2), (1,0,3), (2,1,3), (0,2,3)):
            yield frozenset((nodesOfThisTet[i1], nodesOfThisTet[i2],
                             nodesOfThisTet[i3]))



###################################################################
###################################################################
###################################################################
###################################################################


class VolumeExtrudeOutline(VolumeBase):
    """volume class

    The given surface is projected on the x-y-plane and extruded in z direction
    from zminmx[0] to zminmax[1].

    This volume object keeps a reference of the triangle nodes and coordinates
    of the surface. So don't modify the surface.

    You may supply an open logfile for diagnostic output or None to suppress
    it to the constructor.
    """

    cellSizeFactor=0.5      # cell Size for array relative to 'avgSize'

    def __init__(self, surface, zminmax, logfile = stdout):

        self.triNodes = surface.triNodes
        self.coords = surface.nodeCoords
        self.zminmax = list(zminmax)

        # if logfile == None, no diagnostic output
        if logfile == None:
            self.logfile = QuietLog()
        else:
            self.logfile = logfile
        
        # internal variable, use getboundigbox method! (initialize it)
        self.boundingBox = None

        # variables used for finding points in the volume already set?
        self.point_search_initialized = False
        return


    def set_zminmax(self, *zminmax):
        """change z-range of the volume
        accepts two values: vol.set_zminmax(zmin, zmax)
        or a list or tuple: vol.set_zminmax([zmin, zmax])
        """
        if len(zminmax)==1:
            self.zminmax = map(float, zminmax[0])
        elif len(zminmax)==2:
            self.zminmax = map(float, zminmax)
        else:
            raise ValueError("set_zminmax() expects two values or one tuple!")


    def initialize_point_search(self, smooth_findradius=0.0):
        """call this function once before point_is_inside
        (it is done automatically)

        smoothmode: Additionally consider points to be inside the volume that
        are in fact outside the volume but within the smooth_findradius to the
        closest centre point of any of the elements of the volume.
        """

        # list of centroid-coord-tuples
        self.triCentroids=dict()

        # list of [[xmin,ymin],[xmax,ymax]] for each tri element
        elementBox=dict()

        self.boundingBox = None
        avgSize=0

        # size of the projection on the x-z plane
        # pos facing one direction, neg facing the other
        elSignedSize = dict()
        for element,nodes in self.triNodes.iteritems():

            node_pts = [self.coords[node]
                        for node in nodes]
            v12=vector(node_pts[0], node_pts[1])
            v13=vector(node_pts[0], node_pts[2])
            detA = v12[0]*v13[1] - v12[1]*v13[0]
            elSignedSize[element] = 0.5*detA

            xy=[[],[]]
            x=y=0
            for coords in node_pts:
                for i in range(2): # for min max
                    xy[i].append(coords[i])
                x=x+coords[0]/3.0 # for centroid
                y=y+coords[1]/3.0
            self.triCentroids[element]=[x,y]

            elementBox[element]=[[min(xy[0]),min(xy[1])],
                                 [max(xy[0]),max(xy[1])]]
            
            # update bounding box
            if self.boundingBox:
                for i in range(2): # x,y
                    if elementBox[element][0][i]<self.boundingBox[0][i]: # min
                        self.boundingBox[0][i]=elementBox[element][0][i]
                    if elementBox[element][1][i]>self.boundingBox[1][i]: # max
                        self.boundingBox[1][i]=elementBox[element][1][i]
            else:
                self.boundingBox=[[0,0],[0,0]]
                for i in range(2): # x,y
                    self.boundingBox[0][i]=elementBox[element][0][i]
                    self.boundingBox[1][i]=elementBox[element][1][i]

            avgSize+=self.dist_2D(
                elementBox[element][1],elementBox[element][0])

        avgSize /= len(self.triNodes)

        # in method point_is_inside() tolerance is compared to values 0..1
        self.tolerance = 1E-4

        self.logfile.write(
            'bounding box of the surface: [(%.1f,%.1f),(%.1f,%.1f)]\n'
            % (self.boundingBox[0][0],self.boundingBox[0][1],
               self.boundingBox[1][0],self.boundingBox[1][1]))
        self.logfile.write('avg. element size: %g\n' % avgSize)


        ## -----------------------------------------------------------------
        ## map elements to cells

        self.cellOrigin=self.boundingBox[0]
        cellDimensions=[0,0]
        self.cellSize=avgSize*self.cellSizeFactor
        self.cellNumber=[0,0]

        # x,y,z number of cells
        self.logfile.write('cell array dimensions:')
        for i in range(2):
            cellDimensions[i]=self.boundingBox[1][i]-self.boundingBox[0][i]
            self.logfile.write(" "+str(cellDimensions[i]))

        self.logfile.write('\ncell numbers:')
        for i in range(2):
            self.cellNumber[i]=int(cellDimensions[i]/self.cellSize)+1
            self.logfile.write(" "+str(self.cellNumber[i]))
        
        self.logfile.write(' cells: %d\n' % (
                self.cellNumber[0]*self.cellNumber[1]))

        # create empty array
        for i in range(self.cellNumber[0]):
            if i==0: self.cells=[]
            for j in range(self.cellNumber[1]):
                if j==0: LJ=[]
                LJ.append([])
            self.cells.append(LJ)

        # tris with a too small projected size will be ignored
        zerosize = (self.tolerance*avgSize*0.1)**2

        # map elements
        for element in self.triNodes:

            # ignore too small tris
            if abs(elSignedSize[element]) < zerosize:
                continue

            # get index ranges
            cellIndex=[[0,0],[0,0]] # low,high
            for l in range(2): # low,high - index
                for i in range(2): # dimension (x,y) - index
                    cellIndex[l][i]=int((elementBox[element][l][i]
                                         -self.cellOrigin[i])/self.cellSize)
            # append element ID to cells that overlap the element's bounding box
            for i in range(cellIndex[0][0],cellIndex[1][0]+1):
                for j in range(cellIndex[0][1],cellIndex[1][1]+1):
                    self.cells[i][j].append(element)

        # for smooth mode store radius
        self.smoothmode = smooth_findradius>0
        self.smooth_findradius = float(smooth_findradius)

        # finished initializing
        self.point_search_initialized = True
        return

    def point_is_inside(self, point):
        """test if point is inside the volume
        point is tuple of coordinates
        """

        if len(self.triNodes)==0:
            return False

        if not self.point_search_initialized:
            self.initialize_point_search()

        if (point[2]<self.zminmax[0]
            or point[2]>self.zminmax[1]):
            # point outside z-region
            return False

        # cell check (includes bounding box test from cell index)
        index=[0,0]
        for i in range(2):
            index[i]=int((point[i]-self.cellOrigin[i])/self.cellSize)
            if index[i]<0 or index[i]>=self.cellNumber[i]:
                return False

        # if cell check passed ... go on
        elems_in_cell=self.cells[index[0]][index[1]]
        centroid_to_point_distance=[]
        # sort by distance
        for element in elems_in_cell:
            l=self.dist_2D(point,self.triCentroids[element])
            centroid_to_point_distance.append([l,element])
        centroid_to_point_distance.sort()

        # full check if point inside any tri of elems_in_cell
        inside=0

        if self.smoothmode:
            if (centroid_to_point_distance and
                (centroid_to_point_distance[0][0]<self.smooth_findradius)):
                inside = 1

        for [l,element] in centroid_to_point_distance:
            node_pts = [self.coords[node]
                        for node in self.triNodes[element]]
            # vectors
            v12=vector(node_pts[0], node_pts[1])
            v13=vector(node_pts[0], node_pts[2])
            detA = v12[0]*v13[1] - v12[1]*v13[0]
            adjA = [[v13[1], -v13[0]], [-v12[1], v12[0]]]
            invA = [[x/detA for x in row] for row in adjA]

            v1p=vector(node_pts[0], point)
            a0 = invA[0][0]*v1p[0] + invA[0][1]*v1p[1]
            a1 = invA[1][0]*v1p[0] + invA[1][1]*v1p[1]

#             if a0>=0 and a1>=0 and (a0+a1)<=1:
#                 inside = 1
#                 break

            if ((a0 > -self.tolerance) and
                (a1 > -self.tolerance) and
                ((a0+a1) < (1+self.tolerance))):
                inside = 1
                break

        return inside


    ## ------------------------------------------------------------------------
    ## ------------------------------------------------------------------------
    ## transforming the volume

    def translate(self, move_vector):
        """Translating is done by simply moving the point coordinates for this
        volume. Other volumes or the original mesh from which the volume
        originates are not affected.

        This creates a new node-dictionary, which is refered to from now on.
        It contains only those nodes connected to any element in self.elset.

        The point search initialization is modified and doesn't have to be done
        again.
        """

        all_nodes = set()
        for nodes in self.triNodes.itervalues():
            all_nodes.update(nodes)
        
        new_coords = dict()
        for node_id in all_nodes:
            new_coords[node_id] = vector_plus(self.coords[node_id],move_vector)
        
        self.coords = new_coords
        for i in range(2): # min, max
            self.zminmax[i] += move_vector[2]


        # adapt point search properties
        if self.point_search_initialized and len(self.triNodes)>0:
            
            # tet Centroids
            for centroid in self.triCentroids.itervalues():
                vector_modif_add(centroid, move_vector)

            # BoundingBox
            for i in range(2): # x,y
                self.boundingBox[0][i] += move_vector[i]
                self.boundingBox[1][i] += move_vector[i]

            # search cells: nothing to modify here
            # self.cellOrigin is just an alias for self.boundingBox[0]


    ## ------------------------------------------------------------------------
    ## some service functions: 2D distance

    @staticmethod
    def dist_2D(a,b):
        return sqrt((b[0]-a[0])**2 + (b[1]-a[1])**2)



###################################################################
###################################################################
###################################################################
###################################################################


class VolumeBox(VolumeBase):
    """volume class

    You may supply an open logfile for diagnostic output or None to suppress
    it to the constructor.
    """

    def __init__(self, box, logfile = stdout):
        """if logfile == None, no diagnostic output
        
        box = [[xmin, ymin, zmin], [xmax, ymax, zmax]]
        """
        if logfile == None:
            self.logfile = QuietLog()
        else:
            self.logfile = logfile
        
        # internal variable, use getboundigbox method! (initialize it)
        self.boundingBox = box

        return

    def copy(self):
        """make a (shallow) copy of this volume
        """
        return VolumeBox(
            copy.deepcopy(self.boundigBox),
            logfile=self.logfile)


    def getboundingbox(self):
        """Bounding Box of the whole volume.
        Returns a list of two coordinate tuples (rather lists)
        The first coord tuple states the min value, the last the max.
        I.e. self.getboundingbox()[i][j]:
        i \in {0,1} for min/max value, j \in {0,1,2} coordinate index
        """
        return self.boundingBox

    def get_exterior_surface(self):
        """return a TriangleSurface object with all outer triangles of this
        volume.
        """

        raise Exception("VolumeBox.get_exterior_surface() not implemented yet.")

    def point_is_inside(self, point):
        """test if point is inside the volume
        point is tuple of coordinates
        """
        inside = True
        for i in range(3):
            if (point[i]<self.boundingBox[0][i]
                or point[i]>self.boundingBox[1][i]):
                inside = False
                break
        return inside

    def scale(self, scale_factor, scale_origin):
        """Scale the volume relative to scale_origin by scale_factor.
        """
        for i in range(2):
            orig2old = vector(scale_origin, self.boundingBox[i])
            orig2new = vector_scale(orig2old, scale_factor)
            self.boundingBox[i] = vector_plus(scale_origin, orig2new)

    def translate(self, move_vector):
        """Translate the volume by move_vector.
        """
        for i in range(2):
            self.boundingBox[i] = vector_plus(boundingBox[i], move_vector)
