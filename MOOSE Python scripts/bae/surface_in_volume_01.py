"""surface_in_volume_01.py

class for a surface of triangles embedded in a volume of tets
"""


_version_ = "1.03"
_version_history_ = """\
Versions:
=========

1.01 new   GP: new from surfToContact 2.10
1.02 GP added assignAllTetsGeom(),
        added repairSingleTriHoles(), moved it from surfToContact2.12.py
        added assignAllTetsGeom and renamed assignAllTets to assignAllTetsConn
1.03 GP found that assignAllTetsGeom doesn't work
        fixed assignAllTetsConn determining if face is connected to the border
        treat nodes in self.moebiusStopEdges like freeSurfBorderNodes
        changed assignAllTetsConn: stop looking for neighbouring tets to assign
        to the same side from a tet connected to freeSurfBorderNodes (this may
        be too harsh, let's see)
        if assignAdjacentTets wants to assign the same tet to both sides, new
        heuristics to find the better choice.
"""

## ---------------------------------------------------------------------------
_bugs_ = """\
known bugs/problems:
====================
- needs python 2.4, if in doubt use abaqus python
- method updateSplitNodesOtherTris() not working, because there is no need yet

ideas / todo:
=============

"""

from collections import deque  # deque = double ended queue
import sys

from bae.future_01 import *
from bae.vecmath_01 import *

from bae.surface_01 import TriangleSurface



class LengthOfSetNotOneException(Exception):
    pass
def getOnlyItem(this_set):
    if len(this_set)!=1:
        raise LengthOfSetNotOneException(
            "Set has %d items, wheras getOnlyItem() expects exactly one."
            % len(this_set))
    return tuple(this_set)[0]



## ---------------------------------------------------------------------------
## class for surface in volume

class TriangleSurfaceInVolume(TriangleSurface):
    """adds some functionality concerning the attached tets
    to the TriangleSurface class

    There are two algorithms used for seperation of a tet mesh at this surface:
    ===========================================================================

    

    A. connection based algorithm
    -----------------------------
    
    consists of three steps:
    1. assignAdjacentTets()  ... build dictionaries and find tets directly
          adjacent to surface
    2. findFreeBorder()      ... investigate the border of the surface
    3. assignAllTetsConn()   ... assign all other tets to either side
    4. updateSplitNodesTets(splitNodes, blackList)  ... update list of nodes of
          the mesh to be split at this surface with splitting direction
          

    B. geometry base algorithm (does not work properly!)
    --------------------------
    1. assignAdjacentTets()  ... build dictionaries and find tets directly
          adjacent to surface
    2. findFreeBorder()  ... investigate the border of the surface
    3. assignAllTetsGeom()   ... assign all other tets to either side
    4. updateSplitNodesTets(splitNodes, blackList)  ... update list of nodes of
          the mesh to be split at this surface with splitting direction


    pros and cons A vs. B:
    ----------------------
    B does not work properly!
    A works (if at all) no matter how sharp edges may be there.


    =====================================================================
    Note to the programmer, if you modify this class:
    The methods of this child class should behave as if its base class members
    where members of an external object and therefore only use the
    TriangleSurface methods that are meant as the external interface of
    TriangleSurface. Why? To keep it simple.
    """

    def __init__(self, surfaceName, triList, triNodes, nodeCoords,
                 tetNodes,
                 logfile=sys.stdout, initchecks = True):
        """see class description of TriangleSurface

        tetNodes is a dict {elementID: nodeIDs[4]} of all tets in the mesh
        the TriangleSurfaceInVolume object keeps a reference of tetNodes
        """
        TriangleSurface.__init__(self, surfaceName, triList, triNodes,
                                 nodeCoords, logfile, initchecks)

        # create the connectivity of the surface to the volume
        # self.tetNodes {tetID:corner nodes} only tets with nodes on the surf
        # self.tetSurfNodes {tet: [nodes on surface]}
        #                              where len([nodes on surf])>0
        # self.nodeTets  {node:[tetIDs,]} for tets in self.tetSurfNodes
        self.tetNodes = dict()
        self.tetSurfNodes = dict()
        self.nodeTets = defaultdict(list)

        # loop over all tet elems in the mesh
        for element, nodes in tetNodes.iteritems():
            if len(nodes)<4:
                continue  # this is not a tet
            nodesOnSurf = self.getNodes().intersection(nodes)
            if len(nodesOnSurf)==0:
                continue  # this tet is not connected to the surface
            # self.tetSurfNodes
            self.tetSurfNodes[element] = nodesOnSurf
            # self.NodeTets
            for node in nodesOnSurf:
                self.nodeTets[node].append(element)
            # self.tetNodes
            self.tetNodes[element] = nodes[0:4]

        # self.cuttingFaces will be updated (if nessecary) by addCuttingSurf()
        self.cuttingFaces = set()


    def messageWarn(self, warnString):
        self.logfile.write( "WARNING: " + warnString + '\n')


    def addCuttingSurf(self, othersurface):
        """

        updates self.cuttingFaces
        """

        for otherNodes in othersurface.triNodes.itervalues():
            nodesOnSurf = self.getNodes().intersection(otherNodes)
            if len(nodesOnSurf)>0:
                self.cuttingFaces.add(frozenset(otherNodes))
        return


    #---- some check and repair functions
    def findSingleTriHoles(self):
        """Returns a list of node sets, each set defining a missing tri.

        Missing tris in that sense are three nodes connected only by border
        edges and not already connected by a tri element.

        See also TriangleSurface.findSingleTris(), it finds tri elements with
        only border edges."""

        if not(hasattr(self, "borderEdges")):
            self.calcEdges()

        newTris = set()
        for element, nodesOnSurf in self.tetSurfNodes.iteritems():
            nodesOnBorder = self.getBorderNodes().intersection(nodesOnSurf)
            if len(nodesOnBorder) == 4:
                thisFacesIter = self.tetFacesIter(nodesOnBorder)
            elif len(nodesOnBorder) == 3:
                thisFacesIter = (frozenset(nodesOnBorder), )
            else:
                # this tet does not share a face with the surface
                continue
            for faceNodes in thisFacesIter:
                if faceNodes in newTris:
                    continue
                if ((len(self.borderEdges.intersection(
                                self.triEdgesIter(faceNodes))) == 3)
                    and (faceNodes not in self.getNodesToTri())):
                    # all edges of this face are in self.borderEdges
                    # and there is no tri connecting those edges
                    newTris.add(faceNodes)

        # finished
        return newTris

    def repairSingleTriHoles(self, elemOffset=None):
        """

        returns a tuple:
        1. elemOffset: number of the last inserted triangle element
        2. newTriList: list of all newly inserted triangle element
        """
        if not(elemOffset):
            # just take the last of all known elements
            elemOffset = max(max(self.triNodes), max(self.tetNodes))

        singleTriHoles = self.findSingleTriHoles()
        # look for double edges in this set
        singleTriHolesEdgeToTris = defaultdict(set)
        for thisTriNodes in singleTriHoles:
            for thisEdge in self.triEdgesIter(thisTriNodes):
                singleTriHolesEdgeToTris[thisEdge].add(thisTriNodes)
        # remove all but one new tris from singleTriHoles
        # the best has the highest rating, the smallest angle to neighbours
        for dset in singleTriHolesEdgeToTris.itervalues():
            if len(dset)<=1:   # should never be <1
                continue
            besttri = None
            for thisTriNodes in dset:
                rating = length(self.calcAngleToNeighbours(
                        thisTriNodes, default=0.0))
                if not(besttri):
                    besttri = (rating, thisTriNodes)
                    continue
                if rating>besttri[0]:
                    singleTriHoles.remove(besttri[1])
                    besttri = (rating, thisTriNodes)
                else:
                    singleTriHoles.remove(thisTriNodes)
        # now add all what is left in singleTriHoles to the surface
        newTriList = list()
        for newNodes in singleTriHoles:

            elemOffset += 1
            newTri = elemOffset
            self.insertTri(newTri, newNodes)
            newTriList.append(newTri)

        return (elemOffset, newTriList)
        

    def findCorruptedBorderTris(self):
        """Find triangles on the border that better seem to be deleted.

        If there is a tet with all four nodes on the surface and one of the
        tris lying on the tets faces is connected to three border nodes (at
        least two border edges), that tri is considered a tri that currupts
        the border and should be deleted.
        """

        lostTris = set()
        for element, nodesOnSurf in self.tetSurfNodes.iteritems():
            if len(nodesOnSurf) < 4:
                continue
            nodesOnBorder = self.getBorderNodes().intersection(nodesOnSurf)
            if len(nodesOnBorder) == 4:
                thisFacesIter = self.tetFacesIter(nodesOnBorder)
            elif len(nodesOnBorder) == 3:
                thisFacesIter = (frozenset(nodesOnBorder), )
            else:
                # this tet does not share a face with the surface
                # this is a very questionable case too (a tet connected to the
                # surface by all four nodes but not by a single face!), but we
                # ignore it
                continue
            for nodesOfThisTri in thisFacesIter:
                try:
                    triId = self.getNodesToTri()[nodesOfThisTri]
                except KeyError:
                    continue
                lostTris.add(triId)

        # finished
        return lostTris


    # some iterator function
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
        (see getTriNormalUsc) points inwardly.
        """
        nodesOfThisTet = tuple(nodesOfThisTet)
        for i1,i2,i3 in ((0,1,2), (1,0,3), (2,1,3), (0,2,3)):
            yield frozenset((nodesOfThisTet[i1], nodesOfThisTet[i2],
                             nodesOfThisTet[i3]))

    @staticmethod
    def tetEdgesIter(nodesOfThisTet):
        """Returns an iterator of the edges of the given tet element
        nodesOfThisTet will be converted to a tuple of node ids

        use as in
        myTet = mdb.el_node[5731]  # list of node ids, e.g. [3, 43, 12, 36]
        for edges in tetEdgesIter(myTet):
           ... edges will be frozenset((3, 43, )), next time
               frozenset((43, 3, 36)) and so on.

        An extra comment on the index ordering used in the for loop below: If
        the points order would stay that of the used tuple, each faces normal)
        (see getTriNormalUsc points inwardly.
        """
        nodesOfThisTet = tuple(nodesOfThisTet)
        for i1,i2 in ((0,1), (0,2), (0,3), (1,2), (1,3), (2,3)):
            yield frozenset((nodesOfThisTet[i1], nodesOfThisTet[i2]))


    ## ------------------------------------------------------------------------
    ##    seperate mesh at surface (connection based algorithm)
    ##
    def assignAdjacentTets(self):
        """first step of splitting the mesh

        1. build dictionaries
        self.tetToFace       {surface-Tet-elementID: set(node1, node2, node3)} 
        self.faceToTet       {set(node1, node2, node3): list of elements}
        surfTetsCentroids    {element:centroid[3]} for surf tets


        2. find tets directly adjacent to surface (sharing more than two nodes
           with the surface)

        
        initializes:

        self.assignedTets [set(), set()], set of tets on neg side, pos side
        self.assignedTetFaces  [neg side list, pos side list], contains only
           'early assigned' tets and their faces not on this surface.
           'early assigned' tets are those tets that share three nodes with the
           surface, their assignment to either side at this point is 'earlier'
           than those of the other tets to follow
            - neg/pos side lists are lists cointaining tuples (tet elem, (free
              faces, border flag))
            - the border flag is True, if the edge connecting the free face
              and the corresponding surface tri is on the border
        self.assignedTetFaces [list(), list()],
            each list has tuples (tet elem, free faces)
        self.cohesiveEls { (negElement,posElement) -> 
            [[nodeindices of matching nodes in negElement], [same in posEl.]] }
        """

        # self.tetToFace: {surface-Tet-elementID: set(node1, node2, node3)}
        self.tetToFace = dict()
        # self.faceToTet: {set(node1, node2, node3): list of elements}
        self.faceToTet = defaultdict(list)
        # surfTetsCentroids: {element:centroid[3]} for surf tets
        surfTetsCentroids = dict()
        #
        for element, thisTetNodes in self.tetNodes.iteritems():

            # tetToFace, faceToTet
            self.tetToFace[element] = set()
            for thisFace in self.tetFacesIter(thisTetNodes):
                self.tetToFace[element].add(thisFace)
                self.faceToTet[thisFace].append(element)
            # surfTetsCentroids
            x=y=z=0
            for node in thisTetNodes:
                coords=self.nodeCoords[node]
                x=x+coords[0]/4.0
                y=y+coords[1]/4.0
                z=z+coords[2]/4.0
            surfTetsCentroids[element]=[x,y,z]


        # cohesiveEls: { (negElement,posElement) -> 
        #   [[nodeindices of matching nodes in negElement], [same in posEl.]] }
        # tet pairs with node indices to create cohesive elements
        # node indices follow separation from intersections, node numbers don't
        self.cohesiveEls=dict()

        # assignedTets[0]: set of tets on neg side,
        # assignedTets[1]: ... on pos side
        self.assignedTets = [set(), set()]

        # assignedTetFaces[i]: lists of tuples (tet elem, (free face, onBorder))
        # contains only elements with a face on the splitting surface
        # a free face is a Tet-Faces where surface tets may be connected to,
        # that have not been assigned to either side yet
        # i=0: neg side, i=1: pos side
        self.assignedTetFaces = [list(), list()]

        # find tets directly adjacent to surface (sharing three or four nodes
        # with the surface)
        # update assignedTets, cohesiveEls
        wrongNumOfTetsOnTri = dict()
        twoTetsOnSameSide = list()
        tetsOnBothSide = set()
        for triElem, thisSurfFace in self.triNodes.iteritems():
            oneNode = thisSurfFace[0]  # any one node on the surface triangle
            thisSurfFace = frozenset(thisSurfFace)
            try:
                tetsOnThisTri = self.faceToTet[thisSurfFace]
            except KeyError:
                # no tet connected to this tri
                continue  # no tet connected to this tri

            # each face of the surface should be connected to two tets
            # where not, the tri is being ignored
            if len(tetsOnThisTri) != 2:
                wrongNumOfTetsOnTri[triElem]=tetsOnThisTri
                continue


            elementPair = [None, None]
            nodeIndexList = [list(), list()]
            for tetElem in tetsOnThisTri:
                vectorNodeCent = vector(self.nodeCoords[oneNode],
                                        surfTetsCentroids[tetElem])
                side = dot(self.triNormals[triElem], vectorNodeCent)
                if side == 0:
                    self.messageWarn("Tet element %d has zero thickness."
                                     % tetElem)
                if side>0.:
                    sideidx = 1
                else:
                    sideidx = 0

                # face on surface is not a free face
                # have to compare to all faces (nodes) of this surface
                # not just 
                freeFaces = self.tetToFace[tetElem].difference(
                    self.getNodesToTri())

                # is the edge connecting freeFaces to the surface on the border?
                freeFaceTups = list()
                for thisFreeFace in freeFaces:
                    thisSurfEdge = thisFreeFace.intersection(thisSurfFace)
                    borderFlag = thisSurfEdge in self.borderEdges
                    freeFaceTups.append((thisFreeFace, borderFlag))

                # add new faces
                if tetElem in self.assignedTets[1-sideidx]:
                    # tet already on other side!
                    tetsOnBothSide.add(tetElem)
                    self.assignedTets[1-sideidx].remove(tetElem)

                    # find average normal of all its nodes connected to surface
                    tetNodesOnSurf = set(self.tetNodes[tetElem]).intersection(
                        self.getNodeTris())
                    avgNormal = [0,0,0]
                    for thisNode in tetNodesOnSurf:
                        avgNormal = vector_plus(avgNormal,
                                                self.getNodeNormals(thisNode))

                    # find common centroid of surface tris adjacent to this tet
                    avgCentroid = [0,0,0]
                    avgCentroidCnt = 0
                    for thisFace in self.tetFacesIter(self.tetNodes[tetElem]):
                        if thisFace not in self.getNodesToTri(): continue
                        for thisNode in thisFace:
                            avgCentroid = vector_plus(
                                avgCentroid, self.nodeCoords[thisNode])
                        avgCentroidCnt += 3
                    avgCentroid = vector_scale(avgCentroid, 1./avgCentroidCnt)

                    # vec from avgCentroid to centroid of tet
                    # correct sideidx
                    newside = dot(
                        vector(avgCentroid, surfTetsCentroids[tetElem]),
                        avgNormal)
                    if newside>0.:
                        sideidx = 1
                    else:
                        sideidx = 0
                # any way, assign to the correct (whatever that is) side
                self.assignedTets[sideidx].add(tetElem)
                self.assignedTetFaces[sideidx].append((tetElem, freeFaceTups))
                elementPair[sideidx] = tetElem
                nodeIndexList[sideidx]=[self.tetNodes[tetElem].index(coh_node)
                                        for coh_node in self.triNodes[triElem]]


            # each face of the surface should have one tet on each side
            # if not, this tri is being ignored
            if ((elementPair[0]==None or elementPair[1]==None) and
                (triElem not in wrongNumOfTetsOnTri)):
                twoTetsOnSameSide.append((triElem, tetsOnThisTri))


            # update cohesiveEls-dict
            # cohesiveEls: { (negElement,posElement) -> 
            #    [[nodeindices of matching nodes in negElement],
            #     [same in posEl.]] }
            elementPair = tuple(elementPair)
            if elementPair not in self.cohesiveEls:
                self.cohesiveEls[elementPair]=nodeIndexList

            if (elementPair[1],elementPair[0]) in self.cohesiveEls:
                raise Exception, 'double cohesive element'

        return (wrongNumOfTetsOnTri, twoTetsOnSameSide, tetsOnBothSide)



    def findFreeBorder(self):
        """Investigate the border of the surface.
        Find free borders, where surf tet neighbours connect upper to
        lower side. Those nodes will not be split.

        Add also all nodes in self.moebiusStopEdges to freeSurfBorderNodes
        because those nodes shall not be split either and shall be treated the
        same.

        initializes:
        self.freeSurfBorderNodes   set of surf nodes where surf ends within
                                   model

        uses:
        self.assignedTets with only early assigned tets
                  ... as from self.assignAdjacentTets()

        returns a tuple of error cases:

        borderTetMissing: list of tri ids
           triangles on the border of the surface not having tets on both sides
           That error may also occur if two splitting surfaces partly coincide.

        borderTetDouble: list of tuples (edge, side_idx, tets on this side)
           triangle on the border of the surface having more than one 'early
           assigned' tets (sharing three nodes with the surface) on either
           side. Those tets must be intersecting each other.
           This error may also occur when the tets on intersecting surfaces
           where heavily distorted when the mesh was split by a too large
           amount.

        wrappedTets: list of tet elements
           tets on the border of the surface which are connected to the surface
           with all their nodes.

        danglingTets: list of tet elements
           tets on the border of the surface that do not share a face (three
           nodes) with the surface. That is a contradiction because 'early
           assigned' tets are identified by sharing three nodes with the
           surface.
        
        """

        if not(hasattr(self,"borderEdges")) or not(hasattr(self,"edgeToTri")):
            self.calcEdges()
        
        if not hasattr(self, "assignedTets"):
            # you should rather call assignAdjacentTets before to get
            # diagnostic data and create some output
            self.assignAdjacentTets()

        # dict {edges: list of connected tetElems} ( edge=set(node1,node2) )
        borderEdgeToTet = dict()
        for thisEdge in self.borderEdges:
            tetlist = [tetElem
                       for tetElem in self.tetNodes
                       if len(thisEdge.intersection(self.tetNodes[tetElem]))==2]
            borderEdgeToTet[thisEdge] = tetlist


        # In the following for loop: Try to walk from pos to neg side within
        # borderEdgeToTet[thisEdge] list. Yields freeSurfBorderNodes: set of
        # nodes on the surface border where the surface ends within the model,
        # not at its border or at another surface already splitting the model.
        # Those points will not be split by the current surface.
        self.freeSurfBorderNodes = set()

        # error flags. for explanation see warning messages below
        borderTetMissing = list()
        borderTetDouble = list()
        danglingTets = list()
        wrappedTets = list()
        for thisEdge in self.borderEdges:

            # find assignedTetsOnThisEdge[side_idx]: tets already assigned to
            # either side of the surface:  'early assigned' tets
            # they are those tets that share three nodes with the surface
            edge_ok = True
            assignedTetsOnThisEdge = [0,0]
            for side_idx in range(2):
                assignedTetsOnThisEdge[side_idx] = tuple(
                    set(borderEdgeToTet[thisEdge]).intersection(
                        self.assignedTets[side_idx]))
                if len(assignedTetsOnThisEdge[side_idx])<1:
                    borderTetMissing.append(tuple(self.edgeToTri[thisEdge])[0])
                    # for security reasons add edge nodes to freeSurfBorderNodes
                    self.freeSurfBorderNodes.update(thisEdge)
                    edge_ok = False
                elif len(assignedTetsOnThisEdge[side_idx])>1:
                    borderTetDouble.append(
                        (thisEdge, side_idx, assignedTetsOnThisEdge[side_idx]))
                    # add edge nodes to self.freeSurfBorderNodes
                    # (for security reasons)
                    self.freeSurfBorderNodes.update(thisEdge)
                    edge_ok = False
                else:
                    assignedTetsOnThisEdge[side_idx] = \
                        assignedTetsOnThisEdge[side_idx][0]

            # initialize start at pos side
            # thisTetList contains all tets on the border edge. They are tested
            # for connecting both sides of the surface.
            # nextelem is (initially) the tet connected to pos side
            # nextnode is supposed to be the node of nextelem not on the border
            # edge and connecting nextelem (the start tet) to the next on the
            # walk to the other side
            if edge_ok:
                thisTetList = list(borderEdgeToTet[thisEdge])
                nextelem = assignedTetsOnThisEdge[1]
                thisTetList.remove(nextelem)
                nextnode = set(self.tetNodes[nextelem]).difference(
                    self.getNodeTris().iterkeys())
                if len(nextnode)>1:
                    # more than one nodes of the starting tet are not on this
                    # surf for security reasons add edge nodes to
                    # self.freeSurfBorderNodes
                    danglingTets.append(nextelem)
                    self.freeSurfBorderNodes.update(thisEdge)
                    edge_ok = False
                elif len(nextnode)<1:
                    # all nodes of the starting tet are on this surf
                    # for security reasons add edge nodes to
                    # self.freeSurfBorderNodes
                    wrappedTets.append(nextelem)
                    self.freeSurfBorderNodes.update(thisEdge)
                    edge_ok = False
                else:
                    nextnode = nextnode.pop()

            # Now do the actual walk from pos side to neg side
            if edge_ok:
                isBorderEdge = True
                while nextelem != assignedTetsOnThisEdge[0]:
                    # check if another surface intersects at this face!
                    checkFace = set(thisEdge)
                    checkFace.add(nextnode)
                    checkFace = frozenset(checkFace)
                    if checkFace in self.cuttingFaces:
                        # yes: this face belongs to another surface
                        # => no connection from pos side to neg side
                        # => this edge is no free borderedge
                        isBorderEdge = False
                        break

                    # try to find the next step
                    lastelem = nextelem
                    lastnode = nextnode
                    nextelem = [t for t in thisTetList
                                if lastnode in self.tetNodes[t]]
                    if len(nextelem) == 0:
                        # no next element, end of the walk.
                        # => no connection from pos side to neg side
                        # => this edge is no free borderedge
                        isBorderEdge = False
                        break
                    if len(nextelem) > 1:
                        raise Exception(
                            'more than two edge-tets on non-edge node')
                    nextelem = nextelem[0]
                    thisTetList.remove(nextelem)
                    nextnodeset = set(
                        self.tetNodes[nextelem]).difference(thisEdge)
                    nextnodeset.remove(lastnode)
                    nextnode = getOnlyItem(nextnodeset)
                    continue
                if isBorderEdge:
                    # add edge nodes to self.freeSurfBorderNodes
                    self.freeSurfBorderNodes.update(thisEdge)

        # Add all nodes in self.moebiusStopEdges to freeSurfBorderNodes
        for thisEdge in self.moebiusStopEdges:
            self.freeSurfBorderNodes.update(thisEdge)

        # fini
        return (borderTetMissing, borderTetDouble, danglingTets, wrappedTets)


    def assignAllTetsConn(self):
        """separate tets at surface, assign all other tets to either side.

        connection based algorithm:
        First find all tets connected to surface other than free surface border
        then assign tets to either side by looking at neighbours

        update self.assignedTets

        return:
        - moreThanTwoTets: list of tuples (set of three nodes, set of attached
          tets) for faces with more than two attached tets) Those faces will be
          ignored.
        - notAssignedTets: set of tets (still) not assigned to either side
        - tetsOnBothSide: set of tets which were found on both sides of the
          surface. They are only assigned to the neg side.

        some details:
        - collect all tet elements yet to be assigned to either side in
          notAssignedTets. There are no tets in there which are only connected
          to the surface by free border nodes. But there are still some tets in
          it that don't go on either side: those only connected to border nodes
          of the surface (need not be free border!) and behind another inter-
          secting surface. They will be identified during the walk from already
          assigned tets to its neighbours by this check:
            if onBorder and thisFace in self.cuttingFaces...
        - The walk from already assigned tets to its neighbours stops when the
          already assigned tet to start from is connected to the free border.
          This is to prevend walking to the other side. All tets yet to be
          assigned to either side are connected to at least one neighbour that
          is not connected to a free border at all. (I hope so.)
        """

        if not hasattr(self, "freeSurfBorderNodes"):
            # you should rather call findFreeBorder before to get diagnostic
            # data and create some output
            self.findFreeBorder()


        # notAssignedTets: Tets yet to be assigned to either side
        notAssignedTets = set()
        for element, nodesOfThisTet in self.tetNodes.iteritems():
            nodesOnSurf = [node
                           for node in nodesOfThisTet
                           if node in self.getNodes()
                           and node not in self.freeSurfBorderNodes]
            if len(nodesOnSurf)>0:
                notAssignedTets.add(element)

        # some are already assigned (those sharing three nodes with the surface)
        notAssignedTets -= self.assignedTets[0]
        notAssignedTets -= self.assignedTets[1]


        # find neighbours of already assigned tets
        # update assignedTets, notAssignedTets
        # 
        # Tet element lastTet has not yet investigated faces given in the list
        # neighbourFaces. For each of those faces find out if there is a tet
        # behind we should assign to the current side. (lastTet is in
        # assignedTets[sideidx] already!)
        #
        # A neigbouring tet is considered not yet assigned to either side if
        # it's still in notAssignedTets. It is then moved to
        # assignedTets[sideidx]. Then all of the neighbouring tets faces not
        # connected to lastTet are recursivly investigated.
        # 
        moreThanTwoTets = list()
        tetsOnBothSide = set()
        for sideidx in range(2):
            for lastTet, neighbourFaces in self.assignedTetFaces[sideidx]:
                tetsToFindNeighbours = deque()
                tetsToFindNeighbours.append((lastTet, neighbourFaces))
                while len(tetsToFindNeighbours) > 0:

                    lastTet, neighbourFaces = tetsToFindNeighbours.popleft()
                    for thisFace, onBorder in neighbourFaces:

                        # check if another surface intersects at this face!
                        # if on border stop at this face
                        if onBorder and thisFace in self.cuttingFaces:
                            continue

                        # find neighbour tet on the other side of thisFace
                        newTets = list(self.faceToTet[thisFace])

                        if len(newTets)>2:
                            moreThanTwoTets.append((thisFace, newTets))
                        newTets.remove(lastTet)
                        if len(newTets)!=1:
                            continue # no other tets on this face or error
                        newTet = newTets[0]

                        if newTet not in notAssignedTets:
                            continue

                        # if on border don't proceed to new surface nodes, collect only
                        # tets connected to the old nodes on the surface
                        if onBorder:
                            facesNodesOnSurf = thisFace.intersection(
                                self.getNodeTris())
                            newTetSurfNodes = set(self.tetNodes[newTet]).intersection(
                                self.getNodeTris())
                            if len(newTetSurfNodes-facesNodesOnSurf)>0:
                                continue

                        # ok, found (valid) neighbour tet: newTet
                        notAssignedTets.remove(newTet)
                        if newTet in self.assignedTets[1-sideidx]:
                            # already on other side!
                            tetsOnBothSide.add(newTet)
                        else:
                            self.assignedTets[sideidx].add(newTet)


                        # don't step further from this tet if it is
                        # connected to freeSurfBorderNodes
                        if set(self.tetNodes[newTet]).intersection(
                            self.freeSurfBorderNodes):
                            continue

                        # find free faces and recursivly find neighbouring tets
                        # is one of newFaces connected to the border?
                        newFaces = set(self.tetToFace[newTet])
                        newFaces.remove(thisFace)

                        newFaceTups = list()
                        for thisFreeFace in newFaces:
                            facesNodesOnSurf = thisFreeFace.intersection(
                                self.getNodeTris())
                            borderFlag = all([node in self.getBorderNodes()
                                              for node in facesNodesOnSurf])
                            newFaceTups.append((thisFreeFace, borderFlag))

                        tetsToFindNeighbours.append((newTet, newFaceTups))

        # fini
        return (moreThanTwoTets, notAssignedTets, tetsOnBothSide)


    def assignAllTetsGeom(self):
        """separate tets at surface, assign all other tets to either side.

        This method does not work properly as it was intended. It does not
        recognize if a tet connected to a node on the border of the surface is
        behind another surface and has to be discarded from the set of tets to
        be assigned to either side of this surface. Nevertheless it's left
        here, because it may be of interest some day.

        Geometry based algorithm:
        Assign not yet assigned tets to either side of the surface by comparing
        the vector from the tets nodes on the surface to its centroid with the
        average node normal of the connected surface nodes.

        'Not yet assigned tets' are tets with exactly one or two nodes on the
        surface (those with three or four are already in self.assignedTets)
        and at least one of those nodes on the surface is not in the set
        self.freeSurfBorderNodes.

        updates self.assignedTets
        """

        if not hasattr(self, "freeSurfBorderNodes"):
            # This should be superfluous:
            # You should rather call findFreeBorder before, to generate some
            # diagnostic output from its return values (which are ignored here)
            self.findFreeBorder()

        # notAssignedTets: Tets yet to be assigned to either side
        notAssignedTets = set()
        for element, nodesOfThisTet in self.tetNodes.iteritems():
            nodesOnSurf = [node
                           for node in nodesOfThisTet
                           if node in self.getNodes()
                           and node not in self.freeSurfBorderNodes]
            if len(nodesOnSurf)>0:
                notAssignedTets.add(element)

        # some are already assigned (those sharing three nodes with the surface)
        notAssignedTets -= self.assignedTets[0]
        notAssignedTets -= self.assignedTets[1]

        # separate tets at surface / by average of node normals and centroid
        for thisTet in notAssignedTets:
            normalSum=[0,0,0]
            pointOnSurf=None  # average of at most two nodes on the surface
            x=y=z=0
            for i, node in enumerate(self.tetNodes[thisTet]):
                coords=self.nodeCoords[node]
                # calculation of centroid 
                x=x+coords[0]/4.0
                y=y+coords[1]/4.0
                z=z+coords[2]/4.0

                # surface normal and point on surface
                if node in self.getNodes(): # nodes on surface
                    normalSum=vector_plus(normalSum,self.getNodeNormals()[node])
                    if pointOnSurf:
                        pointOnSurf=vector_scale(
                            vector_plus(pointOnSurf, coords), 0.5)
                    else:
                        pointOnSurf=coords

            centroid=[x,y,z]
            vectorNodeCent=vector(pointOnSurf, centroid)

            side=dot(vectorNodeCent,normalSum)
            if side>0.:
                self.assignedTets[1].add(thisTet)
            else:
                self.assignedTets[0].add(thisTet)

        # fini
        return


    ## ------------------------------------------------------------------------
    ##   update splitNodes
    ##
    ##   with the splittings of the tets attached to this surface
    ##   insert a splitting-entry for each node on the surface
    ## ------------------------------------------------------------------------
    def updateSplitNodesTets(self, splitNodes, blackList):
        """
        updates:
        splitNodes a dict : {node id: [splitting 1, splitting 2, ...]}
            each splitting is a three-tuple (neg. side attached tets, pos. side
            attached tets, averaged splitting normal)

        uses:
        self.freeSurfBorderNodes, self.assignedTets, self.getNodeNormals()
        """

        for nodeA in self.getNodes():
            # check in case loose surf tri nodes have no attached tets (harpoon)
            if (nodeA not in self.nodeTets) or (nodeA in blackList):
                continue

            # nodeA on surface border -> don't seperate
            if nodeA in self.freeSurfBorderNodes:
                continue

            # find attached tets on either side
            attachedTets = [0,0]
            for sideidx in range(2):
                attachedTets[sideidx] =self.assignedTets[sideidx].intersection(
                    self.nodeTets[nodeA])

            # register a new splitting with the splitting normal
            splitNodes[nodeA].append((attachedTets[0], attachedTets[1],
                                      self.getNodeNormals()[nodeA]))



    ## ------------------------------------------------------------------------
    ##   update splitNodes
    ##   
    ##   insert a splitting-entry for each tri element in triNodes
    ## ------------------------------------------------------------------------
    def updateSplitNodesOtherTris(self, splitNodes, blackList, triNodes):
        """
        updates:
        splitNodes a dict : {node id: [splitting 1, splitting 2, ...]}
            each splitting is a three-tuple (neg. side attached tets, pos. side
            attached tets, averaged splitting normal)

        not working correctly, do it later....
        (if there is a need to also correctly split the surfaces)
        """

        coincideTris = list()
        for triElem, elemNodes in triNodes.iteritems():
            if triElem in self.triNodes:
                continue

            # check for surface nodes
            replaceNodesIdx=[(node, idx)
                             for idx, node in enumerate(elemNodes)
                             if node in self.getNodesTris()]
            if len(replaceNodesIdx)==0:
                continue

            freeNodes = list(elemNodes)
            for i,ii in replaceNodesIdx: freeNodes.remove(i)
            if len(freeNodes)==0:
                coincideTris.append(triElem)
                continue
            freeNode = freeNodes[0] # just any node not connected to the surface

            # test if freeNode belongs to any tet on the upper side
            sideidx = None  # which side yet undecided
            for tetElem in self.assignedTets[1]:
                if freeNode in self.tetNodes[tetElem]:
                    sideidx = 1 # upper/pos side
                    break
            if sideidx==None:
                for tetElem in self.assignedTets[0]:
                    if freeNode in self.tetNodes[tetElem]:
                        sideidx = 0 # lower/neg side
                        break
            if sideidx==None:
                messageWarn("tri element %d not connected to surface tets" % triElem)

#             # if tri contains surface nodes -> replace only if element on pos side
#             if sideidx==1:
#                 for node, idx in replaceNodesIdx:

#                     # update global tri element dict
#                     triNodes[triElem][idx] = mapOffsetNode[node]
#                     # update all surface-objects lists
#                     for otherSurf in surfaceList:
#                         if otherSurf==self: continue
#                         if triElem in otherSurf.triNodes:
#                             otherSurf.replaceNode(triElem, node, mapOffsetNode[node])




    def cleanSplittingVars(self):
        del self.assignedTets
        del self.assignedTetFaces
        del self.freeSurfBorderNodes
        del self.tetToFace
        del self.faceToTet

