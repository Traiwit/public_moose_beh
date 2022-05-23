"""surface_in_volume_02.py

class for a surface of triangles embedded in a volume of tets
"""


__version__ = "2.05"
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
2.01 GP moved from bae.surface_01 to bae.surface_02
2.02 GP fixed: two other modifications to the algorithm:
         - take into accont the free border of intersecting surfaces,
         - assignAllTetsConn now first looks for direct neighbours of all tets,
           that are already assigned to either side, then their neighbours and
           so on instead of first asigning one side completely then the other.
2.03 GP changed: messages to logfile and sets to a debugmodel for:
    repairSingleTriHoles, findCorruptedBorderTris
2.05 GP switched to bae.log_01
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
import copy

from bae.future_01 import *
from bae.vecmath_01 import *

from bae.surface_02 import TriangleSurface, checkModelModuleVersion
from bae.log_01 import msg


## ---------------------------------------------------------------------------
## some service stuff

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

    The main purpose of the class is to provide a means to seperate a tet mesh
    at this surface. It uses a connection based algorithm which consists of
    three steps:

     1. assignAdjacentTets()  ... build dictionaries and find tets directly
           adjacent to surface
     2. findFreeBorder()      ... investigate the border of the surface
     3. assignAllTetsConn()   ... assign all other tets to either side
     4. updateSplitNodesTets(splitNodes, blackList) ... update list of nodes of
           the mesh to be split at this surface with splitting direction
    """

#   =====================================================================
#   Note to the programmer, if you modify this class:
#   The methods of this child class should behave as if its base class members
#   where members of an external object and therefore only use the
#   TriangleSurface methods that are meant as the external interface of
#   TriangleSurface. Why? To keep it simple. Otherwise it would be much harder
#   to change the internals of either class.
#   =====================================================================

    def __init__(self, nodeCoords=None, triNodes=None, tetNodes=None,
                 logfile=sys.stdout):
        """see class description of TriangleSurface

        tetNodes (optional) is a dict {elementID: nodeIDs[4]} of at least those
        tets in the mesh that are connected to the triangles. There may be more
        tets in tetNodes. The TriangleSurfaceInVolume object works with a copy
        of only the relevant tets.

        Note: If you specify tetNodes at a later stage, make sure its values
        consist of only four nodes each (the corner nodes of the tets).
        """
        TriangleSurface.__init__(self, nodeCoords, triNodes)

        # create the connectivity of the surface to the volume
        # self.tetNodes {tetID:corner nodes} only tets with nodes on the surf
        # self.tetSurfNodes {tet: [nodes on surface]}
        #                              where len([nodes on surf])>0
        # self.nodeTets  {node:[tetIDs,]} for tets in self.tetSurfNodes
        self.tetNodes = dict()
        self.tetSurfNodes = dict()
        self.nodeTets = defaultdict(list)

        # loop over all tet elems in the mesh
        if tetNodes is not None:
            assert isinstance(tetNodes, dict)
            allSurfNodes = set(self.getNodeTris())
            for element, nodes in tetNodes.iteritems():
                if len(nodes)<4:
                    continue  # this is not a tet
                nodesOnSurf = allSurfNodes.intersection(nodes)
                if len(nodesOnSurf)==0:
                    continue  # this tet is not connected to the surface
                # self.tetSurfNodes
                self.tetSurfNodes[element] = nodesOnSurf
                # self.NodeTets
                for node in nodesOnSurf:
                    self.nodeTets[node].append(element)
                # self.tetNodes
                self.tetNodes[element] = nodes[0:4]

        # self.cuttingFaces and self.surfIntersectionNodes
        # will be updated (if nessecary) by addCuttingSurf()
        self.cuttingFaces = set()
        self.surfIntersectionNodes = set()

    def messageWarn(self, warnString):
        msg("WARNING: %s" % warnString)


    def addCuttingSurf(self, othersurface):
        """updates the following two attributes

         - self.cuttingFaces ... set of faces (each being a frozenset of node
           numbers) on other surfaces intersecting this surface
         - self.surfIntersectionNodes ... set of surf nodes at the intersection
           of this surface with any other surface

        Should be called for all possible intersections before
        self.findFreeBorder().
        """
        allSurfNodes = set(self.getNodeTris())
        for otherNodes in othersurface.triNodes.itervalues():
            nodesOnSurf = allSurfNodes.intersection(otherNodes)
            self.surfIntersectionNodes.update(nodesOnSurf)
            if len(nodesOnSurf)>0:
                self.cuttingFaces.add(frozenset(otherNodes))
        return


    def addCuttingSurfFreeBorder(self, othersurface):
        """updates self.freeSurfBorderNodes with common
        othersurface.freeSurfBorderNodes (free border nodes of the other
        surface at the intersection between the two surfaces)

        Should be called for all possible intersections after
        self.findFreeBorder() since it evaluates
        othersurface.freeSurfBorderNodes and before self.assignAllTetsConn().
        """
        newNodes = othersurface.freeSurfBorderNodes. \
            intersection(self.getNodeTris())
        self.freeSurfBorderNodes.update(newNodes)
        return


    #---- some check and repair functions
    def findSingleTriHoles(self):
        """Returns a list of node sets, each set defining a missing tri.

        Missing tris in that sense are three nodes connected only by border
        edges and not already connected by a tri element.

        See also TriangleSurface.findSingleTris(), it finds tri elements with
        only border edges."""

        allBorderEdges = self.getBorderEdges()
        allBorderNodes = self.getBorderNodes()
        allNodesToTri = self.getNodesToTri()

        newTris = set()
        for element, nodesOnSurf in self.tetSurfNodes.iteritems():
            nodesOnBorder = allBorderNodes.intersection(nodesOnSurf)
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
                if ((len(allBorderEdges.intersection(
                                self.triEdgesIter(faceNodes))) == 3)
                    and (faceNodes not in allNodesToTri)):
                    # all edges of this face are in self.borderEdges
                    # and there is no tri connecting those edges
                    newTris.add(faceNodes)

        # finished
        return newTris

    def repairSingleTriHoles(self, elemOffset=None, debugElset=None):
        """
        @param debugElset: if not None add diagnostic elsets to this dict
        {elset name: elset}.

        @returns: a tuple:
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

        if len(newTriList) and debugElset is not None:
            warn_elset_name = "_FILL-TRI"
            debugElset[warn_elset_name] = set(newTriList)
            msg("Added %d tris to this surface fixing a hole.\n"
                "... Added them to the elset %s for debugging."
                % (len(newTriList), warn_elset_name))

        return (elemOffset, newTriList)


    def findCorruptedBorderTris(self, debugElset=None):
        """Find triangles on the border that better seem to be deleted.

        If there is a tet with all four nodes on the surface and one of the
        tris lying on the tets faces is connected to three border nodes (at
        least two border edges), that tri is considered a tri that currupts
        the border and should be deleted.

        @param debugElset: if not None add diagnostic elsets to this dict
        {elset name: elset}.
        """

        allBorderNodes = self.getBorderNodes()
        allNodesToTri = self.getNodesToTri()

        lostTris = set()
        for element, nodesOnSurf in self.tetSurfNodes.iteritems():
            if len(nodesOnSurf) < 4:
                continue
            nodesOnBorder = allBorderNodes.intersection(nodesOnSurf)
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
                    triId = allNodesToTri[nodesOfThisTri]
                except KeyError:
                    continue
                lostTris.add(triId)

        # finished
        if len(lostTris) and debugElset is not None:
            warn_elset_name = "_BADBORDER"
            debugElset[warn_elset_name] = set(lostTris)
            msg("Removed %d tris from this surface because they were"
                " connected to three border nodes and a tet with all"
                " four nodes on the surface.\n"
                "... Added them to the elset %s for debugging."
                % (len(lostTris), warn_elset_name))

        return lostTris


    # some iterator function
    @staticmethod
    def tetFacesIter(nodesOfThisTet):
        """Returns an iterator of the faces of the given tet element
        nodesOfThisTet gives the node ids, it is converted to a tuple

        use as in

        >>> nodes = mdb.el_node[5731]  # list of node ids, e.g. [3, 43, 12, 36]
        >>> for faces in surf.tetFacesIter(nodes):

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

        >>> nodes = mdb.el_node[5731]  # list of node ids, e.g. [3, 43, 12, 36]
        >>> for edges in tetEdgesIter(nodes):

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
    ##    separate mesh at surface (connection based algorithm)
    ##

    # dict: { ordered tet node indices of a face : side index }
    tetOnTriSideDict = dict.fromkeys([
            (0,1,2), (1,2,0), (2,0,1),
            (1,0,3), (0,3,1), (3,1,0),
            (2,1,3), (1,3,2), (3,2,1),
            (0,2,3), (2,3,0), (3,0,2)], 1)
    tetOnTriSideDict.update(dict.fromkeys([
                (0,2,1), (2,1,0), (1,0,2),
                (1,3,0), (3,0,1), (0,1,3),
                (2,3,1), (3,1,2), (1,2,3),
                (0,3,2), (3,2,0), (2,0,3)], 0))

    def assignAdjacentTets(self, debugElset=None):
        """
        first step of splitting the mesh

         1. build dictionaries
          - self.tetToFace    {surface-Tet-elementID: set(node1, node2, node3)}
          - self.faceToTet    {set(node1, node2, node3): list of elements}

         2. find tets directly adjacent to surface (sharing more than two nodes
            with the surface)


        initializes:

         - self.assignedTets [set(), set()], set of tets on neg side, pos side
         - self.assignedTetFaces  [neg side list, pos side list], contains only
           'early assigned' tets and their faces not on this surface.
           'early assigned' tets are those tets that share three nodes with the
           surface, their assignment to either side at this point is 'earlier'
           than those of the other tets to follow
            - neg/pos side lists are lists cointaining tuples (tet elem, (free
              faces, border flag))
            - freeFaceTups is a list of (free face, border flag)-tuples
            - free face is a set of node numbers, only faces not on the surface
            - the border flag is True, if the edge connecting the free face
              and the corresponding surface tri is on the border
         - self.assignedTetFaces [list(), list()],
           each list has tuples (tet elem, free faces)
         - self.cohesiveEls { (negElement,posElement):
           [[nodeindices of matching nodes in negElement], [same in posEl.]] }

        @param debugElset: if not None add diagnostic elsets to this dict
        {elset name: elset}.

        @returns: a tuple of the following three values:

         - wrongNumOfTetsOnTri ... {tri Id: set of tets} If a surface tri has
           not exactly two tets attached to it face to face (three coincidal
           nodes) this tri is added to this dict with a set of all connected
           tets.

         - twoTetsOnSameSide ... list of (triElem, tetsOnThisTri)-tuples.
           tetsOnThisTri is the set of tets attached to this tri face to face

         - tetsOnBothSide ... a set of tets that could have been assigned to
           both sides of the surface (but in the end it's only assigned to one
           side).
        """

        # reorder all triangle connectivity lists in self.triNodes according
        # to the orientation (counter clock wise when looking from "above")
        allTriNormals = self.getTriNormal()

        # some properties of the tri surface
        allNodesToTri = self.getNodesToTri()
        allBorderEdges = self.getBorderEdges()

        # self.tetToFace: {surface-Tet-elementID: set(node1, node2, node3)}
        self.tetToFace = dict()
        # self.faceToTet: {set(node1, node2, node3): list of elements}
        self.faceToTet = defaultdict(list)
        #
        for element, thisTetNodes in self.tetNodes.iteritems():
            self.tetToFace[element] = set()
            for thisFace in self.tetFacesIter(thisTetNodes):
                self.tetToFace[element].add(thisFace)
                self.faceToTet[thisFace].append(element)

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
        for triElem, thisTriNodes in self.triNodes.iteritems():
            #### needed??? 1 line
            oneNode = thisTriNodes[0]  # any one node on the surface triangle

            thisSurfFace = frozenset(thisTriNodes)
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
                thisTetNodes = self.tetNodes[tetElem]
                triOnTetIndices = [thisTetNodes.index(node)
                                   for node in thisTriNodes]
                sideidx = self.tetOnTriSideDict[tuple(triOnTetIndices)]

                # face on surface is not a free face
                # have to compare to all faces (nodes) of this surface
                # not just 
                freeFaces = self.tetToFace[tetElem].difference(
                    allNodesToTri)

                # The edge connecting freeFaces to the surface: is it on the
                # surface border?
                freeFaceTups = list()
                for thisFreeFace in freeFaces:
                    thisSurfEdge = thisFreeFace.intersection(thisSurfFace)
                    borderFlag = thisSurfEdge in allBorderEdges
                    freeFaceTups.append((thisFreeFace, borderFlag))

                # Check if this tet has already been asigned to the other side.
                if tetElem in self.assignedTets[1-sideidx]:
                    # tet already on other side! Leave it there
                    tetsOnBothSide.add(tetElem)
                    continue

                # add new faces
                self.assignedTetFaces[sideidx].append((tetElem, freeFaceTups))
                # assign tet to the side we just determined
                self.assignedTets[sideidx].add(tetElem)
                elementPair[sideidx] = tetElem
                nodeIndexList[sideidx]=[self.tetNodes[tetElem].index(coh_node)
                                        for coh_node in self.triNodes[triElem]]


            # each face of the surface should have one tet on each side
            # if not, this tri is being ignored
            if ((elementPair[0] is None or elementPair[1] is None)):
                # The case that there are not exactly two tets assigned to this
                # tri element is dealt with earlier by checking the condition
                # len(tetsOnThisTri) != 2. So if there is no tet on either side
                # two tets must be on the same side!
                twoTetsOnSameSide.append((triElem, tetsOnThisTri))
                continue

            # update cohesiveEls-dict
            # cohesiveEls: { (negElement,posElement) ->
            #    [[nodeindices of matching nodes in negElement],
            #     [same in posEl.]] }
            elementPair = tuple(elementPair)
            if elementPair not in self.cohesiveEls:
                self.cohesiveEls[elementPair]=nodeIndexList

            if (elementPair[1],elementPair[0]) in self.cohesiveEls:
                raise Exception('double cohesive element')

        # Warnings
        if debugElset is not None:
            if wrongNumOfTetsOnTri:
                warn_elset_name_tri = "_NUMTET-TRI"
                warn_elset_name_tet = "_NUMTET-TET"
                elsTris = set(); debugElset[warn_elset_name_tri] = elsTris
                elsTets = set(); debugElset[warn_elset_name_tet] = elsTets
                for tri, tets in wrongNumOfTetsOnTri.iteritems():
                    elsTris.add(tri)
                    elsTets.update(tets)
                msg("There are %d triangles in the surface, that do not"
                    " have exactly two tets which are connected to all"
                    " three nodes. They will be ignored.\n"
                    "... The tris are collected in elset %s, the connected"
                    " tets in %s."
                    % (len(wrongNumOfTetsOnTri),
                       warn_elset_name_tri, warn_elset_name_tet))
            if twoTetsOnSameSide:
                warn_elset_name_tri = "_INTERSEC-TET-TRI"
                warn_elset_name_tet = "_INTERSECTING-TET"
                elsTris = set(); debugElset[warn_elset_name_tri] = elsTris
                elsTets = set(); debugElset[warn_elset_name_tet] = elsTets
                for tri, tets in twoTetsOnSameSide:
                    elsTris.add(tri)
                    elsTets.update(tets)
                msg("There are %d triangles in the surface with the two"
                    " tets connected to all three nodes are on the same"
                    " side of the tri. Those tets obviously intersect each"
                    " other. They will be ignored.\n"
                    "... The tris are collected in elset %s, the connected"
                    " tets in %s."
                    % (len(twoTetsOnSameSide),
                       warn_elset_name_tri, warn_elset_name_tet))
            if tetsOnBothSide:
                warn_elset_name = "_BOTHSIDETETS"
                debugElset[warn_elset_name] = set(tetsOnBothSide)
                msg("%d tets are found to be on both sides of this"
                    " surface. This happened during 'early assignment'"
                    " looking at tets directly adjacent to the surface.\n"
                    "... Added them to the elset %s for debugging."
                    % (len(tetsOnBothSide), warn_elset_name))

        return (wrongNumOfTetsOnTri, twoTetsOnSameSide, tetsOnBothSide)


    def findFreeBorder(self):
        """Investigate the border of the surface.
        Find free borders, where surf tet neighbours connect upper to
        lower side. Those nodes will not be split.

        Add also all nodes in self.moebiusStopEdges to freeSurfBorderNodes
        because those nodes shall not be split either and shall be treated the
        same.

        initializes:
         - self.freeSurfBorderNodes: set of surf nodes where surf ends within
           the model

        uses:
         - self.assignedTets with only early assigned tets
           ... as from self.assignAdjacentTets()

        @returns: a tuple of error cases:

         - borderTetMissing: list of tri ids.
           triangles on the border of the surface not having tets on both sides
           That error may also occur if two splitting surfaces partly coincide.

         - borderTetDouble: list of tuples (edge, side_idx, tets on this side).
           triangle on the border of the surface having more than one 'early
           assigned' tets (sharing three nodes with the surface) on either
           side. Those tets must be intersecting each other.
           This error may also occur when the tets on intersecting surfaces
           where heavily distorted when the mesh was split by a too large
           amount.

         - wrappedTets: list of tet elements.
           tets on the border of the surface which are connected to the surface
           with all their nodes.

         - danglingTets: list of tet elements.
           tets on the border of the surface that do not share a face (three
           nodes) with the surface. That is a contradiction because 'early
           assigned' tets are identified by sharing three nodes with the
           surface.
        """

        # some properties of the tri surface
        allBorderEdges = self.getBorderEdges()
        edgeToTri = self.getEdgeToTri()
        allNodeTris = self.getNodeTris()

        if not hasattr(self, "assignedTets"):
            # you should rather call assignAdjacentTets before to get
            # diagnostic data and create some output
            self.assignAdjacentTets()

        # dict {edges: list of connected tetElems} ( edge=set(node1,node2) )
        borderEdgeToTet = dict()
        for thisEdge in allBorderEdges:
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
        for thisEdge in allBorderEdges:

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
                    borderTetMissing.append(tuple(edgeToTri[thisEdge])[0])
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
                else:  # if len(assignedTetsOnThisEdge[side_idx])==1
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
                nextnode = set(self.tetNodes[nextelem]).difference(allNodeTris)
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

        # some properties of the tri surface
        allNodesOnSurface = set(self.getNodeTris())
        allBorderNodes = self.getBorderNodes()

        if not hasattr(self, "freeSurfBorderNodes"):
            # you should rather call findFreeBorder before to get diagnostic
            # data and create some output
            self.findFreeBorder()

        # notAssignedTets: Tets yet to be assigned to either side
        notAssignedTets = set()
        allNodesToSplit = allNodesOnSurface - self.freeSurfBorderNodes
        for element, nodesOfThisTet in self.tetNodes.iteritems():
            if len(allNodesToSplit.intersection(nodesOfThisTet))>0:
                notAssignedTets.add(element)

        # some are already assigned (those sharing three nodes with the surface)
        notAssignedTets -= self.assignedTets[0]
        notAssignedTets -= self.assignedTets[1]

        # tetsToFindNeighbours is a deque (list) with tuples for each tet yet
        # to be processed to find it's neighbours. Each tuple consists of:
        # (element number, list of faces to possible neighbours, side index)
        # initialize tetsToFindNeighbours with a tet connected with
        # at least three nodes to the surface (early assigned tet)
        tetsToFindNeighbours = deque()
        for sideidx in range(2):
            for lastTet, neighbourFaces in self.assignedTetFaces[sideidx]:
                tetsToFindNeighbours.append((lastTet, neighbourFaces, sideidx))

        # following loop: find neighbours of already assigned tets
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
        while len(tetsToFindNeighbours) > 0:

            lastTet, neighbourFaces, sideidx = tetsToFindNeighbours.popleft()
            for thisFace, onBorder in neighbourFaces:
                # onBorder means: this face is connected to the surface
                # only on the surfaces border
                # Or it has been reached via such a face. In that case
                # it must be connected to the initial face on the
                # border by at least one node. So we finds everything
                # connected to those bordernodes on the same side but
                # we don't proceed sideways along the border

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
                    continue  # no other tets on this face or error
                newTet = newTets[0]

                if newTet not in notAssignedTets:
                    continue

#               I think: check that later and first assign this tet to the right side
#               but don't add new face 
#                 # if on border don't proceed to new surface nodes, collect only
#                 # tets connected to the old nodes on the surface
#                 if onBorder:
#                     facesNodesOnSurf = thisFace.intersection(
#                         allNodesOnSurface)
#                     newTetNodesOnSurface = set(self.tetNodes[newTet]).\
#                         intersection(allNodesOnSurface)
#                     if len(newTetNodesOnSurface-facesNodesOnSurf)>0:
#                         continue

                # ok, found (valid) neighbour tet: newTet
                notAssignedTets.remove(newTet)
                if newTet in self.assignedTets[1-sideidx]:
                    # already on other side!
                    tetsOnBothSide.add(newTet)
                else:
                    self.assignedTets[sideidx].add(newTet)

#              this chunk seems to be too cautious in deciding to stop at the free surface
#              border (the border not at the outer border of the whole model)
#              It's been replaced by the one below, which seems to give better results
#                 # don't step further from this tet if it is
#                 # connected to freeSurfBorderNodes
#                 if set(self.tetNodes[newTet]).intersection(
#                     self.freeSurfBorderNodes):
#                     continue

                # don't step further from this tet if it is
                # connected to the surface only on freeSurfBorderNodes
                newTetNodesOnSurface = (set(self.tetNodes[newTet]).
                                        intersection(allNodesOnSurface))
                if len(newTetNodesOnSurface.difference(
                        self.freeSurfBorderNodes)) == 0:
                    continue

                # find free faces and recursivly find neighbouring tets
                # is one of newFaces connected to the border?
                newFaces = set(self.tetToFace[newTet])
                newFaces.remove(thisFace)

                newFaceTups = list()
                for thisFreeFace in newFaces:
                    facesNodesOnSurf = thisFreeFace.intersection(
                        allNodesOnSurface)
                    onBorder_new = (facesNodesOnSurf <= allBorderNodes)
                    newFaceTups.append((thisFreeFace, onBorder_new))

                tetsToFindNeighbours.append((newTet, newFaceTups, sideidx))


        # ignore a tet in notAssignedTets if its nodes that are connected to
        # the surface all are either in self.freeSurfBorderNodes or in
        # self.surfIntersectionNodes. In that case (and since the algorithm did
        # not find them, they are behind an intersecting surface and this
        # surface does not split the nodes on those tets.
        # note that this definition of allNodesToSplit is different from the
        # one earlier in this function
        allNodesToSplit = allNodesOnSurface - (
            self.freeSurfBorderNodes | self.surfIntersectionNodes)
        for tet in list(notAssignedTets):
            nodes = self.tetNodes[tet]
            if len(allNodesToSplit.intersection(nodes))==0:
                notAssignedTets.remove(tet)

        # fini
        return (moreThanTwoTets, notAssignedTets, tetsOnBothSide)

    def assignSomeTetsGeom(self, notAssignedTets):
        """assigne the specified tets to either side based on the position of
        the tets centroid.

        Geometry based algorithm:
        Assign not yet assigned tets to either side of the surface by comparing
        the vector from the tets nodes on the surface to its centroid with the
        average node normal of the connected surface nodes.
        """

        # some properties of the tri surface
        allNodeTris = self.getNodeTris()
        allNodeNormal = self.getNodeNormal()

        # separate tets at surface / by average of node normals and centroid
        for thisTet in notAssignedTets:
            normalSum=[0,0,0]
            pointOnSurf=None  # average of at most two nodes on the surface
            centroid=[0.0, 0.0, 0.0]
            for i, node in enumerate(self.tetNodes[thisTet][:4]):
                coords=self.nodeCoords[node]
                # calculation of centroid
                vector_modif_add(centroid, vector_scale(coords, 0.25))

                # surface normal and point on surface
                if node in allNodeTris: # nodes on surface
                    normalSum=vector_plus(normalSum,allNodeNormal[node])
                    if pointOnSurf:
                        pointOnSurf=vector_scale(
                            vector_plus(pointOnSurf, coords), 0.5)
                    else:
                        pointOnSurf=coords

            vectorNodeCent=vector(pointOnSurf, centroid)

            side=dot(vectorNodeCent,normalSum)
            if side>0.:
                self.assignedTets[1].add(thisTet)
            else:
                self.assignedTets[0].add(thisTet)

        # fini
        return


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

        # some properties of the tri surface
        allNodeTris = self.getNodeTris()
        allBorderNodes = self.getBorderNodes()
        allNodeNormal = self.getNodeNormal()

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
                           if node in allNodeTris
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
                if node in allNodeTris: # nodes on surface
                    normalSum=vector_plus(normalSum,allNodeNormal[node])
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
    def updateSplitNodesTets(self, splitNodes, blackList=set()):
        """
        updates:
         - splitNodes a dict : {node id: [splitting 1, splitting 2, ...]}
           each splitting is a three-tuple (neg. side attached tets, pos. side
           attached tets, averaged splitting normal)

        uses:
         - self.freeSurfBorderNodes, self.assignedTets, self.getNodeNormals()
        """

        allNodeNormal = self.getNodeNormal()
        for nodeA in self.getNodeTris():
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
                                      allNodeNormal[nodeA]))



    ## ------------------------------------------------------------------------
    ##   update splitNodes
    ##   
    ##   insert a splitting-entry for each tri element in triNodes
    ## ------------------------------------------------------------------------
    def updateSplitNodesOtherTris(self, splitNodes, blackList, triNodes):
        """
        updates:

         - splitNodes a dict : {node id: [splitting 1, splitting 2, ...]}.
           each splitting is a three-tuple (neg. side attached tets, pos. side
           attached tets, averaged splitting normal)

        not working correctly, do it later....
        (if there is a need to also correctly split the surfaces)
        """

        allNodeTris = self.getNodeTris()
        coincideTris = list()
        for triElem, elemNodes in triNodes.iteritems():
            if triElem in self.triNodes:
                continue

            # check for surface nodes
            replaceNodesIdx=[(node, idx)
                             for idx, node in enumerate(elemNodes)
                             if node in allNodesTris]
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
            if sideidx is None:
                for tetElem in self.assignedTets[0]:
                    if freeNode in self.tetNodes[tetElem]:
                        sideidx = 0 # lower/neg side
                        break
            if sideidx is None:
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


## ---------------------------------------------------------------------------

class SurfaceInVolFromMeshTris(TriangleSurfaceInVolume):
    """A TriangleSurfaceInVolume created from some tri elements of an abaqus
    model
    """

    def __init__(self, model, elset=None, logfile=sys.stdout):
        """Arguments to the constructor:

        Only those nodes connected to the surface and adjacent tets are copied
        to the self.nodeCoords attribute. But this is not a deep copy: The
        coordinate lists referenced by the dict self.nodeCoords are the same as
        those referenced by the model.nodeCoords dict.

        The node connectivity from model.elNodes however is deep copied: The
        node lists referenced by the dict self.triNodes are others then those
        referenced by model.elNodes.

        Elements in elset that are not found in model.elNodes are ignored, a
        warning is posted.

        @param model: is an abq_model_02.Model object

        @param elset: the surface shall consist of the elements in this elset
              might be anything that model.getUnionSet() accepts as input:
              an elset name, a set of element numbers, a list of elset names
              and element numbers
              elset might also be None, then all tri elements of model are used.

        @param logfile: deprecated argument, only for compatibility reasons.
        No effect.
        """
        TriangleSurface.__init__(self)
        checkModelModuleVersion(model, "%s.%s.__init__"
                                % (__name__, self.__class__.__name__))

        # initialize all tris: self.triNodes
        alltris = set()
        for typ in model.elShapeToTypes['TRI'].intersection(model.typeEl):
            alltris.update(model.typeEl[typ])

        if elset is None:
            elset = alltris.intersection(model.elNodes)
        else:
            elset = model.getUnionSet("elset", elset)
            elset.intersection_update(alltris)

        nodestocopy = set()
        nb_trisnotfound = 0
        for elem in elset:
            try:
                nodes = model.elNodes[elem][:3]
                self.triNodes[elem] = nodes
                nodestocopy.update(nodes)
            except KeyError:
                nb_trisnotfound += 1
                pass
        if nb_trisnotfound>0:
            msg("WARNING from %s.__init__():"
                " Could not find element connectivity for %d"
                " of the %d tri elements in the specified"
                " elset." % (self.__class__.__name__,
                             nb_trisnotfound, len(elset)))
        if len(self.triNodes)==0:
            msg('WARNING: %s initialized with empty element set.'
                % self.__class__.__name__)

        # create the connectivity of the surface to the volume
        # self.tetNodes {tetID:corner nodes} only tets with nodes on the surf
        # self.tetSurfNodes {tet: [nodes on surface]}
        #                              where len([nodes on surf])>0
        # self.nodeTets  {node:[tetIDs,]} for tets in self.tetSurfNodes
        self.tetNodes = dict()
        self.tetSurfNodes = dict()
        self.nodeTets = defaultdict(list)

        allSurfNodes = set(self.getNodeTris())
        alltets = set()
        for typ in model.elShapeToTypes['TET'].intersection(model.typeEl):
            alltets.update(model.typeEl[typ])

        # loop over all tet elems in the mesh
        tetWithLessThenFourNodes = list()
        for element in alltets:
            nodes = model.elNodes[element][0:4]
            if len(nodes)<4:
                tetWithLessThenFourNodes.append(element)
                continue
            nodesOnSurf = allSurfNodes.intersection(nodes)
            if len(nodesOnSurf)==0:
                continue  # this tet is not connected to the surface
            # self.tetSurfNodes
            self.tetSurfNodes[element] = nodesOnSurf
            # self.NodeTets
            for node in nodesOnSurf:
                self.nodeTets[node].append(element)
            # self.tetNodes
            self.tetNodes[element] = nodes
            nodestocopy.update(nodes)

        if len(tetWithLessThenFourNodes):
            msg("WARNING, SEVERE INCONSISTENCY:\n"
                "%d Elements were expected to be tets but had less than four"
                " nodes:")
            for element in tetWithLessThenFourNodes:
                msg(" . element %d, type %s, nodes: %s"
                    % (element, model.elType[element], model.elNodes[element]))

        # copy all connected nodes (no deep copy)
        for node in nodestocopy:
            try:
                coords = model.nodeCoords[node]
                self.nodeCoords[node] = coords
            except KeyError:
                pass  

        # self.cuttingFaces and self.surfIntersectionNodes
        # will be updated (if nessecary) by addCuttingSurf()
        self.cuttingFaces = set()
        self.surfIntersectionNodes = set()

        return

## ---------------------------------------------------------------------------

class SurfaceInVolFromTriangleSurface(SurfaceInVolFromMeshTris):
    """A TriangleSurfaceInVolume created from an already existing
    TriangleSurface-object. It gets the connectivity information for the
    adjacent tet elements from an abq_model_02.Model object.
    """

    def __init__(self, surf, model):
        SurfaceInVolFromMeshTris.__init__(
            self, model, elset=surf.triNodes.keys())

        # the following permanent attributes created by the constructor are
        # simply copied from surf

        # copy unconditionally
        for attr in ['triNormal', 'edgeToTri', 'borderEdges','borderNodes',
                     'nodesToTri', 'nodeNormal', 'nodeTris']:
            setattr(self, attr, copy.deepcopy(getattr(surf, attr)))

        # copy if attribute exists
        for attr in ['nodesToDoubles', 'moebiusStopEdges', 'moebius_flag']:
            if hasattr(surf, attr):
                setattr(self, attr, copy.deepcopy(getattr(surf, attr)))
