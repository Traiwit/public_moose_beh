"""surface_01.py: DEPRECATED, use L{bae.surface_03}

surface class
"""

_version_ = "1.05"

_version_history_ = """\
Versions:
=========

1.01 new   GP: new from surfToContact 2.06, ignore edges with more than two
           triangles in calcNormals
1.02 GP added: some further checks and remove tri from surface functions
     insert tri, getBorderNodes, getNodesToTri
1.03 GP added: in class TriangleSurface added getNodeNormals(),
     self.getNodeTris() and self.getNodes()
     removed removeNoNormal(), added checkNoNormal() instead
1.04 GP changed: checkNoNormal() returns a set not a dict
     calcNormals starts repeatedly, if the surface consists of not connected
     parts, so in the end each triangle has a normal
     added: get_abq_model()
1.05 GP added: set self.moebiusStopEdges created by self.calcNormals()
     get_simple_parts_with_normals(), added check to constructor
"""

## ---------------------------------------------------------------------------
_bugs_ = """\
known bugs/problems:
====================
- needs python 2.4, if in doubt use abaqus python

ideas / todo:
=============

"""

from collections import deque  # deque = double ended queue
from sys import stdout
from math import cos, pi
import copy

from bae.future_01 import *
from bae.vecmath_01 import *
from bae.abq_model_01 import QuietLog, Model, NodesDict
# IMPORTANT: look at the end of the file, there is yet another import

class _SurfaceConnectionError(Exception):
    """raised when there is an error in the connection of the elements
    forming the surface"""
    pass

class _NotImplementedError(Exception):
    """raised when the requested action is not yet implemented, for example
    for a particular element type"""
    pass

class TriangleSurface(object):
    """surface made of plane triangles

    - Keeps a copy of the tri-element list (argument triList) as a set.
    - Creates the dict self.triNodes {element:[node, node,...]} by searching
      all elements of triList in the dict triNodes {element:[node, node,...]},
      which might contain much more elements.
    - Keeps a reference of the {node: [coord1,..]} dictionary
      (argument nodeCoords)

    arguments of the constructor:
    -----------------------------

    surfaceName is just the name used for error messages and warnings
    triList is the list of tri element ids that form this surface
    ...
    logfile (default: stdout) write warnings to that stream
        if logfile == None, no diagnostic output
    initchecks (default: True) check some of the arguments when creating a new
        TriangleSurface object.

    methods meant to be called from outside: (maybe not complete)
    ----------------------------------------

    self.calcNormals() -> yields attribute self.triNormals
    
    self.checkSimplyConnected()

    self.findSingleTris() -> list of tri elements without neighbours

    self.getBorderNodes() -> set of all nodes anywhere on the border.

    self.getNodesToTri() -> {face (set of three nodes): tri element number}

    self.getNodeNormals() -> {node number -> averaged normal of the surface}

    self.getNodeTris() -> dict {node id: [tri ids]}

    self.getNodes() -> set of nodes belonging to the surface
       (this may be superfluous, you can use self.getNodeTris()
        avoid using it, it may be removed, and takes extra time and memory)
    """

    def __init__(self, surfaceName, triList, triNodes, nodeCoords,
                 logfile=stdout, initchecks = True):
        "see class description"
        # if logfile == None, no diagnostic output
        if logfile == None:
            self.logfile = QuietLog()
        else:
            self.logfile = logfile

        # the following are always there
        self.surfaceName = surfaceName
        self.triNodes = dict([(elem, list(triNodes[elem])) for elem in triList])
        self.nodeCoords = nodeCoords

        # the following are None as long as they are not calculated yet
        self.borderNodes = None  # don't query this variable, use getBorderNodes()
        self.nodesToTri = None  # don't query this variable, use getNodesToTri()
        self.nodeNormals = None # don't query this variable, use getNodeNormals()
        self.nodeTris=None      # don't query this variable, use getNodeTris()
        self.nodes=None         # don't query this variable, use getNodes()

        if initchecks:
            for node in self.getNodes():
                try:
                    self.nodeCoords[node]
                except KeyError:
                    raise ValueError("coords for node %d not specified when"
                                     " creating new TriangleSurface object %s."
                                     " Check arguments to constructor!"
                                     % (node, surfaceName))


    @classmethod
    def fromAbqModelSurface(self, abqModel, surfName,
                                  logfile=stdout, initchecks=True):
        """create a TriangleSurface from the surface named surfName in the
        ABAQUS model abqModel.

        The surface must consist of faces of tet elements (C3D4 or C3D10M).
        Only the corner nodes of the tets are considered.

        The elsets abqModel+"_S1", abqModel+"_S2", ... "+_S4" define the
        surface.
        
        The node coords of this surface are a reference (no copy) to those of
        the abaqus model. If the abaqus model moves, so does this surface.
        """

        triElementFaces={'S1':[0,1,2],'S2':[0,1,3],'S3':[1,2,3],'S4':[0,2,3]}

        triNodes = dict()
        elnum = 1
        for faceType, nodeIds in triElementFaces.iteritems():
            for el in abqModel.elset["%s_%s"%(surfName, faceType)]:
                if abqModel.el_type[el] not in ('C3D4', 'C3D10M'):
                    raise _NotImplementedError(
                        "createFromAbqModelSurface() can not create a surface"
                        " from elements of type %s." % abqModel.el_type[el])
                nodes = abqModel.el_node[el]
                triNodes[elnum] = [nodes[i] for i in nodeIds]
                elnum += 1
                
        triList = triNodes.keys()
        triList.sort()
        surf = TriangleSurface(surfName,
                               triList=triList, triNodes=triNodes,
                               nodeCoords=abqModel.node_coord,
                               logfile=logfile, initchecks=initchecks)
        return surf



    def get_abq_model(self, el_offset=0, node_offset=0):
        """create an abq_model from this surface
        
        node numbers start with node_offset+1, element numbers with el_offset+1
        node coords are copied, not referenced, if surfaces node coords change,
        those of the abq_model don't.

        as a not finished yet!, not testet! no docs yet"""
        m_out = Model()
        all_nodes = dict()
        for el, nodes in self.triNodes.iteritems():
            el_offset += 1
            new_nodes = list()
            for node in nodes[:3]:
                try:
                    new_nodes.append(all_nodes[node])
                except KeyError:
                    node_offset += 1
                    all_nodes[node] = node_offset
                    new_nodes.append(node_offset)
            m_out.update_elem(el_offset, 'S3', new_nodes)
        m_out.node_coord = NodesDict()
        for node_old, node_new in all_nodes.iteritems():
            m_out.node_coord[node_new] = list(self.nodeCoords[node_old])

        return m_out


    def get_simple_parts_with_normals(
        self, mintris=1, floodangle=None, skipMoebiusCheck=False):
        """Creates a list of TriangleSurface objects of parts of self.

        Each new surface is simply connected (no T-junction) and all
        triangles are connected. The new surfaces are in a state as if
        calcNormals() and calcEdges() had been called on them. Essentially that
        means: The normals of the triangles are already calculated
        (self.triNormals has been created) and edgeToTri and borderEdges are
        there.

        arguments:
        mintris ... discard all surfaces with less then mintris triangles
        floodangle ... stop at edges with an angle greater than floodangle
                       (in degrees)
        skipMoebiusCheck ... no check if the surface has two sides

        return tuple:
        surflist .. list of TriangleSurface objects
        skipped_parts_cnt .. number of parts, that were skipped because <mintris
        tris_in_skipped_parts .. list of corresponding tri elements
        """

        surflist = list()
        skipped_parts_cnt = 0
        tris_in_skipped_parts = list()

        mintris = max(mintris, 1)

        if not hasattr(self, "edgeToTri"):
            self.calcEdges()

        self.triNormals = dict()
        if len(self.triNodes)==0:
            return surflist, skipped_parts_cnt, tris_in_skipped_parts

        # local flag and set used to sum up
        all_moebius_flag_cnt = 0
        all_moebiusStopEdges = set()

        # set of tris for which a normal has not been calculated yet
        noNormalTris = set(self.triNodes)
        
        while len(noNormalTris)>0:

            # get any one element as a starting point
            startTriId = iter(noNormalTris).next()

            # we want the flag and the set for each part surface alone
            self.moebius_flag = False
            self.moebiusStopEdges = set()

            # calculate normals starting at that point
            trisWithNormalList = self.calcNormalsOnPart(
                startTriId, floodangle, skipMoebiusCheck)
            noNormalTris.difference_update(trisWithNormalList)

            # check mintris
            if len(trisWithNormalList) < mintris:
                skipped_parts_cnt += 1
                tris_in_skipped_parts.extend(trisWithNormalList)
                continue

            # sum up data for self
            if self.moebius_flag:
                all_moebius_flag_cnt += 1
            all_moebiusStopEdges.update(self.moebiusStopEdges)

            # create new surface
            newsurf = TriangleSurface(
                self.surfaceName, trisWithNormalList,
                self.triNodes, self.nodeCoords, self.logfile)

            # create other features of new part surface
            newsurf.calcEdges()
            newsurf.triNormals = dict([
                    (el, self.triNormals[el])
                    for el in trisWithNormalList ])
            newsurf.moebius_flag     = self.moebius_flag
            newsurf.moebiusStopEdges = self.moebiusStopEdges

            # add to list
            surflist.append(newsurf)

        self.moebius_flag     = all_moebius_flag_cnt>0
        self.moebiusStopEdges = all_moebiusStopEdges

        if self.moebius_flag:
            self.logfile.write(
                "WARNING: %d parts of the surface %s are twisted like the"
                " Moebius strip, they have only one side.\n"
                % (all_moebius_flag_cnt, self.surfaceName))

        return surflist, skipped_parts_cnt, tris_in_skipped_parts



    def calcEdges(self):
        """creates object-attribute
        edgeToTri: a dict {set(node1,node1):set of tri-element-ids}
        borderEdges: list of edges ( =set(node1,node2) ) with only one tri
           connected to
        """
        self.edgeToTri = defaultdict(set)
        for element, nodesOfThisTri in self.triNodes.iteritems():
            for thisEdge in self.triEdgesIter(nodesOfThisTri):
                self.edgeToTri[thisEdge].add(element)
        self.borderEdges = set([edge
                           for edge,tris in self.edgeToTri.iteritems()
                           if len(tris)==1])


    def getBorderNodes(self):
        """return the set of all nodes anywhere on the border"""
        if self.borderNodes != None:
            return self.borderNodes
        if not hasattr(self, "borderEdges"):
            self.calcEdges()
        
        self.borderNodes = set()
        for edge in self.borderEdges:
            self.borderNodes.update(edge)
        return self.borderNodes


    def getNodesToTri(self):
        """return a dict {frozenset(three node numbers) : tri id}
        if it's not there yet, create it

        if there are double tris they are placed in self.nodesToDoubles
        {frozenset(node numbers) : list of tri ids}
        after calling this method len(self.nodesToDoubles) may be checked
        """
        if self.nodesToTri != None:
            return self.nodesToTri

        # it's not there already, create it (new)
        self.nodesToTri = dict()
        self.nodesToDoubles = defaultdict(list)
        for element, nodesOfThisTri in self.triNodes.iteritems():
            face = frozenset(nodesOfThisTri)
            if face in self.nodesToTri:
                self.nodesToDoubles[face].append(element)
            else:
                self.nodesToTri[face] = element
        return self.nodesToTri


    def getNodeTris(self):
        """creates or returns the previously calculated
        self.nodeTris: {node id: [tri ids]}
        """
        if self.nodeTris==None:
            # create self.nodeTris: {node id: [tri ids]}
            self.nodeTris = defaultdict(list)
            for element in self.triNodes:
                nodesOfThisTri = self.triNodes[element]
                for node in nodesOfThisTri:
                    self.nodeTris[node].append(element)
        # anyway return the dict
        return self.nodeTris


    def getNodes(self):
        """returns a set of nodes belonging to the surface
        """
        if self.nodes==None:
            self.nodes=set(self.getNodeTris().iterkeys())
        return self.nodes


    def checkSimplyConnected(self):
        """check if the surface triangles are simply connected

        all edges must have exactly one or two triangles connected to it
        """
        if not hasattr(self, "edgeToTri"):
            self.calcEdges()
        return all([(len(tris) in (1,2)) for tris in self.edgeToTri.values()])
            

    def findSingleTris(self):
        """Returns a list of tri elements without neighbours"""
        if not hasattr(self, "borderEdges"):
            self.calcEdges()
        return [triId
                for triId, nodes in self.triNodes.iteritems()
                if len(self.borderEdges.intersection(
                    self.triEdgesIter(nodes))) == 3]

    def calcAngleToNeighbours(self, triNodes, default=None):
        """calculates the cosines of the angles between the normals of the tri
        given by triNodes and of all of its neighbours.

        Their does not have to be a tri on the surface connecting triNodes,
        this is meant to check the angle before actually inserting a new tri.

        The sign of the angle is arbitrary, since node ordering is not taken
        into account.

        The result is always a list of three members. Where there is no
        neighbour, the value of the argument "default" is returned in that
        list.
        """
        if not(hasattr(self, "edgeToTri")):
            self.calcEdges()

        n1 = norm(self.getTriNormalUsc(tuple(triNodes)))

        angles = [default, default, default]
        for edgeIdx, thisEdge in enumerate(self.triEdgesIter(triNodes)):
            nodes = set()
            for thisTri in self.edgeToTri[thisEdge]:
                nodes.update(self.triNodes[thisTri])
            nodes.difference_update(triNodes)
            if len(nodes)!=1: continue
            nodes.update(thisEdge)
            n2 = norm(self.getTriNormalUsc(tuple(nodes)))
            angles[edgeIdx] = dot(n1,n2)

        return angles



    def getTriNormalUsc(self, nodes):
        """Not meant to be called from "outside"

        Computes the unscaled ("Usc") normal of a triangle.
        Nodes come in counter-clockwise succession if you look from "above"
        (positive side), against the direction of the normal.
        """
        vectorAB=vector(self.nodeCoords[nodes[0]],self.nodeCoords[nodes[1]])
        vectorAC=vector(self.nodeCoords[nodes[0]],self.nodeCoords[nodes[2]])
        return cross(vectorAB,vectorAC)

    def getTriNormalSc(self, nodes):
        """Not meant to be called from "outside"

        Computes the scaled ("Sc") normal of a triangle.
        Nodes come in counter-clockwise succession if you look from "above"
        (positive side), against the direction of the normal.
        """
        vectorAB=vector(self.nodeCoords[nodes[0]],self.nodeCoords[nodes[1]])
        vectorAC=vector(self.nodeCoords[nodes[0]],self.nodeCoords[nodes[2]])
        return norm(cross(vectorAB,vectorAC))


    def calcNormalsOnPart(self, startTriId,
                          floodangle=None, skipMoebiusCheck=False):
        """Not meant to be called from "outside"

        Make sure self.edgeToTri exists, otherwise call self.calcEdges()
        before. Make sure len(self.triNodes)>0 and the dict self.triNormals
        exists.

        creates object-attribute triNormals:
        a dict {element:[comp1, comp2, comp3]} with [comp...] being the vector
        components of the *unscaled* normal of that triangle
        
        This function calculates the elements normal vectors of the surface
        triangles in such a way that all normals point to the same side of the
        surface. It starts with the element given by *startTriId* and proceeds
        with all neighbouring tris as long as they are simply connected.

        Additionally the orientation of the tri element is saved by means of
        reordering its nodes in the self.triNodes-dictionary in
        counter-clockwise succession if you look from "above" (positive side).
        (The calculated normal points "up" to the positive side.)

        At an edge with more than two triangles attached, the algorithm stops.
        At triangles that have their normal calculated already
        (newTri in self.triNormals) the algorithm stops as well.
        If the angle in degree between two neighbouring triangles is greater
        than *floodangle*, the algorithm stops at that edge. floodangle==None
        means: don't check angle.

        Simply connected means there are max two triangles connected to each
        edge.

        If skipMoebiusCheck==True, no check is performed, if the surface
        has two sides and the normals are determinable uniquely. Otherwise if
        contradicting results for the same normal are found, self.moebius_flag
        is set to True and the corresponding edge is added to
        self.moebiusStopEdges. Additionally a warning is logged in the log
        file.

        The function returns a list of all triangles it calculated the normals
        for.
        """

        if floodangle != None:
            floodanglecos = cos(float(floodangle)*pi/180.)

        trisWithNormalList = list()
        startTriNodes = self.triNodes[startTriId]
        startNormal = self.getTriNormalSc(startTriNodes)
        self.triNormals[startTriId] = startNormal
        trisWithNormalList.append(startTriId)

        # trisToFindNeighbours is a list of (TriId, TriNodes, Normal)-tuples
        # of those tris for which the neighbours have to be examined yet. The
        # normals of the tris in this list have already been calculated.
        # the TriNodes-item contains the node numbers in counter-clockwise
        # succession if you look from "above" (positive side)
        trisToFindNeighbours = deque()  # a double-ended queue
        trisToFindNeighbours.append((startTriId, startTriNodes, startNormal))

        while len(trisToFindNeighbours) > 0:
            # get the first (oldest) tri in the list
            startTriId, startTriNodes, startNormal = (
                trisToFindNeighbours.popleft())

            for thisEdge, edgeNodes in self.triEdgesOrderedNodesIter(
                startTriNodes):

                trisOnThisEdge = self.edgeToTri[thisEdge]
                if len(trisOnThisEdge)!=2:
                    # more than two tris on this edge: ignore this edge
                    # less than two: edge on border, no other tri connected
                    continue

                # find the other tri connected to this edge (not startTriId)
                for newTri in trisOnThisEdge:
                    if newTri != startTriId: break


                # if no check and normal already calculated, skip the rest
                newTriHasNormalAlready = (newTri in self.triNormals)
                if skipMoebiusCheck and newTriHasNormalAlready:
                    continue

                # find third node not on this edge
                newTriNodes = list(self.triNodes[newTri])
                for i in thisEdge:
                    newTriNodes.remove(i)
                # append nodes on the edge in correct order
                newTriNodes.append(edgeNodes[1])
                newTriNodes.append(edgeNodes[0])
                # calculate new normal
                newNormal = self.getTriNormalSc(newTriNodes)

                # check max angle between neighbouring normals
                if (floodangle != None and
                    dot(startNormal, newNormal)<floodanglecos):
                    continue

                # if normal already calculated, compare and goto next
                if newTriHasNormalAlready:
                    if dot(newNormal, self.triNormals[newTri]) < 0.0:
                        self.moebius_flag = True
                        self.moebiusStopEdges.add(thisEdge)
                        self.logfile.write(
                            "WARNING: The normals flip between element %d and"
                            " %d.\n" % (startTriId, newTri))

                    # in any case keep the already calculated normal and
                    # proceed with the next tri
                    continue

                # update self.triNormals with new normal
                self.triNormals[newTri] = newNormal
                trisWithNormalList.append(newTri)
                # update new node order to self.triNodes
                self.triNodes[newTri] = newTriNodes
                # append this neighbour to trisToFindNeighbours-list
                trisToFindNeighbours.append((newTri, newTriNodes, newNormal))

        # fini
        return trisWithNormalList


    def calcNormals(self, skipMoebiusCheck=False):
        """creates object-attribute triNormals:
        a dict {element:[comp1, comp2, comp3]} with [comp...] being the vector
        components of the *unscaled* normal of that triangle
        
        This function calculates for all tri-elements of the
        surface-elset the elements normal vectors of the surface,
        so that all normals point to the same side of the surface.

        Additionally the orientation of the tri element is saved by means of
        reordering its nodes in the self.triNodes-dictionary in
        counter-clockwise succession if you look from "above" (positive side).
        (The calculated normal points "up" to the positive side.)

        At an edge with more than two triangles attached, the algorithm stops,
        hopefully affected triangle are simply connected to the surface from
        the other side. To check for lost triangles call checkNoNormal()

        Simply connected means there are max two triangles connected to each
        edge.

        if skipMoebiusCheck==True, no check is performed, if the surface
        has two sides and the normals are determinable uniquely. Otherwise a
        warning is logged "surface is like the Moebius strip".
        """

        if not hasattr(self, "edgeToTri"):
            self.calcEdges()

        self.triNormals = dict()
        if len(self.triNodes)==0: return

        self.moebius_flag = False
        self.moebiusStopEdges = set()

        # set of tris for which a normal has not been calculated yet
        noNormalTris = set(self.triNodes)
        
        while len(noNormalTris)>0:

            # get any one element as a starting point
            startTriId = iter(noNormalTris).next()

            # calculate normals starting at that point
            trisWithNormalList = self.calcNormalsOnPart(
                startTriId, floodangle=None,
                skipMoebiusCheck=skipMoebiusCheck)
            noNormalTris.difference_update(trisWithNormalList)

        if self.moebius_flag:
            self.logfile.write(
                "WARNING: The surface %s is twisted like the Moebius"
                "strip.\n" % self.surfaceName)
        return


    def getNodeNormals(self):
        """returns dict
        {node number -> averaged normal of the surface}
        if not yet created do so, store in self.nodeNormals for second use
        """

        if (self.nodeNormals==None) or not(hasattr(self, "triNormals")):

            # if it's not there already, create it (new)
            if not hasattr(self, "triNormals"):
                self.calcNormals()

            # self.nodeNormals
            self.nodeNormals = dict()
            for node in self.getNodes():
                nodeNormal = [0,0,0]
                for element in self.getNodeTris()[node]:
                    nodeNormal = vector_plus(nodeNormal,self.triNormals[element])
                self.nodeNormals[node] = norm(nodeNormal)

        # return the dict anyway
        return self.nodeNormals



    def checkNoNormal(self):
        """search for tris from surface that have not got a normal
        call this only after calcNormals()
        returns a dict of removed triangles (as they are in self.triNodes)
        ... or None if normals have not been calculated yet.
        """
        if not hasattr(self, "triNormals"):
            return None

        obscureTris = set(self.triNodes).difference(self.triNormals)
        return obscureTris



    def removeTri(self, triId):
        """remove a tri from the surface

        update self.triNodes, self.edgeToTri, self.borderEdges
        remove self.borderNodes, self.nodesToTri
        """
        nodes = self.triNodes[triId]
        del self.triNodes[triId]
        if hasattr(self, "edgeToTri"):
            for thisEdge in self.triEdgesIter(nodes):
                self.edgeToTri[thisEdge].remove(triId)
                if len(self.edgeToTri[thisEdge]) == 0:
                    del self.edgeToTri[thisEdge]
                    self.borderEdges.remove(thisEdge)
                elif len(self.edgeToTri[thisEdge]) == 1:
                    self.borderEdges.add(thisEdge)

        # update self.nodeTris, self.nodes
        if self.nodeTris!=None:
            for node in nodes:  # all nodes of the tri to remove
                self.nodeTris[node].remove(triId)
                if len(self.nodeTris[node])==0:
                    del self.nodeTris[node]
                    if self.nodes!=None:
                        self.nodes.remove(node)
        else:
           self.nodes=None

        # delete obscured members
        # they have to be recalculated from scratch (no better idea)
        self.borderNodes = None
        self.nodesToTri = None
        self.nodeNormals = None


    def insertTri(self, newTri, nodes):
        """insert a new tri
        swap node ordering in order to have its normal consistent with its
        neighbours

        update self.triNodes, self.edgeToTri, self.borderEdges
        remove self.borderNodes, self.nodesToTri
        """
        
        if hasattr(self, "triNormals"):
            if not hasattr(self, "edgeToTri"):
                self.calcEdges()

        neighbour = None  # look for neighbour to determine normal

        # update edgeToTri and borderEdges and find a neighbour
        # ... only if edgeToTri already defined
        if hasattr(self, "edgeToTri"):
            for edge in self.triEdgesIter(nodes):
                try:
                    newneighbour = self.edgeToTri[edge]
                    if len(newneighbour)>1:
                        raise _SurfaceConnectionError(
                            "There are already more than one elements"
                            " connected to the edge %s when trying to insert"
                            " new tri element %d on nodes %s."
                            % (str(edge), newTri, str(nodes))) 

                    # if we don't have a neighbour yet, but need one...
                    if (neighbour==None and hasattr(self, "triNormals")):
                        newneighbour = tuple(newneighbour)[0]
                        if newneighbour in self.triNormals:
                            neighbour = tuple(newneighbour)[0]

                    # a border edge getting a new tri attached is no longer
                    # on the border
                    try:
                        self.borderEdges.remove(edge)
                    except ValueError:
                        # The new edge might also not have been there at all
                        # so no error, just pass.
                        # If it was there before it must have been a border
                        # edge but we checked that before by counting edgeToTri
                        pass

                except KeyError:
                    # this edges is not yet in self.edgeToTri, its a new edge
                    # a new edge is an edge on the border
                    self.borderEdges.add(edge)
                    
                # anyway update self.edgeToTri
                self.edgeToTri[edge].add(newTri)
                
        # must have at least one neighbour
        # at least when triNormals is defined
        if hasattr(self, "triNormals"):
            if (not neighbour):
                raise _SurfaceConnectionError(
                    "New tri element %d (nodes: %s) has no neighbour."
                    % (newTri, str(nodes)))
            # sort the nodes in an order such that its normal fits to
            # the neighbour
            newNodes = [0,0,0]
            NodesUnsorted = set(nodes)
            for idx, node in enumerate(self.triNodes[neighbour]):
                if node in NodesUnsorted:
                    newNodes[-idx] = node
                    NodesUnsorted.remove(node)
                else:
                    thirdIdx = -idx
            newNodes[thirdIdx] = NodesUnsorted.pop()

            # calculate normal
            self.triNormals[newTri] = self.getTriNormalUsc(newNodes)
        else:
            # if normals not computed yet, just take the nodes unsorted
            newNodes = list(nodes)

        # update triNodes
        self.triNodes[newTri] = newNodes

        # update nodeTris
        if self.nodeTris!=None:
            for node in newNodes:
                self.nodeTris[node].append(newTri)

        # delete obscured members
        self.borderNodes = None
        self.nodesToTri = None
        self.nodeNormals = None
        self.nodes=None


    def replaceNode(self, triId, oldNode, newNode):
        """Change surface topology: Reconnect triangle to other node.

        Reconnect the corner of the triangle triId which was originally
        connected to oldNode to newNode. That is used when the surface is split
        by surfToContact.

        Update all surface data accordingly. Don't remove the oldNode and don't
        replace the node in other triangles of the surface.
        """
        try:
            idx = self.triNodes[triId].index(oldNode)
        except ValueError:
            raise ValueError("Node %d not in tri element %d in surface %s"
                             " when trying to replace it."
                             % (oldNode, triId, self.surfaceName))
        except IndexError:
            raise IndexError("No tri element %d in surface %s when trying"
                             " to replace node %d on it."
                             % (triId, self.surfaceName, oldNode))
        # update triNodes
        self.triNodes[triId][idx] = newNode

        # update edgeToTri, borderEdges
        if hasattr(self, "edgeToTri"):
            for edge,triList in self.edgeToTri.items():
                if (oldNode in edge) and (triId in triList):
                    triList.remove(triId)
                    if len(triList)==1:
                        self.borderEdges.add(edge)
                    if len(triList)==0:
                        del self.edgeToTri[edge]
                        self.borderEdges.remove(edge)
                    newEdge = set(edge)
                    newEdge.remove(oldNode)
                    newEdge.add(newNode)
                    newEdge = frozenset(newEdge)
                    triList = self.edgeToTri[newEdge]
                    triList.add(triId)
                    if len(triList)==1:
                        self.borderEdges.add(newEdge)
                    if len(triList)==2:
                        self.borderEdges.remove(newEdge)

        # delete obscured members
        self.borderNodes = None
        self.nodesToTri = None
        self.nodeNormals = None
        self.nodeTris=None
        self.nodes=None


    @staticmethod
    def triEdgesOrderedNodesIter(nodesOfThisTri):
        """Returns an iterator of the edges of the given triangle in two forms

        nodesOfThisTri must be a list of node ids

        The return value is a tuple of two elements. Each of those states the
        edge by its node numbers. The first element of the return tuple is a
        frozenset, that means the nodes are not ordered. The second element is
        a tuple of the node ids in the correct order.

        similar to TriangleSurface.triEdgesIter()
        """
        for i1, i2 in ((0,1), (1,2), (2,0)):
            yield (frozenset((nodesOfThisTri[i1], nodesOfThisTri[i2])),
                   (nodesOfThisTri[i1], nodesOfThisTri[i2]))

    @staticmethod
    def triEdgesIter(nodesOfThisTri):
        """Returns an iterator of the edges of the given triangle
        nodesOfThisTri will be converted to a tuple of node ids

        use as in
        myTri = (3, 43, 12)
        for edges in triEdgesIter(myTri):
           ... edges will be frozenset((3, 43)), next time frozenset((43, 12))
           and then frozenset((12, 3))

        see triEdgesOrderedNodesIter if you also need the nodes in the correct
        order.
        """
        nodesOfThisTri = tuple(nodesOfThisTri)
        for i1, i2 in ((0,1), (1,2), (2,0)):
            yield frozenset((nodesOfThisTri[i1], nodesOfThisTri[i2]))


#####################################################################
#####################################################################
#####################################################################
#####################################################################

# Additional import at the end of the file
# This has to be done here, because TriangleSurfaceInVolume needs the
# definition of TriangleSurface.
from bae.surface_in_volume_01 import TriangleSurfaceInVolume


## ---------------------------------------------------------------------------
## convinient constructor from already existing TriangleSurface

def TriangleSurfaceInVolume_fromTriangleSurface(surf, tetNodes):
    """creates a TriangleSurfaceInVolume-object from an already existing
    TriangleSurface-object. Note that the nodeCoords member of the
    TriangleSurface-object already has to contain all nodes for the tets
    as well.

    Don't know if the copy must be deepcopy
    """
    newsurf = TriangleSurfaceInVolume(
        surf.surfaceName,
        surf.triNodes,
        surf.triNodes,
        surf.nodeCoords,
        tetNodes,
        surf.logfile)

    # copy unconditionally
    for attr in ['borderNodes', 'nodesToTri', 'nodeNormals', 'nodeTris',
                 'nodes']:
        setattr(newsurf, attr, copy.deepcopy(getattr(surf, attr)))

    for attr in ['triNormals', 'edgeToTri', 'borderEdges', 'nodesToDoubles',
                 'moebiusStopEdges', 'moebius_flag']:
        if hasattr(surf, attr):
            setattr(newsurf, attr, copy.deepcopy(getattr(surf, attr)))

    return newsurf
