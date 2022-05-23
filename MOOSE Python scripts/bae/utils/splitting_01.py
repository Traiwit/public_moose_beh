"""Functions for splitting a tet mesh at tri-surfaces, inserting cohesive
elements.
"""

__version__ = "1.06"

_version_history_ = """\
1.01 taken from 1_pre/12_interface/addCohsToElsets3.py v3.04
1.02 GP added MeshSplitter
1.03 GP added parameter elemsToAdd as searchset to addCohsToElsets()
1.04 GP incompatible change: renamed MeshSplitter to MeshSplitterTetsOnSurf
    and added new MeshSplitter class as its parent class.
    added new class MeshSplitterTetsAllFaces.
1.05 GP added rejoinSplitMesh
1.06 GP added ignoreLockedWedges-argument for
    MeshSplitterTetsOnSurf.insertWedges. Don't insert locked wedges for
    singular tris. I.e. disregard tris with all nodes on the border of the
    surface (and not on the boundary of the tet mesh). Because the corresponding
    wedges would have bottom and top face attached to the same three nodes and
    could not deform.
"""

todo = """
 - add old surfaceRepair features as separate functions (badBorder, fillHole
   and so on) to MeshSplitterTetsOnSurf
 - add different splitting functions for different purposes to
   MeshSplitterTetsOnSurf:
   . with pwp-wedges with MPCs
   . create contact surface
 - add constant thickness option to MeshSplitter.splitNodes(): so that you can
   really prescribe a certain constant thickness to the faults also at kinks
   and intersections.
   This will not work in all cases because it will not be possible
   geometrically at places. What do we do there?
"""


from collections import defaultdict, deque
from itertools import izip

from bae.abq_model_02 import Model
from bae.surface_03 import TriangleSurface
from bae.misc_01 import findNextRoundNumber
from bae.vecmath_01 import norm, vector_scale, vector_sum, vector_plus
from bae.log_01 import msg, MsgTicker


def addCohsToElsets(
        mesh,
        mSets,
        elsetsNamePattern=[""],
        cohesiveElsetNameTemplate=None,
        joinedElsetNameTemplate="%s_wCOHS",
        separateTetsElsetNameTemplate=None,
        disjoint=True,
        elemsToAdd="COH3D6",
        tetShape="TET"):
    """Add cohesives (or other elements) to element sets.

     - Connectivity based algorithm starting from cohesive nodes
     - Writing results on the fly
     - Cohesives that share three nodes with one individual tet of the current
       elset are considered as belonging to this elset.
     - The created cohesive elsets are disjoint. I.e. a cohesive separating two
       elsets (i.e. if it sits right on the border of the two elsets) will
       belong to the first of the two. (See disjoint option.)

    Possible improvements: sort elsets according to nb of elements and process
    smallest first.

    @param mesh: L{bae.abq_model_02.Model} object containing the mesh including
      tets and cohesives. Elsets in this model are being ignored. This model
      object is not being altered.

    @param mSets: L{bae.abq_model_02.Model} object containing the elsets that
      are supposed to get cohesives. Might be the same object as mesh. This
      model object is not being altered.

    @param elsetsNamePattern: pattern for elsets (the pattern "" stands for all
      elsets)

    @param cohesiveElsetNameTemplate: If you want cohesive-only-elsets then
      specify a name pattern for them.
      E.g. for VOLC (tets) -> VOLC_COHS (corresponding cohs)
      specify "%s_COHS" here.

    @param joinedElsetNameTemplate: If you want joined elsets with both tets
      and cohesives then specify a name pattern for them. If not set to None
      otherwise the default "%s_wCOHS" will be used.

      E.g. for VOLC (tets) -> VOLC_wCohs (tets and corresponding cohs)
      specify "%s_wCohs" here. Specify "%s" if you want the original name.

    @param separateTetsElsetNameTemplate: If you additionally want to keep the
      original tets-only elsets under a different name then specify a name
      pattern for them.
      E.g. for VOLC (w/o cohs) -> VOLC_TETS (w/o cohs) and VOLC (with cohs)
      specify "%s_TETS" here.
      You need those tets-only elset to apply gravity only on tets.
      If you don't want those extra elsets, assign None as value

    @param disjoint: If True then add each cohesive to only one single elset.
      (That's the default for compatibility. And the procedure is faster with
      disjoint=True.) If False then a cohesive might also end up in two elsets
      if it's exactly at the interface between two sequence elsets. This
      might be desired because then the disjoint functionality of makeHistory
      makes the elsets disjoint again in truly chronological order.

      If the original tet-only elsets are intentionally not disjoint (e.g. for
      the "split-elsets-feature" for independent excavation and refill) then
      you should also use disjoint=False in order to assign all cohesives to
      all elsets.

    @param elemsToAdd: Searchset: All elements that shall be investigated
      whether they shall be added to the corresponding elsets.

      Can be a string identifying a particular element type. Or simply an
      iterator of element numbers.

    @returns: L{bae.abq_model_02.Model} object containing all resulting elsets
      as ordered through the ...ElsetNameTemplate arguments.
    """

    if not(cohesiveElsetNameTemplate
           or joinedElsetNameTemplate):
        raise ValueError(
            "Neither joinedElsetNameTemplate nor cohesiveElsetNameTemplate"
            " has been specified. You would not get anything!")

    if isinstance(elemsToAdd, basestring):
        allCohs = mesh.typeEl[elemsToAdd]
        allCohsName = elemsToAdd
    else:
        allCohs = set(elemsToAdd)
        allCohsName = "selected cohesive(?)"

    #allTets = mesh.shapeEl['TET']  # all tets
    allTets = mesh.shapeEl[tetShape]  # all tets
    msg("%d %s elements and %d tets in the mesh."
        % (len(allCohs), allCohsName, len(allTets)))

    # get the elset names
    allElsetNames = mSets.regexpNames("elset", *elsetsNamePattern)
    msg("Found %d elsets to be updated with %ss."
        % (len(allElsetNames), allCohsName))

    # --- Build elem to elset dictionary elemToElset
    # Depending on the disjoint option there are two different outcomes:
    # If disjoint=True then values of elemToElset are single elset names.
    # If elsets intersect in this case, i.e. elements belong to more than one
    # elset then one arbitrary elset is chosen.
    # if disjoint=False then values of elemToElset are lists of elset names.
    if disjoint:
        # values of elemToElset are single elset names
        elemToElset = dict()
        for elsetName in allElsetNames:
            elemToElset.update(
                (elem, elsetName) for elem in mSets.elset[elsetName])
    else:
        # if disjoint=False then values of elemToElset are lists of elset names
        elemToElset = defaultdict(list)
        for elsetName in allElsetNames:
            for elem in mSets.elset[elsetName]:
                elemToElset[elem].append(elsetName)
        elemToElset.default_factory = None
    # final message
    msg("The %s elsets to be updated contain %d tet elements so far."
        % (len(allElsetNames), len(elemToElset)))

    # all nodes on cohesive elements
    cohsNodes = mesh.getNodesFromElset(allCohs)
    msg("Found %d nodes on %s elements." % (len(cohsNodes), allCohsName))
    
    # apply mpc: replace slave nodes in cohsNodes by master nodes
    if (elemsToAdd=="C3D6P" or elemsToAdd=="DC3D6") and len(mesh.mpc)>0:
        mpcLookup = dict( (nds[0],nds[1]) for type,nds in mesh.mpc )
        msg("Within this nodeset, %d slave nodes will be replaced by master nodes"
            %len(cohsNodes & set(mpcLookup)) )
        cohsNodes = set( mpcLookup.get(ndl,ndl) for ndl in cohsNodes )
    else:
        mpcLookup = dict()

    # {cohs node: set(connected tet elements)}
    nodeToTets = defaultdict(set)
    allSeqTets = mSets.getUnionSet("elset",allElsetNames).intersection(allTets)
    msg("Building database for node to tet elem connections.")
    ticker = MsgTicker("processed tet element %%d/%d" % len(allSeqTets))
    for elem in allSeqTets:
        for node in set(mesh.elNodes[elem]).intersection(cohsNodes):
            nodeToTets[node].add(elem)
        ticker.tick()
    del ticker
    nodeToTets.default_factory = None  # transform to ordinary dict

    # assign each cohesive to an elset, ... preparations
    cohsElset = defaultdict(set)
    def findElsetForCohs(cohs, nodes, connectedTets):
        """assign cohesive element to right elset."""

        elem = connectedTets.pop()
        # extra check. More than one tet connected to either side.
        # Should never fail else mesh or algorithm not ok
        if len(connectedTets) and (nodes[0:3]!=nodes[3:6]):
            connectedTets.add(elem)
            msg("WARNING: On one side of cohesive element %d there are more"
                " than one connected tets: %s. Check mesh or script and"
                " result. Choosing one of the tets to proceed. Results not"
                " reliable!" % (cohs, list(connectedTets)))

        # assign cohesive element to right elset(s)
        if disjoint:
            elsetName = elemToElset[elem]
            cohsElset[elsetName].add(cohs)
        else:
            for elsetName in elemToElset[elem]:
                cohsElset[elsetName].add(cohs)
        return

    # assign each cohesive to an elset, ... actual loop
    msg("Assigning %s to elsets." % allCohsName)
    ticker = MsgTicker("processing %s element %%d/%d"
                       % (allCohsName, len(allCohs)))
    for cohs in allCohs:
        ticker.tick()
        nodes = [ mpcLookup.get(ndl,ndl) for ndl in mesh.elNodes[cohs]]
        

        # find tet connected to cohesive lower side (node 1,2,3)
        try:
            connectedTets = nodeToTets[nodes[0]].intersection(
                nodeToTets[nodes[1]], nodeToTets[nodes[2]])
        except KeyError:
            connectedTets = set()

        # assign to elset
        if connectedTets:
            findElsetForCohs(cohs, nodes, connectedTets)
            if disjoint:
                continue

        # find tet connected to cohesive upper side (node 4,5,6)
        try:
            connectedTets = nodeToTets[nodes[3]].intersection(
                nodeToTets[nodes[4]], nodeToTets[nodes[5]])
        except KeyError:
            connectedTets = set()

        # assign to elset
        if connectedTets:
            findElsetForCohs(cohs, nodes, connectedTets)

    del ticker
    msg("Finished assigning %s to elsets." % allCohsName)

    # generate output
    mOut = Model()
    msg("Assemblying elsets for output.")
    ticker = MsgTicker("processing elset %s/%d\n" % ("%d", len(allElsetNames)))
    for elsetName in sorted(allElsetNames):
        ticker.tick()
        elset = mSets.elset[elsetName]
        foundCohs = cohsElset[elsetName]

        # assembling the elset for output
        if separateTetsElsetNameTemplate:
            newElsetName = separateTetsElsetNameTemplate % elsetName
            mOut.elset[newElsetName] = elset
        if cohesiveElsetNameTemplate:
            newElsetName = cohesiveElsetNameTemplate % elsetName
            mOut.elset[newElsetName] = foundCohs
        if joinedElsetNameTemplate:
            newElsetName = joinedElsetNameTemplate % elsetName
            newset = elset.union(foundCohs)
            mOut.elset[newElsetName] = newset

    del ticker

    return mOut


class MeshSplitter(object):
    """Base class for mesh splitter.
    """
    # some iterator function
    @staticmethod
    def tetFacesIter(nodesOfThisTet):
        """Returns an iterator of the faces of the given tet element
        nodesOfThisTet gives the node ids, it is converted to a tuple

        use as in

        >>> myTet = mdb.el_node[5731]  # list of node ids, e.g. [3, 43, 12, 36]
        >>> for faces in surf.tetFacesIter(myTet):

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

        >>> myTet = mdb.el_node[5731]  # list of node ids, e.g. [3, 43, 12, 36]
        >>> for edges in tetEdgesIter(myTet):

        ... edges will be frozenset((3, 43, )), next time
        frozenset((43, 3, 36)) and so on.

        An extra comment on the index ordering used in the for loop below: If
        the points order would stay that of the used tuple, each faces normal)
        (see getTriNormalUsc points inwardly.
        """
        nodesOfThisTet = tuple(nodesOfThisTet)
        for i1,i2 in ((0,1), (0,2), (0,3), (1,2), (1,3), (2,3)):
            yield frozenset((nodesOfThisTet[i1], nodesOfThisTet[i2]))

    # decide on which side of a triangle a tet is connected
    # tetOnTriSideDict: { ordered tet node indices of a face : direction },
    # pos side = 1, neg side = -1
    # use like this:
    #  >>> triOnTetIndices = tuple(tetNodes.index(node)
    #  >>>                         for node in triNodes)
    #  >>> direction = self.tetOnTriSideDict[triOnTetIndices]
    tetOnTriSideDict = dict.fromkeys([
            (0,1,2), (1,2,0), (2,0,1),
            (1,0,3), (0,3,1), (3,1,0),
            (2,1,3), (1,3,2), (3,2,1),
            (0,2,3), (2,3,0), (3,0,2)], 1)
    tetOnTriSideDict.update(dict.fromkeys([
                (0,2,1), (2,1,0), (1,0,2),
                (1,3,0), (3,0,1), (0,1,3),
                (2,3,1), (3,1,2), (1,2,3),
                (0,3,2), (3,2,0), (2,0,3)], -1))


class MeshSplitterTetsOnSurf(MeshSplitter):
    """Tool to split a tet-mesh on triangle surfaces.

    Usage:
     >>> model = Model().read("MYPRJ_meshL_wTris.inp")
     >>> surfaceNames = model.regexpNames("elset", "FAULT")
     >>> splitter = MeshSplitterTetsOnSurf(model, surfaceNames)
     >>> wrongTris = splitter.filterOpenTris()
     >>> nodeSides = splitter.findSplittings()
     >>> newNodeToOld, triSplitNodes = splitter.splitNodes(nodeSides)
     >>> triWedgeTuples = splitter.insertWedges(triSplitNodes)

    @ivar mesh: L{bae.mesh_01.Mesh}-object or L{bae.abq_model_02.Model} object
      containing the mesh including tets, tris and possibly elsets for the
      surfaces. Supplied as argument to the constructor.
    @ivar surf: L{TriangleSurface} object consisting of all splitting
      triangles.
    @ivar tetElems: set of all tet elements in self.mesh
    """

    def __init__(self, mesh, surfaceElsets):
        """
        provides self.L{tetElems}, self.L{surf}

        @param mesh: L{bae.mesh_01.Mesh}-object or L{bae.abq_model_02.Model}-
          object containing the mesh including tets, tris and possibly elsets
          for the surfaces.
        @param surfaceElsets: A dict {elset name: set of element numbers}
          identifying the triangle surfaces that shall split the tet mesh.
          In case of the mesh-parameter being a L{bae.abq_model_02.Model}-
          object surfaceElsets can also be a list of elset names of those
          surfaces.
        """

        # store mesh object
        self.mesh = mesh

        # process surfaceElsets, create TriangleSurface-object
        if isinstance(mesh, Model) and not isinstance(surfaceElsets, dict):
            surfaceElsets = dict((eln, mesh.elset[eln])
                                 for eln in surfaceElsets)
        triElems = set(x for elems in surfaceElsets.itervalues() for x in elems)
        self.surf = TriangleSurface.fromMeshTris(mesh, triElems)
        triElems = set(self.surf.elNodes)

        # some diagnostic outputs on the elements in the mesh
        tetElems = mesh.shapeEl["TET_L"]

        # search for unrecognized element types/shapes
        if isinstance(mesh, Model):
            unknowntypes = defaultdict(int)
            for el in set(mesh.elNodes).difference(tetElems, triElems):
                unknowntypes[mesh.elType[el]] += 1
            unknowntypes = sorted(unknowntypes.iteritems(), key=lambda x: x[1])
            if unknowntypes:
                msg("WARNING: Some elements are not being considered, (type,"
                    " number): %s"
                    % (", ".join("%s: %d" % x for x in unknowntypes)))
            del unknowntypes
        else:
            unknownshapes = defaultdict(int)
            for el in set(mesh.elNodes).difference(tetElems, triElems):
                unknownshapes[mesh.elShape[el]] += 1
            unknownshapes = sorted(unknownshapes.iteritems(),key=lambda x:x[1])
            if unknownshapes:
                msg("WARNING: Some elements are not being considered, (shape,"
                    " number): %s"
                    % (", ".join("%s: %d" % x for x in unknownshapes)))
            del unknownshapes

        # store permanently
        self.tetElems = tetElems

    def filterOpenTris(self):
        """Identify surface tri elements that don't have a tet connected
        to each side.

        It could be that only one tet is connected to one side, the other side
        is free. It could be that no tet is attached at all.

        Or the mesh could be inconsistent with more than one tet attached to
        one side.

        Those tris will be removed from the splitting surface before the
        actual splitting. Otherwise they would result in cohesive elements
        connected to node number zero.

        @Return: a list of tri elements with only one tet attached to it.

        @Note: It's just the number of attached tets that are counted. If two
        tets are connected to the same side this algoritm does not detect it.
        """
        msg("Searching for badly connected tris. Steps (a)-(b)-(c)")
        elNodes = self.mesh.elNodes  # abbreviation

        msg(" (a) Identifying all nodes on the surface...")
        surfNodes = self.surf.getConnectedNodes()
        msg("... found %d nodes on the surface." % len(surfNodes))

        surfFaces = set(frozenset(nodes)
                        for nodes in self.surf.elNodes.itervalues())
        msg(" (a) Identified %d faces on the surface" % len(surfFaces))

        faceToTets = defaultdict(list)
        ticker = MsgTicker("Building tri to tet connectivity data (%%d/%d)"
                           % len(self.tetElems))
        for tet in self.tetElems:
            ticker.tick()
            nodes = surfNodes.intersection(elNodes[tet])
            if len(nodes)>=3:
                nodes = frozenset(nodes)
                if len(nodes)==3 and nodes in surfFaces:
                    faceToTets[nodes].append(tet)
                elif len(nodes)==4:
                    for face in self.tetFacesIter(nodes):
                        if face in surfFaces:
                            faceToTets[face].append(tet)
        del ticker
        faceToTets.default_factory = None  # finalize faceToTets

        ticker = MsgTicker(
            " (b) Searching for tris with incomplete tet connectivity (%%d/%d)"
            % len(surfNodes), printLast=True)
        wrongTris = list()
        for tri, triNodes in self.surf.elNodes.iteritems():
            ticker.tick()
            try:
                attachedTets = faceToTets[frozenset(triNodes)]
            except KeyError:
                wrongTris.append(tri)
                continue
            if len(attachedTets)!=2:
                wrongTris.append(tri)
                continue
            sideIdxSum = sum(
                self.tetOnTriSideDict[tuple(
                    tetNodes.index(node) for node in triNodes)]
                for tetNodes in (elNodes[tet] for tet in attachedTets))
            if sideIdxSum!=0:
                wrongTris.append(tri)
                continue
        del ticker
        msg(" (b) Found %d tris which did not have a tet connected on each"
            " side." % len(wrongTris))
        if len(wrongTris):
            msg(" ... Here are at least some of them:\n%s" % wrongTris[:100])

            # update splitting surface
            self.surf.removeElems(wrongTris)
            msg(" (c) Removed %d badly connected tris from the splitting"
                " surface." % len(wrongTris))
        else:
            msg(" (c) No need to remove badly connected tris from the"
                " splitting surface.")
        return wrongTris

    def findSplittings(self):
        """Identify nodes possibly to be split and chunks of attached tet
        elements that will stay on the same node after splitting.

        @returns: a dict {node: list of (tetsOnThisSide, borderTris)-tuples}
            with tetsOnThisSide being a list of tet elements on this node
            directly connected without tri element in between
            and borderTris being a list of tri elements confining the
            tets-chunk
        """

        # just an abbreviation
        elNodes = self.mesh.elNodes

        # all nodes on the splitting surface / on tri elements
        msg("Identifying all nodes on the surface...")
        surfNodes = self.surf.getConnectedNodes()

        msg("Identifying tet elements connected to tri elements.")
        # set of tet elements connected to the surface by at least one node
        connectedTets = set(
            tet for tet in self.tetElems
            if surfNodes.intersection(elNodes[tet]))
        msg("We have %d nodes on the tri elements and %d tet elements"
            " attached to them." % (len(surfNodes), len(connectedTets)))

        # {face: tri element nb}
        faceToTri = dict(
            (frozenset(elNodes[el]), el)
            for el in self.surf.elNodes)

        # nodeTets: {node id: set of attached tet elements}
        # Take only elements that are also attached to the cutting surface.
        nodeTets = defaultdict(set)
        ticker = MsgTicker(
            "Preparing node to tet connectivity data for element %%d/%d"
            % len(connectedTets))
        for element in connectedTets:
            ticker.tick()
            for node in set(elNodes[element]).intersection(surfNodes):
                nodeTets[node].add(element)
        del ticker
        nodeTets.default_factory = None  # finalize nodeTets

        # Only elements that are also attached to a tri are considered.
        # tetToFace: {surface-Tet-elementID: set(node1, node2, node3)}
        # faceToTet: {set(node1, node2, node3): list of elements}
        tetToFace = dict()
        faceToTet = defaultdict(list)
        ticker = MsgTicker(
            "Preparing face connectivity data for element %%d/%d"
            % len(connectedTets))
        for element in connectedTets:
            ticker.tick()
            thisTetToFace = set()
            for thisFace in self.tetFacesIter(elNodes[element]):
                thisTetToFace.add(thisFace)
                faceToTet[thisFace].append(element)
            tetToFace[element] = thisTetToFace
        del ticker
        faceToTet.default_factory = None
        if any(len(tets)>2 for tets in faceToTet.itervalues()):
            msg("WARNING: There are more than two tets attached to a single"
                " face.")

        # nodeSides will be filled in the following loop:
        # {node: list of (tetsOnThisSide, borderTris)-tuples}
        # tetsOnThisSide: tet elements on this node directly connected without
        #                 tri element in between
        # borderTris: tri elements confining the tets-chunk
        nodeSides = dict()

        ticker = MsgTicker("Calculating splittings of node %%d/%d"
                           % len(surfNodes))
        for node in surfNodes:
            ticker. tick()
            tetsToVisit = set(nodeTets[node])  # make a copy

            #
            sidesList = list()
            nodeSides[node] = sidesList
            while tetsToVisit:

                # initialize new side
                tetsOnThisSide = list()
                borderTris = list()

                # tetsToFindNeighbours is a deque (list) with tuples for each
                # tet yet to be processed to find it's neighbours. Each tuple
                # consists of:
                # (element number, list of faces to possible neighbours)
                # initialize tetsToFindNeighbours with an arbitrary choice from
                # tetsToVisit
                lastTet = tetsToVisit.pop()
                tetsToFindNeighbours = deque((
                    (lastTet,
                     set(self.tetFacesIter(self.mesh.elNodes[lastTet]))), ))

                # following loop: find neighbours of lastTet without crossing
                # a triangle element. When reaching a triangle element store
                # its number and orientation towards the tet.
                #
                # Tet element lastTet has not yet investigated faces given in
                # the list neighbourFaces. For each of those faces find out if
                # there is a tet behind that we should assign to the current
                # side.
                # (lastTet is not in tetsToVisit anymore!)
                #
                # Neigbouring tets not in tetsToVisit anymore have been
                # reached by other paths before. So don't look at them a
                # second time.
                while len(tetsToFindNeighbours) > 0:

                    lastTet, neighbourFaces = tetsToFindNeighbours.popleft()
                    tetsOnThisSide.append(lastTet)
                    for thisFace in neighbourFaces:
                        # is this face on the current node?
                        if node not in thisFace:
                            # no, connection away from the node, ignore
                            continue

                        try:
                            splittingTri = faceToTri[thisFace]
                        except KeyError:
                            # ordinary connection, ... go on to the next tet

                            # Find any (there should be no more than one) tet
                            # attached to this face that has not been assigned
                            # to either side yet. If at all this should be the
                            # neighbour tet on the other side of thisFace
                            # opposite to lastTet.
                            # Side note: lastTet does not belong to
                            # tetsToVisit anymore
                            newTetsSet = tetsToVisit.intersection(
                                faceToTet[thisFace])
                            try:
                                newTet = newTetsSet.pop()
                            except KeyError:
                                # pop from an empty set
                                # ... no other tet on this face
                                continue

                            # find free faces and recursivly find neighbour tets
                            newFaces = set(tetToFace[newTet])
                            newFaces.remove(thisFace)
                            tetsToFindNeighbours.append((newTet, newFaces))
                            tetsToVisit.remove(newTet)

                        else:
                            # stop at this face
                            # ... store the tri
                            direction = self.tetOnTriSideDict[
                                tuple(elNodes[lastTet].index(node)
                                      for node in elNodes[splittingTri])]
                            borderTris.append((splittingTri, direction))
                            continue

                # finished collecting all connected tets
                sidesList.append((tetsOnThisSide, borderTris))

        # end of loop: for node in surfNodes
        del ticker
        return nodeSides

    def splitNodes(self, nodeSides, nodeDistOffset=0.0, nodeNbOffset=100000):
        """Create new nodes as split-offs and attach tet elements to the new
        nodes.

        Updates self.mesh.nodeCoords.

        @param nodeSides: As from L{findSplittings}
        @param nodeDistOffset: Distance from old nodes to new split-off nodes
        @param nodeNbOffset: New node numbers will start at the next multiple
           of this value after the highest current node number. Plus one.
           For example 100000 might give 3200001 or 100001 as first new node.

        @returns: (newNodeToOld, triSplitNodes)-tuple.
        newNodeToOld is a dict that only contains really new nodes. If you
        want the old node for an assumed new node regardless if it has been
        split actually then use the term "newNodeToOld.get(oldNode, oldNode)".
        triSplitNodes is a dict {triElemNb: wedgeElemNodes},
        wedgeElemNodes is a list of six node numbers that could be used for
        a linear wedge element connecting the tet element sides that have been
        split by the algorithm. In Abaqus node numbering order.
        """

        # abbreviation
        nodeCoords = self.mesh.nodeCoords
        elNodes = self.mesh.elNodes

        # connection data gathered in the following main loop
        newNodeToOld = dict()
        triSplitNodes = dict(  # {tri elem nb: nodes for splitting wedge}
            (tri, [0]*6) for tri in self.surf.elNodes)

        triNormal = self.surf.getTriNormal()

        # for next new node number
        lastNewNode = findNextRoundNumber(
            max(self.mesh.nodeCoords), nodeNbOffset)

        ticker = MsgTicker("splitting off new nodes %%d/%d" % len(nodeSides))
        for oldNode in sorted(nodeSides):
            ticker.tick()
            sides = nodeSides[oldNode]
            oldCoords = nodeCoords[oldNode]

            # loop over compact tet chunks == over new nodes
            for tetsOnThisSide, borderTris in sides:

                # initialize new node
                if len(sides)==1:
                    # no splitting here, keep old node
                    newNode = oldNode
                else:
                    lastNewNode += 1
                    newNode = lastNewNode

                    # relation new to old node
                    newNodeToOld[newNode] = oldNode

                    # get direction to move node
                    try:
                        direction = norm(vector_sum(*(
                            vector_scale(triNormal[splittingTri], dirFact)
                            for splittingTri, dirFact in borderTris)))
                    except IndexError:
                        direction = [1.,0.,0.]
                        msg("WARNING: Unresolved error in"
                            " bae.utils.splitting_01.MeshSplitterTetsOnSurf.splitNodes."
                            " borderTris empty on node to split. This needs debugging."
                            " But the result is likely correct, nevertheless!")
                        msg("Node %d with %d borderTris, tetsOnThisSide: %s" %
                            (oldNode, len(borderTris), tetsOnThisSide))

                    # define new node
                    newCoords = vector_plus(oldCoords, vector_scale(
                        direction, nodeDistOffset))
                    nodeCoords[newNode] = newCoords

                    # attach tets to new node
                    for tet in tetsOnThisSide:
                        nodes = elNodes[tet]
                        nodeIdx = nodes.index(oldNode)
                        nodes[nodeIdx] = newNode

                # store new node in triSplitNodes
                for splittingTri, dirFact in borderTris:
                    nodeIdx = elNodes[splittingTri].index(oldNode)
                    if dirFact==1:
                        nodeIdx += 3  # pos. side goes on upper side of wedge
                    triSplitNodes[splittingTri][nodeIdx] = newNode

        del ticker
        return newNodeToOld, triSplitNodes

    def insertWedges(self, triSplitNodes, elType="WEDGE_L",
                     elemNbOffset=100000, ignoreLockedWedges=True):
        """Inserts linear wedge elements in the places indicated by the
        triSplitNodes result item of L{splitNodes}.

        Updates self.mesh.elNodes (and elShape, shapeEl, elType and typeEl).

        @param triSplitNodes: as from L{splitNodes}:
           a dict {triElemNb: wedgeElemNodes},
           wedgeElemNodes is a list of six node numbers that could be used for
           a linear wedge element connecting the tet element sides that have
           been split by the algorithm. In Abaqus node numbering order.

        @param elType: element shape or -type of the elements to be inserted.
           If L{self.mesh<mesh>} is a L{bae.mesh_01.Mesh}-object it would
           be "WEDGE_L". If L{self.mesh<mesh>} is a L{bae.abq_model_02.Model}
           object it would typically be "COH3D6".

           If elType=="COH3D6P" then add cohesive elements with 9 nodes: 1..3
           lower side, 4..6 upper side, 7..9 middle nodes. The middle nodes
           being the nodes of the old triangle elements.

        @param elemNbOffset: Element numbers for the new wedges will start at
           the next multiple of this value after the highest current element
           number. Plus one.
           For example 100000 might give 3200001 or 100001 as first new number.

        @param ignoreLockedWedges: disregard wedges with top and bottom face
           attached to the same three nodes, i.e. nodes[:3]==nodes[3:]

        @return: a list of (tri elem, new wedge elem)-tuples. (element numbers)
        """

        if ignoreLockedWedges:
            try:
                shape = self.mesh.elTypeToShapeSep[elType]
            except (AttributeError, KeyError):
                shape = elType
            N = self.mesh.elShapeToNbCornerNodes[shape] / 2

            triSplitNodes = dict(
                (tri, nodes)
                for tri, nodes in triSplitNodes.iteritems()
                if nodes[:N]!=nodes[N:])

        offset = findNextRoundNumber(max(self.mesh.elNodes), elemNbOffset) + 1
        triWedgeTuples = [(tri, offset+i)
                          for i, tri in enumerate(sorted(triSplitNodes))]
        if elType=="COH3D6P":
            self.mesh.updateElems(dict(
                (wedge, triSplitNodes[tri]+self.mesh.elNodes[tri])
                for tri, wedge in triWedgeTuples), elType)
        else:
            self.mesh.updateElems(dict(
                (wedge, triSplitNodes[tri])
                for tri, wedge in triWedgeTuples), elType)

        return triWedgeTuples


def splitMesh(mesh, surfaceElsets, removeTris=True,
              nodeNbOffset=100000, elemNbOffset=100000,
              nodeDistOffset=0.0, elType="COH3D6",
              outputFileName=None):
    """Split mesh on surfaces given by triangle elements.

    @param mesh: might be a L{bae.abq_model_02.Model} or L{bae.mesh_01.Mesh}
       instance or a file name for an Abaqus input file to be read.

       If mesh is not an abq_model_02.Model-object or a file name then no
       elsets will be created. This might be of limited practical use.
    @param surfaceElsets: supplies the triangle elements and associated names
       at which the mesh is to be split. It might be a dict {elset name:
       set of element numbers}.
       Or --if the mesh argument is a L{bae.abq_model_02.Model} instance-- it
       might be an elset name pattern that will be passed to
       L{bae.abq_model_02.Model.regexpNames} to get the actual elset names.
       Or --under the same condition-- it might be a list of elset names.
    @param removeTris: if True then remove original surface tris from the model
       This implies a cleanup that removes all nodes not attached to any
       element.
    @param nodeNbOffset: New node numbers will start at the next multiple
       of this value after the highest current node number. Plus one.
       For example 100000 might give 3200001 or 100001 as first new node.
    @param elemNbOffset: Element numbers for the new wedges will start at
       the next multiple of this value after the highest current element
       number. Plus one.
       For example 100000 might give 3200001 or 100001 as first new number.
    @param nodeDistOffset: Distance from old nodes to new split-off nodes.

    @param elType: element shape or -type of the elements to be inserted.
       If the mesh-argument is a L{bae.mesh_01.Mesh}-object it would
       typically be "WEDGE_L".
       If the mesh-argument is a L{bae.abq_model_02.Model}-object it would
       typically be "COH3D6".

       If elType=="COH3D6P" then add cohesive elements with 9 nodes: 1..3
       lower side, 4..6 upper side, 7..9 middle nodes. The middle nodes
       being the nodes of the old triangle elements.

    @param outputFileName: If not specified and mesh is given as a file name
       then a default output file name will be derived from the
       input name (mesh argument) and elset name pattern (surfaceElsets
       argument).

       If not specified and mesh is not given as a file name then no output
       file will be written. If mesh is a L{bae.abq_model_02.Model}-object
       then this object will be modfified and can be exported afterwards.
    """

    # some preparation: read mesh, identify surface elsets, output name
    if isinstance(mesh, basestring):
        if not(outputFileName):
            outputFileName = (
                "%s_%s.inp" % (mesh.split(".inp")[0],surfaceElsets))
        mesh = Model().read(mesh)

    if isinstance(mesh, Model):
        if isinstance(surfaceElsets, basestring):
            surfaceElsets = mesh.regexpNames("elset", surfaceElsets)
        if not isinstance(surfaceElsets, dict):
            surfaceElsets = dict(
                (name, mesh.elset[name])
                for name in surfaceElsets)

    # start splitting...
    splitter = MeshSplitterTetsOnSurf(mesh, surfaceElsets)
    # diagnostic outputs
    msg("Number of elements: %d total, %d surface tris, %d tets"
        % (len(mesh.elNodes), len(splitter.surf.elNodes),
           len(splitter.tetElems)))
    if isinstance(mesh, Model):
        msg("Number of elsets: %d total, %d for splitting."
            % (len(mesh.elset), len(surfaceElsets)))
    else:
        msg("Number of elsets for splitting: %d."
            % len(surfaceElsets))
    msg('Mesh bounding box: %s' % mesh.getBoundingBox())

    # remove badly connected tris (open tris not having one tet on each side)
    splitter.filterOpenTris()

    msg("Identifying nodes to be split and chunks of attached tet elements.")
    nodeSides = splitter.findSplittings()
    msg("Investigated %d nodes. Found %d nodes to be split."
        % (len(nodeSides),
           sum(1 for sides in nodeSides.itervalues() if len(sides)>1)))

    msg("Splitting nodes.")
    newNodeToOld, triSplitNodes = splitter.splitNodes(
        nodeSides, nodeDistOffset=nodeDistOffset, nodeNbOffset=nodeNbOffset)
    if not newNodeToOld:
        raise ValueError("ERROR: Nothing to split here.")
    msg("Created %d new nodes (%d...%d) and found %d places to insert"
        " cohesives."
        % (len(newNodeToOld), min(newNodeToOld), max(newNodeToOld),
           len(triSplitNodes)))
    oneSideTris = [tri for tri, nodes in triSplitNodes.iteritems()
                   if 0 in nodes]
    if oneSideTris:
        msg("ERROR: This should not happen because 'open tris' were filtered"
            " out before: There are %d tris with no tet on one side:\n%s"
            % (len(oneSideTris), oneSideTris))
        msg("Cleaning this up... (You'd better check the resulting mesh.)")

        msg("... Collecting nodes to collapse.")
        wrongNewNodes = set()
        for tri in oneSideTris:
            wrongNewNodes.update(triSplitNodes[tri])
        wrongNewNodes.remove(0)

        msg("... Collecting attached tets.")
        wrongConnectedTets = [
            tet for tet in splitter.tetElems
            if wrongNewNodes.intersection(splitter.mesh.elNodes[tet])]

        msg("... Reattaching tet nodes.")
        for el in wrongConnectedTets:
            splitter.mesh.elNodes[el] = [
                newNodeToOld.get(node, node)
                for node in splitter.mesh.elNodes[el]]

        msg("... Cleaning up triSplitNodes, newNodeToOld.")
        for tri in oneSideTris:
            del triSplitNodes[tri]
        for node in wrongNewNodes:
            del newNodeToOld[node]

    msg("Creating cohesive elements.")
    triWedgeTuples = splitter.insertWedges(
        triSplitNodes, elType=elType, elemNbOffset=elemNbOffset)
    msg("Created %d cohesive elements." % len(triWedgeTuples))

    if isinstance(mesh, Model):
        msg("Creating distinct elsets for the cohesives.")
        triToWedge = dict(triWedgeTuples)
        for surfaceElsetName, triElems in surfaceElsets.iteritems():
            newElsetName = "COHSV_%s" % surfaceElsetName
            newElset = set(triToWedge.get(tri) for tri in triElems)
            newElset.discard(None)  # remove Nones from not converted tris
            mesh.elset[newElsetName] = newElset
        msg("Created %d new elsets." % len(surfaceElsets))

    if removeTris:
        if isinstance(mesh, Model):
            nbRemovedTris = 0
            empty = set()
            allTris = (mesh.shapeEl.get("TRI_L", empty)
                       | mesh.shapeEl.get("TRI_Q", empty)
                       | mesh.shapeEl.get("TRI", empty))
            for elsetName, elset in surfaceElsets.iteritems():
                elsetRemove = elset.intersection(allTris)
                if len(elsetRemove) == len(elset):
                    del mesh.elset[elsetName]
                mesh.removeElems(elsetRemove)
                nbRemovedTris += len(elsetRemove)
            msg("Removed %d surface tri elsets and elements from the model,"
                " leaving %d tri elements." % (nbRemovedTris, len(allTris)))

        msg("Removing unused nodes from the mesh...")
        allNodes = mesh.getConnectedNodes()
        removeNodes = set(mesh.nodeCoords).difference(allNodes)
        for node in removeNodes:
            del mesh.nodeCoords[node]
        msg("... Finished removing unused nodes.")

    if outputFileName:
        msg("Writing result to %s" % outputFileName)
        mesh.write(outputFileName)
        msg("fini.")


def rejoinSplitMesh(model, cohsElset=None):
    """Remove the given cohesive elements from the mesh and re-join the mesh
    accordingly. Modifies the model in place, no return value.

    It's very likely that unconnected nodes are left behind. They are not being
    deleted by this method. For clean up use
    L{bae.mesh_01.Mesh.getConnectedNodes}():
     >>> rejoinSplitMesh(mesh)
     >>> allNodes = mesh.getConnectedNodes()
     >>> removeNodes = set(mesh.nodeCoords).difference(allNodes)
     >>> for node in removeNodes:
     ...    del mesh.nodeCoords[node]

    @param model: A L{bae.abq_model_02.Model} or L{bae.mesh_01.Mesh}
       instance.
    @param cohsElset: Cohesive elements to remove. May be a set of element
       numbers or an element set name or a list of element set names, anything
       that Model.getUnionSet accepts as input.

       If not given or None then take model.shapeEl["WEDGE"] as default.

    @Note: The same algorithm is also implemented in
       L{bae.surface_03.ElemFaceSurface.addFreeSurfFromSplitMesh}.
    """
    elNodes = model.elNodes  # abbreviation

    # determine / check cohsElset
    if cohsElset is None:
        cohsElset = model.shapeEl["WEDGE"]
    elif isinstance(cohsElset, set):
        pass
    elif isinstance(model, Model):
        cohsElset = model.getUnionSet("elset", cohsElset)

    # determine / check volElems
    volElems = model.shapeEl["TET"]
    if not volElems:
        msg("WARNING: rejoinSplitMesh did not find any elements of shape TET"
            " to re-join.")
        return

    msg("rejoinSplitMesh: Determining cohsNodes.", debugLevel=1)
    cohsNodes = (elNodes[cohs] for cohs in cohsElset)

    msg("rejoinSplitMesh: Preparing oldToNew.")
    # create node pairs, sorted ascending
    # Some notes on the algorithm:
    # - cohsNodePairs will be iterated from small node numbers to large
    # - assignments oldToNew will assign from n2/large/old to n1/small/new
    #   node numbers
    # - Therefore in most cases the assignment destination (n1/new) in
    #   oldToNew will not have to be changed later in the process of
    #   assembling it because we won't find n2->n1 with smaller n2 anymore.
    #   Because of the ascending order of cohsNodePairs and the node tuples
    #   in it.
    # - There is one exception: If a new n2->n1 appears with n2/large/old
    #   already reassigned earlier, i.e. we have an earlier assignment
    #   n2->n1' already stored in oldToNew. In this case we store n1->n1'
    #   (or the inverse n1'->n1 such that we maintain the rule to assign
    #   from large to small). Now we also have to check all earlier
    #   oldToNew contents for assignments to a destination that is just
    #   about to be reassigned. I.e. find assignments n2"->n1 and replace
    #   their destinations n1 by n1'.
    cohsNodePairs = sorted(set(
        tuple(sorted(twoNodes))
        for twoNodes in (
            (nodes[i+len(nodes)/2], nodes[i])
            for nodes in cohsNodes
            for i in range(len(nodes)/2))
        if twoNodes[0]!=twoNodes[1]))
    msg("rejoinSplitMesh: ... found %d node pairs..." % len(cohsNodePairs),
        debugLevel=1)

    # higher node number shall be replaced by lower
    oldToNew = dict()
    newToOld = defaultdict(list)
    for n1, n2 in cohsNodePairs:
        n1 = oldToNew.get(n1, n1)
        try:
            # check if the assignment destination n1/new has already been
            # reassigned earlier.
            n2 = oldToNew[n2]
        except KeyError:
            pass
        else:
            # higher node already assigned earlier to oldToNew
            n1, n2 = sorted((n1, n2))
            # n2 might be already assigned as lower number (destination in
            # oldToNew). Check all earlier assignment destinations:
            try:
                nn2List = newToOld[n2]
            except KeyError:
                pass
            else:
                for nn2 in nn2List:
                    oldToNew[nn2] = n1
                del newToOld[n2]
                newToOld[n1].extend(nn2List)
        # in any case, store result: assign higher node to lower node
        oldToNew[n2] = n1
        newToOld[n1].append(n2)
    msg("rejoinSplitMesh: ... found %d pairs of duplicate nodes."
        % len(oldToNew), debugLevel=1)

    # reattach volume elements using oldToNew
    msg("rejoinSplitMesh: Reattaching volume elements.", debugLevel=1)
    elNodes.update(
        (thisElem, [oldToNew.get(node, node) for node in elNodes[thisElem]])
        for thisElem in volElems)

    # remove cohs elements
    msg("rejoinSplitMesh: Removing %d cohesive elements." % len(cohsElset),
        debugLevel=1)
    model.removeElems(cohsElset)
    msg("rejoinSplitMesh: ...fini.", debugLevel=1)


class MeshSplitterTetsAllFaces(MeshSplitter):
    """Tool to split a tet-mesh on all tet faces

    Usage:
     >>> model = Model().read("MYPRJ_meshL.inp")
     >>> splitter = MeshSplitterTetsAllFaces(model)
     >>> newNodeToOld, faceSplitNodes = splitter.splitNodes()
     >>> faceWedgeTuples = splitter.insertWedges(faceSplitNodes)
     >>> # clean up nodes not connected anymore
     >>> allNodes = model.getConnectedNodes()
     >>> removeNodes = set(model.nodeCoords).difference(allNodes)
     >>> for node in removeNodes:
     ...    del model.nodeCoords[node]

    @ivar mesh: L{bae.mesh_01.Mesh}-object or L{bae.abq_model_02.Model} object
      containing the mesh including tets, tris and possibly elsets for the
      surfaces. Supplied as argument to the constructor.
    """

    def __init__(self, mesh):
        """
        @param mesh: L{bae.mesh_01.Mesh}-object or L{bae.abq_model_02.Model}-
          object containing the mesh including tets, tris and possibly elsets
          for the surfaces.
       """

        # store mesh object
        self.mesh = mesh

    def splitNodes(self, nodeNbOffset=100000):
        """Create new nodes as split-offs and attach tet elements to the new
        nodes.

        Updates self.mesh.nodeCoords.

        @param nodeNbOffset: New node numbers will start at the next multiple
           of this value after the highest current node number. Plus one.
           For example 100000 might give 3200001 or 100001 as first new node.

        @returns: (newNodeToOld, faceSplitNodes)-tuple
        newNodeToOld is a dict that contains all new nodes.
        faceSplitNodes is a dict {face: wedgeElemNodes}
        face is a frozenset of three node numbers. Those are the *old* nodes
        that the tet has been connected to before this function call.
        wedgeElemNodes is a list of six node numbers that could be used for
        a linear wedge element connecting the tet element sides that have been
        split by the algorithm. In Abaqus node numbering order.
        """

        # some diagnostic outputs on the elements in the mesh
        tetElems = self.mesh.shapeEl["TET_L"]
        msg("Number of elements: %d total, found %d tets."
            % (len(self.mesh.elNodes), len(tetElems)))

        # just an abbreviation
        elNodes = self.mesh.elNodes
        nodeCoords = self.mesh.nodeCoords

        # Only elements that are also attached to a tri are considered.
        # faceToTet: {set(node1, node2, node3): list of elements}
        faceToTet = defaultdict(list)
        ticker = MsgTicker(
            "Preparing face connectivity data for element %%d/%d"
            % len(tetElems))
        for element in tetElems:
            ticker.tick()
            for thisFace in self.tetFacesIter(elNodes[element]):
                faceToTet[thisFace].append(element)
        del ticker
        faceToTet.default_factory = None

        wrongFaces = [face for face, tets in faceToTet.iteritems()
                      if len(tets)>2]
        if wrongFaces:
            msg("WARNING: There %d faces with more than two tets attached."
                " They will be ignored" % len(wrongFaces))
            for face in wrongFaces:
                del faceToTet[face]

        singleFaces = [face for face, tets in faceToTet.iteritems()
                       if len(tets)<2]
        for face in singleFaces:
            del faceToTet[face]
        msg("Found %d exteriour faces. %d internal faces left that will be"
            " split." % (len(singleFaces), len(faceToTet)))

        # initialize result variables
        newNodeToOld = dict()
        faceSplitNodes = dict(  # {tri elem nb: nodes for splitting wedge}
            (face, [0]*6) for face in faceToTet)

        # for next new node number
        nextNewNode = findNextRoundNumber(
            max(self.mesh.nodeCoords), nodeNbOffset) + 1

        # pos. side tet nodes go on upper side of wedge
        dirFactToWedgeNodeIdxOffset = {1:3, -1:0}

        ticker = MsgTicker("Splitting off new nodes, element %%d/%d"
                           % len(tetElems))
        for element in sorted(tetElems):
            ticker.tick()

            # create new node numbers
            oldNodes = elNodes[element]
            newNodes = [nextNewNode+i for i in range(len(oldNodes))]
            nextNewNode += len(oldNodes)

            # update newNodeToOld
            newNodeToOld.update(
                (newNode, oldNode)
                for newNode, oldNode in izip(newNodes, oldNodes))

            # update faceSplitNodes
            for face in self.tetFacesIter(elNodes[element]):
                if face not in faceToTet:
                    continue

                wedgeNodes = faceSplitNodes[face]
                directedFace = sorted(face)
                wedgeNodeIdxOffset = dirFactToWedgeNodeIdxOffset[
                    self.tetOnTriSideDict[
                        tuple(oldNodes.index(node) for node in directedFace)]]
                for i, oldNode in enumerate(directedFace):
                    newNode = newNodes[oldNodes.index(oldNode)]
                    wedgeNodes[i+wedgeNodeIdxOffset] = newNode

            # attach tet to new nodes
            elNodes[element] = newNodes

        del ticker

        msg("Updating node coordinates...")
        nodeCoords.update(
            (newNode, nodeCoords[oldNode])
            for newNode, oldNode in newNodeToOld.iteritems())

        msg("Created %d new nodes and splitted %d faces."
            % (len(newNodeToOld), len(faceSplitNodes)))
        return (newNodeToOld, faceSplitNodes)

    def insertWedges(self, faceSplitNodes, elType="COH3D6",
                     elemNbOffset=100000):
        """Inserts linear wedge elements in the places indicated by the
        faceSplitNodes result item of L{splitNodes}.

        Updates self.mesh.elNodes (and associates elShape and shapeEl).

        @param faceSplitNodes: as from L{splitNodes}

        @param elType: element shape or -type of the elements to be inserted.
           If L{self.mesh<mesh>} is a L{bae.mesh_01.Mesh}-object it would
           be "WEDGE_L". If L{self.mesh<mesh>} is a L{bae.abq_model_02.Model}
           object it would typically be "COH3D6".

        @param elemNbOffset: Element numbers for the new wedges will start at
           the next multiple of this value after the highest current element
           number. Plus one.
           For example 100000 might give 3200001 or 100001 as first new number.

        @return: a list of (face, new wedge elem number)-tuples.
           face is a frozenset of three node numbers. Those are the *old* nodes
           that the tet has been connected to before the call to L{splitNodes}.
        """

        offset = findNextRoundNumber(max(self.mesh.elNodes), elemNbOffset) + 1
        faceWedgeTuples = [(face, offset+i)
                           for i, face in enumerate(faceSplitNodes)]
        self.mesh.updateElems(dict(
            (wedge, faceSplitNodes[face])
            for face, wedge in faceWedgeTuples), elType)

        return faceWedgeTuples
