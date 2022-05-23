"""That's basically what used to be called findAtMesh
"""

## ---------------------------------------------------------------------------
__version__ = "1.05"
_version_history_ = """\
1.00 GP new: from findAtmesh 4.03
1.01 GP fixed bugs, split off makeDisJoint function and removed makeDisJoint
    and smoothElsets option from mapElsets (incompatible parameter change)
    added mapWholeVolume
1.02 GP added: mapElems
1.03 GP changed: use KDTree.boxSearch in kdtreePointsInVolMesh
        and also some (incompatible) parameter name changes (no worries)
1.04 GP accelerated kdtreePointsInVolMesh
1.05 GP added MeshToVolMapperCurruptSurface
"""

todo = """
 - complete the version using scipy, using selectattrversion
"""

from collections import defaultdict, deque

from bae.misc_01 import selectattrversion
from bae.abq_model_02 import Model
from bae.log_01 import msg, MsgTicker

numpy = None
scipy = None
# try:
#     import numpy
# except ImportError:
#     numpy = None
#     scipy = None
# else:
#     try:
#         import scipy
#     except ImportError:
#         scipy = None

if scipy:
    from scipy.spatial import KDTree
else:
    from bae.misc_01 import KDTree


class MeshToVolMapper(object):
    """That's basically the findAtMesh4 functionality
    Usage: Suppose MyPrj2014_G01_tetL.inp is the mesh for which we seek elsets,
    volumeMesh.inp and volumeElsets.inp contain the volumes (called "targets"
    in the dark age of total confusion):
     >>> from bae.utils.mapping_01 import MeshToVolMapper
     >>> mapper = MeshToVolMapper("MyPrj2014_G01_tetL.inp")
     >>> mappedElsets = mapper.mapElsets(
     ...     volModel=["volumeMesh.inp", "volumeElsets.inp"] )
     >>> mout = Model()
     >>> mout.elset.update( mappedElsets )
     >>> mout.write("mappedElsets.inp")

    Or even simpler if you've got a whole-volume-model, i.e. the volume is
    defined by *all* elements of volumeMesh.inp:
     >>> from bae.utils.mapping_01 import MeshToVolMapper
     >>> mapper = MeshToVolMapper("MyPrj2014_G01_tetL.inp")
     >>> mout = Model()
     >>> mout.elset["VOLUME"] = mapper.mapWholeVolume("volumeMesh.inp")
     >>> mout.write("mappedElsets.inp")

    Or possibly use L{bae.misc_01.withPickle}():
     >>> from bae.utils.mapping_01 import MeshToVolMapper
     >>> from bae.misc_01 import withPickle
     >>> def getMapper():
     >>>     return MeshToVolMapper("MyPrj2014_G01_tetL.inp")
     >>> mapper = withPickle("MyPrj2014_G01_mapper.pickle", getMapper)
     >>> ...

    @Note: See L{bae.misc_01.withPickle}(..., keepInMem=True) and
    L{bae.misc_01.needsUpdate} for a convenient, efficient and flexible way to
    map different volumes into one and the same mesh.
    """

    def __init__(self, mesh, elems=None):
        """
        @param mesh: can be a L{bae.mesh_01.Mesh}-object or of a derived class
           (e.g. L{bae.abq_model_02.Model}) or a filename (a string).
        @param elems: Anything that L{bae.mesh_01.Mesh.getElCentroids} or
           L{bae.abq_model_02.Model.getElCentroids} accepts as elems parameter.
           States the elements to be tested whether they are within the volume.
           If None then all elements are considered.
        """
        if isinstance(mesh, basestring):
            mesh = Model().read(mesh)
            msg("Finished reading mesh with %d elements and %d nodes."
                % (len(mesh.elNodes), len(mesh.nodeCoords)))

        if len(mesh.nodeCoords)==0 or len(mesh.elNodes)==0:
            raise Exception("The mesh does not contain elements and nodes.")

        msg("Initializing the mapper, ...calculating element centroids...")
        elemPtsDict = mesh.getElCentroids(elems)
        meshElemNbList = sorted(elemPtsDict)
        elemCentroidsList = [elemPtsDict[elem] for elem in meshElemNbList]

        if scipy:
            self.elemLabels = numpy.array(meshElemNbList, dtype=float)
            self.elemCentroids = numpy.array(elemCentroidsList, dtype=float)
            self.kdt = KDTree(self.elemCentroids, leafsize=10)
        else:
            # initialize KDTree on mesh centroids
            self.meshElemNbList = meshElemNbList
            self.elemCentroidsList = elemCentroidsList
            msg("... initializing KDTree ...")
            self.kdt = KDTree(self.elemCentroidsList)
            msg(" done.")


    def kdtreePointsInVolMesh(self, volmesh, elems):
        """For each element in elems find all points in the kdtree that are inside
        the element. Treat all elements in volmesh as tets.

        Could remove all non-tet elements from elems initially. Not implemented
        (yet).
        """

        tolerance = 1E-4
        avg_n = 4
        ptIdsFound = set()
        ticker = MsgTicker("... processing element %s/%d in the volume mesh"
                           % ("%s", len(elems)))
        for cnt, elem in enumerate(elems):

            ticker.msg(cnt+1)
            # corner coordinates
            nodes = volmesh.elNodes[elem]
            coords = [volmesh.nodeCoords[node] for node in nodes[:avg_n]]

            # bounding box for this element
            coordsT = zip(*coords)
            elemBox = [ map(min, coordsT), map(max, coordsT) ]
            ptIdsToCheck = self.kdt.boxSearch(elemBox)

            # ignore points that are already identified
            ptIdsToCheck = set(ptIdsToCheck).difference(ptIdsFound)
            if not ptIdsToCheck:
                continue

            # prepare final test of pts in ptIdsToCheck (i.e. in elemBox)
            x = coords
            # coordinates relative to the first node
            # dx = [[xxi-x0i for xxi, x0i in zip(xx, x[0])]
            #       for xx in x[1:]]
            dx = [[x[1][0]-x[0][0], x[1][1]-x[0][1], x[1][2]-x[0][2]],
                  [x[2][0]-x[0][0], x[2][1]-x[0][1], x[2][2]-x[0][2]],
                  [x[3][0]-x[0][0], x[3][1]-x[0][1], x[3][2]-x[0][2]]]
            # determinant of the "Jakobian"
            # cylcePos = [[0,1,2], [1,2,0], [2,0,1]]
            # detJ = sum(
            #     dx[0][ip[0]]*dx[1][ip[1]]*dx[2][ip[2]]
            #     - dx[2][ip[0]]*dx[1][ip[1]]*dx[0][ip[2]]
            #     for ip in cylcePos)
            detJ = (dx[0][0]*dx[1][1]*dx[2][2]
                    + dx[0][1]*dx[1][2]*dx[2][0]
                    + dx[0][2]*dx[1][0]*dx[2][1]
                    - dx[2][0]*dx[1][1]*dx[0][2]
                    - dx[2][1]*dx[1][2]*dx[0][0]
                    - dx[2][2]*dx[1][0]*dx[0][1])
            if detJ <= 0.0:
                # negative or zero volume, ignore this element
                continue

            # inverse * determinant
            # invJ = [ [a[ip2[1]]*b[ip2[2]] - a[ip2[2]]*b[ip2[1]]
            #           for ip2 in cylcePos]
            #          for a, b in ((dx[ip1[1]], dx[ip1[2]])
            #                       for ip1 in cylcePos)]
            invJ = [[dx[1][1]*dx[2][2] - dx[1][2]*dx[2][1],
                     dx[1][2]*dx[2][0] - dx[1][0]*dx[2][2],
                     dx[1][0]*dx[2][1] - dx[1][1]*dx[2][0]],
                    [dx[2][1]*dx[0][2] - dx[2][2]*dx[0][1],
                     dx[2][2]*dx[0][0] - dx[2][0]*dx[0][2],
                     dx[2][0]*dx[0][1] - dx[2][1]*dx[0][0]],
                    [dx[0][1]*dx[1][2] - dx[0][2]*dx[1][1],
                     dx[0][2]*dx[1][0] - dx[0][0]*dx[1][2],
                     dx[0][0]*dx[1][1] - dx[0][1]*dx[1][0]]]

            # final test of pts in elemBox: in tet?
            for ptIdx in ptIdsToCheck:
                point = self.kdt.pointList[ptIdx]

                # coordinates relative to the first node
                # dp = [pi - x0i for pi, x0i in zip(point, x[0])]
                dp = [point[0]-x[0][0], point[1]-x[0][1], point[2]-x[0][2]]

                # parametric coordinates: xi_p := invJ * dp
                # xip = [float(sum(invJ[i][j]*dp[j] for j in range(3)))/detJ
                #        for i in range(3)]
                xip = [float(invJ[0][0]*dp[0]+invJ[0][1]*dp[1]+invJ[0][2]*dp[2])/detJ,
                       float(invJ[1][0]*dp[0]+invJ[1][1]*dp[1]+invJ[1][2]*dp[2])/detJ,
                       float(invJ[2][0]*dp[0]+invJ[2][1]*dp[1]+invJ[2][2]*dp[2])/detJ]

                shapeFuncs = [1-sum(xip), xip[0], xip[1], xip[2]]
                if all(xi>-tolerance for xi in shapeFuncs):
                    ptIdsFound.add(ptIdx)
                    continue
            # this element: done, save results

        del ticker
        return ptIdsFound

    def mapElsets(self, volModel, volElsetNames=None):
        """Map volumes represented by individual elsets in the volume-mesh.

        @param volModel: List of names of ABAQUS input file for the volume.
        All files are merged to one model. Or just one filename. Or an
        abq_model_02.Model instance.

        @param volElsetNames: List of elset names in the model representing the
        volume ("target") made up from volModel. Each given elset creates
        one volume and finally one elset in the model made up from
        meshFileNames. If not specified: use all elsets in this model.

        @returns: An iterator of (elset name, set of mesh element numbers)-
        tuples, one tuple for each elset. Can easily be cast to the popular
        elset dictioniary: dict(myMapper.mapElsets(volModel))
        """

        #--- get volModel  ... model containing the volume ("target") elsets
        if isinstance(volModel, basestring):
            volModel = Model().read(
                volModel, recognizedBlocks="NODE,ELEMENT,ELSET")
        elif isinstance(volModel, Model):
            pass
        else:
            v2 = Model()
            for fi in volModel:
                v2.read(fi, recognizedBlocks="NODE,ELEMENT,ELSET")
            volModel = v2; del v2
        if len(volModel.nodeCoords)==0 or len(volModel.elNodes)==0:
            raise Exception("The model for the volume(s) does not contain"
                            " elements and nodes.")
        msg("Finished reading model for the volume(s) with %d elements and"
            " %d nodes and %d elsets."
            % (len(volModel.elNodes), len(volModel.nodeCoords),
               len(volModel.elset)))

        #--- determine elsets in the volModel to be considered
        if (volElsetNames is None or len(volElsetNames)==0):
            volElsetNames = sorted(volModel.elset)
            if len(volElsetNames)==0:
                raise Exception("The model for the volume(s) does not contain"
                                " element sets.")
            msg("Will find elements in the volume of all %d element sets of"
                " the volume model." % len(volElsetNames))
        else:
            # remove non-existing elsets but retain old order
            newList = [eln for eln in volElsetNames if eln in volModel.elset]
            if len(newList)==0:
                raise Exception(
                    "The model for the volume(s) does not contain any"
                    " of the element sets specified in volElsetNames.")
            notFoundElsets = [eln for eln in volElsetNames
                              if eln not in volModel.elset]
            volElsetNames = newList
            if len(notFoundElsets)>0:
                msg("WARNING: %d elsets from volElsetNames are not defined in"
                    " the model for the volume(s).")
                if len(notFoundElsets)>10:
                    msg("Here are the first 10 of them:\n%s"
                        % notFoundElsets[:10])
                else:
                    msg("Those are the missing element sets:\n%s"
                        % notFoundElsets)

        #--- filter tet elements in volModel
        allTetsInVol = volModel.shapeEl.get("TET_L", set()).union(
            volModel.shapeEl.get("TET_Q", set()))

        for elCnt, elsetName in enumerate(volElsetNames):
            msg("Processing elset %s (%d/%d)"
                % (elsetName, elCnt+1, len(volElsetNames)))

            # intersect mesh
            elems = allTetsInVol.intersection(volModel.elset[elsetName])
            foundElset = self.mapElems(volModel, elems)

            yield (elsetName, foundElset)

    def mapWholeVolume(self, volMesh):
        """Map volume represented by the whole volume-mesh.

        @param volMesh: name of ABAQUS input file or abq_model_02.Model object

        @returns: a set of mesh element numbers
        """

        #--- get volMesh  ... model containing the volume ("target") elsets
        if isinstance(volMesh, basestring):
            volMesh = Model().read(volMesh)

        if len(volMesh.nodeCoords)==0 or len(volMesh.elNodes)==0:
            raise Exception("The model for the volume(s) does not contain"
                            " elements and nodes.")

        msg("Model for the volume contains %d elements and %d nodes."
            % (len(volMesh.elNodes), len(volMesh.nodeCoords)))


        #--- filter tet elements in volMesh
        allTetsInVol = volMesh.shapeEl.get("TET_L", set()).union(
            volMesh.shapeEl.get("TET_Q", set()))

        # intersect mesh
        return self.mapElems(volMesh, allTetsInVol)

    def mapElems(self, volMesh, elems):
        """Map volume represented by the given elements
        @param volMesh: abq_model_02.Model object, could be a mesh_01.Mesh as
            well. Should only contain tet elements.
        @param elems: Iterable of element numbers of elements defined in
            volMesh.
        @returns: a set of mesh element numbers

        @Note: To ensure that only tet elements are in the given elset do:
         >>> allTetsInVol = volMesh.shapeEl.get("TET_L", set()).union(
         >>>     volMesh.shapeEl.get("TET_Q", set()))
         >>> elems.intersection_update(allTetsInVol)
        """
        ptIdsFound = self.kdtreePointsInVolMesh(volMesh, elems)
        foundElset = set(self.meshElemNbList[idx] for idx in ptIdsFound)
        return foundElset


class MeshToVolMapperCurruptSurface(object):
    """class for mapping corrupt Rhino meshes into an FE-mesh.
    CONSTRUCTION SITE!!! this needs cleanup, interface will change!!

    Current usage, THIS WILL CHANGE!!!
    ==================================

     >>> from bae.utils.mapping_01 import MeshToVolMapperCurruptSurface as M
     >>> from bae.abq_model_02 import Model
     >>> mapper = M("MyPrj2014_G01_tetL.inp", elems="SEARCHSET")
     >>> borderElsets = Model().read("MyBordersElsets.inp", "ELSET").elset
     >>> mout = Model()
     >>> for eln, els in borderElsets.iteritems():
     >>>     mappedElset = mapper.mapVolumeSurf(els)
     >>>     mout.elset["eln"] = mappedElset
     >>> mout.write("mappedElsets.inp")

    The following is where we want to go. Not implemented yet.

    Usage: Suppose MyPrj2014_G01_tetL.inp is the mesh for which we seek elsets,
    volumeMesh.inp and volumeElsets.inp contain the volumes (called "targets"
    in the dark age of total confusion):
     >>> from bae.utils.mapping_01 import MeshToVolMapperCurruptSurface as M
     >>> from bae.surface_03 import TriangleSurface
     >>> mapper = M("MyPrj2014_G01_tetL.inp")
     >>> surf = TriangleSurface.fromSTL("WreckedSurf.stl")
     >>> mappedElsets = mapper.mapVolumeSurf(surf)
     >>> mout = Model()
     >>> mout.elset.update( mappedElsets )
     >>> mout.write("mappedElsets.inp")
    """

    def __init__(self, mesh, elems=None):
        """
        The mesh (the searchset) must not contain any holes. I.e. you should
        not take the quadratic mesh that might contain cohesives!

        @param mesh: can be a L{bae.mesh_01.Mesh}-object or of a derived class
           (e.g. L{bae.abq_model_02.Model}) or a filename (a string).
        @param elems: "Searchset", elements in mesh to consider.
           Anything that L{bae.mesh_01.Mesh.getElCentroids} or
           L{bae.abq_model_02.Model.getElCentroids} accepts as elems parameter.
           States the elements to be tested whether they are within the volume.
           If None then all elements are considered.
        """
        if isinstance(mesh, basestring):
            mesh = Model().read(mesh)
            msg("Finished reading mesh with %d elements and %d nodes."
                % (len(mesh.elNodes), len(mesh.nodeCoords)))

        if len(mesh.nodeCoords)==0 or len(mesh.elNodes)==0:
            raise Exception("The mesh does not contain elements and nodes.")

        if elems:
            self.mesh = mesh.copySubmodelFromElsets(elems)
        else:
            self.mesh = mesh

        # some diagnostic outputs on the elements in the mesh
        tetElems = self.mesh.shapeEl["TET_L"]

        # search for unrecognized element types/shapes
        if isinstance(self.mesh, Model):
            unknowntypes = defaultdict(int)
            for el in set(self.mesh.elNodes).difference(tetElems):
                unknowntypes[self.mesh.elType[el]] += 1
            unknowntypes = sorted(unknowntypes.iteritems(), key=lambda x: x[1])
            if unknowntypes:
                msg("WARNING: Some elements are not being considered, (type,"
                    " number): %s"
                    % (", ".join("%s: %d" % x for x in unknowntypes)))
            del unknowntypes
        else:
            unknownshapes = defaultdict(int)
            for el in set(self.mesh.elNodes).difference(tetElems):
                unknownshapes[self.mesh.elShape[el]] += 1
            unknownshapes = sorted(unknownshapes.iteritems(),key=lambda x:x[1])
            if unknownshapes:
                msg("WARNING: Some elements are not being considered, (shape,"
                    " number): %s"
                    % (", ".join("%s: %d" % x for x in unknownshapes)))
            del unknownshapes

        # store permanently
        self.tetElems = tetElems

        # just an abbreviation
        elNodes = self.mesh.elNodes

        msg("Initializing the mapper, ...determine face tet connectivity...")
        # faceToTet: {set(node1, node2, node3): list of elements}
        faceToTet = defaultdict(list)
        ticker = MsgTicker(
            "Preparing face connectivity data for element %%d/%d"
            % len(self.tetElems))
        for element in self.tetElems:
            ticker.tick()
            for thisFace in self.tetFacesIter(elNodes[element]):
                faceToTet[thisFace].append(element)
        del ticker
        faceToTet.default_factory = None
        if any(len(tets)>2 for tets in faceToTet.itervalues()):
            msg("WARNING: There are more than two tets attached to a single"
                " face.")

        # store permanently
        self.faceToTet = faceToTet

        msg("Initializing the mapper, ...determine outer/free surface...")
        self.tetsOnFreeSurface = [
            tets[0] for face, tets in self.faceToTet.iteritems()
            if len(tets)==1]

        msg("Initializing the mapper, ...   ...")


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

    # def mapVolumeSurf(self, surf):
    #         @param surf: L{bae.surface_03.TriangleSurface}-object

    def mapVolumeSurf(self, borderELset):
        """CONSTRUCTION SITE!!! this needs cleanup, interface will change!!
        Map volume represented by the given triangle surface.

        @param borderELset: set of element numbers defining the border of the
           region to be identified

        @returns: a set of mesh element numbers
        """

        # initialize containers for the main iteration
        tetsNotTouched = set(self.tetElems)
        neighbourTets = deque(self.tetsOnFreeSurface)

        # iterate over elems yet to be checked
        # following loop: find neighbours of currentTet.
        #
        # Tet element currentTet has not yet investigated faces given in
        # the list neighbourFaces. For each of those faces find out if
        # there is a tet behind that we have not yet visited.
        # (currentTet is not in tetsToVisit anymore!)
        #
        # Neigbouring tets not in tetsToVisit anymore have been
        # reached by other paths before. So don't look at them a
        # second time.
        while len(neighbourTets) > 0:

            currentTet = neighbourTets.popleft()
            tetsNotTouched.remove(currentTet)
            neighbourFaces = list(
                self.tetFacesIter(self.mesh.elNodes[currentTet]))

            # check tet if it's on the border
            ############## CONSTRUCTION SITE
            # if on border: continue

            for thisFace in neighbourFaces:

                    # Find any (there should be no more than one) tet
                    # attached to this face that has not been visited yet.
                    # If at all this should be the neighbour tet on the
                    # other side of thisFace opposite to currentTet.
                    # Side note: currentTet does not belong to
                    # tetsNotTouched anymore
                    newTetsSet = tetsNotTouched.intersection(
                        self.faceToTet[thisFace])
                    try:
                        newTet = newTetsSet.pop()
                    except KeyError:
                        # pop from an empty set
                        # ... no other tet on this face
                        continue

                    neighbourTets.append(newTet)

        # return anything that is not "outside"
        return tetsNotTouched

    def mapVolumeSurfaces(self, surfaces, disjoint=True):
        """Map volumes represented by the triangle surfaces. Default is that
        resulting elsets are disjoint.

        @param surfaces: an iterator of
           L{bae.surface_03.TriangleSurface}-objects
        @param disjoint: If True then resulting elsets are made disjoint: The
           first surface takes all, following surfaces can only take what's
           left from previous surfaces.

        @returns: an iterator of elsets: sets of mesh element numbers
        """

        elemsOutside = set(self.tetsOnFreeSurface)
        # iterate over elems yet to be checked


        # return anything that is not "outside"
        return set(self.mesh.elNodes).difference(elemsOutside)


def smoothElset(elNodes, elset):
    """Try to smooth holes in otherwise rather smooth surfaces of the
    given elsets.

    Add those elements from elNodes that have all four corner nodes already
    connected to the given elset to this particular elset. This is meant to
    even out single-element-holes in otherwise rather smooth surfaces of
    elsets.

    Usage in scripts:
     >>> from bae.abq_model_02 import Model
     >>> from bae.utils.mapping_01 import smoothElset
     >>> mesh = Model().read("G03_TETL.inp")   # nodes and elements
     >>> mSets = Model().read("SETS_MAT.inp", recognizedBlocks="ELSET")
     >>> for els in mSets.elset.itervalues()):
     >>>     smoothElset(mesh.elNodes, els)
     >>> mSets.write("SETS_MAT_2.inp")

    This particular task can be performed more efficiently using the different
    algorithm from smoothElset.py.

    @param elNodes: dict {elem: list of node numbers}
    @param elset: set of element numbers. Will be updated with the elements
       intended to smooth the elset.
    """

    # border nodes for each elset
    msg("Preparing border node set.")
    insideNodes = set(node for elem in elset
                      for node in elNodes[elem])
    outsideElems = set(elNodes).difference(elset)
    outsideNodes = set(node for elem in outsideElems
                       for node in elNodes[elem])
    borderNodes = insideNodes.intersection(outsideNodes)
    del outsideNodes, insideNodes
    msg("Found %d nodes on elset borders, %d elements outside original elset."
        % (len(borderNodes), len(outsideElems)))

    # {border node: set(connected tet elements not in elset)}
    # consider only node 1..4 (corner nodes), don't check element type
    nodeToTets = defaultdict(set)
    msg("Building database for node to elem connections.")
    ticker = MsgTicker("... processed tet element %%d/%d"
                       % len(outsideElems))
    for elem in outsideElems:
        for node in set(elNodes[elem][0:4]).intersection(borderNodes):
            nodeToTets[node].add(elem)
        ticker.tick()
    del ticker
    nodeToTets.default_factory = None  # transform to ordinary dict

    # add tet elements with all 4 nodes connected to the current elset
    msg("Adding smoothing elements to elset.")
    nbAddedElems = 0
    elemsOnBorder = set(
        elem for node in borderNodes for elem in nodeToTets[node])
    for elem in elemsOnBorder:
        nbNodesOnBorder = len(borderNodes.intersection(elNodes[elem][0:4]))
        if nbNodesOnBorder == 4:
            elset.add(elem)
            nbAddedElems += 1
    msg("Added %d elements to elset." % nbAddedElems)


def makeDisJoint(elsetsDict, orderedList):
    """Make elsets disjoint by scraping elsets later in the list from
    elements of any elset found earlier in the list. I.e. the first elset in
    the list will be left untouched.

    The dictionary will be modified in place.

    Usage in scripts:
     >>> from bae.abq_model_02 import Model
     >>> from bae.utils.mapping_01 import smoothElset, makeDisJoint
     >>> mesh = Model().read("G03_TETL.inp")   # nodes and elements
     >>> mSets = Model().read("SETS_MAT.inp", recognizedBlocks="ELSET")
     >>> for els in mSets.elset.itervalues()):
     >>>     smoothElset(mesh.elNodes, els)
     >>> orderedElsetList = sorted(mesh.elset)
     >>> makeDisJoint(mesh.elset, orderedElsetList)
     >>> mSets.write("SETS_MAT_2.inp")

    @param elsetsDict: dict {elset name: set of element numbers} or
        L{bae.abq_model_02.container.ElsetsDict}-object
    @param orderedList: list of elset names stating the order, see
        description above.
    @return: None.

    @Note: If you also want to use smoothElset() then makeDisJoint() should
    be called last.
    """

    allFoundsoFar = set()

    for elCnt, elsetName in enumerate(orderedList):
        elset = elsetsDict[elsetName]
        elset.difference_update(allFoundsoFar)
        allFoundsoFar.update(elset)

    return None
