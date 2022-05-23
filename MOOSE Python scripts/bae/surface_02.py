"""surface_02.py
Will soon be DEPRECATED, preferably use L{bae.surface_03}

surface class
"""

_version_ = "2.12"

_version_history_ = """\
Versions:
=========

2.00 GP new: derived from surface_01.py version 1.05
2.01 GP changed: messages to logfile and sets to a debugmodel for:
         getSimpleSurfacesFromMeshTris, findSingleTris, checkSimplyConnected
     fixed: joinParts option of getSimpleSurfacesFromMeshTris fixed the actual
         joining of surface parts
2.02 GP added: asAbqModel()
2.03 GP modif: exportAsRaw() accepts an open file object and a new objectName
        argument.
2.04 GP added: exportAsVtk()
2.05 GP switched to bae.log_01, added SurfaceFromCohesives class
2.06 GP added: ???
2.07 GP added: insertSurface, colour and transparency options to exportAsIv
2.08 GP added: getRhinoMeshGeom
2.09 GP added: ElemFaceSurface
2.10 GP: getRhinoMeshGeom now stores element and node numbers as texture
         coordinates
2.11 GP: addFreeSurfFromElset now accepts elset=None argument
2.12 AF: changed triElementFaces to have outward face normal directions
"""

## ---------------------------------------------------------------------------
_bugs_ = """\
known bugs/problems:
====================
- needs python 2.4, if in doubt use abaqus python
- exportAsVtk() issues a warning from the underlying pyvtk library that can be
  ignored. It complains about the mesh having no data attribute, which
  apparently is fine for paraview. Suppressing this warning is not easily
  possible. Would have to modify the pyvtk library.

ideas / todo:
=============
- getSimpleSurfacesFromMeshTris(): at edges with more than two tris instead of
  splitting it up in any case accept a connection with an angle of less then
  5 deg if all other intersections are > 30 deg.
"""

from collections import defaultdict, deque  # deque = double ended queue
from itertools import izip
from sys import stdout
from math import cos, pi, sqrt
import struct

from bae.vecmath_01 import vector, vector_plus, vector_modif_add,\
    vector_scale, norm, dot, cross, length, dist
from bae.abq_model_02 import Model, checkModelModuleVersion, \
    container as abqModelInternal
from bae.mesh_01 import Mesh
from bae.misc_01 import BoundingBox
from pyvtk import UnstructuredGrid, VtkData
from bae.log_01 import msg
# IMPORTANT: look at the end of the file, there is yet another import

# import Rhino if possible
try:
    import Rhino, System
except ImportError:
    Rhino = None


class _SurfaceConnectionError(Exception):
    """raised when there is an error in the connection of the elements
    forming the surface"""
    pass

## ---------------------------------------------------------------------------

class Surface(object):
    pass

## ---------------------------------------------------------------------------

class TriangleSurface(Surface):
    """
    A surface made of triangles.

    @ivar triNodes: is a dictionary: {element number:  list of node numbers}.
      It contains only the elements that actually form the surface.
    @ivar nodeCoords: is a reference to the dictionary: {node number:
      coordinates of the node} as given to the constructor.
      If you move the original nodes, the surface moves also.
    """
    def __init__(self, nodeCoords=None, triNodes=None, logfile=None):
        """
        @param nodeCoords: optional, initialize self.nodeCoords with this
        @param triNodes: optional, initialize self.triNodes with this
        @param logfile: deprecated argument, only for compatibility reasons.
        No effect.

        @note: The triangles are assumed to be simply connected
        (no T-junctions, all connected, no separate parts).
        """

        self.logfile = stdout
        del logfile

        # the following are always there: self.nodeCoords, self.triNodes
        if nodeCoords is None:
            self.nodeCoords = dict()
        else:
            assert isinstance(nodeCoords, dict)
            self.nodeCoords = nodeCoords

        if triNodes is None:
            self.triNodes = dict()
        else:
            assert isinstance(triNodes, dict)
            self.triNodes = triNodes

        # the following are None as long as they are not calculated yet
        # don't query these variables, use the corresponding get...() method,
        # ie. for self.edgeToTri use self.getEdgeToTri()
        # if any operation renders those data invalid, set the correspondig
        # variable (eg. self.edgeToTri) to None and it will be recalculated
        # automatically
        self.edgeToTri = None
        self.borderEdges = None
        self.borderNodes = None
        self.nodesToTri = None
        self.nodeNormal = None
        self.nodeTris = None
        self.triNormal = None

        # other flags
        self.pointSearchInitialized = False

    def __str__(self):
        """Give a short description of the surface.

        >>> surf = SurfaceFromSTL("example.stl")
        >>> print surf
        TriangleSurface (SurfaceFromSTL) 214 nodes, 201 triangles, ...
        ... min: [0.0,0.4,1.2], max: [5.3,8.4,3.2]

        """
        return ("TriangleSurface (%s) %d nodes, %d triangles, %s"
                % (self.__class__.__name__, len(self.getNodeTris()),
                   len(self.triNodes), self.getBoundingBox()))


    def insertSurface(self, otherSurface, mergeNodes={}):
        """Add the contents of another surface to self.

        This has been taken from abq_model_02.Model.insertModel().

        Node and element number of the otherSurface will be conserved if there
        is no node resp. element in self falling in between the those new
        numbers. I.e. if there is node 1, 2, 10, 11 in self and nodes 3, 4, 5 in
        otherSurface then they will retain their numbers upon import to self and
        will be added to self as nodes 3, 4, 5. If there is again node 1, 2,
        10, 11 in self and now nodes with numbers 7, 12, 17 in otherSurface then
        they will get new numbers and become the new nodes 12, 13, 14 in self.
        (The same applies to element numbers accordingly.)

        Nodes that are not defined in otherSurface.nodeCoords but used in the
        element definitions otherSurface.triNodes will keep their old numbers.
        This is a feature to enable import of elements connected to already
        defined nodes with this function. You might want to check beforehand
        if all nodes are well defined if you don't need this feature.

        Node coordinates are deep copied, i.e. if you later change individual
        node coordinates in otherSurface then the node coordinates of self are
        not affected. Same with the element connectivity triNodes.

        @param otherSurface: Surface-instance to be added to self.

        @param mergeNodes: A dictionary {nodeOther : nodeThis}. nodeOther
          is the node number in otherSurface of the node to be replaced by
          the nodeThis-node in self.

          I.e. say mergeNodes={1000:10, ...}. Node 10 is a node in self. All
          occurences of node 1000 in otherSurface are replaced by 10 during the
          process of adding otherSurface to self. The node 1000 from
          otherSurface won't be imported to self.

        @Return: A tuple (nodesOldToNew, elemsOldToNew):
          - nodesOldToNew: {node number in otherSurface, new node num}
          - elemsOldToNew: {elNum in otherSurface, new elNum}
        """

        otherNodesToAdd = sorted( node for node in otherSurface.nodeCoords
                                  if node not in mergeNodes )

        # determine new node numbers
        if len(otherNodesToAdd)==0:
            newNodes = []
        elif set(self.nodeCoords).intersection(
            range(otherNodesToAdd[0], otherNodesToAdd[-1]+1)):
            # conflicting node number ranges
            firstNode = max(self.nodeCoords)+1
            newNodes = range(firstNode, firstNode+len(otherNodesToAdd))
        else:
            newNodes = otherNodesToAdd

        # insert nodes
        self.nodeCoords.update(
            (newN, otherSurface.nodeCoords[oldN][:])
            for newN, oldN in izip(newNodes, otherNodesToAdd)
            )

        # initialize nodesOldToNew: {node number in otherSurface, new node num}
        nodesOldToNew = dict(mergeNodes)
        nodesOldToNew.update(izip(otherNodesToAdd, newNodes))
        del otherNodesToAdd, newNodes

        # determine new element numbers
        oldElems = sorted(otherSurface.triNodes)
        if len(oldElems) and set(self.triNodes).intersection(
                range(oldElems[0], oldElems[-1]+1)):
            # conflicting element number ranges
            firstElem = max(self.triNodes)+1
            newElems = range(firstElem, firstElem+len(oldElems))
        else:
            newElems = oldElems

        # initialize elemsOldToNew: {elNum in otherSurface, new elNum}
        elemsOldToNew = dict( izip(oldElems, newElems) )
        del oldElems, newElems

        # insert elements
        triNodesDict = dict(
            ( elemsOldToNew[elem],
              [ nodesOldToNew.get(node, node) for node in nodes] )
            for elem, nodes in otherSurface.triNodes.iteritems())
        # update Surface object attributes
        dict.update(self.triNodes, triNodesDict)

        # The following attributes are reset to None.
        # They are supposed to be None as long as they are not calculated yet
        # don't query these variables, use the corresponding get...() method,
        # ie. for self.edgeToTri use self.getEdgeToTri()
        # This operation renders this data invalid / incomplete therefore
        # the corresponding variable (eg. self.edgeToTri) is reset to None.
        # It will then be recalculated automatically
        self.edgeToTri = None
        self.borderEdges = None
        self.borderNodes = None
        self.nodesToTri = None
        self.nodeNormal = None
        self.nodeTris = None
        self.triNormal = None
        self.pointSearchInitialized = False
        
        # return some of the transfer parameters
        return (nodesOldToNew, elemsOldToNew)

        

    def getAbqModel(self, elOffset=0, nodeOffset=0, elType="S3"):
        """create an abq_model from this surface
        
        node numbers start with nodeOffset+1, element numbers with elOffset+1
        node coords are copied, not referenced, if surfaces node coords change,
        those of the abq_model don't.
        """
        m_out = Model()
        all_nodes = dict()
        for el, nodes in self.triNodes.iteritems():
            elOffset += 1
            new_nodes = list()
            for node in nodes[:3]:
                try:
                    new_nodes.append(all_nodes[node])
                except KeyError:
                    nodeOffset += 1
                    all_nodes[node] = nodeOffset
                    new_nodes.append(nodeOffset)
            m_out.updateElem(elOffset, elType, new_nodes)
        for node_old, node_new in all_nodes.iteritems():
            m_out.nodeCoords[node_new] = list(self.nodeCoords[node_old])

        return m_out


    def asAbqModel(self, elType="S3"):
        """Returns self as an instance of abq_model_02.Model, so it can perform
        any actions specific to this class like reading and writing.

        Node numbers and element numbers are left unchanged.
        Node coords and element connectivity (i.e. Model.elNodes) are
        referenced, not copied. If surfaces node coords change, those of the
        abq_model do as well.

        This method is intended as a replacement for self.getAbqModel() when
        it's rather a conversion than a duplication that you are after. It's a
        good idea to delete the original TriangleSurface object after calling
        this function.

        For different element and node numbers use the function
        abq_model_02.Model.renumber().

        Example:
         >>> from bae.surface_02 import SurfaceFromSTL
         >>>
         >>> # read a surface and remove disconnected single tris
         >>> surf = SurfaceFromSTL("mysurface.stl")
         >>> toRemove = surf.findSingleTris()
         >>> for tri in toRemove:
         >>>     surf.removeTri(tri)
         >>>
         >>> # convert to an Abaqus model
         >>> model = surf.asAbqModel()
         >>> del surf
         >>>
         >>> # renumber nodes and elements and export as Abaqus input file
         >>> model.renumber(nodeStart=20000, elemStart=10000)
         >>> model.write("mysurface.inp")
        """
        m_out = Model()
        m_out.nodeCoords = abqModelInternal.NodesDict(self.nodeCoords)
        m_out.updateElems(self.triNodes, elType)

        return m_out


    def exportAsRaw(self, outputFile, objectName="Object1", numberFormat="%g"):
        """Export the triangle surface to the specified file in .raw format.
        This is suitable for import as mesh into Rhino.

        @param outputFile: If outputFile has got a write method (i.e. if
        outputFile is a file object or something similar) the output is being
        written to that file and the file is left open. Otherwise outputFile
        is treated as a file name and a corresponding file is being opened
        for write and closed after this procedure. If a file with that name
        existed before it is silently overwritten.
        @param objectName: Each mesh object in the raw file may be given a
        name that will become the mesh objects name in Rhino. This parameter
        accepts a string for this name.
        @param numberFormat: You may specify a format specifier that'll be used
        for each number (float) in the resulting raw file output. This can be
        used to increase the precision. E.g. "%1.3f" will print three decimal
        places.
        """
        try:
            writemethod = outputFile.write
        except AttributeError:
            output = open(outputFile, "w")
            writemethod = output.write

            
        writemethod("%s\n" % objectName)
        elems = self.triNodes.keys()
        elems.sort()
        for el in elems:
            nodes = self.triNodes[el]
            coords = [self.nodeCoords[node][ix]
                      for node in nodes
                      for ix in xrange(3)]
            line = " ".join(numberFormat%x for x in coords)
            writemethod("%s\n" % line)

        try:
            output.close()
        except NameError:
            pass
        return


    def exportAsSTL(self, fileName):
        """Export the triangle surface to the new file with the specified file
        name in .stl (stereolithography) format. This is suitable for import as
        mesh into Rhino or as geometry into icem.

        The binary format is exported. The normals are not exported.

        It might still be advisable to call self.getTriNormal() in advance.
        This will ensure that the vertices (nodes) are ordered in a consistent
        way. Not necessarily with normals pointing outwards as STL requires,
        though.

        @Note: Existing files with the same name are silently overwritten.
        """

        msg("Writing binary stl file %s" % fileName)

        # Write binary stl file:
        # Check wikipedia for the binary layout of the file.
        # Uses the struct library to convert into binary form.
        fp=open(fileName,'wb')
        
        header = struct.pack('80s', "binary stl")
        fp.write(header)

        length = struct.pack('I', len(self.triNodes))
        fp.write(length)

        emptynormal = struct.pack("<3f", 0.0, 0.0, 0.0)
        emptyattrib = struct.pack("2s", "")
        ptpacker = struct.Struct("<3f")

        for nodes in self.triNodes.itervalues():
            assert(len(nodes)==3)
            fp.write(emptynormal)
            for node in nodes:
                coords = self.nodeCoords[node]
                fp.write(ptpacker.pack(*coords))
            fp.write(emptyattrib)

        fp.close()
        msg("Finished writing binary stl file %s" % fileName)
        return


    def exportAsIv(self, fileName, colour=None, transparency=0):
        """Export the triangle surface to the new file with the specified file
        name in OpenInventor (.iv) format. This is suitable for import as shapes
        into voxler.
        """

        msg("Writing OpenInventor ascii file %s" % fileName)

        # find all connected nodes
        allnodes = set()
        for nodes in self.triNodes.itervalues():
            allnodes.update(nodes)
        allnodes = sorted(allnodes)

        # point list
        pointList = list()
        nodeToPtIdx = dict()
        for cnt, node in enumerate(allnodes):
            coords = self.nodeCoords[node]
            pointList.append(coords)
            nodeToPtIdx[node] = cnt

        # cell list
        triangles = [
            [nodeToPtIdx[node] for node in nodes]
            for nodes in self.triNodes.itervalues()
            ]

        output = open(fileName, "w")
        output.write("#Inventor V2.0 ascii\n")
        output.write("Separator {\n")

        # write material
        if colour is not None or transparency>0:
            if colour is None:
                colour = (1,1,1)
            else:
                colour = tuple(colour)
            output.write("  Material {\n")
            output.write("    ambientColor 0 0 0\n")
            output.write("    diffuseColor %g %g %g\n" % colour)
            output.write("    specularColor 0 0 0\n")
            output.write("    emissiveColor 0 0 0\n")
            output.write("    shininess 1\n")
            output.write("    transparency %g\n" % transparency)
            output.write("  }\n")
            output.write("  ShapeHints {\n")
            output.write("    vertexOrdering COUNTERCLOCKWISE\n")
            output.write("    shapeType UNKNOWN_SHAPE_TYPE\n")
            output.write("    faceType CONVEX\n")
            output.write("    creaseAngle 0\n")
            output.write("  }\n")

        # write point / node data
        output.write("  Coordinate3 {\n    point [\n")
        for coords in pointList:
            output.write("%g %g %g,\n" % tuple(coords))
        output.write("    ]\n  }\n")

        # write face / triangle data
        output.write("  IndexedFaceSet {\n    coordIndex [\n")
        for pts in triangles:
            output.write("%d,%d,%d,-1,\n" % tuple(pts))
        output.write("    ]\n  }\n")

        # fini
        output.write("}\n")
        output.close()
        msg("Finished writing OpenInventor file %s"% fileName)
        return

    def exportAsVtk(self, fileName, objectName="Object1", format='binary'):
        """Export the triangle surface to the specified file in lagacy .vtk
        format. This is suitable for paraview and others. Not for voxler.

        @param fileName: ... guess! 
        @param objectName: Name parameter in the vtk file. This parameter
        accepts a string.
        @param format: "binary" or "ascii"
        @Note: Existing files with that name are being overwritten silently.
        """
        # point list
        pointList = list()
        nodeToPtIdx = dict()
        for cnt, (node, coords) in enumerate(self.nodeCoords.iteritems()):
            pointList.append(coords)
            nodeToPtIdx[node] = cnt

        # cell list
        triangles = [
            [nodeToPtIdx[node] for node in nodes]
            for nodes in self.triNodes.itervalues()]

        # export to vtk
        vtkGrid = UnstructuredGrid( pointList, triangle=triangles )
        v = VtkData(vtkGrid, objectName)
        v.tofile(fileName, format)

        # fini
        msg("Finished writing vtk file %s"% fileName)

    # getRhinoMeshGeom --- only if used with Rhino
    if Rhino:
        def getRhinoMeshGeom(self, welded=False):
            """DEPRECATED: Use bert.lib.fearo_02.getModelFromMeshSurface
            instead!

            Return a Rhino.Geometry.Mesh object of self

            Example:
             >>> surf = TriangleSurface(....)
             >>> rhMesh = surf.getRhinoMeshGeom()
             >>> guid = scriptcontext.doc.Objects.AddMesh(rhMesh)
             >>> rs.ObjectName(guid, "new surface") # Rhino object name
             >>> scriptcontext.doc.Views.Redraw()

            @param welded: if True then create a mesh object with each vertex /
                node stored only once: a welded mesh. Otherwise each node /
                vertex is stored once for each facet / triangle / cohesive
                element. Node and element numbers will only be stored if this
                argument is false.
            """
            mesh = Rhino.Geometry.Mesh()
            # abbreviations (for speed)
            Point3d = Rhino.Geometry.Point3d
            Point2f = Rhino.Geometry.Point2f

            if welded:

                # find all connected nodes
                allnodes = set()
                for nodes in self.triNodes.itervalues():
                    allnodes.update(nodes)
                allnodes = sorted(allnodes)

                # vertices
                nodeToPtIdx = dict()
                for cnt, node in enumerate(allnodes):
                    mesh.Vertices.Add(Point3d(*self.nodeCoords[node]))
                    nodeToPtIdx[node] = cnt

                # faces
                faces = [ Rhino.Geometry.MeshFace(*[nodeToPtIdx[node]
                                                    for node in nodes])
                          for nodes in self.triNodes.itervalues() ]
                mesh.Faces.AddFaces(faces)

            else:

                nodeIdx = 0
                faces = list()
                tcs = System.Array.CreateInstance(Point2f, len(self.triNodes)*3)
                for elem, nodes in self.triNodes.iteritems():
                    nbnodes = len(nodes)
                    facenodes = range(nodeIdx, nodeIdx+nbnodes)

                    for i, node in izip(facenodes, nodes):
                        # vertices
                        coords = self.nodeCoords[node]
                        pt = Point3d(*coords)

                        mesh.Vertices.Add(pt)
                        # element and node number in texture coordinates
                        tcs[i] = Rhino.Geometry.Point2f(elem, node)

                    # add face to faces list
                    faces.append( Rhino.Geometry.MeshFace(*facenodes) )

                    nodeIdx += nbnodes
                mesh.Faces.AddFaces(faces)
                mesh.TextureCoordinates.SetTextureCoordinates(tcs)

            # result: mesh geom object
            return mesh


    def getEdgeToTri(self):
        """returns a dict {set(node1,node1):set of tri-element-ids}"""
        if self.edgeToTri == None:
            self.edgeToTri = defaultdict(set)
            for element, nodesOfThisTri in self.triNodes.iteritems():
                for thisEdge in self.triEdgesIter(nodesOfThisTri):
                    self.edgeToTri[thisEdge].add(element)
        return self.edgeToTri

    def getBorderEdges(self):
        """returns a set of edges ( =set(node1,node2) ) with only one tri
        connected to any other of this surface"""
        if self.borderEdges==None or self.edgeToTri==None:
            edgeToTri = self.getEdgeToTri()
            self.borderEdges = set([edge
                                    for edge,tris in edgeToTri.iteritems()
                                    if len(tris)==1])
        return self.borderEdges


    def getBorderNodes(self):
        """return the set of all nodes anywhere on the border"""
        if self.borderNodes==None or self.borderEdges==None:
            self.getBorderEdges()
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
        if self.nodesToTri == None:
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
        """returns a dict {node id: [tri ids]}

        also recommended if you need a list of all connected nodes:
        >>> allNodes = set(self.getNodeTris())
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


    def getTriNormal(self, skipMoebiusCheck=True):
        """returns a dict: {triId : normal vector} and reorders the triangle
        connectivity lists in self.triNodes according to the orientation.

        The normal vectors of the surface triangles will all point to the same
        side of the surface. They are normed, i.e. have length 1.

        Additionally the orientation of the tri element is saved by means of
        reordering its nodes in the self.triNodes-dictionary in
        counter-clockwise succession if you look from "above" (positive side).
        (The calculated normal points "up" to the positive side.)

        At an edge with more than two triangles attached, the algorithm stops.
        Then it starts again with an arbitrary triangle that has not yet got a
        normal. Until all normals are calculated for each triangle.

        Simply connected means there are max two triangles connected to each
        edge.

        if skipMoebiusCheck==True, no check is performed, if the surface
        has two sides and the normals are determinable uniquely. Otherwise a
        warning is logged "surface is like the Moebius strip".
        """
        if self.triNormal==None:

            self.triNormal=dict()
            self.moebius_flag = False
            self.moebiusStopEdges = set()
            if len(self.triNodes)==0: return

            # set of tris for which a normal has not been calculated yet
            noNormalTris = set(self.triNodes)

            while len(noNormalTris)>0:

                # get any one element as a starting point
                startTriId = iter(noNormalTris).next()

                # calculate normals starting at that point
                trisWithNormalList = self.calcNormalsOnPart(
                    startTriId, splitangle=None,
                    skipMoebiusCheck=skipMoebiusCheck)
                noNormalTris.difference_update(trisWithNormalList)

            if self.moebius_flag:
                msg("WARNING: This surface is twisted like the Moebius strip.")
            
        return self.triNormal


    def getNodeNormal(self):
        """
        Returns dict {node number -> averaged normal of the surface}
        
        The dictionary is stored internally, for the case that the function is
        called again.
        """

        if ((self.nodeNormal==None) or (self.triNormal==None)):
            # if not yet created do so, store in self.nodeNormal for second use

            triNormal = self.getTriNormal()
            nodeTris = self.getNodeTris()

            # self.nodeNormal
            self.nodeNormal = dict()
            for node in nodeTris.iterkeys():
                nodeNormal = [0,0,0]
                for element in nodeTris[node]:
                    nodeNormal = vector_plus(nodeNormal,triNormal[element])
                self.nodeNormal[node] = norm(nodeNormal)

        # return the dict anyway
        return self.nodeNormal


    def checkSimplyConnected(self, debugElset=None, debugNset=None):
        """
        check if the surface triangles are simply connected

        all edges must have exactly one or two triangles connected to it

        @param debugElset: if not None add diagnostic elsets to this dict
        {elset name: elset}.

        @param debugNset: if not None add diagnostic nsets to this dict
        {nset name: nset}.

        @returns: True if all triangles are simply connected, otherwise false
        """
        wrongTris = set()
        wrongNodes = set()
        for edge, tris in self.getEdgeToTri().iteritems():
            if len(tris) > 2:
                wrongTris.update(tris)
                wrongNodes.update(edge)

        if len(wrongTris):
            if debugElset!=None or debugNset!=None:
                msg("Found tris that are not 'simply connected'.")
            if debugElset!=None:
                warn_elset_name = "_BADCONNECT"
                debugElset[warn_elset_name] = wrongTris
                msg("All tris on an edge with more than two tris on it"
                    " (a 'surface junction') are collected in elset %s."
                    % warn_elset_name)
            if debugNset!=None:
                warn_nset_name = "_BADCONNECT"
                debugNset[warn_nset_name] = wrongNodes
                msg("All nodes on an edge with more than two tris on it"
                    " (a 'surface junction') are collected in nset %s."
                    % warn_nset_name)

        return len(wrongTris)==0

            
    def findSingleTris(self, debugElset=None):
        """Returns a list of tri elements without neighbours

        @param debugElset: if not None add diagnostic elsets to this dict
        {elset name: elset}.
        """
        borderEdges = self.getBorderEdges()

        singleTris = [
            triId for triId, nodes in self.triNodes.iteritems()
            if len(borderEdges.intersection(self.triEdgesIter(nodes))) == 3]

        if len(singleTris)>0 and debugElset!=None:
            warn_elset_name = "_SINGLE-TRI"
            debugElset[warn_elset_name] = set(singleTris)
            msg("Removed %d tris from this surface because they had three edges"
                " on the border of the surface, so they were not connected at"
                " all.\n"
                "Added them to the elset %s for debugging."
                % (len(singleTris), warn_elset_name))
        
        return singleTris


    def calcAngleToNeighbours(self, triNodes, default=None):
        """calculates the cosines of the angles between the normals of the tri
        given by triNodes and of all of its neighbours.

        There does not have to be a tri on the surface connecting triNodes,
        this is meant to check the angle before actually inserting a new tri.

        The sign of the angle is arbitrary, since node ordering is not taken
        into account.

        The result is always a list of three members. Where there is no
        neighbour, the value of the argument "default" is returned in that
        list.
        """
        edgeToTri = self.getEdgeToTri()

        n1 = norm(self.getTriNormalUsc(tuple(triNodes)))

        angles = [default, default, default]
        for edgeIdx, thisEdge in enumerate(self.triEdgesIter(triNodes)):
            nodes = set()
            for thisTri in edgeToTri[thisEdge]:
                nodes.update(self.triNodes[thisTri])
            nodes.difference_update(triNodes)
            if len(nodes)!=1: continue
            nodes.update(thisEdge)
            n2 = norm(self.getTriNormalUsc(tuple(nodes)))
            angles[edgeIdx] = dot(n1,n2)

        return angles


    def checkNoNormal(self, debugElset=None):
        """search for tris from surface that have not got a normal
        call this only after getTriNormal()
        returns a set of removed triangles (as they are in self.triNodes)
        ... or None if normals have not been calculated yet.

        @param debugElset: if not None add diagnostic elsets to this dict
        {elset name: elset}.
        """
        if self.triNormal == None:
            return None

        obscureTris = set(self.triNodes).difference(self.triNormal)

        if len(obscureTris) and debugElset!=None:
            warn_elset_name = "_IGNORED-TRI"
            debugElset[warn_elset_name] = obscureTris
            msg("The normals of %d tri elements could not be calculated due to"
                " poor connectivity.\n"
                "They have been removed from the surface and can be found in"
                " elset %s."
                % (len(obscureTris), warn_elset_name))

        return obscureTris


    def getBoundingBox(self):
        """
        Returns the bounding box of the surface

        Example:
          >>> # bounding box of that surface
          >>> bb = surf.getBoundingBox()
          >>> print "lower left front corner:", bb[0]
          >>> print "upper right back corner:", bb[1]

        @Returns: a BoundingBox object which is basically list of two
        coordinate tuples (rather lists).
        The first coord tuple states the min value, the last the max.
        I.e. self.getBoundingBox()[i][j]:
        i \in {0,1} for min/max value, j \in {0,1,2} coordinate index
        """

        # nodes set (contains unique nodes)
        nodes = set(self.getNodeTris())

        # update bounding box
        boundingBox = BoundingBox()
        boundingBox.update(nodeCoords=self.nodeCoords, elNodes=nodes)
        return boundingBox


    def removeTri(self, triId):
        """remove a tri from the surface

        update self.triNodes, self.edgeToTri, self.borderEdges
        remove self.borderNodes, self.nodesToTri
        """

        # update self.triNodes
        nodes = self.triNodes[triId]
        del self.triNodes[triId]

        # update self.edgeToTri
        if self.edgeToTri!=None:
            removedEdges = list()
            for thisEdge in self.triEdgesIter(nodes):
                self.edgeToTri[thisEdge].remove(triId)
                removedEdges.append(thisEdge)
                if len(self.edgeToTri[thisEdge]) == 0:
                    del self.edgeToTri[thisEdge]

        # update self.borderEdges
        if self.borderEdges!=None:
            if self.edgeToTri==None:
                self.borderEdges=None
            else:
                for thisEdge in removedEdges:
                    if len(self.edgeToTri[thisEdge]) == 0:
                        self.borderEdges.remove(thisEdge)
                    elif len(self.edgeToTri[thisEdge]) == 1:
                        self.borderEdges.add(thisEdge)

        # update self.nodeTris
        if self.nodeTris!=None:
            for node in nodes:  # all nodes of the tri to remove
                self.nodeTris[node].remove(triId)
                if len(self.nodeTris[node])==0:
                    del self.nodeTris[node]

        # update self.triNormal
        if self.triNormal != None:
            del self.triNormal[triId]

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

        # if normals have been computed we want to update self.triNormal
        # for that we need the neighbour and for that self.edgeToTri
        if self.triNormal != None:
            self.getEdgeToTri()

        # update edgeToTri and
        # find a neighbour to determine normal
        # ... only if edgeToTri already defined
        newBorderEdges = list()
        neighbour = None
        if self.edgeToTri != None:
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
                    if (neighbour==None and self.triNormal!=None):
                        newneighbour = tuple(newneighbour)[0]
                        if newneighbour in self.triNormal:
                            neighbour = newneighbour

                except KeyError:
                    # this edges is not yet in self.edgeToTri, its a new edge
                    # a new edge is an edge on the border
                    newBorderEdges.append(edge)
                    
                # anyway update self.edgeToTri
                self.edgeToTri[edge].add(newTri)

        # update borderEdges
        if self.borderEdges != None:
            if self.edgeToTri == None:
                self.borderEdges = None
            else:
                for edge in self.triEdgesIter(nodes):
                    try:
                        # a border edge getting a new tri attached is no longer
                        # on the border if it has been on the border before
                        self.borderEdges.remove(edge)
                    except ValueError:
                        pass
                self.borderEdges.update(newBorderEdges)

        # must have at least one neighbour
        # at least when triNormal is defined
        if self.triNormal != None:
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
            self.triNormal[newTri] = self.getTriNormalSc(newNodes)
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
        self.nodeNormal = None


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
            raise ValueError("Node %d not in tri element %d in this surface"
                             " when trying to replace it."
                             % (oldNode, triId))
        except IndexError:
            raise IndexError("No tri element %d in this surface when trying"
                             " to replace node %d on it."
                             % (triId, oldNode))
        # update triNodes
        self.triNodes[triId][idx] = newNode

        # update edgeToTri
        if self.edgeToTri != None:
            borderEdgesAdd = list()
            borderEdgesRemove = list()
            for edge,triList in self.edgeToTri.items():
                if (oldNode in edge) and (triId in triList):
                    triList.remove(triId)
                    if len(triList)==1:
                        borderEdgesAdd.append(edge)
                    if len(triList)==0:
                        del self.edgeToTri[edge]
                        borderEdgesRemove.append(edge)
                    newEdge = set(edge)
                    newEdge.remove(oldNode)
                    newEdge.add(newNode)
                    newEdge = frozenset(newEdge)
                    triList = self.edgeToTri[newEdge]
                    triList.add(triId)
                    if len(triList)==1:
                        borderEdgesAdd.append(newEdge)
                    if len(triList)==2:
                        borderEdgesRemove.append(newEdge)

        # update borderEdges
        if self.borderEdges != None:
            if self.edgeToTri == None:
                self.borderEdges = None
            else:
                self.borderEdges.update(borderEdgesAdd)
                self.borderEdges.difference_update(borderEdgesRemove)

        # update self.triNormal
        if self.triNormal != None:
            self.triNormal[triId] = self.getTriNormalSc(self.triNodes[triId])

        # delete obscured members
        self.borderNodes = None
        self.nodesToTri = None
        self.nodeNormal = None
        self.nodeTris=None

    #-- is point below surface?

    cellSizeFactor=0.5      # cell Size for array relative to 'avgSize'
    # in method pointIsBelow() tolerance is compared to values 0..1
    tolerance = 1E-4


    def translate(self, move_vector):
        """Move the surface in space by move_vector

        Example:

        >>> # move this surface up by 10
        >>> surf.translate(move_vector=[0,0,10])

        Translating is done by simply moving the point coordinates for this
        surface. Other surfaces or the original mesh from which the surface
        originates are not affected.

        This creates a new node-dictionary, which is refered to from now on.
        It contains only those nodes connected to any element in self.triNodes.

        The point search initialization is modified and doesn't have to be done
        again.
        """

        all_nodes = set(self.getNodeTris())
        
        new_coords = dict()
        for node_id in all_nodes:
            new_coords[node_id] = vector_plus(
                self.nodeCoords[node_id],move_vector)
        
        self.nodeCoords = new_coords

        # adapt point search properties
        if self.pointSearchInitialized and len(self.triNodes)>0:
            
            # tet Centroids
            for centroid in self.triCentroids.itervalues():
                centroid[0] += move_vector[0]
                centroid[1] += move_vector[1]

            # BoundingBox
            self.boundingBox[0][0] += move_vector[0]
            self.boundingBox[0][1] += move_vector[1]
            self.boundingBox[1][0] += move_vector[0]
            self.boundingBox[1][1] += move_vector[1]


    def scale(self, scale_factor, scale_origin):
        """Scale the surface by a certain factor relative to scale_origin

        scale_factor may be a single number in which case it is applied to all
        three dimensions or a list of three numbers with one factor for each
        dimension.

        Example:

        >>> # scale this surface by 2 in each direction
        >>> surf.scale(scale_factor=2, scale_origin=[0,0,0])

        Scaling is done by simply moving the point coordinates for this
        surface. Other surfaces or the original mesh from which the surface
        originates are not affected.

        This creates a new node-dictionary, which is refered to from now on.
        It contains only those nodes connected to any element in self.triNodes.

        After scaling point search has to be initialized again.
        """

        all_nodes = set(self.getNodeTris())

        if isinstance(scale_factor, (float, int)):
            scale_factor = [scale_factor,scale_factor,scale_factor]
        
        new_coords = dict()
        for node_id in all_nodes:
            orig2old = vector(scale_origin, self.nodeCoords[node_id])
            orig2new = [orig2old[0]*scale_factor[0],
                        orig2old[1]*scale_factor[1],
                        orig2old[2]*scale_factor[2]]
            new_coords[node_id] = vector_plus(scale_origin, orig2new)
        
        self.nodeCoords = new_coords

        # after scaling point search has to be initialized again
        self.pointSearchInitialized = False


    def initializePointSearch(self, smooth_findradius=0.0):
        """
        Usually you don't need to call this function.

        You may call this function once before pointIsBelow or projectPoints,
        otherwise it is done automatically.

        smoothmode: Additionally consider points to be below the surface that
        are in fact outside but within the smooth_findradius to the
        closest centre point of any of the triangles of the surface (in the
        x-y-plane).
        """

        # list of centroid-coord-tuples
        self.triCentroids=dict()

        # list of [[xmin,ymin],[xmax,ymax]] for each tri element
        elementBox=dict()

        boundingBox = BoundingBox(dim=2)
        avgSize=0

        # size of the projection on the x-z plane
        # pos facing one direction, neg facing the other
        elSignedSize = dict()
        for element,nodes in self.triNodes.iteritems():

            node_pts = [self.nodeCoords[node]
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
            try:
                for i in xrange(2): # x,y,z
                    if elementBox[element][0][i]<boundingBox[0][i]: # min
                        boundingBox[0][i]=elementBox[element][0][i]
                    if elementBox[element][1][i]>boundingBox[1][i]: # max
                        boundingBox[1][i]=elementBox[element][1][i]
            except IndexError:
                boundingBox[0]=elementBox[element][0][:2]
                boundingBox[1]=elementBox[element][1][:2]

            avgSize+=self.dist2D(
                elementBox[element][1],elementBox[element][0])

        avgSize /= len(self.triNodes)

        msg('bounding box of the surface: %s' % boundingBox)
        msg('avg. element size: %g' % avgSize)


        ## -----------------------------------------------------------------
        ## map elements to cells

        cellSize=avgSize*self.cellSizeFactor

        # x,y,z number of cells
        cellDimensions = [boundingBox[1][i]-boundingBox[0][i] for i in range(2)]
        msg('cell array dimensions: %s' % " ".join(map(str, cellDimensions)))

        cellNumber = [int(x/cellSize)+1 for x in cellDimensions]
        msg('cell numbers: %s cells: %d'
            % (" ".join(map(str, cellNumber)), cellNumber[0]*cellNumber[1]))

        # create empty array
        for i in range(cellNumber[0]):
            if i==0: cells=[]
            for j in range(cellNumber[1]):
                if j==0: LJ=[]
                LJ.append([])
            cells.append(LJ)

        # tris with a too small projected size will be ignored
        zerosize = (self.tolerance*avgSize*0.1)**2

        # map elements
        for element in self.triNodes:

            # ignore too small tris
            if abs(elSignedSize[element]) < zerosize:
                continue

            # get index ranges
            cellIndex=[[0,0],[0,0]] # low,high
            cellOrigin = boundingBox[0]
            for l in range(2): # low,high - index
                for i in range(2): # dimension (x,y) - index
                    cellIndex[l][i]=int((elementBox[element][l][i]
                                         -cellOrigin[i])/cellSize)
            # append element ID to cells that overlap the element's bounding box
            for i in range(cellIndex[0][0],cellIndex[1][0]+1):
                for j in range(cellIndex[0][1],cellIndex[1][1]+1):
                    cells[i][j].append(element)

        # permanently store objects as object attributes
        self.boundingBox = boundingBox
        self.cellNumber = cellNumber
        self.cellSize = cellSize
        self.cells = cells

        # for smooth mode store radius
        self.smoothmode = smooth_findradius>0
        self.smooth_findradius = float(smooth_findradius)

        # finished initializing
        self.pointSearchInitialized = True
        return

    def pointIsBelow(self, point, returnElems=False):
        """
        Test if point is exactly below or above the surface in the sense of
        a vertical projection on the x-y plane. The function does not decide
        whether a point is actually below or above the surface. The z
        coordinates are ignored completely.

        @param point: a tuple of coordinates (of which the z coordinate is
        neglected)

        @param returnElems: If True the function returns a dictionary
        {triNb:(a0, a1)}.
        The vertical projection of the point on the tri element triNb is
        pt_0 + a0*vec_01 + a1*vec_02, with pt_0 being
        self.nodeCoords[self.triNodes[triNb][0]], vec_01 being
        the vector from node self.triNodes[triNb][0] to self.triNodes[triNb][1]
        and vec_02 from self.triNodes[triNb][0] to self.triNodes[triNb][2]:
        
        >>> from bae.vecmath_01 import *
        >>> point = [...some coords...]
        >>> surface = TriangleSurface(....)
        >>> triA0A1 = surface.pointIsBelow(point, returnElems=True)
        ... for triNb, (a0,a1) in triA0A1.iteritems():
        ...     coords = [surface.nodeCoords[node]
        ...               for node in surface.triNodes[triNb]]:
        ...                  # newPt = pt_0 + a0*vec_01 + a1*vec_02
        ...                   newPt = list(coords[0])
        ...                   vector_modif_add(newPt, vector_scale(vec_01,a0))
        ...                   vector_modif_add(newPt, vector_scale(vec_02,a1))

        Elements that are only found by means of the smoothmode (see
        self.initializePointSearch) will have (a0, a1) = (0,0).

        @Note: On the first invocation of pointIsBelow() the search algorithm
        is initialized by automatically calling self.initializePointSearch().
        If you now change the coordinates or the connectivity of the surface
        subsequent calls to pointIsBelow() will not work properly. To force a
        new initialization based on the updated data either call
        self.initializePointSearch() manually or set
        self.pointSearchInitialized = False.

        @Bug: tris with a too small projected size will be ignored
        """
        if returnElems:
            notFound = dict()
        else:
            notFound = False

        if len(self.triNodes)==0:
            return notFound

        if not self.pointSearchInitialized:
            self.initializePointSearch()

        # cell check (includes bounding box test from cell index)
        index=[0,0]
        cellOrigin = self.boundingBox[0]
        for i in range(2):
            index[i]=int((point[i]-cellOrigin[i])/self.cellSize)
            if index[i]<0 or index[i]>=self.cellNumber[i]:
                return notFound

        # if cell check passed ... go on
        elems_in_cell=self.cells[index[0]][index[1]]
        centroid_to_point_distance=[]
        # sort by distance
        for element in elems_in_cell:
            l=self.dist2D(point,self.triCentroids[element])
            centroid_to_point_distance.append([l,element])
        centroid_to_point_distance.sort()

        # full check if point inside any tri of elems_in_cell
        inside = dict()

        if (self.smoothmode and centroid_to_point_distance):
            for l, element in centroid_to_point_distance:
                if l>=self.smooth_findradius:
                    break
                inside[element] = (0,0)
                if not returnElems:
                    break

        for [l,element] in centroid_to_point_distance:
            node_pts = [self.nodeCoords[node]
                        for node in self.triNodes[element]]
# Explanation of the following algorithm:

#       2
#      ^
#     /  .
#  b /     .
#   /        .
# 0 ----------> 1
#        a   

# How to find out, if a point is inside the triangle 0 - 1 - 2:
# or, more precise, if the projection of that point is inside 
# the triangle:

# Vectors starting on 0 end ending on the dotted line are:
#  v. = delta a + (1-delta) b              with 0 < delta < 1
# Vectors shorter than v. but pointing in the same direction are:
#  v = gamma v.                          with 0 < gamma < 1
#    = gamma ( delta a + (1 - delta) b)
#    = gamma delta a + gamma ( 1 - delta) b
#
# Define: 
# alpha := gamma delta                     with 0 < alpha < 1
#  beta := gamma (1 - delta)               with 0 < beta  < 1
#
# So that we have:
#  v = alpha a + beta b

# Looking at which values alpha and beta can take:
# alpha + beta = gamma (delta + ( 1 - delta) ) = gamma < 1
# <=> alpha + beta < 1   =>  alpha < 1 AND beta < 1

# Algorithm:
# Given a point p, which could be inside the triangle 0 - 1 - 2, by
# its position vector x_p. Find alpha, beta by:

# 1) Compute v = x_p - x_0, a vector pointing from 0 to the point.
# 2) Compute alpha and beta from the first 2 equations of:
#    v = alpha a + beta b
# 3) Check if alpha + beta < 1 and if alpha > 0 and if beta > 0.

            # vectors
            v12=vector(node_pts[0], node_pts[1])
            v13=vector(node_pts[0], node_pts[2])
            detA = float(v12[0]*v13[1] - v12[1]*v13[0])
            adjA = [[v13[1], -v13[0]], [-v12[1], v12[0]]]
            invA = [[x/detA for x in row] for row in adjA]

            v1p=vector(node_pts[0], point)
            a0 = invA[0][0]*v1p[0] + invA[0][1]*v1p[1]
            a1 = invA[1][0]*v1p[0] + invA[1][1]*v1p[1]

            # if a0>=0 and a1>=0 and (a0+a1)<=1:
            #     inside = 1
            #     break
            if ((a0 > -self.tolerance) and
                (a1 > -self.tolerance) and
                ((a0+a1) < (1+self.tolerance))):
                inside[element] = (a0, a1)
                if not returnElems:
                    break

        if returnElems:
            return inside
        else:
            return bool(len(inside))

    def projectPoints(self, points):
        """projects arbitrary points vertically on the surface

        Returns a list of as many items as there are items in points. Each of
        those items is a list of projected points. It may be empty if the point
        is not below or above the surface or it may contain one ore many points
        on the surface.

        If two of the points are closer than self.tolerance only the first is
        recognized.

        Note: Uses self.pointIsBelow(). See its documentation in case the
        coordinates of the surface change between calls.
        """
        result = list()
        for point in points:
            triA0A1 = self.pointIsBelow(point, returnElems=True)
            projectedPts = list()
            for triNb, (a0,a1) in triA0A1.iteritems():
                coords = [self.nodeCoords[node]
                          for node in self.triNodes[triNb]]
                vec_01 = vector(coords[0], coords[1])
                vec_02 = vector(coords[0], coords[2])                

                # newPt = pt_0 + a0*vec_01 + a1*vec_02
                newPt = list(coords[0])
                vector_modif_add(newPt, vector_scale(vec_01,a0))
                vector_modif_add(newPt, vector_scale(vec_02,a1))

                tolDist = self.tolerance*(
                    0.5*(length(vec_01)+length(vec_02)))
                for oldPt in projectedPts:
                    if dist(newPt, oldPt)<tolDist:
                        newPt = None
                        break
                if newPt:
                   projectedPts.append(newPt)
            result.append(projectedPts)
        return result


    #-- some service functions: 2D distance

    @staticmethod
    def dist2D(a,b):
        "2D distance only considering x and y coordinates"
        return sqrt((b[0]-a[0])**2 + (b[1]-a[1])**2)


    #-- functions not meant to be called from "outside"

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
        """ Returns an iterator of the edges of the given triangle

        nodesOfThisTri will be converted to a tuple of node ids
        
        >>> myTri = (3, 43, 12)
        >>> for edges in triEdgesIter(myTri):
        >>> # edges will be:
        >>> # frozenset((3, 43)), frozenset((43, 12)), frozenset((12, 3))

        @note: See triEdgesOrderedNodesIter if you also need the nodes in
        the correct order.

        """
        nodesOfThisTri = tuple(nodesOfThisTri)
        for i1, i2 in ((0,1), (1,2), (2,0)):
            yield frozenset((nodesOfThisTri[i1], nodesOfThisTri[i2]))

    def getTriNormalUsc(self, nodes):
        """Not meant to be called from "outside"

        Computes the unscaled ("Usc") normal of a triangle. I.e. the length of
        the normal vector is arbitrary.

        Nodes come in counter-clockwise succession if you look from "above"
        (positive side), against the direction of the normal.
        """
        vectorAB=vector(self.nodeCoords[nodes[0]],self.nodeCoords[nodes[1]])
        vectorAC=vector(self.nodeCoords[nodes[0]],self.nodeCoords[nodes[2]])
        return cross(vectorAB,vectorAC)

    def getTriNormalSc(self, nodes):
        """Not meant to be called from "outside"

        Computes the scaled ("Sc") normal of a triangle. I.e. this normal
        vector is of length 1.

        Nodes come in counter-clockwise succession if you look from "above"
        (positive side), against the direction of the normal.
        """
        vectorAB=vector(self.nodeCoords[nodes[0]],self.nodeCoords[nodes[1]])
        vectorAC=vector(self.nodeCoords[nodes[0]],self.nodeCoords[nodes[2]])
        return norm(cross(vectorAB,vectorAC))

    def calcNormalsOnPart(self, startTriId,
                          splitangle=None, skipMoebiusCheck=False):
        """Not meant to be called from "outside"

        Make sure len(self.triNodes)>0, and self.moebiusStopEdges is a set.
        self.triNormal must already be a dict!

        creates object-attribute triNormal:
        a dict {element:[comp1, comp2, comp3]} with [comp...] being the vector
        components of the *scaled* normal of that triangle. I.e. each normal
        vector has length one.
        
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
        (newTri in self.triNormal) the algorithm stops as well.
        If the angle in degree between two neighbouring triangles is greater
        than *splitangle*, the algorithm stops at that edge. splitangle==None
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

        if splitangle != None:
            splitanglecos = cos(float(splitangle)*pi/180.)

        edgeToTri = self.getEdgeToTri()
        assert len(self.triNodes)>0

        trisWithNormalList = list()
        startTriNodes = self.triNodes[startTriId]
        startNormal = self.getTriNormalSc(startTriNodes)
        self.triNormal[startTriId] = startNormal
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

                trisOnThisEdge = edgeToTri[thisEdge]
                if len(trisOnThisEdge)!=2:
                    # more than two tris on this edge: ignore this edge
                    # less than two: edge on border, no other tri connected
                    continue

                # find the other tri connected to this edge (not startTriId)
                for newTri in trisOnThisEdge:
                    if newTri != startTriId: break


                # if no check and normal already calculated, skip the rest
                newTriHasNormalAlready = (newTri in self.triNormal)
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
                if (splitangle != None and
                    dot(startNormal, newNormal)<splitanglecos):
                    continue

                # if normal already calculated, compare and goto next
                if newTriHasNormalAlready:
                    if dot(newNormal, self.triNormal[newTri]) < 0.0:
                        self.moebius_flag = True
                        self.moebiusStopEdges.add(thisEdge)
                        msg("WARNING: The normals flip between element %d and"
                            " %d." % (startTriId, newTri))

                    # in any case keep the already calculated normal and
                    # proceed with the next tri
                    continue

                # update self.triNormal with new normal
                self.triNormal[newTri] = newNormal
                trisWithNormalList.append(newTri)
                # update new node order to self.triNodes
                self.triNodes[newTri] = newTriNodes
                # append this neighbour to trisToFindNeighbours-list
                trisToFindNeighbours.append((newTri, newTriNodes, newNormal))

        # fini
        return trisWithNormalList


## ---------------------------------------------------------------------------

class SurfaceFromMeshTris(TriangleSurface):
    """a surface created from some tri elements of an abaqus model
    """

    def __init__(self, model, elset=None, logfile=None):
        """
        All kinds of triangle elements in the model (or in the given elsets)
        are considered, non triangle elements are ignored.

        Only those nodes connected to the surface are copied to the
        self.nodeCoords attribute. But this is not a deep copy: The coordinate
        lists referenced by the dict self.nodeCoords are the same as those
        referenced by the model.nodeCoords dict.

        The node connectivity from model.elNodes however is deep copied: The
        node lists referenced by the dict self.triNodes are others then those
        referenced by model.elNodes.

        Elements in elset that are not found in model.elNodes are ignored, a
        warning is posted.

        @param model: an abq_model_02.Model object or a filename (a string) or
        an open file of an abaqus input file.

        @param elset: The surface shall consist of the elements in this elset
        might be anything that model.getUnionSet() accepts as input:
        an elset name, a set of element numbers, a list of elset names
        and element numbers
        elset might also be None, then all tri elements of model are used.

        @param logfile: deprecated argument, only for compatibility reasons.
        No effect.
        """

        TriangleSurface.__init__(self)

        if isinstance(model, (basestring, file)):
            # basestring = any kind of string
            m = Model()
            m.read(model)
            model = m
        else:
            checkModelModuleVersion(model, "%s.%s.__init__"
                                    % (__name__, self.__class__.__name__))

        alltris = set()
        for typ in model.elShapeToTypes['TRI'].intersection(model.typeEl):
            alltris.update(model.typeEl[typ])

        if elset==None:
            elset = alltris.intersection(model.elNodes)
        else:
            elset = model.getUnionSet("elset", elset)
            elset.intersection_update(alltris)

        nb_trisnotfound = 0
        nodestocopy = set()
        for elem in elset:
            try:
                nodes = model.elNodes[elem][:3]
                self.triNodes[elem] = nodes
                nodestocopy.update(nodes)
            except KeyError:
                nb_trisnotfound += 1
                pass
        if nb_trisnotfound>0:
            msg("WARNING from SurfaceFromMeshTris.__init__():"
                " Could not find element connectivity for %d"
                " of the %d tri elements in the specified"
                " elset." % (nb_trisnotfound, len(elset)))
        
        # copy all connected nodes (no deep copy)
        for node in nodestocopy:
            try:
                coords = model.nodeCoords[node]
                self.nodeCoords[node] = coords
            except KeyError:
                pass  

        if len(self.triNodes)==0:
            msg('WARNING: SurfaceFromMeshTris initialized with'
                ' empty element set.')
        return


## ---------------------------------------------------------------------------

class SurfaceFromAbqModelSurface(TriangleSurface):
    """A TriangleSurface created from the surface named surfName in the
    ABAQUS model abqModel.

    The surface must consist of faces of tet elements (C3D4 or C3D10M).
    Only the corner nodes of the tets are considered.

    The elsets abqModel+"_S1", abqModel+"_S2", ... "+_S4" define the
    surface.
    
    The node coords of this surface are a reference (no copy) to those of
    the abaqus model. If the abaqus model moves, so does this surface.
    """

    def __init__(self, abqModel, surfName, logfile=None):
        """create a TriangleSurface from the surface named surfName in the
        ABAQUS model abqModel.

        The surface must consist of faces of tet elements (C3D4 or C3D10M).
        Only the corner nodes of the tets are considered.

        The elsets abqModel+"_S1", abqModel+"_S2", ... "+_S4" define the
        surface.
        
        The node coords of this surface are a reference (no copy) to those of
        the abaqus model. If the abaqus model moves, so does this surface.

        @param logfile: deprecated argument, only for compatibility reasons.
        No effect.
        """

        TriangleSurface.__init__(self)
        checkModelModuleVersion(abqModel, "%s.%s.__init__"
                                % (__name__, self.__class__.__name__))

        triElementFaces={'S1':[0,2,1],'S2':[0,1,3],'S3':[1,2,3],'S4':[0,3,2]}

        nodestocopy = set()
        elnum = 1
        for faceType in set(abqModel.surface[surfName]).intersection(
            triElementFaces):
            nodeIds = triElementFaces[faceType]
            for el in abqModel.elset[
                abqModel.surface[surfName][faceType]]:
                if abqModel.elType[el] not in ('C3D4', 'C3D10M'):
                    raise NotImplementedError(
                        "%s.%s.__init__() can not create a surface from"
                        " elements of type %s."
                        % (__name__, self.__class__.__name__,
                           abqModel.elType[el]))
                nodes = [abqModel.elNodes[el][i] for i in nodeIds]
                self.triNodes[elnum] = nodes
                nodestocopy.update(nodes)
                elnum += 1
                
        for node in nodestocopy:
            try:
                coords = abqModel.nodeCoords[node]
                self.nodeCoords[node] = coords
            except KeyError:
                pass  

        if len(self.triNodes)==0:
            msg('WARNING: SurfaceFromMeshTris initialized with'
                ' empty element set.')
        return


## ---------------------------------------------------------------------------

class SurfaceFromDXF(TriangleSurface):
    """
    A surface created from some tri elements read from a DXF file

    Derived from dxv2inp.py

    @Note: This does not read DXF files created by rhino (default settings).
    Don't know why.
    """

    def __init__(self, filename, logfile=None, faceID = '3DFACE'):
        """
        @param filename: File name of a DXF file to read surface data from.

        @param faceID: Face identifier in the DXF file.

        @param logfile: deprecated argument, only for compatibility reasons.
        No effect.
        """

        TriangleSurface.__init__(self)

        inputFile = open(filename)
        lineNb = 0

        nodes = []
        elements = []
        while 1:
            line = inputFile.readline().strip()
            lineNb += 1
            if not line: break

            if line.split()[0]==faceID:

                line = inputFile.readline() # 2 lines junk
                line = inputFile.readline()
                lineNb += 2

                points = [0,0,0]
                for i in range(4):

                    coords = [0,0,0]
                    for j in range(3):

                        line = inputFile.readline() # point/xyz counter
                        line = inputFile.readline() # coord value
                        lineNb += 2
                        coords[j] = float(line)

                    try:
                        nodesIndex = nodes.index(coords) # node ID
                    except ValueError:
                        nodesIndex = len(nodes)
                        nodes.append(coords)

                    try:
                        if nodesIndex not in points:
                            points[i]=nodesIndex # element connectivity
                    except IndexError:
                        msg("WARNING: Found more than three different points"
                            " for triangle %d on line %d of the DXF input"
                            " file. Ignoring this point."
                            % (len(elements)+1, lineNb))

                elements.append(points)

        inputFile.close()

        self.nodeCoords.update((cnt+1, coords)
                               for cnt, coords in enumerate(nodes))
        self.triNodes.update((cnt+1, [node+1 for node in nodes])
                              for cnt, nodes in enumerate(elements))
        return


## ---------------------------------------------------------------------------

class SurfaceFromRaw(TriangleSurface):
    """
    A surface created from some tri elements read from a .raw file like Rhino
    creates them.

    Derived from raw2inp.py
    """

    def __init__(self, filename, logfile=None):
        """
        @param filename: File name of a .raw file to read surface data from.

        @param logfile: deprecated argument, only for compatibility reasons.
        No effect.
        """

        TriangleSurface.__init__(self)

        inputFile = open(filename)
        lineNb = 0

        nodes = []
        elements = []

        # ignore first line
        line = inputFile.readline().strip()
        msg("Reading an object called %s from %s"
            % (line, filename))

        while 1:
            line = inputFile.readline().strip()
            lineNb += 1
            if not line: break

            try:
                raw=map(float, line.split())
            except ValueError:
                # if something is not convertible to a float then
                # it must be another object, ignore this line
                msg("Reading an object called %s from %s"
                    % (line, filename))
                continue

            points = [0,0,0]
            for i in range(3):

                coords = raw[i*3:(i+1)*3]

                try:
                    nodesIndex = nodes.index(coords) # node ID
                except ValueError:
                    nodesIndex = len(nodes)
                    nodes.append(coords)

                points[i]=nodesIndex # element connectivity

            elements.append(points)

        inputFile.close()

        self.nodeCoords.update((cnt+1, coords)
                               for cnt, coords in enumerate(nodes))
        self.triNodes.update((cnt+1, [node+1 for node in nodes])
                              for cnt, nodes in enumerate(elements))
        return


## ---------------------------------------------------------------------------

class SurfaceFromSTL(TriangleSurface):
    """
    A surface created from some tri elements read from a .stl (stereo
    lithography or what) file.

    Derived from raw2inp.py and something from the internet
    """

    def __init__(self, filename, logfile=None):
        """
        @param filename: File name of a .raw file to read surface data from.

        @param logfile: deprecated argument, only for compatibility reasons.
        No effect.
        """

        TriangleSurface.__init__(self)

        coordsToNodes = dict()
        nodes = []
        elements = []

        # read start of file to determine if its a binay stl file or an ascii
        # stl file
        fp=open(filename,'rb')
        h=fp.read(80)
        type=h[0:5]
        fp.close()

        if type=='solid':
            msg("reading text file %s" % filename)

            # read text stl match keywords to grab the points to build the model
            fp=open(filename,'r')

            for line in fp:
                words=line.split()
                if len(words)>0:
                    if words[0]=='solid':
                        msg("Reading an object called %s from %s"
                            % (words[1], filename))
                    if words[0]=='facet':
                        triangle=[]
                        # ignoring the normal information
                        # normal=(eval(words[2]),eval(words[3]),eval(words[4]))

                    if words[0]=='vertex':
                        coords = tuple(float(x) for x in words[1:4])
                        try:
                            node = coordsToNodes[coords]
                        except KeyError:
                            node = len(nodes)
                            nodes.append(coords)
                            coordsToNodes[coords] = node
                        triangle.append(node)

                    if words[0]=='endloop':
                        # make sure we got the correct number of values before
                        # storing
                        if len(triangle)==3:
                            elements.append(triangle)
                        else:
                            msg("WARNING: found facet with %d vertices in %s"
                                % (len(triangle), filename))
            fp.close()

        else:
            msg("reading binary stl file %s" % filename)

            # Load binary stl file:
            # Check wikipedia for the binary layout of the file.
            # We use the struct library to read in and convert binary data
            # into a format we can use.
            fp=open(filename,'rb')
            h=fp.read(80)

            l=struct.unpack('I',fp.read(4))[0]
            def parse_pt(bytes):
                if len(bytes)==12:
                    return tuple(struct.unpack('f',bytes[slice(i, i+4)])[0]
                            for i in (0, 4, 8))
                else:
                    raise EOFError("Read less than 12 bytes for a coordinate"
                                   " tripel in %s." % filename)

            count=0
            while True:
                try:
                    p=fp.read(12)
                    n=parse_pt(p)

                    triangle=[]
                    for i in range(3):
                        p=fp.read(12)
                        coords = parse_pt(p)
                        try:
                            node = coordsToNodes[coords]
                        except KeyError:
                            node = len(nodes)
                            nodes.append(coords)
                            coordsToNodes[coords] = node
                        triangle.append(node)

                    elements.append(triangle)
                    count+=1
                    fp.read(2)

                except EOFError:
                    break
            fp.close()

        # copy data into triangle object
        self.nodeCoords.update((cnt+1, list(coords))
                               for cnt, coords in enumerate(nodes))
        self.triNodes.update((cnt+1, [node+1 for node in nodes])
                              for cnt, nodes in enumerate(elements))
        return

## ---------------------------------------------------------------------------

class SurfaceFromCohesives(TriangleSurface):
    """a surface created from some cohesive elements COH3D6 of an abaqus model
    """

    def __init__(self, model, elset=None, elType="COH3D6", logfile=None):
        """
        All kinds of triangle elements in the model (or in the given elsets)
        are considered, non triangle elements are ignored.

        Only those nodes connected to the surface are copied to the
        self.nodeCoords attribute. But this is not a deep copy: The coordinate
        lists referenced by the dict self.nodeCoords are the same as those
        referenced by the model.nodeCoords dict.

        The node connectivity from model.elNodes however is deep copied: The
        node lists referenced by the dict self.triNodes are others then those
        referenced by model.elNodes.

        Elements in elset that are not found in model.elNodes are ignored, a
        warning is posted.

        @param model: an abq_model_02.Model object or a filename (a string) or
        an open file of an abaqus input file.

        @param elset: The surface shall consist of the elements in this elset
        might be anything that model.getUnionSet() accepts as input:
        an elset name, a set of element numbers, a list of elset names
        and element numbers
        elset might also be None, then all COH3D6 elements of model are used.

        @param logfile: deprecated argument, only for compatibility reasons.
        No effect.
        """

        TriangleSurface.__init__(self)

        if isinstance(model, (basestring, file)):
            # basestring = any kind of string
            m = Model()
            m.read(model)
            model = m
        else:
            checkModelModuleVersion(model, "%s.%s.__init__"
                                    % (__name__, self.__class__.__name__))

        if elset==None:
            elset = model.typeEl[elType]
        else:
            elset = model.getUnionSet("elset", elset)
            elset.intersection_update(model.typeEl[elType])

        nb_trisnotfound = 0
        nodestocopy = set()
        for elem in elset:
            try:
                nodes = model.elNodes[elem][:3]
                self.triNodes[elem] = nodes
                nodestocopy.update(nodes)
            except KeyError:
                nb_trisnotfound += 1
                pass
        if nb_trisnotfound>0:
            msg("WARNING from SurfaceFromCohesives.__init__():"
                " Could not find element connectivity for %d"
                " of the %d tri elements in the specified"
                " elset." % (nb_trisnotfound, len(elset)))
        
        # copy all connected nodes (no deep copy)
        for node in nodestocopy:
            try:
                coords = model.nodeCoords[node]
                self.nodeCoords[node] = coords
            except KeyError:
                pass  

        if len(self.triNodes)==0:
            msg('WARNING: SurfaceFromCohesives initialized with'
                ' empty element set.')
        return

## ---------------------------------------------------------------------------

class SurfaceFromVRML(TriangleSurface):
    """Work in progress or possibly stalled... Does not work yet at all.
    A surface created from some tri elements read from a VRML (.wrl) file like
    Abaqus/CAE creates them.
    """

    def __init__(self, filename, logfile=None):
        """
        @param filename: File name of a .wrl file to read surface data from.

        @param logfile: deprecated argument, only for compatibility reasons.
        No effect.

        @Note: This is a quite crude implementation of the VRML syntax which
        is sufficient to read output from Abaqus/CAE. For other cases it might
        or might not work.
        """

        TriangleSurface.__init__(self)

        inputFile = open(filename)
        headerline = inputFile.next()
        if not(headerline=="#VRML V2.0 utf8\n"):
            msg("SurfaceFromVRML: Can't read file %s. Not a recognized VRML"
                " dialect.")
            return

        lineNb = 0
        nodes = []
        elements = []

        for line in inputFile:
            line = line.strip()
            lineNb += 1

            # empty or comment line
            if len(line)==0 or line[0]=="#":
                continue

            if line.startswith("NavigationInfo"):
                # ignore this block
                # go on from here....
                # look here:
                # //REDBACK/gero/redback-work/abaqus/perseverance/2010/Analysis/PERSEDEC2010_COUPLE_LR02_YR2011_M10-M03_RESTART/PERSE_YR2011_M10-M03_dU3_020.wrl
                # http://www.web3d.org/x3d/specifications/vrml/vrml97/index.htm
                # http://www.wiley.com/legacy/compbooks/vrml2sbk/ch02/02fig01.htm
                #
                # alternatively we could use some parsers to parse the VRML
                # file:
                # ZestyParser: http://pypi.python.org/pypi/ZestyParser
                # simpleparse: http://simpleparse.sourceforge.net/ http://www.ibm.com/developerworks/linux/library/l-simple.html http://en.wikipedia.org/wiki/Extended_Backus%E2%80%93Naur_Form
                # ...: http://effbot.org/zone/simple-top-down-parsing.htm
                # list of python parsers:
                # http://wiki.python.org/moin/LanguageParsing
                # http://nedbatchelder.com/text/python-parsers.html
                pass
        return

## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------

class ElemFaceSurface(Surface):
    """
    A surface made of faces of a mesh. This is very much related to Abaqus
    meshes: It uses Abaqus' element face numbering.

    In other words: This is a Abaqus continuum element surface.

    @ivar faceEl: {face key: set of element numbers} face key like "S1", "S2"
    @ivar mesh: a abq_model_02.Model instance defining the elements referenced
        in faceEl

    @Note: Might be generalized to use a mesh_01.Mesh instance in the future.
    """
    # Todo:
    #   - union, difference, intersection functions returning a new surface
    #     object

    def __init__(self, faceEl=None, mesh=None):
        """Can also serve as a copy constructor if only the first argument is
        supplied and of type ElemFaceSurface. In this case a deep copy of the
        faceEl item is performed, for self.mesh we simply store the same
        reference.

        Whereas if you supply faceEl as a dictionary then it is just stored as
        is!

        @param faceEl: Optional, initialize self.faceEl with this.
           If this argument is of type ElemFaceSurface and if it's the only
           argument then this function serves as a copy constructor and
           performes a deep copy of the faceEl attribute.
        @param mesh: optional, initialize self.mesh with this
        """
        if isinstance(faceEl, ElemFaceSurface) and mesh is None:
            # copy constructor
            # if the first argument is another ElemFaceSurface-instance
            self.mesh = faceEl.mesh
            self.faceEl = dict(
                (typ, set(elset))
                for typ, elset in faceEl.faceEl.iteritems() )
        elif faceEl is None:
            self.faceEl = dict()
            self.mesh = mesh
        else:
            self.faceEl = faceEl
            self.mesh = mesh

    def countFaces(self):
        return sum(map(len, self.faceEl.itervalues()))

    def addFreeSurfFromElset(self, model, elset=None):
        """Adds to self the perimeter surface / exterior surface
        of the given elset.

        @param model: A L{bae.abq_model_02.Model} or L{bae.mesh_01.Mesh}
           instance.
        @param elset: may be an element set name or a list of element set names,
           anything that Model.getUnionSet accepts as input. If not given or
           None then take the free surface of the whole model.
        @returns: self, so that surf=ElemFaceSurface().addFreeSurfFromElset(...)
           is possible.
        """

        if self.mesh is None:
            self.mesh = model
        elif self.mesh is not model:
            raise ValueError(
                "ElemFaceSurface.addFreeSurfFromElset(): Given model argument"
                " must be identical to self.mesh and it's not.")
        if elset is None:
            elset = set(model.elNodes)
        elif isinstance(model, Model):
            # if we have an abq_model_02.Model as self.mesh==model
            # then we can use model.getUnionSet()
            elset = model.getUnionSet("elset", elset)
        else:
            # if self.mesh==model is a mesh_01.Mesh object then the elset
            # argument must be a set of element numbers
            assert(isinstance(model, Mesh))
            if not isinstance(elset, set):
                # convert lists or so to a set
                newelset = set(elset)
            if not isinstance(iter(newelset).next(), int):
                # items of newelset must be element numbers of type integer
                raise ValueError(
                    "elset argument to ElemFaceSurface.addFreeSurfFromElset()"
                    " must be an iterable of integers, instead its a %s:\n%s"
                    % (type(elset), elset))
            elset = newelset

        faceToElem = defaultdict(list)
        for thisElem in elset:
            for faceNum, face in enumerate(model.elemFacesIter(thisElem)):
                faceToElem[face].append([thisElem,faceNum])

        freeSurfElems = defaultdict(set)
        for face, elems in faceToElem.iteritems():
            if len(elems)==1:
                freeSurfElems['S%d' %(elems[0][1]+1)].add(elems[0][0])
            else:
                pass

        self.faceEl.update(freeSurfElems)
        return self

    def updateFromModel(self, model, surfaceName):
        """Returns a dict of elements that define the surface surf.

        @param model: a abq_model_02.Model instance
        @param surfaceName: Surface name of type string

        @returns: self, so that surf=ElemFaceSurface().addFreeSurfFromElset(...)
           is possible.

        @note: This stores references of the element sets. If you modify this
          data you may also modify the elset the surface referes to.
        """

        if self.mesh is None:
            self.mesh = model
        elif self.mesh is not model:
            raise ValueError(
                "ElemFaceSurface.updateFromModel(): Given model argument"
                " must be identical to self.mesh and it's not.")

        try:
            surfElemSets = model.surface[surfaceName]
        except KeyError:
            raise KeyError("Surface %s is not defined in model." % surfaceName)

        for faceId, elsetName in surfElemSets.iteritems():
            self.faceEl[faceId] = model.elset[elsetName]
        return self

    def storeInAbqModel(self, model, name):
        """Add this surface to abq_model with a specified name, i.e. add
        corresponding elsets and model.surface - dictionary

        @param model: abq_model_02.Model instance
        @param name: name of the new surface
        """
        if name in model.surface:
            raise ValueError("Surface with the specified name %s already exists"
                             " in the model." % name)
        surfdict = dict()
        for typ in sorted(self.faceEl):
            elsetName = "%s_%s" % (name, typ)
            if elsetName in model.elset:
                cnt = 0
                template = "%s_%%d" % elsetName
                while elsetName in model.elset:
                    cnt += 1
                    elsetName = template % cnt
            model.elset[elsetName] = self.faceEl[typ]
            surfdict[typ] = elsetName
        model.surface[name] = surfdict

    def getTriangleSurface(self, quadTetTo4Tris=False):
        """Return a L{TriangleSurface} created from self.

        The surface must consist of faces of tet elements (C3D4 or C3D10M).
        Only the corner nodes of the tets are considered if the parameter
        quadTetTo4Tris is not set to True.

        The node coords of this surface are a reference (no copy) to those of
        the abaqus model self.mesh. If self.mesh moves, so does this surface.

        @param quadTetTo4Tris: If this is set to true then create four tris
        for each quadratic tet face.

        @Note: This implementation has to be cleaned up. Currently you should
        only use the quadTetTo4Tris parameter if the mesh consists of quadratic
        tets only.
        """

        newsurf = TriangleSurface()

        if quadTetTo4Tris:

            # collect tri nodes
            newTriNodes = list()
            # k = nodesIdx4ShellWedges[i][j] gives for the i-th wedge elem of
            # four to be attached to any one tet face the j-th node index if
            # it's on face S1.
            # For any of the four tet faces (S1, S2, S3 or S4) you'll get the
            # node index from
            # n = model.elShapeFaceNodes["TET_Q"][][k] and the node number from
            # model.elNodes[elem][n]
            nodesIdx4ShellWedges = [[0,5,3], [5,2,4], [4,1,3], [5,4,3]]
            # for quad tet elements: edge (two node indexes) to mid node index
            tetMesh = self.mesh

            for typ in sorted(self.faceEl):
                faceIdx = int(typ[1:]) - 1

                tetNodeIds4Face = tetMesh.elShapeFaceNodes["TET_Q"][faceIdx]
                nodesIdx4ShellWedgesForCurrentTetFace = [
                    [tetNodeIds4Face[i] for i in nodeIdsFace0]
                    for nodeIdsFace0 in nodesIdx4ShellWedges ]

                for tetElem in self.faceEl[typ]:
                    tetNodes = tetMesh.elNodes[tetElem]

                    # update newTriNodes list
                    newTriNodes.extend(
                        [tetNodes[i] for i in nodeIds]
                        for nodeIds in nodesIdx4ShellWedgesForCurrentTetFace )

            # store tri nodes in newsurf
            newsurf.triNodes.update(
                (i+1, nodes)
                for i, nodes in enumerate(newTriNodes) )
            nodestocopy = set( i for nodes in newTriNodes for i in nodes )
            del newTriNodes

        else: # if not quadTetTo4Tris

            triElementFaces={'S1':[0,2,1],'S2':[0,1,3],'S3':[1,2,3],'S4':[0,3,2]} 

            nodestocopy = set()
            elnum = 1
            for faceType in set(self.faceEl).intersection(
                triElementFaces):
                nodeIds = triElementFaces[faceType]
                for el in self.faceEl[faceType]:
                    if self.mesh.elShape[el] not in ('TET_L', 'TET_Q'):
                        raise NotImplementedError(
                            "%s.%s.getTriangleSurface() can not create a"
                            " surface from elements of shape %s."
                            % (__name__, self.__class__.__name__,
                               self.mesh.elShape[el]))
                    nodes = [self.mesh.elNodes[el][i] for i in nodeIds]
                    newsurf.triNodes[elnum] = nodes
                    nodestocopy.update(nodes)
                    elnum += 1

        # store nodes in newsurf
        for node in nodestocopy:
            try:
                coords = self.mesh.nodeCoords[node]
                newsurf.nodeCoords[node] = coords
            except KeyError:
                pass

        if len(newsurf.triNodes)==0:
            msg('WARNING: SurfaceFromMeshTris initialized with'
                ' empty element set.')
        return newsurf

    def update(self, other):
        """
        Constructs a boolean add operator for Abaqus continuum element surfaces.
        
        @param other: ElemFaceSurface object to be added.
        """
        assert isinstance(other, ElemFaceSurface)

        if self.mesh is None:
            self.mesh = other.mesh
        elif self.mesh is not other.mesh:
            raise ValueError(
                "ElemFaceSurface.update(): mesh of other surface"
                " must be identical to self.mesh and it's not.")
        for face, elems in other.faceEl.iteritems():
            try:
                oldElems = self.faceEl[face]
            except KeyError:
                oldElems = set()
                self.faceEl[face] = oldElems
            oldElems.update(elems)

    def difference_update(self, other):
        """
        Constructs a boolean subtract operator for Abaqus continuum element
        surfaces.

        @param other: ElemFaceSurface object to be subtracted.
        """
        assert isinstance(other, ElemFaceSurface)

        if self.mesh is None:
            self.mesh = other.mesh
        elif self.mesh is not other.mesh:
            raise ValueError(
                "ElemFaceSurface.difference_update(): mesh of other surface"
                " must be identical to self.mesh and it's not.")
        for face, elems in other.faceEl.iteritems():
            try:
                oldElems = self.faceEl[face]
            except KeyError:
                continue
            else:
                oldElems.difference_update(elems)

    def intersection_update(self, other):
        """Constructs a boolean intersect operator for Abaqus continuum element
        surfaces. I.e. remove all surface parts from self that are not found in
        other.
        
        @param other: ElemFaceSurface object to be intersected with self.
        """
        assert isinstance(other, ElemFaceSurface)

        if self.mesh is None:
            self.mesh = other.mesh
        elif self.mesh is not other.mesh:
            raise ValueError(
                "ElemFaceSurface.difference_update(): mesh of other surface"
                " must be identical to self.mesh and it's not.")
        for face, elems in other.faceEl.iteritems():
            try:
                oldElems = self.faceEl[face]
            except KeyError:
                del self.faceEl[face]
            else:
                oldElems.intersection_update(elems)

## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------

def getSimpleSurfacesFromMeshTris(
    abqModel, elset=None,
    mintris=1, splitangle=None, skipMoebiusCheck=False, joinParts=False,
    logfile=None, debugElset=None):

    """Creates a list of TriangleSurface objects of the triangles in elset
    in the abqModel.

    Each surface is simply connected (no T-junction) and all
    triangles are connected. The normals of the triangles in those surfaces are
    already calculated (self.triNormal has been created) and edgeToTri and
    borderEdges are there.


    @param abqModel: bae.abq_model_02.Model object

    @param elset: Elset name or element numbers or combination
              (Everything acceptable as input to abqModel.getUnionSet().)
              if None or not specified, take all tri elements out of the model
    @param mintris: Discard all surfaces with less then mintris triangles

    @param splitangle: Stop at edges with an angle greater than
              splitangle (in degrees)
    @param skipMoebiusCheck: No check if the surface has two sides

    @param joinParts: If False don't try to join parts of the surface, otherwise
    parts that are connected over an edge with an angle smaller than
    _joinPartsAngleDeg are possibly reconnected. (conditions apply)

    @param debugElset: if not None add diagnostic elsets to this dict
    {elset name: elset}.

    @param logfile: deprecated argument, only for compatibility reasons.
    No effect.

    @return: tuple:
     - surflist .. list of TriangleSurface objects
     - skipped_parts_cnt .. number of parts, that were skipped because < mintris
     - tris_in_skipped_parts .. list of corresponding tri elements
    """

    checkModelModuleVersion(abqModel, "%s.getSimpleSurfacesFromMeshTris"
                            % __name__)

    surflist = list()
    skipped_parts_cnt = 0
    tris_in_skipped_parts = list()

    mintris = max(mintris, 1)

    # create a TriangleSurface object of the given elset
    # This is a misuse and probably not really necessary but it's done for
    # historical reasons. And it violates the assumption that the
    # TriangleSurface objects are simply connected.
    sAll = SurfaceFromMeshTris(abqModel, elset)

    # local flag and set used to sum up
    all_moebius_flag_cnt = 0
    all_moebiusStopEdges = set()

    # store edges with more than two tris attached to it
    multiTriEdges = set()
    # for tris on multiTriEdges: {tri : part index in trisWithNormalList}
    mteTriPartIdx = dict()

    # set of tris for which a normal has not been calculated yet
    noNormalTris = set(sAll.triNodes)

    sAll.triNormal = dict()

    while len(noNormalTris)>0:

        # get any one element as a starting point
        startTriId = iter(noNormalTris).next()

        # we want the flag and the set for each part surface alone
        sAll.moebius_flag = False
        sAll.moebiusStopEdges = set()

        # calculate normals starting at that point
        trisWithNormalList, newMTEs, newMteTris = getConnectedTris(
            sAll, startTriId, splitangle, skipMoebiusCheck)
        multiTriEdges.update(newMTEs)
        mteTriPartIdx.update([(ttt, len(surflist)) for ttt in newMteTris])

#         trisWithNormalList = sAll.calcNormalsOnPart(
#             startTriId, splitangle, skipMoebiusCheck)
        noNormalTris.difference_update(trisWithNormalList)

        # sum up data for sAll
        if sAll.moebius_flag:
            all_moebius_flag_cnt += 1
        all_moebiusStopEdges.update(sAll.moebiusStopEdges)

        # create new surface
        newsurf = SurfaceFromMeshTris(abqModel, trisWithNormalList)

        # create other features of new part surface
        newsurf.triNormal = dict([
                (el, sAll.triNormal[el])
                for el in trisWithNormalList ])
        newsurf.moebius_flag     = sAll.moebius_flag
        newsurf.moebiusStopEdges = sAll.moebiusStopEdges

        # add to list
        surflist.append(newsurf)

    msg("This triangle elset has been split into %d surface parts."
        % len(surflist))

    ## now try to join parts together again
    if joinParts:
        msg("Try to reconnect parts at edges with more than two triangles."
            " There are %d such edges. Minimum join angle: %g."
            % (len(multiTriEdges), _joinPartsAngleDeg))

        joinPartsAngleCos = cos(_joinPartsAngleDeg*pi/180)
        edgeToTri = sAll.getEdgeToTri()

        ## possibleConnections is a dict { connection key : value list }
        # connection key is a set of two part indexes in trisWithNormalList
        # value list components:
        #  . sum of edge connection values = values of the angle weighting
        #    function _joinPartsAngleWeightingFunction()
        #  . list of contradicting connections
        #  . min number of tris in parts
        #  . twist normals factor: -1 when the normals have to be twisted and
        #    1 if not
        possibleConnections = dict()
        for edgeCnt, thisEdge in enumerate(multiTriEdges):

            # ignore all tris on this edge that belong to a part with more than
            # one tris on this edge
            # {connected tri : part index in trisWithNormalList}
            connTriToPart = dict([
                    (tri, mteTriPartIdx[tri])
                    for tri in edgeToTri[thisEdge]])
            connPartToNbTri = defaultdict(int) # {part index  }
            for pn in connTriToPart.values():
                connPartToNbTri[pn] += 1

            # list of tris connected to this edge
            connectedTris = [tri for tri in edgeToTri[thisEdge]
                             if connPartToNbTri[connTriToPart[tri]]==1]
            if len(connectedTris)<2:
                continue

            # determine the orientation of the tris with respect to the
            # connecting edge
            orderedEdge = tuple(thisEdge)
            orderedEdgeTwisted = tuple(orderedEdge[::-1])
            triTwistNormalsFactor = list()
            for tri in connectedTris:
                triOrderedEdges = set([
                        ee[1] for ee in sAll.triEdgesOrderedNodesIter(
                            sAll.triNodes[tri])])
                if orderedEdge in triOrderedEdges:
                    triTwistNormalsFactor.append(1.0)
                elif orderedEdgeTwisted in triOrderedEdges:
                    triTwistNormalsFactor.append(-1.0)
                else:
                    raise Exception(
                        "Multiple connected edge %s not on tri %d with nodes"
                        " %s anymore."
                        % (orderedEdge, tri, sAll.triNodes[tri]))

            # all possible part connections
            for i1,t1 in enumerate(connectedTris):
                part1 = mteTriPartIdx[t1]
                normal1 = vector_scale(sAll.triNormal[t1],
                                       triTwistNormalsFactor[i1])
                for i2 in range(i1+1, len(connectedTris)):
                    t2 = connectedTris[i2]
                    part2 = mteTriPartIdx[t2]
                    twistNormalsFactor2 = -triTwistNormalsFactor[i2]
                    normal2 = vector_scale(sAll.triNormal[t2],
                                           twistNormalsFactor2)
                    connTwistNormalsFactor = (triTwistNormalsFactor[i1]
                                              * twistNormalsFactor2)

                    # ignore all connections with an angle greater than
                    # _joinPartsAngleDeg
                    t1t2AngleCos = dot(normal1, normal2)
                    if t1t2AngleCos<joinPartsAngleCos:
                        continue

                    # get or initialize the value list for this connection
                    thisConn = frozenset((part1, part2))
                    if thisConn in possibleConnections:
                        thisConnVal = possibleConnections[thisConn]
                    else:
                        thisConnVal = [0.0, set(), 0, connTwistNormalsFactor]
                        possibleConnections[thisConn] = thisConnVal

                    # Check 4th component of connection value list:
                    # twist normals factor: -1 when the normals have to be
                    # twisted and 1 if not. 0 if found to be contradicting.
                    # For such a connection we don't need further information
                    if thisConnVal[3]*connTwistNormalsFactor < 0.5:
                        thisConnVal[3] = 0
                        continue

                    # Update 3rd component of connection value list:
                    # min number of tris in parts
                    thisConnVal[2] = min(thisConnVal[2],
                                         len(surflist[i1].triNodes),
                                         len(surflist[i2].triNodes))

                    # Update 1st component of connection value list:
                    # sum of edge connection values = values of the angle
                    # weighting function
                    thisConnVal[0] += _joinPartsAngleWeightingFunction(
                        t1t2AngleCos)

                    # Update 2nd component of connection value list:
                    # list of contradicting connections
                    otherConnTris = set(connectedTris)
                    otherConnTris.remove(t1)
                    otherConnTris.remove(t2)
                    for t3 in otherConnTris:
                        thisConnVal[1].add(
                            frozenset((connTriToPart[t1], connTriToPart[t3])))
                        thisConnVal[1].add(
                            frozenset((connTriToPart[t2], connTriToPart[t3])))
        ## finished possibleConnections dict
        msg("Found %d possible part connections in the first place."
            % len(possibleConnections))

        # remove impossible connections
        # 4th component: twist normals factor == 0 if found to be contradicting
        for key in possibleConnections.keys():
            if possibleConnections[key][3] < 0.5:
                del possibleConnections[key]

        # calculate rating, 1st component (sum of edge connection values)
        # rating = rating of self*nbPossibleConnects minus sum of contradicting
        #          ratings
        nbPossibleConnects = len(possibleConnections)
        ratingConnect = list()
        for key, connVal in possibleConnections.iteritems():
            connVal[1].intersection_update(possibleConnections)
            rating = connVal[0] * (nbPossibleConnects-1)
            for contraKey in connVal[1]:
                rating -= possibleConnections[contraKey][0]
            ratingConnect.append((rating,key))

        # delete contradicting connections from top to bottom of rating list
        ratingConnect.sort()
        for rating, key in ratingConnect:
            try:
                for contraKey in possibleConnections[key][1]:
                    del possibleConnections[contraKey]
            except KeyError:
                # this connection has already been deleted
                pass

        # sort possibleConnections, result in list joinConnections:
        # transform frozenset keys to tuples (high number, low number)
        joinConnections = [list(key) for key in possibleConnections]
        for partList in joinConnections:
            partList.sort(reverse=True)
        joinConnections.sort(reverse=True)

        # Remove duplicates, each surface is at most joined into one other
        # (the second, target surface having a lower number).
        # Because joinConnections is sorted, if a surface could be joined with
        # several others, the lowest surface number is chosen.
        joinConnections = dict(joinConnections)
        joinConnections = list(joinConnections.iteritems())
        joinConnections.sort(reverse=True)
        
        # join the surfaces
        deleteParts = set()
        for partList in joinConnections:
            deleteParts.add(partList[0])

            msg("Merging surface part %d into part %d."
                % (partList[0]+1, partList[1]+1))
            
            # merging partList[0] into partList[1]
            surf = surflist[partList[1]]
            tomerge = surflist[partList[0]]
            surf.triNodes.update(tomerge.triNodes)
            surf.nodeCoords.update(tomerge.nodeCoords)
            surf.triNormal.update(tomerge.triNormal)
            surf.moebius_flag = surf.moebius_flag or tomerge.moebius_flag
            surf.moebiusStopEdges.update(tomerge.moebiusStopEdges)
            if surf.moebius_flag and tomerge.moebius_flag:
                all_moebius_flag_cnt -= 1

        # delete the merged in parts
        surflist = [surf
                    for i, surf in enumerate(surflist)
                    if i not in deleteParts]

        msg("Joined parts of the surface %d times, leaving a total of %d"
            " parts." % (len(joinConnections), len(surflist)))


    # check mintris
    newsurflist = list()
    for surf in surflist:
        if len(surf.triNodes) >= mintris:
            newsurflist.append(surf)
        else:
            skipped_parts_cnt += 1
            tris_in_skipped_parts.extend(surf.triNodes)
    surflist = newsurflist

    if all_moebius_flag_cnt>0:
        msg("WARNING from getSimpleSurfacesFromMeshTris(): %d surfaces are"
            " twisted like the Moebius strip, they have only one side."
            % all_moebius_flag_cnt)

    return surflist, skipped_parts_cnt, tris_in_skipped_parts


def getConnectedTris(surf, startTriId,
                     splitangle=None, skipMoebiusCheck=False):
    """Not meant to be called from "outside"
    getSimpleSurfacesFromMeshTris() uses it.
    """
    extraDoc = """\
    This function is meant to replace TriangleSurface.calcNormalsOnPart()
    ...and it's derived from that

    Make sure len(surf.triNodes)>0, and surf.moebiusStopEdges is a set.
    surf.triNormal must already be a dict!

    creates object-attribute triNormal:
    a dict {element:[comp1, comp2, comp3]} with [comp...] being the vector
    components of the *unscaled* normal of that triangle

    This function calculates the elements normal vectors of the surface
    triangles in such a way that all normals point to the same side of the
    surface. It starts with the element given by *startTriId* and proceeds
    with all neighbouring tris as long as they are simply connected.

    Additionally the orientation of the tri element is saved by means of
    reordering its nodes in the surf.triNodes-dictionary in
    counter-clockwise succession if you look from "above" (positive side).
    (The calculated normal points "up" to the positive side.)

    At an edge with more than two triangles attached, the algorithm stops.
    At triangles that have their normal calculated already
    (newTri in surf.triNormal) the algorithm stops as well.
    If the angle in degree between two neighbouring triangles is greater
    than *splitangle*, the algorithm stops at that edge. splitangle==None
    means: don't check angle.

    Simply connected means there are max two triangles connected to each
    edge.

    If skipMoebiusCheck==True, no check is performed, if the surface
    has two sides and the normals are determinable uniquely. Otherwise if
    contradicting results for the same normal are found, surf.moebius_flag
    is set to True and the corresponding edge is added to
    surf.moebiusStopEdges. Additionally a warning is logged in the log
    file.

    The function returns a tuple:
     1. a list of all triangles that are connected directly or indirectly to
        startTriId
     2. a set of edges with more than two tris attached to it
     3. a set of triangles belonging to this part connected to any edge in the
        former set of multiply connected edges
    """

    if splitangle != None:
        splitanglecos = cos(float(splitangle)*pi/180.)

    edgeToTri = surf.getEdgeToTri()
    assert len(surf.triNodes)>0

    trisWithNormalList = list()
    startTriNodes = surf.triNodes[startTriId]
    startNormal = surf.getTriNormalSc(startTriNodes)
    surf.triNormal[startTriId] = startNormal
    trisWithNormalList.append(startTriId)

    # trisToFindNeighbours is a list of (TriId, TriNodes, Normal)-tuples
    # of those tris for which the neighbours have to be examined yet. The
    # normals of the tris in this list have already been calculated.
    # the TriNodes-item contains the node numbers in counter-clockwise
    # succession if you look from "above" (positive side)
    trisToFindNeighbours = deque()  # a double-ended queue
    trisToFindNeighbours.append((startTriId, startTriNodes, startNormal))

    # multiTriEdges is a set of edges with more than two tris attached to it
    multiTriEdges = set()
    # mteTris is a set of tris belonging to this part and connected to an edge
    # in multiTriEdges
    mteTris = set()

    while len(trisToFindNeighbours) > 0:
        # get the first (oldest) tri in the list
        startTriId, startTriNodes, startNormal = (
            trisToFindNeighbours.popleft())

        for thisEdge, edgeNodes in surf.triEdgesOrderedNodesIter(
            startTriNodes):

            trisOnThisEdge = edgeToTri[thisEdge]
            if len(trisOnThisEdge)>2:
                # more than two tris on this edge:
                # ignore this edge at this point, store for later processing
                multiTriEdges.add(thisEdge)
                mteTris.add(startTriId)
                continue
            elif len(trisOnThisEdge)<2:
                # less than two: edge on border, no other tri connected
                continue

            # find the other tri connected to this edge (not startTriId)
            for newTri in trisOnThisEdge:
                if newTri != startTriId: break


            # if no check and normal already calculated, skip the rest
            newTriHasNormalAlready = (newTri in surf.triNormal)
            if skipMoebiusCheck and newTriHasNormalAlready:
                continue

            # find third node not on this edge
            newTriNodes = list(surf.triNodes[newTri])
            for i in thisEdge:
                newTriNodes.remove(i)
            # append nodes on the edge in correct order
            newTriNodes.append(edgeNodes[1])
            newTriNodes.append(edgeNodes[0])
            # calculate new normal
            newNormal = surf.getTriNormalSc(newTriNodes)

            # check max angle between neighbouring normals
            if (splitangle != None and
                dot(startNormal, newNormal)<splitanglecos):
                continue

            # if normal already calculated, compare and goto next
            if newTriHasNormalAlready:
                if dot(newNormal, surf.triNormal[newTri]) < 0.0:
                    surf.moebius_flag = True
                    surf.moebiusStopEdges.add(thisEdge)
                    msg("WARNING: The normals flip between element %d and"
                        " %d." % (startTriId, newTri))

                # in any case keep the already calculated normal and
                # proceed with the next tri
                continue

            # update surf.triNormal with new normal
            surf.triNormal[newTri] = newNormal
            trisWithNormalList.append(newTri)
            # update new node order to surf.triNodes
            surf.triNodes[newTri] = newTriNodes
            # append this neighbour to trisToFindNeighbours-list
            trisToFindNeighbours.append((newTri, newTriNodes, newNormal))

    # fini
    return (trisWithNormalList, multiTriEdges, mteTris)

_joinPartsAngleDeg = 15.0
def _joinPartsAngleWeightingFunction(angleCos):
    """Not meant to be called from "outside"
    angle weighting function used by getSimpleSurfacesFromMeshTris():
    at alpha = 0 => weighting function f = 1.0
    at alpha = 1 => weighting function f = 0.0
    => f = sum_i a_i cos^i alpha
    a_i = aList
    """
    # alpha max = 60 deg
#     cosMaxAngle = cos(60.0 *pi/180.0)
#     aList = [ -3.52478665,  24.41584112, -61.15783254,  64.43440953, -23.16763147]
    # alpha max = 15 deg
    cosMaxAngle = cos(15.0 *pi/180.0)
    aList = [ 1673924.77719293, -6884242.83380361, 10617036.49685724, -7277187.95816989, 1870470.51792333]

    if angleCos<cosMaxAngle:
        return 0.0
    else:
        return sum([a_i * angleCos**i
                    for i, a_i in enumerate(aList)])


#####################################################################
# Additional import at the end of the file
# This has to be done here, because TriangleSurfaceInVolume needs the
# definition of TriangleSurface.
from bae.surface_in_volume_02 import TriangleSurfaceInVolume, \
    SurfaceInVolFromMeshTris, SurfaceInVolFromTriangleSurface, \
    LengthOfSetNotOneException
#####################################################################
