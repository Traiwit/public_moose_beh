r"""gs_mesh_01.py

This module defines a class to create a ground support mesh as abaqus model

Usage
=====
  >>> from bae.generatemesh.gs_mesh_01 import GSMesh4545
  >>> m = GSMesh4545()
  >>> help(m) for more information


Example
=======

>>> from bae.abq_model_02 import Model
>>> from bae.generatemesh.tunnel_01 import StraightTunnel, XSectCircBack
>>> from bae.generatemesh.gs_mesh_01 import \
>>>     GSMesh4545, GSBoltsOnMesh, GSscsOnMesh, GSscsFree
>>>
>>> from math import sqrt
>>> #---------------------------------------------------------
>>>
>>> model = Model()
>>>
>>> #-- tunnel geometry
>>> # tunnel cross section
>>> # originally width=4.4, height=5.6
>>> # 2010NOV16, David: add 0.3 to each
>>> xSect = XSectCircBack(width=4.7, height=5.9, backRadius=3)
>>>
>>> # tunnel length:
>>> # DFN blocks are about 7-8m diameter, we should model approx 3 of them
>>> # => model length (and diamater as well) 25m
>>> # tunnel heading:
>>> # correct direction vector (with inclined slope): [0.789,-0.614,0.031]
>>> # projected on the xy-plane: [0.789,-0.614,0.000] = 127.89 deg
>>> tunnelContour = StraightTunnel(xSect, length=25, heading=127.89)
>>> tunnelOrigin = [3352.767674424111,4325.757436761255,4152.956413536747]
>>>
>>> #-- GS mesh
>>> # according to David: D=3mm, meshSize=50mm, sigma_yield=500MPa
>>> # But there is also: D=3.6mm, meshSize=75mm
>>> # The mesh is in a 100mm shotcrete layer => offset -50mm
>>> meshContour = tunnelContour.getOffsetContour(offset=-0.05)
>>> # It starts 1.7m (Fausto said 1.5 to 2m) above the ground
>>> meshWidth = 2.0*meshContour.sFromBottom(1.7)
>>> gsMesh = GSMesh4545(meshContour,
>>>                     approxMeshSize=0.075, meshWidth=meshWidth,
>>>                     elsetName="MESH")
>>>
>>> #-- connected bolts
>>> boltLength = 3.9
>>> bolts = GSBoltsOnMesh(beamMeshSize=0.05)
>>> gsMesh.connectBolts(bolts)
>>> # bolts D=25mm, really 2.25m in rock, total length=2.4m, in old submodel: 2m
>>> # Fausto: bolt pattern 1.0m or 1.2m diagonal
>>> boltPatternZmin = 0.5 # lowest bolt row (dont go deeper by all means)
>>> horizontalBoltAtZ = 1.5 # all bolts below point downwards, all above upwards
>>> boltPatternSpacing = 1.0  # diagonal to the tunnel
>>> boltPatternDiam = sqrt(2)*boltPatternSpacing # along the tunnel
>>>
>>> def boltAnglesFunction(s):
>>>     angle = (s/meshContour.sFromBottom(horizontalBoltAtZ))*90.0
>>>     angleToTunnel = 90.0
>>>     return (angle, angleToTunnel)
>>>
>>> # first of two sets with a row at the very top (s=0)
>>> (boltsSCoords, boltsSIndexList) = bolts.addOrthoPattern(
>>>     yStart=0.5*boltPatternDiam, yWidth=boltPatternDiam,
>>>     sStart=0.25*boltPatternDiam, sWidth=boltPatternDiam,
>>>     sMax=bolts.contour.sFromBottom(boltPatternZmin),
>>>     anglesFunction=boltAnglesFunction,
>>>     boltLength=boltLength, excess=0.06)
>>> print "created first set of bolts rows at s=...\n%s" % boltsSCoords
>>>
>>> # second of two sets with the first row at s=0.5*boltPatternDiam
>>> (boltsSCoords, boltsSIndexList) = bolts.addOrthoPattern(
>>>     yStart=boltPatternDiam, yWidth=boltPatternDiam,
>>>     sStart=0.75*boltPatternDiam, sWidth=boltPatternDiam,
>>>     sMax=bolts.contour.sFromBottom(boltPatternZmin),
>>>     anglesFunction=boltAnglesFunction,
>>>     boltLength=boltLength, excess=0.06)
>>> print "created second set of bolts rows at s=...\n%s" % boltsSCoords
>>>
>>> #-- connected shotcrete
>>> # shotcrete all over and to the bottom
>>> sEnd = meshContour.sFromBottom(0.0)
>>> newShotcrete = GSscsOnMesh(sStart=0.0, sEnd=sEnd, elsetName="SHOTCRETE_OLD")
>>> gsMesh.connectShotcrete(newShotcrete)
>>>
>>> #-- Old GS with mesh, bolts and shotcrete: finish
>>> gsMesh.addToModel(model)
>>>
>>>
>>> # #-- New GS: Shotcrete layer not connected to anything
>>> # scfrContour = tunnelContour.getOffsetContour(offset=-0.15)
>>> # sStart = scfrContour.sFromBottom(2.1)
>>> # sEnd = scfrContour.sFromBottom(0.0)
>>> # newShotcrete = GSscsFree(
>>> #     scfrContour, sStart=sStart, sEnd=sEnd,
>>> #     meshSize=0.2, elsetName="SHOTCRETE_NEW")
>>> # newShotcrete.addToModel(model)
>>>
>>> # renumbering
>>> model.renumber(nodeStart=5000000, elemStart=5000000)
>>> model.translate(tunnelOrigin)
>>>
>>> # intermediate output
>>> outputName = "AMP_GS11_cheap-02.inp"
>>> print "... writing input file %s" % outputName
>>> model.write(outputName)
>>>
>>> print ("Now please select and save the nset footnodes and the elset warped"
>>>        " with cae and the write_nset/write_elset functions from the"
>>>        " bae.cae_01 module. When the files nset_footnodes.inp and"
>>>        " elset_warped.inp are there...")
>>> raw_input("Then please press ENTER")
>>>
>>> #-- some extra corrections
>>>
>>> # pull bottom nodes of shotcrete to z=0
>>> print "pulling the footnodes to the floor"
>>> m2=Model().read("nset_footnodes.inp")
>>> for node in m2.nset["footnodes"]:
>>>     model.nodeCoords[node][2] = 0.0
>>>
>>> # split warped quad shell elements into tri elements
>>> print "splitting the wraped elements into triangles"
>>> m2=Model().read("elset_warped.inp")
>>> warpedElems = m2.elset["warped"]
>>> nextElem = max(model.elNodes) + 1
>>> newTriElems = dict()
>>> for el in warpedElems:
>>>     nodes = model.elNodes[el]
>>>     model.removeElem(el, removeNodes=False, updateElsets=True)
>>>     newTriElems[nextElem] = [nodes[i] for i in [0,1,3]]
>>>     nextElem += 1
>>>     newTriElems[nextElem] = [nodes[i] for i in [1,2,3]]
>>>     nextElem += 1
>>> model.updateElems(newTriElems, "S3")
>>> model.elset["SHOTCRETE_OLD"].update(newTriElems)
>>>
>>> #-- fini
>>> outputName = "AMP_GS11_cheap-02.inp"
>>> print "... writing input file %s" % outputName
>>> model.write(outputName)

"""

__version__ = "1.1"

_version_history_ = """\
Versions:
=========

gs_mesh_01
1.0 : GP: derived from Perilya GS project scripts used for the Antamina GS
          project in 2010
1.1 GP switch to model.insertElem, model.insertNode
"""

from math import sqrt, sin, cos, pi
from itertools import izip, combinations
from collections import defaultdict
from bae.vecmath_01 import vector_plus, vector_scale
from bae.generatemesh.tunnel_01 import Tunnel
from bae.log_01 import msg


def _getLastNodeNum(model, lastNodeNum=0):
    try:
        node_num = max(model.nodeCoords.iterkeys())
    except ValueError:
        node_num = lastNodeNum
    return node_num

def _getNextElem(model, nextElem=None):
    # check nextElem
    if nextElem is None:
        if len(model.elNodes):
            nextElem = max(model.elNodes.iterkeys())+1
        else:
            nextElem = 1
    return nextElem


################################################################

class GSMesh(object):
    """Base class for GS mesh applied to a certain tunnel contour.
    """

    def __init__(self, contour, meshWidth, elsetName="MESH",
                 numNodesBetweenCrossings=0, beamElType='B31'):
        """
        Constructor

        @param contour: the contour of the mesh of type Tunnel
        @param meshWidth: The length of the mesh along the circumefernce of the
        tunnel cross section.
        @param numNodesBetweenCrossings: for the fem mesh, ignored at the moment
        """
        self.contour = contour
        self.elsetName = elsetName
        self.length = float(self.contour.length)
        self.halfMeshWidth = 0.5*meshWidth
        self.numNodesBetweenCrossings = numNodesBetweenCrossings
        self.beamElType = beamElType
        self.meshbolts = dict()
        self.shotcrete = None

    def connectBolts(self, meshbolts, boltPlatesElType="S3",
                     boltPlatesElsetName="BOLTPLATES"):
        assert isinstance(meshbolts, GSBoltsOnMesh)
        self.meshbolts = meshbolts
        self.boltPlatesElType = boltPlatesElType
        self.boltPlatesElsetName = boltPlatesElsetName
        meshbolts.contour = self.contour
        meshbolts.connectedMesh = self

    def connectShotcrete(self, shotcrete):
        assert isinstance(shotcrete, GSscsOnMesh)
        self.shotcrete = shotcrete
        shotcrete.contour = self.contour
        # note that the derived classes overload this function to
        # calculate numFullRows


class GSMesh4545(GSMesh):
    """Mesh with a diagonal orientation.

    If you just want mesh:
    >>> mesh = GSMesh4545(...)
    >>> mesh.addToModel(model)

    Mesh and Bolts:
    >>> from bae.abq_model_02 import Model
    >>> from bae.generatemesh.tunnel_01 import StraightTunnel, XSectCircBack
    >>> from bae.generatemesh.gs_mesh_01 import GSMesh4545, GSBoltsOnMesh
    >>>
    >>> model = Model()
    >>>
    >>> xSect = XSectCircBack(width=4.4, height=5.6, backRadius=6)
    >>> tunnelContour = StraightTunnel(xSect, length=15, heading=45.0)
    >>>
    >>> #-- GS mesh 1
    >>> meshContour = tunnelContour.getOffsetContour(offset=-0.05)
    >>> gsMesh = GSMesh4545(meshContour,
    >>>                     approxMeshSize=0.15, meshWidth=12)
    >>>
    >>> #-- connected bolts
    >>> bolts1 = GSBoltsOnMesh(beamMeshSize=0.2)
    >>> bolts1.addRow(0, 0, length=4, spacing=5.0, excess=0.05)
    >>> bolts1.addRow(10, 20, length=4, spacing=5.0, excess=0.05)
    >>> bolts1.addRow(31, 30, length=4, spacing=5.0, excess=0.05)
    >>> bolts1.addRow(-10, -20, length=4, spacing=5.0, excess=0.05)
    >>> bolts1.addRow(-31, -30, length=4, spacing=5.0, excess=0.05)
    >>> gsMesh.connectBolts(bolts1)
    >>> gsMesh.addToModel(model)
    >>>
    >>> #-- fini
    >>> model.write("GS_mesh.inp")

    @ivar meshDy: meshspacing in the y direction along the tunnel
    @ivar meshDs: meshspacing in the s direction along the contour
    """
    def __init__(self, contour, approxMeshSize, meshWidth, elsetName="MESH",
                 numNodesBetweenCrossings=0, beamElType='B31'):
        """
        Constructor

        @param contour: the contour of the mesh of type Tunnel
        @param approxMeshSize: the approximate mesh size of the ground support
        mesh. This will be adapted to fit in an integer number of mesh pieces
        on the whole length of the tunnel. This is the edge length of the
        square mesh openings. The mesh nodes are 1.4 times this length away
        from each other in the direction of the tunnel. See the meshDy
        class instance attribute.
        @param meshWidth: The length of the mesh along the circumference of the
        tunnel cross section.
        @param elsetName: for the beam elements
        @param numNodesBetweenCrossings: for the fem mesh, being ignored at the
        moment.
        """
        GSMesh.__init__(self, contour, meshWidth, elsetName,
                        numNodesBetweenCrossings, beamElType)
        self.numMeshIntervalls = int(round(self.length/(
            approxMeshSize*sqrt(2.0))))
        self.meshDy = self.length / self.numMeshIntervalls
        self.meshDs = self.meshDy
        self.numFullRows = int(self.halfMeshWidth/self.meshDy)+1

    def _newFullRow(self, model, s):
        """Add a new row of nodes to model at coordinate s and return the new
        node numbers"""
        new_nodes = list()
        node_num = _getLastNodeNum(model)
        for i in range(self.numMeshIntervalls+1):
            node_num += 1
            model.nodeCoords[node_num] = self.contour.ys2xyzContour(
                y=i*self.meshDy, s=s)
            new_nodes.append(node_num)
        return new_nodes

    def _newMidRow(self, model, s):
        """Add a new row of mid nodes to model at coordinate s and return the
        new numbers"""
        new_nodes = list()
        node_num = _getLastNodeNum(model)
        for i in range(self.numMeshIntervalls):
            node_num += 1
            model.nodeCoords[node_num] = self.contour.ys2xyzContour(
                y=(0.5+i)*self.meshDy, s=s)
            new_nodes.append(node_num)
        return new_nodes

    def _insertOneSide(self, model, top_nodes, otherside_nodes, dir_mesh_ds,
                       nextElem):
        """
        dir_mesh_ds>0 -> direction pos s, dir_mesh_ds<0 -> direction neg s

        otherside_nodes must be None on the first call. On the second call it
        must be the the second item of the resulting tuple of the first call.
        """

        if isinstance(self.shotcrete, GSscsOnMesh):
            numFullRows = max(self.numFullRows, self.shotcrete.lastFullRow)
            scLastFullRow = self.shotcrete.lastFullRow
            scTopRow = self.shotcrete.topRow
        else:
            numFullRows = self.numFullRows
            scLastFullRow = -10  # not applicable
            scTopRow = 1E6
        s_list = [(i-1, dir_mesh_ds*i) for i in range(1, numFullRows)]

        last_nodes = top_nodes
        beforelast_nodes = otherside_nodes

        elNodesMesh = dict()
        elNodesPlates = dict()
        elNodesSc4 = dict()
        elNodesSc3 = dict()
        for s_ix, s in s_list:

            #-- next mid row, nodes
            next_nodes = self._newMidRow(model, s-0.5*dir_mesh_ds)

            # mesh beam elements
            if s_ix < self.numFullRows-1:
                for i, nnode in enumerate(next_nodes):
                    elNodesMesh[nextElem] = [last_nodes[i],nnode]
                    elNodesMesh[nextElem+1] = [last_nodes[i+1],nnode]
                    nextElem+=2

            # directional mesh loop index needed to look up bolts
            if dir_mesh_ds<0:
                dir_s_ix = -2*s_ix
            else:
                dir_s_ix = 2*s_ix

            if otherside_nodes is None:
                # bolts with index 0 are treated on the second side
                dir_s_ix = None  # flag not to insert bolt or shotcrete shell
                otherside_nodes = next_nodes

            # bolts
            try:
                (bolt_angle, bolt_length, bolts_dy,
                 startpos, angleToTunnel, excess) = self.meshbolts[dir_s_ix]
            except KeyError:
                pass
            else:
                # bolt positions in the mesh
                num_bolt_rings = int(  # number of rings along the tunnel
                    (self.length-startpos-self.meshDy)/bolts_dy) + 1
                bolt_ypos = [startpos+bolts_dy*i for i in range(num_bolt_rings)]
                bolt_mesh_iy = [int(round(y/self.meshDy-0.5))
                                for y in bolt_ypos]

                for b_iy in bolt_mesh_iy:

                    # create the bolt
                    bolt_node, nextElem = self.meshbolts._createBolt(
                        model, (b_iy+0.5)*self.meshDy, s-dir_mesh_ds,
                        bolt_angle, bolt_length,
                        nextElem=nextElem,
                        secondNodeDist=excess,
                        angle_to_tunnel=angleToTunnel)

                    # connect it to the mesh
                    for node1, node2 in [
                        [beforelast_nodes[b_iy], last_nodes[b_iy]],
                        [last_nodes[b_iy], next_nodes[b_iy]],
                        [next_nodes[b_iy], last_nodes[b_iy+1]],
                        [last_nodes[b_iy+1], beforelast_nodes[b_iy]]]:
                        elNodesPlates[nextElem] = [bolt_node,node1,node2]
                        nextElem += 1

            # shotcrete shell elements, first row with tri elements
            if (s_ix==scTopRow and abs(self.shotcrete.sStart)>1e-6):
                if dir_mesh_ds>0:
                    for i in xrange(0, len(next_nodes)):
                        elNodesSc3[nextElem] = [
                            last_nodes[i], next_nodes[i], last_nodes[i+1]]
                        nextElem+=1
                else:
                    for i in xrange(0, len(next_nodes)):
                        elNodesSc3[nextElem] = [
                            next_nodes[i], last_nodes[i], last_nodes[i+1]]
                        nextElem+=1

            # shotcrete shell elements, row with quad elements
            # insert top row on second call with beforelast_nodes!=None
            if (s_ix<scLastFullRow-1 and (s_ix>scTopRow or (
                    s_ix==scTopRow==0 and beforelast_nodes is not None))):
                for i in xrange(len(next_nodes)):
                    if dir_mesh_ds>0:
                        elNodesSc4[nextElem] = [
                            beforelast_nodes[i], last_nodes[i],
                            next_nodes[i], last_nodes[i+1]]
                    else:
                        elNodesSc4[nextElem] = [
                            next_nodes[i], last_nodes[i],
                            beforelast_nodes[i], last_nodes[i+1]]
                    nextElem+=1

            # swap over to the next row of nodes
            beforelast_nodes = last_nodes
            last_nodes = next_nodes

            #-- next full row
            next_nodes = self._newFullRow(model, s)

            # mesh beam elements
            if s_ix < self.numFullRows-1:
                for i, lnode in enumerate(last_nodes):
                    elNodesMesh[nextElem] = [lnode, next_nodes[i]]
                    elNodesMesh[nextElem+1] = [lnode, next_nodes[i+1]]
                    nextElem+=2

            # directional mesh loop index needed to look up bolts
            if dir_mesh_ds<0:
                dir_s_ix = -2*s_ix-1
            else:
                dir_s_ix = 2*s_ix+1

            # bolts
            try:
                (bolt_angle, bolt_length, bolts_dy,
                 startpos, angleToTunnel, excess) = self.meshbolts[dir_s_ix]
            except KeyError:
                pass
            else:
                # bolt positions in the mesh
                num_bolt_rings = int(  # number of rings along the tunnel
                    (self.length-startpos-2.0*self.meshDy)/bolts_dy) + 1
                bolt_ypos = [startpos+bolts_dy*i for i in range(num_bolt_rings)]
                bolt_mesh_iy = [int(round(y/self.meshDy))
                                for y in bolt_ypos]

                for b_iy in bolt_mesh_iy:

                    # create the bolt
                    bolt_node, nextElem = self.meshbolts._createBolt(
                        model, b_iy*self.meshDy, s-0.5*dir_mesh_ds,
                        bolt_angle, bolt_length,
                        nextElem=nextElem,
                        secondNodeDist=excess,
                        angle_to_tunnel=angleToTunnel)

                    # connect it to the mesh
                    # note: this bolt might be on the border and only connect
                    # to three nodes
                    try:
                        node1 = beforelast_nodes[b_iy]
                    except IndexError:
                        node1 = None
                    for nodelist, offset in [
                        [last_nodes, -1],
                        [next_nodes, 0],
                        [last_nodes, 0],
                        [beforelast_nodes, 0]]:
                        try:
                            node2=nodelist[b_iy+offset]
                        except IndexError:
                            node2=None
                        if (node1 is not None) and (node2 is not None):
                            elNodesPlates[nextElem] = [bolt_node,node1,node2]
                            nextElem += 1
                        node1 = node2

            # shotcrete shell elements, row with quad elements
            if (s_ix<scLastFullRow-1 and s_ix>=scTopRow):
                if dir_mesh_ds>0:
                    elNodesSc3[nextElem] = [
                        beforelast_nodes[0], next_nodes[0], last_nodes[0]]
                    nextElem+=1
                    for i in xrange(1, len(next_nodes)-1):
                        elNodesSc4[nextElem] = [
                            beforelast_nodes[i], last_nodes[i-1],
                            next_nodes[i], last_nodes[i]]
                        nextElem+=1
                    elNodesSc3[nextElem] = [
                        beforelast_nodes[-1], last_nodes[-1], next_nodes[-1]]
                    nextElem+=1
                else:
                    elNodesSc3[nextElem] = [
                        next_nodes[0], beforelast_nodes[0], last_nodes[0]]
                    nextElem+=1
                    for i in xrange(1, len(next_nodes)-1):
                        elNodesSc4[nextElem] = [
                            next_nodes[i], last_nodes[i-1],
                            beforelast_nodes[i], last_nodes[i]]
                        nextElem+=1
                    elNodesSc3[nextElem] = [
                        next_nodes[-1], last_nodes[-1], beforelast_nodes[-1]]
                    nextElem+=1

            # shotcrete shell elements, last row: tri elements
            if s_ix==scLastFullRow-2:
                if dir_mesh_ds>0:
                    for i in xrange(0, len(next_nodes)-1):
                        elNodesSc3[nextElem] = [
                            last_nodes[i], next_nodes[i], next_nodes[i+1]]
                        nextElem+=1
                else:
                    for i in xrange(0, len(next_nodes)-1):
                        elNodesSc3[nextElem] = [
                            next_nodes[i], last_nodes[i], next_nodes[i+1]]
                        nextElem+=1

            # swap over to the next row of nodes
            beforelast_nodes = last_nodes
            last_nodes = next_nodes

        # evaluate elNodes, really insert the elements...
        model.updateElems(elNodes=elNodesMesh, elType=self.beamElType)
        model.forceElset(self.elsetName).update(elNodesMesh.iterkeys())
        if len(elNodesPlates):
            model.updateElems(
                elNodes=elNodesPlates, elType=self.boltPlatesElType)
            model.forceElset(self.boltPlatesElsetName).update(
                elNodesPlates.iterkeys())
        if len(elNodesSc4):
            model.updateElems(elNodes=elNodesSc4,
                              elType=self.shotcrete.elTypeQuad)
            model.updateElems(elNodes=elNodesSc3,
                              elType=self.shotcrete.elTypeTri)
            elset = model.forceElset(self.shotcrete.elsetName)
            elset.update(elNodesSc4.iterkeys())
            elset.update(elNodesSc3.iterkeys())

        return (nextElem, otherside_nodes)

    def addToModel(self, model, nextElem=None):
        """Add this mesh to the given abaqus model.

        @param model: an abq_model_02.Model object
        @param nextElem: Start element numbering with number nextElem. If None
        or not specified take the highest existing number plus one.

        @Note: At the moment this function assumes that there are no nodes in
        the model so far or that any existing node numbers are sequential
        starting from 1.
        """
        nextElem = _getNextElem(model, nextElem)

        # initialize first row
        top_nodes = self._newFullRow(model, s=0)

        # pos s side
        nextElem, otherside_nodes = self._insertOneSide(
            model, top_nodes, None, self.meshDy, nextElem)

        # neg s side
        nextElem, otherside_nodes = self._insertOneSide(
            model, top_nodes, otherside_nodes, -self.meshDy, nextElem)

    def connectShotcrete(self, shotcrete):
        GSMesh.connectShotcrete(self, shotcrete)

        shotcrete.topRow = int(
            abs(shotcrete.sStart)/self.meshDy)
        shotcrete.lastFullRow = int(
            abs(shotcrete.sEnd)/self.meshDy)+1

    def getNearestMeshIndex(self, s):
        """Returns the index of the mesh loop for the specified s coordinate.
        Zero identifies the top row, pos side starts with one,
        neg side with -1."""
        return int(round(2.0*s/self.meshDs))


class GSMesh090(GSMesh):
    """Mesh with horizontal and vertical wires, aligned to the tunnel axis
    and perpenticular.

    Example:

    >>> from bae.abq_model_02 import Model
    >>> from bae.generatemesh.tunnel_01 import StraightTunnel, XSectCircBack
    >>> from bae.generatemesh.gs_mesh_01 import GSMesh090, GSBoltsOnMesh
    >>>
    >>> model = Model()
    >>>
    >>> xSect = XSectCircBack(width=4.4, height=5.6, backRadius=6)
    >>> tunnelContour = StraightTunnel(xSect, length=15, heading=45.0)
    >>>
    >>> #-- GS mesh 2
    >>> mesh2Contour = tunnelContour.getOffsetContour(offset=-0.15)
    >>> gsMesh2 = GSMesh090(mesh2Contour, meshSize=0.1, meshWidth=12)
    >>>
    >>> #-- connected bolts
    >>> bolts2 = GSBoltsOnMesh(beamMeshSize=0.3)
    >>> bolts2.addRow(1, 0, length=2.5, spacing=4.0, excess=0.15)
    >>> bolts2.addRow(-7, -70, length=2.5, spacing=4.0, excess=0.15)
    >>> gsMesh2.connectBolts(bolts2)
    >>> gsMesh2.addToModel(model)
    >>>
    >>> #-- fini
    >>> model.write("GS_mesh.inp")

    @ivar meshSize: meshspacing
    @ivar numMeshIntervalls: number of mesh loops along the tunnel
    @ivar numFullRows: number of mesh loops from the top down along the contour
    """
    def __init__(self, contour, meshSize, meshWidth, elsetName="MESH",
                 numNodesBetweenCrossings=0, beamElType='B31'):
        """
        Constructor

        @param contour: the contour of the mesh of type Tunnel
        @param meshSize: mesh size of the ground support mesh
        @param meshWidth: The length of the mesh along the circumefernce of the
        tunnel cross section.
        @param numNodesBetweenCrossings: for the fem mesh, being ignored at the
        moment.
        """
        GSMesh.__init__(self, contour, meshWidth, elsetName,
                        numNodesBetweenCrossings, beamElType)
        self.meshSize = meshSize
        self.numMeshIntervalls = int(round(self.length/meshSize))
        self.numFullRows = int(self.halfMeshWidth/self.meshSize)+1

    def _newMeshRow(self, model, s, firstNode):
        """Add a new row of nodes to model at coordinate s and return the new
        node numbers"""
        new_nodes = model.insertNodes(
            ( self.contour.ys2xyzContour(y=i*self.meshSize, s=s)
              for i in range(self.numMeshIntervalls+1) ),
            firstNode=firstNode)
        return new_nodes

    def _createMeshOneSide(self, model, top_nodes, dir_mesh_ds,
                           firstNode, firstElem):
        "dir_mesh_ds>0 -> direction pos s, dir_mesh_ds<0 -> direction neg s"

        last_nodes = top_nodes
        s_list = [(i, dir_mesh_ds*i) for i in range(1, self.numFullRows)]

        elNodesMesh = list()
        elNodesPlates = list()
        for (s_ix, s) in s_list:

            next_nodes = self._newMeshRow(model, s, firstNode)

            # elements between rows
            for i, nnode in enumerate(next_nodes):
                elNodesMesh.append([last_nodes[i],nnode])

            # elements in the new row
            last_node = next_nodes[0]
            for new_node in next_nodes[1:]:
                elNodesMesh.append([last_node,new_node])
                last_node = new_node

            # bolts
            if dir_mesh_ds<0:
                dir_s_ix = -s_ix
            else:
                dir_s_ix = s_ix

            try:
                (bolt_angle, bolt_length, bolts_dy,
                 startpos, angleToTunnel, excess) = self.meshbolts[dir_s_ix]
            except KeyError:
                pass
            else:

                # bolt positions in the mesh
                num_bolt_rings = int(  # number of rings along the tunnel
                    (self.length-startpos-self.meshSize)/bolts_dy) + 1
                bolt_ypos = [startpos+bolts_dy*i for i in range(num_bolt_rings)]
                bolt_mesh_iy = [int(round(y/self.meshSize-0.5))
                                for y in bolt_ypos]

                for b_iy in bolt_mesh_iy:

                    # create the bolt
                    bolt_node, nextElem = self.meshbolts._createBolt(
                        model, (b_iy+0.5)*self.meshSize, s-0.5*dir_mesh_ds,
                        bolt_angle, bolt_length,
                        nextElem=nextElem,
                        secondNodeDist=excess,
                        angle_to_tunnel=angleToTunnel)

                    # connect it to the mesh
                    for node1, node2 in [
                        [last_nodes[b_iy+1], last_nodes[b_iy]],
                        [last_nodes[b_iy], next_nodes[b_iy]],
                        [next_nodes[b_iy], next_nodes[b_iy+1]],
                        [next_nodes[b_iy+1], last_nodes[b_iy+1]]]:
                        elNodesPlates.append([bolt_node,node1,node2])

            last_nodes = next_nodes

        # evaluate elNodes, really insert the elements...
        newElems = model.insertElems(elNodesMesh, self.beamElType, firstElem)
        model.forceElset(self.elsetName).update(newElems)
        if len(elNodesPlates):
            newElems = model.insertElems(elNodesPlates, self.boltPlatesElType, firstElem)
            model.forceElset(self.boltPlatesElsetName).update(newElems)
        return

    def addToModel(self, model, nextElem=None, firstNode=1, firstElem=None):
        """Add this shotcrete layer to the given abaqus model.

        @param model: an abq_model_02.Model object
        @param nextElem: deprecated, for compatibility only; use firstElem
        @param firstNode: Start node numbering with this number
        @param firstElem: Start element numbering with this number. Default=1.
        If not specified, nextElem will be taken.

        @Note: Nodes and elements are inserted, i.e. they don't replace
        existing ones. That implies that new node and element numbers might
        not be sequential if there are already nodes or elements in the model
        with numbers higher than firstNode, firstElem.
        """

        # for compatibility
        if firstElem is None:
            firstElem = nextElem
        if firstElem is None:
            firstElem = 1

        nextElem = _getNextElem(model, nextElem)

        # initialize first row
        top_nodes = self._newMeshRow(model, 0, firstNode)
        last_node = top_nodes[0]
        elNodes = list()
        for new_node in top_nodes[1:]:
            elNodes.append([last_node,new_node])
            last_node = new_node
        newElems = model.insertElems(elNodes, self.beamElType, firstElem)
        model.forceElset(self.elsetName).update(newElems)

        # pos s side
        nextElem = self._createMeshOneSide(
            model, top_nodes, self.meshSize, nextElem)

        # neg s side
        nextElem = self._createMeshOneSide(
            model, top_nodes, -self.meshSize, nextElem)

    def connectShotcrete(self, shotcrete):
        GSMesh.connectShotcrete(self, shotcrete)

        shotcrete.numFullRows = int(
            abs(shotcrete.sEnd - shotcrete.sStart)/self.meshSize)+1

    def getNearestMeshIndex(self, s):
        """Returns the index of the mesh loop for the specified s coordinate.
        There is no zero index, pos side starts with one, neg side with -1."""
        res = int(round((s+0.5)/self.meshSize))
        if res<1:
            res -= 1
        return res


#######################################################################

class GSBolts(object):
    """Base class for sets / rings of GS bolts.
    """

    def __init__(self, beamMeshSize, beamElType='B31',
                 boltsName='BOLTS', boltsInRockName='BOLTSINROCK'):
        self.beamMeshSize = beamMeshSize
        self.beamElType = beamElType
        self.boltsName = boltsName
        self.boltsInRockName = boltsInRockName

    def _createBolt(self, model, y_xs, s_xs, angle, length, nextElem,
                    firstnode=None, secondNodeDist=None,
                    angle_to_tunnel=90.0):
        """creates the nodes and elements for one bolt

        If firstnode is not None, take that node as first node of the bolt,
        otherwise create a new node as first node

        All bolt elements and nodes are added to elset and nset 'BOLTS', all but
        the first element are added to the elset 'BOLTSINROCK'

        @param secondNodeDist: meshsize between the first and the second node,
        (if None defaults to self.beamMeshsize)
        @param angle_to_tunnel: .. =90: perpendicular to the drive, =0:
        parallel to drive

        @Returns: a tuple (node number of the first node, nextElem at the end)

        @Note: This function assumes the tunnel to stretch exactly horizontally.
        """

        assert hasattr(self, "contour"), \
            "Meshbolts not connected to mesh. Use mesh.connectBolts()!"
        if secondNodeDist is None:
            secondNodeDist=self.beamMeshSize

        node_num = _getLastNodeNum(model)
        num_nodes = int(round(length/self.beamMeshSize))+1

        # radial postions relative to the first node
        rad_posits = [0, secondNodeDist]
        rad_dist = float(length-secondNodeDist)/(num_nodes-2)
        for i in range(num_nodes-2):
            rad_posits.append(secondNodeDist+(i+1)*rad_dist)

        # start coords
        (x,y,z) = self.contour.ys2xyzContour(y_xs,s_xs)

        # direction
        dir_xs = sin(angle*pi/180)*sin(angle_to_tunnel*pi/180)
        dir_ys = cos(angle*pi/180)*sin(angle_to_tunnel*pi/180)
        dir_along = cos(angle_to_tunnel*pi/180)
        t_dir = self.contour.orientation
        dir_vec = [t_dir[1]*dir_xs + t_dir[0]*dir_along,
                   -t_dir[0]*dir_xs + t_dir[1]*dir_along,
                   dir_ys]

        if firstnode:
            # idx_list starts at 1, firstnode already in the model
            idx_list = range(num_nodes)[1:]
        else:
            # idx_list starts at 0, firstnode will be created in the loop
            firstnode = node_num+1
            idx_list = range(num_nodes)

        elNodes = dict()
        for i in idx_list:

            # add new node
            node_num += 1
            model.nodeCoords[node_num] = vector_plus(
                [x,y,z], vector_scale(dir_vec, rad_posits[i]))
            model.forceNset(self.boltsName).add(node_num)
            if i>0:
                # create element between this node and its predecessor
                if i==1:
                    node1 = firstnode
                else:
                    model.forceElset(self.boltsInRockName).add(nextElem)
                    node1 = node_num-1

                elNodes[nextElem] = [node1,node_num]
                model.forceElset(self.boltsName).add(nextElem)
                nextElem += 1

        # really insert the elements
        model.updateElems(elNodes=elNodes, elType=self.beamElType)
        return (firstnode, nextElem)


class GSBoltsFree(GSBolts):
    """An installation of bolts not connected to mesh or shotcrete.

    Example (with mesh and connected bolts as well):

    >>> from bae.abq_model_02 import Model
    >>> from bae.generatemesh.tunnel_01 import StraightTunnel, XSectCircBack
    >>> from bae.generatemesh.gs_mesh_01 import \
    >>>     GSMesh4545, GSBoltsFree, GSBoltsOnMesh
    >>>
    >>> model = Model()
    >>>
    >>> xSect = XSectCircBack(width=4.4, height=5.6, backRadius=6)
    >>> tunnelContour = StraightTunnel(xSect, length=15, heading=45.0)
    >>>
    >>> #-- GS mesh 1
    >>> meshContour = tunnelContour.getOffsetContour(offset=-0.05)
    >>> gsMesh = GSMesh4545(meshContour,
    >>>                     approxMeshSize=0.15, meshWidth=12)
    >>>
    >>> #-- connected bolts
    >>> bolts1 = GSBoltsOnMesh(beamMeshSize=0.2)
    >>> bolts1.addRow(0, 0, length=4, spacing=5.0, excess=0.05)
    >>> bolts1.addRow(10, 20, length=4, spacing=5.0, excess=0.05)
    >>> gsMesh.connectBolts(bolts1)
    >>> gsMesh.addToModel(model)
    >>>
    >>> #-- bolts not connected to the mesh
    >>> freeBoltsContour = tunnelContour.getOffsetContour(offset=-0.1)
    >>> freeBolts = GSBoltsFree(freeBoltsContour, beamMeshSize=0.2)
    >>> # s, angle, length, ...
    >>> freeBolts.addRow(7, 110, length=4.0, spacing=2.0, excess=0.1)
    >>> freeBolts.addRow(-7, -110, length=4.0, spacing=2.0, excess=0.1)
    >>> freeBolts.addToModel(model)
    >>>
    >>> #-- fini
    >>> model.write("GS_mesh.inp")

    @ivar boltlist: a list of (s, angle, length, x_spacing) tuples, one for
    each row of bolts, see the arguments of addRow

    @Note: Might not work properly if the tunnel doesn't stretch exactly
    horizontally. Desired behaviour needs to be defined and checked.
    """

    def __init__(self, contour, *args, **kwargs):
        """Constructor

        @param contour: the contour of the mesh of type Tunnel
        @keyword beamMeshSize:
        @keyword beamElType: defaults to 'B31'
        """

        GSBolts.__init__(self, *args, **kwargs)
        assert isinstance(contour, Tunnel)
        self.contour = contour
        self.boltlist = list()

    def addRow(self, s, angle, length, spacing,
               startpos=None, angleToTunnel=90.0, excess=None):
        """
        Add one row of bolts to the installation.
        @param s: s coordinate on the cross section of the contour
        @param angle: 0 is vertical, pos. on right side (pos. x_xs)
        @param length: of the bolts
        @param spacing: along the tunnel (y-coordinate, perpendicular to cross
        section
        @param startpos: y-coordinate of the first bolt. Second will be
        (startpos+spacing), then (startpos+2*spacing), ...
        Defaults to 0.5*spacing.
        @param angleToTunnel:
        @param excess: length of the first element that doesn't belog to the
        elset BOLTSINROCK. This should be approx. the distance of self.contour
        to the rock. Defaults to self.beamMeshSize.
        """
        if startpos is None:
            startpos = 0.5*spacing
        if excess is None:
            excess=self.beamMeshSize
        self.boltlist.append((s, angle, length, spacing,
                              startpos, angleToTunnel, excess))

    def addToModel(self, model, nextElem=None):
        """Add this bolts intallation to the given abaqus model.

        @param model: an abq_model_02.Model object
        @param nextElem: Start element numbering with number nextElem. If None
        or not specified take the highest existing number plus one.

        @Note: At the moment this function assumes that there are no nodes in
        the model so far or that any existing node numbers are sequential
        starting from 1.
        """
        nextElem = _getNextElem(model, nextElem)

        # iterate through self.boltlist
        for (s, bolt_angle, bolt_length, bolts_dy, firstrow_x,
             angleToTunnel, excess) in self.boltlist:

            # bolt positions in the mesh, number in x direction
            num_bolt_rings = int(self.contour.length/bolts_dy)

            bolt_ypos = [bolts_dy*i+firstrow_x for i in range(num_bolt_rings)]

            for y in bolt_ypos:

                # create the bolt
                bolt_node, nextElem = self._createBolt(
                    model, y, s, bolt_angle, bolt_length,
                    nextElem=nextElem,
                    secondNodeDist=excess,
                    angle_to_tunnel=angleToTunnel)


class GSBoltsOnMesh(GSBolts, dict):
    """An installation of bolts connected to GS mesh.

    There is no addToModel method. Those bolts need to be added to a mesh
    and the addToModel method of the mesh also creates the bolts.

    @Note: Might not work properly if the tunnel doesn't stretch exactly
    horizontally. Desired behaviour needs to be defined and checked.
    """

    def __init__(self, *args, **kwargs):
        """Constructor

        @keyword beamMeshSize:
        @keyword beamElType: defaults to 'B31'
        @keyword boltsName: elset name for all bolt elements, defaults to 'BOLTS'
        @keyword boltsInRockName: elset name for all bolt elements but the
        innermost, defaults to 'BOLTSINROCK'.
        @keyword boltPlatesName: elset name for the bolt plates, defaults to
        'BOLTPLATES' Must be given as keyword argument.
        """
        try:
            self.boltPlatesName = kwargs['boltPlatesName']
            del kwargs['boltPlatesName']
        except KeyError:
            self.boltPlatesName = 'BOLTPLATES'

        GSBolts.__init__(self, *args, **kwargs)
        self.contour = None
        self.connectedMesh = None

    def addRow(self, sIndex, angle, length, spacing,
               startpos=None, angleToTunnel=90.0, excess=None):
        """
        Add one row of bolts to the installation.

        @param sIndex: index of the mesh square in s direction into which
        this bolts are to inserted. s_index=1 means bolt at s=meshSize/2 on the
        right side, s_index=10 means bolt at s=meshSize*9.5, s_index=-1 means
        bolt at s=-meshSize/2 on the left side.
        @param angle: 0 is vertical, pos. on right side (pos. x_xs)
        @param length: of the bolts
        @param spacing: along the tunnel (perpendicular to cross section)
        @param startpos: y-coordinate of the first bolt. Second will be
        (startpos+spacing), then (startpos+2*spacing), ...
        Defaults to 0.5*spacing.
        @param angleToTunnel: Angle between the bolt and the tunnel axis. The
        default 90 deg means bolts perpendicular to tunnel axis.
        @param excess: Length of the first element that doesn't belog to the
        elset BOLTSINROCK. This should be approx. the distance of self.contour
        to the rock. Defaults to self.beamMeshSize.
        """
        if startpos is None:
            startpos = 0.5*spacing
        if excess is None:
            excess=self.beamMeshSize
        self[sIndex] = (angle, length, spacing,
                        startpos, angleToTunnel, excess)

    def addOrthoPattern(self,
                        yStart, yWidth, sStart, sWidth, sMax,
                        anglesFunction, boltLength, excess=None):
        """
        Add an orthogonal pattern of bolts to the GS.

        You have to have connected the mesh to this bolts instance beforehand:
         >>> bolts = GSBoltsOnMesh(...)
         >>> gsMesh.connectBolts(bolts)
         >>> bolts.addOrthoPattern(...)

        @param anglesFunction: A function being supplied with the current
        s-coordinate of the connection point to the mesh returning a tuple
        of two angles for the bolts orientation: First item is the angle in the
        x-y-plane counting from zero at the (vertical) y-axis towards the
        (horizontal, pointing right) x-axis. Second item is the angle towards
        the tunnel long axis x (which will be 90 deg in most cases). Both
        angles are to be supplied in degrees.
        @param excess: Length of the first element that doesn't belog to the
        elset BOLTSINROCK. This should be approx. the distance of self.contour
        to the rock. Defaults to self.beamMeshSize.

        @returns: A tuple of two lists: The s-coordinates and the mesh loop
        index list.
        """
        assert (self.connectedMesh is not None
                and self.contour is not None), \
                "call mesh.connectBolts(bolts) before bolts.addOrthoPattern()"

        # pos side bolts
        numBoltRowsPos = int((sMax - sStart) / sWidth) + 1
        boltsSCoords = [sWidth*i+sStart for i in range(numBoltRowsPos)]

        # neg side bolts
        numBoltRowsNeg = int((sMax + sStart) / sWidth)
        boltsSCoords.extend(sStart-sWidth*(i+1) for i in range(numBoltRowsNeg))

        # prepare list of indexes, angles, add bolt rows
        boltsSIndexList = [self.connectedMesh.getNearestMeshIndex(s)
                           for s in boltsSCoords]
        boltsAnglesList = [anglesFunction(s) for s in boltsSCoords]
        for bi, (angle,angleToTunnel) in zip(boltsSIndexList, boltsAnglesList):
            self.addRow(bi, angle, length=boltLength,
                        spacing=yWidth, startpos=yStart, excess=excess)

        return (boltsSCoords, boltsSIndexList)


#######################################################################

class GSShotCreteShells(object):
    """Base class for shotcrete (fibrecrete) represented by shell elements.
    """
    pass


class GSscsFree(GSShotCreteShells, GSMesh090):
    """Shotcrete (fibrecrete) represented by shell elements not connected to a
    GS mesh.

    Example:

    >>> from bae.abq_model_02 import Model
    >>> from bae.generatemesh.tunnel_01 import StraightTunnel, XSectCircBack
    >>> from bae.generatemesh.gs_mesh_01 import GSscsFree
    >>>
    >>> model = Model()
    >>>
    >>> xSect = XSectCircBack(width=4.4, height=5.6, backRadius=6)
    >>> tunnelContour = StraightTunnel(xSect, length=15, heading=45.0)
    >>>
    >>> #-- Shotcrete layer not connected to anything
    >>> scfrContour = tunnelContour.getOffsetContour(offset=-0.30)
    >>> shotcreteFree = GSscsFree(
    >>>     scfrContour, sStart=1.5, sEnd=4.5,
    >>>     meshSize=0.5, elsetName="SHOTCRETE_FREE")
    >>> shotcreteFree.addToModel(model)
    >>>
    >>> model.write("GS_mesh.inp")

    It should be possible to connect bolts to this shotcrete-layer just like
    You would do with GSMesh090. Not tested, though.

    >>> shotcreteFree.connectBolts(meshbolts)

    @ivar sectionStartNodes: Two lists for the right (s>=0) and left side (s<0)
    of this shotcrete section. Lists the first column of nodes (e.g. the y_min
    side) from top down (in ascending s-coordinate order). Can be queried for
    connecting purposes. Not implemented yet.
    @ivar sectionEndNodes: Two lists similar to sectionStartNodes for the last
    column (e.g. y_max side).
    """
    def __init__(self, contour,
                 sStart, sEnd, meshSize,
                 elsetName='SHOTCRETE',
                 elType='S4R',
                 wallConnectorLength=None,
                 wallConnectorElsetName="SHOTCRETE_CONNECTOR"):
        """
        Constructor:

        @param contour: defines the shells mid-surface
        @param sStart: s coordinate of the top end of the shotcrete (s=0 is the
        top of the tunnel). Note that the shotcrete layer is symmetrical left
        and right side of the tunnel.
        @param sEnd: s coordinate of the bottom of the shotcrete.
        @param meshSize: mesh size along the contour and along the tunnel
        @param elsetName: all new shell elements go in this elset
        @param elType: shell element type
        @param wallConnectorLength: If not None then each node of the shotcrete
        shells gets a connector element of type JOIN and a further truss
        element reaching perpendicular to the wall. The truss element should be
        completely inside in the rock and is intended to be embedded.
        wallConnectorLength states the length of the connector element, i.e.
        should be a bit more than then distance of the shotcrete midsurface to
        the wall. This length must not be zero in the current implementation,
        it is used as the length for the truss as well.

        UPDATE:  ... now it's two beam elements, each with the given length

        @param wallConnectorElsetName: All connector elements go in this elset.
        The truss elements go into an elset with the same name plus "_EMBED".

        @note: needs further explanation and possibly further development as
        well ... docs might be outdated.
        """
        self.contour = contour
        self.length = float(self.contour.length)
        self.elType = elType
        self.elsetName = elsetName

        self.wallConnectorLength = wallConnectorLength
        self.wallConnectorElsetName = wallConnectorElsetName

        self.sStart = float(sStart)
        self.meshSize = float(meshSize)

        self.numFullRows = int(abs(float(sEnd)-self.sStart)/self.meshSize)+1

        self.numMeshIntervalls = int(round(self.length/self.meshSize))

        self.meshbolts = dict()

        self.sectionStartNodes = [list(), list()]
        self.sectionEndNodes = [list(), list()]

    def _newMeshRow(self, model, s, nextElem):
        """Add a new row of nodes to model at coordinate s and connector
        and truss elements from those nodes pointing out perpendicular to
        the surface.
        If self.wallConnectorLength is None, the connectors and trusses are
        not added.

        @returns: list of new node numbers, nextElem
        """
        node_num = _getLastNodeNum(model)

        yList = [i*self.meshSize for i in range(self.numMeshIntervalls+1)]
        points1 = [self.contour.ys2xyzContour(y, s)
                   for y in yList]
        new_nodes = model.insertNodes(points1, firstNode=node_num)

        if self.wallConnectorLength is not None:
            # nodes for the connector elements
            offsets = [vector_scale(self.contour.ys2normal(y, s),
                                    self.wallConnectorLength)
                       for y in yList]
            points2 = [vector_plus(point1, offs)
                       for point1, offs in izip(points1, offsets)]
            nodes2 = model.insertNodes(points2, firstNode=node_num)
            points3 = [vector_plus(point2, offs)
                       for point2, offs in izip(points2, offsets)]
            nodes3 = model.insertNodes(points3, firstNode=node_num)

            # connector elements
            newElems = model.insertElems(
                elNodes=zip(new_nodes, nodes2),
                elType="CONN3D2", firstElem=nextElem)
            nextElem = newElems[-1]+1
            model.forceElset(self.wallConnectorElsetName).update(newElems)

            # second elements (for embedding)
            newElems = model.insertElems(
                elNodes=zip(nodes2, nodes3),
                elType="T3D2", firstElem=nextElem)
            nextElem = newElems[-1]+1
            model.forceElset("%s_EMBED"%self.wallConnectorElsetName).update(newElems)

        return new_nodes, nextElem

    def _createMeshOneSide(self, model, top_nodes, dir_mesh_ds, nextElem):
        "dir_mesh_ds>0 -> direction pos s, dir_mesh_ds<0 -> direction neg s"

        last_nodes = top_nodes
        if dir_mesh_ds<0:
            sStart = -self.sStart
        else:
            sStart = self.sStart
        s_list = [(i, sStart+dir_mesh_ds*i)
                  for i in range(1, self.numFullRows)]

        elNodesSC = dict()
        elNodesPlates = dict()
        for (s_ix, s) in s_list:

            next_nodes, nextElem = self._newMeshRow(model, s, nextElem)

            if dir_mesh_ds>0:
                for i in range(len(next_nodes)-1):
                    elNodesSC[nextElem] = [last_nodes[i],next_nodes[i],
                                           next_nodes[i+1],last_nodes[i+1]]
                    nextElem += 1
            else:
                for i in range(len(next_nodes)-1):
                    elNodesSC[nextElem] = [last_nodes[i],last_nodes[i+1],
                                           next_nodes[i+1],next_nodes[i]]
                    nextElem += 1

            # bolts
            if dir_mesh_ds<0:
                dir_s_ix = -s_ix
            else:
                dir_s_ix = s_ix

            try:
                (bolt_angle, bolt_length, bolts_dy,
                 startpos, angleToTunnel, excess) = self.meshbolts[dir_s_ix]
            except KeyError:
                pass
            else:

                # bolt positions in the mesh
                num_bolt_rings = int(  # number of rings along the tunnel
                    (self.length-startpos-self.meshSize)/bolts_dy) + 1
                bolt_ypos = [startpos+bolts_dy*i for i in range(num_bolt_rings)]
                bolt_mesh_iy = [int(round(y/self.meshSize-0.5))
                                for y in bolt_ypos]

                for b_iy in bolt_mesh_iy:

                    # create the bolt
                    bolt_node, nextElem = self.meshbolts._createBolt(
                        model, (b_iy+0.5)*self.meshSize, s-0.5*dir_mesh_ds,
                        bolt_angle, bolt_length,
                        nextElem=nextElem,
                        secondNodeDist=excess,
                        angle_to_tunnel=angleToTunnel)

                    # connect it to the mesh
                    for node1, node2 in [
                        [last_nodes[b_iy+1], last_nodes[b_iy]],
                        [last_nodes[b_iy], next_nodes[b_iy]],
                        [next_nodes[b_iy], next_nodes[b_iy+1]],
                        [next_nodes[b_iy+1], last_nodes[b_iy+1]]]:
                        elNodesPlates[nextElem] = [bolt_node,node1,node2]
                        nextElem += 1

            last_nodes = next_nodes

        # evaluate elNodes, really insert the elements...
        model.updateElems(elNodes=elNodesSC, elType=self.elType)
        model.forceElset(self.elsetName).update(elNodesSC.iterkeys())
        if len(elNodesPlates):
            model.updateElems(elNodes=elNodesPlates, elType=self.boltPlatesElType)
            model.forceElset(self.boltPlatesElsetName).update(elNodesPlates.iterkeys())
        return nextElem

    def addToModel(self, model, nextElem=None):
        """Add this mesh to the given abaqus model.

        @param model: an abq_model_02.Model object
        @param nextElem: Start element numbering with number nextElem. If None
        or not specified take the highest existing number plus one.

        @Note: At the moment this function assumes that there are no nodes in
        the model so far or that any existing node numbers are sequential
        starting from 1.
        """
        nextElem = _getNextElem(model, nextElem)

        # initialize first row
        top_nodes, nextElem = self._newMeshRow(model, self.sStart, nextElem)

        # pos s side
        nextElem = self._createMeshOneSide(
            model, top_nodes, self.meshSize, nextElem)

        if abs(self.sStart) > 1E-6:
            # mesh doesn't start at the top, initialize new top row
            top_nodes, nextElem = self._newMeshRow(model,-self.sStart,nextElem)

        # neg s side
        nextElem = self._createMeshOneSide(
            model, top_nodes, -self.meshSize, nextElem)


class GSscsOnMesh(GSShotCreteShells):
    """Shotcrete (fibrecrete) represented by shell elements on top of GS mesh.
    The shotcrete shells are connected by nodes to the mesh.
    """
    def __init__(self, sStart, sEnd,
                 elsetName='SHOTCRETE',
                 elTypeQuad='S4R', elTypeTri='S3'):
        self.sStart = float(sStart)
        self.sEnd = float(sEnd)
        self.elsetName = elsetName
        self.elTypeQuad = elTypeQuad
        self.elTypeTri = elTypeTri

