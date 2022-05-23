"""dyna_model_01.py

This module is for reading(*), writing and manipulating LS-Dyna input files.

Usage
=====
  >>> from bae.dyna_model_01 import DynaModel
  >>> m = DynaModel()
  >>> help(m) for more information

(*) reading not implemented yet.
"""

__version__ = "1.03"

_version_history_ = """
Versions:
=========

0.01 : new
0.02 : added loadRemovePartSet
1.03 GP renumbered version from 0.xx to 1.xx
        added StressDistribution, StressGeostatic and
        DynaModel.updateInitStressSolidData()
"""

_todo_ = """
"""
from collections import defaultdict

from bae.mesh_01 import Mesh
# from bae.log_01 import msg
from bae.misc_01 import AutoKeyDict, groupIter

#---------------------------------------------------------
# model class to hold all data in a LS-Dyna keyword file
#---------------------------------------------------------

class DynaModel(Mesh):
    """
    Reads and writes ABAQUS input files, provides the data as member
    attributes and has some functionality to manipulate this data.

    Usage:

      >>> from bae.abq_model_02 import Model
      >>> m = Model()
      >>> m.read(inp_name)  # specify filename *with* extension .inp

    @ivar nodeCoords: dict {node number: [x,y,z]},
    @ivar elNodes: dict {element number: [node numbers]}, B{do not modify
        directly}, use L{updateElems}, L{removeElems}.
    @ivar part: dict {part id: L{PartData} }
    @ivar section: dict {section id: L{SectionSolidData} } to hold section
        definitions
    @ivar material: dict {material id: L{MaterialData} } for material
        properties
    @ivar setPart: dict {set id: set of part ids}
    @ivar table: dict {table/curve id: L{TableData} } for curves and tables
    @ivar loadStiffenPartSet: dict {part set id: curve id} to define load
        stiffen associations between part sets an load curves
    @ivar loadRemovePartSet: dict {part set id: (T1, T2)} to define removal
        of the part set between times T1 and T2
    @ivar initStressSolidData: L{InitStressSolidData}-object, basically a dict
        {element number: InitStressSolidDataValue-object}. See the convenience
        function L{updateInitStressSolidData}, for example to make use of
        L{StressGeostatic}
    """

    def __init__(self, oldmodel=None):
        """
        @param oldmodel: can be a L{abq_model_02.Model} or L{mesh_01.Mesh}
           or another DynaModel object in order to work as a copy constructor.

           Note: No deep copy is performed. The values of the nodeCoords,
           elNodes, parts dictionaries are identical to the ones in oldmodel,
           i.e. node coordinates-lists are identical as are the element
           conectivity lists.
        """

        # data members, they are always there; from base class
        self.nodeCoords = AutoKeyDict(factory=list)
        self.elNodes = AutoKeyDict(factory=list)
        self.elShape = dict()
        self.shapeEl = defaultdict(set)
        
        # data members, they are always there; special to DynaModel
        self.part = AutoKeyDict(factory=PartData)
        self.section = AutoKeyDict(factory=SectionSolidData)
        self.material = AutoKeyDict(factory=MaterialData)
        self.setPart = AutoKeyDict(factory=set)
        self.table = AutoKeyDict(factory=TableData)
        self.loadStiffenPartSet = dict()
        self.loadRemovePartSet = dict()
        self.initStressSolidData = InitStressSolidData()

        # self.elShape = internal.FakeElShape(self.elType, self.elTypeToShapeSep)
        # self.shapeEl = internal.FakeShapeEl(
        #     self.typeEl, self.elShapeToTypes, self.elTypeToShapeSep)

        if isinstance(oldmodel, Mesh):
            self.nodeCoords.update(oldmodel.nodeCoords)
            self.elNodes.update(oldmodel.elNodes)
            self.elShape.update(oldmodel.elShape.iteritems())
            self.shapeEl.update(oldmodel.shapeEl.iteritems())

            if isinstance(oldmodel, DynaModel):
                self.part.update(oldmodel.part)
                self.section.update(oldmodel.section)
                self.material.update(oldmodel.material)
            else:
                # oldmodel is a mesh_01.Mesh or abq_model_02.Model object

                # add dummy material with matId=1
                matId = self.material.append(
                    # type,    rho,    E,   nu
                    "ELASTIC", 2700.0, 1E9, 0.2)

                # one part for each element shape
                for shape, elems in oldmodel.shapeEl.iteritems():

                    # correct node ordering
                    if shape=="TET_Q":
                        # indexes in elNodes: Abq: 7, 8, 9 -> dyna: 9, 7, 8
                        for elem in elems:
                            nodes = self.elNodes[elem]
                            self.elNodes[elem] = [
                                nodes[i] for i in [0,1,2,3,4,5,6,9,7,8]]
                    elif shape=="WEDGE_L":
                        # indexes in elNodes:
                        # Abq: 0,1,2,3,4,5 -> dyna: 0,4,1,3,5,2
                        for elem in elems:
                            nodes = self.elNodes[elem]
                            self.elNodes[elem] = [
                                nodes[i] for i in [0,4,1,3,5,2]]

                    # create section and part
                    secId = self.section.append(
                        *self._shapeToSectionSolidData[shape])
                    self.part.append(
                        name="part_%s"%shape,
                        secId=secId, matId=matId, elems=set(elems))

    _shapeToSectionSolidData = {
        # ... : (elform, aet)
        "TET_L" : (10, 0),
        "TET_Q" : (17, 0), # node ordering changed!
        "WEDGE_L" : (21, 0), # node ordering changed!
        "HEX_L" : (1, 0),
        }

    def checkParts(self):
        """
        Check:
         - If all elements listed in parts actually exist (are defined).
         - If all defined parts are referred to by at least one element.
         - That all parts are disjoint.
         - That all parts listed in any set_part actually exist
        """
        emptyParts = [ pid
                       for pid, partData in self.part.iteritems()
                       if len(partData.elems)==0 ]
        assert not(emptyParts), \
            ("There are %d parts with no elements:\n%s."
             % (len(emptyParts), ",\n".join(
             "part %s (pid %d) is empty" % (self.part[pid].name, pid)
             for pid in emptyParts)))
        notDefinedElems = [
            elem
            for partData in self.part.itervalues()
            for elem in partData.elems
            if elem not in self.elNodes ]
        assert len(notDefinedElems)==0, \
            ("%d elements listed in parts are not defined."
             % len(notDefinedElems))

        parts = [self.part[pid] for pid in sorted(self.part)]
        intersections = [
            (i, part1, part2)
            for i, part1 in enumerate(parts)
            for part2 in parts[i+1:]
            if part1.elems.intersection(part2.elems) ]
        assert not(intersections),\
            ("%d parts intersect:\n%s." % ",\n".join(
                "part %s (pid %d) and part %s intersect"
                % (part1.name, sorted(self.part)[i], part2.name)
                for (i, part1, part2) in intersections))

        notDefinedParts = [
            pid
            for setPart in self.setPart.itervalues()
            for pid in setPart
            if pid not in self.part]
        assert len(notDefinedParts)==0, \
            ("%d parts listed in set_parts are not defined."
             % len(notDefinedParts))
        
    # #### to be implemented!
    # def removeElems(self, elems, removeNodes=False):
    #     """removes an element from the mesh
    #     updates elNodes, elShape, shapeEl dicts and all elsets

    #     @param elems: A list or a set of element numbers (also accepts a dict,
    #     in that case takes its keys.)

    #     @param removeNodes: if True remove unused nodes. This is not very
    #     efficient if many elements are beeing deleted, in this case rather
    #     clean up after all elements have been removed using
    #     getNodesFromElset():
    #       >>> mesh.removeElems(elems=myBigElset)
    #       >>> allNodes = mesh.getNodesFromElset()
    #       >>> removeNodes = set(mesh.nodeCoords).difference(allNodes)
    #       >>> for node in removeNodes:
    #       ...    del mesh.nodeCoords[node] 

    #     @note: Silently ignores missing elements.
    #     @note: It is save to do the following (it is not only save but the
    #     recommended way to accomplish those tasks):
    #       >>> mesh.removeElems(elems=mesh.elNodes)  # remove all elements
    #       >>> mesh.removeElems(elems=mesh.shapeEl['TRI_L']) # remove one type
    #     """

    #     # special treatment if you want to remove all, wouldn't work otherwise
    #     if elems is self.elNodes:
    #         self.elNodes.clear()
    #         self.shapeEl.clear()
    #         self.elShape.clear()
    #         return

    #     # node list of element nodes
    #     if removeNodes:
    #         removeNodesSet = set()
    #         for elNum in elems:
    #             try:
    #                 removeNodesSet.update(self.elNodes[elNum])
    #             except KeyError:
    #                 # elNum not in self.elNodes
    #                 pass

    #     # remove elements from self.elNodes and self.elShape
    #     for elNum in elems:
    #         try:
    #             del self.elShape[elNum]
    #         except KeyError:
    #             pass

    #         try:
    #             del self.elNodes[elNum]
    #         except KeyError:
    #             pass

    #     # remove elements from self.shapeEl
    #     # use a copy (self.shapeEl.items()) in the loop
    #     for thisshape, els in self.shapeEl.items():
    #         els.difference_update(elems)
    #         if len(els)==0:
    #             del self.shapeEl[thisshape]

    #     # remove nodes
    #     if removeNodes:
    #         for nodes in self.elNodes.itervalues():
    #             removeNodesSet.difference_update(nodes)
    #             if len(removeNodesSet)==0:
    #                 break
    #         for node in removeNodesSet:
    #             del self.nodeCoords[node]

    def splitParts(self, elems, name,
                   oldPids=None,
                   section=None, material=None):
        """
        Split off elements from existing parts. New parts are being created
        with optional new section and/or material ids.

        @param elems: iterable (list, set) of element ids
        @param name: heading (=alphanumeric part name) of the new parts.
           May contain the string {oldName} (including the braces) to be
           replaced by the name/heading of the part from which the new part
           is being split off.
        @param oldPids: iterable (list, set) of part ids. Only consider parts
           contained in this list for the split. If not given consider all
           parts.
        @param section: Section id for the new parts. If not given then retain
           the section id of the old part.
        @param material: Material id for the new parts. If not given then keep
           the material id of the old part.

        @returns: a set containing all parts that have been split off. Those
           will generally be the newly created parts. But old parts may be
           included as well if these parts are completely contained in elems.
        """

        elems = set(elems)
        newparts = set()

        def updatePart(partData, section, material):
            if section:
                partData.secId = section
            if material:
                partData.matId = material

        if oldPids:
            oldPids = sorted(oldPids)
        else:
            oldPids = sorted(self.part)

        for pid in oldPids:
            partData = self.part[pid]
            newElems = partData.elems.intersection(elems)
            if not newElems:
                continue

            newName = name.format(oldName=partData.name)

            # remove elems from old part
            partData.elems.difference_update(newElems)
            if not(len(partData.elems)):
                # rest of old is empty, so replace old part by new
                newparts.add(pid)
                updatePart(partData, section, material)
                partData.elems = newElems
                partData.name = newName
            else:
                # create new part
                old = self.part[pid]
                newpid = self.part.append(
                    newName, old.secId, old.matId, newElems)
                newparts.add(newpid)
                updatePart(self.part[newpid], section, material)

                # update set_part
                for parts in self.setPart.itervalues():
                    if pid in parts:
                        parts.add(newpid)

        return newparts

    def elemsOfSetPart(self, setPart):
        """return a set of all elements of a given set_part"""
        elems = set(elem
                    for part in self.setPart[setPart]
                    for elem in self.part[part].elems)
        return elems

    def updateInitStressSolidData(self, elset, **kwargs):
        """Assign a L{InitStressSolidDataValue} object to all elements in
        elset (overwriting existing assignments).

        For more arguments not listed below see
        L{InitStressSolidDataValue.__init__}.

        @param elset: Iterable of element numbers

        @kwarg numbers: a tuple of (nint, nhisv, nthint, nthhsv)

          nint is the number of integration points, nhisv the number of history
          varibles per integration point, nthint the number of thermal
          integration points per element, nthhsv the number of thermal history
          variables per integration point. The last two must be zero currently.

        @kwarg initStress: Might be a list of nint 6-tuple of stress
           components: [sig11, sig22, sig33, sig12, sig23, sig13]. Or an
           object of a subclass of L{StressDistribution}, e.g.
           L{StressGeostatic}.
        """

        initStress = kwargs["initStress"]
        if isinstance(initStress, StressDistribution):
            nint = kwargs["numbers"][0]
            elGaussPoints = self.getElGaussPoints(elset)
            if not all(len(gps)==nint for gps in elGaussPoints.itervalues()):
                elShapes = sorted(set(self.elShape[eid] for eid in elset))
                raise ValueError(
                    "Number of Gauss-points for elShapes %s does not"
                    " correspond to nint=%d" % (", ".join(elShapes), nint))
            for eid in elset:
                kwargs["initStress"] = [
                    initStress.getStress(pt) for pt in elGaussPoints[eid]]
                self.initStressSolidData[eid] = InitStressSolidDataValue(
                    **kwargs)
        else:
            values = InitStressSolidDataValue(**kwargs)
            self.initStressSolidData.update( (e, values) for e in elset )


    @staticmethod
    def floatTo10Char(x):
        if x>=1e100:
            # exponential with 3 digit exponent: 1.234e+123
            return "%.4g" % x
        elif x>9999999999:
            # exponential with 2 digit exponent: 1.2345e+12
            return "%.5g" % x
        elif x>999999999:
            # (rounded) integer, no decimal point: 1234567890
            return "%.10g" % x  # "%.0f" would yield the same
        elif x>=1:
            # float with decimal point: 123.567891
            return "%.9g" % x
        elif x>0.001:
            # float starting "0.": 0.99999999 .. 0.00100001
            return "%.8g" % x   # note: "%.8f" causes trailing zeros
        elif x>=1e-99:
            # exponential with 2 digit exponent: 1.2345e-12
            return "%.5g" % x
        elif x>=0:
            # exponential with 3 digit exponent: 1.234e-123
            return "%.4g" % x
        elif x>=-9.99e-100:
            # exponential with 3 digit exponent: -1.23e-123
            return "%.3g" % x
        elif x>=-0.001:
            # exponential with 2 digit exponent: -1.234e-12
            return "%.4g" % x
        elif x>=-1:
            # float starting "-0.": -0.2345678 .. -0.0010001
            return "%.7g" % x
        elif x>=-99999999:
            # float with decimal point: -123.45678
            return "%.8g" % x
        elif x>=-999999999:
            # (rounded) integer, no decimal point: -123456789
            return "%.0f" % x
        elif x>-1e100:
            # exponential with 2 digit exponent: -1.234e+12
            return "%.4g" % x
        else:
            # exponential with 3 digit exponent: -1.23e+123
            return "%.3g" % x

    @staticmethod
    def floatTo20Char(x):
        # currently implemented: cheap version
        # because 14 or 13 significant digits seem sufficient
        if x>0:
            return "%.14g" % x
        else:
            return "%.13g" % x

    def write(self, outputFile):
        """Writes the model data to an LS-Dyna keyword file.

        Examples, Usage
        ===============

        >>> m = DynaModel()
        >>> m.write("mymodel.k")

        >>> fo = open("mymodel.k", "w")
        >>> fo.write("** Some extra things to write to this file.")
        >>> m.write(fo)
        >>> fo.close()

        Interface
        =========

        @param outputFile: May be an open file object or a filename (a string).

         If you supply an open file object it will be left open by this method.
         There will be no header (seems to be optional anyway) and no trailing
         *END statement. You have to add an *END manually afterwards.

         If you supply a string as filename a new file will be created silently
         overwriting an existing file with the same name. A header and an *END
         command will be added to the file. The file will be closed at the end.
        """
        if type(outputFile) == str:
            outputFile = open(outputFile, "w")
            close_file = True
        else:
            assert hasattr(outputFile, "write")
            close_file = False

        # write header
        if close_file:
            import sys, time
            outputFile.write(
                "$ LS-DYNA Keyword file created by %s\n"
                "$ Created on %s\n"
                % (sys.argv[0], time.asctime()))

        # write node definitions
        outputFile.write("*NODE\n"
                         "$ nid,x,y,z,tc,rc\n")
        for nid in sorted(self.nodeCoords):
            coords = ",".join(
                self.floatTo20Char(x) for x in self.nodeCoords[nid])
            outputFile.write("%d,%s\n" % (nid, coords))

        # write element definitions
        outputFile.write("*ELEMENT_SOLID (ten nodes format)\n"
                         "$ first line: eid, pid; second line: node numbers\n")
        for pid in sorted(self.part):
            for elem in sorted(self.part[pid].elems):
                # first line: eid, pid
                outputFile.write("%d,%d\n" % (elem, pid))
                # second line: node numbers
                outputFile.write(",".join(
                    str(node) for node in self.elNodes[elem]))
                outputFile.write("\n")

        # write section definitions
        outputFile.write("*SECTION_SOLID\n"
                         "$ secId, elform, aet\n")
        for secId in sorted(self.section):
            secData = self.section[secId]
            outputFile.write(
                "%d,%d,%d\n" % (secId, secData.elform, secData.aet))

        # write material definitions
        for matId in sorted(self.material):
            matData = self.material[matId]
            matDataList = [str(matId)] + [
                self.floatTo10Char(x) for x in matData.parameter]
            matParams = "\n".join(
                ",".join(line)
                for line in groupIter(matDataList, maxcounts=8))
            outputFile.write("*MAT_%s\n%s\n"
                             % (matData.type_, matParams))

        # write part definitions
        for pid in sorted(self.part):
            partData = self.part[pid]
            outputFile.write(
                "*PART\n"
                "$ pid, secId, matId, eosid, hgid, grav, adpopt, tmid\n")
            outputFile.write("%s\n%s\n" % (partData.name, ",".join(
                str(x) for x in [pid, partData.secId, partData.matId,
                                 0,0,0,0,0])))

        # write set_part definitions
        for setId in sorted(self.setPart):
            parts = self.setPart[setId]
            outputFile.write(
                "*SET_PART\n"
                "$ setId, DA1, DA2, DA3, DA4, SOLVER, parts on next lines\n")
            outputFile.write(
                "%d,0.0,0.0,0.0,0.0,MECH\n%s\n" % (setId, "\n".join(
                    ",".join(line)
                    for line in groupIter(parts, maxcounts=8, format="%d") )))

        # write table/curve definitions
        for tableId in sorted(self.table):
            data = self.table[tableId]
            if data.type_=="curve":
                outputFile.write(
                    "*DEFINE_CURVE\n"
                    "$ id, opt, sfx, sfy, offx, offy, type\n")
                outputFile.write("%d,0,1.0,1.0,0.0,0.0,0\n" % tableId)
                outputFile.write("$ x, y\n")
                for line in data.tab:
                    outputFile.write(
                        "%s,%s\n" % tuple(map(self.floatTo20Char, line)))
            else:
                raise NotImplementedError(
                    "Table type %s not implemented" % data.type_)

        # write LOAD_STIFFEN_PART_SET
        if len(self.loadStiffenPartSet):
            outputFile.write(
                "*LOAD_STIFFEN_PART_SET\n"
                "$ partSetId, load curve id\n")
            for psid in sorted(self.loadStiffenPartSet):
                outputFile.write(
                    "%d,%d\n" % (psid, self.loadStiffenPartSet[psid]))

        # write LOAD_REMOVE_PART_SET
        if len(self.loadRemovePartSet):
            outputFile.write(
                "*LOAD_REMOVE_PART_SET\n"
                "$ partSetId, start time, end time\n")
            for psid in sorted(self.loadRemovePartSet):
                t0, t1 = self.loadRemovePartSet[psid]
                outputFile.write(
                    "%d,%s,%s\n"
                    % (psid, self.floatTo10Char(t0), self.floatTo10Char(t1)))

        # write INITIAL_STRESS_SOLID
        if len(self.initStressSolidData):
            outputFile.write(
                "*INITIAL_STRESS_SOLID\n"
                "$# eid, nint, nhisv, large, iveflg, ialegp, nthint, nthhsv\n"
                "$# sig11, sig22, sig33, sig12, sig23\n")
            for eid in sorted(self.initStressSolidData):
                dataVal = self.initStressSolidData[eid]
                nint, nhisv, nthint, nthhsv = dataVal.numbers
                if nint!=1:
                    raise NotImplementedError(
                        "Only nint=1 implemented. Don't know what to do with"
                        " nint=%s. Element %d." % (nint, eid))
                if nthint!=0 or nthhsv!=0:
                    raise NotImplementedError(
                        "Only nthint=0 or nthhsv=0 implemented. Don't know"
                        " what to do with nthint=%s and nthhsv=%s. Element %d."
                        % (nthint, nthhsv, eid))
                values = (dataVal.initStress[0]
                          + [dataVal.initEps[0]]
                          + dataVal.initHisv[0])
                outputFile.write(
                    "%s\n%s\n" % (
                        # eid, nint,nhisv, large,iveflg,ialegp, nthint,nthhsv
                        ",".join("%d"%x for x in (
                            eid, nint, nhisv, 1, 0, 0, nthint, nthhsv)),
                        "\n".join(
                            ",".join(line) for line in groupIter(
                                map(self.floatTo20Char, values), maxcounts=5))
                    ))

        if close_file:
            # write *END
            outputFile.write("*END\n")

            outputFile.close()
        return

#---------------------------------------------------------
# classes for specific data items
#---------------------------------------------------------


class PartData(object):
    """Data container for a part

    @ivar name: part name string (heading)
    @ivar secId: integer section id
    @ivar matId: integer material id
    @ivar elems: set of element ids
    """

    def __init__(self, name="", secId=0, matId=0, elems=None):
        self.name = name
        self.secId = secId
        self.matId = matId
        if elems:
            self.elems = set(elems)
        else:
            self.elems = set()

    def __str__(self):
        return ("PartData for %s: secId %d, matId %d, #elems %d"
                % (self.name, self.secId, self.matId, len(self.elems)))


class SectionSolidData(object):
    """Data container for a section definition (element type)

    @ivar elform: element type: 0 = 8-node brick
    @ivar aet: usually zero.
    """

    def __init__(self, elform=0, aet=0):
        self.elform = elform
        self.aet = aet

    def __str__(self):
        return ("SectionSolidData: elform %d, aet (ambient element type) %d"
                % (self.elform, self.aet))


class MaterialData(object):
    """Data container for the material properties of one material

    @ivar type_: Material type like "ELASTIC"
    @ivar parameter: list of parameter values (floats)
    """
    
    def __init__(self, type_, *parameter):
        if not(isinstance(type_, basestring)):
            raise ValueError("Material type must be a string but got %s"
                             % type_)
        self.type_ = type_
        self.parameter = list(parameter)


class TableData(object):
    """Data for curves and tables.

    @ivar type_: "curve"
    @ivar tab: list of [x, y]-tuples if type_="curve"

    @Note: Currently only type_="curve" implemented.
    @Note: object holds a reference of the data passed to __init__(), no copy
       is made.
    """
    def __init__(self, type_, tab):
        self.type_ = type_.lower()
        if self.type_=="curve":
            self.tab = tab
        else:
            raise NotImplementedError("Table type %s not implemented" % type_)

class InitStressSolidDataValue(object):
    """Data for a particular set of *INITIAL_STRESS_SOLID data values
    """
    def __init__(self, numbers, initStress, initEps, initHisv):
        """
        @param numbers: a tuple of (nint, nhisv, nthint, nthhsv)

          nint is the number of integration points, nhisv the number of history
          varibles per integration point, nthint the number of thermal
          integration points per element, nthhsv the number of thermal history
          variables per integration point. The last two must be zero currently.

        @param initStress: list of nint 6-tuples of stress values or an
           object of a subclass of L{StressDistribution}, e.g.
           L{StressGeostatic}.
        @param initEps: list of nint floats, the effective plastic strains
        @param initHisv: list of nint tuple of nhisv floats to initialize the
          history variables
        """
        self.numbers = numbers
        self.initStress = initStress
        self.initEps = initEps
        self.initHisv = initHisv

class InitStressSolidData(dict):
    """All the data for the *INITIAL_STRESS_SOLID command
    basically a dict {element number: InitStressSolidDataValue-object}.
    Maybe with some extra housekeeping attributes in the future...
    """
    pass


class StressDistribution(object):
    pass


class StressGeostatic(StressDistribution):
    """Object to compute geostatic stress distribution.

    Usage:
     >>> initStress = StressGeostatic(
     >>>     rho=2700, g=9.81, zeroLevel=213.0, k1=1.2, k2=1.3)
     >>> model.updateInitStressSolidData(
     >>>     ringsElset, numbers=(nint, nhisv, nthint, nthhsv),
     >>>     initStress=initStress, initEps=0.0, 
    """
    def __init__(self, **kwargs):
        """Constructor.

        Possible sets of arguments:
         - rho, g, zeroLevel, k1, k2

        @kwarg rho: rock density
        @kwarg g: gravity constant, usally 9.81
        @kwarg zeroLevel: z-coordinate where we have S=0
        @kwarg k1: S_11 = k1*S_33
        @kwarg k2: S_22 = k2*S_33
        """
        if all((x in kwargs) for x in ("rho", "g", "zeroLevel", "k1", "k2")):
            self.zeroLevel = kwargs["zeroLevel"]
            rhoG = kwargs["rho"] * kwargs["g"]
            k1, k2 = kwargs["k1"], kwargs["k2"]
            self.gradS = [rhoG*k1, rhoG*k2, rhoG, 0.0, 0.0, 0.0]
        else:
            raise NotImplementedError(
                "Arguments %s not sufficient for StressGeostatic."
                % (",".join(kwargs)))

    def getStress(self, point):
        """Return a list of stress components [sig11, sig22, sig33, sig12,
        sig23, sig13] for the given point.
        @param point: List of x,y,z coordinates.
        """

        zR = point[2]-self.zeroLevel
        return [s*zR for s in self.gradS]
