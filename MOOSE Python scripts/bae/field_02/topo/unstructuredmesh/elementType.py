"""Module holding the L{ElementType} class.
"""

# An ElementType-object holds all relevant (and known/determined) attributes
# of a certain element type like nodesPerElem, midSideNodes,
# gaussPtElCoords, Abaqus name, ... everything that was stored in separate
# dicts in bae.mesh_01 and bae.abq_model_02.
# Have sensible default values for everything not explicitly specified or
# have properties (i.e. functions that determine a sensible default on the
# fly in case some data is provided after initialisation of the ElementType
# object.

import numpy as np

class ElementType(object):
    """Objects identifying element types, element shapes, element connectivity.
    Always query particular attributes of the ElementType object like
    nodesPerElem. Never test for something like the abqName and try to infer
    some other attribute that you are actually interested in. I.e. a query
    elType.abqName=="C3D4" should likely be replaced by
    elType.nodesPerElem==4 and/or elType.shape=="TET".

    Initialization should always be done with all information available, all
    other attributes will be derived from that, i.e. they get default values
    by the __init__ method.

    @ivar abqName: Abaqus element type label like C3D4, S3
    @ivar vtkCellType: VTK element type number.
        Linear types: 3=line, 5=tri, 9=quad, 10=tet, 12=brick, 13=wedge.
        Quad types: 21=line, 22=tri, 23=quad, 24=tet, 25=hex
    @ivar shape: "LINE", "TRI", "QUAD", "TET", "WEDGE", "HEX"

    @ivar order: 1 or 2
    @ivar nodesPerElem: an int

    @ivar nodesOnFace: Indexes for the nodes on each face. For linear elements
       with only corner nodes this should be identical to faceNodes. For
       quadratic elements nodesOnFace also contains the mid side nodes of each
       face. nodesOnFace[i,j] is the index in elNodes of the j-th node on the
       i-th face.
       The index ordering of each face is such that it fits an element type
       that can represent the face with normal pointing inwardly.

    @ivar faceNodes: Node indexes for each face. The face is represented by its
       corner nodes only, for quadratic elements the mid-side nodes are left
       out.
       faceNodes[i,j] is the index in elNodes of the j-th node on the i-th
       face. The index ordering of each face is such that each faces normal
       points inwardly.
    """

    #--- conversion to Abaqus element type
    # note: the first Abaqus elType in the list is the default
    shapeToAbqName = {
        "LINE": ("B31","T3D2","CONN3D2"),
        "TRI": ("S3","S3R","DS3","M3D3","CPE3",
                "CPS6M","CPE6","STRI65","DS6"),
        "QUAD": ("S4","S4R","S4RS","CPS4R","CPE4","CPE4R",
                 "CPS8","CPS8R","CPE8","CPE8R","CPE8RP"),
        "TET": ("C3D4","C3D4H","C3D4P","C3D4T","DC3D4",
                "C3D10","C3D10M","C3D10H","C3D10MH","C3D10MP"),
        "WEDGE": ("C3D6","C3D6R","C3D6P","DC3D6","C3D15","SC6R",
                  "COH3D6","COH3D6P",),
        "HEX": ("C3D8","C3D8P","C3D8R","C3D20","C3D20R","SC8R"),
    }
    abqNameToShape = dict(
        (t, s) for s, tt in shapeToAbqName.iteritems() for t in tt)

    # note: the first Abaqus elType in the list is the default
    vtkCellTypeToAbqName = {
        3: ("B31","T3D2", "CONN3D2"),
        5: ("S3","M3D3","CPE3","S3R","DS3",),
        9: ("S4","S4R","S4RS","CPS4R","CPE4","CPE4R",),
        10: ("C3D4","C3D4H","C3D4P","C3D4T",
             "DC3D4",),
        13: ("C3D6","C3D6R","C3D6P","DC3D6","SC6R",  # linear wedge
             "COH3D6","COH3D6P",),
        12: ("C3D8","C3D8P","C3D8R","SC8R",),

        # Note that for "B32" the node ordering does not fit to the
        # getCentroid()-algorithm
        21: ("B32",),  # quad line
        22: ("DS6","CPS6M","CPE6","STRI65"),  # quadratic triangle
        # quadratic quadrilateral
        23: ("CPS8","CPS8R","CPE8","CPE8R", "CPE8RP"),
        24: ("C3D10","C3D10M","C3D10H","C3D10MH","C3D10MP"),  # quad tet
        # ??: ("C3D15",),  # quadratic wedge not implemented in vtk
        25: ("C3D20","C3D20P","C3D20R"),  # quad hex
    }
    abqNameToVtkCellType = dict(
        (at, vt) for vt, att in vtkCellTypeToAbqName.iteritems() for at in att)

    #--- addional info generally based on Abaqus element type
    abqNameToNodesPerElem = {
        # lines
        "B31": 2, "T3D2": 2, "CONN3D2": 2,

        # tris
        "S3": 3, "S3R": 3, "DS3": 3, "M3D3": 3, "CPE3": 3,
        "CPS6M": 6, "CPE6": 6, "STRI65": 6, "DS6": 6,

        # quads
        "S4": 4, "S4R": 4, "S4RS": 4, "CPS4R": 4, "CPE4": 4, "CPE4R": 4,
        "CPS8": 8, "CPS8R": 8, "CPE8": 8, "CPE8R": 8, "CPE8RP": 8,

        # tets
        "C3D4": 4,"C3D4H": 4, "C3D4P": 4, "C3D4T": 4, "DC3D4": 4,
        "C3D10": 10, "C3D10M": 10, "C3D10H": 10, "C3D10MH": 10, "C3D10MP": 10,

        # wedges
        "C3D6": 6, "C3D6R": 6, "C3D6P": 6, "DC3D6": 6, "SC6R": 6,
        "COH3D6": 6, "COH3D6P": 6,
        "C3D15": 15,

        # hex / bricks
        "C3D8": 8, "C3D8P": 8, "C3D8R": 8, "SC8R": 8,
        "C3D20": 20, "C3D20R": 20,
    }

    # order of the shape functions, 1 for (bi-/tri-) linear, 2 for quadratic
    abqNameToOrder = dict(
        (t, o) for o, tt in [
            (1, ("B31","T3D2", "CONN3D2")),
            (1, ("CPE3","S3","S3R","DS3","M3D3")),
            (1, ("CPS4R","CPE4","CPE4R","S4","S4R","S4RS")),
            (1, ("C3D4","C3D4H","C3D4P","C3D4T")),
            (1, ("DC3D4",)),
            (1, ("COH3D6","COH3D6P",
                 "C3D6","C3D6R","C3D6P","DC3D6","SC6R")),
            (1, ("C3D8","C3D8P","C3D8R","SC8R")),

            # Note that for "B32" the node ordering does not fit to the
            # getCentroid()-algorithm
            (2, ("B32",)),
            (2, ("CPS6M","CPE6","STRI65","DS6")),
            (2, ("CPS8","CPS8R","CPE8","CPE8R", "CPE8RP")),
            (2, ("C3D10","C3D10M","C3D10H","C3D10MH","C3D10MP")),
            (2, ("C3D15",)),
            (2, ("C3D20","C3D20P","C3D20R")),
        ] for t in tt)

    #: Node indexes for each face of particular element types.
    #: The index ordering of each face is such that each faces normal
    #: points inwardly.
    abqNameToNodeIdsOnFace = [
        # tri_L
        (("S3","S3R","DS3","M3D3","CPE3"), ((0,1), (1,2), (0,2))),
        # tri_Q
        (("CPS6M","CPE6","STRI65","DS6"), ((0,1,3), (1,2,4), (0,2,5))),
        # quad_L
        (("S4","S4R","S4RS","CPS4R","CPE4","CPE4R"),
         ((0,1), (1,2), (2,3), (3,0))),
        # quad_Q
        (("CPS8","CPS8R","CPE8","CPE8R","CPE8RP"),
         ((0,1,4), (1,2,5), (2,3,6), (3,0,7))),
        # tet_L
        (("C3D4","C3D4H","C3D4P","C3D4T","DC3D4"),
         ((0,1,2), (1,0,3), (2,1,3), (0,2,3))),
        # tet_Q
        (("C3D10","C3D10M","C3D10H","C3D10MH","C3D10MP"),
         ((0,1,2,4,5,6), (1,0,3,4,7,8), (2,1,3,5,8,9), (0,2,3,6,9,7))),
        # wedge_L
        (("C3D6","C3D6R","C3D6P","DC3D6","SC6R","COH3D6","COH3D6P",),
         ((0,1,2), (3,5,4), (0,3,4,1), (1,4,5,2), (0,2,5,3))),
        # wedge_Q
        (("C3D15",),
         ((0,1,2,6,7,8), (3,5,4,9,10,11), (0,3,4,1,12,9,13,8),
          (1,4,5,2,13,10,14,7), (0,2,5,3,8,14,11,12))),
        # hex_L
        (("C3D8","C3D8P","C3D8R","SC8R"),
         ((0,1,2,3), (4,7,6,5), (0,4,5,1), (1,5,6,2), (2,6,7,3),
          (0,3,7,4))),
        # hex_Q
        (("C3D20","C3D20R"),
         ((0,1,2,3,8,9,10,11), (4,7,6,5,12,13,14,15),
          (0,4,5,1,16,11,17,8), (1,5,6,2,17,13,18,9),
          (2,6,7,3,18,14,19,10), (0,3,7,4,11,19,15,16))),
        ]
    # convert to np.array (do it here separately from the dict initialization
    # then we have only one instance of each ndarray referenced multiple times
    # in the final dict
    abqNameToNodeIdsOnFace = [
        (nn, np.array(fn)) for nn, fn in abqNameToNodeIdsOnFace]
    # convert to dict {abqName: faceNodes}
    abqNameToNodeIdsOnFace = dict(
        (n, faceNodes)
        for names, faceNodes in abqNameToNodeIdsOnFace for n in names)


    #: Node indexes for each face of particular element types.
    #: The index ordering of each face is such that each faces normal
    #: points inwardly.
    abqNameToFaceNodes = [
        # tri
        (("S3","S3R","DS3","M3D3","CPE3",
          "CPS6M","CPE6","STRI65","DS6"), ((0,1), (1,2), (0,2))),
        # quad
        (("S4","S4R","S4RS","CPS4R","CPE4","CPE4R",
          "CPS8","CPS8R","CPE8","CPE8R","CPE8RP"),
         ((0,1), (1,2), (2,3), (3,0))),
        # tet
        (("C3D4","C3D4H","C3D4P","C3D4T","DC3D4",
          "C3D10","C3D10M","C3D10H","C3D10MH","C3D10MP"),
         ((0,1,2), (1,0,3), (2,1,3), (0,2,3))),
        # wedge
        (("C3D6","C3D6R","C3D6P","DC3D6","SC6R","C3D15","COH3D6","COH3D6P",),
         ((0,1,2), (3,5,4), (0,3,4,1), (1,4,5,2), (0,2,5,3))),
        # hex
        (("C3D8","C3D8P","C3D8R","SC8R","C3D20","C3D20R",),
         ((0,1,2,3), (4,7,6,5), (0,4,5,1), (1,5,6,2), (2,6,7,3),
          (0,3,7,4))),
        ]
    # convert to np.array (do it here separately from the dict initialization
    # then we have only one instance of each ndarray referenced multiple times
    # in the final dict
    abqNameToFaceNodes = [(nn, np.array(fn)) for nn, fn in abqNameToFaceNodes]
    # convert to dict {abqName: faceNodes}
    abqNameToFaceNodes = dict(
        (n, faceNodes)
        for names, faceNodes in abqNameToFaceNodes for n in names)


    def __init__(self, **kwargs):
        """Constructor. At most one of abqName, vtkCellType or shape must be
        given as keyword argument. An arbitrary number of additional keyword
        arguments can be passed and will overwrite corresponding attributes.

        Example:
         >>> elType = ElementType(abqName="C3D4")  # all defaults
         >>> print elType.order
         1
         >>> elType = ElementType(abqName="C3D4", order=2)  # cheating
         >>> print elType.order
         2

        @kwargs abqName: Abaqus element type label like C3D4, S3
        @kwargs vtkCellType: VTK element type number.
            Linear types: 3=line, 5=tri, 9=quad, 10=tet, 12=brick, 13=wedge.
            Quad types: 21=line, 22=tri, 23=quad, 24=tet, 25=hex
        @kwargs shape: "LINE", "TRI", "QUAD", "TET", "WEDGE", "HEX"
        """

        # check that at most one of abqName, vtkCellType or shape is given
        keyargs = set(("abqName", "vtkCellType", "shape")).intersection(kwargs)
        if len(keyargs) > 1:
            raise ValueError(
                "ElementType.__init__ got %d arguments including %s. From the"
                " latter list only one is allowed because they would possibly"
                " contradict each other." % (len(kwargs), keyargs))

        # initialize defaults based on abqName, vtkCellType or shape
        if "abqName" in kwargs:
            self.initFromAbqName(kwargs["abqName"])
        elif "vtkCellType" in kwargs:
            self.abqName = self.vtkCellTypeToAbqName[kwargs["vtkCellType"]][0]
            self.initFromAbqName(self.abqName)
        elif "shape" in kwargs:
            self.abqName = self.shapeToAbqName[kwargs["shape"]][0]
            self.initFromAbqName(self.abqName)

        # store given arguments
        for n, v in kwargs.iteritems():
            setattr(self, n,v)

    def initFromAbqName(self, abqName):
        self.shape = self.abqNameToShape[abqName]
        self.vtkCellType = self.abqNameToVtkCellType[abqName]

        self.nodesPerElem = self.abqNameToNodesPerElem[abqName]
        self.order = self.abqNameToOrder[abqName]
        self.nodesOnFace = self.abqNameToNodeIdsOnFace[abqName]
        self.faceNodes = self.abqNameToFaceNodes[abqName]

# some predefined default types
elTypeLinTet = ElementType(abqName="C3D4")
elTypeLinTri = ElementType(abqName="S3")
