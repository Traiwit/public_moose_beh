"""
I/O for VTK <https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf>.

Don't use this module directly. Consider this code and the module interface as
highly volatile!

taken from U{meshio<https://github.com/nschloe/meshio>}
"""
from functools import reduce
import collections

import numpy

from ..topo.structuredgrid import StructuredGrid
from ..topo.unstructuredmesh.meshbase import \
    HomoMesh, MeshTopoBase, MeshNodes, MeshElems, MeshGaussPts
from ..topo.unstructuredmesh.trimesh import TriMesh
from ..topo.unstructuredmesh.elementType import elTypeLinTet, elTypeLinTri

from bae.log_01 import msg

# from .._common import (
#     _meshio_to_vtk_order,
#     _vtk_to_meshio_order,
#     meshio_to_vtk_type,
#     vtk_to_meshio_type,
# )
# from .._exceptions import ReadError, WriteError
# from .._files import open_file
# from .._helpers import register
# from .._mesh import CellBlock, Mesh

# directly imported from meshio._exceptions:
class ReadError(Exception):
    pass

class WriteError(Exception):
    pass


# directly imported from meshio.__common:
def _vtk_to_meshio_order(vtk_type, numnodes, dtype=int):
    # meshio uses the same node ordering as VTK for most cell types. However,
    # for the linear wedge, the ordering of the gmsh Prism [1] is adopted since
    # this is found in most codes (Abaqus, Ansys, Nastran,...). In the vtkWedge
    # [2], the normal of the (0,1,2) triangle points outwards, while in gmsh
    # this normal points inwards.
    # [1] http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
    # [2] https://vtk.org/doc/nightly/html/classvtkWedge.html
    if vtk_type == 13:
        return numpy.array([0, 2, 1, 3, 5, 4], dtype=dtype)
    else:
        return numpy.arange(0, numnodes, dtype=dtype)


def _meshio_to_vtk_order(vtkCellType, numnodes, dtype=int):
    if vtkCellType == 13:  # wedge
        return numpy.array([0, 2, 1, 3, 5, 4], dtype=dtype)
    else:
        return numpy.arange(0, numnodes, dtype=dtype)

# directly imported from meshio._mesh:
CellBlock = collections.namedtuple("CellBlock", ["type", "data"])


vtk_type_to_numnodes = numpy.array(
    [
        0,  # empty
        1,  # vertex
        -1,  # poly_vertex
        2,  # line
        -1,  # poly_line
        3,  # triangle
        -1,  # triangle_strip
        -1,  # polygon
        -1,  # pixel
        4,  # quad
        4,  # tetra
        -1,  # voxel
        8,  # hexahedron
        6,  # wedge
        5,  # pyramid
        10,  # penta_prism
        12,  # hexa_prism
        -1,
        -1,
        -1,
        -1,
        3,  # line3
        6,  # triangle6
        8,  # quad8
        10,  # tetra10
        20,  # hexahedron20
        15,  # wedge15
        13,  # pyramid13
        9,  # quad9
        27,  # hexahedron27
        6,  # quad6
        12,  # wedge12
        18,  # wedge18
        24,  # hexahedron24
        7,  # triangle7
        4,  # line4
    ]
)


# These are all VTK data types.
# One sometimes finds 'vtktypeint64', but this is ill-formed.
vtk_to_numpy_dtype_name = {
    "bit": "bool",
    "unsigned_char": "uint8",
    "char": "int8",
    "unsigned_short": "uint16",
    "short": "int16",
    "unsigned_int": "uint32",
    "int": "int32",
    "unsigned_long": "uint64",
    "long": "int64",
    "float": "float32",
    "double": "float64",
}

numpy_to_vtk_dtype = {v: k for k, v in vtk_to_numpy_dtype_name.iteritems()}

# supported vtk dataset types
vtk_dataset_types = [
    "UNSTRUCTURED_GRID",
    "STRUCTURED_POINTS",
    "STRUCTURED_GRID",
    "RECTILINEAR_GRID",
]
# additional infos per dataset type
vtk_dataset_infos = {
    "UNSTRUCTURED_GRID": [],
    "STRUCTURED_POINTS": [
        "DIMENSIONS",
        "ORIGIN",
        "SPACING",
        "ASPECT_RATIO",  # alternative for SPACING in version 1.0 and 2.0
    ],
    "STRUCTURED_GRID": ["DIMENSIONS"],
    "RECTILINEAR_GRID": [
        "DIMENSIONS",
        "X_COORDINATES",
        "Y_COORDINATES",
        "Z_COORDINATES",
    ],
}

# all main sections in vtk
vtk_sections = [
    "METADATA",
    "DATASET",
    "POINTS",
    "CELLS",
    "CELL_TYPES",
    "POINT_DATA",
    "CELL_DATA",
    "LOOKUP_TABLE",
    "COLOR_SCALARS",
]


class Info:
    """Info Container for the VTK reader.

    @ivar points: Nx3 array of N point coordinates. For an unstructured grid.

    @ivar c: cell connectivity array, a flat array of integers containing the
        number of points/vertices for the first cell, then the point indexes for
        the first cell, then for the second cell and so on.
    @ivar ct: cell type array, contains the vtk cell type number for each cell
    @ivar cell_connectivity: Created py translate_cells() for unstructured
        grids only: A list of CellBlock items. Each one of those items
        -cellBlock- looks like that: cellBlock.type is the vtk cell type (an
        integer) and cellBlock.data is an NxM array of node indexes with N the
        number of cells (of this cell type) and M the number of nodes per cell.

    @ivar point_data: a dictionary {field name: field data ndarray} for all
        point data fields.
    @ivar cell_data_raw: a dictionary {field name: field data ndarray} for all
        cell data fields.
    @ivar cell_data: Created py translate_cells() for unstructured grids only:
        A dictionary {field name: list of data arrays}. The list of data arrays
        corresponds to the list of cell blocks for distinct cell types in the
        cell_connectivity attribute.
    @ivar field_data:

    @ivar section: 
        One of the problem in reading VTK files are POINT_DATA and CELL_DATA fields.
        They can contain a number of SCALARS+LOOKUP_TABLE tables, without giving and
        indication of how many there are. Hence, SCALARS must be treated like a
        first-class section.  To associate it with POINT/CELL_DATA, we store the
        'active' section in this variable.
    """

    def __init__(self):
        self.points = None
        self.field_data = collections.OrderedDict()
        self.cell_data_raw = collections.OrderedDict()
        self.point_data = collections.OrderedDict()
        self.dataset = {}
        self.c = None
        self.ct = None
        self.active = None
        self.is_ascii = False
        self.split = []
        self.num_items = 0
        # One of the problem in reading VTK files are POINT_DATA and CELL_DATA fields.
        # They can contain a number of SCALARS+LOOKUP_TABLE tables, without giving and
        # indication of how many there are. Hence, SCALARS must be treated like a
        # first-class section.  To associate it with POINT/CELL_DATA, we store the
        # `active` section in this variable.
        self.section = None


def read(inputPath):
    """Reads a VTK vtk file.
    """
    with open(inputPath, "rb") as f:
        out = read_buffer(f)
    return out


def read_buffer(f):
    """Reads an input stream, i.e. an open file-like object.

    @returns: all data in an L{Info} object/
    """
    
    # initialize output data
    info = Info()

    # skip header and read title
    f.readline()
    info.title = f.readline().strip()

    data_type = f.readline().decode("utf-8").strip().upper()
    if data_type not in ["ASCII", "BINARY"]:
        raise ReadError("Unknown VTK data type '{}'.".format(data_type))
    info.is_ascii = data_type == "ASCII"

    while True:
        line = f.readline().decode("utf-8")
        if not line:
            # EOF
            break

        line = line.strip()
        if len(line) == 0:
            # skip empty lines
            continue

        info.split = line.split()
        info.section = info.split[0].upper()

        if info.section in vtk_sections:
            _read_section(f, info)
        else:
            _read_subsection(f, info)

    # _check_mesh(info)  # not needed anymore, prepareTopo checks the data
    # translate_cells(info)  # prepareTopo now calls translate_cells
    topo = prepareTopo(info)
    return info, topo

def _read_section(f, info):
    if info.section == "METADATA":
        _skip_meta(f)

    elif info.section == "DATASET":
        info.active = "DATASET"
        info.dataset["type"] = info.split[1].upper()
        if info.dataset["type"] not in vtk_dataset_types:
            raise ReadError(
                "Only VTK '{}' supported (not {}).".format(
                    "', '".join(vtk_dataset_types), info.dataset["type"]
                )
            )

    elif info.section == "POINTS":
        info.active = "POINTS"
        info.num_points = int(info.split[1])
        data_type = info.split[2].lower()
        info.points = _read_points(f, data_type, info.is_ascii, info.num_points)

    elif info.section == "CELLS":
        info.active = "CELLS"
        info.num_items = int(info.split[2])
        info.c = _read_cells(f, info.is_ascii, info.num_items)

    elif info.section == "CELL_TYPES":
        info.active = "CELL_TYPES"
        info.num_items = int(info.split[1])
        info.ct = _read_cell_types(f, info.is_ascii, info.num_items)

    elif info.section == "POINT_DATA":
        info.active = "POINT_DATA"
        info.num_items = int(info.split[1])

    elif info.section == "CELL_DATA":
        info.active = "CELL_DATA"
        info.num_items = int(info.split[1])

    elif info.section == "LOOKUP_TABLE":
        info.num_items = int(info.split[2])
        data = numpy.fromfile(f, count=info.num_items * 4, sep=" ", dtype=float)
        rgba = data.reshape((info.num_items, 4))  # noqa F841

    elif info.section == "COLOR_SCALARS":
        nValues = int(info.split[2])
        # re-use num_items from active POINT/CELL_DATA
        num_items = info.num_items
        dtype = numpy.ubyte
        if info.is_ascii:
            dtype = float
        data = numpy.fromfile(f, count=num_items * nValues, dtype=dtype)


def _read_subsection(f, info):
    if info.active == "POINT_DATA":
        d = info.point_data
    elif info.active == "CELL_DATA":
        d = info.cell_data_raw
    elif info.active == "DATASET":
        d = info.dataset
    else:
        d = info.field_data

    if info.section in vtk_dataset_infos[info.dataset["type"]]:
        if info.section[1:] == "_COORDINATES":
            info.num_points = int(info.split[1])
            data_type = info.split[2].lower()
            d[info.section] = _read_coords(f, data_type, info.is_ascii, info.num_points)
        else:
            if info.section == "DIMENSIONS":
                d[info.section] = list(map(int, info.split[1:]))
            else:
                d[info.section] = list(map(float, info.split[1:]))
            if len(d[info.section]) != 3:
                raise ReadError(
                    "Wrong number of info in section '{}'. Need 3, got {}.".format(
                        info.section, len(d[info.section])
                    )
                )
    elif info.section == "SCALARS":
        d.update(_read_scalar_field(f, info.num_items, info.split, info.is_ascii))
    elif info.section == "VECTORS":
        d.update(_read_field(f, info.num_items, info.split, [3], info.is_ascii))
    elif info.section == "TENSORS":
        d.update(_read_field(f, info.num_items, info.split, [3, 3], info.is_ascii))
    elif info.section == "FIELD":
        d.update(_read_fields(f, int(info.split[2]), info.is_ascii))
    else:
        raise ReadError("Unknown section '{}'.".format(info.section))

def prepareTopo(info):
    """Create an appropriate field_02.topo object from the vtk dataset
    """
    if info.dataset["type"] == "UNSTRUCTURED_GRID":
        if info.c is None:
            raise ReadError("Required section CELLS not found.")
        if info.ct is None:
            raise ReadError("Required section CELL_TYPES not found.")
        translate_cells(info)
        distinctCellTypes = numpy.unique(info.ct)

        if len(distinctCellTypes)==1 and distinctCellTypes[0]==5:
            # only tri elements
            topo = TriMesh(
                nodeCoords=info.points,
                elNodes=info.cell_connectivity[0].data,
                elType=elTypeLinTri)

        elif len(distinctCellTypes)==1 and distinctCellTypes[0]==10:
            # only four node (linear) tet elements
            topo = TetMesh(
                nodeCoords=info.points,
                elNodes=info.cell_connectivity[0].data,
                elType=elTypeLinTet)

        elif len(distinctCellTypes)==1:
            # only elements of one kind
            raise NotImplementedError("General cell type not yet implemented.")
            # ... HomoMesh() what about nodesPerElem?
            # ??? nodesPerElem = vtk_type_to_numnodes[distinctCellTypes[0]]

        else:
            # different kinds of elements
            raise NotImplementedError(
                "Inhomogeneous cell types not yet implemented.")
        
    elif info.dataset["type"] == "STRUCTURED_POINTS":
        try:
            gridPtNb = info.dataset["DIMENSIONS"]
        except KeyError:
            raise ReadError("Required parameter DIMENSIONS not found for"
                            " STRUCTURED_POINTS.")
        try:
            origin = info.dataset["ORIGIN"]
        except KeyError:
            raise ReadError("Required parameter ORIGIN not found for"
                            " STRUCTURED_POINTS.")
        spacing = (
            info.dataset.get("SPACING") or info.dataset.get("ASPECT_RATIO"))
        if not spacing:
            raise ReadError("Neither SPACING nor ASPECT_RATIO found for"
                            " STRUCTURED_POINTS.")
        topo = StructuredGrid(
            origin=origin, spacing=spacing, gridPtNb=gridPtNb)
        
    elif info.dataset["type"] == "RECTILINEAR_GRID":
        raise NotImplementedError("reading dataset type RECTILINEAR_GRID not"
                                  " implemented yet.")
        if not all(x in info.dataset for x in (
                "X_COORDINATES", "Y_COORDINATES", "Z_COORDINATES")):
            raise ReadError("Required parameter X/Y/Z_COORDINATES not found"
                            " for RECTILINEAR_GRID.")

    return topo
    
def _check_mesh(info):
    """Just checking if all (?) required data is present for the topo/grid/mesh
    """
    # Originally this function also created points and cells data from
    # structured grid / -points data
    if info.dataset["type"] == "UNSTRUCTURED_GRID":
        if info.c is None:
            raise ReadError("Required section CELLS not found.")
        if info.ct is None:
            raise ReadError("Required section CELL_TYPES not found.")
    elif info.dataset["type"] == "STRUCTURED_POINTS":
        if not "DIMENSIONS" in info.dataset:
            raise ReadError("Required parameter DIMENSIONS not found for"
                            " STRUCTURED_POINTS.")
        if not "ORIGIN" in info.dataset:
            raise ReadError("Required parameter ORIGIN not found for"
                            " STRUCTURED_POINTS.")
        if not "SPACING" in info.dataset and not "ASPECT_RATIO" in info.dataset:
            raise ReadError("Neither SPACING nor ASPECT_RATIO found for"
                            " STRUCTURED_POINTS.")
    elif info.dataset["type"] == "RECTILINEAR_GRID":
        if not all(x in info.dataset for x in (
                "X_COORDINATES", "Y_COORDINATES", "Z_COORDINATES")):
            raise ReadError("Required parameter X/Y/Z_COORDINATES not found"
                            " for RECTILINEAR_GRID.")

    return
    # We don't need unstructured data (points, cells) from structured grid
    # data. The following chunk is deactivated and for reference purposes only.
    if info.dataset["type"] in ("STRUCTURED_POINTS", "RECTILINEAR_GRID"):
        info.points = _generate_points(info)
    if info.dataset["type"] in (
            "STRUCTURED_POINTS", "RECTILINEAR_GRID", "STRUCTURED_GRID"):
        info.c, info.ct = _generate_cells(dim=info.dataset["DIMENSIONS"])


def _generate_cells(dim):
    """Calculates and returns two arrays -cells and cell_types- from structured
    grid data. Cells is a flattened array containing for each cell the nb of
    vertices, followed by the indexes of the vertices. Cell_types is an array
    of single ints for each cell with all the same values: 3 for line cells
    (1-D), 9 for quad cells (2-D) and 12 for a hex cells (3D).

    Will only be called for info.dataset["type"] in ("STRUCTURED_POINTS",
    "RECTILINEAR_GRID", "STRUCTURED_GRID")
    """
    ele_dim = [d - 1 for d in dim if d > 1]
    ele_no = numpy.prod(ele_dim, dtype=int)
    spatial_dim = len(ele_dim)

    if spatial_dim == 1:
        # cells are lines in 1D
        cells = numpy.empty((ele_no, 3), dtype=int)
        cells[:, 0] = 2
        cells[:, 1] = numpy.arange(ele_no, dtype=int)
        cells[:, 2] = cells[:, 1] + 1
        cell_types = numpy.full(ele_no, 3, dtype=int)
    elif spatial_dim == 2:
        # cells are quad in 2D
        cells = numpy.empty((ele_no, 5), dtype=int)
        cells[:, 0] = 4
        cells[:, 1] = numpy.arange(0, ele_no, dtype=int)
        cells[:, 1] += numpy.arange(0, ele_no, dtype=int) // ele_dim[0]
        cells[:, 2] = cells[:, 1] + 1
        cells[:, 3] = cells[:, 1] + 2 + ele_dim[0]
        cells[:, 4] = cells[:, 3] - 1
        cell_types = numpy.full(ele_no, 9, dtype=int)
    else:
        # cells are hex in 3D
        cells = numpy.empty((ele_no, 9), dtype=int)
        cells[:, 0] = 8
        cells[:, 1] = numpy.arange(ele_no)
        cells[:, 1] += (ele_dim[0] + ele_dim[1] + 1) * (
            numpy.arange(ele_no) // (ele_dim[0] * ele_dim[1])
        )
        cells[:, 1] += (numpy.arange(ele_no) % (ele_dim[0] * ele_dim[1])) // ele_dim[0]
        cells[:, 2] = cells[:, 1] + 1
        cells[:, 3] = cells[:, 1] + 2 + ele_dim[0]
        cells[:, 4] = cells[:, 3] - 1
        cells[:, 5] = cells[:, 1] + (1 + ele_dim[0]) * (1 + ele_dim[1])
        cells[:, 6] = cells[:, 5] + 1
        cells[:, 7] = cells[:, 5] + 2 + ele_dim[0]
        cells[:, 8] = cells[:, 7] - 1
        cell_types = numpy.full(ele_no, 12, dtype=int)

    return cells.reshape(-1), cell_types


def _generate_points(info):
    """generates an array of point coordinates from structured grid data.

    Should only be called for info.dataset["type"] in ("STRUCTURED_POINTS",
    "RECTILINEAR_GRID")
    """

    if info.dataset["type"] == "STRUCTURED_POINTS":
        dim = info.dataset["DIMENSIONS"]
        ori = info.dataset["ORIGIN"]
        try:
            spa = info.dataset["SPACING"]
        except KeyError:
            spa = info.dataset["ASPECT_RATIO"]
        axis = [
            numpy.linspace(ori[i], ori[i] + (dim[i] - 1.0) * spa[i], dim[i])
            for i in range(3)
        ]
    elif info.dataset["type"] == "RECTILINEAR_GRID":
        axis = [
            info.dataset["X_COORDINATES"],
            info.dataset["Y_COORDINATES"],
            info.dataset["Z_COORDINATES"],
        ]
    elif info.dataset["type"] == "STRUCTURED_GRID":
        # structured grid apparrently has no points ?
        return None

    x_dim = len(axis[0])
    y_dim = len(axis[1])
    z_dim = len(axis[2])
    pnt_no = x_dim * y_dim * z_dim
    x_id, y_id, z_id = numpy.mgrid[0:x_dim, 0:y_dim, 0:z_dim]
    points = numpy.empty((pnt_no, 3), dtype=axis[0].dtype)
    # VTK sorts points and cells in Fortran order
    points[:, 0] = axis[0][x_id.reshape(-1, order="F")]
    points[:, 1] = axis[1][y_id.reshape(-1, order="F")]
    points[:, 2] = axis[2][z_id.reshape(-1, order="F")]
    return points


def _read_coords(f, data_type, is_ascii, num_points):
    dtype = numpy.dtype(vtk_to_numpy_dtype_name[data_type])
    if is_ascii:
        coords = numpy.fromfile(f, count=num_points, sep=" ", dtype=dtype)
    else:
        # Binary data is big endian, see
        # <https://www.vtk.org/Wiki/VTK/Writing_VTK_files_using_python#.22legacy.22>.
        dtype = dtype.newbyteorder(">")
        coords = numpy.fromfile(f, count=num_points, dtype=dtype)
        # moved skipping empty lines
        # line = f.readline().decode("utf-8")
        # if line != "\n":
        #     raise ReadError()
    return coords


def _read_points(f, data_type, is_ascii, num_points):
    """Read the point coordinates e.g. for an unstructured mesh from file

    @param f: input stream
    @returns: Nx3 array of N point coordinates
    """
    ###GP: need to check: Are points always 3D in vtk?
    dtype = numpy.dtype(vtk_to_numpy_dtype_name[data_type])
    if is_ascii:
        points = numpy.fromfile(f, count=num_points * 3, sep=" ", dtype=dtype)
    else:
        # Binary data is big endian, see
        # <https://www.vtk.org/Wiki/VTK/Writing_VTK_files_using_python#.22legacy.22>.
        dtype = dtype.newbyteorder(">")
        points = numpy.fromfile(f, count=num_points * 3, dtype=dtype)
        # moved skipping empty lines
        # line = f.readline().decode("utf-8")
        # if line != "\n":
        #     raise ReadError()
    return points.reshape((num_points, 3))


def _read_cells(f, is_ascii, num_items):
    if is_ascii:
        c = numpy.fromfile(f, count=num_items, sep=" ", dtype=int)
    else:
        c = numpy.fromfile(f, count=num_items, dtype=">i4")
        # moved skipping empty lines
        # line = f.readline().decode("utf-8")
        # if line != "\n":
        #     raise ReadError()
    return c


def _read_cell_types(f, is_ascii, num_items):
    if is_ascii:
        ct = numpy.fromfile(f, count=int(num_items), sep=" ", dtype=int)
    else:
        # binary
        ct = numpy.fromfile(f, count=int(num_items), dtype=">i4")
        # moved skipping empty lines
        # line = f.readline().decode("utf-8")
        # # Sometimes, there's no newline at the end
        # if line.strip() != "":
        #     raise ReadError()
    return ct


def _read_scalar_field(f, num_data, split, is_ascii):
    data_name = split[1]
    data_type = split[2].lower()
    try:
        num_comp = int(split[3])
    except IndexError:
        num_comp = 1

    # The standard says:
    # > The parameter numComp must range between (1,4) inclusive; [...]
    if not (0 < num_comp < 5):
        raise ReadError("The parameter numComp must range between (1,4) inclusive")

    dtype = numpy.dtype(vtk_to_numpy_dtype_name[data_type])
    lt, _ = f.readline().decode("utf-8").split()
    if lt.upper() != "LOOKUP_TABLE":
        raise ReadError()

    if is_ascii:
        data = numpy.fromfile(f, count=num_data, sep=" ", dtype=dtype)
    else:
        # Binary data is big endian, see
        # <https://www.vtk.org/Wiki/VTK/Writing_VTK_files_using_python#.22legacy.22>.
        dtype = dtype.newbyteorder(">")
        data = numpy.fromfile(f, count=num_data, dtype=dtype)
        # moved skipping empty lines
        # line = f.readline().decode("utf-8")
        # if line and line != "\n":
        #     raise ReadError(
        #         "Expected empty line, but got line of length %d:\n%s"
        #         % (len(line), line))

    # Note: We generally use ordered dicts for fields. As we return only one
    # field in this case it's ok to use an ordinary dict.
    return {data_name: data}


def _read_field(f, num_data, split, shape, is_ascii):
    data_name = split[1]
    data_type = split[2].lower()

    dtype = numpy.dtype(vtk_to_numpy_dtype_name[data_type])
    # <https://stackoverflow.com/q/2104782/353337>
    k = reduce((lambda x, y: x * y), shape)
    
    if is_ascii:
        data = numpy.fromfile(f, count=k * num_data, sep=" ", dtype=dtype)
    else:
        # Binary data is big endian, see
        # <https://www.vtk.org/Wiki/VTK/Writing_VTK_files_using_python#.22legacy.22>.
        dtype = dtype.newbyteorder(">")
        data = numpy.fromfile(f, count=k * num_data, dtype=dtype)
        # moved skipping empty lines
        # line = f.readline().decode("utf-8")
        # if line and line != "\n":
        #     raise ReadError(
        #         "Expected empty line, but got line of length %d:\n%s"
        #         % (len(line), line))

    data = data.reshape(-1, *shape)

    # Note: We generally use ordered dicts for fields. As we return only one
    # field in this case it's ok to use an ordinary dict.
    return {data_name: data}


def _read_fields(f, num_fields, is_ascii):
    """reads FIELD dataset attribute for cell data or point data datasets.

    @param f: input stream
    @param num_fields: number of arrays in this "field" dataset attribute
    @param is_ascii:
    """
    data = collections.OrderedDict()
    for _ in range(num_fields):
        # skip empty lines
        while 1:
            line = f.readline().decode("utf-8")
            if not line:
                raise ReadError("Reached end-of-file unexpectedly.")
            line = line.split()
            if line:
                break

        if line[0] == "METADATA":
            _skip_meta(f)
            name, shape0, shape1, data_type = f.readline().decode("utf-8").split()
        else:
            name, shape0, shape1, data_type = line

        shape0 = int(shape0)
        shape1 = int(shape1)
        dtype = numpy.dtype(vtk_to_numpy_dtype_name[data_type.lower()])

        if is_ascii:
            dat = numpy.fromfile(f, count=shape0 * shape1, sep=" ", dtype=dtype)
        else:
            # Binary data is big endian, see
            # <https://www.vtk.org/Wiki/VTK/Writing_VTK_files_using_python#.22legacy.22>.
            dtype = dtype.newbyteorder(">")
            dat = numpy.fromfile(f, count=shape0 * shape1, dtype=dtype)
            # moved skipping empty lines
            # line = f.readline().decode("utf-8")
            # if line and line != "\n":
            #     raise ReadError(
            #         "Expected empty line, but got line of length %d:\n%s"
            #         % (len(line), line))

        if shape0 == 9:
            dat = dat.reshape((shape1, 3, 3))
        elif shape0 != 1:
            dat = dat.reshape((shape1, shape0))

        data[name] = dat

    return data


def _skip_meta(f):
    # skip possible metadata
    # https://vtk.org/doc/nightly/html/IOLegacyInformationFormat.html
    while True:
        line = f.readline().decode("utf-8").strip()
        if not line:
            # end of metadata is a blank line
            break


def translate_cells(info):
    """For unstructured grid data sets restructures the cell connectivity data
    and the cell data.

    Adds cell_connectivity and cell_data items to the given info object.

    cell_connectivity is a list of CellBlock items. Each one of those items
    -cellBlock- looks like that: cellBlock.type is the vtk cell type (an
    integer) and cellBlock.data is an NxM array of node indexes with N the
    number of cells (of this cell type) and M the number of nodes per cell.

    cell_data is a dictionary {field name: list of data arrays}. The list of
    data arrays corresponds to the list of cell types.
    """
    # https://www.vtk.org/doc/nightly/html/vtkCellType_8h_source.html
    # Translate it into the cells array.
    # conn_raw / info.ct is a one-dimensional vector with
    # (num_points0, p0, p1, ... ,pk, numpoints1, p10, p11, ..., p1k, ...
    
    conn_raw = info.c
    types = info.ct
    cell_data_raw = info.cell_data_raw
    
    cells = []
    cell_data = collections.OrderedDict()

    # Deduct offsets from the cell types. This is much faster than manually going
    # through the conn_raw array. Slight disadvantage: This doesn't work for
    # cells with a custom number of points like polygon --vtk cell type 7.
    # But we're not using polygons.
    numnodes = vtk_type_to_numnodes[types]
    if not numpy.all(numnodes > 0):
        raise ReadError("File contains cells of types that we cannot handle.")

    # offsets gives the first index in conn_raw (connectivity array) for each cell
    offsets = numpy.cumsum(numnodes + 1) - (numnodes + 1)

    if not numpy.all(numnodes == conn_raw[offsets]):
        raise ReadError(
            "The nb of vertices in cell data is not consistent with the nb"
            " of vertices deducted from the cell types.")

    # array of indexes at which the cell type changes
    # i.e.  types = [1,1,1,1,4,4,4,4,2,2,2] ==> b = [0, 4, 8, 11]
    b = numpy.concatenate(
        [[0], numpy.where(types[:-1] != types[1:])[0] + 1, [len(types)]]
    )

    # iterate over blocks of cells with one cell type
    for start, end in zip(b[:-1], b[1:]):
        cell_type = types[start]
        # meshio_type = vtk_to_meshio_type[types[start]]
        n = conn_raw[offsets[start]]
        cell_idx = _vtk_to_meshio_order(types[start], n, dtype=offsets.dtype)
        indices = numpy.add.outer(offsets[start:end], cell_idx + 1)
        cells.append(CellBlock(cell_type, conn_raw[indices]))
        # cells.append(CellBlock(meshio_type, conn_raw[indices]))
        for name, d in cell_data_raw.items():
            if name not in cell_data:
                cell_data[name] = []
            cell_data[name].append(d[start:end])

    info.cell_connectivity = cells
    info.cell_data = cell_data

#{ writing vtk files
def write(filename, topo, fields, iTime=0, description="Field data", binary=True):
    """
    @param filename:

    @param description: String for the second line of the vtk data file
        (title field).

    @param topo: The L{FieldsCollection.topo} object, subclass of
        L{MeshTopoBase} or L{StructuredGrid}
    @param fields: A dictionary {fieldName: field data} like
        L{FieldsCollection.fields}.
    @param iTime: the time index to be used for fields
    @param binary: True for binary or False for ascii
    """
    def pad(array):
        # pad array with "constant" values. ((0,0),(0,1)) means: for the first
        # axis add "0" values before and "0" after, for the second axis add
        # "0' values before and "1" after.
        return numpy.pad(array, ((0, 0), (0, 1)), "constant")

    with open(filename, "wb") as f:

        # version changed for compatibility, was: 4.2
        f.write(b"# vtk DataFile Version 2.0\n")
        # f.write("written by meshio v{}\n".format(__version__).encode("utf-8"))
        f.write("{}\n".format(description).encode("utf-8"))
        f.write(("BINARY\n" if binary else "ASCII\n").encode("utf-8"))
        
        if isinstance(topo, MeshTopoBase):

            # topo is an unstructured mesh
            if topo.mesh.nodeCoords.shape[1] == 2:
                msg("WARNING: VTK requires 3D points, but 2D points given."
                    " Appending 0 third component.")
                points = pad(topo.mesh.nodeCoords)
            else:
                points = topo.mesh.nodeCoords

            f.write(b"DATASET UNSTRUCTURED_GRID\n")

            # write points and cells
            _write_points(f, points, binary)
            _write_cells(f, topo.mesh, binary)

            # write point data
            if isinstance(topo, MeshNodes):
                # nodal fields
                f.write("POINT_DATA {}\n".format(len(topo)).encode("utf-8"))
                _write_field_data(f, fields, iTime, binary)

            # write cell data
            elif isinstance(topo, MeshElems):
                f.write("CELL_DATA {}\n".format(len(topo)).encode("utf-8"))
                _write_field_data(f, fields, iTime, binary)

            else:
                raise NotImplementedError(
                    "VTK export for topo of type {} not implemented"
                    .format(type(topo)))
        
        elif isinstance(topo, StructuredGrid):
        
            f.write(b"DATASET STRUCTURED_POINTS\n")
            f.write("DIMENSIONS {} {} {}\n"
                    .format(*topo.gridPtNb).encode("utf-8"))
            f.write("ORIGIN {} {} {}\n"
                    .format(*topo.origin).encode("utf-8"))
            f.write("SPACING {} {} {}\n"
                    .format(*topo.spacing).encode("utf-8"))
            f.write("POINT_DATA {}\n".format(len(topo)).encode("utf-8"))
            _write_field_data(f, fields, iTime, binary)
        
        else:
            raise NotImplementedError(
                "Writing topo type {} to vtk not imlemented, yet."
                .format(type(topo)))


def _write_points(f, points, binary):
    f.write(
        "POINTS {} {}\n".format(
            len(points), numpy_to_vtk_dtype[points.dtype.name]
        ).encode("utf-8")
    )

    if binary:
        # Binary data must be big endian, see
        # <https://www.vtk.org/Wiki/VTK/Writing_VTK_files_using_python#.22legacy.22>.
        # if points.dtype.byteorder == "<" or (
        #     points.dtype.byteorder == "=" and sys.byteorder == "little"
        # ):
        #     logging.warn("Converting to new byte order")
        points.astype(points.dtype.newbyteorder(">")).tofile(f, sep="")
    else:
        # ascii
        points.tofile(f, sep=" ")
    f.write(b"\n")


def _write_cells(f, mesh, binary):
    """write cell connectivity array for unstructured meshes to vtk file
    """
    if isinstance(mesh, HomoMesh):
        total_num_cells = len(mesh.elNodes)
        total_num_idx = mesh.elType.nodesPerElem * total_num_cells
        # For each cell, the number of nodes is stored
        total_num_idx += total_num_cells
        f.write("CELLS {} {}\n".format(total_num_cells, total_num_idx).encode("utf-8"))
        n = len(mesh.elNodes)  # ==mesh.elNodes.shape[0]
        m = mesh.elType.nodesPerElem  # ==mesh.elNodes.shape[1]
        cell_idx = _meshio_to_vtk_order(mesh.elType.vtkCellType, m)
        if binary:
            dtype = numpy.dtype(">i4")
            # One must force endianness here:
            # <https://github.com/numpy/numpy/issues/15088>
            # prepend a column with the value m
            numpy.column_stack(
                [
                    numpy.full(n, m, dtype=dtype),
                    mesh.elNodes[:, cell_idx].astype(dtype),
                ],
            ).astype(dtype).tofile(f, sep="")
            f.write(b"\n")
        else:
            # ascii
            # prepend a column with the value m
            numpy.column_stack(
                [
                    numpy.full(n, m, dtype=mesh.elNodes.dtype),
                    mesh.elNodes[:, cell_idx],
                ]
            ).tofile(f, sep="\n")
            f.write(b"\n")

        # write cell types
        f.write("CELL_TYPES {}\n".format(total_num_cells).encode("utf-8"))
        if binary:
            numpy.full(n, mesh.elType.vtkCellType, dtype=dtype).tofile(
                f, sep="")
            f.write(b"\n")
        else:
            # ascii
            numpy.full(n, mesh.elType.vtkCellType).tofile(f, sep="\n")
            f.write(b"\n")

    elif isinstance(mesh, HeteroMesh):

        raise NotImplementedError(
            "VTK export of heteromesh topo not implemented yet.")

        total_num_cells = sum([len(c.data) for c in cells])
        total_num_idx = sum([c.data.size for c in cells])
        # For each cell, the number of nodes is stored
        total_num_idx += total_num_cells
        f.write("CELLS {} {}\n".format(total_num_cells, total_num_idx).encode("utf-8"))
        if binary:
            for c in cells:
                n = c.data.shape[1]
                cell_idx = _meshio_to_vtk_order(c.type, n)
                dtype = numpy.dtype(">i4")
                # One must force endianness here:
                # <https://github.com/numpy/numpy/issues/15088>
                numpy.column_stack(
                    [
                        numpy.full(c.data.shape[0], n, dtype=dtype),
                        c.data[:, cell_idx].astype(dtype),
                    ],
                ).astype(dtype).tofile(f, sep="")
            f.write(b"\n")
        else:
            # ascii
            for c in cells:
                n = c.data.shape[1]
                cell_idx = _meshio_to_vtk_order(c.type, n)
                # prepend a column with the value n
                numpy.column_stack(
                    [
                        numpy.full(c.data.shape[0], n, dtype=c.data.dtype),
                        c.data[:, cell_idx],
                    ]
                ).tofile(f, sep="\n")
                f.write(b"\n")

        # write cell types
        f.write("CELL_TYPES {}\n".format(total_num_cells).encode("utf-8"))
        if binary:
            for c in cells:
                key_ = c.type[:7] if c.type[:7] == "polygon" else c.type
                vtk_type = meshio_to_vtk_type[key_]
                numpy.full(len(c.data), vtk_type, dtype=numpy.dtype(">i4")).tofile(
                    f, sep=""
                )
            f.write(b"\n")
        else:
            # ascii
            for c in cells:
                key_ = c.type[:7] if c.type[:7] == "polygon" else c.type
                numpy.full(len(c.data), meshio_to_vtk_type[key_]).tofile(f, sep="\n")
                f.write(b"\n")

    else:
        raise NotImplementedError("Unrecognized mesh type {}".format(type(mesh)))


def _write_field_data(f, fields, iTime, binary):
    f.write(("FIELD FieldData {}\n".format(len(fields))).encode("utf-8"))
    for name, values in fields.items():
        if len(values.shape) == 2:  # scalar value has only space and time dim
            num_tuples = values.shape[0]
            num_components = 1
        elif len(values.shape) == 3:  # vector value
            num_tuples = values.shape[0]
            num_components = values.shape[2]
        elif len(values.shape) == 4:  # tensor value
            num_tuples = values.shape[0]
            num_components = values.shape[2]*values.shape[3]
        else:
            raise ValueError(
                "_write_field_data called with a field of shape %s." % values.shape)

        if " " in name:
            raise WriteError(
                "VTK doesn't support spaces in field names ('{}').".format(name)
            )

        f.write(
            (
                "{} {} {} {}\n".format(
                    name,
                    num_components,
                    num_tuples,
                    numpy_to_vtk_dtype[values.dtype.name],
                )
            ).encode("utf-8")
        )
        if binary:
            values[:,iTime,...].astype(values.dtype.newbyteorder(">")).tofile(f, sep="")
        else:
            # ascii
            values[:,iTime,...].tofile(f, sep=" ")
            # numpy.savetxt(f, points)
        f.write(b"\n")
