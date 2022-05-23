"""write vtk data files
"""

from warnings import warn
class VtkWarn(UserWarning):
    pass


class BaseElementType(object):

    #--- must be specified in subclass
    nb_ips = "number of integration points"
    nb_nodes = "number of nodes"
    vtk_celltype = "corresponding vtk cell type"

    #--- may be specified in subclass (default value given here)
    # if this member attribute is specified in a subclass it is a list of the
    # vtk node indices in the order of the corresponding abaqus element node
    # index:
    # vtk_abqnodeidx[index in vtk file] := abaqus node index (zero based!)
    vtk_abqnodeidx = None
    pass

class ElementTypeC3D10M(BaseElementType):
    description = "10 node quadratic tet"
    nb_ips = 4
    nb_nodes = 10
    vtk_celltype = 24

class ElementTypeCOH3D6(BaseElementType):
    description = "6 node wedge cohesive element"
    nb_ips = 3
    nb_nodes = 6
    vtk_celltype = 13
    vtk_abqnodeidx = [0, 2, 1, 3, 5, 4]

class ElementTypeB31(BaseElementType):
    description = "2 node linear beam"
    nb_ips = 1
    nb_nodes = 2
    vtk_celltype = 3

ii = len("ElementType")
supportedTypes = dict([
        (cname[ii:], locals()[cname])
        for cname in dir()
        if cname[:ii]=="ElementType"
        ])   
del ii
# print "supported types:", supportedTypes.keys()


class FieldData(object):
    """class to temporarily hold the data for one data field"""
    def __init__(self, fieldName, fieldType, dataList):
        self.name = fieldName
        self.type = fieldType
        self.dataList = dataList
    def getFieldName(self):
        return self.name


class VtkFile(file):
    """Class for *writing* vtk files. (Simple Legacy Format)"""
    def __init__(self, output_filename):
        file.__init__(self, output_filename , 'w')

        # write vtk header
        self.write('# vtk DataFile Version 3.0\n')
        self.write('vtk output\n')
        self.write('ASCII\n')
        self.write('DATASET UNSTRUCTURED_GRID\n')


    def addMeshData(self, model, nbDim = 3):
        """write mesh data

        model may be a abq_model_02.Model object or something that has similar
        nodeCoords, elNodes and elType members

        nbDim is the number of point coordinates (default:3)
        """
        
        self.nb_elements = len(model.elNodes)
        self.nb_nodes = len(model.nodeCoords)

        # extract undeformed node coordinates
        # modelNodeVtkNode: {abaqus node number : vtk node index}
        # modelNodes[ modelNodeVtkNode[num] ] = num
        self.modelNodeVtkNode = dict()
        self.modelNodes = range(self.nb_nodes)

        self.write('\nPOINTS %s float\n'%(self.nb_nodes))
        vtkNodeCnt = 0
        nodelist = model.nodeCoords.keys()
        nodelist.sort()
        for node in nodelist:
            coords = model.nodeCoords[node]
            self.modelNodeVtkNode[node] = vtkNodeCnt
            self.modelNodes[vtkNodeCnt] = node
            # list of undeformed node positions
            node_pos = " ".join(["%s"%x for x in coords[:nbDim]])
            self.write(node_pos+'\n')
            vtkNodeCnt += 1

        # extract element/cell connectivity
        # modelElemVtkCell: {abaqus element number : vtk cell index}
        # modelElems[ modelElemVtkCell[num] ] = num
        # cellList: list of lists [number nodes, vtk node indices]
        self.modelElemVtkCell = dict()
        self.modelElems = range(self.nb_elements)
        cellList = list()
        modelElemType = dict()

        vtkElemCnt = 0
        cellListSize = 0
        elemlist = model.elNodes.keys()
        elemlist.sort()
        for elt_i in elemlist:
            elt_type = model.elType[elt_i]
            try:
                elType = supportedTypes[elt_type]
            except KeyError:
                msg = 'Type %s not supported yet!' %(elt_type)
                raise KeyError(msg)

            modelElemType[elt_i] = elType
            self.modelElemVtkCell[elt_i] = vtkElemCnt
            self.modelElems[vtkElemCnt] = elt_i
            if elType.vtk_abqnodeidx:
                connectivity = [model.elNodes[elt_i][i] for i in elType.vtk_abqnodeidx]
            else:
                connectivity = model.elNodes[elt_i]
            cellList.append([elType.nb_nodes]
                            + [self.modelNodeVtkNode[k] for k in connectivity])
            cellListSize += 1+elType.nb_nodes
            vtkElemCnt += 1

        # write CELLS list
        self.write('\nCELLS %s %s\n'%(self.nb_elements, cellListSize))
        for line in cellList:
            self.write(' '.join(map(str, line)) + '\n')

        # write CELL_TYPES
        self.write('\nCELL_TYPES %s\n'%(self.nb_elements))
        for elt_i in self.modelElems:
            self.write('%s\n' % modelElemType[elt_i].vtk_celltype)

        # attribute to store all the data
        self.pointData = list()
        self.cellData = list()

    def addFieldData(self, fieldName, fieldType, pointOrCell, dataDict):
        """write one field to the vtk file
        
        fieldType may be "scalar", "vector" or "tensor"
        pointOrCell may be "point" or "cell" indicating whether the data is
           point (nodal) or cell (element) data.
        dataDict is a dict {node/element nb : value}
        """

        if pointOrCell=="point":
            dataList = [dataDict.get(node, None) for node in self.modelNodes]
            self.pointData.append(FieldData(fieldName, fieldType, dataList))
            
        elif pointOrCell=="cell":
            dataList = [dataDict.get(elem, None) for elem in self.modelElems]
            self.cellData.append(FieldData(fieldName, fieldType, dataList))
        else:
            raise ValueError(
                "VtkFile.addDataList(): invalid pointOrCell argument '%s'."
                " Must be 'point' or 'cell'."
                % pointOrCell)

    def writeFieldData(self):
        """this function has to be called once all fields have been added with
        the addFieldData function.
        """
        for vtk_keyword, nb_pos, fieldData in (
            ("POINT_DATA", self.nb_nodes, self.pointData),
            ("CELL_DATA", self.nb_elements, self.cellData)):

            if len(fieldData):
                fieldData.sort(None, key=FieldData.getFieldName)
                self.write('\n%s %s\n'%(vtk_keyword,nb_pos))
                for field in fieldData:
                    if field.type == "scalar":
                        self.write('SCALARS %s float 1\n'%self.getValidVtkName(field.name))
                        self.write('LOOKUP_TABLE default\n')
                        for data in field.dataList:
                            if data == None:
                                data = 0
                            self.write('%s\n'%data)
                    elif field.type=="vector":
                        self.write('VECTORS %s float\n'%(self.getValidVtkName(field.name)))
                        for data in field.dataList:
                            if data == None:
                                data = [0,0,0]
                            elif len(data)<3:
                                raise Exception("In field %s of type %s len<3:\n%s"
                                                % (field.name, field.type, data))
                            self.write(' '.join(map(str,data))+'\n')
                    elif field.type=="tensor":
                        self.write('TENSORS %s float\n'%(self.getValidVtkName(field.name)))
                        for data in field.dataList:
                            if data == None:
                                data = [0,0,0,0,0,0]
                            elif len(data)<6:
                                raise Exception("In field %s of type %s len<6:\n%s"
                                                % (field.name, field.type, data))
                            self.write(' '.join([str(data[i])
                                                 for i in [0,3,4,3,1,5,4,5,2]])+'\n')
                    else:
                        warn(VtkWarn("Don't know how to write data of type %s yet!"
                                     %(field.type)))



    @staticmethod
    def getValidVtkName(field_keyword):
        ''' convert keyword in string acceptable for VTK (no spaces) '''
        return field_keyword.replace(' ' ,'_')


#------------------------------------------------------------------------------

class DatasetAttribute(list):
    """Class for holding one set of ("dataset attribute-") data.
    Base class for more specific types.

    @ivar pointsCells: POINT DATA -> "point"; CELL DATA -> "cell"
    """
    __slots__ = ["pointsCells"]

class DatasetAttributeScalar(DatasetAttribute):
    """Class for holding one set of ("dataset attribute-") data of type SCALAR.

    @ivar pointsCells: POINT DATA -> "point"; CELL DATA -> "cell"
    @ivar dataName: the name of the attribute
    @ivar dataType: may be bit, unsigned_char, char, unsigned_short, short,
    unsigned_int, int, unsigned_long, long, float, or double.
    @ivar numComp: number of components (1...4)
    @ivar lookupTable: name of the lookup table or "default"
    """
    __slots__ = ["pointsCells", "dataName", "dataType", "numComp",
                 "lookupTable"]

class VtkReader(object):
    """Read dat from an vtk file (Simple Legacy Format).
    Does not work generally, only for our special case: Only points and scalar
    attributes.

    Docs on the vtk file format: search the web for vtk file formats pdf
    or see D:\gero\doc\vtk-file-formats.pdf ...

    @ivar points: list of coordinate tuples (lists)
    @ivar point_data: dict {field name: list of values}
    @ivar cell_data: dict {field name: list of values}
    """

    _dataTypeToConvFunc = {
        "float": float,
        "double": float,
        }

    def __init__(self):
        # list of coordinate tuples (lists)
        self.points = list()
        
        # dict {field name: list of values}
        self.point_data = dict()
        self.cell_data = dict()
        
    def readSpecial(self, fileName):
        vtkFile = open(fileName)

        vtkFile.next() # file identifier and format version 
        vtkFile.next() # header string

        # format ASCII or BINARY, can only read ASCII
        format = vtkFile.next().strip().upper()
        assert format=="ASCII"

        for part_header in vtkFile:
            if part_header=="":
                continue # ignore empty lines
            part_header = part_header.strip().upper().split()

            if part_header[0] == "DATASET":
                dataset_type = part_header[1]
                if dataset_type == "UNSTRUCTURED_GRID":
                    dataset_item_type = ""
                    while dataset_item_type=="": # ignore empty lines
                        dataset_item_type = vtkFile.next()
                    dataset_item_type = dataset_item_type.strip(
                        ).upper().split()

                    if dataset_item_type[0] == "POINTS":
                        nb_pts = int(dataset_item_type[1])
                        data_type = dataset_item_type[2]
                        assert data_type=="FLOAT"
                        self.points = [
                            map(float, vtkFile.next().strip().split())
                            for i in xrange(nb_pts)]
                    else:
                        ValueError(
                            "dataset item type %s (after DATASET"
                            " UNSTRUCTURED_GRID) not implemented."
                            % dataset_item_type[0])
                else:
                    ValueError("DATASET type %s not implemented."
                               % dataset_type)


            elif part_header[0]=="POINT_DATA" or part_header[0]=="CELL_DATA":
                if part_header[0]=="POINT_DATA":
                    point_cells = "point"
                    attr = self.point_data
                elif part_header[0]=="CELL_DATA":
                    point_cells = "cell"
                    attr = self.cell_data
                    raise ValueError("cell data not implemented.")

                nb_values = int(part_header[1])
                for data_item_type in vtkFile:
                    if data_item_type=="":
                        continue # ignore empty lines
                    data_item_type = data_item_type.strip().split()
                    data_item_type[0] = data_item_type[0].upper()
                    data_name = data_item_type[1]
                    if data_item_type[0] == "SCALARS":

                        # evaluate data type (third value on the line)
                        data_type = data_item_type[2].lower()
                        try:
                            converter = self._dataTypeToConvFunc[data_type]
                        except KeyError:
                            raise ValueError("dataType %s not implemented yet."
                                             % data_type)

                        # evaluate number of components (optional 4. value on
                        # the line)
                        try:
                            nb_components = data_item_type[3]
                        except IndexError:
                            nb_components = 1

                        # read the lookup table
                        lookup_table = vtkFile.next().split()
                        assert lookup_table[0].upper() == "LOOKUP_TABLE"
                        lookup_table = lookup_table[1]

                        # read the data
                        if nb_components == 1:
                            new_data = DatasetAttributeScalar(
                                converter(vtkFile.next())
                                for i in xrange(nb_values))
                        else:
                            new_data = DatasetAttributeScalar(
                                map(converter, vtkFile.next().strip().split())
                                for i in xrange(nb_values))

                        # assign other attributes
                        new_data.pointsCells = point_cells
                        new_data.dataName = data_name
                        new_data.dataType = data_type
                        new_data.numComp = nb_components
                        new_data.lookupTable = lookup_table

                        # store new data in the dictionary
                        if data_name in attr:
                            raise ValueError(
                                "Data %s appears multiple times." % data_name)
                        attr[data_name] = new_data
                        del new_data

                    else:
                        raise ValueError(
                            "data item type %s not implemented yet."
                            % data_item_type[0])

            else:
                raise ValueError("Unexpected line: %s" % " ".join(part_header))
