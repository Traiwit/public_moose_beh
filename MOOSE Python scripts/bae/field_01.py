"""Module that provides means to store field data.

@Note: This module provides dynamically created classes / types for storing the
field data. They have some issues. See the documentation of the base class
L{Field<bae.field_01.Field>}.

"""

__version__ = "1.11"

_version_history_ = r"""
1.0 GP: unversioned versions
1.1 GP fixed: interpolation yields None if element or node not in the field
1.2 GP added: FieldPosStructuredPoints
       changed: you can now create different field classes with the same name.
              Each time you call Field.classFromPosType() you get a new class.
1.3 GP/AF added: createFieldClass() synonym for Field.classFromPosType()
1.4 GP added: FieldPosMeshElemIP.getNodeField(), now uses bae.log_01
1.5 GP added: datatype str (useful for vtk_02.FramesFieldsOnUnstructuredPoints
          (csv files) and possibly not much more.
1.6 GP added: interpolate methods that use bae.mesh_01.InterpolMeshToPoints
       instances for interpolation.
1.7 GP added pickle support. Made createFieldClass the initial function
       definition and Field.classFromPosType assigned to it.
       Added createFieldObject
1.8 GP added eigenVV
1.9 GP added Field.fromComponents and componentSuffixes to vector and tensor
       classes
1.10 GP added FieldPosStructuredPoints.interpolateToGridPwc
1.11 GP added FieldPosStructuredPoints.fromBoxMask,
     make FieldPosStructuredPoints.interpolateToPoints() accept a
     MeshStructuredPoints-objects as point parameter and process it efficiently
"""

import cPickle as pickle
from collections import defaultdict
from itertools import izip
from array import array

from math import ceil, floor
try:
    import numpy as np
except ImportError:
    np = None

from bae.mesh_01 import MeshStructuredPoints
from bae.log_01 import msg


def createFieldClass(fieldName, position, dataType, addBaseClasses=()):
    """
    Dynamically generate a new class for holding field data.

    This is basically a convenience function for creating the new class
    from the appropriate base classes. I.e.
     >>> Field_U = Field.classFromPosType("U", "node", "vector")
     ... # or ...
     >>> Field_U = createFieldClass("U", "node", "vector")

    is basically a replacement for
     >>> class Field_U(FieldPosMeshNode, FieldTypeVector):
     >>>     fieldName = "U"

    Each time you call Field.classFromPosType() / createFieldClass() you
    get a new class. Different classes of the same name just have the same
    __name__ attribute but are different classes besides that.

    @param fieldName: identifies the subclass / type. The class name
    self.__name__ will be "Field_%s"%fieldName.
    @param position: can be "node", "element", "elemIP", "point"
    or "structPt"
    @param dataType: "scalar", "vector", "tensor", "str"
    @param addBaseClasses: A tuple of other base classes of the field class
    to be. This can be used add additional methods and attributes to this
    class. It is in fact being used for fields to be read from the odb.

    @Note: The __name__ attribute of the class will be "Field_%s"%fieldName
    regardless of the identifier you assign it to. I.e.
     >>> MyField = Field.classFromPosType("U", "node", "vector")
     ... # ... or ...
     >>> MyField = createFieldClass("U", "node", "vector")
     >>> f = MyField()
     >>> print f.__class__
     <class '__main__.Field_U'>

    whereas
     >>> class MyField(FieldPosMeshNode, FieldTypeVector):
     >>>     fieldName = "U"
     >>> f = MyField()
     >>> print f.__class__
     <class '__main__.MyField'>
    """

    classValues = dict()
    classValues["fieldName"] = fieldName
    classValues["__slots__"] = []

    baseClasses = (
        (_positionToClass[position], _typeToClass[dataType])
        + addBaseClasses)

    cls = type("Field_%s"%fieldName, baseClasses, classValues)

    return cls

def createFieldObject(
        fieldName, position, dataType, addBaseClasses=(),
        initArgs=[]):
    """Convenience function to directly create a field object of a dynamically
    generated class. In case you just want one object of the new class.
    And this is required for the pickle-magic.

    Arguments as for L{createFieldClass} plus initArgs: a tuple or list of
    arguments passed to the constructor of the new object.
    """
    cls = createFieldClass(fieldName, position, dataType, addBaseClasses)
    obj = cls(*initArgs)
    return obj


class Field(object):
    """Base class for all classes that store field information.

    The child class that is actually used is derived from dict or list and so
    the contents can be addressed using the []-operator. The key or index is
    a node number, element number, point index depending on the actual type.

    Objects of this type can be initialized with an iterable as suitable for
    the corresponding base class: For position "node", "element" and "elemIP"
    this base class is dict. For position "point" and "structPt" this base
    class is list. I.e.
     >>> Field_B = Field.classFromPosType("B", "point", "scalar")
     >>> fr1_B = Field_B([1,1,1,0,0,0,1])

    An object of this type holds the data for a specifig time point. This class
    and its subclasses are completely timeless and unaware of time.

    Usage: Create a subclass for a particular field, e.g. displacement "U"
    using the method classFromPosType(). Then create objects of this class e.g.
    for each time step:
     >>> Field_U = Field.classFromPosType("U", "node", "vector")
     >>> fr1_U = Field_U()
     >>> fr1_U[1] = [0.0, 0.0, 0.0]  # no displacement for node 1
     >>> fr1_U[7] = [0.1, 1.2, 1.8]  # displacement of node 7

    @cvar fieldName: The fieldName argument of OdbReader.getOdbFieldClass().
    @cvar position: specifies what key or index of this dict/list stands for.
    Can be "node", "element", "elemIP", "point" or "structPt".
    @cvar dataType: "scalar", "vector", "tensor", "str"

    @Note: You will usually use dynamically created classes / types based on
    class Field. Objects / instances of those types can not easily be pickled
    using the ordinary pickle.dump() method. Use the pickle_dump() and
    pickle_load() methods instead.

    @Note: Field.classFromPosType is synonymous to createFieldClass.
    I.e.
     >>> from bae.field_01 import Field
     >>> Field_U = Field.classFromPosType("U", "node", "vector")

    is identical to
     >>> from bae.field_01 import createFieldClass
     >>> Field_U = createFieldClass("U", "node", "vector")

    @Note: dataType "str" might not work in all circumstances
    """

    @staticmethod
    def classFromPosType(fieldName, position, dataType, addBaseClasses=()):
        """Dynamically generates a new class for holding field data.

        @param fieldName: Used to identify the field for example in data
        containers. The subclass / type name (self.__name__) will be
        "Field_%s"%fieldName.
        @param position: can be "node", "element", "elemIP", "point"
        or "structPt"
        @param dataType: "scalar", "vector", "tensor", "str"
        @param addBaseClasses: Optional and typically not required: A tuple of
        other base classes of the field class to be. This can be used add
        additional methods and attributes to this class.
        It is being used for fields to be read from the odb by the (deprecated)
        L{bae.odb_02} module.

        @Note: The __name__ attribute of the class will be "Field_%s"%fieldName
        regardless of the identifier you assign it to. I.e.
         >>> MyField = Field.classFromPosType("U", "node", "vector")
         >>> f = MyField()
         >>> print f.__class__
         <class '__main__.Field_U'>

        whereas
         >>> class MyField(FieldPosMeshNode, FieldTypeVector):
         >>>     fieldName = "U"
         >>> f = MyField()
         >>> print f.__class__
         <class '__main__.MyField'>

        @Note: Field.classFromPosType() is an alias for L{createFieldClass}.

        @Note: This is basically a convenience function for creating the new
        class from the appropriate base classes. I.e.
         >>> Field_U = Field.classFromPosType("U", "node", "vector")

        is basically a replacement for
         >>> class Field_U(FieldPosMeshNode, FieldTypeVector):
         >>>     fieldName = "U"

        @Note: Each time you call Field.classFromPosType() / createFieldClass()
        you get a new class. Different classes of the same name just have the
        same __name__ attribute but are different classes besides that.
        """
        return createFieldClass(fieldName, position, dataType, addBaseClasses)

    @classmethod
    def cloneType(cls, fieldName=None, position=None, dataType=None):
        """Create a copy of the type of self.
        This is the preferred way of getting a new type/class with only one or
        some class attributes changed.

        For Example to get a new field with changed name:
         >>> fld = odb.getFieldFromOdb("SDV1")
         >>> renamed = fld.cloneType(fieldName="PST")(fld)

        @Note: Mutable class attributes (e.g. lists) are only shallow-copied!
        """
        newtype = type(cls.fieldName, cls.__bases__, dict(cls.__dict__))
        if fieldName:
            newtype.fieldName = fieldName
        if position:
            newtype.position = position
        if dataType:
            newtype.dataType = dataType
        return newtype

    @classmethod
    def fromComponents(cls, fieldName, dataType, data):
        """Return a tensor field that is joined from individual components
        stored in the provided data dictionary.

        Grabs the fields (dict values) identified by the dictionary keys
        <fieldName>_11, <fieldName>_22, ... and returns an equivalent tensor
        field.

        @param dataType: must be "vector" or "tensor"
        """

        # determine component names
        if dataType=="vector":
            class DataTypeCls(FieldTypeVector):
                pass
        elif dataType=="tensor":
            class DataTypeCls(FieldTypeTensor):
                pass
        else:
            raise ValueError(
                "dataType for Field.fromComponents must be vector or tensor.")
        DataTypeCls.fieldName = fieldName
        compFieldNames = DataTypeCls.getComponentNamesList()
        del DataTypeCls

        # determine position and type of field data container (list or dict)
        firstFld = data[compFieldNames[0]]
        position = firstFld.position
        islist = isinstance(firstFld, list)
        del firstFld

        # create results class
        ResClass = Field.classFromPosType(fieldName, position, dataType)

        # copy data
        if islist:
            valIterator = izip(*(data[comp] for comp in compFieldNames))
            res = ResClass(list(xx) for xx in valIterator)
        else:
            raise NotImplementedError(
                "Field.fromComponents only implemented for position point and"
                " structPt so far.")
        return res


    def __getstate__(self):

        # possibly store the keys in state and get an iterator for the values
        if isinstance(self, list):
            state = [None]
            valIter = iter(self)
        elif isinstance(self, dict):
            valIter = self.itervalues()
            # store the keys as signed long integers
            state = [None, array("l", self.iterkeys()).tostring()]

        # serialize vectors and tensors
        if self.dataType in ("vector", "tensor"):
            valIter = (x for xx in valIter for x in xx)

        # now store the values
        if self.dataType=="str":
            state[0] = tuple(valIter)
        elif self.dataType in ("scalar", "vector", "tensor"):
            state[0] = array("d", valIter).tostring()
        else:
            raise NotImplementedError(
                "ERROR: dataType %s not implemented in Field.__getstate__.")

        return state

    def __setstate__(self, state):

        # unpack dataType str or else array of floats
        if self.dataType=="str":
            valIter = iter(state[0])
        else:
            arr = array("d")
            arr.fromstring(state[0])
            valIter = iter(arr)

        # pack vectors and tensors from flat/serial representation in array
        if self.dataType in ("vector", "tensor"):
            valIter = (list(x)
                       for x in izip(*[valIter]*self.nbComponents))

        # store data
        if isinstance(self, list):
            list.extend(self, valIter)
        elif isinstance(self, dict):
            arr = array("l")
            arr.fromstring(state[1])
            keyIter = iter(arr)
            dict.update(self, izip(keyIter, valIter))

    def __reduce__(self):
        """Needed for pickling. This is magic...

        See: https://docs.python.org/2/library/pickle.html#the-pickle-protocol
        """
        return (
            createFieldObject,    # function to (re-)create the object
            (self.fieldName, self.position, self.dataType), # arguments to above
            self.__getstate__())  # arguments for self.__setstate__()

    def pickle_dump(self, file, protocol=0):
        """pickle (export) the field to the given open file.

        Arguments are the same as for the dump function and method of the
        builtin pickle modules.

        Usage:
         >>> myField.pickle_dump(open("myfield.pickle", "wb"))
         >>> restoredField = Field.pickle_load(open("myfield.pickle", "rb"))
        """

        if isinstance(self, list):
            data = list(self)
        elif isinstance(self, dict):
            data = self.items()

        pickle.dump((self.fieldName, self.position, self.dataType),
                    file, protocol)
        pickle.dump(data, file, protocol)

    @staticmethod
    def pickle_load(file):
        """ unpickle (import) a field from the given open file.

        Usage:
         >>> myField.pickle_dump(open("myfield.pickle", "wb"))
         >>> restoredField = Field.pickle_load(open("myfield.pickle", "rb"))
        """

        # get type info and data from the pickle file
        (fieldName, position, dataType) = pickle.load(file)
        data = pickle.load(file)

        # create the new instance and fill with the data
        cls = Field.classFromPosType(fieldName, position, dataType)
        self=cls()
        if isinstance(self, list):
            list.extend(self, data)
        elif isinstance(self, dict):
            dict.update(self, data)
        return self


class FieldPosMeshNode(Field, dict):
    """field with values at the nodes"""

    position = "node"

    def interpolateToPoints(self, mesh, elemCoords):
        """Interpolate the field values onto arbitrary points within the
        corresponding mesh.

        Also note alternative self.@L{interpolate}() method. This method
        interpolateToPoints() might be considered deprecated, not sure about
        that.

        Usage:
         >>> myPoints = [[1.0,0.0,3.1], [0.7,0.0,3.5], [0.1,1.5,3.4]]
         >>> elemCoords = mesh.getElemCoords(myPoints)
         >>> values = myField.interpolateToPoints(mesh, elemCoords)
         >>> import csv   # now output to a csv file
         >>> output = csv.writer(open("myoutput.csv", "wb"))
         >>> for pt, val in zip(myPoints, values):
         >>>     output.writerow(pt+[val])  # suppose val is a scalar

        @param elemCoords: a sequence of (element number, element coordinates)
          tuples, each describing one point. This data is provided by
          Mesh.getElemCoords().

          elemCoords may contain tuples with None as element number, as
          Mesh.getElemCoords() provides for points outside the mesh. Those
          items yield None on output.

        @Returns: An iterator over all values for each of the points in
        elemCoords

        @Note: An iterator as returned by this function can only be used once.
        If you need the values several times, convert to (i.e. store in) a
        list:
         >>> values = list(myField.interpolateToPoints(mesh, elemCoords))
        """
        getWeightedAverageFunc = self.getWeightedAverage4Item
        elemSFuncs = mesh.getElemShapeFuncs(elemCoords)
        for element, shapeFuncs in elemSFuncs:
            if element is None:
                yield None
            else:
                nodes = mesh.elNodes[element]
                yield getWeightedAverageFunc(nodes, shapeFuncs)

    def interpolate(self, interpolation, externalPtValue=None):
        """Interpolate the field values onto points according to the
        specified InterpolMeshToPoints object.

        Usage:
         >>> myPoints = [[1.0,0.0,3.1], [0.7,0.0,3.5], [0.1,1.5,3.4]]
         >>> myGrid = MeshUnstructuredPoints(myPoints)
         >>> interpolation = InterpolMeshToPoints(myGrid)
         >>> values = myField.interpolate(interpolation)
         >>> import csv   # now output to a csv file
         >>> output = csv.writer(open("myoutput.csv", "wb"))
         >>> for pt, val in zip(myPoints, values):
         >>>     output.writerow(pt+[val])  # suppose val is a scalar

        @param interpolation: a L{bae.mesh_01.InterpolMeshToPoints}
          instance. interpolation.elemCoords may contain tuples with
          None as element number, as Mesh.getElemCoords() provides for points
          outside the mesh. Those items yield None on output.

        @param externalPtValue: value to be returned for points that are
          outside the domain of the interpolation. Defaults to None.

        @Returns: An iterator over all values for each of the points in
          interpolation.

        @Note: An iterator as returned by this function can only be used once.
        If you need the values several times, convert to (i.e. store in) a
        list:
         >>> values = list(myField.interpolate(interpolation))
        """
        if not hasattr(interpolation, "getElemShapeFuncs"):
            raise NotImplementedError(
                "Don't know what to do with Interpolation object of type <%s>."
                "Maybe use Field.interpolate2()? interpolation:\n%s"
                % (type(interpolation), interpolation))

        # interpolation of type InterpolMeshToPoints or similar
        # resulting in an iterator of values for single points
        getWeightedAverageFunc = self.getWeightedAverage4Item
        elNodes = interpolation.mesh.elNodes
        elemSFuncs = interpolation.getElemShapeFuncs()
        for element, shapeFuncs in elemSFuncs:
            if element is None:
                yield externalPtValue
            else:
                nodes = elNodes[element]
                yield getWeightedAverageFunc(nodes, shapeFuncs)

    def interpolate2(self, interpolation, externalPtValue=None):
        """Interpolate the field values according to the specified
        interpolation object.

        Usage:
         >>> myPoints = [[1.0,0.0,3.1], [0.7,0.0,3.5], [0.1,1.5,3.4]]
         >>> myGrid = MeshUnstructuredPoints(myPoints)
         >>> interpolation = InterpolMeshToPoints(myGrid)
         >>> fld = myField.interpolate2(interpolation)
         >>> 

        @param interpolation: Can be a L{bae.mesh_01.InterpolMeshToPoints}
          instance. Then, interpolation.elemCoords may contain tuples with
          None as element number, as Mesh.getElemCoords() provides for points
          outside the mesh. Those items yield None on output.

          Can be a L{bae.mesh_01.InterpolMeshFacesToSurface} instance. In this
          case there is actually no interpolation. The task is rather to select
          only the values for the required nodes.

        @param externalPtValue: value to be returned for points that are
          outside the domain of the interpolation. Defaults to None.

        @Returns: A new L{Field} instance.

        @Note: If interpolation is a L{bae.mesh_01.InterpolMeshFacesToSurface}
          instance then a L{Field} object is returned that might just be a
          subset of self. The values are shallow-copied from the source, vector
          values will refer to the very same list of component values as in
          self.
        """
        if hasattr(interpolation, "getElemShapeFuncs"):
            # interpolation of type InterpolMeshToPoints or similar
            getWeightedAverageFunc = self.getWeightedAverage4Item
            elNodes = interpolation.mesh.elNodes
            elemSFuncs = interpolation.getElemShapeFuncs()
            if (hasattr(interpolation, "grid")
                and isinstance(interpolation.grid, MeshStructuredPoints)):
                position = "structPt"
            else:
                position = "point"
            FldType = createFieldClass(self.fieldName, position, self.dataType)
            res = FldType(
                externalPtValue if element is None
                else getWeightedAverageFunc(elNodes[element], shapeFuncs)
                for element, shapeFuncs in elemSFuncs)
        
        elif hasattr(interpolation, "targetTopo"):
            # interpolation of type InterpolMeshFacesToSurface
            # resulting in an object of type FieldPosMeshNode...
            targetNodes = interpolation.targetTopo.getConnectedNodes()
            FldType = type(self)
            res = FldType(
                (key, val)
                for key, val in self.iteritems()
                if key in targetNodes)

        else:
            raise NotImplementedError(
                "Don't know what to do with Interpolation object of type <%s>:"
                "\n%s" % (type(interpolation), interpolation))

        return res


class FieldPosMeshElement(Field, dict):
    """Field with only one value per element.
    """

    position = "element"

    def interpolateToPoints(self, mesh, elemCoords):
        """Interpolate the field values onto arbitrary points within the
        corresponding mesh. See the description of the L{corresponding method
        for FieldPosMeshNode<bae.field_01.FieldPosMeshNode.interpolateToPoints>}.

        Also note alternative self.@L{interpolate}() method. This method
        interpolateToPoints() might be considered deprecated, not sure about
        that.
        """
        for element, shapeFuncs in elemCoords:
            if element is None:
                yield None
            else:
                yield self[element]

    def interpolate(self, interpolation, externalPtValue=None):
        """Interpolate the field values onto arbitrary points according to the
        specified InterpolMeshToPoints object. See the description of the
        L{corresponding method for
        FieldPosMeshNode<bae.field_01.FieldPosMeshNode.interpolate>}.
        """
        for element, shapeFuncs in interpolation.elemCoords:
            if element is None:
                yield externalPtValue
            else:
                yield self[element]

    def interpolate2(self, interpolation, externalPtValue=None):
        """Interpolate the field values according to the specified
        interpolation object.

        See the description of the L{corresponding method for
        FieldPosMeshNode<bae.field_01.FieldPosMeshNode.interpolate>}.
        """
        if hasattr(interpolation, "elemCoords"):
            # interpolation of type InterpolMeshToPoints or similar
            if (hasattr(interpolation, "grid")
                and isinstance(interpolation.grid, MeshStructuredPoints)):
                position = "structPt"
            else:
                position = "point"
            FldType = createFieldClass(self.fieldName, position, self.dataType)
            res = FldType(
                externalPtValue if element is None
                else self[element]
                for element, shapeFuncs in interpolation.elemCoords)

        elif hasattr(interpolation, "volRef"):
            # interpolation of type InterpolMeshFacesToSurface
            # resulting in an object of type FieldPosMeshElement
            FldType = type(self)
            res = FldType(
                (targetElem, self[sourceElem])
                for targetElem, (sourceElem, faceId, gpIW)
                in interpolation.volRef.iteritems())

        else:
            raise NotImplementedError(
                "Don't know what to do with Interpolation object of type <%s>:"
                "\n%s" % (type(interpolation), interpolation))

        return res


class FieldPosMeshElemIP(Field, dict):
    """field with values at the integration points"""
    position = "elemIP"

    def interpolateToPoints(self, mesh, elemCoords):
        """Interpolate the field values onto arbitrary points within the
        corresponding mesh. See the description of the L{corresponding method
        for FieldPosMeshNode<bae.field_01.FieldPosMeshNode.interpolateToPoints>}.

        This interpolation is linear inside each element and non-contiuous
        across element boundaries.

        Also note alternative self.@L{interpolate}() method. This method
        interpolateToPoints() might be considered deprecated, not sure about
        that.
        """
        getWeightedAverageFunc = self.getWeightedAverageIP4Elem
        elemIPFuncs = mesh.getElemIPInterpolFuncs(elemCoords)
        for element, weights in elemIPFuncs:
            if element is None:
                yield None
            else:
                yield getWeightedAverageFunc(element, weights)

    def interpolate(self, interpolation, externalPtValue=None):
        """Interpolate the field values onto arbitrary points according to the
        specified InterpolMeshToPoints object. See the description of the
        L{corresponding method for
        FieldPosMeshNode<bae.field_01.FieldPosMeshNode.interpolate>}.

        This interpolation is linear inside each element and non-contiuous
        across element boundaries.
        """
        getWeightedAverageFunc = self.getWeightedAverageIP4Elem
        elemIPFuncs = interpolation.getElemIPInterpolFuncs()
        for element, weights in elemIPFuncs:
            if element is None:
                yield externalPtValue
            else:
                yield getWeightedAverageFunc(element, weights)


    def interpolate2(self, interpolation, externalPtValue=None):
        """Interpolate the field values according to the specified
        interpolation object.

        See the description of the L{corresponding method for
        FieldPosMeshNode<bae.field_01.FieldPosMeshNode.interpolate>}.
        """
        getWeightedAverageFunc = self.getWeightedAverageIP4Elem  # abbreviation
        if hasattr(interpolation, "getElemIPInterpolFuncs"):
            # interpolation of type InterpolMeshToPoints or similar
            elemIPFuncs = interpolation.getElemIPInterpolFuncs()
            if (hasattr(interpolation, "grid")
                and isinstance(interpolation.grid, MeshStructuredPoints)):
                position = "structPt"
            else:
                position = "point"
            FldType = createFieldClass(self.fieldName, position, self.dataType)
            res = FldType(
                externalPtValue if element is None
                else getWeightedAverageFunc(element, weights)
                for element, weights in elemIPFuncs)

        elif hasattr(interpolation, "volRef"):
            # interpolation of type InterpolMeshFacesToSurface
            # resulting in an object of type FieldPosMeshElemIP
            FldType = type(self)
            res = FldType(
                (targetElem,
                 [getWeightedAverageFunc(sourceElem, w) for w in weights_all])
                for targetElem, (sourceElem, faceId, weights_all)
                in interpolation.volRef.iteritems())

        else:
            raise NotImplementedError(
                "Don't know what to do with Interpolation object of type <%s>:"
                "\n%s" % (type(interpolation), interpolation))

        return res

    def interpolateToPointsPiecewiseConst(self, mesh, elemCoords):
        """Interpolate the field values onto arbitrary points within the
        corresponding mesh. See the description of the L{corresponding method
        for FieldPosMeshNode<bae.field_01.FieldPosMeshNode.interpolateToPoints>}.

        This interpolation is piecewise constant. In general it is
        non-contiuous. It yields the value of the N-th integration point where
        N is the index of the largest element coordinate / closest corner.

        Also note alternative self.@L{interpolatePiecewiseConst}() method.
        This method interpolateToPointsPiecewiseConst() might be considered
        deprecated, not sure about that.
        """
        shapeNotImplemented = defaultdict(int)
        for element, coords in elemCoords:
            if element is None:
                yield None
            else:
                # determine element shape
                shape = self.elShape[element]
                if shape=="TET_L":
                    idx = 0
                    yield self[element][idx]
                elif shape=="TET_Q":
                    # index of the largest element coordinate / closest corner
                    idx = coords.index(max(coords))
                    yield self[element][idx]
                else:
                    shapeNotImplemented[shape] += 1
                    yield None

        # output possible error message
        for shape, cnt in shapeNotImplemented.iteritems():
            msg("WARNING: Element shape %s is not implemented in"
                " FieldPosMeshElemIP.interpolateToPointsPiecewiseConst()."
                " This affects %d elements.\n" % (shape, cnt))
        return

    def interpolatePiecewiseConst(self, interpolation):
        """Interpolate the field values onto arbitrary points according to the
        specified InterpolMeshToPoints object. See the description of the
        L{corresponding method for
        FieldPosMeshNode<bae.field_01.FieldPosMeshNode.interpolate>}.

        This interpolation is piecewise constant. In general it is
        non-contiuous. It yields the value of the N-th integration point where
        N is the index of the largest element coordinate / closest corner.
        """
        shapeNotImplemented = defaultdict(int)
        for element, coords in interpolation.elemCoords:
            if element is None:
                yield None
            else:
                # determine element shape
                shape = self.elShape[element]
                if shape=="TET_L":
                    idx = 0
                    yield self[element][idx]
                elif shape=="TET_Q":
                    # index of the largest element coordinate / closest corner
                    # Note: in pure python this method is faster than:
                    # max(enumerate(coords), key=lambda x:x[1])[0]
                    # max(ids, key=coords.__getitem__)
                    # max((c,i) for i,c in enumerate(coords))[1]
                    # max((c,i) for c,i in izip(coords, ids))[1]
                    idx = coords.index(max(coords))
                    yield self[element][idx]
                else:
                    shapeNotImplemented[shape] += 1
                    yield None

        # output possible error message
        for shape, cnt in shapeNotImplemented.iteritems():
            msg("WARNING: Element shape %s is not implemented in"
                " FieldPosMeshElemIP.interpolatePiecewiseConst()."
                " This affects %d elements.\n" % (shape, cnt))
        return

    def getElementMeanField(self):
        """Return a field object with position="element" that contains the mean
        of the values at the integration points. The fieldName of the new field
        class will be the old fieldName plus "_EM" for element mean."""
        fieldName = "%s_EM" % self.fieldName
        FieldClass = Field.classFromPosType(
            fieldName, "element", self.dataType)
        getMeanFunc = self.getMeanIP4Values
        newField = FieldClass(
            (elem, getMeanFunc(valueIP))
            for elem, valueIP in self.iteritems())
        return newField

    def getNodeField(self, nodeElIpWt, OutputFieldClass=None):
        """Return a field object with position="node" that contains the field
        continuously extrapolated to the nodes. The values at the integration
        points are being multiplied with the corresponding weight factor and
        the sum then yields the nodal value.

        @param nodeElIpWt: Data structure containing weight factors supplied by
           L{bae.mesh_01.Mesh.getNodalAveragingIPCoeffs()} of the form
           { nodelabel : [ [elemlabel, [weight for IP1, weight for IP2, ....]],
           [elemLabel, [weights ...]], ],   ... }

        @param OutputFieldClass: If you call this function multiple times for
           the same kind of field (same field from different time frames for
           instance) it's more economic to create the class for the resulting
           nodal field only once, for example like this:
            >>> FieldClassIP = Field.classFromPosType("U", "elemIP", "vector")
            >>> # ... get some data of type FieldClassIP (fields on mymesh)
            >>> #     and store it in frameDataList
            >>> FieldClassNodal = Field.classFromPosType(
            >>>     "U_ANV", "node", "vector")
            >>> nodeElIpWt = mymesh.getNodalAveragingIPCoeffs()
            >>> for fieldIP in frameDataList:
            >>>     fieldNodal = fieldIP.getNodeField(
            >>>         nodeElIpWt, OutputFieldClass=FieldClassNodal)
            >>>     # ... further process the nodal field

           If not specified will be derived from self.
           The fieldName of the new field class will be the old fieldName plus
           "_ANV" for average nodal value. This can also be used for subsequent
           calls to this functions like so:
            >>> nodeElIpWt = mymesh.getNodalAveragingIPCoeffs()
            >>> FieldClassNodal = None
            >>> for fieldIP in frameDataList:
            >>>     fieldNodal = fieldIP.getNodeField(
            >>>         nodeElIpWt, OutputFieldClass=FieldClassNodal)
            >>>     if FieldClassNodal is None:
            >>>         FieldClassNodal = type(fieldNodal)
            >>>     # ... further process the nodal field

        @Note: To get the effect of the "don't average on region boundaries"
        option in Abaqus/CAE call Mesh.getNodalAveragingIPCoeffs() several
        times region by region and merge the results into one dictionary to
        be supplied to this function. Well this is not the whole story as you
        would obviously need two values for one node where two regions
        intersect...
        """
        if OutputFieldClass is None:
            fieldName = "%s_ANV" % self.fieldName
            OutputFieldClass = Field.classFromPosType(
                fieldName, "node", self.dataType)
        else:
            assert OutputFieldClass.position == "node",\
                "Incompatible OutputFieldClass"
            assert OutputFieldClass.dataType == self.dataType,\
                "Incompatible OutputFieldClass"

        # this is the faster version without intermediate diagnostic output
        msg("Extrapolating integration point values to nodes.")
        field = OutputFieldClass(
            (node, self.weightedSumT(
                (value,weight)
                for elem, weightsList in elemWeightsList
                for value, weight in izip(self[elem], weightsList)) )
            for node, elemWeightsList in nodeElIpWt.iteritems()
            )
        # # this is the slower version with intermediate diagnostic output
        # field = OutputFieldClass()
        # msg("Extrapolating integration point values to nodes.")
        # ticker = MsgTicker(" Already processed %s/%d."%("%d",len(elems)))
        # for cnt, (node, elemWeightsList) in enumerate(nodeElIpWt.iteritems()):
        #     ticker.msg(cnt+1)
        #     field[node] = self.weightedSumT(
        #         (value,weight)
        #         for elem, weightsList in elemWeightsList
        #         for value, weight in izip(self[elem], weightsList))
        # del ticker

        return field


class FieldPosPoint(Field, list):
    """In contrast to the other positions this class defines a non-continuous
    field only defined at certain points. There is no interpolation."""
    position = "point"


class FieldPosStructuredPoints(Field, list):
    """field with values at points in a structured points grid"""
    position = "structPt"

    @classmethod
    def fromBoxMask(cls, mesh, box, inValue, outValue):
        """Create a field with a particular value inside a given box and
        another value outside.

        Example:
         >>> mesh = getMesh(...)  # a MeshStructuredPoints-object
         >>> box = mesh.getBoundingBox()
         >>> box.scale(0.5)
         >>> heightField = createFieldObject("height", "structPt", "scalar",
         >>>     initArgs=[(xyz[2] for xyz in mesh.getPointsIter()),])
         >>> heightInBox = type(heightField).fromBoxMask(
         >>>     mesh, box, heightField, -9999)

        @param mesh: L{MeshStructuredPoints} object for the field self
        @param box: [[xmin, ymin, zmin], [xmax, ymax, zmax]]
           or L{bae.misc_01.BoundingBox} object.
        @param inValue: float or other field defined on the same grid (mesh) to
           be taken for points in or on the border of the box
        @param outValue: float or other field defined on the same grid (mesh) to
           be taken for points outside the box
        """

        # define getInVal-function: return constant given value or list item
        if isinstance(inValue, list):
            getInVal = inValue.__getitem__
        else:
            getInVal = lambda i: inValue
        # define getOutVal-function: as above
        if isinstance(outValue, list):
            getOutVal = outValue.__getitem__
        else:
            getOutVal = lambda i: outValue

        # determine index range of the box
        ximin = [
            int(ceil(float(box[0][j] - mesh.origin[j])/mesh.spacing[j]))
            for j in 0,1,2]
        ximax = [
            int(floor(float(box[1][j] - mesh.origin[j])/mesh.spacing[j]))
            for j in 0,1,2]

        # prepare result
        def getVal(i, j, k, idx, getInVal, getOutVal):
            if (ximin[0]<=i<=ximax[0]
                and ximin[1]<=j<=ximax[1]
                and ximin[2]<=k<=ximax[2]):
                return getInVal(idx)
            else:
                return getOutVal(idx)
        res = cls(
            getVal(i,j,k, idx, getInVal, getOutVal)
            for idx, (i,j,k) in enumerate(mesh.getGridIdxIter()))
        return res

    def interpolateToPoints(self, mesh, points):
        """Interpolate the field values onto arbitrary points within the
        corresponding mesh / grid.

        Usage:
         >>> myPoints = [[1.0,0.0,3.1], [0.7,0.0,3.5], [0.1,1.5,3.4]]
         >>> values = myField.interpolateToPoints(mesh, myPoints)
         >>> import csv   # now output to a csv file
         >>> output = csv.writer(open("myoutput.csv", "wb"))
         >>> for pt, val in zip(myPoints, values):
         >>>     output.writerow(pt+[val])  # suppose val is a scalar

        @param mesh: the topography of this field, i.e. the grid
        (L{MeshStructuredPoints} object) on which the field values (self) are
        defined.
        @param points: A list of point coordinates/L{MeshUnstructuredPoints}
        or a L{MeshStructuredPoints}-object.

        @Returns: An iterator over all values for each of the points. Points
        found outside the grid yield None on output.

        @Note: An iterator as returned by this function can only be used once.
        If you need the values several times, convert to (i.e. store in) a
        list:
         >>> values = list(myField.interpolateToPoints(mesh, points))
        """
        getWeightedAverageFunc = self.getWeightedAverage4Item
        outOfBounds = len(self)+9999

        if type(points) == MeshStructuredPoints:

            otherOrigin = [(points.origin[i]-mesh.origin[i])/mesh.spacing[i]
                           for i in xrange(3) ]
            otherSpacing = [points.spacing[i]/mesh.spacing[i]
                            for i in xrange(3) ]

            # idxFloat[[x-ids], [y-ids], [z-ids]]
            idsFloat = [
                [otherOrigin[i]+j*otherSpacing[i]
                 for j in range(points.gridPtNb[i])]
                for i in range(3)]

            # three lists of (inFlag, id0, xi)-tuples, one list for each dimension
            # inFlag: point from points in mesh?
            # id0: first index of the mesh-interval point is in
            # xi: 0..1 position of point in the mesh-interval
            id0Xis = [[(xi>=0 and int(xi)<=mesh.gridPtNb[i]-2, int(xi), xi-int(xi))
                       for xi in ids]
                      for i, ids in enumerate(idsFloat)]

            # iterate over i, j, k and yield result for each point
            for kIn, k, zeta in id0Xis[2]:
                for jIn, j, eta in id0Xis[1]:
                    for iIn, i, xi in id0Xis[0]:

                        ptIdx0 = (
                            # if we're outside the mesh, the first index will
                            # cause an IndexError and getWeightedAverageFunc()
                            # will return None
                            not(iIn and jIn and kIn) and outOfBounds or
                            i*mesh.strides[0]+j*mesh.strides[1]
                            +k*mesh.strides[2])
                        ptIds = [
                    ptIdx0,
                    ptIdx0+mesh.strides[0],
                    ptIdx0                +mesh.strides[1],
                    ptIdx0+mesh.strides[0]+mesh.strides[1],
                    ptIdx0                                +mesh.strides[2],
                    ptIdx0+mesh.strides[0]                +mesh.strides[2],
                    ptIdx0                +mesh.strides[1]+mesh.strides[2],
                    ptIdx0+mesh.strides[0]+mesh.strides[1]+mesh.strides[2] ]
                        weights = [
                            (1-xi)*(1-eta)*(1-zeta),
                            xi    *(1-eta)*(1-zeta),
                            (1-xi)*eta    *(1-zeta),
                            xi    *eta    *(1-zeta),
                            (1-xi)*(1-eta)*zeta    ,
                            xi    *(1-eta)*zeta    ,
                            (1-xi)*eta    *zeta    ,
                            xi    *eta    *zeta     ]
                        yield getWeightedAverageFunc(ptIds, weights)

        elif isinstance(points, list):

            # list of point coordinates/L{MeshUnstructuredPoints}
            ptIdsWeights = mesh.getPtIdsWeights(points)
            for ptIds, weights in ptIdsWeights:
                if ptIds is None:
                    yield None
                else:
                    yield getWeightedAverageFunc(ptIds, weights)
        else:
            raise NotImplementedError(
                "Points argument of type <%s> not implemented." % type(points))

    def interpolateToPointsPiecewiseConst(self, mesh, points):
        """Interpolate the field values piecewise constant onto a second
        structured grid.

        This interpolation is piecewise constant. In general the result is
        non-contiuous. It yields the value of the closest grid point of mesh.

        @param mesh: the topography of this field, i.e. the grid
        (L{MeshStructuredPoints} object) on which the field values (self) are
        defined.
        @param points: a list of point coordinates/L{MeshUnstructuredPoints}
        or a L{MeshStructuredPoints}-object.

        @Return: another Field object.

        @Note: The case points is a MeshStructuredPoints object is quite well
        optimized and should be fast.
        @Note: The closest grid point may be arbitrarily far away. It's just
        the closest, not necessarily close.
        """
        if type(points) == MeshStructuredPoints:

            otherOrigin = [(points.origin[i]-mesh.origin[i])/mesh.spacing[i]
                           for i in xrange(3) ]
            otherSpacing = [points.spacing[i]/mesh.spacing[i]
                            for i in xrange(3) ]
            closestIds = [
                [   # index between 0..gridPtNb; int(round(...) finds closest
                    min(mesh.gridPtNb[i]-1, max(0, int(round(idxFloat))))
                    # idxFloat: position of points-pt in mesh
                    for idxFloat in (otherOrigin[i]+j*otherSpacing[i]
                                     for j in range(points.gridPtNb[i]))]
                for i in range(3)]
            result = type(self)(
                self[i*mesh.strides[0]+j*mesh.strides[1]+k*mesh.strides[2]]
                for k in closestIds[2]
                for j in closestIds[1]
                for i in closestIds[0])
            return result

        elif isinstance(points, list):

            # e.g. a MeshUnstructuredPoints object
             raise NotImplementedError(
                "Expected MeshStructuredPoints object as points argument"
                " but got <%s>." % type(points))
        else:
            raise NotImplementedError(
                "Expected MeshStructuredPoints object as points argument"
                " but got <%s>." % type(points))


_positionToClass = {
    "node": FieldPosMeshNode,
    "element": FieldPosMeshElement,
    "elemIP": FieldPosMeshElemIP,
    "point": FieldPosPoint,
    "structPt": FieldPosStructuredPoints,
    }


class FieldTypeScalar(Field):
    dataType = "scalar"

    @staticmethod
    def weightedSum(values, weights):
        """sum_i value_i*weight_i
        See also L{weightedSumT} which does the same with "transposed"
        arguments:
        weightedSum(values, weights) == weightedSumT(izip(values, weights))

        @param values: iterable of float values
        @param weights: iterable of float weights
        """
        return sum(v*w for v,w in izip(values, weights))

    @staticmethod
    def weightedSumT(valueWeights):
        """sum_i value_i*weight_i. Takes the "transposed" argument of
        L{weightedSum}:
        weightedSum(values, weights) == weightedSumT(izip(values, weights))

        @param valueWeights: an iterable of (value, weight)-tuples
        @returns: sum_i value_i*weight_i, a float
        """
        return sum(v*w for v,w in valueWeights)

    @staticmethod
    def weightedSums(valuesWeightsIter):
        """calculate the weighted sum for a list (iterable) of values:

        sum_i values[i]*weights[i]
        for values, weights in valuesWeightsIter

        @Returns: an iterator of the weighted sums.

        @Note: the sum of all weights for one value in the list should be 1:
         >>> all (sum(weights) for values, weights in valuesWeightsIter == 1)
        """
        return ( sum(v*w for v,w in izip(values, weights))
                 for values, weights in valuesWeightsIter )

    def getWeightedAverage4Item(self, labels, weights):
        """Return the weighted average of those values identified by the labels
        argument.

        Intended for interpolation. E.g. for a field defined at the nodes,
        labels contains all node numbers for a particular element and weights
        contains the values of the corresponding shape functions evaluated at a
        particular point inside the element. This function will then return the
        interpolated value at the point.

        Suitable for nodal and per element data, not for integration point data.
        """
        assert len(labels)==len(weights)
        try:
            return sum(self[i]*weight for i, weight in izip(labels, weights))
        except (KeyError, IndexError):
            return None

    def getWeightedAverageIP4Elem(self, elem, weights):
        """Return the weighted average of the values at the integration point
        of the element elem.

        Intended for interpolation. weights contains the values of the shape
        functions evaluated at a particular point inside the element. This
        function will then return the interpolated value at the point.

        Only suitable for elemIP data
        """
        assert (self.position in ("elemIP", "elemSP"))
        try:
            values = self[elem]
        except KeyError:
            return None

        assert len(values)==len(weights), (
            "Mismatch between number of values (%d) and number of weights (%s)"
            " in FieldTypeScalar.getWeightedAverageIP4Elem() for element %d."
            " This occurs when the last frame is special when the simulation"
            " crashed."
            % (len(values), len(weights), elem) )
        return sum(val*weight for val, weight in izip(values, weights))

    def getMeanIP4Elem(self, elem):
        """For the given element return the arithmetic mean of all integration
        point values.

        Only suitable for elemIP data"""
        assert (self.position in ("elemIP", "elemSP"))
        try:
            values = self[elem]
        except KeyError:
            return None

        denom = 1.0/len(values)
        return sum(values)*denom

    def getMeanIP4Values(self, values):
        """Return the arithmetic mean of the given integration point values.

        Only suitable for elemIP data"""
        assert (self.position in ("elemIP", "elemSP"))
        denom = 1.0/len(values)
        return sum(values)*denom

    @classmethod
    def getComponentNamesList(cls):
        """Return a list with just one item: self.fieldName
        """
        return [cls.fieldName]


class FieldTypeVector(Field):
    """each value of the dict/list is a list of three vector components"""
    dataType = "vector"
    nbComponents = 3
    componentSuffixes = ["1", "2", "3"]

    @staticmethod
    def _vecadd(r, (v,w)):
        r[0] += v[0]*w
        r[1] += v[1]*w
        r[2] += v[2]*w
        return r

    @staticmethod
    def weightedSum(values, weights):
        """sum_i value_i*weight_i, value_i being a vector.
        See also L{weightedSumT} which does the same with "transposed"
        arguments:
        weightedSum(values, weights) == weightedSumT(izip(values, weights))

        @param values: iterable of vector values
        @param weights: iterable of float weights
        """
        return reduce(FieldTypeVector._vecadd,
                      izip(values,weights),
                      [0.0, 0.0, 0.0])

    @staticmethod
    def weightedSumT(valueWeights):
        """sum_i value_i*weight_i, value_i being a vector.
        Takes the "transposed" argument of L{weightedSum}:
        weightedSum(values, weights) == weightedSumT(izip(values, weights))

        @param valueWeights: an iterable of (value, weight)-tuples
        @returns: sum_i value_i*weight_i, a vector
        """
        return reduce(FieldTypeVector._vecadd,
                      valueWeights,
                      [0.0, 0.0, 0.0])

    @staticmethod
    def weightedSums(valuesWeightsIter):
        """calculate the weighted sum for a list (iterable) of values:

        sum_i values[i]*weights[i]
        for values, weights in valuesWeightsIter

        @Returns: an iterator of the weighted sums.

        @Note: the sum of all weights for one value in the list should be 1:
         >>> all (sum(weights) for values, weights in valuesWeightsIter == 1)
        """

        return ( reduce(FieldTypeVector._vecadd,
                        izip(values,weights),
                        [0.0, 0.0, 0.0])
                 for values, weights in valuesWeightsIter )

    def getWeightedAverage4Item(self, labels, weights):
        """Return the weighted average of those values identified by the lables
        argument.

        Intended for interpolation. E.g. for a field defined at the nodes,
        labels contains all node numbers for a particular element and weights
        contains the values of the corresponding shape functions evaluated at a
        particular point inside the element. This function will then return the
        interpolated value at the point.

        Suitable for nodal and per element data, not for integration point data.
        """
        assert len(labels)==len(weights)
        try:
            return [sum(self[i][j]*weight
                        for i, weight in izip(labels, weights))
                    for j in xrange(3)]
        except KeyError:
            return None


    def getWeightedAverageIP4Elem(self, elem, weights):
        """Return the weighted average of the values at the integration point
        of the element elem.

        Intended for interpolation. weights contains the values of the shape
        functions evaluated at a particular point inside the element. This
        function will then return the interpolated value at the point.

        Only suitable for elemIP data
        """
        assert (self.position in ("elemIP", "elemSP"))
        try:
            values = self[elem]
        except KeyError:
            return None

        assert len(values)==len(weights)
        return [sum(val[j]*weight
                    for val, weight in izip(values, weights))
                for j in xrange(3)]

    def getMeanIP4Elem(self, elem):
        """For the given element return the arithmetic mean of all integration
        point values.

        Only suitable for elemIP data"""
        assert (self.position in ("elemIP", "elemSP"))
        try:
            values = self[elem]
        except KeyError:
            return None

        denom = 1.0/len(values)
        return [sum(val[j] for val in values)*denom
                for j in xrange(3)]

    def getMeanIP4Values(self, values):
        """Return the arithmetic mean of the given integration point values.

        Only suitable for elemIP data"""
        assert (self.position in ("elemIP", "elemSP"))
        denom = 1.0/len(values)
        return [sum(val[j] for val in values)*denom
                for j in xrange(3)]

    @classmethod
    def getComponentNamesList(cls):
        """Return a list with the names of the field components. That will be
        self.fieldName plus an index 1 to 3
        """
        name = cls.fieldName
        return ["_".join((name, sfx)) for sfx in cls.componentSuffixes]


class FieldTypeTensor(Field):
    """each value of the dict/list is a list of six tensor components in the
    following order: 11, 22, 33, 12, 13, 23"""
    dataType = "tensor"
    nbComponents = 6
    componentSuffixes = ["11", "22", "33", "12", "13", "23"]

    @staticmethod
    def _vecadd(r, (v,w)):
        r[0] += v[0]*w
        r[1] += v[1]*w
        r[2] += v[2]*w
        r[3] += v[3]*w
        r[4] += v[4]*w
        r[5] += v[5]*w
        return r

    @staticmethod
    def weightedSum(values, weights):
        """sum_i value_i*weight_i, value_i being a tensor.
        See also L{weightedSumT} which does the same with "transposed"
        arguments:
        weightedSum(values, weights) == weightedSumT(izip(values, weights))

        @param values: iterable of tensor values
        @param weights: iterable of float weights
        """
        return reduce(FieldTypeVector._vecadd,
                      izip(values,weights),
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    @staticmethod
    def weightedSumT(valueWeights):
        """sum_i value_i*weight_i, value_i being a tensor.
        Takes the "transposed" argument of L{weightedSum}:
        weightedSum(values, weights) == weightedSumT(izip(values, weights))

        @param valueWeights: an iterable of (value, weight)-tuples
        @returns: sum_i value_i*weight_i, a vector
        """
        return reduce(FieldTypeVector._vecadd,
                      valueWeights,
                      [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    @staticmethod
    def weightedSums(valuesWeightsIter):
        """calculate the weighted sum for a list (iterable) of values:

        sum_i values[i]*weights[i]
        for values, weights in valuesWeightsIter

        @Returns: an iterator of the weighted sums.

        @Note: the sum of all weights for one value in the list should be 1:
         >>> all (sum(weights) for values, weights in valuesWeightsIter == 1)
        """

        return ( reduce(FieldTypeVector._vecadd,
                        izip(values,weights),
                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
                 for values, weights in valuesWeightsIter )

    def getWeightedAverage4Item(self, labels, weights):
        """Return the weighted average of those values identified by the lables
        argument.

        Intended for interpolation. E.g. for a field defined at the nodes,
        labels contains all node numbers for a particular element and weights
        contains the values of the corresponding shape functions evaluated at a
        particular point inside the element. This function will then return the
        interpolated value at the point.

        Suitable for nodal and per element data, not for integration point data.
        """
        assert len(labels)==len(weights)
        try:
            return [sum(self[i][j]*weight
                        for i, weight in izip(labels, weights))
                    for j in xrange(6)]
        except KeyError:
            return None


    def getWeightedAverageIP4Elem(self, elem, weights):
        """Return the weighted average of the values at the integration point
        of the element elem.

        Intended for interpolation. weights contains the values of the shape
        functions evaluated at a particular point inside the element. This
        function will then return the interpolated value at the point.

        Only suitable for elemIP data
        """
        assert (self.position in ("elemIP", "elemSP"))
        try:
            values = self[elem]
        except KeyError:
            return None

        assert len(values)==len(weights)
        return [sum(val[j]*weight
                    for val, weight in izip(values, weights))
                for j in xrange(6)]

    def getMeanIP4Elem(self, elem):
        """For the given element return the arithmetic mean of all integration
        point values.

        Only suitable for elemIP data"""
        assert (self.position in ("elemIP", "elemSP"))
        try:
            values = self[elem]
        except KeyError:
            return None

        denom = 1.0/len(values)
        return [sum(val[j] for val in values)*denom
                for j in xrange(6)]

    def getMeanIP4Values(self, values):
        """Return the arithmetic mean of the given integration point values.

        Only suitable for elemIP data"""
        assert (self.position in ("elemIP", "elemSP"))
        denom = 1.0/len(values)
        return [sum(val[j] for val in values)*denom
                for j in xrange(6)]

    @classmethod
    def getComponentNamesList(cls):
        """Return a list with the names of the field components. That will be
        self.fieldName plus an index 1 to 6
        """
        name = cls.fieldName
        return ["_".join((name, sfx)) for sfx in cls.componentSuffixes]


class FieldTypeString(Field):
    dataType = "str"

    @classmethod
    def getComponentNamesList(cls):
        """Return a list with just one item: self.fieldName
        """
        return [cls.fieldName]



_typeToClass = {
    "scalar": FieldTypeScalar,
    "vector": FieldTypeVector,
    "tensor": FieldTypeTensor,
    "str": FieldTypeString,
    }


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

def eigenVV(tensorFld, convertStress=False):
    """Calculate eigenvalues and eigenvectors of the given (symetrical) tensor
    field.

    The tensor is a list of six tensor components in the following order:
    11, 22, 33, 12, 13, 23

    @param convertStress: If True then assume a stress vector that is to be
    converted to geotechnical convention. I.e. compressive stress is pos.

    @returns: list of eigenvalues, list of unit eigenvectors. Largest eigenvalue
    first.
    """

    if not np or "eigh" not in vars(np.linalg):
        raise ImportError(
            "eigenVV works with numpy.linalg.eigh. Doesn't seem to be"
            " available.")

    position = tensorFld.position
    if position not in ("point", "structPt"):
        raise NotImplementedError(
            "eigenVV currently only accepts fields of position point or"
            " structPt. Instead got %s" % position)

    N = len(tensorFld)
    msg("Calculating Eigenvalues and -vectors of %d tensors." % N)

    M = N
    m = 1
    while 1:
        try:
            mat = np.zeros((M,3,3), float)
            break
        except MemoryError:
            m *= 2
            M = int(ceil(N / m)+0.1)  # +0.1 for safety
            pass
    if m>1:
        # make it even smaller for safety reasons
        m *= 2
        M = int(ceil(N / m)+0.1)  # +0.1 for safety
    msg("Assigned space for N=%d. Will need %d iterations." % (M, m))

    # preparing result fields
    eigVals = [createFieldObject("EVal%s" % i, position, "scalar")
               for i in "123"]
    eigVecs = [createFieldObject("EVec%s" % i, position, "vector")
               for i in "123"]

    # components in the tensor are: 11, 22, 33, 12, 13, 23
    ids = np.array([[0, 3, 4],
                    [3, 1, 5],
                    [4, 5, 2]])
    for i in range(m):
        if m>1:
            msg("Processing iteration %d/%d..." % (i+1,m))

        # possibly smaller vector on last iteration
        msg("Fetching data for %d points" % M)
        if (i+1)*M>N:
            MM = M
            M = N-(i*M)
            mat[:,:,:] = np.array(tensorFld[i*MM:i*MM+M], float)[:,ids]
        else:
            mat[:,:,:] = np.array(tensorFld[i*M:(i+1)*M], float)[:,ids]

        msg("computing eigenvalues and -vectors.")
        # v[..., :, i] is the normalized eigenvector corresponding to the
        # eigenvalue w[..., i]
        if np.__version__ >= "1.8":
            w,v = np.linalg.eigh(mat)
        else:
            w = np.zeros((M, 3), float)
            v = np.zeros((M, 3, 3), float)
            for j in range(M):
                w[j,:], v[j,:,:] = np.linalg.eigh(mat[j])

        # sorting eigenvalues
        msg("sorting eigenvalues and vectors")
        ids0 = np.arange(M).reshape((M,1))
        ids1 = np.arange(3).reshape((1,3,1))
        if convertStress:
            ids2 = np.argsort(w, axis=-1)
        else:
            ids2 = np.argsort(-w, axis=-1)  # sort in descending order
        # also convert to float (not: double precision)
        w = w[ids0,ids2]
        v = v[ids0[:,np.newaxis,:],ids1,ids2[:,np.newaxis,:]]
        # w = w[ids0,ids2].astype("f")
        # v = v[ids0[:,np.newaxis,:],ids1,ids2[:,np.newaxis,:]].astype("f")

        # normalize eigenvectors
        msg("normalizing %d eigenvectors" % v.shape[0])
        v /= np.sqrt((v**2).sum(axis=1))[:,np.newaxis,:]

        for i in range(3):
            if convertStress:
                eigVals[i].extend(-w[:,i])
            else:
                eigVals[i].extend(w[:,i])
            eigVecs[i].extend(v[:,:,i])
    del mat
    return eigVals, eigVecs
