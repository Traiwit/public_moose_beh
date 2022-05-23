r"""Module to read field values and other data from an Abaqus output database
(odb).

Provides class L{OdbReaderInterpolateToPoints} for reading values
interpolated to arbitrary positions.
"""

### IMPORTANT NOTE:
# The versioning information is stored in bae/odb_03/__init__.py

from array import array

from bae.mesh_01 import \
    getMesh, MeshBaseType, \
    InterpolMeshToPoints, InterpolMeshToCentroids, InterpolMeshToFaceCentroids
from bae.field_01 import createFieldClass
from bae.log_01 import msg

from bae.odb_04 import defaultNotDefinedKey, defaultNotDefinedValue
from bae.odb_04.odbReader import OdbReader

## ------------------------------------------------------------------------
#{ odb reader classes


class OdbReaderInterpolateToPoints(object):
    """For getting data at arbitrary points from fields in the odb.

    Usage:
     >>> from bae.odb_04 import OdbReader, OdbReaderInterpolateToPoints
     >>> from bae.mesh_01 import getMesh
     >>>
     >>> gridData = dict(firstPoint=[...], lastPoint=[...], ...)
     >>> grid = getMesh(**gridData)
     >>>
     >>> odbReader = OdbReader(odbPath, version="6.12-2")
     >>> odbReader.openOdb()
     >>> odb = OdbReaderInterpolateToPoints(odbReader)
     >>> odb.initOutputPoints(grid, interpolPickleName="....interpol")
     >>> fldU = odb.getFieldFromOdb(odbStep, odbFrame, "U_1")

    An odbAccessServer is running in a separate child process that opens the
    odb and does all the communication with the Abaqus API.
    This object is a like a proxy that uses the odbAccessServer to read data
    from the odb and then further processes it.

    @ivar odbReader: A L{OdbReader} instance.

    @ivar interpolation: Depending on whether L{initOutputPoints},
    L{initOutputAtCentroids} or L{initOutputAtFaceCentroids} has been called
    it's of type
    L{InterpolMeshToPoints<bae.mesh_01.InterpolMeshToPoints>},
    L{InterpolMeshToCentroids<bae.mesh_01.InterpolMeshToCentroids>} or
    L{InterpolMeshToFaceCentroids<bae.mesh_01.InterpolMeshToFaceCentroids>}.

    @ivar elemList: only if L{initOutputAtCentroids} has been called this
    is a reference to self.interpolation.elemList supplied for convenience.
    It's a (sorted) list of all element labels.
    """

    # conversion for interpolationType arguments string to integer
    _interpolationTypes = {
        "default": 1,
        "const": 2,
        }

    def __init__(self, odbReader):
        """
        @param odbReader: An L{OdbReader} object. Remember to also
        L{open<OdbReader.openOdb>} the odb before accessing it through
        methods of this class, i.e. before L{getFieldFromOdb} or
        L{getMatFieldFromOdb}
        """
        if not isinstance(odbReader, OdbReader):
            raise ValueError(
                "The odb argument to OdbReaderInterpolateToPoints must be of"
                " type OdbReader. Instead we got %s." % (type(odbReader)))
        self.odbReader = odbReader

    def _getAbqModel(self, recognizedBlocks):
        """Read the mesh from the odb.

        If the odb is not open, yet then open it and after getting the
        mesh close it again.

        This is a service function for the various initOutput functions
        that sometimes --but not always-- require the mesh from the odb.
        At that early stage the odb might still be considered closed by the
        calling process because the desired frame is not yet available in the
        odb.
        """
        if self.odbReader.odbIsOpen:
            return self.odbReader.getAbqModel(recognizedBlocks)
        else:
            msg("OdbReaderInterpolateToPoints.initOutputPoints needs the"
                " mesh from the odb which is not opened, yet. Opening it"
                " now.", debugLevel=5)
            self.odbReader.openOdb()
            model = self.odbReader.getAbqModel(recognizedBlocks)
            self.odbReader.closeOdb()
            return model

    def initOutputPoints(self, *args, **kwargs):
        """Initialize output points. Calculate the element coordinates for each
        output point and possibly store them for subsequent use in a pickle
        file. Or load this data from the so created pickle file.

        Creates the instance variable self.interpolation, an
        L{InterpolMeshToPoints<bae.mesh_01.InterpolMeshToPoints>}-object.

        @kwarg grid: A
          L{MeshStructuredPoints<bae.mesh_01.MeshStructuredPoints>} or
          L{MeshUnstructuredPoints<bae.mesh_01.MeshUnstructuredPoints>}
          object identifying the output points for the interpolation. This is
          optional if interpolPickleName is specified and the grid can be
          loaded from the corresponding file.

          If grid is not supplied then args and all kwargs except
          interpolPickleName are passed to L{bae.mesh_01.getMesh} to create
          the grid.

          If the grid is not supplied in any way it must be stored in the
          pickle file given by interpolPickleName.

        @kwarg interpolPickleName: A file containing the position of the
          grid points in the mesh (i.e. element numbers and
          element coordinates) will be saved and used on subsequent
          invocations.

          This is the file name base without the ".interpol" extension. See the
          pickleName argument of L{bae.mesh_01.InterpolMeshToPoints} for further
          details.

          If the mesh or the requested points have changed remove
          this file. If you don't specify this argument the point position will
          be calculated and no pickle file will be generated for subsequent
          invocations.
        """

        # determine interpolPickleName
        try:
            interpolPickleName = kwargs["interpolPickleName"]
        except KeyError:
            msg("OdbReaderInterpolateToPoints.initOutputPoints: Not using"
                " pickle file for the interpolation.", debugLevel=20)
            interpolPickleName = None
        else:
            msg("OdbReaderInterpolateToPoints.initOutputPoints: Using pickle"
                " file %s" % interpolPickleName, debugLevel=20)
            del kwargs["interpolPickleName"]

        # determine grid
        if "grid" in kwargs:
            # keyword argument grid
            grid = kwargs["grid"]
        elif args and isinstance(args[0], MeshBaseType):
            # first positional argument is a grid
            grid = args[0]
        elif args or kwargs:
            # create new grid from arguments
            grid = getMesh(*args, **kwargs)
        else:
            # no grid found
            grid = None

        # provide interpolation data to the remap kernel
        self.interpolation = InterpolMeshToPoints(
            grid=grid,
            getMeshCallBack=self._getAbqModel,
            gmCallBackArgs=dict(recognizedBlocks="NODE,ELEMENT"),
            pickleName=interpolPickleName
            )
        # the following chunk seems to be too expensive in normal operations
        # please uncomment if you need it for debugging
        # msg("OdbReaderInterpolateToPoints.initOutputPoints: %d points, found"
        #     " in following elements:\n%s"
        #     % (len(self.interpolation),
        #        [x for x,y in self.interpolation.elemCoords]),
        #     debugLevel=20)

        # store remap weights
        (nodalweights, gaussPtsweights) = self.interpolation.writeRemapParam(
            notDefinedKey=defaultNotDefinedKey)
        self.odbReader.comm.storeRemapWeights(
            len(self.interpolation), nodalweights, gaussPtsweights,
            defaultNotDefinedKey)
        msg("OdbReaderInterpolateToPoints.initOutputPoints:"
            " stored weights for %d points." % self.odbReader.comm.nbPoints,
            debugLevel=10)

        # prepare filtersets for the odbAccess-server
        self.setFilterElset(set(e for e, xyz in self.interpolation.elemCoords
                                if e is not None))

    def initOutputAtCentroids(self, model=None, elems=None, boundingBox=None):
        """Initialize output at element centroids.

        Creates the instance variable self.interpolation, an
        L{InterpolMeshToCentroids<bae.mesh_01.InterpolMeshToCentroids>}-object.

        And creates the instance variable elemList: a (sorted) list of all
        element labels as given by the parameter elems.

        Fields returned by subsequent calls to L{getFieldFromOdb} or
        L{getMatFieldFromOdb} are of type L{bae.field_01.Field} with attribute
        position="point". This is essentially an ordered list of values.
        For convenience self.elemList is supplied as a reference to
        self.interpolation.elemList. self.interpolation is of type
        L{bae.mesh_01.InterpolMeshToCentroids}.

        @param model: L{bae.mesh_01.Mesh} object. If not given then it will be
           loaded from the odb.
        @param elems: Set of element labels. If model is not given or a
           L{bae.abq_model_02.Model} object then elems can also be an elset
           name or a list of elset names -anything that
           L{bae.abq_model_02.Model.getUnionSet}() accepts as items argument.
        @param boundingBox: L{BoundingBox<bae.misc_01.BoundingBox>}-object (or
           [[xmin,ymin,zmin],[xmax,ymax,zmax]]). Consider only elements whose
           centroid lies within this box. If None then all elements will be
           considered.
        """
        if not model:
            # get topographic data (node connectivity, element centroids)
            recognizedBlocks = set(("NODE", "ELEMENT"))
            if isinstance(elems, (basestring, list)):
                recognizedBlocks.add("ELSET")
            model = self._getAbqModel(recognizedBlocks=recognizedBlocks)

        # provide interpolation data to the remap kernel
        self.interpolation = InterpolMeshToCentroids(
            model=model, elems=elems, boundingBox=boundingBox)

        # store a reference of interpolation.elemList in self
        self.elemList = self.interpolation.elemList

        # store remap weights
        (nodalweights, gaussPtsweights) = self.interpolation.writeRemapParam(
            notDefinedKey=defaultNotDefinedKey)
        self.odbReader.comm.storeRemapWeights(
            len(self.interpolation), nodalweights, gaussPtsweights,
            defaultNotDefinedKey)
        msg("OdbReaderInterpolateToPoints.initOutputAtCentroids:"
            " stored weights for %d points." % self.odbReader.comm.nbPoints,
            debugLevel=10)

        # prepare filtersets for the odbAccess-server
        self.setFilterElset(self.elemList)

    def initOutputAtFaceCentroids(self, *args, **kwargs):
        """Initialize output at face centroids.

        Creates the instance variable self.interpolation, an
        L{InterpolMeshToFaceCentroids<bae.mesh_01.InterpolMeshToFaceCentroids>}
        -object.

        Usage:
         >>> odb = OdbReaderInterpolateToPoints(...)
         >>> model = odb.getAbqModel()
         >>> surf = ElemFaceSurface().updateFromModel(model, "PITSURF")
         >>> odb.initOutputAtFaceCentroids(surf)
         >>> ...

        Instead of using ElemFaceSurface the following is possible, too:
         >>> ...
         >>> surfName = "PITSURF"
         >>> faceEl = dict((faceId, model.elset[elsetName])
         >>>               for faceId, elsetName
         >>>               in model.surface[surfName].iteritems())
         >>> odb.initOutputAtFaceCentroids(model, faceEl)
         >>> ...

        Fields returned by subsequent calls to L{getFieldFromOdb} or
        L{getMatFieldFromOdb} are of type L{bae.field_01.Field} with attribute
        position="point". This is essentially an ordered list of values. To
        determine the corresponding element labels and face ids check
        self.interpolation of type L{bae.mesh_01.InterpolMeshToFaceCentroids}.

        @kwarg model: a L{bae.abq_model_02.Model}; requires faceEl argument,
           precludes surface argument. DEPRECATED, better supply surface.
        @kwarg faceEl: a "faceEl"-dictionary:
           {"S1": set((4346, 4768, 6821, ...)), ...}; requires model argument,
           precludes surface argument. DEPRECATED, better supply surface.
        @kwarg surface: a L{bae.surface_03.ElemFaceSurface} object; precludes
           model and faceEl arguments.
        @Note: model and faceEl arguments are provided for compatibility
           reasons only.
        """

        # provide interpolation data to the remap kernel
        self.interpolation = InterpolMeshToFaceCentroids(*args, **kwargs)

        # store remap weights
        (nodalweights, gaussPtsweights) = self.interpolation.writeRemapParam(
            notDefinedKey=defaultNotDefinedKey)
        self.odbReader.comm.storeRemapWeights(
            len(self.interpolation), nodalweights, gaussPtsweights,
            defaultNotDefinedKey)
        msg("OdbReaderInterpolateToPoints.initOutputAtFaceCentroids:"
            " stored weights for %d points." % self.odbReader.comm.nbPoints,
            debugLevel=10)

        # prepare filtersets for the odbAccess-server
        self.setFilterElset(set(e for e, xyz in self.interpolation.elemCoords))

    def setFilterElset(self, elems):
        """Define elset and nset filter to be switched on and off with
        L{activateFilter} and L{deactivateFilter}.

        @param elems: an iterable of element labels. A copy will be made, it's
        save to pass iterators or temporary/mutable objects.
        """
        # using sorted to make sure we have a copy
        self.filterElset = sorted(elems)
        try:  # make sure there is no self.filterElsetName
            del self.filterElsetName
        except AttributeError:
            pass
        msg("OdbReaderInterpolateToPoints.setFilterElset:"
            " createFilterElset will be supplied with (first 50): %s"
            % self.filterElset[:50], debugLevel=20)

        # converting into sorted list only for diagnostic output
        self.filterNset = sorted(
            self.interpolation.mesh.getConnectedNodes(self.filterElset))
        try:  # make sure there is no self.filterNsetName
            del self.filterNsetName
        except AttributeError:
            pass
        msg("OdbReaderInterpolateToPoints.setFilterElset:"
            " createFilterNset will be supplied with (first 50): %s"
            % self.filterNset[:50], debugLevel=20)

    def activateFilter(self):
        """Activate the filters (elset and nset) to only read the values for
        the elements and nodes required for the output points.

        On the first call this creates a nset and an elset in the odb-API. The
        names of those sets will be stored in instance variables. On subsequent
        calls the old sets will then be re-used.

        When the requested data has been read the filters should be deactivated
        by calling self.L{deactivateFilter} to avoid them staying activated
        when it's not intended.
        """
        try:
            self.odbReader.setFilterElset(self.filterElsetName)
        except AttributeError:
            self.filterElsetName = (
                self.odbReader.createFilterElset(self.filterElset))
        try:
            self.odbReader.setFilterNset(self.filterNsetName)
        except AttributeError:
            self.filterNsetName = (
                self.odbReader.createFilterNset(self.filterNset))

    def deactivateFilter(self):
        """Deactivate the filters set by L{activateFilter}
        """
        self.odbReader.setFilterElset()
        self.odbReader.setFilterNset()

    def getFieldFromOdb(self, odbStep, odbFrame, fieldName,
                        notDefinedValue=None, interpolationType="default"):
        """
        @param odbStep: String specifying the odb step. E.g. "Step-2".
        @param odbFrame: Integer specifying the odb frame number / index.

        @param fieldName: examples: "S", "V", "S_MISES", "S_MIN_PRINCIPAL",
          "U_1", "S_33", "SDV4", "SDV_PST"

        @param notDefinedValue: value to be returned for points outside of the
          mesh.
        @param interpolationType: a string:
          - "default": depending on the element type, for C3D10M: linear for
            integration point values, quadratic or nodal values
          - "const": piecewise constant: take the value of the closest node or
            integration point.

        @raises ValueError: if the specified frame is not in the odb.
        """
        # Calls two different methods on the odbAccessServer to
        # accomplish the task of getting field data from the odb.

        # First it calls L{getFieldProps} to determine some properties of the
        # field. Then it initializes the resulting field values object (a
        # L{bae.field_01.Field} object) and then calls the getFieldData.

        if self.odbReader.comm.nbPoints<=0:
            raise RuntimeError(
                "OdbReaderInterpolateToPoints.getFieldFromOdb: no output"
                " points defined. Maybe initOutputPoints has not been called?")

        self.activateFilter()

        # set frame, can raise ValueError
        assert isinstance(odbStep, basestring)
        assert isinstance(odbFrame, int)
        self.odbReader._callmethod("setFrame", odbStep, odbFrame)

        # field, invariant and component will be derived from the name
        (fieldName, position, dataType, odbFieldName, component, invariant)\
            = self.odbReader.getFieldProps(fieldName)

        # create appropriate bae.field_01.Field
        FieldClass = createFieldClass(fieldName, "point", dataType)
        field = FieldClass()

        # store notDefinedValue
        if notDefinedValue is None:
            notDefinedValue = defaultNotDefinedValue[dataType]
        if dataType=="scalar":
            notDefinedValue = [notDefinedValue]
        else:
            if dataType == "tensor":
                nbComp = 6
            elif dataType == "vector":
                nbComp = 3
            # note: nbComp not defined -> unknown dataType
            if nbComp!=len(notDefinedValue):
                raise ValueError(
                    "Incompatible 'notDefinedValue' argument to"
                    " getFieldFromOdb: expected %s with %d components"
                    " but got '%s'."
                    % (dataType, nbComp, notDefinedValue))
        notDefinedArr = array("f", notDefinedValue)

        # get the actual data from the odb
        self.odbReader._startRequest(
            "getFieldData", odbFieldName, component, invariant)
        self.odbReader.comm.updateStructPtFieldValues(
            position, field, notDefinedArr, defaultNotDefinedKey,
            self._interpolationTypes[interpolationType])
        self.odbReader._finalizeRequest()
        self.deactivateFilter()

        msg("OdbReaderInterpolateToPoints.getFieldFromOdb: read %d field"
            " values. ending." % len(field), debugLevel=10)
        return field

    def getMatFieldFromOdb(self,
                           nameFilter=None, firstMatCode=None, noMatCode=-1,
                           matNameToMatNumber=None):
        """Read material section name from the odb at grid point
        positions.

        @param nameFilter: A callable or None for material name conversion.
          For example if there are region variants that shall be treated as
          one. I.e. if SEDIMENT, SEDIMENT_CAVE and SEDIMENT_UCUT shall all be
          treated as SEDIMENT.

          Must be a function otherwise it'll be ignored. It will be called as a
          conversion function taking the original material name from the odb as
          first argument and a dummy section type "Unknown type" as second
          argument. It's assumed to return the desired converted material name.
          (The section type was intended to read something like "solid",
          "beam", "shell" but this functionality is deprecated. The second
          argument is only left there for compatibility and should be ignored.)

          Example: To take everything up to the first "_" or "-" use as
          matRegionsNameFilter:
           >>> lambda matname,secType: matname.split("_")[0].split("-")[0]

          That's equivalent to:
           >>> def convert(matname,secType):
           >>>     return matname.split("_")[0].split("-")[0]
           >>> odbToVtkGridBox(..., matRegionsNameFilter=convert, ...)

        @param firstMatCode: Material code number of the first material in
          the list of material names. Or if you prefer: Offset of the material
          code numbers in the matNumberField field and the index of the
          corresponding material in the matNamesList.

          I.e. if firstMatCode=1 then all points with assigned material number
          of 1 are of the first material type in the list matNamesList which is
          matNamesList[0]. All points having material number 5 are of material
          type matNamesList[4].

          Specify firstMatCode=0 if you want the material numbers in
          matNumberField to directly specify the index in matNumberField. The
          default value of 1 seems to be more appropriate for the usual use
          case of exporting matNumberField to a vtk file and matNamesList to a
          seperate csv file.

          If not specified then it defaults to 1 or the next after what's in
          the matNameToMatNumber argument.

        @param noMatCode: Material number to represent points that are not
          being assinged any material to (because they sit outside the mesh).

        @param matNameToMatNumber: Optional dictionary determinining the
          material name to material number assignment. Keys are the final
          material names as returned by the nameFilter function if applicable.

          Typically the keys are the uppercase material names as specified
          in the Abaqus input file of the job.

          Note: Should a material name appear that is not in the supplied
          dictionary yet it will be added.

          Note: If matNameToMatNumber is specified then the resulting
          matNamesList only contains new material names that were not
          specified in the matNameToMatNumber argument.

          Note: All keys (mat names) are silently converted to uppercase
          if needed. This is done in place, i.e. the dictionary provided
          by the calling function is modified as a side effect.

        @returns: (matNamesList, matNumberField)-tuple
          matNamesList is the list of material names.
          matNumberField is the field (i.e. a list) of material

        @Note: Only materials that are found at any of the given points will
          be listed in the resulting matNamesList. (See optional parameter
          matNameToMatNumber for an exception to this.) Material names from
          other regions of the model are not considered.
        """

        # determine or check firstMatCode
        if firstMatCode is None:
            # if firstMatCode not given then determine from matNameToMatNumber
            if matNameToMatNumber:
                firstMatCode = max(matNameToMatNumber.itervalues()) + 1
            else:
                firstMatCode = 1
        elif (matNameToMatNumber and
                firstMatCode<=max(matNameToMatNumber.itervalues())):
            # if firstMatCode is given and interferes with matNameToMatNumber
            # then issue a warning
            lastName, lastNumber = max(
                matNameToMatNumber.iteritems(), key=lambda (x,y): y)
            msg("WARNING"
                " from OdbReaderInterpolateToPoints.getMatFieldFromOdb():"
                " firstMatCode %d <= last specified mat number from"
                " matNameToMatNumber %d (%s). To prevent conflicts correcting"
                " firstMatCode to %d."
                % (firstMatCode, lastNumber, lastName, lastNumber+1))
            firstMatCode = lastNumber+1
            del lastName, lastNumber

        # check that given mat names are upper case
        if matNameToMatNumber:
            for matName, matNb in matNameToMatNumber.items():
                nameUpper = matName.upper()
                if nameUpper != matName:
                    matNameToMatNumber[nameUpper] = matNameToMatNumber[matName]
                    del matNameToMatNumber[matName]

        # get the actual data from the odb
        self.activateFilter()
        msg("OdbReaderInterpolateToPoints.getMatFieldFromOdb: starting"
            " getMatNames on the odbAccessServer.", debugLevel=10)
        self.odbReader._startRequest("getMatNames", firstMatCode, 0)
        matNumbersPrelim = list()
        self.odbReader.comm.updateStructPtFieldValuesConstInt(
            matNumbersPrelim, noMatCode, defaultNotDefinedKey)
        msg("OdbReaderInterpolateToPoints.getMatFieldFromOdb: Got interpolated"
            " matNumbersPrelim for %d pts: %s"
            % (len(matNumbersPrelim), matNumbersPrelim[:10]), debugLevel=10)

        # get matNamesList from the pipe
        (nbMatNames,) = self.odbReader._pipeRead("i")
        msg("OdbReaderInterpolateToPoints.getMatFieldFromOdb: will read"
            " %d matNames from the pipe." % nbMatNames, debugLevel=10)
        matNamesList = list(self.odbReader._pipeRead("s"*nbMatNames))
        self.odbReader._finalizeRequest()
        msg("OdbReaderInterpolateToPoints.getMatFieldFromOdb: read %d mat"
            " names." % len(matNamesList), debugLevel=10)
        self.deactivateFilter()

        # extract material region info from section names and apply nameFilter
        nextNewNumber = firstMatCode
        oldToNewMatNumber = dict()
        newMatNamesList = list()
        if matNameToMatNumber is None:
            matNameToMatNumber = dict()
        for oldNumber, matName in enumerate(matNamesList, start=firstMatCode):
            if callable(nameFilter):
                matName = nameFilter(matName, "Unknown type")
            try:
                newNumber = matNameToMatNumber[matName]
            except KeyError:
                newNumber = nextNewNumber
                nextNewNumber += 1
                matNameToMatNumber[matName] = newNumber
                newMatNamesList.append(matName)
            oldToNewMatNumber[oldNumber] = newNumber
        matNamesList = newMatNamesList
        del newMatNamesList

        # create field object converting matNumbers
        MatNumberFieldClass = createFieldClass("Material", "structPt", "scalar")
        matNumberField = MatNumberFieldClass(
            oldToNewMatNumber.get(x, noMatCode)
            for x in matNumbersPrelim)

        msg("Identified %d different new material codes for %d grid points."
            " (Not counting the %d preassigned material codes.)"
            % (len(matNamesList), len(matNumberField),
               len(matNameToMatNumber)-len(matNamesList)))
        return (matNamesList, matNumberField)

#} end of odb reader classes
## ------------------------------------------------------------------------
