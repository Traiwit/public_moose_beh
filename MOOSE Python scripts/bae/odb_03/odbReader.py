r"""Module to read field values and other data from an Abaqus output database
(odb).

Provides classes L{OdbReader} for reading nodal or Gauss point data as is and
L{OdbReaderInterpolateToPoints} for reading values interpolated to arbitrary
positions.
"""

### IMPORTANT NOTE:
# The versioning information is stored in bae/odb_03/__init__.py

import re
from array import array

from bae.mesh_01 import \
    getMesh, MeshBaseType, \
    InterpolMeshToPoints, InterpolMeshToCentroids, InterpolMeshToFaceCentroids
import bae.abq_model_02 as abq_model
from bae.field_01 import createFieldClass
from bae.log_01 import msg

from bae.odb_03 import defaultNotDefinedKey, defaultNotDefinedValue
from bae.odb_03.odbReaderBase import OdbReaderBase


## ------------------------------------------------------------------------
#{ odb reader classes

class OdbReader(OdbReaderBase):
    """For getting data the odb.

    An odbAccessServer is running in a separate child process that opens the
    odb and does all the communication with the Abaqus API.
    This object is a like a proxy that uses the odbAccessServer to read data
    from the odb and then further processes it.
    """

    def getAbqModel(self, recognizedBlocks=None):
        r"""get an L{abq_model_02.Model<bae.abq_model_02.Model>}-object from
        the odb model data,

        reads nodes, elements and elsets

        assumes there is only one part instance in the model

        ideas for the future:
         - read nsets, surfaces, all the rest

        @param recognizedBlocks: If specified only get the specified property
        from the odb. Not implemented so far. Reads "NODE", "ELEMENT", "ELSET".
        Can be specified as tuple of strings or comma separated string.

        @Note: elsets are being read from the first part instance and from the
        root assembly as well. If the same elset name occurs twice the elset
        in the root assembly is ignored.
        @Note: Look in abq_model_02.Model.read how recognizedBlocks is treated
        there and adapt this here accordingly.
        @Note: Can only handle 3D models: Node coords are always float-triples.
        """
        model = abq_model.Model()

        # preprocess the recognizedBlocks argument
        recognizedBlocks = model._checkRecognizedBlocks(recognizedBlocks)

        if "NODE" in recognizedBlocks:
            # get the nodes
            self._startRequest("getNodeCoordsFromOdb")
            self.comm.updateNodeCoordsFromPipe(model.nodeCoords)
            self._finalizeRequest()

        if "ELEMENT" in recognizedBlocks:
            # get the elements
            self._startRequest("getElementsFromOdb")
            self.comm.updateElNodesFromPipe(
                model.elNodes, model.elType, model.typeEl)
            self._finalizeRequest()

        if "ELSET" in recognizedBlocks:
            # get the elements
            self._startRequest("getElsetsFromOdb")
            while 1:
                (elsetName,) = self._pipeRead("s")
                if not elsetName:
                    break
                (nbElems,) = self._pipeRead("i")
                elset = set(self._pipeReadArray("i", nbElems))
                model.elset[elsetName] = elset
            self._finalizeRequest()

        msg("OdbReader.getAbqModel: Got model from the odb: %s" % model,
            debugLevel=10)

        return model

    def createFilterElset(self, elset=None):
        """Set a elset as filter for all subsequent read operations for
        element based fields.
        @param elset: set or other iterable of element labels (ints) Supply an
           empty list to disable the filter.
        @returns: the Abaqus-API internal set name. This set name can later be
           used for L{setFilterElset} to re-activate this filter set.
        """
        if elset:
            elset = sorted(elset)
        else:
            elset = []
        elsetName = self._callmethod(
            "createFilterElset", elset, resultTypeCodes="s")
        return elsetName

    def createFilterNset(self, nset=None):
        """Set a nset as filter for all subsequent read operations for
        nodal fields.
        @param nset: set or other iterable of node labels (ints). Supply an
           empty list to disable the filter.
        @returns: the Abaqus-API internal set name. This set name can later be
           used for L{setFilterNset} to re-activate this filter set.
        """
        if nset:
            nset = sorted(nset)
        else:
            nset = []
        nsetName = self._callmethod(
            "createFilterNset", sorted(nset), resultTypeCodes="s")
        return nsetName

    def setFilterElset(self, elsetName=""):
        """Set a elset as filter for all subsequent read operations for
        element based fields.
        @param elsetName: name of an elset already in the odb or created
          earlier. Supply an empty string to disable the filter.
        @returns: nb of elements in the filter elset; 0 if the filter is
          disabled; -1 if the specified elset name does not exist.
        """
        if not elsetName:
            elsetName = ""
        elif not isinstance(elsetName, basestring):
            raise ValueError("setFilterElset expects a string as elsetName"
                             " argument but got %s." % elsetName)
        nbItems = self._callmethod(
            "setFilterElset", elsetName, resultTypeCodes="i")
        return nbItems

    def setFilterNset(self, nsetName=""):
        """Set a nset as filter for all subsequent read operations for
        nodal fields.
        @param nsetName: name of an nset already in the odb or created
          earlier. Supply an empty string to disable the filter.
        @returns: nb of nodes in the filter nset; 0 if the filter is
          disabled; -1 if the specified nset name does not exist.
        """
        if not nsetName:
            nsetName = ""
        elif not isinstance(nsetName, basestring):
            raise ValueError("setFilterNset expects a string as nset name but"
                             " got %s." % nsetName)
            
        nbItems = self._callmethod(
            "setFilterNset", nsetName, resultTypeCodes="i")
        return nbItems

    @staticmethod
    def writeMatNamesToCsv(outputFileName, matNamesList, firstMatCode=1):
        """Convenience function to write the material number - material regions
        relations to a csv file.
        """
        import csv
        output = csv.writer(open(outputFileName, "wb"))
        matList = [ [cnt+firstMatCode, matName]
                    for cnt, matName in enumerate(matNamesList)
                    if matName]
        output.writerow(["material number", "material name"])
        output.writerows(matList)

    @property
    def postData(self):
        """A dictionary being fed from the corresponding ..._postData.py
        database.

        Usually contains some seqElsets... dictionaries like seqElsetsExcav in
        the exampe below. And a frameNames dictionary.

        Example: Suppose in frame 3 of step 2 labelled "YR2012_M02" the elsets
        "Pit_2012_02" and "Dev_2012_02" will be excavated and in frame 4
        "YR2012_M03" only "Dev_2012_03":
         >>> odb = OdbReader("/datapth/myodb.odb")
         >>> frameNames = odb.postData["frameNames"]
         >>> print frameNames[2][3:5]
         ... ['YR2012_M02','YR2012_M03']
         >>> seqSets = odb.postData["seqElsetsExcav"]
         >>> print seqSets[2][3:5]
         ... [ ('Pit_2012_02','Dev_2012_02'),(Dev_2012_03,) ]

        The data will only be read when the portData property is being accessed
        for the first time.

        Raises IOError if the XXX_postData.py file can not be accessed. XXX is
        the name of the odb without the a trailing ".odb"-suffix.
        """
        if not hasattr(self, "_postData"):
            self._postData = dict()
            if self.odbPath:
                # deduct postData-filename from odb name
                if self.odbPath[-4:].lower() == ".odb":
                    postDataFileName = self.odbPath[:-4]+"_postData.py"
                else:
                    postDataFileName = self.odbPath+"_postData.py"
                # load/parse postData from external file
                execfile(postDataFileName, self._postData)

        return self._postData

    def getFieldFromOdb(self, stepNb, frameNb, fieldName):
        """
        Example:
         >>> odb = OdbReader("myodb.odb")
         >>> odb.setFilterElset("ALL_TETS")
         >>> fld = odb.getFieldFromOdb(2, 11, "S")
         >>> print fld.intPtSectPt, fld.dataType, fld.position
         [(1,-1), (2,-1), (3,-1), (4,-1)] tensor elemIP
         >>> print fld[101234]  # get stress tensors for elem 101234
         [[1.434E6, 3.123E4, ... 6 components for Gauss pt 1...], [6 comps
         ... for Gauss pt 2], [...], [...]]
         >>> odb.setFilterElset("BEAMS")
         >>> fld = odb.getFieldFromOdb(2, 11, "S_11")
         >>> print fld.intPtSectPt, fld.dataType, fld.position
         [(1,3), (1,7), (1,11), (1,15)] scalar elemIP
         >>> print fld[9012]  # get stress component 11 for elem 9012
         [3.1415E6, 3.333E5, ... 4 components for 4 section points]

        @param stepNb: integer number of the step in the odb. If the first
           step in the odb is called "Step-3" specify stepNb=1 for this.
        @param frameNb: integer frame number
        @param fieldName: examples: "S", "V", "S_MISES", "S_MIN_PRINCIPAL",
           "U_1", "S_33", "SDV4", "SDV_PST"

        @returns: A L{Field<bae.field_01.Field>} object.

           If possible the field has an attribute intPtSectPt: a list of (Gauss
           pt number, section point number) tuples. The items in this list
           correspond to the items in each value of the returned field, see
           examples provided.

           This attribute is missing if elements of different types with
           different numbers of integration points and section points are being
           requested.

           If you need this attribute and don't get it yet then use
           L{setFilterElset} or L{createFilterElset} to constrain the set of
           elements returned by a specific call to this function.

        @raises ValueError: if the specified frame is not in the odb.

        @Note: For tensor fields of structural elements the full tensor might
           not be available and cause an error. In this case take individual
           components one by one. The reason is that Abaqus has different
           types of tensors (TENSOR_3D_FULL, TENSOR_3D_SURFACE,
           TENSOR_3D_PLANAR) that are not fully implemented in this module.
        @Note: For beam elements (B31 in A/Explicit) it's suggested to only
           consider S_11. Gero doesn't understand the rest.
        @Note: For shell elements it might be necessary to query the
           components as S_11, S_22, S_33 in order to get S_11, S_22, S_12
           due to limitations of this module.
        """
        # Calls two different methods on the odbAccessServer to
        # accomplish the task of getting field data from the odb.

        # First it calls L{getFieldProps} to determine some properties of the
        # field. Then it initializes the resulting field values object (a
        # L{bae.field_01.Field} object) and then calls the getFieldData.

        # set frame, can raise ValueError
        self._callmethod("setFrame", stepNb, frameNb)

        # field, invariant and component will be derived from the name
        (fieldName, position, dataType, odbFieldName, component, invariant)\
            = self.getFieldProps(fieldName)

        # create appropriate bae.field_01.Field
        FieldClass = createFieldClass(fieldName, position, dataType)
        field = FieldClass()

        # get the actual data from the odb
        self._startRequest("getFieldData", odbFieldName, component, invariant)
        self.comm.updateFieldValues(position, field)
        self._finalizeRequest()
        msg("OdbReader.getFieldFromOdb: read %d field values. ending."
            % len(field), debugLevel=10)
        return field

    def getMatFieldFromOdb(self,
                           nameFilter=None, firstMatCode=None,
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
                matNameToMatNumber.iteritems(), key=lambda x,y: y)
            msg("WARNING"
                " from OdbReader.getMatFieldFromOdb():"
                " firstMatCode %d <= last specified mat number from"
                " matNameToMatNumber %d (%s). Automatically added mat numbers"
                " starting at %d may conflict with the given"
                " matNameToMatNumber dict without notice."
                % (firstMatCode, lastNumber, lastName))
            del lastName, lastNumber

        # check that given mat names are upper case
        if matNameToMatNumber:
            for matName, matNb in matNameToMatNumber.items():
                nameUpper = matName.upper()
                if nameUpper != matName:
                    matNameToMatNumber[nameUpper] = matNameToMatNumber[matName]
                    del matNameToMatNumber[matName]

        # get the actual data from the odb
        msg("OdbReader.getMatFieldFromOdb: starting getMatNames on the"
            " odbAccessServer.", debugLevel=10)
        self._startRequest("getMatNames", firstMatCode, 0)
        matNumbersPrelim = dict()
        self.comm.updateFieldValuesConstInt(matNumbersPrelim)
        msg("OdbReader.getMatFieldFromOdb: Got matNumbersPrelim for %d pts."
            % len(matNumbersPrelim), debugLevel=10)

        # get matNamesList from the pipe
        (nbMatNames,) = self._pipeRead("i")
        msg("OdbReader.getMatFieldFromOdb: will read %d matNames from the"
            " pipe." % nbMatNames, debugLevel=10)
        matNamesList = list(self._pipeRead("s"*nbMatNames))
        self._finalizeRequest()
        msg("OdbReader.getMatFieldFromOdb: read %d mat numbers and %d mat"
            " names."
            % (len(matNumbersPrelim), len(matNamesList)), debugLevel=10)

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

        # create appropriate bae.field_01.Field object converting matNumbers
        MatNumberFieldClass = createFieldClass("Material", "element", "scalar")
        try:
            matNumberField = MatNumberFieldClass(
                (e, oldToNewMatNumber[x])
                for (e, x) in matNumbersPrelim.iteritems())
        except KeyError:
            raise KeyError(
                "ERROR when translating material numbers: Not all mat numbers"
                " %s have a corresponding entry in the list of mat names. This"
                " list has %d entries and starts with %d."
                % (sorted(set(matNumbersPrelim.itervalues())),
                   len(oldToNewMatNumber), firstMatCode))

        msg("Identified %d different new material codes for %d elements."
            " (Not counting the %d preassigned material codes.)"
            % (len(matNamesList), len(matNumberField),
               len(matNameToMatNumber)-len(matNamesList)))
        return (matNamesList, matNumberField)


class OdbReaderInterpolateToPoints(OdbReader):
    """For getting data at arbitrary points from fields in the odb.

    Usage:
     >>> from bae.odb_03 import OdbReaderInterpolateToPoints
     >>> from bae.mesh_01 import getMesh
     >>>
     >>> gridData = dict(firstPoint=[...], lastPoint=[...], ...)
     >>> grid = getMesh(**gridData)
     >>>
     >>> odb = OdbReaderInterpolateToPoints(odbPath, version="6.12-2")
     >>> odb.initOutputPoints(grid, interpolPickleName="....interpol")
     >>> fldU = odb.getFieldFromOdb(stepNb, frameNb, "U_1")

    An odbAccessServer is running in a separate child process that opens the
    odb and does all the communication with the Abaqus API.
    This object is a like a proxy that uses the odbAccessServer to read data
    from the odb and then further processes it.

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

        def getAbqModel(recognizedBlocks):
            if self.odbIsOpen:
                return self.getAbqModel(recognizedBlocks)
            else:
                msg("OdbReaderInterpolateToPoints.initOutputPoints needs the"
                    " mesh from the odb which is not opened, yet. Opening it"
                    " now.", debugLevel=5)
                self.openOdb()
                model = self.getAbqModel(recognizedBlocks)
                self.closeOdb()
                return model
            
        # provide interpolation data to the remap kernel
        self.interpolation = InterpolMeshToPoints(
            grid=grid,
            getMeshCallBack=getAbqModel,
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
        self.comm.storeRemapWeights(
            len(self.interpolation), nodalweights, gaussPtsweights,
            defaultNotDefinedKey)
        msg("OdbReaderInterpolateToPoints.initOutputPoints:"
            " stored weights for %d points." % self.comm.nbPoints,
            debugLevel=10)

        # create filtersets in the odbAccess-server
        filterElset = set(e for e, xyz in self.interpolation.elemCoords
                          if e is not None)
        if self.odbIsOpen:
            msg("OdbReaderInterpolateToPoints.initOutputPoints:"
                " createFilterElset will be supplied with (first 50): %s"
                % sorted(filterElset)[:50], debugLevel=20)
            self.createFilterElset(filterElset)
        else:
            self.deferredElset = filterElset
            
        filterNset = self.interpolation.mesh.getConnectedNodes(filterElset)
        if self.odbIsOpen:
            msg("OdbReaderInterpolateToPoints.initOutputPoints:"
                " createFilterNset will be supplied with (first 50): %s"
                % sorted(filterNset)[:50], debugLevel=20)
            self.createFilterNset(filterNset)
        else:
            self.deferredNset = filterNset

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
            model = self.getAbqModel(recognizedBlocks=recognizedBlocks)

        # provide interpolation data to the remap kernel
        self.interpolation = InterpolMeshToCentroids(
            model=model, elems=elems, boundingBox=boundingBox)

        # store a reference of interpolation.elemList in self
        self.elemList = self.interpolation.elemList

        # store remap weights
        (nodalweights, gaussPtsweights) = self.interpolation.writeRemapParam(
            notDefinedKey=defaultNotDefinedKey)
        self.comm.storeRemapWeights(
            len(self.interpolation), nodalweights, gaussPtsweights,
            defaultNotDefinedKey)
        msg("OdbReaderInterpolateToPoints.initOutputAtCentroids:"
            " stored weights for %d points." % self.comm.nbPoints,
            debugLevel=10)

        # create filtersets in the odbAccess-server
        filterElset = sorted(self.interpolation.elemList)
        msg("OdbReaderInterpolateToPoints.initOutputAtCentroids:"
            " createFilterElset will be supplied with (first 50): %s"
            % filterElset[:50], debugLevel=20)
        self.createFilterElset(filterElset)
        filterNset = sorted(
            self.interpolation.mesh.getConnectedNodes(filterElset))
        msg("OdbReaderInterpolateToPoints.initOutputAtCentroids:"
            " createFilterNset will be supplied with (first 50): %s"
            % filterNset[:50], debugLevel=20)
        self.createFilterNset(filterNset)

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
        self.comm.storeRemapWeights(
            len(self.interpolation), nodalweights, gaussPtsweights,
            defaultNotDefinedKey)
        msg("OdbReaderInterpolateToPoints.initOutputAtFaceCentroids:"
            " stored weights for %d points." % self.comm.nbPoints,
            debugLevel=10)

        # create filtersets in the odbAccess-server
        filterElset = sorted(e for e, xyz in self.interpolation.elemCoords)
        msg("OdbReaderInterpolateToPoints.initOutputAtFaceCentroids:"
            " createFilterElset will be supplied with (first 50): %s"
            % filterElset[:50], debugLevel=20)
        self.createFilterElset(filterElset)
        filterNset = sorted(
            self.interpolation.mesh.getConnectedNodes(filterElset))
        msg("OdbReaderInterpolateToPoints.initOutputAtFaceCentroids:"
            " createFilterNset will be supplied with (first 50): %s"
            % filterNset[:50], debugLevel=20)
        self.createFilterNset(filterNset)

    #{ open/close
    def openOdb(self):
        """Open the odb.
        """
        OdbReader.openOdb(self)

        try:
            filterElset = self.deferredElset
        except AttributeError:
            pass
        else:
            msg("OdbReaderInterpolateToPoints.initOutputPoints:"
                " createFilterElset will be supplied with (first 50): %s"
                % sorted(filterElset)[:50], debugLevel=20)
            self.createFilterElset(filterElset)

        try:
            filterNset = self.deferredNset
        except AttributeError:
            pass
        else:
            msg("OdbReaderInterpolateToPoints.initOutputPoints:"
                " createFilterNset will be supplied with (first 50): %s"
                % sorted(filterNset)[:50], debugLevel=20)
            self.createFilterNset(filterNset)
    #} end of open/close
        
    def getFieldFromOdb(self, stepNb, frameNb, fieldName,
                        notDefinedValue=None, interpolationType="default"):
        """
        @param stepNb: Integer specifying the odb step. 1= first step in the odb
        @param frameNb: Integer specifying the odb frame number / index.

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

        if self.comm.nbPoints<=0:
            raise RuntimeError(
                "OdbReaderInterpolateToPoints.getFieldFromOdb: no output"
                " points defined. Maybe initOutputPoints has not been called?")

        # set frame, can raise ValueError
        self._callmethod("setFrame", stepNb, frameNb)

        # field, invariant and component will be derived from the name
        (fieldName, position, dataType, odbFieldName, component, invariant)\
            = self.getFieldProps(fieldName)

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
        self._startRequest("getFieldData", odbFieldName, component, invariant)
        self.comm.updateStructPtFieldValues(
            position, field, notDefinedArr, defaultNotDefinedKey,
            self._interpolationTypes[interpolationType])
        self._finalizeRequest()
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
                matNameToMatNumber.iteritems(), key=lambda x,y: y)
            msg("WARNING"
                " from OdbReaderInterpolateToPoints.getMatFieldFromOdb():"
                " firstMatCode %d <= last specified mat number from"
                " matNameToMatNumber %d (%s). Automatically added mat numbers"
                " starting at %d may conflict with the given"
                " matNameToMatNumber dict without notice."
                % (firstMatCode, lastNumber, lastName))
            del lastName, lastNumber

        # check that given mat names are upper case
        if matNameToMatNumber:
            for matName, matNb in matNameToMatNumber.items():
                nameUpper = matName.upper()
                if nameUpper != matName:
                    matNameToMatNumber[nameUpper] = matNameToMatNumber[matName]
                    del matNameToMatNumber[matName]

        # get the actual data from the odb
        msg("OdbReaderInterpolateToPoints.getMatFieldFromOdb: starting"
            " getMatNames on the odbAccessServer.", debugLevel=10)
        self._startRequest("getMatNames", firstMatCode, 0)
        matNumbersPrelim = list()
        self.comm.updateStructPtFieldValuesConstInt(
            matNumbersPrelim, noMatCode, defaultNotDefinedKey)
        msg("OdbReaderInterpolateToPoints.getMatFieldFromOdb: Got interpolated"
            " matNumbersPrelim for %d pts: %s"
            % (len(matNumbersPrelim), matNumbersPrelim[:10]), debugLevel=10)

        # get matNamesList from the pipe
        (nbMatNames,) = self._pipeRead("i")
        msg("OdbReaderInterpolateToPoints.getMatFieldFromOdb: will read"
            " %d matNames from the pipe." % nbMatNames, debugLevel=10)
        matNamesList = list(self._pipeRead("s"*nbMatNames))
        self._finalizeRequest()
        msg("OdbReaderInterpolateToPoints.getMatFieldFromOdb: read %d mat"
            " names." % len(matNamesList), debugLevel=10)

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
