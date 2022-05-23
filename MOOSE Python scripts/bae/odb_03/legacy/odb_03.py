r"""Module to handle odbs

This is the multiprocessing remote host version still under construction!

first of all: read data easily

All scripts using this module should be run with "python myscript.py" in order
to be able to use numpy.
"""

__version__ = "1.01"
_version_history_ = r"""\
Versions:

1.01 new     - based on odb_02 version 1.08
"""

import sys, os, time
from subprocess import Popen
from struct import Struct
from array import array
import cPickle as pickle

from bae.field_01 import Field
from bae.mesh_01 import Mesh, MeshStructuredPoints
from bae.abq_model_02 import Model
from bae.future_01 import defaultdict

from bae.log_01 import log, msg

#---------------------------------------------------------------------
#--- read data from the odb
#... wrapper class for reading from the odb 

class OdbReader(object):
    r"""Class to read data from the odb.
    """

    def __init__(self, odb, abaqusexe="abaqus", proxyTermTimeout=5):
        """Constructor

        @param odb: file name of the odb.
        @param abaqusexe: abaqus executable that might differentiate different
           versions. E.g. "abq6102". Might also contain the complete path.
        @param proxyTermTimeout: time period in sec to wait for the proxy to
           terminate orderly. After that period it is being killed by SIGKILL.
        """

        self.proxyTermTimeout = proxyTermTimeout

        # find the path to the proxy script to be run with abaqus python
        proxyScript = os.path.join(
            os.path.dirname(                 # the directory where we found ...
                sys.modules[self.__module__] # the current module
                .__file__ ),                 # the path to the current module
            "odb_proxy_03.py")

        # launch the proxy script
        (self._fd_inp, remote_out) = os.pipe()
        (remote_inp, self._fd_out) = os.pipe()
        self.proxy = Popen([abaqusexe, "python", proxyScript,
                            str(remote_inp), str(remote_out)])

        # translation structure for message length
        self.lenStruct = Struct("i")

        # open odb
        assert isinstance(odb, basestring)
        self.simpleCmd("openOdb", odb=odb)

    def _waitForProxyTerm(self):
        rc = None
        endtime = time.clock() + self.proxyTermTimeout
        while time.clock() < endtime:
            rc = self.proxy.poll()
            if rc is not None:
                break
            time.sleep(0.05)
        if rc is None:
            msg("Proxy still not terminated after waiting for %g seconds."
                % self.proxyTermTimeout)
        return rc

    def __del__(self):
        msg("Sending terminate to proxy.", debugLevel=10)
        self.sendCmd("terminate")
        
        rc = self._waitForProxyTerm()
        if rc is None:
            msg("Calling terminate method of proxy process.")
            self.proxy.terminate()
            rc = self._waitForProxyTerm()
        if rc is None:
            msg("Calling kill method of proxy process.")
            self.proxy.kill()
            rc = self._waitForProxyTerm()
        self.proxy.wait()

    def sendStr(self, string):
        """send string to proxy script."""
        os.write(self._fd_out, self.lenStruct.pack(len(string)))
        os.write(self._fd_out, string)

    def receiveStr(self):
        """receive message from proxy script."""
        length = self.lenStruct.unpack(
            os.read(self._fd_inp, self.lenStruct.size))[0]
        string = ""
        while length>0:
            addstr = os.read(self._fd_inp, length)
            length -= len(addstr)
            string += addstr
        return string

    def sendCmd(self, cmd, **kwargs):
        """Send the command cmd with keyword arguments to the odb reader
        proxy script. The handling of any data comming back is not done by this
        method. See self.simpleCmd() for an alternative that also deals with
        simple return values.

        Usage:
         >>> self.sendCmd("openOdb", odb="myodb.odb")
         >>> self.receiveStr() # fetch (and ignore) the return value

        @Note: The odb-proxy always sends back at least the return value of the
        function. You have to fetch this with the self.recieveStr() method
        """
        self.sendStr(cmd)
        msg("sent cmd %s" % cmd, debugLevel=10)
        self.sendStr(pickle.dumps(kwargs, 2))
        msg("sent kwargs %s" % kwargs, debugLevel=10)
        return

    def simpleCmd(self, cmd, **kwargs):
        """Send the command cmd with keyword arguments to the odb reader
        proxy script.

        Usage:
         >>> result = self.simpleCmd("openOdb", odb="myodb.odb")

        @Returns: The return value of the corresponding proxy method.
        """
        self.sendStr(cmd)
        msg("sent cmd %s" % cmd, debugLevel=10)
        self.sendStr(pickle.dumps(kwargs, 2))
        msg("sent kwargs %s" % kwargs, debugLevel=10)
        result = pickle.loads(self.receiveStr())
        msg("received return value of type %s" % type(result), debugLevel=10)
        return result

    #-------------------------------------------------------------


    # def getMesh(self):
    #     r"""get a mesh_01.Mesh-object from the odb model data,

    #     reads nodes and elements

    #     assumes there is only one part instance in the model
    #     """
    #     self.sendCmd("getMesh")
    #     msg("Requested mesh info from the odb")

    #     mesh = Mesh()
    #     nodeStruct = Struct("Lddd")
    #     stride = nodeStruct.size
    #     nodesStr = self.receiveStr()
    #     msg("Received data for %d nodes, transfering them to the mesh object."
    #         % (len(nodesStr)/stride))
    #     for i in range(0, len(nodesStr), stride):
    #         nodeLabel, x,y,z = nodeStruct.unpack(nodesStr[i:i+stride])
    #         mesh.nodeCoords[nodeLabel] = [x,y,z]
            
    #     msg("Now waiting for connectivity data from the odb.")
    #     elTypesList = pickle.loads(self.receiveStr())
    #     msg("Recieved %d element types. Reading connectivity data.")
    #     elStruct = Struct("3L")
    #     stride1 = elStruct.size
    #     stride2 = Struct("L").size
    #     connStr = self.receiveStr()
        
    #     msg("Received connectivity data for elements, transfering them to"
    #         " the mesh object.")
    #     i = 0
    #     while i<len(connStr):
    #         elemLabel, elType, nbNodes = elStruct.unpack(connStr[i:i+stride1])
    #         i += stride1
    #         stride3 = nbNodes*stride2
    #         elNodes = list(array("L", connStr[i:i+stride3]))
    #         i += stride3
    #         mesh.elNodes[elemLabel] = elNodes
    #         mesh.elShape    Model.elTypeToShapeSep

    def getAbqModel(self):
        r"""get a mesh_01.Mesh-object from the odb model data,

        reads nodes and elements

        assumes there is only one part instance in the model

        ideas for future enhancements:
         - add the possibility to selectively read only elements, nodes (see
           recognizedBlocks argument of abq_model_02.Model.read() / write())
         - read elsets, nsets, surfaces, all the rest
        """

        self.sendCmd("getMesh")
        msg("Requested mesh info from the odb")

        model = Model()
        nodeStruct = Struct("Lddd")
        stride = nodeStruct.size
        nodesStr = self.receiveStr()
        msg("Received data for %d nodes, transfering them to the model object."
            % (len(nodesStr)/stride))
        for i in range(0, len(nodesStr), stride):
            nodeLabel, x,y,z = nodeStruct.unpack(nodesStr[i:i+stride])
            model.nodeCoords[nodeLabel] = [x,y,z]
            
        msg("Now waiting for connectivity data from the odb.")
        elTypesList = pickle.loads(self.receiveStr())
        msg("Recieved %d element types. Reading connectivity data."
            % len(elTypesList))
        elStruct = Struct("3L")
        stride1 = elStruct.size
        stride2 = Struct("L").size
        connStr = self.receiveStr()
        
        msg("Received connectivity data for elements, transfering them to"
            " the model object.")
        typElNodes = defaultdict(dict)
        i = 0
        while i<len(connStr):
            elemLabel, elType, nbNodes = elStruct.unpack(connStr[i:i+stride1])
            i += stride1
            stride3 = nbNodes*stride2
            elNodes = list(array("L", connStr[i:i+stride3]))
            i += stride3
            
            typElNodes[elType][elemLabel] = nbNodes

        for elType, elNodes in typElNodes.iteritems():
            model.updateElems(elNodes, elType)
        msg("Finished transfering %d elements to the new model."
            % len(model.elNodes))

        # fetch return value of the proxy getMesh method
        self.receiveStr()
        return model

    def odbFrame(self, stepNo, frameNo):
        """Returns an OdbFrame instance identifying the specified frame in the odb.
        This can be used as argument to getOdbField for instance.
        """
        stepName = "Step-%s" % stepNo
        objId = self.simpleCmd("getFrame", stepName=stepName, frameNo=frameNo)
        return OdbFrame(self, objId)

    def odbFrameIterator(self, framesList=None):
        r"""Iterator over frames identified by step nb and frame nb, yields a
        (step number, frame number, OdbFrame object) - tuple in each iteration.

        Usage:
         >>> framesList = [
         ...    (2, [10, 11, 12]),
         ...    (3, [1, 2, 3, 4, 5]),
         ...    ]
         >>> for stepNb, frameNb, frame in self.getFrameIterator2(framesList):
         ...     print ("odb frame %d, simulation time %g"
         ...            % (frame.incrementNumber, frame.frameValue))

        @param framesList: list (or other iterable) of (step number, frame
          number list) - tuples. If None or not specified: yield all frames up
          to the end. 

          The step number is the number in the step name of the step in the
          odb. I.e. step number 3 identifies odb step "Step-3". In the ordinary
          two steps analysis the possible values are 1 and 2. Note that if the
          first step in the odb is called step-4 the corresponding step number
          is 4.

          The frame number list contains frame numbers as in the odb. Instead
          of a frame number list there can be None which means: yield all
          frames for this step. Negative values can be used to count from the
          end like usual list indices. Frame numbers that do not exist in the
          odb are ignored without any warning.

        @Returns: a (stepNo, frameNb, odbFrame) tuple on each iteration.
        """

        if framesList is None:
            framesList = [[int(s[5:]), range(f)]
                          for s,f in self.simpleCmd("getStepsList")]

        for stepNo, frameIdList in framesList:
            stepName = "Step-%s" % stepNo
            self.simpleCmd("initFrameIter",
                           stepName=stepName, frameIdList=frameIdList)
            while 1:
                frameNo, odbFrame = self.simpleCmd("frameIterNext")
                if frameNo is None:
                    break
                # important: no local reference to OdbFrame object
                # otherwise there might be communication problems when it's
                # being deleted while the next frameIterNext cmd is running
                yield (stepNo, frameNo, OdbFrame(self, odbFrame))


    def getOdbField(self, odbFrame, fieldName):
        """
        """
        pass

    #-------------------------------------------------
    #{ methods for reading and interpolation to grid

    def initOutputGrid(self, **kwargs):
        """
        You may specify either of the following argument combinations:
         - firstPoint, lastPoint and spacing
         - gridPtNbs, origin and spacing

        Specify a boxName argument for speed up on subsequent runs.

        Creates the following instance variables:
         - grid: A mesh_01.MeshStructuredPoints object representing the output
              grid
         - mesh: mesh object of only those elements that include one of the
              grid points
         - elemCoords: list of (element number, elem coord)-tuples for each
              grid point
         - filterElset: odbset containing all elements of which output data
              will be read from the odb.
         - filterNset: odbset containing all nodes of which output data
              will be read from the odb.

        @kwarg gridPtNbs: A triple of integers specifying the number of grid
          points in x, y, z direction
        @kwarg origin: first point (min x, y, z) (synonym to firstPoint)
        @kwarg spacing: A triple of doubles: Point spacing in x, y, z direction
          or just one number the same spacing in each direction.
        @kwarg firstPoint: first point (min x, y, z) (synonym to origin)
        @kwarg lastPoint: last point (max x, y, z)

        @kwarg boxName: A file "%s_coordinates.pickle"%boxName containing the
          position of the grid points in the mesh (i.e. element numbers and
          element coordinates) will be saved and used on subsequent
          invocations. If the mesh or the requested points have changed remove
          this file. If you don't specify this argument the point position will
          be calculated and no pickle file will be generated for subsequent
          invocations.

        @returns: self.grid
        """

        def getOrigin(kwargs):
            try:
                return kwargs["origin"]
            except KeyError:
                pass

            try:
                return kwargs["firstPoint"]
            except KeyError:
                raise KeyError("Neither origin nor firstPoint in keyw args.")

        def getGridPtNb(kwargs, spacing):

            # try "lastPoint"
            try:
                return [int((x1-x0)/dx)+1 for x0, x1, dx
                            in zip(origin, kwargs["lastPoint"], spacing)]
            except KeyError:
                pass

            # try "gridPtNbs"
            try: 
                return kwargs["gridPtNbs"]
            except KeyError:
                raise KeyError("Neither lastPoint nor gridPtNbs in keyw args.")

        try:
            # "origin" or "firstPoint" must be there in any case
            origin = getOrigin(kwargs)

            # "spacing" must be there
            spacing = kwargs["spacing"]
            if isinstance(spacing, (int, float)):
                spacing = [spacing,spacing,spacing]

            # "lastPoint" or "gridPtNbs" must be there
            gridPtNb = getGridPtNb(kwargs, spacing)

        except KeyError:
            raise ValueError(
                "Unsufficient arguments to define the output grid.")

        # create grid object
        msg("Creating grid points list.")
        self.grid = MeshStructuredPoints(
            gridPtNb=gridPtNb, origin=origin, spacing=spacing)
        gridPoints = list(self.grid.getPointsIter())
        msg("... finished grid with %d points" % len(gridPoints))
        if not len(gridPoints):
            raise Exception("No points in the grid.")

        # calculate or look up position of the grid points in the mesh
        boxName = kwargs.get("boxName", None)
        if boxName:
            meshinfopickleName = "%s_coordinates.pickle" % boxName
            try:
                meshPickleFile = open(meshinfopickleName, "rb")
            except IOError:
                # get meshinfo from the odb
                msg("Did not find mesh information pickle file %s."
                    % meshinfopickleName)
                getMeshinfoFromPickle = False
                writeMeshinfoPickle = True
            else:
                getMeshinfoFromPickle = True
                writeMeshinfoPickle = False

        else:
            getMeshinfoFromPickle = False
            writeMeshinfoPickle = False

        # actually read the meshinfo from the pickle file
        if getMeshinfoFromPickle:
            msg("Reading mesh information pickle file %s."%meshinfopickleName)
            msg("Loading point positions for interpolation.")
            self.mesh = pickle.load(meshPickleFile)
            msg("Loaded mesh with %d elements and %d nodes."
                % (len(self.mesh.elNodes), len(self.mesh.nodeCoords)))
            self.elemCoords = pickle.load(meshPickleFile)
            msg("Loaded point positions for %d points." % len(self.elemCoords))
            del meshPickleFile

        # actually calculate position of the grid points in the mesh
        else:
            msg("Reading mesh from the odb.")
            model = self.getAbqModel(recognizedBlocks = "NODE,ELEMENT")
            msg("Finished reading mesh information from the odb.")

            # initializing point search only within the bounding box of points
            model.initializePointSearch(box=self.grid.getBoundingBox())

            # list of (element number, elem coord)-tuples for each pt in gridPoints
            self.elemCoords = model.getElemCoords(gridPoints)
            msg("Calculated element coordinates for %d points, %d points are"
                " within the mesh." % (len(gridPoints), len([
                            None for elem, coords in self.elemCoords
                            if elem is not None])))
            allElems = set(elem for elem, coords in self.elemCoords)
            self.mesh = model.partMeshFromElements(allElems)
            msg("Created a submesh with %d elements and %d nodes."
                % (len(self.mesh.elNodes), len(self.mesh.nodeCoords)))
            del model

        assert len(gridPoints)==len(self.elemCoords)

        # write position of the grid points in the mesh to pickle file
        if writeMeshinfoPickle:
            msg("Writing mesh information pickle file %s."
                % meshinfopickleName)
            meshPickleFile = open(meshinfopickleName, "wb")
            pickle.dump(self.mesh, meshPickleFile)
            msg("Saved submesh.")
            pickle.dump(self.elemCoords, meshPickleFile)
            msg("Saved point positions for %d points." % len(self.elemCoords))
            del meshPickleFile

        # create a filterset for the output data
        msg("Creating filterset for the output data.")
        self.filterElset = self.odbSetFromElements(self.mesh.elNodes.iterkeys())
        self.filterNset = self.odbSetFromNodes(self.mesh.nodeCoords.iterkeys())
        msg("len(filterElset): %d"%len(self.filterElset.elements), debugLevel=10)
        msg("len(filterNset): %d"%len(self.filterNset.nodes), debugLevel=10)

        return self.grid

    #} end of methods for reading / interpolating


#---------------------------------------------------------------------
#--- representations of odb object (proxy objects)

class OdbObject(object):
    """A reference to an object stored in the odb proxy, possibly an Abaqus
    odb object.

    @Note: IMPORTANT! Don't delete objects of this type in the middle of other
    communication. Mind the scope of automatic variables. I.e. make sure those
    object will onlyu be deleted when no other OdbReader-command is still in
    execution.

    Objects of this class make sure that the corresponding proxy object will be
    deleted automatically when this object is being deleted. This involves some
    communication with the proxy that must not interfere with another ongoing
    interaction.
    """
    def __init__(self, odbReader, objId):
        """
        @param odbReader: ... to be added
        @param objId: ... to be added
        """
        self.odbReader = odbReader
        self.objId = objId
    
    def __del__(self):
        """Make sure that the odb object in the proxy is deleted as well
        """
        self.odbReader.simpleCmd("removeOdbObject", odbId=self.objId)

class OdbFrame(OdbObject):
    """A reference to a frame in the odb
    """
    pass
