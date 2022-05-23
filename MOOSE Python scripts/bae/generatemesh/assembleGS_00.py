"""Some utility functions for assemblying ground support structures.

The API is not finished and maybe not very consistent. Therefore this
preliminary version should be copied to your working directory and imported
from there, i.e.:

>>> from assembleGS_00 import assembleBolts, embedFeature

Some stuff might go to gs_mesh_02 or elsewhere or the name and the whole
purpose might change... 
"""

__version__ = "0.02"

_version_history_ = """\
Versions:
=========

0.01 : GP: developed for the CULLIGS2012 project. Derived from earlier stuff:
       Yieldlok, Antamina, Perilya
0.02 : GP incompatible interface change for embedFeature: elsetName argument
       now optional.
"""

import cPickle as pickle
import csv
from math import pi

from bae.vecmath_01 import vector, length, norm, trans_vert2dir
from bae.misc_01 import Container, RowIterFromCsvFile
from bae.generatemesh.gs_mesh_02 import adjustEmbeddingPos
from bae.abq_model_02 import Model

from bae.log_01 import log, msg, MsgTicker

def assembleBolts(
    boltPosFileName,
    layerToBoltTemplate,

    outputFile,
    boltdataFileName=None,
    boltdataPickleName=None,

    embeddingHostMesh=None,

    boltFirstNode=1,
    boltFirstElem=1,
    ):
    r"""
    Assemble bolts of (possibly) different types according to a csv file
    stating location, orientation and type identifier (Rhino layer) into one
    Abaqus input file.

    @param boltPosFileName: csv file with the position and direction of the
       bolts. Columns must be "layer", "x0","y0","z0", "dx","dy","dz".
       Might be created by rs_exportLinesToCsv.py in Rhino.

    @param layerToBoltTemplate: dict { layer name: bolt model object }
       The bolt model object will be copied, moved and rotated into place.

    @param outputFile: anything that abq_model_02.Model.write accepts as
       outputFile argument. An open file object for example.
    @param boltdataFileName: If not None then write some bolt data to the
       newly created csv file with that name.
    @param boltdataPickleName: If not None then write the boltDataFile list
       to the newly created pickle file with that name.

    @param embeddingHostMesh: bolt footnode positions will be adjusted such
       that they lie inside any element of this host mesh.

    @param boltFirstNode: first node number of the bolts assembly
    @param boltFirstElem: first element number of the bolts assembly
    """

    model = Model()  # resulting GS model
    boltDataList = list() # some extra data for postprocessing and stuff

    #--- read the bolt types (layer), positions and directions
    tab = RowIterFromCsvFile(
        boltPosFileName, columns=["layer",
                                  (["x0","y0","z0"], float, "pt0"),
                                  (["dx","dy","dz"], float, "dir"),
                                  ] )

    msg("Reading bolts table and creating bolts.")
    adjustedPtsIds = list()
    maxAdjustDist = 0.0
    
    ticker = MsgTicker("Processing bolt nb %d")
    for cntBolts, row in enumerate(tab.getContainerIter()):
        ticker.msg(cntBolts+1)

        # adjust foot node position to be within a rock tet element
        # log.setDebugLevel(maxDebugLevel=0)
        footpointList = [row.pt0]
        adjustedPtsIds.extend(
            adjustEmbeddingPos(footpointList, embeddingHostMesh) )
        footpoint = footpointList[0]
        adjustDist = length(vector(row.pt0, footpoint))
        maxAdjustDist = max(adjustDist, maxAdjustDist)
        # log.setDebugLevel(maxDebugLevel=maxDebugLevel)

        # create this bolt
        boltTemplate = layerToBoltTemplate[row.layer]
        bolt = boltTemplate.copy()

        bolt.renumber(nodeStart=boltFirstNode, elemStart=boltFirstElem)
        bolt.rotate([0,0,0], trans_vert2dir(norm(row.dir)))
        bolt.translate(footpoint)

        # insert into common model
        (nodesNew, elemsNew) = model.insertModel(bolt)

        # store some bolt data
        boltDataList.append(Container(
            footpoint = footpoint,
            elastElems = [elemsNew[elem] for elem in bolt.elemsElast],
            nodeCollar = nodesNew[bolt.nodeCollar],
            nodeBoltEnd = nodesNew[bolt.nodeBoltEnd],
            nodesElast = [nodesNew[node] for node in bolt.nodesElast],
            nodesEmbedded = [nodesNew[node] for node in bolt.nodesEmbedded],
            ))

    del ticker
    cntBolts += 1
    msg("Found %d bolts." % cntBolts)

    if len(adjustedPtsIds):
        msg("WARNING: The positions of %d bolts have been corrected in order"
            " to be located inside a host element. Maximum adjustment: %g."
            % (len(adjustedPtsIds), maxAdjustDist))

    #--- write model to Abaqus input file
    model.write(outputFile)

    #--- create bolt data csv file
    if boltdataFileName is not None:
        boltDataFile = csv.writer(open(boltdataFileName, "wb"))

        nbElastElemsPerBolt = len(boltDataList[0].elastElems)
        boltDataFile.writerow(['x','y','z', 'footnode',]
                            + ['boltElems']*nbElastElemsPerBolt
                            )
        for bolt in boltDataList:
            boltDataFile.writerow(
                bolt.footpoint
                + [ bolt.nodeCollar, ]
                +  bolt.elastElems)
        del boltDataFile
        msg("wrote bolt data to %s" % boltdataFileName)


    #--- create bolt data pickle file
    if boltdataPickleName is not None:
        pickle.dump(boltDataList, open(boltdataPickleName, "wb"))
        msg("Wrote boltDataList pickle to %s" % boltdataPickleName)

    return


def embedFeature(
    sourceFile,
    outputFile,

    elsetName=None,
    firstNode=None,
    firstElem=None,

    embeddingHostMesh=None,

    embeddedElsets=None,
    ):
    r"""
    Copy a certain feature (like the fibrecrete structure) to another Abaqus
    input file given by the argument outputFile. Element and node labels will
    be adapted.

    @param sourceFile: Abaqus input file (name) name containing the features
       to be added to the model.
    @param outputFile: anything that abq_model_02.Model.write accepts as
       outputFile argument. An open file object for example.

       Please always specify this argument as keyword argument in case more
       arguments are being added in future versions!

    @param elsetName: All new features will be collected in an elset of this
       name. (Optional)
    @param firstNode: defaults to firstElem if given or else to 1
    @param firstElem: defaults to firstNode if given or else to 1

    @param embeddingHostMesh: Node positions will be adjusted such that they
       lie inside any element of this host mesh.

    @param embeddedElsets: A set of elset names.

       If embeddedElsets and elsetName are specified then elsetName will be
       added to it.

       If embeddedElsets is specified but elsetName is not then all elsets
       that are being written to the outputFile will be added to it.
    """

    if firstNode is None and firstElem is None:
        firstNode = 1
        firstElem = 1
    elif firstNode is None:
        firstNode = firstElem
    elif firstElem is None:
        firstElem = firstNode

    if elsetName:
        extraText = " for elset %s" % elsetName
    else:
        extraText = ""
    msg("Now inserting a new feature from %s%s." % (sourceFile, extraText))

    # read new features
    extraModel = Model().read(sourceFile)
    extraModel.renumber(nodeStart=firstNode, elemStart=firstElem)

    # adjust position
    if embeddingHostMesh is not None:
        adjustedPtsIds, maxAdjustDist = adjustEmbeddingPos(
            extraModel.nodeCoords, embeddingHostMesh,
            returnValues=["adjustedPtsIds", "maxAdjustDist"])
        adjusted = len(adjustedPtsIds)
        if adjusted:
            msg("%d nodes had to be adjusted by max. %g"
                % (adjusted, maxAdjustDist))

    # resulting model
    model = Model()
    nodesOldToNew, elemsOldToNew = model.insertModel(extraModel)
    if elsetName:
        model.forceElset(elsetName).update(set(extraModel.elNodes))

    if isinstance(embeddedElsets, set):
        if elsetName:
            embeddedElsets.add(elsetName)
        else:
            embeddedElsets.update(extraModel.elset)

    #--- write model to Abaqus input file
    model.write(outputFile)
    msg("Finished inserting a new feature for elset %s starting at"
        " element number %d." % (elsetName, elemsOldToNew[firstElem]))

    return
