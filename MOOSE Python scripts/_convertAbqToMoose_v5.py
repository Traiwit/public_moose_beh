"""
convertes the abaqus input files used to run an Abq/Standard (Hydro) simulation
to a Moose input file.

Using from Abq/Input:
    - mesh file containing NODE and ELEMENT (type DC3D4 and DC3D6)
    - set files containing NSET and ELSET (and SURFACE?!)
    - mpc file  containing MPC

and writing:
    - new *.inp Abaqus file consisting of
        > NODE (renumbered nodeLabels)
        > ELEMENT (renumbered elemLabels, replaced elTypes DC3D4->C3D4, DC3D6->C3D6)
        > NSET and ELSET (using renumberedLabels)
    - an *.i Moose file consisting of
        > [Mesh] ... []-block
        > [Constraints] ... []-block
        
        
        
version 5 Trai: add 0, 0 at the end of MPCs file, so the new MOOSE MPCbe can read it  



suggestiong for version 6: creating list of surfaces        
"""
from collections import defaultdict

from bae.log_01 import log as logfile, msg
from bae.abq_model_02 import Model
from bae.abq_model_02.container import MpcsList


### ###########################################################################
### #####  P A R A M E T E R S
### ###########################################################################

### ----- which abaqus-input filenames?
inpFilenamesAndBlocks = [
    ### MESH-related
    ("GNKT2020066_G01HTv1_TETL_W2FF.inp",      "NODE,ELEMENT"),
    ("GNKT2020066_G01HTv1_SETS_ALL_SETS.inp",  "NSET,ELSET"),
    ("GNKT2020066_G01HTv1_MPCS.inp",           "MPC"),
    ("GNKT2020066_G01HTv1_SETS_MPC_NSETS.inp", "NSET"),
    ### GEODOMAINS-related
   #("GNKT2020066_G01_SETS_MATVOL_D01_bELSETS.inp",    "ELSET"),
   #("GNKT2020066_G01HTv1_SETS_INTERFACE_ELSETS.inp",  "ELSET"),
    ### SEQUENCE-related
   #("GNKT2020066_G01HTv1_SETS_Q01HTshort_TOPOL_cSRFSETS.inp",     "NSET,ELSET"),
    ("GNKT2020066_G01HTv1_SETS_Q01HTshort3_DRAIN_cSRFSETS.inp",     "ELSET,SURFACE"),
   #("GNKT2020066_G01HTv1_SETS_Q01HTshort_DRAIN_mdlchgELSETS.inp", "ELSET")
    ]

### ----- which abaqus-filename (as new input for Moose)?
outFilename = "GNKT2020066_G01HTv1_D01_Q01_v7_MOOSEMESH.inp"

### ----- which moose-filename (containing [Mesh] and [Constraints]-bocks)?
mpccsvFilename = outFilename.replace(".inp","_MPCs.csv").replace("_MOOSEMESH","_MOOSEMPCS")
lkpFilenameNd = outFilename.replace("_MOOSEMESH.inp","_lkpNdLabels.csv")
lkpFilenameEl = outFilename.replace("_MOOSEMESH.inp","_lkpElLabels.csv")

### ----- which parameters for Moose's LinearNodalConstraint
# lncWeight = "1"
# lncPenalty = "1e10"


### ###########################################################################
### #####  M A I N
### ###########################################################################

msg("= some preperations...")
strHeader = [
    "**",
    "** This file was created using the following abqInput files:",
    "**",
    ]
lenFilename = 8
lenBlocks = 6
lenContent = 15
allBlocks = set()
for fn,kw in inpFilenamesAndBlocks:
    lenFilename = max(lenFilename,len(fn))
    lenBlocks = max(lenBlocks,len(kw))
    allBlocks.update(set(kw.split(",")))
allBlocks = list(allBlocks)
msg("  DEBUG: lenFilename = %d"%lenFilename)
msg("  DEBUG: lenBlocks   = %d"%lenBlocks)
msg("  DEBUG: lenContent  = %d"%lenContent)
msg("  DEBUG: allBlocks   = %s"%str(allBlocks))
strFormat = "**"+" | %%-%ds"%lenFilename+" | %%-%ds"%lenBlocks+" | %%-%ds"%lenContent
msg("  DEBUG: strFormat   = '%s'"%strFormat)
strSep = "**"+" + "+"-"*lenFilename+" + "+"-"*lenBlocks+" + "+"-"*lenContent
msg("  DEBUG: strSep      = '%s'"%strSep)
strHeader += [
    strFormat%("filename","blocks","content summary"),
    strSep,
    ]
msg("  done\n")

msg("= reading all data...")
mdl = Model()
for fn,kw in inpFilenamesAndBlocks:
    tmp = Model()
    tmp.read(fn,kw)
    strModelAdd = strFormat%(fn,kw,str(tmp).replace("Abaqus model data: ",""))
    msg("  DEBUG: %s"%strModelAdd)
    strHeader += [strModelAdd,]
    #mdl.insertModel(tmp)
    if tmp.nodeCoords:
        mdl.nodeCoords.update(tmp.nodeCoords)
        msg("  DEBUG: added nodeCoords")
    for elt in tmp.typeEl.keys():
        mdl.updateElems( dict((el,tmp.elNodes[el]) for el in tmp.typeEl[elt]), elt.replace("DC3D","C3D"))
        msg("  DEBUG: added elNodes of type '%s'->'%s'"%(elt,elt.replace("DC3D","C3D")))
    if tmp.elset:
        mdl.elset.update(tmp.elset)
        msg("  DEBUG: added elsets")
    if tmp.nset:
        mdl.nset.update(tmp.nset)
        msg("  DEBUG: added nsets")
    if tmp.surface:
        mdl.surface.update(tmp.surface)
        msg("  DEBUG: added surfaces")
    if tmp.mpc:
        mdl.mpc.extend(tmp.mpc)
        msg("  DEBUG: added mpcs")
    del tmp
strHeader += [
    strSep,
    "** | --> %s"%str(mdl),
    "** ",
    ]
msg("summary:\n"+"\n".join(strHeader))
msg("  done\n")

msg("= renumbering...")
lkpNodes,lkpElems = mdl.renumber(nodeStart=0,elemStart=1)  # default: nodeStart=1,elemStart=1
msg("  ...done with nodes,elements,nsets,elsets")
mdl.mpc = MpcsList([ (type,[lkpNodes[nd] for nd in nodes]) for type,nodes in mdl.mpc ])
msg("  ...done with mpcs")
msg("  done\n")

msg("= writing output...")
mdl.elset.writeCompressed = True                # default: True
mdl.elset.writeCompressedWithCommaOne = True    # default: False, Moose requires True
mdl.nset.writeCompressed = True                 # default: True
mdl.nset.writeCompressedWithCommaOne = True     # default: False, Moose requires True
mdl.write(
    outFilename,
    list(set(allBlocks)-set(["MPC",])),
    header="\n".join(strHeader),
    withSummary=True
    )
msg("  ... written renumbered Abq/Input to '%s'..."%outFilename)

# mooseConstraints = defaultdict(set)
# for type,nodes in mdl.mpc:
    # slaveNd = nodes[0]
    # masterNd = nodes[1]
    # mooseConstraints[masterNd].add(slaveNd)
# cnt = 0
# strConstraints = []
# for primNd,secdNds in mooseConstraints.iteritems():
    # strConstraints.extend([
        # "  [pwp_constr_%d]"%cnt,
        # "    type = LinearNodalConstraint",
        # "    variable = porepressure",
        # "    weights = %s"%lncWeight,
        # "    penalty = %s"%lncPenalty,
        # "    primary = '%d'"%primNd,
        # "    secondary_node_ids = '"+" ".join(["%d"%nd for nd in secdNds])+"'",
        # "  []",
        # ])
    # cnt += 1
fw = open(mpccsvFilename,"wb")
fw.write("\n".join([
        "slaveNd,masterNd",
        ]+[
        "%d,%d"%(nodes[0],nodes[1])
        for type,nodes in mdl.mpc
        ])+"\n"
        '0,0'"\n"
    )


fw.close()
del fw
msg("  ... written renumbered MPCs to '%s'..."%mpccsvFilename)


fw = open(lkpFilenameNd,"wb")
fw.write("\n".join([
    "mooseNdLabel,abqNdLabel",
    ]+[
    "%d,%d"%(n1,n2) for n1,n2 in sorted([
        (newNd,oldNd) for oldNd,newNd in lkpNodes.iteritems()
        ])
    ])+"\n")
fw.close()
del fw
msg("  ... written nodeLookup (Abq<->Moose) to '%s'..."%lkpFilenameNd)

fw = open(lkpFilenameEl,"wb")
fw.write("\n".join([
    "mooseElLabel,abqElLabel",
    ]+[
    "%d,%d"%(e1,e2) for e1,e2 in sorted([
        (newEl,oldEl) for oldEl,newEl in lkpElems.iteritems()
        ])
    ])+"\n")
fw.close()
del fw
msg("  ... written elemLookup (Abq<->Moose) to '%s'..."%lkpFilenameEl)
msg("  done\n")
