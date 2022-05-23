"""
convertes the abaqus input files used to assign matprops in an Abq/Standard 
(Hydro) simulation to a Moose input csv-file.

rowNumber in csv-file is the Moose-elLabel

csv-file content:
kw_xx,kw_yy,kw_zz, Aw, kw_max

todo for future, preferred:
option1: kw_xx,kw_yy,kw_zz,kw_xy,kw_xz,kw_yz, rmdmax, kw_max
option2: kw_inplane,kw_normal,n_nearing,n_plunge, rmdmax, kw_max
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
    ### GEODOMAINS-related
    ("GNKT2020066_G01_SETS_MATVOL_D01_bELSETS.inp",    "ELSET"),
    ("GNKT2020066_G01HTv1_SETS_INTERFACE_ELSETS.inp",  "ELSET"),
    ### MATERIALS-related
    ("GNKT2020066_MATPROPS_M01H.inp",      "MATERIAL"),
    ### SECTIONASSIGNMENT-related
    # ("GNKT2020066_MATSECTS_D01H.inp",      "SOLIDSECTION,COHESIVESECTION"),
    ("GNKT2020066_MATSECTS_D01H.inp",      "SOLIDSECTION,ELSET"),
    ]

### ----- which abaqus-filename (as new input for Moose)?
outFilename = "GNKT2020066_G01HTv1_D01_Q01_v5_MOOSEMESH.inp" # not written to, only to create dependent filenames

### ----- which moose-filename (containing [Mesh] and [Constraints]-bocks)?
lkpFilenameEl = outFilename.replace("_MOOSEMESH.inp","_lkpElLabels.csv")
matFilenamePerm = outFilename.replace("_MOOSEMESH.inp","_MATPROPS_PERM.csv")
matFilenamePoro = outFilename.replace("_MOOSEMESH.inp","_MATPROPS_PORO.csv")

### ----- which indices in *MATERIAL *USERMATERIAL,TYPE=THERMAL?
idxUmatPorosity = 2
idxUmatKw0 = 4
idxUmatKwmax = 5
idxUmatAw = 6


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
    # tmp = Model()
    # tmp.read(fn,kw)
    # strModelAdd = strFormat%(fn,kw,str(tmp).replace("Abaqus model data: ",""))
    # msg("  DEBUG: %s"%strModelAdd)
    # strHeader += [strModelAdd,]
    # #mdl.insertModel(tmp)
    # if tmp.elset:
        # mdl.elset.update(tmp.elset)
        # msg("  DEBUG: added elsets")
    # if tmp.nset:
        # mdl.nset.update(tmp.nset)
        # msg("  DEBUG: added nsets")
    # if tmp.surface:
        # mdl.surface.update(tmp.surface)
        # msg("  DEBUG: added surfaces")
    # if tmp.properties:
        # mdl.properties.update(tmp.properties)
        # msg("  DEBUG: added sections")
    # if tmp.material:
        # mdl.material.update(tmp.material)
        # msg("  DEBUG: added materials")
    # del tmp
    mdl.read(fn,kw)
    strModelAdd = strFormat%(fn,kw,str(mdl).replace("Abaqus model data: ",""))
    strHeader += [strModelAdd,]
strHeader += [
    strSep,
    "** | --> %s"%str(mdl),
    "** ",
    ]
msg("summary:\n"+"\n".join(strHeader))
msg("  done\n")


msg("= creating elOldToNew-dict from file '%s'..."%lkpFilenameEl)
lkpOldToNew = dict()
fr = open(lkpFilenameEl, "rb")
for line in fr:
    ln = line.strip().split(",")
    try:
        elMoose,elAbq = map(int,ln)
    except:
        pass
    else:
        lkpOldToNew[elAbq] = elMoose
fr.close()
del fr
msg("  done\n")


msg("= creating dict {elMooseLabel: matName, ...} ...")
elMaterials = defaultdict(set)
for elsName,sprop in mdl.properties.iteritems():
    for elOld in mdl.elset[elsName]:
        elMaterials[lkpOldToNew[elOld]].add(sprop['MATERIAL'])
for elNew in elMaterials.iterkeys():
    nAssignments = len(elMaterials[elNew])
    if not nAssignments==1:
        raise Exception("unidentified material assignment for mooseElLabel %d: %s"%(elNew,str(elMaterials[elNew])))
    elMaterials[elNew] = list(elMaterials[elNew])[-1]
msg("  done\n")


msg("= creating dict {matName: (kw0_xx,kw0_yy,kw0_zz,kw0_xy,kw0_xz,kw0_yz,RMDmax,kw_max), ...} ...")
lkpMatprops = dict()
for matName in mdl.material.keys():
    porosity = mdl.material[matName]['Umat'][idxUmatPorosity]
    kw0      = mdl.material[matName]['Umat'][idxUmatKw0]
    Aw       = mdl.material[matName]['Umat'][idxUmatAw]
    kwmax    = mdl.material[matName]['Umat'][idxUmatKwmax]
    lkpMatprops[matName] = dict(
        perm=",".join([ "%.4e"%val for val in [kw0,kw0,kw0,0,0,0,Aw,kwmax]]),
        poro="%.4e"%porosity,
        )
msg("  done\n")


msg("= writing matprops (perm and poro) to csv-files...")
elLabels = sorted(list(elMaterials.keys()))
fwPerm = open(matFilenamePerm,"wb")
fwPoro = open(matFilenamePoro,"wb")
cntEl = 0
for newEl in elLabels:
    cntEl += 1
    assert(cntEl==newEl)
    strMatpropsDict = lkpMatprops[elMaterials[newEl]]
    fwPerm.write(strMatpropsDict["perm"]+"\n")
    fwPoro.write(strMatpropsDict["poro"]+"\n")
fwPerm.close()
fwPoro.close()
del fwPerm,fwPoro
msg("  done\n")
