from bae.log_01 import log as logfile, msg
from bae.abq_model_02 import Model
from bae.abq_model_02.container import MpcsList

### ###########################################################################
### #####  P A R A M E T E R S 
### ###########################################################################

mpcFilename = "GNKT2020066_G01HTv1_MPCS.inp"
outFilename = "GNKT2020066_G01HTv1_MPCS+1.inp"

    
### ###########################################################################
### #####  M A I N  
### ###########################################################################

msg("= reading all data...")
mdl = Model()
mdl.read(mpcFilename,"MPC")
msg("  done\n")

msg("= renumbering...")
newMPCs = MpcsList([ (type,[nd+1 for nd in nodes]) for type,nodes in mdl.mpc ])
msg("  ...done with mpcs")
mdl.mpc = newMPCs
msg("  done\n")

msg("= writing output to '%s'..."%outFilename)
mdl.write(
    outFilename,
    "MPC",
    withSummary=True
    )
msg("  done\n")
