"""
make sure you are using bae.abq_model_02.container.py with rev.6168 or later
"""

from bae.abq_model_02 import Model

for fr in [
        # "GNKT2020066_G01HTv1_SETS_Q01H_cELSETS.inp",
        "GNKT2020066_G01HTv1_D01_Q01_MOOSE_2_reaaranged.inp",
        ]:
    mdl = Model()
    mdl.read(fr)
    mdl.elset.writeCompressed = True                # default was True already
    mdl.elset.writeCompressedWithCommaOne = True    # default changed from False to True
    mdl.nset.writeCompressed = True                
    mdl.nset.writeCompressedWithCommaOne = True    
    fw = fr.replace(".inp","_GCO.inp")  # GCO means: GenerateCommaOne
    mdl.write(fw,
        header="** added ',1' to all GENERATE datalines in sets \n** from original file '%s'"%(fr),
        withSummary=True
        )
    del mdl
