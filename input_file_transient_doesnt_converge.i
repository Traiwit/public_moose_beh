### mpiexec -n 4 ~/moose_beh/moose_beh-opt -i /home/moose/XXXXXX_HR02_G01HTv2_Q02HfullSeq/XXXXXX2023052_R01_HR02_Q02HfullSeq_part1_TR_PS_test.i --t --color off --n-threads=8 > XXXXXX2023052_R01_HR02_Q02HfullSeq_part1_TR_PS_SEVDA_newflux.log 2>&1
### this file 'G01HTv2_Q02HfullSeq\XXXXXX2023052_R01_HR02_Q02HfullSeq_part1.i' was created 
###     by user 'arnd' on 'NIMBUS (win32)'
###     on '2024-03-15 12:15:10'
###     via '4_createMooseInputFiles_update_readfileNEW.py'
###     using 'XXXXXX2023052_R01_HR02_Q02HfullSeq_DATALIST2.csv'
### 
### this file covers the following 150-3 mining steps:
###     # NOT: 'F001_Y2018_M12', 'F002_Y2019_M06', 'F003_Y2019_M10', 
###     # BUT:
###     [ 'F004_Y2019_M11', 'F005_Y2019_M12', 'F006_Y2020_M01',
###      ...
###      'F150_Y2032_M01']
### 
### ################################################################### #######
### Global Parameters
### ################################################################### #######

### gravity vector
gravity = '0 0 -9.81'  # in m/s^2

### bounding box and zLevels
zLevelPwpBC = 10460.000

### moose porous flow dictator (to be used in every block)
[GlobalParams]
   PorousFlowDictator = uobj_pfdictator
[]

[UserObjects]
  [uobj_pfdictator]
    type = PorousFlowDictator
    porous_flow_vars = porepressure_L1
    number_fluid_phases = 1
    number_fluid_components = 1
    execute_on = 'TIMESTEP_END'   # the default
    # execute_on = 'INITIAL NONLINEAR TIMESTEP_END'
  []
[]


### ################################################################### #######
### MESH (NODE, ELEMENTS) AND
### SETS (ELSETS, SRFSETS, NSETS) AS WELL AS 
### MPCS
### ################################################################### #######

[Mesh]
  construct_node_list_from_side_list = true
  construct_side_list_from_node_list = false
  build_all_side_lowerd_mesh = false
  allow_renumbering = false
  add_subdomain_ids = '3 4'
  ### include abaqus mesh
  [mesh_G01HTv2]
    type = FileMeshGenerator
    file = 'XXXXXX2023052_R01_HR02_G01HTv2_Q02HfullSeq_part0_F000.e'
    allow_renumbering = false
  []

  [subdomains]
    type = ParsedSubdomainMeshGenerator
    input = mesh_G01HTv2
    combinatorial_geometry = 'x >= -10000 &  y >= -10000  & z >= -10000 '
    block_id = 0
  []

  [separate]
    type = MeshRepairGenerator
    input = 'subdomains'
    separate_blocks_by_element_types = true
  []

[]


[UserObjects]
  [1]
    type = CoupledVarThresholdElementSubdomainModifier
    coupled_var = 'state_M0'
    block = '1 3'  # tet element block
    criterion_type = EQUAL
    threshold = 1 
    subdomain_id = 3
    complement_subdomain_id = 1
    moving_boundary_name = moving_boundary
    execute_on = 'TIMESTEP_BEGIN'
    apply_initial_conditions = true
  []

  # [2]
  #   type = CoupledVarThresholdElementSubdomainModifier
  #   coupled_var = 'state_M0'
  #   block = '2 4'  # wedge element block
  #   criterion_type = EQUAL
  #   threshold = 1
  #   subdomain_id = 4
  #   complement_subdomain_id = 2
  #   execute_on = 'TIMESTEP_BEGIN TIMESTEP_END'
  #   apply_initial_conditions = false
  # []
[]

[MPCbe]
  csv_file = 'XXXXXX2023052_G02HTv2_MPCbe.csv'
  header = true
  offset = -1
  variable = porepressure_L1
  # penalty = 1e10
  # penalty = 1e9
  penalty = 1e6
[]


### ################################################################### #######
### MATERIAL PROPERTIES
### ################################################################### #######
### LIST OF REQUIRED MATPROPS (ACCORDING TO [KERNELS] DEFINED IN 
### SECTION 'HYDROGEOLOGICAL PROBLEM DEFINTION' BELOW):
### 
###   - FLUID PROPERTIES (IE., FLUID DENSITY, FLUID VISCOSITY, ...)
###   - ROCK PERMEABILITY TENSOR
###   - RELATIVE PERMEABILITY
###   - SATURATION
###   - POROSITY
### 
[Materials]
  ###
  ### MATPROPS --- FLUID PROPERTIES
  ###
  [mat_groundwater]
    ### this changes all time-units from seconds to years 
    ### CAUTION: within this input-file, define all quantities related to time
    ###          accordingly, except in [FluidProperties], but everywhere else
    type = PorousFlowSingleComponentFluid
    fp = the_simple_fluid
    phase = 0
    time_unit = years
  []
  ###
  ### MATPROPS --- ROCK PERMEABILITY TENSOR
  ###
  [mat_permeability_tensor]
    ### choose one of the following options:

    ### option-BEH-isotropic:
    type = PorousFlowPermeabilityPropDepIsotropicFromCsv3
    read_prop_user_object = uobj_read_perm_fromCsv
    damage = prop_damage
    state = prop_state
    state_values = '1 2 3'
    perm_state_values = '1e-15 5e-11 1e-12'
  []
  
  ###
  ### MATPROPS --- RELATIVE PERMEABILITY
  ###
  [mat_relative_permeability]
    ### choose one of the following options:

    # ### option-FLAC3: (FLAC, see https://mooseframework.inl.gov/source/materials/PorousFlowRelativePermeabilityFLAC.html)
    type = PorousFlowRelativePermeabilityFLACBE
    phase = 0
    m = 2
    state = prop_state # Kr at excav elements = 1
    cap_relperm = true
    cap_relperm_value = 0.1
  []
  
  ###
  ### MATPROPS --- SATURATION
  ###
  [mat_saturation]
    type = PorousFlow1PhaseP
    porepressure = porepressure_L1
    # capillary_pressure = uobj_saturation_optionConst
    capillary_pressure = uobj_saturation_optionVG
  []
  
  ###
  ### MATPROPS --- POROSITY
  ###
  # [mat_porosity]
  #   type = PorousFlowPorosityConst
  #   porosity = 0.25
  #   constant_on = NONE
  # []
  [mat_porosity]
    type = PorousFlowPorosityFromCSV
    read_prop_user_object = uobj_read_poro_fromCsv
    state = prop_state
    state_values = '1 2'
    poro_state_values = '1 0.25'
  []
  
  ###
  ### MATPROPS --- FLUID TEMPERATURE (NEEDED FOR SOME FLUID PROPERTIES)
  ###
  [mat_temperature]
    type = PorousFlowTemperature
  []
  
  ###
  ### MATPROPS --- FLUID MASS FRACTION (NEEDED FOR SOME FLUID PROPERTIES)
  ###
  [mat_massfrac]
    type = PorousFlowMassFraction
  []
  
  ###
  ### MATPROPS --- DAMAGE AND STATE
  ###
  [mat_damage]
    type = GenericConstantArrayBEH
    prop_name = prop_damage
    read_prop_user_object = uobj_read_damage_fromCsv
    csv_col_index = 1
  []
  [mat_state]
    type = GenericConstantArrayBEH
    prop_name = prop_state
    read_prop_user_object = uobj_read_state_fromCsv
    csv_col_index = 1
  []
[]

### -------------------------------------------------
### QUANTITIES ON WHICH THE MATPROPS ABOVE DEPEND ON:
### -------------------------------------------------
###
### MATPROP-DEPENDENCY --- THE SIMPLE FLUID (GROUNDWATER PROPERTIES)
###
[FluidProperties]
  [the_simple_fluid]
    type = SimpleFluidProperties
  []
[]

###
### MATPROP-DEPENDENCY --- CAPILLARY PRESSURE (NEEDED FOR SATURATION)
###
[UserObjects]

  [uobj_saturation_optionVG]
    ### ### vanGenuchten (https://mooseframework.inl.gov/source/userobjects/PorousFlowCapillaryPressureVG.html)
    type = PorousFlowCapillaryPressureVG
    alpha = 1e-6
    m = 0.6
    log_extension = True
    sat_lr = 0.2
  []
[]

###
### MATPROP-DEPENDENCY --- PERMEABILITY-CSV-READER
###
[UserObjects]
  [uobj_read_perm_fromCsv]
    type = PropertyReadFile
    prop_file_name = 'XXXXXX2023052_G01HTv2_MATERIALS_M01H_D02H_PERM.csv'
    nprop = 3     # that is (kiso,RMDmax,Kmax)
    # nprop = 8   # that is (kxx,kyy,kzz,kxy,kxz,kyz,RMDmax,Kmax)
    read_type = element
    execute_on = 'TIMESTEP_BEGIN'
  []
[]

###
### MATPROP-DEPENDENCY --- POROSITY-CSV-READER
###
[UserObjects]
  [uobj_read_poro_fromCsv]
    type = PropertyReadFile
    prop_file_name = 'XXXXXX2023052_G01HTv2_MATERIALS_M01H_D02H_PORO.csv'
    nprop = 1
    read_type = element
    execute_on = 'TIMESTEP_BEGIN'
  []
[]

###
### MATPROP-DEPENDENCY --- DAMAGE AND STATE CSV-READERS
###
[UserObjects]
  [uobj_read_damage_fromCsv]
    type = PropertyReadFileTimeBased
    prop_file_name = 'LOGP_FROM_R01/EXCHG_0004_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0005_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0006_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0007_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0008_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0009_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0010_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0011_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0012_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0013_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0014_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0015_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0016_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0017_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0018_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0019_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0020_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0021_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0022_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0023_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0024_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0025_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0026_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0027_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0028_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0029_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0030_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0031_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0032_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0033_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0034_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0035_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0036_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0037_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0038_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0039_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0040_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0041_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0042_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0043_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0044_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0045_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0046_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0047_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0048_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0049_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0050_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0051_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0052_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0053_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0054_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0055_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0056_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0057_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0058_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0059_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0060_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0061_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0062_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0063_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0064_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0065_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0066_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0067_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0068_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0069_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0070_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0071_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0072_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0073_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0074_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0075_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0076_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0077_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0078_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0079_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0080_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0081_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0082_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0083_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0084_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0085_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0086_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0087_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0088_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0089_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0090_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0091_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0092_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0093_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0094_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0095_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0096_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0097_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0098_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0099_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0100_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0101_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0102_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0103_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0104_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0105_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0106_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0107_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0108_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0109_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0110_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0111_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0112_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0113_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0114_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0115_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0116_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0117_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0118_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0119_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0120_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0121_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0122_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0123_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0124_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0125_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0126_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0127_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0128_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0129_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0130_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0131_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0132_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0133_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0134_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0135_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0136_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0137_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0138_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0139_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0140_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0141_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0142_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0143_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0144_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0145_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0146_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0147_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0148_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0149_LOGP.inp
                      LOGP_FROM_R01/EXCHG_0150_LOGP.inp'
    # time_per_file = '1.000 2.000 3.000 4.000 5.000 6.000 7.000 8.000 9.000 10.000 11.000 12.000 13.000 14.000 15.000 16.000 17.000 18.000 19.000 20.000 21.000 22.000 23.000 24.000 25.000 26.000 27.000 28.000 29.000 30.000 31.000 32.000 33.000 34.000 35.000 36.000 37.000 38.000 39.000 40.000 41.000 42.000 43.000 44.000 45.000 46.000 47.000 48.000 49.000 50.000 51.000 52.000 53.000 54.000 55.000 56.000 57.000 58.000 59.000 60.000 61.000 62.000 63.000 64.000 65.000 66.000 67.000 68.000 69.000 70.000 71.000 72.000 73.000 74.000 75.000 76.000 77.000 78.000 79.000 80.000 81.000 82.000 83.000 84.000 85.000 86.000 87.000 88.000 89.000 90.000 91.000 92.000 93.000 94.000 95.000 96.000 97.000 98.000 99.000 100.000 101.000 102.000 103.000 104.000 105.000 106.000 107.000 108.000 109.000 110.000 111.000 112.000 113.000 114.000 115.000 116.000 117.000 118.000 119.000 120.000 121.000 122.000 123.000 124.000 125.000 126.000 127.000 128.000 129.000 130.000 131.000 132.000 133.000 134.000 135.000 136.000 137.000 138.000 139.000 140.000 141.000 142.000 143.000 144.000 145.000 146.000 147.000 148.000 149.000 150.000'
    time_per_file = '1.833	2	2.083	2.167	2.25	2.333	2.417	2.5	2.583	2.667	2.75	2.833	2.917	3	3.083	3.167	3.25	3.333	3.417	3.5	3.583	3.667	3.75	3.833	3.917	4	4.083	4.167	4.25	4.333	4.417	4.5	4.583	4.667	4.75	4.833	4.917	5	5.083	5.167	5.25	5.333	5.417	5.5	5.583	5.667	5.75	5.833	5.917	6	6.083	6.167	6.25	6.333	6.417	6.5	6.583	6.667	6.75	6.833	6.917	7	7.083	7.167	7.25	7.333	7.417	7.5	7.583	7.667	7.75	7.833	7.917	8	8.083	8.167	8.25	8.333	8.417	8.5	8.583	8.667	8.75	8.833	8.917	9	9.083	9.167	9.25	9.333	9.417	9.5	9.583	9.667	9.75	9.833	9.917	10	10.083	10.167	10.25	10.333	10.417	10.5	10.583	10.667	10.75	10.833	10.917	11	11.083	11.167	11.25	11.333	11.417	11.5	11.583	11.667	11.75	11.833	11.917	12	12.083	12.167	12.25	12.333	12.417	12.5	12.583	12.667	12.75	12.833	12.917	13	13.083	13.167	13.25	13.333	13.417	13.5	13.583	13.667	13.75	13.833	13.917	14	14.083    '
    nprop = 2
    read_type = element
    execute_on = 'TIMESTEP_BEGIN'
  []
  [uobj_read_state_fromCsv]
    type = PropertyReadFileTimeBased
    prop_file_name = 'STATE_FROM_Q02H/EXCHG_0004_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0005_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0006_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0007_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0008_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0009_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0010_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0011_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0012_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0013_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0014_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0015_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0016_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0017_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0018_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0019_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0020_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0021_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0022_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0023_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0024_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0025_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0026_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0027_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0028_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0029_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0030_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0031_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0032_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0033_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0034_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0035_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0036_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0037_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0038_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0039_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0040_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0041_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0042_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0043_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0044_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0045_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0046_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0047_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0048_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0049_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0050_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0051_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0052_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0053_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0054_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0055_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0056_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0057_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0058_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0059_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0060_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0061_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0062_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0063_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0064_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0065_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0066_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0067_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0068_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0069_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0070_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0071_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0072_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0073_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0074_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0075_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0076_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0077_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0078_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0079_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0080_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0081_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0082_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0083_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0084_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0085_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0086_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0087_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0088_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0089_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0090_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0091_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0092_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0093_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0094_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0095_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0096_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0097_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0098_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0099_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0100_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0101_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0102_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0103_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0104_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0105_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0106_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0107_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0108_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0109_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0110_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0111_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0112_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0113_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0114_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0115_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0116_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0117_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0118_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0119_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0120_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0121_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0122_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0123_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0124_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0125_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0126_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0127_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0128_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0129_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0130_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0131_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0132_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0133_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0134_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0135_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0136_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0137_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0138_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0139_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0140_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0141_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0142_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0143_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0144_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0145_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0146_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0147_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0148_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0149_STATE.inp
                      STATE_FROM_Q02H/EXCHG_0150_STATE.inp'
    time_per_file = '1.833	2	2.083	2.167	2.25	2.333	2.417	2.5	2.583	2.667	2.75	2.833	2.917	3	3.083	3.167	3.25	3.333	3.417	3.5	3.583	3.667	3.75	3.833	3.917	4	4.083	4.167	4.25	4.333	4.417	4.5	4.583	4.667	4.75	4.833	4.917	5	5.083	5.167	5.25	5.333	5.417	5.5	5.583	5.667	5.75	5.833	5.917	6	6.083	6.167	6.25	6.333	6.417	6.5	6.583	6.667	6.75	6.833	6.917	7	7.083	7.167	7.25	7.333	7.417	7.5	7.583	7.667	7.75	7.833	7.917	8	8.083	8.167	8.25	8.333	8.417	8.5	8.583	8.667	8.75	8.833	8.917	9	9.083	9.167	9.25	9.333	9.417	9.5	9.583	9.667	9.75	9.833	9.917	10	10.083	10.167	10.25	10.333	10.417	10.5	10.583	10.667	10.75	10.833	10.917	11	11.083	11.167	11.25	11.333	11.417	11.5	11.583	11.667	11.75	11.833	11.917	12	12.083	12.167	12.25	12.333	12.417	12.5	12.583	12.667	12.75	12.833	12.917	13	13.083	13.167	13.25	13.333	13.417	13.5	13.583	13.667	13.75	13.833	13.917	14	14.083    '
    nprop = 2
    read_type = element
    execute_on = 'TIMESTEP_BEGIN'
  []
[]


### ################################################################### #######
### HYDROGEOLOGICAL PROBLEM DEFINTION
### ################################################################### #######

[Problem]
  type = 'FEProblem'   # (see https://mooseframework.inl.gov/source/problems/FEProblem.html)
  
  ### checking parameters:
  boundary_restricted_elem_integrity_check = true   # default: true
  boundary_restricted_node_integrity_check = true   # default: true
  check_uo_aux_state = false                        # default: false
  error_on_jacobian_nonzero_reallocation = false    # default: false
  fv_bcs_integrity_check = true                     # default: true
  kernel_coverage_check = false                      # default: true
  material_coverage_check = true                    # default: true
  material_dependency_check = false                  # default: true
  skip_nl_system_check = false                      # default: false
[]

[Kernels]
  [PorousFlowMassTimeDerivative]
    type = PorousFlowMassTimeDerivative
    ### see https://mooseframework.inl.gov/source/kernels/PorousFlowMassTimeDerivative.html
    variable = porepressure_L1
    fluid_component = 0
    multiply_by_density = true
  []
  
  [fluid_flux_fullupwind]
    type = PorousFlowAdvectiveFlux   # (see https://mooseframework.inl.gov/source/kernels/PorousFlowAdvectiveFlux.html)
    variable = porepressure_L1
    gravity = '${gravity}'
    fluid_component = 0
  []
[]

[Variables]
  [porepressure_L1]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxVariables]
  [saturation_M0]
    type = MooseVariableConstMonomial
  []
  [porosity_M0]
    type = MooseVariableConstMonomial
  []
  [charElLength_M0]
    type = MooseVariableConstMonomial
  []
  [darcy_M0_x]
    type = MooseVariableConstMonomial
  []
  [darcy_M0_y]
    type = MooseVariableConstMonomial
  []
  [darcy_M0_z]
    type = MooseVariableConstMonomial
  []
  # [darcy_M0_normal_srf_boundary]
  #   type = MooseVariableConstMonomial
  # []
  [permeability_M0_xx]
    type = MooseVariableConstMonomial
  []
  [damage_M0]
    type = MooseVariableConstMonomial
  []
  [state_M0]
    type = MooseVariableConstMonomial
  []
  [Rel_permeability_M0]
    type = MooseVariableConstMonomial
  []
[]

[AuxKernels]
  [Rel_permeability_M0]
    type = PorousFlowPropertyAux
    variable = Rel_permeability_M0
    property = relperm
    # execute_on = 'LINEAR TIMESTEP_END'      
    # execute_on = 'NONLINEAR TIMESTEP_END'   
    # execute_on = 'TIMESTEP_BEGIN TIMESTEP_END'
    execute_on = 'INITIAL TIMESTEP_BEGIN LINEAR TIMESTEP_END' 
  []
  [saturation_M0]
    type = PorousFlowPropertyAux
    variable = saturation_M0
    property = saturation
    # execute_on = 'LINEAR TIMESTEP_END'             # update saturation each linear interation (moose-default) 
    # execute_on = 'TIMESTEP_BEGIN TIMESTEP_END'     # keep saturation constant within a timestep
    # execute_on = 'INITIAL NONLINEAR TIMESTEP_END'  # update saturation each nonlinear interation
    execute_on = 'INITIAL TIMESTEP_BEGIN LINEAR TIMESTEP_END'       # update saturation each linear interation
  []
  [porosity_M0]
    type = PorousFlowPropertyAux
    variable = porosity_M0
    property = porosity
    # execute_on = 'LINEAR TIMESTEP_END'      
    # execute_on = 'NONLINEAR TIMESTEP_END'   
    # execute_on = 'TIMESTEP_BEGIN TIMESTEP_END'
    execute_on = 'INITIAL TIMESTEP_BEGIN LINEAR TIMESTEP_END'
  []
  [charElLength_M0]
    type = ElementLengthAuxsqrtBEH   # THIS IS OUR OWN DEVELOPMENT
    variable = charElLength_M0
    method = min
    execute_on = 'INITIAL TIMESTEP_BEGIN LINEAR'
  []
  
  [darcy_M0_x]
    # type = PorousFlowDarcyVelocityComponentYearToSec    # DEPRECATED
    type = PorousFlowDarcyVelocityComponent               # THIS IS IN M/YEAR
    variable = darcy_M0_x
    component = x
    gravity = '${gravity}'
    execute_on = 'INITIAL LINEAR TIMESTEP_END'
  []
  [darcy_M0_y]
    # type = PorousFlowDarcyVelocityComponentYearToSec    # DEPRECATED
    type = PorousFlowDarcyVelocityComponent               # THIS IS IN M/YEAR
    variable = darcy_M0_y
    component = y
    gravity = '${gravity}'
    execute_on = 'INITIAL LINEAR TIMESTEP_END'
  []
  [darcy_M0_z]
    # type = PorousFlowDarcyVelocityComponentYearToSec    # DEPRECATED
    type = PorousFlowDarcyVelocityComponent               # THIS IS IN M/YEAR
    variable = darcy_M0_z
    component = z
    gravity = '${gravity}'
    execute_on = 'INITIAL LINEAR TIMESTEP_END'
  []

  [permeability_M0_xx]
    type = PorousFlowPropertyAux
    variable = permeability_M0_xx
    property = permeability
    row = 0
    column = 0
    execute_on = 'INITIAL TIMESTEP_BEGIN LINEAR TIMESTEP_END'
  []
  [damage_M0]
    type = MaterialRealAux
    variable = damage_M0
    property = prop_damage
    # execute_on = 'LINEAR TIMESTEP_END'         # (moose-default)
    # execute_on = 'TIMESTEP_BEGIN TIMESTEP_END' # damage does not change within a timestep
    execute_on = 'INITIAL TIMESTEP_BEGIN LINEAR TIMESTEP_END'        # damage does not change within a timestep
  []
  [state_M0]
    type = MaterialRealAux
    variable = state_M0
    property = prop_state
    # execute_on = 'LINEAR TIMESTEP_END'         # (moose-default)
    # execute_on = 'TIMESTEP_BEGIN TIMESTEP_END' # state does not change within a timestep
    execute_on = 'INITIAL TIMESTEP_BEGIN LINEAR TIMESTEP_END'        # state does not change within a timestep
  []
[]


### ################################################################### #######
### CONTROLS (turn on/off objects wrt. time)
### ################################################################### #######
### CONTROLS FOR 150-3 STEPS:


[Controls]

  [controls0004]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0004_UG'
    # conditional_function = 'if(t>=1.833 & t<2.000,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>1.750 & t<=1.833,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0005]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0005_UG'
    # conditional_function = 'if(t>=2.000 & t<2.083,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>1.833 & t<=2.000,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0006]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0006_UG'
    # conditional_function = 'if(t>=2.083 & t<2.167,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>2.000 & t<=2.083,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0007]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0007_UG'
    # conditional_function = 'if(t>=2.167 & t<2.250,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>2.083 & t<=2.167,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0008]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0008_UG'
    # conditional_function = 'if(t>=2.250 & t<2.333,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>2.167 & t<=2.250,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0009]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0009_UG'
    # conditional_function = 'if(t>=2.333 & t<2.417,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>2.250 & t<=2.333,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0010]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0010_UG'
    # conditional_function = 'if(t>=2.417 & t<2.500,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>2.333 & t<=2.417,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0011]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0011_UG'
    # conditional_function = 'if(t>=2.500 & t<2.583,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>2.417 & t<=2.500,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0012]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0012_UG'
    # conditional_function = 'if(t>=2.583 & t<2.667,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>2.500 & t<=2.583,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0013]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0013_UG'
    # conditional_function = 'if(t>=2.667 & t<2.750,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>2.583 & t<=2.667,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0014]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0014_UG'
    # conditional_function = 'if(t>=2.750 & t<2.833,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>2.667 & t<=2.750,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0015]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0015_UG'
    # conditional_function = 'if(t>=2.833 & t<2.917,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>2.750 & t<=2.833,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0016]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0016_UG'
    # conditional_function = 'if(t>=2.917 & t<3.000,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>2.833 & t<=2.917,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0017]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0017_UG'
    # conditional_function = 'if(t>=3.000 & t<3.083,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>2.917 & t<=3.000,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0018]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0018_UG'
    # conditional_function = 'if(t>=3.083 & t<3.167,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>3.000 & t<=3.083,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0019]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0019_UG'
    # conditional_function = 'if(t>=3.167 & t<3.250,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>3.083 & t<=3.167,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0020]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0020_UG'
    # conditional_function = 'if(t>=3.250 & t<3.333,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>3.167 & t<=3.250,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0021]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0021_UG'
    # conditional_function = 'if(t>=3.333 & t<3.417,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>3.250 & t<=3.333,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0022]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0022_UG'
    # conditional_function = 'if(t>=3.417 & t<3.500,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>3.333 & t<=3.417,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0023]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0023_UG'
    # conditional_function = 'if(t>=3.500 & t<3.583,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>3.417 & t<=3.500,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0024]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0024_UG'
    # conditional_function = 'if(t>=3.583 & t<3.667,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>3.500 & t<=3.583,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0025]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0025_UG'
    # conditional_function = 'if(t>=3.667 & t<3.750,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>3.583 & t<=3.667,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0026]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0026_UG'
    # conditional_function = 'if(t>=3.750 & t<3.833,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>3.667 & t<=3.750,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0027]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0027_UG'
    # conditional_function = 'if(t>=3.833 & t<3.917,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>3.750 & t<=3.833,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0028]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0028_UG'
    # conditional_function = 'if(t>=3.917 & t<4.000,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>3.833 & t<=3.917,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0029]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0029_UG'
    # conditional_function = 'if(t>=4.000 & t<4.083,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>3.917 & t<=4.000,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0030]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0030_UG'
    # conditional_function = 'if(t>=4.083 & t<4.167,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>4.000 & t<=4.083,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0031]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0031_UG'
    # conditional_function = 'if(t>=4.167 & t<4.250,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>4.083 & t<=4.167,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0032]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0032_UG'
    # conditional_function = 'if(t>=4.250 & t<4.333,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>4.167 & t<=4.250,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0033]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0033_UG'
    # conditional_function = 'if(t>=4.333 & t<4.417,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>4.250 & t<=4.333,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0034]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0034_UG'
    # conditional_function = 'if(t>=4.417 & t<4.500,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>4.333 & t<=4.417,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0035]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0035_UG'
    # conditional_function = 'if(t>=4.500 & t<4.583,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>4.417 & t<=4.500,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0036]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0036_UG'
    # conditional_function = 'if(t>=4.583 & t<4.667,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>4.500 & t<=4.583,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0037]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0037_UG'
    # conditional_function = 'if(t>=4.667 & t<4.750,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>4.583 & t<=4.667,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0038]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0038_UG'
    # conditional_function = 'if(t>=4.750 & t<4.833,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>4.667 & t<=4.750,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0039]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0039_UG'
    # conditional_function = 'if(t>=4.833 & t<4.917,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>4.750 & t<=4.833,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0040]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0040_UG'
    # conditional_function = 'if(t>=4.917 & t<5.000,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>4.833 & t<=4.917,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0041]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0041_UG'
    # conditional_function = 'if(t>=5.000 & t<5.083,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>4.917 & t<=5.000,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0042]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0042_UG'
    # conditional_function = 'if(t>=5.083 & t<5.167,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>5.000 & t<=5.083,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0043]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0043_UG'
    # conditional_function = 'if(t>=5.167 & t<5.250,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>5.083 & t<=5.167,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0044]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0044_UG'
    # conditional_function = 'if(t>=5.250 & t<5.333,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>5.167 & t<=5.250,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0045]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0045_UG'
    # conditional_function = 'if(t>=5.333 & t<5.417,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>5.250 & t<=5.333,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0046]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0046_UG'
    # conditional_function = 'if(t>=5.417 & t<5.500,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>5.333 & t<=5.417,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0047]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0047_UG'
    # conditional_function = 'if(t>=5.500 & t<5.583,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>5.417 & t<=5.500,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0048]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0048_UG'
    # conditional_function = 'if(t>=5.583 & t<5.667,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>5.500 & t<=5.583,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0049]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0049_UG'
    # conditional_function = 'if(t>=5.667 & t<5.750,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>5.583 & t<=5.667,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0050]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0050_UG'
    # conditional_function = 'if(t>=5.750 & t<5.833,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>5.667 & t<=5.750,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0051]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0051_UG'
    # conditional_function = 'if(t>=5.833 & t<5.917,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>5.750 & t<=5.833,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0052]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0052_UG'
    # conditional_function = 'if(t>=5.917 & t<6.000,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>5.833 & t<=5.917,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0053]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0053_UG'
    # conditional_function = 'if(t>=6.000 & t<6.083,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>5.917 & t<=6.000,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0054]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0054_UG'
    # conditional_function = 'if(t>=6.083 & t<6.167,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>6.000 & t<=6.083,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0055]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0055_UG'
    # conditional_function = 'if(t>=6.167 & t<6.250,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>6.083 & t<=6.167,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0056]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0056_UG'
    # conditional_function = 'if(t>=6.250 & t<6.333,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>6.167 & t<=6.250,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0057]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0057_UG'
    # conditional_function = 'if(t>=6.333 & t<6.417,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>6.250 & t<=6.333,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0058]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0058_UG'
    # conditional_function = 'if(t>=6.417 & t<6.500,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>6.333 & t<=6.417,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0059]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0059_UG'
    # conditional_function = 'if(t>=6.500 & t<6.583,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>6.417 & t<=6.500,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0060]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0060_UG'
    # conditional_function = 'if(t>=6.583 & t<6.667,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>6.500 & t<=6.583,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0061]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0061_UG'
    # conditional_function = 'if(t>=6.667 & t<6.750,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>6.583 & t<=6.667,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0062]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0062_UG'
    # conditional_function = 'if(t>=6.750 & t<6.833,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>6.667 & t<=6.750,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0063]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0063_UG'
    # conditional_function = 'if(t>=6.833 & t<6.917,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>6.750 & t<=6.833,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0064]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0064_UG'
    # conditional_function = 'if(t>=6.917 & t<7.000,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>6.833 & t<=6.917,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0065]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0065_UG'
    # conditional_function = 'if(t>=7.000 & t<7.083,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>6.917 & t<=7.000,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0066]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0066_UG'
    # conditional_function = 'if(t>=7.083 & t<7.167,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>7.000 & t<=7.083,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0067]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0067_UG'
    # conditional_function = 'if(t>=7.167 & t<7.250,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>7.083 & t<=7.167,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0068]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0068_UG'
    # conditional_function = 'if(t>=7.250 & t<7.333,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>7.167 & t<=7.250,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0069]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0069_UG'
    # conditional_function = 'if(t>=7.333 & t<7.417,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>7.250 & t<=7.333,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0070]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0070_UG'
    # conditional_function = 'if(t>=7.417 & t<7.500,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>7.333 & t<=7.417,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0071]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0071_UG'
    # conditional_function = 'if(t>=7.500 & t<7.583,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>7.417 & t<=7.500,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0072]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0072_UG'
    # conditional_function = 'if(t>=7.583 & t<7.667,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>7.500 & t<=7.583,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0073]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0073_UG'
    # conditional_function = 'if(t>=7.667 & t<7.750,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>7.583 & t<=7.667,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0074]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0074_UG'
    # conditional_function = 'if(t>=7.750 & t<7.833,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>7.667 & t<=7.750,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0075]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0075_UG'
    # conditional_function = 'if(t>=7.833 & t<7.917,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>7.750 & t<=7.833,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0076]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0076_UG'
    # conditional_function = 'if(t>=7.917 & t<8.000,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>7.833 & t<=7.917,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0077]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0077_UG'
    # conditional_function = 'if(t>=8.000 & t<8.083,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>7.917 & t<=8.000,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0078]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0078_UG'
    # conditional_function = 'if(t>=8.083 & t<8.167,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>8.000 & t<=8.083,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0079]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0079_UG'
    # conditional_function = 'if(t>=8.167 & t<8.250,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>8.083 & t<=8.167,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0080]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0080_UG'
    # conditional_function = 'if(t>=8.250 & t<8.333,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>8.167 & t<=8.250,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0081]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0081_UG'
    # conditional_function = 'if(t>=8.333 & t<8.417,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>8.250 & t<=8.333,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0082]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0082_UG'
    # conditional_function = 'if(t>=8.417 & t<8.500,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>8.333 & t<=8.417,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0083]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0083_UG'
    # conditional_function = 'if(t>=8.500 & t<8.583,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>8.417 & t<=8.500,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0084]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0084_UG'
    # conditional_function = 'if(t>=8.583 & t<8.667,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>8.500 & t<=8.583,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0085]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0085_UG'
    # conditional_function = 'if(t>=8.667 & t<8.750,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>8.583 & t<=8.667,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0086]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0086_UG'
    # conditional_function = 'if(t>=8.750 & t<8.833,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>8.667 & t<=8.750,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0087]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0087_UG'
    # conditional_function = 'if(t>=8.833 & t<8.917,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>8.750 & t<=8.833,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0088]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0088_UG'
    # conditional_function = 'if(t>=8.917 & t<9.000,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>8.833 & t<=8.917,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0089]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0089_UG'
    # conditional_function = 'if(t>=9.000 & t<9.083,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>8.917 & t<=9.000,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0090]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0090_UG'
    # conditional_function = 'if(t>=9.083 & t<9.167,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>9.000 & t<=9.083,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0091]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0091_UG'
    # conditional_function = 'if(t>=9.167 & t<9.250,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>9.083 & t<=9.167,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0092]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0092_UG'
    # conditional_function = 'if(t>=9.250 & t<9.333,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>9.167 & t<=9.250,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0093]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0093_UG'
    # conditional_function = 'if(t>=9.333 & t<9.417,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>9.250 & t<=9.333,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0094]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0094_UG'
    # conditional_function = 'if(t>=9.417 & t<9.500,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>9.333 & t<=9.417,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0095]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0095_UG'
    # conditional_function = 'if(t>=9.500 & t<9.583,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>9.417 & t<=9.500,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0096]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0096_UG'
    # conditional_function = 'if(t>=9.583 & t<9.667,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>9.500 & t<=9.583,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0097]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0097_UG'
    # conditional_function = 'if(t>=9.667 & t<9.750,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>9.583 & t<=9.667,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0098]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0098_UG'
    # conditional_function = 'if(t>=9.750 & t<9.833,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>9.667 & t<=9.750,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0099]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0099_UG'
    # conditional_function = 'if(t>=9.833 & t<9.917,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>9.750 & t<=9.833,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0100]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0100_UG'
    # conditional_function = 'if(t>=9.917 & t<10.000,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>9.833 & t<=9.917,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0101]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0101_UG'
    # conditional_function = 'if(t>=10.000 & t<10.083,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>9.917 & t<=10.000,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0102]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0102_UG'
    # conditional_function = 'if(t>=10.083 & t<10.167,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>10.000 & t<=10.083,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0103]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0103_UG'
    # conditional_function = 'if(t>=10.167 & t<10.250,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>10.083 & t<=10.167,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0104]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0104_UG'
    # conditional_function = 'if(t>=10.250 & t<10.333,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>10.167 & t<=10.250,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0105]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0105_UG'
    # conditional_function = 'if(t>=10.333 & t<10.417,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>10.250 & t<=10.333,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0106]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0106_UG'
    # conditional_function = 'if(t>=10.417 & t<10.500,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>10.333 & t<=10.417,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0107]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0107_UG'
    # conditional_function = 'if(t>=10.500 & t<10.583,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>10.417 & t<=10.500,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0108]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0108_UG'
    # conditional_function = 'if(t>=10.583 & t<10.667,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>10.500 & t<=10.583,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0109]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0109_UG'
    # conditional_function = 'if(t>=10.667 & t<10.750,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>10.583 & t<=10.667,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0110]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0110_UG'
    # conditional_function = 'if(t>=10.750 & t<10.833,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>10.667 & t<=10.750,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0111]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0111_UG'
    # conditional_function = 'if(t>=10.833 & t<10.917,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>10.750 & t<=10.833,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0112]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0112_UG'
    # conditional_function = 'if(t>=10.917 & t<11.000,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>10.833 & t<=10.917,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0113]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0113_UG'
    # conditional_function = 'if(t>=11.000 & t<11.083,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>10.917 & t<=11.000,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0114]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0114_UG'
    # conditional_function = 'if(t>=11.083 & t<11.167,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>11.000 & t<=11.083,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0115]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0115_UG'
    # conditional_function = 'if(t>=11.167 & t<11.250,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>11.083 & t<=11.167,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0116]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0116_UG'
    # conditional_function = 'if(t>=11.250 & t<11.333,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>11.167 & t<=11.250,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0117]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0117_UG'
    # conditional_function = 'if(t>=11.333 & t<11.417,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>11.250 & t<=11.333,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0118]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0118_UG'
    # conditional_function = 'if(t>=11.417 & t<11.500,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>11.333 & t<=11.417,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0119]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0119_UG'
    # conditional_function = 'if(t>=11.500 & t<11.583,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>11.417 & t<=11.500,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0120]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0120_UG'
    # conditional_function = 'if(t>=11.583 & t<11.667,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>11.500 & t<=11.583,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0121]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0121_UG'
    # conditional_function = 'if(t>=11.667 & t<11.750,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>11.583 & t<=11.667,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0122]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0122_UG'
    # conditional_function = 'if(t>=11.750 & t<11.833,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>11.667 & t<=11.750,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0123]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0123_UG'
    # conditional_function = 'if(t>=11.833 & t<11.917,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>11.750 & t<=11.833,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0124]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0124_UG'
    # conditional_function = 'if(t>=11.917 & t<12.000,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>11.833 & t<=11.917,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0125]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0125_UG'
    # conditional_function = 'if(t>=12.000 & t<12.083,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>11.917 & t<=12.000,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0126]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0126_UG'
    # conditional_function = 'if(t>=12.083 & t<12.167,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>12.000 & t<=12.083,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0127]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0127_UG'
    # conditional_function = 'if(t>=12.167 & t<12.250,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>12.083 & t<=12.167,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0128]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0128_UG'
    # conditional_function = 'if(t>=12.250 & t<12.333,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>12.167 & t<=12.250,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0129]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0129_UG'
    # conditional_function = 'if(t>=12.333 & t<12.417,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>12.250 & t<=12.333,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0130]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0130_UG'
    # conditional_function = 'if(t>=12.417 & t<12.500,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>12.333 & t<=12.417,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0131]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0131_UG'
    # conditional_function = 'if(t>=12.500 & t<12.583,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>12.417 & t<=12.500,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0132]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0132_UG'
    # conditional_function = 'if(t>=12.583 & t<12.667,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>12.500 & t<=12.583,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0133]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0133_UG'
    # conditional_function = 'if(t>=12.667 & t<12.750,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>12.583 & t<=12.667,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0134]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0134_UG'
    # conditional_function = 'if(t>=12.750 & t<12.833,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>12.667 & t<=12.750,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0135]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0135_UG'
    # conditional_function = 'if(t>=12.833 & t<12.917,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>12.750 & t<=12.833,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0136]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0136_UG'
    # conditional_function = 'if(t>=12.917 & t<13.000,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>12.833 & t<=12.917,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0137]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0137_UG'
    # conditional_function = 'if(t>=13.000 & t<13.083,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>12.917 & t<=13.000,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0138]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0138_UG'
    # conditional_function = 'if(t>=13.083 & t<13.167,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>13.000 & t<=13.083,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0139]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0139_UG'
    # conditional_function = 'if(t>=13.167 & t<13.250,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>13.083 & t<=13.167,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0140]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0140_UG'
    # conditional_function = 'if(t>=13.250 & t<13.333,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>13.167 & t<=13.250,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0141]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0141_UG'
    # conditional_function = 'if(t>=13.333 & t<13.417,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>13.250 & t<=13.333,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0142]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0142_UG'
    # conditional_function = 'if(t>=13.417 & t<13.500,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>13.333 & t<=13.417,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0143]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0143_UG'
    # conditional_function = 'if(t>=13.500 & t<13.583,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>13.417 & t<=13.500,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0144]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0144_UG'
    # conditional_function = 'if(t>=13.583 & t<13.667,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>13.500 & t<=13.583,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0145]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0145_UG'
    # conditional_function = 'if(t>=13.667 & t<13.750,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>13.583 & t<=13.667,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0146]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0146_UG'
    # conditional_function = 'if(t>=13.750 & t<13.833,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>13.667 & t<=13.750,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0147]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0147_UG'
    # conditional_function = 'if(t>=13.833 & t<13.917,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>13.750 & t<=13.833,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0148]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0148_UG'
    # conditional_function = 'if(t>=13.917 & t<14.000,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>13.833 & t<=13.917,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0149]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0149_UG'
    # conditional_function = 'if(t>=14.000 & t<14.083,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>13.917 & t<=14.000,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
  [controls0150]
    type = ConditionalFunctionEnableControl
    enable_objects = '    BCs::bc_drain0150_UG'
    # conditional_function = 'if(t>=14.083 & t<14.166,1,0)'
    # execute_on = 'INITIAL TIMESTEP_BEGIN TIMESTEP_END'
    conditional_function = 'if(t>14.000 & t<=14.083,1,0)'
    execute_on = 'INITIAL TIMESTEP_BEGIN NONLINEAR'
  []
[]

### ################################################################### #######
### INITIAL CONDITIONS (ie., the results of previous frame)
### ################################################################### #######

[ICs]
  [ic_porepressure_L1]
    type = FunctionIC
    variable = porepressure_L1
    # function = '0.'
    # function = 'fcn_pwp_initialCondition'
    # function = 'fcn_pwp_initialCondition_csv'
    function = fcn_pwp_from_prev_results
  []
[]
[Functions]
  [fcn_pwp_initialCondition]
    type = ParsedFunction
    expression = '9810.*(${zLevelPwpBC}-z)'
  []
  # [fcn_pwp_initialCondition_csv]
  #   type = PiecewiseConstantFromCSV
  #   read_prop_user_object = uobj_pwp_read_AllNodes
  #   read_type = node
  #   column_number = 1
  # []
  [fcn_pwp_from_prev_results]
    type = SolutionFunction
    from_variable = porepressure_L1
    solution = uobj_prev_results
  []
[]
[UserObjects]
  # [uobj_pwp_read_AllNodes]
  #   type = PropertyReadFile
  #   prop_file_name = 'XXXXXX2023052_G01HTv2_IC_atALLNODES_rhoGdZminus10m.csv'
  #   nprop = 2
  #   read_type = node
  # []
  [uobj_prev_results]
    type = SolutionUserObject
    # mesh = 'outputExodus/XXXXXX2023052_R01_HR02_G01HTv2_Q02HfullSeq_part0_F003.e'
    mesh = 'XXXXXX2023052_R01_HR02_G01HTv2_Q02HfullSeq_part0_F004.e' # this is correct part0_0003 is F004.e
    timestep = 'LATEST'
  []
[]


### ################################################################### #######
### BOUNDARY CONDITIONS
### ################################################################### #######

## BOUNDARY CONDITIONS --- FARFIELD RECHARGE (VIA WALLS)
[BCs]
  [bc_walls]
    type = FunctionDirichletBC
    variable = porepressure_L1
    boundary = 'BOUNDSRF_WALLS'
    # function = fcn_pwp_walls
    function = fcn_pwp_initialCondition
    # function = fcn_pwp_initialCondition_csv
  []


  [bc_topol0003_LGFlux]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'TOPOLSRFNOTTSF_F003TOF150'
    # pt_vals = '-1e9 -200 0 1e9'
    # multipliers = '-200 -200 0 1e9'
    pt_vals = '-1e9  0 1e9'
    multipliers = '-1e9 0 1e9'
    flux_function = 1
    PT_shift = -539550.0   # in Pa, ie., -9810.N/m3*55.0m
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []

  [bc_topol0003_RiverLake]
    type = PorousFlowSink
    variable = porepressure_L1
    boundary = 'TOPOLSRFTSF_F003TOF150'
    # function = fcn_pwp_tailingsBC
    flux_function = -5000
  []

  [bc_topol0003_PitLake]
    type = FunctionDirichletBC
    variable = porepressure_L1
    boundary = 'DRAINSRFPITLOW_F003TOF150'
    function = fcn_pwp_pitlakeBC
  []

  [bc_drain0003_Pit]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFPITUPP_F003TOF150'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []  
[]

### BOUNDARY CONDITIONS FOR 150-3 STEPS
###   --- RECHARGE VIA TOPOLOGY (LinGradFlux and/or Rainfall and/or Rivers+Lakes)
###   --- DRAINAGE (PIT and/or UG)
[BCs]
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0004: F004_Y2019_M11
  ### NOT DEFINING [bc_topol0004_LGFlux]
  ### NOT DEFINING [bc_topol0004_Rainfall]
  ### NOT DEFINING [bc_topol0004_RiverLake]
  ### NOT DEFINING [bc_drain0004_Pit]
  [bc_drain0004_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F004_Y2019_M11'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0005: F005_Y2019_M12
  ### NOT DEFINING [bc_topol0005_LGFlux]
  ### NOT DEFINING [bc_topol0005_Rainfall]
  ### NOT DEFINING [bc_topol0005_RiverLake]
  ### NOT DEFINING [bc_drain0005_Pit]
  [bc_drain0005_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F005_Y2019_M12'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0006: F006_Y2020_M01
  ### NOT DEFINING [bc_topol0006_LGFlux]
  ### NOT DEFINING [bc_topol0006_Rainfall]
  ### NOT DEFINING [bc_topol0006_RiverLake]
  ### NOT DEFINING [bc_drain0006_Pit]
  [bc_drain0006_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F006_Y2020_M01'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0007: F007_Y2020_M02
  ### NOT DEFINING [bc_topol0007_LGFlux]
  ### NOT DEFINING [bc_topol0007_Rainfall]
  ### NOT DEFINING [bc_topol0007_RiverLake]
  ### NOT DEFINING [bc_drain0007_Pit]
  [bc_drain0007_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F007_Y2020_M02'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0008: F008_Y2020_M03
  ### NOT DEFINING [bc_topol0008_LGFlux]
  ### NOT DEFINING [bc_topol0008_Rainfall]
  ### NOT DEFINING [bc_topol0008_RiverLake]
  ### NOT DEFINING [bc_drain0008_Pit]
  [bc_drain0008_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F008_Y2020_M03'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0009: F009_Y2020_M04
  ### NOT DEFINING [bc_topol0009_LGFlux]
  ### NOT DEFINING [bc_topol0009_Rainfall]
  ### NOT DEFINING [bc_topol0009_RiverLake]
  ### NOT DEFINING [bc_drain0009_Pit]
  [bc_drain0009_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F009_Y2020_M04'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0010: F010_Y2020_M05
  ### NOT DEFINING [bc_topol0010_LGFlux]
  ### NOT DEFINING [bc_topol0010_Rainfall]
  ### NOT DEFINING [bc_topol0010_RiverLake]
  ### NOT DEFINING [bc_drain0010_Pit]
  [bc_drain0010_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F010_Y2020_M05'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0011: F011_Y2020_M06
  ### NOT DEFINING [bc_topol0011_LGFlux]
  ### NOT DEFINING [bc_topol0011_Rainfall]
  ### NOT DEFINING [bc_topol0011_RiverLake]
  ### NOT DEFINING [bc_drain0011_Pit]
  [bc_drain0011_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F011_Y2020_M06'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0012: F012_Y2020_M07
  ### NOT DEFINING [bc_topol0012_LGFlux]
  ### NOT DEFINING [bc_topol0012_Rainfall]
  ### NOT DEFINING [bc_topol0012_RiverLake]
  ### NOT DEFINING [bc_drain0012_Pit]
  [bc_drain0012_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F012_Y2020_M07'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0013: F013_Y2020_M08
  ### NOT DEFINING [bc_topol0013_LGFlux]
  ### NOT DEFINING [bc_topol0013_Rainfall]
  ### NOT DEFINING [bc_topol0013_RiverLake]
  ### NOT DEFINING [bc_drain0013_Pit]
  [bc_drain0013_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F013_Y2020_M08'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0014: F014_Y2020_M09
  ### NOT DEFINING [bc_topol0014_LGFlux]
  ### NOT DEFINING [bc_topol0014_Rainfall]
  ### NOT DEFINING [bc_topol0014_RiverLake]
  ### NOT DEFINING [bc_drain0014_Pit]
  [bc_drain0014_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F014_Y2020_M09'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0015: F015_Y2020_M10
  ### NOT DEFINING [bc_topol0015_LGFlux]
  ### NOT DEFINING [bc_topol0015_Rainfall]
  ### NOT DEFINING [bc_topol0015_RiverLake]
  ### NOT DEFINING [bc_drain0015_Pit]
  [bc_drain0015_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F015_Y2020_M10'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0016: F016_Y2020_M11
  ### NOT DEFINING [bc_topol0016_LGFlux]
  ### NOT DEFINING [bc_topol0016_Rainfall]
  ### NOT DEFINING [bc_topol0016_RiverLake]
  ### NOT DEFINING [bc_drain0016_Pit]
  [bc_drain0016_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F016_Y2020_M11'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0017: F017_Y2020_M12
  ### NOT DEFINING [bc_topol0017_LGFlux]
  ### NOT DEFINING [bc_topol0017_Rainfall]
  ### NOT DEFINING [bc_topol0017_RiverLake]
  ### NOT DEFINING [bc_drain0017_Pit]
  [bc_drain0017_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F017_Y2020_M12'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0018: F018_Y2021_M01
  ### NOT DEFINING [bc_topol0018_LGFlux]
  ### NOT DEFINING [bc_topol0018_Rainfall]
  ### NOT DEFINING [bc_topol0018_RiverLake]
  ### NOT DEFINING [bc_drain0018_Pit]
  [bc_drain0018_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F018_Y2021_M01'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0019: F019_Y2021_M02
  ### NOT DEFINING [bc_topol0019_LGFlux]
  ### NOT DEFINING [bc_topol0019_Rainfall]
  ### NOT DEFINING [bc_topol0019_RiverLake]
  ### NOT DEFINING [bc_drain0019_Pit]
  [bc_drain0019_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F019_Y2021_M02'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0020: F020_Y2021_M03
  ### NOT DEFINING [bc_topol0020_LGFlux]
  ### NOT DEFINING [bc_topol0020_Rainfall]
  ### NOT DEFINING [bc_topol0020_RiverLake]
  ### NOT DEFINING [bc_drain0020_Pit]
  [bc_drain0020_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F020_Y2021_M03'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0021: F021_Y2021_M04
  ### NOT DEFINING [bc_topol0021_LGFlux]
  ### NOT DEFINING [bc_topol0021_Rainfall]
  ### NOT DEFINING [bc_topol0021_RiverLake]
  ### NOT DEFINING [bc_drain0021_Pit]
  [bc_drain0021_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F021_Y2021_M04'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0022: F022_Y2021_M05
  ### NOT DEFINING [bc_topol0022_LGFlux]
  ### NOT DEFINING [bc_topol0022_Rainfall]
  ### NOT DEFINING [bc_topol0022_RiverLake]
  ### NOT DEFINING [bc_drain0022_Pit]
  [bc_drain0022_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F022_Y2021_M05'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0023: F023_Y2021_M06
  ### NOT DEFINING [bc_topol0023_LGFlux]
  ### NOT DEFINING [bc_topol0023_Rainfall]
  ### NOT DEFINING [bc_topol0023_RiverLake]
  ### NOT DEFINING [bc_drain0023_Pit]
  [bc_drain0023_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F023_Y2021_M06'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0024: F024_Y2021_M07
  ### NOT DEFINING [bc_topol0024_LGFlux]
  ### NOT DEFINING [bc_topol0024_Rainfall]
  ### NOT DEFINING [bc_topol0024_RiverLake]
  ### NOT DEFINING [bc_drain0024_Pit]
  [bc_drain0024_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F024_Y2021_M07'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0025: F025_Y2021_M08
  ### NOT DEFINING [bc_topol0025_LGFlux]
  ### NOT DEFINING [bc_topol0025_Rainfall]
  ### NOT DEFINING [bc_topol0025_RiverLake]
  ### NOT DEFINING [bc_drain0025_Pit]
  [bc_drain0025_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F025_Y2021_M08'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0026: F026_Y2021_M09
  ### NOT DEFINING [bc_topol0026_LGFlux]
  ### NOT DEFINING [bc_topol0026_Rainfall]
  ### NOT DEFINING [bc_topol0026_RiverLake]
  ### NOT DEFINING [bc_drain0026_Pit]
  [bc_drain0026_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F026_Y2021_M09'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0027: F027_Y2021_M10
  ### NOT DEFINING [bc_topol0027_LGFlux]
  ### NOT DEFINING [bc_topol0027_Rainfall]
  ### NOT DEFINING [bc_topol0027_RiverLake]
  ### NOT DEFINING [bc_drain0027_Pit]
  [bc_drain0027_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F027_Y2021_M10'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0028: F028_Y2021_M11
  ### NOT DEFINING [bc_topol0028_LGFlux]
  ### NOT DEFINING [bc_topol0028_Rainfall]
  ### NOT DEFINING [bc_topol0028_RiverLake]
  ### NOT DEFINING [bc_drain0028_Pit]
  [bc_drain0028_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F028_Y2021_M11'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0029: F029_Y2021_M12
  ### NOT DEFINING [bc_topol0029_LGFlux]
  ### NOT DEFINING [bc_topol0029_Rainfall]
  ### NOT DEFINING [bc_topol0029_RiverLake]
  ### NOT DEFINING [bc_drain0029_Pit]
  [bc_drain0029_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F029_Y2021_M12'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0030: F030_Y2022_M01
  ### NOT DEFINING [bc_topol0030_LGFlux]
  ### NOT DEFINING [bc_topol0030_Rainfall]
  ### NOT DEFINING [bc_topol0030_RiverLake]
  ### NOT DEFINING [bc_drain0030_Pit]
  [bc_drain0030_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F030_Y2022_M01'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0031: F031_Y2022_M02
  ### NOT DEFINING [bc_topol0031_LGFlux]
  ### NOT DEFINING [bc_topol0031_Rainfall]
  ### NOT DEFINING [bc_topol0031_RiverLake]
  ### NOT DEFINING [bc_drain0031_Pit]
  [bc_drain0031_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F031_Y2022_M02'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0032: F032_Y2022_M03
  ### NOT DEFINING [bc_topol0032_LGFlux]
  ### NOT DEFINING [bc_topol0032_Rainfall]
  ### NOT DEFINING [bc_topol0032_RiverLake]
  ### NOT DEFINING [bc_drain0032_Pit]
  [bc_drain0032_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F032_Y2022_M03'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0033: F033_Y2022_M04
  ### NOT DEFINING [bc_topol0033_LGFlux]
  ### NOT DEFINING [bc_topol0033_Rainfall]
  ### NOT DEFINING [bc_topol0033_RiverLake]
  ### NOT DEFINING [bc_drain0033_Pit]
  [bc_drain0033_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F033_Y2022_M04'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0034: F034_Y2022_M05
  ### NOT DEFINING [bc_topol0034_LGFlux]
  ### NOT DEFINING [bc_topol0034_Rainfall]
  ### NOT DEFINING [bc_topol0034_RiverLake]
  ### NOT DEFINING [bc_drain0034_Pit]
  [bc_drain0034_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F034_Y2022_M05'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0035: F035_Y2022_M06
  ### NOT DEFINING [bc_topol0035_LGFlux]
  ### NOT DEFINING [bc_topol0035_Rainfall]
  ### NOT DEFINING [bc_topol0035_RiverLake]
  ### NOT DEFINING [bc_drain0035_Pit]
  [bc_drain0035_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F035_Y2022_M06'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0036: F036_Y2022_M07
  ### NOT DEFINING [bc_topol0036_LGFlux]
  ### NOT DEFINING [bc_topol0036_Rainfall]
  ### NOT DEFINING [bc_topol0036_RiverLake]
  ### NOT DEFINING [bc_drain0036_Pit]
  [bc_drain0036_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F036_Y2022_M07'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0037: F037_Y2022_M08
  ### NOT DEFINING [bc_topol0037_LGFlux]
  ### NOT DEFINING [bc_topol0037_Rainfall]
  ### NOT DEFINING [bc_topol0037_RiverLake]
  ### NOT DEFINING [bc_drain0037_Pit]
  [bc_drain0037_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F037_Y2022_M08'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0038: F038_Y2022_M09
  ### NOT DEFINING [bc_topol0038_LGFlux]
  ### NOT DEFINING [bc_topol0038_Rainfall]
  ### NOT DEFINING [bc_topol0038_RiverLake]
  ### NOT DEFINING [bc_drain0038_Pit]
  [bc_drain0038_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F038_Y2022_M09'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0039: F039_Y2022_M10
  ### NOT DEFINING [bc_topol0039_LGFlux]
  ### NOT DEFINING [bc_topol0039_Rainfall]
  ### NOT DEFINING [bc_topol0039_RiverLake]
  ### NOT DEFINING [bc_drain0039_Pit]
  [bc_drain0039_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F039_Y2022_M10'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0040: F040_Y2022_M11
  ### NOT DEFINING [bc_topol0040_LGFlux]
  ### NOT DEFINING [bc_topol0040_Rainfall]
  ### NOT DEFINING [bc_topol0040_RiverLake]
  ### NOT DEFINING [bc_drain0040_Pit]
  [bc_drain0040_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F040_Y2022_M11'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0041: F041_Y2022_M12
  ### NOT DEFINING [bc_topol0041_LGFlux]
  ### NOT DEFINING [bc_topol0041_Rainfall]
  ### NOT DEFINING [bc_topol0041_RiverLake]
  ### NOT DEFINING [bc_drain0041_Pit]
  [bc_drain0041_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F041_Y2022_M12'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0042: F042_Y2023_M01
  ### NOT DEFINING [bc_topol0042_LGFlux]
  ### NOT DEFINING [bc_topol0042_Rainfall]
  ### NOT DEFINING [bc_topol0042_RiverLake]
  ### NOT DEFINING [bc_drain0042_Pit]
  [bc_drain0042_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F042_Y2023_M01'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0043: F043_Y2023_M02
  ### NOT DEFINING [bc_topol0043_LGFlux]
  ### NOT DEFINING [bc_topol0043_Rainfall]
  ### NOT DEFINING [bc_topol0043_RiverLake]
  ### NOT DEFINING [bc_drain0043_Pit]
  [bc_drain0043_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F043_Y2023_M02'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0044: F044_Y2023_M03
  ### NOT DEFINING [bc_topol0044_LGFlux]
  ### NOT DEFINING [bc_topol0044_Rainfall]
  ### NOT DEFINING [bc_topol0044_RiverLake]
  ### NOT DEFINING [bc_drain0044_Pit]
  [bc_drain0044_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F044_Y2023_M03'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0045: F045_Y2023_M04
  ### NOT DEFINING [bc_topol0045_LGFlux]
  ### NOT DEFINING [bc_topol0045_Rainfall]
  ### NOT DEFINING [bc_topol0045_RiverLake]
  ### NOT DEFINING [bc_drain0045_Pit]
  [bc_drain0045_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F045_Y2023_M04'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0046: F046_Y2023_M05
  ### NOT DEFINING [bc_topol0046_LGFlux]
  ### NOT DEFINING [bc_topol0046_Rainfall]
  ### NOT DEFINING [bc_topol0046_RiverLake]
  ### NOT DEFINING [bc_drain0046_Pit]
  [bc_drain0046_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F046_Y2023_M05'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0047: F047_Y2023_M06
  ### NOT DEFINING [bc_topol0047_LGFlux]
  ### NOT DEFINING [bc_topol0047_Rainfall]
  ### NOT DEFINING [bc_topol0047_RiverLake]
  ### NOT DEFINING [bc_drain0047_Pit]
  [bc_drain0047_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F047_Y2023_M06'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0048: F048_Y2023_M07
  ### NOT DEFINING [bc_topol0048_LGFlux]
  ### NOT DEFINING [bc_topol0048_Rainfall]
  ### NOT DEFINING [bc_topol0048_RiverLake]
  ### NOT DEFINING [bc_drain0048_Pit]
  [bc_drain0048_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F048_Y2023_M07'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0049: F049_Y2023_M08
  ### NOT DEFINING [bc_topol0049_LGFlux]
  ### NOT DEFINING [bc_topol0049_Rainfall]
  ### NOT DEFINING [bc_topol0049_RiverLake]
  ### NOT DEFINING [bc_drain0049_Pit]
  [bc_drain0049_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F049_Y2023_M08'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0050: F050_Y2023_M09
  ### NOT DEFINING [bc_topol0050_LGFlux]
  ### NOT DEFINING [bc_topol0050_Rainfall]
  ### NOT DEFINING [bc_topol0050_RiverLake]
  ### NOT DEFINING [bc_drain0050_Pit]
  [bc_drain0050_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F050_Y2023_M09'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0051: F051_Y2023_M10
  ### NOT DEFINING [bc_topol0051_LGFlux]
  ### NOT DEFINING [bc_topol0051_Rainfall]
  ### NOT DEFINING [bc_topol0051_RiverLake]
  ### NOT DEFINING [bc_drain0051_Pit]
  [bc_drain0051_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F051_Y2023_M10'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0052: F052_Y2023_M11
  ### NOT DEFINING [bc_topol0052_LGFlux]
  ### NOT DEFINING [bc_topol0052_Rainfall]
  ### NOT DEFINING [bc_topol0052_RiverLake]
  ### NOT DEFINING [bc_drain0052_Pit]
  [bc_drain0052_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F052_Y2023_M11'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0053: F053_Y2023_M12
  ### NOT DEFINING [bc_topol0053_LGFlux]
  ### NOT DEFINING [bc_topol0053_Rainfall]
  ### NOT DEFINING [bc_topol0053_RiverLake]
  ### NOT DEFINING [bc_drain0053_Pit]
  [bc_drain0053_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F053_Y2023_M12'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0054: F054_Y2024_M01
  ### NOT DEFINING [bc_topol0054_LGFlux]
  ### NOT DEFINING [bc_topol0054_Rainfall]
  ### NOT DEFINING [bc_topol0054_RiverLake]
  ### NOT DEFINING [bc_drain0054_Pit]
  [bc_drain0054_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F054_Y2024_M01'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0055: F055_Y2024_M02
  ### NOT DEFINING [bc_topol0055_LGFlux]
  ### NOT DEFINING [bc_topol0055_Rainfall]
  ### NOT DEFINING [bc_topol0055_RiverLake]
  ### NOT DEFINING [bc_drain0055_Pit]
  [bc_drain0055_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F055_Y2024_M02'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0056: F056_Y2024_M03
  ### NOT DEFINING [bc_topol0056_LGFlux]
  ### NOT DEFINING [bc_topol0056_Rainfall]
  ### NOT DEFINING [bc_topol0056_RiverLake]
  ### NOT DEFINING [bc_drain0056_Pit]
  [bc_drain0056_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F056_Y2024_M03'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0057: F057_Y2024_M04
  ### NOT DEFINING [bc_topol0057_LGFlux]
  ### NOT DEFINING [bc_topol0057_Rainfall]
  ### NOT DEFINING [bc_topol0057_RiverLake]
  ### NOT DEFINING [bc_drain0057_Pit]
  [bc_drain0057_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F057_Y2024_M04'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0058: F058_Y2024_M05
  ### NOT DEFINING [bc_topol0058_LGFlux]
  ### NOT DEFINING [bc_topol0058_Rainfall]
  ### NOT DEFINING [bc_topol0058_RiverLake]
  ### NOT DEFINING [bc_drain0058_Pit]
  [bc_drain0058_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F058_Y2024_M05'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0059: F059_Y2024_M06
  ### NOT DEFINING [bc_topol0059_LGFlux]
  ### NOT DEFINING [bc_topol0059_Rainfall]
  ### NOT DEFINING [bc_topol0059_RiverLake]
  ### NOT DEFINING [bc_drain0059_Pit]
  [bc_drain0059_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F059_Y2024_M06'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0060: F060_Y2024_M07
  ### NOT DEFINING [bc_topol0060_LGFlux]
  ### NOT DEFINING [bc_topol0060_Rainfall]
  ### NOT DEFINING [bc_topol0060_RiverLake]
  ### NOT DEFINING [bc_drain0060_Pit]
  [bc_drain0060_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F060_Y2024_M07'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0061: F061_Y2024_M08
  ### NOT DEFINING [bc_topol0061_LGFlux]
  ### NOT DEFINING [bc_topol0061_Rainfall]
  ### NOT DEFINING [bc_topol0061_RiverLake]
  ### NOT DEFINING [bc_drain0061_Pit]
  [bc_drain0061_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F061_Y2024_M08'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0062: F062_Y2024_M09
  ### NOT DEFINING [bc_topol0062_LGFlux]
  ### NOT DEFINING [bc_topol0062_Rainfall]
  ### NOT DEFINING [bc_topol0062_RiverLake]
  ### NOT DEFINING [bc_drain0062_Pit]
  [bc_drain0062_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F062_Y2024_M09'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0063: F063_Y2024_M10
  ### NOT DEFINING [bc_topol0063_LGFlux]
  ### NOT DEFINING [bc_topol0063_Rainfall]
  ### NOT DEFINING [bc_topol0063_RiverLake]
  ### NOT DEFINING [bc_drain0063_Pit]
  [bc_drain0063_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F063_Y2024_M10'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0064: F064_Y2024_M11
  ### NOT DEFINING [bc_topol0064_LGFlux]
  ### NOT DEFINING [bc_topol0064_Rainfall]
  ### NOT DEFINING [bc_topol0064_RiverLake]
  ### NOT DEFINING [bc_drain0064_Pit]
  [bc_drain0064_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F064_Y2024_M11'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0065: F065_Y2024_M12
  ### NOT DEFINING [bc_topol0065_LGFlux]
  ### NOT DEFINING [bc_topol0065_Rainfall]
  ### NOT DEFINING [bc_topol0065_RiverLake]
  ### NOT DEFINING [bc_drain0065_Pit]
  [bc_drain0065_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F065_Y2024_M12'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0066: F066_Y2025_M01
  ### NOT DEFINING [bc_topol0066_LGFlux]
  ### NOT DEFINING [bc_topol0066_Rainfall]
  ### NOT DEFINING [bc_topol0066_RiverLake]
  ### NOT DEFINING [bc_drain0066_Pit]
  [bc_drain0066_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F066_Y2025_M01'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0067: F067_Y2025_M02
  ### NOT DEFINING [bc_topol0067_LGFlux]
  ### NOT DEFINING [bc_topol0067_Rainfall]
  ### NOT DEFINING [bc_topol0067_RiverLake]
  ### NOT DEFINING [bc_drain0067_Pit]
  [bc_drain0067_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F067_Y2025_M02'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0068: F068_Y2025_M03
  ### NOT DEFINING [bc_topol0068_LGFlux]
  ### NOT DEFINING [bc_topol0068_Rainfall]
  ### NOT DEFINING [bc_topol0068_RiverLake]
  ### NOT DEFINING [bc_drain0068_Pit]
  [bc_drain0068_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F068_Y2025_M03'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0069: F069_Y2025_M04
  ### NOT DEFINING [bc_topol0069_LGFlux]
  ### NOT DEFINING [bc_topol0069_Rainfall]
  ### NOT DEFINING [bc_topol0069_RiverLake]
  ### NOT DEFINING [bc_drain0069_Pit]
  [bc_drain0069_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F069_Y2025_M04'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0070: F070_Y2025_M05
  ### NOT DEFINING [bc_topol0070_LGFlux]
  ### NOT DEFINING [bc_topol0070_Rainfall]
  ### NOT DEFINING [bc_topol0070_RiverLake]
  ### NOT DEFINING [bc_drain0070_Pit]
  [bc_drain0070_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F070_Y2025_M05'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0071: F071_Y2025_M06
  ### NOT DEFINING [bc_topol0071_LGFlux]
  ### NOT DEFINING [bc_topol0071_Rainfall]
  ### NOT DEFINING [bc_topol0071_RiverLake]
  ### NOT DEFINING [bc_drain0071_Pit]
  [bc_drain0071_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F071_Y2025_M06'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0072: F072_Y2025_M07
  ### NOT DEFINING [bc_topol0072_LGFlux]
  ### NOT DEFINING [bc_topol0072_Rainfall]
  ### NOT DEFINING [bc_topol0072_RiverLake]
  ### NOT DEFINING [bc_drain0072_Pit]
  [bc_drain0072_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F072_Y2025_M07'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0073: F073_Y2025_M08
  ### NOT DEFINING [bc_topol0073_LGFlux]
  ### NOT DEFINING [bc_topol0073_Rainfall]
  ### NOT DEFINING [bc_topol0073_RiverLake]
  ### NOT DEFINING [bc_drain0073_Pit]
  [bc_drain0073_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F073_Y2025_M08'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0074: F074_Y2025_M09
  ### NOT DEFINING [bc_topol0074_LGFlux]
  ### NOT DEFINING [bc_topol0074_Rainfall]
  ### NOT DEFINING [bc_topol0074_RiverLake]
  ### NOT DEFINING [bc_drain0074_Pit]
  [bc_drain0074_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F074_Y2025_M09'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0075: F075_Y2025_M10
  ### NOT DEFINING [bc_topol0075_LGFlux]
  ### NOT DEFINING [bc_topol0075_Rainfall]
  ### NOT DEFINING [bc_topol0075_RiverLake]
  ### NOT DEFINING [bc_drain0075_Pit]
  [bc_drain0075_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F075_Y2025_M10'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0076: F076_Y2025_M11
  ### NOT DEFINING [bc_topol0076_LGFlux]
  ### NOT DEFINING [bc_topol0076_Rainfall]
  ### NOT DEFINING [bc_topol0076_RiverLake]
  ### NOT DEFINING [bc_drain0076_Pit]
  [bc_drain0076_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F076_Y2025_M11'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0077: F077_Y2025_M12
  ### NOT DEFINING [bc_topol0077_LGFlux]
  ### NOT DEFINING [bc_topol0077_Rainfall]
  ### NOT DEFINING [bc_topol0077_RiverLake]
  ### NOT DEFINING [bc_drain0077_Pit]
  [bc_drain0077_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F077_Y2025_M12'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0078: F078_Y2026_M01
  ### NOT DEFINING [bc_topol0078_LGFlux]
  ### NOT DEFINING [bc_topol0078_Rainfall]
  ### NOT DEFINING [bc_topol0078_RiverLake]
  ### NOT DEFINING [bc_drain0078_Pit]
  [bc_drain0078_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F078_Y2026_M01'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0079: F079_Y2026_M02
  ### NOT DEFINING [bc_topol0079_LGFlux]
  ### NOT DEFINING [bc_topol0079_Rainfall]
  ### NOT DEFINING [bc_topol0079_RiverLake]
  ### NOT DEFINING [bc_drain0079_Pit]
  [bc_drain0079_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F079_Y2026_M02'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0080: F080_Y2026_M03
  ### NOT DEFINING [bc_topol0080_LGFlux]
  ### NOT DEFINING [bc_topol0080_Rainfall]
  ### NOT DEFINING [bc_topol0080_RiverLake]
  ### NOT DEFINING [bc_drain0080_Pit]
  [bc_drain0080_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F080_Y2026_M03'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0081: F081_Y2026_M04
  ### NOT DEFINING [bc_topol0081_LGFlux]
  ### NOT DEFINING [bc_topol0081_Rainfall]
  ### NOT DEFINING [bc_topol0081_RiverLake]
  ### NOT DEFINING [bc_drain0081_Pit]
  [bc_drain0081_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F081_Y2026_M04'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0082: F082_Y2026_M05
  ### NOT DEFINING [bc_topol0082_LGFlux]
  ### NOT DEFINING [bc_topol0082_Rainfall]
  ### NOT DEFINING [bc_topol0082_RiverLake]
  ### NOT DEFINING [bc_drain0082_Pit]
  [bc_drain0082_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F082_Y2026_M05'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0083: F083_Y2026_M06
  ### NOT DEFINING [bc_topol0083_LGFlux]
  ### NOT DEFINING [bc_topol0083_Rainfall]
  ### NOT DEFINING [bc_topol0083_RiverLake]
  ### NOT DEFINING [bc_drain0083_Pit]
  [bc_drain0083_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F083_Y2026_M06'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0084: F084_Y2026_M07
  ### NOT DEFINING [bc_topol0084_LGFlux]
  ### NOT DEFINING [bc_topol0084_Rainfall]
  ### NOT DEFINING [bc_topol0084_RiverLake]
  ### NOT DEFINING [bc_drain0084_Pit]
  [bc_drain0084_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F084_Y2026_M07'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0085: F085_Y2026_M08
  ### NOT DEFINING [bc_topol0085_LGFlux]
  ### NOT DEFINING [bc_topol0085_Rainfall]
  ### NOT DEFINING [bc_topol0085_RiverLake]
  ### NOT DEFINING [bc_drain0085_Pit]
  [bc_drain0085_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F085_Y2026_M08'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0086: F086_Y2026_M09
  ### NOT DEFINING [bc_topol0086_LGFlux]
  ### NOT DEFINING [bc_topol0086_Rainfall]
  ### NOT DEFINING [bc_topol0086_RiverLake]
  ### NOT DEFINING [bc_drain0086_Pit]
  [bc_drain0086_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F086_Y2026_M09'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0087: F087_Y2026_M10
  ### NOT DEFINING [bc_topol0087_LGFlux]
  ### NOT DEFINING [bc_topol0087_Rainfall]
  ### NOT DEFINING [bc_topol0087_RiverLake]
  ### NOT DEFINING [bc_drain0087_Pit]
  [bc_drain0087_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F087_Y2026_M10'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0088: F088_Y2026_M11
  ### NOT DEFINING [bc_topol0088_LGFlux]
  ### NOT DEFINING [bc_topol0088_Rainfall]
  ### NOT DEFINING [bc_topol0088_RiverLake]
  ### NOT DEFINING [bc_drain0088_Pit]
  [bc_drain0088_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F088_Y2026_M11'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0089: F089_Y2026_M12
  ### NOT DEFINING [bc_topol0089_LGFlux]
  ### NOT DEFINING [bc_topol0089_Rainfall]
  ### NOT DEFINING [bc_topol0089_RiverLake]
  ### NOT DEFINING [bc_drain0089_Pit]
  [bc_drain0089_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F089_Y2026_M12'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0090: F090_Y2027_M01
  ### NOT DEFINING [bc_topol0090_LGFlux]
  ### NOT DEFINING [bc_topol0090_Rainfall]
  ### NOT DEFINING [bc_topol0090_RiverLake]
  ### NOT DEFINING [bc_drain0090_Pit]
  [bc_drain0090_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F090_Y2027_M01'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0091: F091_Y2027_M02
  ### NOT DEFINING [bc_topol0091_LGFlux]
  ### NOT DEFINING [bc_topol0091_Rainfall]
  ### NOT DEFINING [bc_topol0091_RiverLake]
  ### NOT DEFINING [bc_drain0091_Pit]
  [bc_drain0091_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F091_Y2027_M02'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0092: F092_Y2027_M03
  ### NOT DEFINING [bc_topol0092_LGFlux]
  ### NOT DEFINING [bc_topol0092_Rainfall]
  ### NOT DEFINING [bc_topol0092_RiverLake]
  ### NOT DEFINING [bc_drain0092_Pit]
  [bc_drain0092_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F092_Y2027_M03'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0093: F093_Y2027_M04
  ### NOT DEFINING [bc_topol0093_LGFlux]
  ### NOT DEFINING [bc_topol0093_Rainfall]
  ### NOT DEFINING [bc_topol0093_RiverLake]
  ### NOT DEFINING [bc_drain0093_Pit]
  [bc_drain0093_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F093_Y2027_M04'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0094: F094_Y2027_M05
  ### NOT DEFINING [bc_topol0094_LGFlux]
  ### NOT DEFINING [bc_topol0094_Rainfall]
  ### NOT DEFINING [bc_topol0094_RiverLake]
  ### NOT DEFINING [bc_drain0094_Pit]
  [bc_drain0094_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F094_Y2027_M05'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0095: F095_Y2027_M06
  ### NOT DEFINING [bc_topol0095_LGFlux]
  ### NOT DEFINING [bc_topol0095_Rainfall]
  ### NOT DEFINING [bc_topol0095_RiverLake]
  ### NOT DEFINING [bc_drain0095_Pit]
  [bc_drain0095_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F095_Y2027_M06'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0096: F096_Y2027_M07
  ### NOT DEFINING [bc_topol0096_LGFlux]
  ### NOT DEFINING [bc_topol0096_Rainfall]
  ### NOT DEFINING [bc_topol0096_RiverLake]
  ### NOT DEFINING [bc_drain0096_Pit]
  [bc_drain0096_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F096_Y2027_M07'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0097: F097_Y2027_M08
  ### NOT DEFINING [bc_topol0097_LGFlux]
  ### NOT DEFINING [bc_topol0097_Rainfall]
  ### NOT DEFINING [bc_topol0097_RiverLake]
  ### NOT DEFINING [bc_drain0097_Pit]
  [bc_drain0097_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F097_Y2027_M08'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0098: F098_Y2027_M09
  ### NOT DEFINING [bc_topol0098_LGFlux]
  ### NOT DEFINING [bc_topol0098_Rainfall]
  ### NOT DEFINING [bc_topol0098_RiverLake]
  ### NOT DEFINING [bc_drain0098_Pit]
  [bc_drain0098_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F098_Y2027_M09'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0099: F099_Y2027_M10
  ### NOT DEFINING [bc_topol0099_LGFlux]
  ### NOT DEFINING [bc_topol0099_Rainfall]
  ### NOT DEFINING [bc_topol0099_RiverLake]
  ### NOT DEFINING [bc_drain0099_Pit]
  [bc_drain0099_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F099_Y2027_M10'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0100: F100_Y2027_M11
  ### NOT DEFINING [bc_topol0100_LGFlux]
  ### NOT DEFINING [bc_topol0100_Rainfall]
  ### NOT DEFINING [bc_topol0100_RiverLake]
  ### NOT DEFINING [bc_drain0100_Pit]
  [bc_drain0100_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F100_Y2027_M11'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0101: F101_Y2027_M12
  ### NOT DEFINING [bc_topol0101_LGFlux]
  ### NOT DEFINING [bc_topol0101_Rainfall]
  ### NOT DEFINING [bc_topol0101_RiverLake]
  ### NOT DEFINING [bc_drain0101_Pit]
  [bc_drain0101_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F101_Y2027_M12'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0102: F102_Y2028_M01
  ### NOT DEFINING [bc_topol0102_LGFlux]
  ### NOT DEFINING [bc_topol0102_Rainfall]
  ### NOT DEFINING [bc_topol0102_RiverLake]
  ### NOT DEFINING [bc_drain0102_Pit]
  [bc_drain0102_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F102_Y2028_M01'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0103: F103_Y2028_M02
  ### NOT DEFINING [bc_topol0103_LGFlux]
  ### NOT DEFINING [bc_topol0103_Rainfall]
  ### NOT DEFINING [bc_topol0103_RiverLake]
  ### NOT DEFINING [bc_drain0103_Pit]
  [bc_drain0103_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F103_Y2028_M02'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0104: F104_Y2028_M03
  ### NOT DEFINING [bc_topol0104_LGFlux]
  ### NOT DEFINING [bc_topol0104_Rainfall]
  ### NOT DEFINING [bc_topol0104_RiverLake]
  ### NOT DEFINING [bc_drain0104_Pit]
  [bc_drain0104_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F104_Y2028_M03'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0105: F105_Y2028_M04
  ### NOT DEFINING [bc_topol0105_LGFlux]
  ### NOT DEFINING [bc_topol0105_Rainfall]
  ### NOT DEFINING [bc_topol0105_RiverLake]
  ### NOT DEFINING [bc_drain0105_Pit]
  [bc_drain0105_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F105_Y2028_M04'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0106: F106_Y2028_M05
  ### NOT DEFINING [bc_topol0106_LGFlux]
  ### NOT DEFINING [bc_topol0106_Rainfall]
  ### NOT DEFINING [bc_topol0106_RiverLake]
  ### NOT DEFINING [bc_drain0106_Pit]
  [bc_drain0106_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F106_Y2028_M05'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0107: F107_Y2028_M06
  ### NOT DEFINING [bc_topol0107_LGFlux]
  ### NOT DEFINING [bc_topol0107_Rainfall]
  ### NOT DEFINING [bc_topol0107_RiverLake]
  ### NOT DEFINING [bc_drain0107_Pit]
  [bc_drain0107_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F107_Y2028_M06'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0108: F108_Y2028_M07
  ### NOT DEFINING [bc_topol0108_LGFlux]
  ### NOT DEFINING [bc_topol0108_Rainfall]
  ### NOT DEFINING [bc_topol0108_RiverLake]
  ### NOT DEFINING [bc_drain0108_Pit]
  [bc_drain0108_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F108_Y2028_M07'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0109: F109_Y2028_M08
  ### NOT DEFINING [bc_topol0109_LGFlux]
  ### NOT DEFINING [bc_topol0109_Rainfall]
  ### NOT DEFINING [bc_topol0109_RiverLake]
  ### NOT DEFINING [bc_drain0109_Pit]
  [bc_drain0109_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F109_Y2028_M08'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0110: F110_Y2028_M09
  ### NOT DEFINING [bc_topol0110_LGFlux]
  ### NOT DEFINING [bc_topol0110_Rainfall]
  ### NOT DEFINING [bc_topol0110_RiverLake]
  ### NOT DEFINING [bc_drain0110_Pit]
  [bc_drain0110_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F110_Y2028_M09'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0111: F111_Y2028_M10
  ### NOT DEFINING [bc_topol0111_LGFlux]
  ### NOT DEFINING [bc_topol0111_Rainfall]
  ### NOT DEFINING [bc_topol0111_RiverLake]
  ### NOT DEFINING [bc_drain0111_Pit]
  [bc_drain0111_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F111_Y2028_M10'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0112: F112_Y2028_M11
  ### NOT DEFINING [bc_topol0112_LGFlux]
  ### NOT DEFINING [bc_topol0112_Rainfall]
  ### NOT DEFINING [bc_topol0112_RiverLake]
  ### NOT DEFINING [bc_drain0112_Pit]
  [bc_drain0112_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F112_Y2028_M11'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0113: F113_Y2028_M12
  ### NOT DEFINING [bc_topol0113_LGFlux]
  ### NOT DEFINING [bc_topol0113_Rainfall]
  ### NOT DEFINING [bc_topol0113_RiverLake]
  ### NOT DEFINING [bc_drain0113_Pit]
  [bc_drain0113_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F113_Y2028_M12'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0114: F114_Y2029_M01
  ### NOT DEFINING [bc_topol0114_LGFlux]
  ### NOT DEFINING [bc_topol0114_Rainfall]
  ### NOT DEFINING [bc_topol0114_RiverLake]
  ### NOT DEFINING [bc_drain0114_Pit]
  [bc_drain0114_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F114_Y2029_M01'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0115: F115_Y2029_M02
  ### NOT DEFINING [bc_topol0115_LGFlux]
  ### NOT DEFINING [bc_topol0115_Rainfall]
  ### NOT DEFINING [bc_topol0115_RiverLake]
  ### NOT DEFINING [bc_drain0115_Pit]
  [bc_drain0115_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F115_Y2029_M02'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0116: F116_Y2029_M03
  ### NOT DEFINING [bc_topol0116_LGFlux]
  ### NOT DEFINING [bc_topol0116_Rainfall]
  ### NOT DEFINING [bc_topol0116_RiverLake]
  ### NOT DEFINING [bc_drain0116_Pit]
  [bc_drain0116_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F116_Y2029_M03'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0117: F117_Y2029_M04
  ### NOT DEFINING [bc_topol0117_LGFlux]
  ### NOT DEFINING [bc_topol0117_Rainfall]
  ### NOT DEFINING [bc_topol0117_RiverLake]
  ### NOT DEFINING [bc_drain0117_Pit]
  [bc_drain0117_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F117_Y2029_M04'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0118: F118_Y2029_M05
  ### NOT DEFINING [bc_topol0118_LGFlux]
  ### NOT DEFINING [bc_topol0118_Rainfall]
  ### NOT DEFINING [bc_topol0118_RiverLake]
  ### NOT DEFINING [bc_drain0118_Pit]
  [bc_drain0118_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F118_Y2029_M05'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0119: F119_Y2029_M06
  ### NOT DEFINING [bc_topol0119_LGFlux]
  ### NOT DEFINING [bc_topol0119_Rainfall]
  ### NOT DEFINING [bc_topol0119_RiverLake]
  ### NOT DEFINING [bc_drain0119_Pit]
  [bc_drain0119_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F119_Y2029_M06'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0120: F120_Y2029_M07
  ### NOT DEFINING [bc_topol0120_LGFlux]
  ### NOT DEFINING [bc_topol0120_Rainfall]
  ### NOT DEFINING [bc_topol0120_RiverLake]
  ### NOT DEFINING [bc_drain0120_Pit]
  [bc_drain0120_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F120_Y2029_M07'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0121: F121_Y2029_M08
  ### NOT DEFINING [bc_topol0121_LGFlux]
  ### NOT DEFINING [bc_topol0121_Rainfall]
  ### NOT DEFINING [bc_topol0121_RiverLake]
  ### NOT DEFINING [bc_drain0121_Pit]
  [bc_drain0121_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F121_Y2029_M08'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0122: F122_Y2029_M09
  ### NOT DEFINING [bc_topol0122_LGFlux]
  ### NOT DEFINING [bc_topol0122_Rainfall]
  ### NOT DEFINING [bc_topol0122_RiverLake]
  ### NOT DEFINING [bc_drain0122_Pit]
  [bc_drain0122_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F122_Y2029_M09'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0123: F123_Y2029_M10
  ### NOT DEFINING [bc_topol0123_LGFlux]
  ### NOT DEFINING [bc_topol0123_Rainfall]
  ### NOT DEFINING [bc_topol0123_RiverLake]
  ### NOT DEFINING [bc_drain0123_Pit]
  [bc_drain0123_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F123_Y2029_M10'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0124: F124_Y2029_M11
  ### NOT DEFINING [bc_topol0124_LGFlux]
  ### NOT DEFINING [bc_topol0124_Rainfall]
  ### NOT DEFINING [bc_topol0124_RiverLake]
  ### NOT DEFINING [bc_drain0124_Pit]
  [bc_drain0124_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F124_Y2029_M11'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0125: F125_Y2029_M12
  ### NOT DEFINING [bc_topol0125_LGFlux]
  ### NOT DEFINING [bc_topol0125_Rainfall]
  ### NOT DEFINING [bc_topol0125_RiverLake]
  ### NOT DEFINING [bc_drain0125_Pit]
  [bc_drain0125_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F125_Y2029_M12'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0126: F126_Y2030_M01
  ### NOT DEFINING [bc_topol0126_LGFlux]
  ### NOT DEFINING [bc_topol0126_Rainfall]
  ### NOT DEFINING [bc_topol0126_RiverLake]
  ### NOT DEFINING [bc_drain0126_Pit]
  [bc_drain0126_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F126_Y2030_M01'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0127: F127_Y2030_M02
  ### NOT DEFINING [bc_topol0127_LGFlux]
  ### NOT DEFINING [bc_topol0127_Rainfall]
  ### NOT DEFINING [bc_topol0127_RiverLake]
  ### NOT DEFINING [bc_drain0127_Pit]
  [bc_drain0127_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F127_Y2030_M02'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0128: F128_Y2030_M03
  ### NOT DEFINING [bc_topol0128_LGFlux]
  ### NOT DEFINING [bc_topol0128_Rainfall]
  ### NOT DEFINING [bc_topol0128_RiverLake]
  ### NOT DEFINING [bc_drain0128_Pit]
  [bc_drain0128_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F128_Y2030_M03'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0129: F129_Y2030_M04
  ### NOT DEFINING [bc_topol0129_LGFlux]
  ### NOT DEFINING [bc_topol0129_Rainfall]
  ### NOT DEFINING [bc_topol0129_RiverLake]
  ### NOT DEFINING [bc_drain0129_Pit]
  [bc_drain0129_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F129_Y2030_M04'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0130: F130_Y2030_M05
  ### NOT DEFINING [bc_topol0130_LGFlux]
  ### NOT DEFINING [bc_topol0130_Rainfall]
  ### NOT DEFINING [bc_topol0130_RiverLake]
  ### NOT DEFINING [bc_drain0130_Pit]
  [bc_drain0130_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F130_Y2030_M05'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0131: F131_Y2030_M06
  ### NOT DEFINING [bc_topol0131_LGFlux]
  ### NOT DEFINING [bc_topol0131_Rainfall]
  ### NOT DEFINING [bc_topol0131_RiverLake]
  ### NOT DEFINING [bc_drain0131_Pit]
  [bc_drain0131_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F131_Y2030_M06'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0132: F132_Y2030_M07
  ### NOT DEFINING [bc_topol0132_LGFlux]
  ### NOT DEFINING [bc_topol0132_Rainfall]
  ### NOT DEFINING [bc_topol0132_RiverLake]
  ### NOT DEFINING [bc_drain0132_Pit]
  [bc_drain0132_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F132_Y2030_M07'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0133: F133_Y2030_M08
  ### NOT DEFINING [bc_topol0133_LGFlux]
  ### NOT DEFINING [bc_topol0133_Rainfall]
  ### NOT DEFINING [bc_topol0133_RiverLake]
  ### NOT DEFINING [bc_drain0133_Pit]
  [bc_drain0133_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F133_Y2030_M08'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0134: F134_Y2030_M09
  ### NOT DEFINING [bc_topol0134_LGFlux]
  ### NOT DEFINING [bc_topol0134_Rainfall]
  ### NOT DEFINING [bc_topol0134_RiverLake]
  ### NOT DEFINING [bc_drain0134_Pit]
  [bc_drain0134_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F134_Y2030_M09'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0135: F135_Y2030_M10
  ### NOT DEFINING [bc_topol0135_LGFlux]
  ### NOT DEFINING [bc_topol0135_Rainfall]
  ### NOT DEFINING [bc_topol0135_RiverLake]
  ### NOT DEFINING [bc_drain0135_Pit]
  [bc_drain0135_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F135_Y2030_M10'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0136: F136_Y2030_M11
  ### NOT DEFINING [bc_topol0136_LGFlux]
  ### NOT DEFINING [bc_topol0136_Rainfall]
  ### NOT DEFINING [bc_topol0136_RiverLake]
  ### NOT DEFINING [bc_drain0136_Pit]
  [bc_drain0136_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F136_Y2030_M11'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0137: F137_Y2030_M12
  ### NOT DEFINING [bc_topol0137_LGFlux]
  ### NOT DEFINING [bc_topol0137_Rainfall]
  ### NOT DEFINING [bc_topol0137_RiverLake]
  ### NOT DEFINING [bc_drain0137_Pit]
  [bc_drain0137_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F137_Y2030_M12'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0138: F138_Y2031_M01
  ### NOT DEFINING [bc_topol0138_LGFlux]
  ### NOT DEFINING [bc_topol0138_Rainfall]
  ### NOT DEFINING [bc_topol0138_RiverLake]
  ### NOT DEFINING [bc_drain0138_Pit]
  [bc_drain0138_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F138_Y2031_M01'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0139: F139_Y2031_M02
  ### NOT DEFINING [bc_topol0139_LGFlux]
  ### NOT DEFINING [bc_topol0139_Rainfall]
  ### NOT DEFINING [bc_topol0139_RiverLake]
  ### NOT DEFINING [bc_drain0139_Pit]
  [bc_drain0139_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F139_Y2031_M02'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0140: F140_Y2031_M03
  ### NOT DEFINING [bc_topol0140_LGFlux]
  ### NOT DEFINING [bc_topol0140_Rainfall]
  ### NOT DEFINING [bc_topol0140_RiverLake]
  ### NOT DEFINING [bc_drain0140_Pit]
  [bc_drain0140_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F140_Y2031_M03'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0141: F141_Y2031_M04
  ### NOT DEFINING [bc_topol0141_LGFlux]
  ### NOT DEFINING [bc_topol0141_Rainfall]
  ### NOT DEFINING [bc_topol0141_RiverLake]
  ### NOT DEFINING [bc_drain0141_Pit]
  [bc_drain0141_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F141_Y2031_M04'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0142: F142_Y2031_M05
  ### NOT DEFINING [bc_topol0142_LGFlux]
  ### NOT DEFINING [bc_topol0142_Rainfall]
  ### NOT DEFINING [bc_topol0142_RiverLake]
  ### NOT DEFINING [bc_drain0142_Pit]
  [bc_drain0142_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F142_Y2031_M05'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0143: F143_Y2031_M06
  ### NOT DEFINING [bc_topol0143_LGFlux]
  ### NOT DEFINING [bc_topol0143_Rainfall]
  ### NOT DEFINING [bc_topol0143_RiverLake]
  ### NOT DEFINING [bc_drain0143_Pit]
  [bc_drain0143_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F143_Y2031_M06'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0144: F144_Y2031_M07
  ### NOT DEFINING [bc_topol0144_LGFlux]
  ### NOT DEFINING [bc_topol0144_Rainfall]
  ### NOT DEFINING [bc_topol0144_RiverLake]
  ### NOT DEFINING [bc_drain0144_Pit]
  [bc_drain0144_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F144_Y2031_M07'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0145: F145_Y2031_M08
  ### NOT DEFINING [bc_topol0145_LGFlux]
  ### NOT DEFINING [bc_topol0145_Rainfall]
  ### NOT DEFINING [bc_topol0145_RiverLake]
  ### NOT DEFINING [bc_drain0145_Pit]
  [bc_drain0145_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F145_Y2031_M08'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0146: F146_Y2031_M09
  ### NOT DEFINING [bc_topol0146_LGFlux]
  ### NOT DEFINING [bc_topol0146_Rainfall]
  ### NOT DEFINING [bc_topol0146_RiverLake]
  ### NOT DEFINING [bc_drain0146_Pit]
  [bc_drain0146_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F146_Y2031_M09'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0147: F147_Y2031_M10
  ### NOT DEFINING [bc_topol0147_LGFlux]
  ### NOT DEFINING [bc_topol0147_Rainfall]
  ### NOT DEFINING [bc_topol0147_RiverLake]
  ### NOT DEFINING [bc_drain0147_Pit]
  [bc_drain0147_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F147_Y2031_M10'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0148: F148_Y2031_M11
  ### NOT DEFINING [bc_topol0148_LGFlux]
  ### NOT DEFINING [bc_topol0148_Rainfall]
  ### NOT DEFINING [bc_topol0148_RiverLake]
  ### NOT DEFINING [bc_drain0148_Pit]
  [bc_drain0148_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F148_Y2031_M11'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0149: F149_Y2031_M12
  ### NOT DEFINING [bc_topol0149_LGFlux]
  ### NOT DEFINING [bc_topol0149_Rainfall]
  ### NOT DEFINING [bc_topol0149_RiverLake]
  ### NOT DEFINING [bc_drain0149_Pit]
  [bc_drain0149_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F149_Y2031_M12'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
  
  ### --- TOPO-RECHARGE AND DRAINAGE FOR 0150: F150_Y2032_M01
  ### NOT DEFINING [bc_topol0150_LGFlux]
  ### NOT DEFINING [bc_topol0150_Rainfall]
  ### NOT DEFINING [bc_topol0150_RiverLake]
  ### NOT DEFINING [bc_drain0150_Pit]
  [bc_drain0150_UG]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure_L1
    boundary = 'DRAINSRFUG_F150_Y2032_M01'
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1
    PT_shift = 0.
    fluid_phase = 0
    use_mobility = true   # default is false
    use_relperm = false    # default is false
  []
[]
[Functions]
  [fcn_pwp_tailingsBC]
    type = ParsedFunction
    expression = '10.'
  []
  [fcn_pwp_pitlakeBC]
    type = ParsedFunction
    expression = '9810.*(10340.-z)'
  []
[]


### ################################################################### #######
### PRECONDITIONING AND EXECUTIONER
### ################################################################### #######

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient   # always use 'Transient' (not 'Steady'), see [Kernels]
  solve_type = NEWTON
  petsc_options = '-snes_converged_reason'
  petsc_options_iname = '-pc_type -pc_hypre_type -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'hypre    boomeramg 0.8'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e4
  nl_max_its = 30
  l_tol = 1e-7
  l_max_its = 100
  line_search = 'bt'
  ###
  ### time stepping:
  # dt = 1.0
  # dtmin = 0.083
  # dtmax = 1.0
  start_time = 1.75
  end_time = 3
  [TimeStepper]
    ### option-1: read time steps from csv-file
    # type = CSVTimeSequenceStepper
    # file_name = 'COWA2023_HR01v2_part1.csv'
    # column_name = time
    ### option-2: define time steps here, directly
    type = TimeSequenceStepper
    time_sequence = '1.000 1.500 1.750 1.833 2.000 2.083 2.167 2.250 2.333 2.417 2.500 2.583 2.667 2.750 2.833 2.917 3.000 3.083 3.167 3.250 3.333 3.417 3.500 3.583 3.667 3.750 3.833 3.917 4.000 4.083 4.167 4.250 4.333 4.417 4.500 4.583 4.667 4.750 4.833 4.917 5.000 5.083 5.167 5.250 5.333 5.417 5.500 5.583 5.667 5.750 5.833 5.917 6.000 6.083 6.167 6.250 6.333 6.417 6.500 6.583 6.667 6.750 6.833 6.917 7.000 7.083 7.167 7.250 7.333 7.417 7.500 7.583 7.667 7.750 7.833 7.917 8.000 8.083 8.167 8.250 8.333 8.417 8.500 8.583 8.667 8.750 8.833 8.917 9.000 9.083 9.167 9.250 9.333 9.417 9.500 9.583 9.667 9.750 9.833 9.917 10.000 10.083 10.167 10.250 10.333 10.417 10.500 10.583 10.667 10.750 10.833 10.917 11.000 11.083 11.167 11.250 11.333 11.417 11.500 11.583 11.667 11.750 11.833 11.917 12.000 12.083 12.167 12.250 12.333 12.417 12.500 12.583 12.667 12.750 12.833 12.917 13.000 13.083 13.167 13.250 13.333 13.417 13.500 13.583 13.667 13.750 13.833 13.917 14.000 14.083'
  []
[]


### ################################################################### #######
### OUTPUT
### ################################################################### #######

[Debug]
  show_var_residual_norms = true
  show_actions = false
  show_action_dependencies = false
[]

[VectorPostprocessors]
  [nodal_pwp]
    type = NodalValueSampler
    variable = porepressure_L1
    sort_by = id
    unique_node_execute = true
  []
[]

[Outputs]
  perf_graph = true
  print_linear_residuals = true
  # print_linear_residuals = false

  [output_exodus_single_file]
    type = ExodusBE
    execute_on = 'TIMESTEP_END'
    sequence = true
    file_base = 'outputExodusTEST_newflux/XXXXXX2023052_R01_HR02_Q02HfullSeq_TR_PS_test03_single'
    times = '1.833	2	2.083	2.167	2.25	2.333	2.417	2.5	2.583	2.667	2.75	2.833	2.917	3	3.083	3.167	3.25	3.333	3.417	3.5	3.583	3.667	3.75	3.833	3.917	4	4.083	4.167	4.25	4.333	4.417	4.5	4.583	4.667	4.75	4.833	4.917	5	5.083	5.167	5.25	5.333	5.417	5.5	5.583	5.667	5.75	5.833	5.917	6	6.083	6.167	6.25	6.333	6.417	6.5	6.583	6.667	6.75	6.833	6.917	7	7.083	7.167	7.25	7.333	7.417	7.5	7.583	7.667	7.75	7.833	7.917	8	8.083	8.167	8.25	8.333	8.417	8.5	8.583	8.667	8.75	8.833	8.917	9	9.083	9.167	9.25	9.333	9.417	9.5	9.583	9.667	9.75	9.833	9.917	10	10.083	10.167	10.25	10.333	10.417	10.5	10.583	10.667	10.75	10.833	10.917	11	11.083	11.167	11.25	11.333	11.417	11.5	11.583	11.667	11.75	11.833	11.917	12	12.083	12.167	12.25	12.333	12.417	12.5	12.583	12.667	12.75	12.833	12.917	13	13.083	13.167	13.25	13.333	13.417	13.5	13.583	13.667	13.75	13.833	13.917	14	14.083  '
    frames = '4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39	40	41	42	43	44	45	46	47	48	49	50	51	52	53	54	55	56	57	58	59	60	61	62	63	64	65	66	67	68	69	70	71	72	73	74	75	76	77	78	79	80	81	82	83	84	85	86	87	88	89	90	91	92	93	94	95	96	97	98	99	100	101	102	103	104	105	106	107	108	109	110	111	112	113	114	115	116	117	118	119	120	121	122	123	124	125	126	127	128	129	130	131	132	133	134	135	136	137	138	139	140	141	142	143	144	145	146	147	148	149	150  '
  []

[]



