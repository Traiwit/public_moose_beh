
[GlobalParams]
  PorousFlowDictator = dictator
[]

[Mesh]
  [mesh]
  type = FileMeshGenerator
   file = HTTRANS2020_G01HTv2_TETL_W2FF.inp
   allow_renumbering = false
   []

  [side_back]
     type = ParsedGenerateSideset
     input = mesh
     combinatorial_geometry = 'z=0'
     new_sideset_name = 'back'
   []

   [side_front]
      type = ParsedGenerateSideset
      input = side_back
      combinatorial_geometry = 'z=30'
      new_sideset_name = 'front'
    []

    [side_top]
       type = ParsedGenerateSideset
       input = side_front
       combinatorial_geometry = 'y=10'
       new_sideset_name = 'top'
     []

     [side_bottom]
        type = ParsedGenerateSideset
        input = side_top
        combinatorial_geometry = 'y=0'
        new_sideset_name = 'bottom'
      []

      [side_left]
         type = ParsedGenerateSideset
         input = side_bottom
         combinatorial_geometry = 'x=0'
         new_sideset_name = 'left'
       []

       [side_right]
          type = ParsedGenerateSideset
          input = side_left
          combinatorial_geometry = 'x = 60'
          new_sideset_name = 'right'
        []
[]

#############################################################################

[Modules]
[./FluidProperties]
[./the_simple_fluid]
type = SimpleFluidProperties
[../]
[../]
[]
[PorousFlowFullySaturated]
porepressure = porepressure
coupling_type = Hydro
gravity = '0 0 -9.81'
fp = the_simple_fluid
time_unit = years
add_darcy_aux = false
[]

[Variables]
  [./porepressure]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./darcy_x]
    order = FIRST
    family = MONOMIAL
  [../]
  [./darcy_y]
    order = FIRST
    family = MONOMIAL
  [../]
  [./darcy_z]
    order = FIRST
    family = MONOMIAL
  [../]
  [./perm_map]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./poro_map]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./min]
  order = CONSTANT
  family = MONOMIAL
[../]
[]

[AuxKernels]

  [./min]
   type = ElementLengthAuxsqrtBEH
   variable = min
   method = min
   execute_on = initial
 [../]
  [./darcy_x]
    type = PorousFlowDarcyVelocityComponent
    variable = darcy_x
    gravity = '0 0 -9.81'
    component = x
    [../]
    [./darcy_y]
      type = PorousFlowDarcyVelocityComponent
      variable = darcy_y
      gravity = '0 0 -9.81'
      component = y
    [../]
    [./darcy_z]
      type = PorousFlowDarcyVelocityComponent
      variable = darcy_z
      gravity = '0 0 -9.81'
      component = z
    [../]
  [./perm_map]
  type = PorousFlowPropertyAux
  variable = 'perm_map'
  property = permeability
  [../]
  [./poro_map]
  type = PorousFlowPropertyAux
  variable = 'poro_map'
  property = porosity
  [../]
  []

#######################################################################
[BCs]

     [./water_grad_left]
       type = FunctionDirichletBC
       variable = porepressure
       boundary = 'left'
       function = ini_pp
     [../]


          [./water_grad_right]
            type = PorousFlowPiecewiseLinearSinkBEH
          variable = porepressure
          boundary = 'right'
          pt_vals = '0 1e9'
          multipliers = '0 1e9'
          flux_function =  10
          v = min
          use_mobility = true
          fluid_phase = 0
          [../]

[]
#############################################################

[UserObjects]
  [./poro]
    type = PropertyReadFile
    prop_file_name = 'porosity_test.csv'
    nprop = 1
    read_type = element
  [../]

  [./perm_read]
    type = PropertyReadFile
    prop_file_name = 'perm_hydro_fault_V2.csv'
    nprop = 8
    read_type = element
  [../]

  [./damage_read_1]
    type = PropertyReadFile
    prop_file_name = 'damage_01.csv'
    nprop = 1
    read_type = element
  [../]
[]

##################################################
[Functions]
  [./ini_pp]
      type = ParsedFunction
      value =  9.81*1000*(25-z)
  [../]
[]
##################################################################
[Materials]
  [./porosity]
    type = PorousFlowPorosityConstCSV
    read_prop_user_object = poro
  [../]
  [./permeability_1]
    type = PorousFlowPermeabilityConstCSV
    read_prop_user_object = perm_read
    damage = damage
  [../]
  [./damage1]
    type = GenericConstantArrayBEH
    prop_name = damage
    read_prop_user_object = damage_read_1
  [../]
[]

[Problem]
  material_dependency_check = false
[]
########################################################################

########################################################################
# to producing pwp v node for coupling
########################################################################
[VectorPostprocessors]
  [./nodal_pwp]
    type = NodalValueSampler
    variable = 'porepressure'
    sort_by = id
    execute_on = FINAL
    unique_node_execute = true
  [../]
[]
##########################################################################
[Preconditioning]
[./SMP]
type = SMP
full = true
[../]
[]
[Executioner]
type = Steady
solve_type = NEWTON
petsc_options = '-snes_converged_reason'
petsc_options_iname = '-pc_type -pc_hypre_type'
petsc_options_value = 'hypre    boomeramg'
nl_rel_tol = 5e-5
# nl_abs_tol = 5000.0
l_tol = 1e-06
l_max_its = 150
nl_max_its = 50
[]
[Debug]
show_var_residual_norms = true
[]
##########################################################################
[Outputs]

  [./pwp_list]
    type = CSV
    file_base = /home/moose/projects/moose_beh/TEST_MOOSE_model/hydro_fault_test_pwp_list
    execute_on = FINAL
  []

  [./out]
    type = Exodus
    file_base = /home/moose/projects/moose_beh/TEST_MOOSE_model/hydro_fault_test
    execute_on = FINAL
  [../]
[]


##########################################################################
##### MPC
##########################################################################
[MPCbe2]
  csv_file = HTTRANS2020_G01HTv2_TETL_W2FF_MOOSEMPCS_MPCs.csv
header = true
[]
