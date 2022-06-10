[Mesh]
allow_renumbering = false
[mesh1]
type = FileMeshGenerator
allow_renumbering = false
file= GNKT2020066_G01HTv1_D01_Q01_v7_MOOSEMESH.inp
[]
[side_top]
type = ParsedGenerateSideset
input = mesh1
combinatorial_geometry = 'y<=-4449.01'
new_sideset_name = 'top'
[]
[side_bottom]
type = ParsedGenerateSideset
input = side_top
combinatorial_geometry = 'y>=5700'
new_sideset_name = 'bottom'
[]
[side_left]
type = ParsedGenerateSideset
input = side_bottom
combinatorial_geometry = 'x>=5110'
new_sideset_name = 'left'
[]
[side_right]
type = ParsedGenerateSideset
input = side_left
combinatorial_geometry = 'x<=-4135.03'
new_sideset_name = 'right'
[]
[]
[Functions]
[./ini_pp]
type = ParsedFunction
value  = 9.81*1000*(120-z)
[../]
[]
[ICs]
[./porepressure]
type = FunctionIC
variable = porepressure
function = ini_pp
[../]
[]
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
nl_rel_tol = 1e-13
nl_abs_tol = 10000.0
l_tol = 1e-06
l_max_its = 150
nl_max_its = 75
[]
[Debug]
show_var_residual_norms = true
show_actions = true
[]
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
[GlobalParams]
PorousFlowDictator = dictator
[]
[Problem]
material_dependency_check = false
[]
[Variables]
[./porepressure]
order = FIRST
family = LAGRANGE
[../]
[]
[AuxVariables]
[./min]
order = CONSTANT
family = MONOMIAL
[../]
[./damage]
order = CONSTANT
family = MONOMIAL
[../]
[./perm_map]
order = CONSTANT
family = MONOMIAL
[../]
[/poro_map]
order = CONSTANT
family = MONOMIAL
[../]
[./darcy_x]
order = CONSTANT
family = MONOMIAL
[../]
[./darcy_y]
order = CONSTANT
family = MONOMIAL
[../]
[./darcy_z]
order = CONSTANT
family = MONOMIAL
[../]
[]
[AuxKernels]
[./min]
type = ElementLengthAuxsqrtBEH
variable = min
method = min
execute_on = INITIAL
[../]
[./damage]
type = MaterialRealAux
variable = damage
property = damage
execute_on = INITIAL
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
[./darcy_x]
type = PorousFlowDarcyVelocityComponentYearToSec
variable = darcy_x
gravity = '0 0 -9.81'
component = x
[../]
[./darcy_y]
type = PorousFlowDarcyVelocityComponentYearToSec
variable = darcy_y
gravity = '0 0 -9.81'
component = y
[../]
[./darcy_z]
type = PorousFlowDarcyVelocityComponentYearToSec
variable = darcy_z
gravity = '0 0 -9.81'
component = z
[../]
[]
[Materials]
[./porosity]
type = PorousFlowPorosityConstCSV
read_prop_user_object = poro_read
[../]
[./permeability_1]
type = PorousFlowPermeabilityConstCSV
 read_prop_user_object = perm_read
 damage = damage
 [../]
[./damage1]
type = GenericConstantArrayBEH
prop_name = damage
read_prop_user_object = damage_read
 [../]
[]
[MPCbe]
csv_file = GNKT2020066_G01HTv1_D01_Q01_v7_MOOSEMPCS_MPCs_dense.csv
header = true
[]
