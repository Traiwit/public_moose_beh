
[BCs]
[./water_grad_left]
type = FunctionDirichletBC
variable = porepressure
boundary = 'left'
function = ini_pp
[../]
[./water_grad_right]
type = FunctionDirichletBC
variable = porepressure
boundary = 'right'
function = ini_pp
[../]
[./water_grad_top]
type = FunctionDirichletBC
variable = porepressure
boundary = 'top'
function = ini_pp
[../]
[./water_grad_bottom]
type = FunctionDirichletBC
variable = porepressure
boundary = 'bottom'
function = ini_pp
[../]
[]
[UserObjects]
[./perm_read]
type = PropertyReadFile
prop_file_name = 'GNKT2020066_G01HTv1_D01_Q01_v5_MATPROPS_PERM.csv'
nprop = 8 #kxx kyy kzz kxy kxz kyz Aw Kmax
read_type = element
[../]
[./poro_read]
type = PropertyReadFile
prop_file_name = 'GNKT2020066_G01HTv1_D01_Q01_v5_MATPROPS_PORO.csv'
nprop = 1
read_type = element
[../]
[./damage_read]
type = PropertyReadFile
prop_file_name = 'Damage_0000.csv'
nprop = 1
 read_type = element
[../]
[]
[VectorPostprocessors]
[./nodal_pwp]
type = NodalValueSampler
variable = 'porepressure'
sort_by = id
execute_on = FINAL
unique_node_execute = true
[../]
[]
[Outputs]
perf_graph = true
[./run_DRAINSRF_F009_UG5]
type = Exodus
execute_on = 'final'
file_base = /home/moose/projects/GNKT_short_9steps/GNKT_short_9steps_EQUAL
[../]
[./pwp_list]
type = CSV
file_base = /home/moose/projects/GNKT_short_9steps/GNKT_short_9steps_EQUAL_pwp_list
execute_on = FINAL
[]
[]
