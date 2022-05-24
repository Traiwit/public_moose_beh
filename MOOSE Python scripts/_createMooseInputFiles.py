#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 13:40:58 2022

@author: moose
"""

#input parameters

import csv

surf_list = list(csv.reader(open('GNKT_9steps_surfacelist.csv')))   #read in the csv file

damage_file_list = list(csv.reader(open('GNKT_9steps_damagelist.csv')))

frame_num = len(surf_list)

MESH_file = 'GNKT2020066_G01HTv1_D01_Q01_v7_MOOSEMESH.inp'
PERM_file = 'GNKT2020066_G01HTv1_D01_Q01_v5_MATPROPS_PERM.csv' # perm file name
PORO_file = 'GNKT2020066_G01HTv1_D01_Q01_v5_MATPROPS_PORO.csv' # poro file name
MPC_file = 'GNKT2020066_G01HTv1_D01_Q01_v7_MOOSEMPCS_MPCs.csv' #MPC file name
project_tag = 'GNKT_short_9steps' #project name
app_name = 'moose_beh-opt'
directory = '/home/moose/projects' #project directory

# parameters
y_pos = 5700
y_neg = -4449.01
x_pos = 5110
x_neg = -4135.03
z_zero_level = 120
flux_fx = 10


## solver related
mpi = 4   # number of MPI cores
nl_rel_tol = 5e-14
nl_abs_tol = 5e3
l_tol = 1e-6
l_max_its = 150
nl_max_its = 50


#%% Surface_files




for block in range(1, frame_num):
            block_pre = block - 1
            surf_name = (surf_list[block][0]) # a =
            damage_file = (damage_file_list[block][0])

            with open(project_tag+'_'+surf_name+'_surface.i', 'w') as f:
                f.write( '[BCs] \n'
                           '[./water_grad_left]\n'
                            'type = FunctionDirichletBC\n'
                             'variable = porepressure\n'
                             'boundary = \'left\'\n'
                             'function = ini_pp\n'
                           '[../]\n'
                           '[./water_grad_right]\n'
                             'type = FunctionDirichletBC\n'
                             'variable = porepressure\n'
                             'boundary = \'right\'\n'
                             'function = ini_pp\n'
                           '[../]\n'
                           '[./water_grad_top]\n'
                             'type = FunctionDirichletBC\n'
                             'variable = porepressure\n'
                             'boundary = \'top\'\n'
                             'function = ini_pp\n'
                           '[../]\n'
                           '[./water_grad_bottom]\n'
                             'type = FunctionDirichletBC\n'
                             'variable = porepressure\n'
                             'boundary = \'bottom\'\n'
                             'function = ini_pp\n'
                           '[../]\n'
                           )

            with open(project_tag+'_'+surf_name+'_surface.i', 'a') as f:

                f.write( '[./drain_'+'%s] \n' %surf_name +
                 ' type = PorousFlowPiecewiseLinearSinkBEH \n'
                 ' variable = porepressure \n'
                 ' boundary = \'%s\'\n'  %surf_name +
                 ' pt_vals = \'0 1e15\' \n'
                 ' multipliers = \'0 1e15\' \n'
                 ' flux_function = %s \n' %flux_fx +
                 ' v = min\n'
                 ' PT_shift = 0 \n'
                 ' use_mobility = true\n'
                 ' fluid_phase = 0\n'
                 '[../]\n'
                       )
            with open(project_tag+'_'+surf_name+'_surface.i', 'a') as f:
                   f.write( '[] \n'
                           )




            with open(project_tag+'_'+surf_name+'_surface.i', 'a') as f:
                   f.write('[UserObjects] \n'
                              '[./perm_read] \n'
                                'type = PropertyReadFile \n'
                                'prop_file_name = \'%s\' \n' %PERM_file +
                                'nprop = 8 #kxx kyy kzz kxy kxz kyz Aw Kmax \n'
                                'read_type = element \n'
                              '[../] \n'

                              '[./poro_read] \n'
                                'type = PropertyReadFile \n'
                                'prop_file_name = \'%s\' \n' %PORO_file +
                                'nprop = 1 \n'
                                'read_type = element \n'
                              '[../] \n'

                                '[./damage_read]\n'
                                'type = PropertyReadFile\n'
                                'prop_file_name = \'%s\' \n' %damage_file +
                                'nprop = 1\n'
                                ' read_type = element\n'
                                '[../]\n'

                            '[] \n'
                            )

            with open(project_tag+'_'+surf_name+'_surface.i', 'a') as f:
                   f.write('[Outputs] \n'
                           'file_base = '+'%s/%s/%s_%s \n' %(directory,project_tag,project_tag,surf_name) +
                           'perf_graph = true \n'
                            '[./run_'+'%s] \n' %surf_name +
                              'type = Exodus \n'
                              'execute_on = \'final\' \n'
                            '[../] \n'
                            '[] \n'
                            )



#%% Base_files

for block in range(1, frame_num):
            block_pre = block - 1
            surf_name = (surf_list[block][0])
            pre = (surf_list[block_pre][0])

            with open(project_tag+'_'+surf_name+'_base.i', 'w') as f:

                   f.write( '[Mesh]\n'
                    'allow_renumbering = false\n'
                    '[mesh1]\n'
                    'type = FileMeshGenerator\n'
                    'allow_renumbering = false\n'
                    'file= %s_%s.e  \n' %(project_tag, pre) +
                    '[] \n'
                    '[side_top]\n'
                    'type = ParsedGenerateSideset\n'
                    'input = mesh1\n'
                    'combinatorial_geometry = \'y<=%s\'  \n' %y_neg +
                    'new_sideset_name = \'top\'  \n'
                    '[] \n'
                    '[side_bottom]\n'
                    'type = ParsedGenerateSideset\n'
                    'input = side_top\n'
                    'combinatorial_geometry = \'y>=%s\'  \n'  %y_pos +
                    'new_sideset_name = \'bottom\'  \n'
                    '[] \n'
                    '[side_left]\n'
                    'type = ParsedGenerateSideset\n'
                    'input = side_bottom\n'
                    'combinatorial_geometry = \'x>=%s\'  \n' %x_pos +
                    'new_sideset_name = \'left\'  \n'
                    '[]\n'
                    '[side_right]\n'
                    'type = ParsedGenerateSideset\n'
                    'input = side_left \n'
                    'combinatorial_geometry = \'x<=%s\'  \n' %x_neg +
                    'new_sideset_name = \'right\'  \n'
                    '[]\n'
                    '[]\n'
                    '[UserObjects]\n'
                    '[initial_mesh]\n'
                    'type = SolutionUserObject\n'
                    'execute_on = INITIAL\n'
                    'mesh= %s_%s.e  \n' %(project_tag, pre) +
                    'timestep = LATEST\n'
                    'system_variables = porepressure\n'
                    '[]\n'
                    '[]\n' )

            with open(project_tag+'_'+surf_name+'_base.i', 'a') as f:

               f.write( '[Functions] \n'
                  '[./ini_pp] \n'
                    'type = ParsedFunction\n'
                    'value  = 9.81*1000*(%s-z)\n' %z_zero_level +
                  '[../]\n'

                  '[equl_pp]\n'
                    'type = SolutionFunction\n'
                    'from_variable = porepressure\n'
                    'solution = initial_mesh\n'
                  '[]\n'
                '[]\n'
                #
                '[ICs]\n'
                  '[./porepressure]\n'
                    'type = FunctionIC\n'
                    'variable = porepressure\n'
                    'function = equl_pp\n'
                  '[../]\n'
                '[]\n'

                '[Preconditioning]\n'
                  '[./SMP]\n'
                    'type = SMP\n'
                    'full = true\n'
                  '[../]\n'
                '[]\n'
                '[Executioner]\n'
                  'type = Steady\n'
                  'solve_type = NEWTON\n'
                  'petsc_options = \'-snes_converged_reason\'  \n'
                  'petsc_options_iname = \'-pc_type -pc_hypre_type\'  \n'
                  'petsc_options_value = \'hypre    boomeramg\'  \n'

                  'nl_rel_tol = %s\n' %nl_rel_tol +
                  'nl_abs_tol = %s\n' %nl_abs_tol +
                  'l_tol = %s\n'      %l_tol +
                  'l_max_its = %s\n'  %l_max_its +
                  'nl_max_its = %s\n' %nl_max_its +
                  '[]\n'
                  '[Debug]\n'
                     'show_var_residual_norms = true \n'
                  '[]\n'

                '[Modules]\n'
                  '[./FluidProperties]\n'
                    '[./the_simple_fluid]\n'
                      'type = SimpleFluidProperties\n'
                    '[../]\n'
                  '[../]\n'
                '[]\n'

                '[PorousFlowFullySaturated]\n'
                  'porepressure = porepressure\n'
                  'coupling_type = Hydro\n'
                  'gravity = \'0 0 -9.81\'  \n'
                  'fp = the_simple_fluid\n'
                  'time_unit = years\n'
                  'add_darcy_aux = false\n'
                '[]\n'

                '[GlobalParams]\n'
                  'PorousFlowDictator = dictator\n'
                '[]\n'

                '[Problem]\n'
                  'material_dependency_check = false\n'
                '[]\n'

                '[Variables]\n'
                  '[./porepressure]\n'
                    'order = FIRST\n'
                    'family = LAGRANGE\n'
                  '[../]\n'
                '[]\n'

                '[AuxVariables]\n'
                  '[./min]\n'
                    'order = CONSTANT\n'
                    'family = MONOMIAL\n'
                  '[../]\n'
                   '[./damage]\n'
                    'order = CONSTANT\n'
                    'family = MONOMIAL\n'
                  '[../]\n'
                  '[./perm_map]\n'
                    'order = CONSTANT\n'
                    'family = MONOMIAL\n'
                  '[../]\n'
                  '[/poro_map]\n'
                    'order = CONSTANT\n'
                    'family = MONOMIAL\n'
                  '[../]\n'
                  '[./darcy_x]\n'
                    'order = CONSTANT\n'
                    'family = MONOMIAL\n'
                  '[../]\n'
                  '[./darcy_y]\n'
                    'order = CONSTANT\n'
                    'family = MONOMIAL\n'
                  '[../]\n'
                  '[./darcy_z]\n'
                    'order = CONSTANT\n'
                    'family = MONOMIAL\n'
                  '[../]\n'
                '[]\n'

                '[AuxKernels]\n'
                    '[./min]\n'
                      'type = ElementLengthAuxsqrtBEH\n'
                      'variable = min\n'
                      'method = min\n'
                      'execute_on = INITIAL\n'
                    '[../]\n'
                    '[./damage]\n'
                      'type = MaterialRealAux\n'
                      'variable = damage\n'
                      'property = damage\n'
                      'execute_on = INITIAL\n'
                    '[../]\n'
                    '[./perm_map]\n'
                      'type = PorousFlowPropertyAux\n'
                      'variable = \'perm_map\'  \n'
                      'property = permeability\n'
                    '[../]\n'
                    '[./poro_map]\n'
                      'type = PorousFlowPropertyAux\n'
                      'variable = \'poro_map\'  \n'
                      'property = porosity \n'
                    '[../]\n'
                    '[./darcy_x]\n'
                      'type = PorousFlowDarcyVelocityComponentYearToSec\n'
                      'variable = darcy_x\n'
                      'gravity = \'0 0 -9.81\'  \n'
                      'component = x\n'
                    '[../]\n'
                    '[./darcy_y]\n'
                      'type = PorousFlowDarcyVelocityComponentYearToSec\n'
                      'variable = darcy_y\n'
                      'gravity = \'0 0 -9.81\'  \n'
                      'component = y\n'
                    '[../]\n'
                    '[./darcy_z]\n'
                      'type = PorousFlowDarcyVelocityComponentYearToSec\n'
                      'variable = darcy_z\n'
                      'gravity = \'0 0 -9.81\'  \n'
                      'component = z\n'
                    '[../]\n'
                  '[]\n'


                    '[Materials]\n'
                      '[./porosity]\n'
                        'type = PorousFlowPorosityConstCSV\n'
                        'read_prop_user_object = poro_read\n'
                      '[../]\n'

                      '[./permeability_1]\n'
                        'type = PorousFlowPermeabilityConstCSV\n'
                       ' read_prop_user_object = perm_read\n'
                       ' damage = damage\n'
                     ' [../]\n'

                           '[./damage1]\n'
                           'type = GenericConstantArrayBEH\n'
                           'prop_name = damage\n'
                           'read_prop_user_object = damage_read\n'
                           ' [../]\n'
                    '[] \n'

                    '[MPCbe2]\n'
                    'csv_file = %s \n' %MPC_file +
                    'header = true\n'
                    '[]\n'
                    )


               #%% producing Equal model

               damage_0_file = damage_file_list[0][0];
               equal_name = surf_list[0][0];

               with open(project_tag+'_'+equal_name+'_base.i', 'w') as f:
                    f.write('[Mesh]\n'
                    'allow_renumbering = false\n'
                    '[mesh1]\n'
                    'type = FileMeshGenerator\n'
                    'allow_renumbering = false\n'
                    'file= %s  \n' %(MESH_file) +
                    '[] \n'
                    '[side_top]\n'
                    'type = ParsedGenerateSideset\n'
                    'input = mesh1\n'
                    'combinatorial_geometry = \'y<=%s\'  \n' %y_neg +
                    'new_sideset_name = \'top\'  \n'
                    '[] \n'
                    '[side_bottom]\n'
                    'type = ParsedGenerateSideset\n'
                    'input = side_top\n'
                    'combinatorial_geometry = \'y>=%s\'  \n'  %y_pos +
                    'new_sideset_name = \'bottom\'  \n'
                    '[] \n'
                    '[side_left]\n'
                    'type = ParsedGenerateSideset\n'
                    'input = side_bottom\n'
                    'combinatorial_geometry = \'x>=%s\'  \n' %x_pos +
                    'new_sideset_name = \'left\'  \n'
                    '[]\n'
                    '[side_right]\n'
                    'type = ParsedGenerateSideset\n'
                    'input = side_left \n'
                    'combinatorial_geometry = \'x<=%s\'  \n' %x_neg +
                    'new_sideset_name = \'right\'  \n'
                    '[]\n'
                    '[]\n'

                    '[Functions] \n'
                  '[./ini_pp] \n'
                    'type = ParsedFunction\n'
                    'value  = 9.81*1000*(%s-z)\n' %z_zero_level +
                  '[../]\n'

                '[]\n'
                #
                '[ICs]\n'
                  '[./porepressure]\n'
                    'type = FunctionIC\n'
                    'variable = porepressure\n'
                    'function = ini_pp\n'
                  '[../]\n'
                '[]\n'

                '[Preconditioning]\n'
                  '[./SMP]\n'
                    'type = SMP\n'
                    'full = true\n'
                  '[../]\n'
                '[]\n'
                '[Executioner]\n'
                  'type = Steady\n'
                  'solve_type = NEWTON\n'
                  'petsc_options = \'-snes_converged_reason\'  \n'
                  'petsc_options_iname = \'-pc_type -pc_hypre_type\'  \n'
                  'petsc_options_value = \'hypre    boomeramg\'  \n'

                  'nl_rel_tol = %s\n' %nl_rel_tol +
                  'nl_abs_tol = %s\n' %nl_abs_tol +
                  'l_tol = %s\n'      %l_tol +
                  'l_max_its = %s\n'  %l_max_its +
                  'nl_max_its = %s\n' %nl_max_its +
                  '[]\n'
                  '[Debug]\n'
                     'show_var_residual_norms = true \n'
                  '[]\n'

                '[Modules]\n'
                  '[./FluidProperties]\n'
                    '[./the_simple_fluid]\n'
                      'type = SimpleFluidProperties\n'
                    '[../]\n'
                  '[../]\n'
                '[]\n'

                '[PorousFlowFullySaturated]\n'
                  'porepressure = porepressure\n'
                  'coupling_type = Hydro\n'
                  'gravity = \'0 0 -9.81\'  \n'
                  'fp = the_simple_fluid\n'
                  'time_unit = years\n'
                  'add_darcy_aux = false\n'
                '[]\n'

                '[GlobalParams]\n'
                  'PorousFlowDictator = dictator\n'
                '[]\n'

                '[Problem]\n'
                  'material_dependency_check = false\n'
                '[]\n'

                '[Variables]\n'
                  '[./porepressure]\n'
                    'order = FIRST\n'
                    'family = LAGRANGE\n'
                  '[../]\n'
                '[]\n'

                '[AuxVariables]\n'
                  '[./min]\n'
                    'order = CONSTANT\n'
                    'family = MONOMIAL\n'
                  '[../]\n'
                   '[./damage]\n'
                    'order = CONSTANT\n'
                    'family = MONOMIAL\n'
                  '[../]\n'
                  '[./perm_map]\n'
                    'order = CONSTANT\n'
                    'family = MONOMIAL\n'
                  '[../]\n'
                  '[/poro_map]\n'
                    'order = CONSTANT\n'
                    'family = MONOMIAL\n'
                  '[../]\n'
                  '[./darcy_x]\n'
                    'order = CONSTANT\n'
                    'family = MONOMIAL\n'
                  '[../]\n'
                  '[./darcy_y]\n'
                    'order = CONSTANT\n'
                    'family = MONOMIAL\n'
                  '[../]\n'
                  '[./darcy_z]\n'
                    'order = CONSTANT\n'
                    'family = MONOMIAL\n'
                  '[../]\n'
                '[]\n'

                '[AuxKernels]\n'
                    '[./min]\n'
                      'type = ElementLengthAuxsqrtBEH\n'
                      'variable = min\n'
                      'method = min\n'
                      'execute_on = INITIAL\n'
                    '[../]\n'
                    '[./damage]\n'
                      'type = MaterialRealAux\n'
                      'variable = damage\n'
                      'property = damage\n'
                      'execute_on = INITIAL\n'
                    '[../]\n'
                    '[./perm_map]\n'
                      'type = PorousFlowPropertyAux\n'
                      'variable = \'perm_map\'  \n'
                      'property = permeability\n'
                    '[../]\n'
                    '[./poro_map]\n'
                      'type = PorousFlowPropertyAux\n'
                      'variable = \'poro_map\'  \n'
                      'property = porosity \n'
                    '[../]\n'
                    '[./darcy_x]\n'
                      'type = PorousFlowDarcyVelocityComponentYearToSec\n'
                      'variable = darcy_x\n'
                      'gravity = \'0 0 -9.81\'  \n'
                      'component = x\n'
                    '[../]\n'
                    '[./darcy_y]\n'
                      'type = PorousFlowDarcyVelocityComponentYearToSec\n'
                      'variable = darcy_y\n'
                      'gravity = \'0 0 -9.81\'  \n'
                      'component = y\n'
                    '[../]\n'
                    '[./darcy_z]\n'
                      'type = PorousFlowDarcyVelocityComponentYearToSec\n'
                      'variable = darcy_z\n'
                      'gravity = \'0 0 -9.81\'  \n'
                      'component = z\n'
                    '[../]\n'
                  '[]\n'


                    '[Materials]\n'
                      '[./porosity]\n'
                        'type = PorousFlowPorosityConstCSV\n'
                        'read_prop_user_object = poro_read\n'
                      '[../]\n'

                      '[./permeability_1]\n'
                        'type = PorousFlowPermeabilityConstCSV\n'
                       ' read_prop_user_object = perm_read\n'
                       ' damage = damage\n'
                     ' [../]\n'

                           '[./damage1]\n'
                           'type = GenericConstantArrayBEH\n'
                           'prop_name = damage\n'
                           'read_prop_user_object = damage_read\n'
                           ' [../]\n'

                    '[] \n'

                      '[MPCbe2]\n'
                    'csv_file = %s \n' %MPC_file +
                    'header = true \n'
                    '[]\n'
                                )

                    with open(project_tag+'_'+equal_name+'_surface.i', 'w') as f:
                        f.write( '[BCs] \n'
                           '[./water_grad_left]\n'
                            'type = FunctionDirichletBC\n'
                             'variable = porepressure\n'
                             'boundary = \'left\'\n'
                             'function = ini_pp\n'
                           '[../]\n'
                           '[./water_grad_right]\n'
                             'type = FunctionDirichletBC\n'
                             'variable = porepressure\n'
                             'boundary = \'right\'\n'
                             'function = ini_pp\n'
                           '[../]\n'
                           '[./water_grad_top]\n'
                             'type = FunctionDirichletBC\n'
                             'variable = porepressure\n'
                             'boundary = \'top\'\n'
                             'function = ini_pp\n'
                           '[../]\n'
                           '[./water_grad_bottom]\n'
                             'type = FunctionDirichletBC\n'
                             'variable = porepressure\n'
                             'boundary = \'bottom\'\n'
                             'function = ini_pp\n'
                           '[../]\n'
                           '[]\n'

                           '[UserObjects] \n'
                              '[./perm_read] \n'
                                'type = PropertyReadFile \n'
                                'prop_file_name = \'%s\' \n' %PERM_file +
                                'nprop = 8 #kxx kyy kzz kxy kxz kyz Aw Kmax \n'
                                'read_type = element \n'
                              '[../] \n'

                              '[./poro_read] \n'
                                'type = PropertyReadFile \n'
                                'prop_file_name = \'%s\' \n' %PORO_file +
                                'nprop = 1 \n'
                                'read_type = element \n'
                              '[../] \n'

                                '[./damage_read]\n'
                                'type = PropertyReadFile\n'
                                'prop_file_name = \'%s\' \n' %damage_0_file +
                                'nprop = 1\n'
                                ' read_type = element\n'
                                '[../]\n'

                            '[] \n'

                            '[Outputs] \n'
                           'file_base = '+'%s/%s/%s_%s \n' %(directory,project_tag,project_tag,equal_name) +
                           'perf_graph = true \n'
                            '[./run_'+'%s] \n' %equal_name +
                              'type = Exodus \n'
                              'execute_on = \'final\' \n'
                            '[../] \n'
                            '[] \n'
                           )






               #%%



with open('run_command.sh', 'w') as f:

                f.write('mpiexec -n %s ./%s -i ~/projects/%s/%s_%s_base.i ~/projects/%s/%s_%s_surface.i --n-threads=4 --t \n' %(mpi,app_name,project_tag,project_tag,equal_name,project_tag,project_tag,equal_name) +
                            'echo done \n'
                            )

                for block in range(1, frame_num):
                    surf_name = (surf_list[block][0])

                    f.write(

                            'mpiexec -n %s ./%s -i ~/projects/%s/%s_%s_base.i ~/projects/%s/%s_%s_surface.i --n-threads=4 --t \n' %(mpi,app_name,project_tag,project_tag,surf_name,project_tag,project_tag,surf_name) +
                            'echo done \n'

                            )
