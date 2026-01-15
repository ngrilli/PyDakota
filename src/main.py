# Nicolò Grilli
# University of Bristol
# 10 Aprile 2022

# optimization code for FEM simulations
# main function

import numpy as np
import pathlib
import os
from moose_simulation import MOOSE_Simulation
from experiment_comp import Experiment_Comp
from parameters_optimizer import Parameters_Optimizer

###################################
#   parameters that must be set   #

# name of the simulation input file
input_file_name = 'CalibrationSC.i'

# name of the template input file
template_file_name = 'CalibrationSC.template'

# simulation output file with x-y data
output_file_name = 'CalibrationSC_out.csv'

# file with experimental data
experiment_file = 'StrainStress.csv'

# name of the png figure file
figure_file_name = 'StrainStress.png'

# name of the experimental data on the x and y axes
# must correspond to the header of the column in the experiment_file
name_of_x_experimental_data = 'Strain'
name_of_y_experimental_data = 'Stress (MPa)'

# names of the corresponding simulated data
# must correspond to the header of the column in the output_file_name
name_of_x_simulated_data = 'time'
name_of_y_simulated_data = 'stress_on_load_surface'

# name of the simulation command
# that you would normally run from a terminal open
# in the folder with all the simulation files
name_of_simulation_cmd = '~/projects/c_pfor_am/c_pfor_am-opt -i CalibrationSC.i'

# it is possible to define prefactors for simulated
# and experimental data, for instance to convert simulation time
# into simulation strain
prefactor_x_experimental = 1.0
prefactor_y_experimental = 1.0
prefactor_x_simulated = 0.2230 / 60.0
prefactor_y_simulated = 1.0

# type of scipy solver
type_of_solver = 'Nelder-Mead'

# name of the parameters to be substituted in the template file
# and initial values and bounds for the optimization
# in the template file parameters for optimization appear as {name_of_parameter}
# and are substituted by a number into the simulation input file
# at each iteration of the algorithm
opt_variables = ['init_dislo','hard_rate']
init_variables = np.array([20.0,0.04])
bnds_variables = np.array([(0, None), (0, None)])

# calibration for cyclic data using backstress and effective stress
cyclic = False

# end parameters that must be set #
###################################

# create lists of data names
name_exp = [name_of_x_experimental_data,name_of_y_experimental_data]
name_sim = [name_of_x_simulated_data,name_of_y_simulated_data]

# create lists of prefactors
prefactors_exp_xy = np.array([prefactor_x_experimental,prefactor_y_experimental])
prefactors_sim_xy = np.array([prefactor_x_simulated,prefactor_y_simulated])

# create a MOOSE simulation
ms = MOOSE_Simulation(input_file_name,template_file_name,name_of_simulation_cmd,opt_variables)

# create an experiment
exp = Experiment_Comp(output_file_name,experiment_file,name_exp,name_sim,prefactors_exp_xy,prefactors_sim_xy,figure_file_name)

# create a parameters optimizer
po = Parameters_Optimizer(ms,exp,init_variables,bnds_variables,type_of_solver)

# launch optimization
final_result = po.find_optimal_parameters()

# print optimization results to terminal
print(final_result)














