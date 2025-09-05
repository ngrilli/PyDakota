# Nicolò Grilli
# Wan Wan Mohammad
# University of Bristol
# 10 Aprile 2022

# parameter optimization procedure

from scipy.optimize import minimize
import numpy as np
from moose_simulation import MOOSE_Simulation
from experiment_comp import Experiment_Comp

class Parameters_Optimizer:

    def __init__(self,simulation,experiment,init_variables,bnds_variables,type_of_solver):
        self.simulation = simulation # this is a MOOSE_Simulation object
        self.experiment = experiment # this is a Experiment_Comp object
        self.init_variables = init_variables
        self.bnds_variables = bnds_variables
        self.type_of_solver = type_of_solver
        
    # find optimal parameters
    def find_optimal_parameters(self):
        risultato = minimize(self.calc_residual,self.init_variables,method=self.type_of_solver,options={'disp': True}) # add bounds
        return risultato
        
    # calculate residual
    # this function is used by scipy minimize
    # x is the array of parameters to optimize
    # that will change at each iteration
    def calc_residual(self,x):
        self.simulation.create_input_file(x)
        self.simulation.launch_simulation()
        temp_residual = self.experiment.calc_residual()
        self.experiment.plot_exp_sim()
        return temp_residual
    
    
    





	

