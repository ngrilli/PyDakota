# Nicolò Grilli
# Michael Salvini
# University of Bristol
# 10 Aprile 2022

# an experiment and how it compares with simulation

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class Experiment_Comp:

    def __init__(self,output_file_name,experiment_file,name_exp,name_sim,prefactors_exp_xy,prefactors_sim_xy,figure_file_name):
        self.output_file_name = output_file_name
        self.experiment_file = experiment_file
        self.name_exp = name_exp
        self.name_sim = name_sim
        self.prefactors_exp_xy = prefactors_exp_xy
        self.prefactors_sim_xy = prefactors_sim_xy
        self.figure_file_name = figure_file_name
        
    # calculate the difference between experiment and simulation
    def calc_residual(self):
        # read experimental file
        exp_data = pd.read_csv(self.experiment_file)
        # name of postprocessor file in MOOSE format
        sim_data = pd.read_csv(self.output_file_name)
        # create arrays with x and y data
        exp_x_data = self.prefactors_exp_xy[0] * exp_data[self.name_exp[0]].to_numpy()
        exp_y_data = self.prefactors_exp_xy[1] * exp_data[self.name_exp[1]].to_numpy()
        sim_x_data = self.prefactors_sim_xy[0] * sim_data[self.name_sim[0]].to_numpy()
        sim_y_data = self.prefactors_sim_xy[1] * sim_data[self.name_sim[1]].to_numpy()
        # interpolate experimental y over simulated x array
        # which is usually longer
        exp_y_interp_onto_sim = np.interp(sim_x_data, exp_x_data, exp_y_data)
        # store exp and sim data into class attributes
        # to plot them later
        self.exp_x_data = exp_x_data
        self.exp_y_data = exp_y_data
        self.sim_x_data = sim_x_data
        self.sim_y_data = sim_y_data
        return np.linalg.norm(exp_y_interp_onto_sim-sim_y_data)
        
    # plot and save comparison between experiment and simulation
    def plot_exp_sim(self):
        fig, ax = plt.subplots()
        ax.plot(self.sim_x_data,self.sim_y_data,linewidth=3,label='simulation')
        ax.plot(self.exp_x_data,self.exp_y_data,linewidth=3,label='experiment')
        ax.set_xlabel(self.name_exp[0],fontsize=14)
        ax.set_ylabel(self.name_exp[1],fontsize=14)
        ax.tick_params(axis="both",labelsize=14)
        fig.tight_layout()
        ax.legend(fontsize=14)
        fig.savefig(self.figure_file_name,dpi=200)
        return 1



	

