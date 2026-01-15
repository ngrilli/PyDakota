# Wan Wan
# Nicolò Grilli
# University of Bristol
# 15 Gennaio 2026

# class for cyclic data and calibration using backstress and effective stress

import numpy as np
from experiment_comp import Experiment_Comp

class Cyclic(Experiment_Comp)

    def __init__(self, cycle_period=1.0, strain_rate=0.001):
        super().__init__()
        self.cycle_period = cycle_period # full cycle
        self.strain_rate = strain_rate
        
    # calculate max and min stress in each cycle
    def max_stress(self):

    def triangular_strain(self, strain_rate):
        amplitude = 
        phase = (t % self.cycle_period) / self.cycle_period
    strain = np.zeros_like(phase)

    for i, p in enumerate(phase):
        if p < 0.25:
            strain[i] = 4 * amplitude * p
        elif p < 0.75:
            strain[i] = 2 * amplitude - 4 * amplitude * p
        else:
            strain[i] = -4 * amplitude + 4 * amplitude * p

    return strain

    def calc_residual(self): # override
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
        # calculate total number of cycles: assume x data is time
        self.number_of_cycles = sim_x_data.max()
        # calculate max stress (positive/negative)
        self.maxstress()
        return 1
