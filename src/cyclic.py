# Wan Wan
# Nicolò Grilli
# University of Bristol
# 15 Gennaio 2026

# class for cyclic data and calibration using backstress and effective stress

import numpy as np
from experiment_comp import Experiment_Comp

class Cyclic(Experiment_Comp)

    def __init__(self, cycle_period=1.0, strain_rate=0.001, elastic_strain_range=30):
        super().__init__()
        self.cycle_period = cycle_period # full cycle
        self.strain_rate = strain_rate
        self.elastic_strain_range = elastic_strain_range
        
    # calculate max and min stress in each cycle
    def max_stress(self):
        self.max_stress_array = np.zeros(shape=(self.number_of_cycles))
        self.min_stress_array = np.zeros(shape=(self.number_of_cycles))
        
        return 1

    # calculate strain with triangular time dependence
    def triangular_strain(self, strain_rate):
        self.strain = np.zeros(shape=(self.number_of_cycles))
        amplitude_cumulative = strain_rate * self.cycle_period
        for i in len(self.sim_x_data):
            phase = (self.sim_x_data[i] % self.cycle_period) / self.cycle_period
            if (phase < 0.25):
                self.strain[i] = amplitude_cumulative * phase
            elif (phase < 0.75):
                self.strain[i] = (amplitude_cumulative / 2.0) - amplitude_cumulative * phase
            else:
                self.strain[i] = amplitude_cumulative * phase - amplitude_cumulative
        return 1

    # Least squares fit: y = m x + b. Returns m youngs modulus, b intercept
    def linear_fit(self, x, y):
        A = np.vstack([x, np.ones_like(x)]).T
        m, b = np.linalg.lstsq(A, y, rcond=None)[0]
        return float(m), float(b)

    # calculate the difference between experiment and simulation: override
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
        # store exp and sim data into data structures
        self.exp_x_data = exp_x_data
        self.exp_y_data = exp_y_data
        self.sim_x_data = sim_x_data
        self.sim_y_data = sim_y_data
        self.exp_y_interp_onto_sim = exp_y_interp_onto_sim
        # calculate total number of cycles: assume x data is time
        self.number_of_cycles = np.floor(sim_x_data.max() / self.cycle_period)
        # calculate max stress (positive/negative)
        self.max_stress()
        return 1
