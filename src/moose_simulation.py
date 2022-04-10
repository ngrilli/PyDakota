# Nicolò Grilli
# University of Bristol
# 10 Aprile 2022

# a MOOSE simulation

import numpy as np
import subprocess

class MOOSE_Simulation:

    def __init__(self,input_file_name,template_file_name,name_of_simulation_cmd,opt_variables):
        self.input_file_name = input_file_name
        self.template_file_name = template_file_name
        self.name_of_simulation_cmd = name_of_simulation_cmd
        self.opt_variables = opt_variables
        
    def launch_simulation(self):
        proc = subprocess.Popen(self.name_of_simulation_cmd,shell=True)
        # waits for the process to finish
        # no need for further checks that the simulation is still running
        proc.wait()
        return 1
        
    # generate input file from template file
    # x is the array with the variables at the current iteration
    # of the optimization process
    def create_input_file(self,x):
        # Read the template file
        with open(self.template_file_name,'r') as file:
            filedata = file.read()
        # Replace the target strings
        for i in range(len(self.opt_variables)):
            temp_string_to_replace = '{' + self.opt_variables[i] + '}'
            filedata = filedata.replace(temp_string_to_replace,str(x[i]))
        # Write the modified input file
        with open(self.input_file_name,'w') as file:
            file.write(filedata)	
        return 1

