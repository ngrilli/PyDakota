Wan Wan Mohammad
Nicolo Grilli
University of Bristol
14 March 2026

this extends optimisation framework to calibrate cyclic crystal plasticity 
against experimental low-cycle fatigue (LCF) data in cyclic.py

This code compares cycle-by-cycle stress evolution,backstress, and effective stress
between simulation and experiment. 

Required Input files:
1) experimental cyclic data.csv with test times, cycle number, strain, stress columns
2) experimental processed data.csv with gamma_cum, backstress, effective stress column 
3) simulation output.csv with time, macroscopic stress (avg_szz), strain(strain_eng) and cumulative slip (gamma_avg) columns

cyclic post-processing: reads experimental cyclic data -> read simulation output -> detect tensile peak stress per cycle
-> tracks cycle-by-cycle stress evolution -> compute backstress, effective stress, cum slip

simulation and calibrated output:
1) resolved_stresses_gamma_simulation.csv
2) 2 plots (max stress vs cycle number, backstress & effective stress vs cum slip)



