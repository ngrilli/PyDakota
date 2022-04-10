Nicolo Grilli
Michael Salvini
University of Bristol
10 Aprile 2022

Code to optimize x-y data variables obtained
from an external solver against experimental data

Installation: git clone the repository

Usage:

Put the simulation files in an arbitrary folder,
it must contain a template file from which the input or parameter file will be generated
and all other necessary files for your simulation

Set properly the parameters in /src/main.py

The template file will contain strings like {name_of_my_parameter}
that will be substituted by numbers during the optimization iterations.

An arbitrary number of these parameters can be specified in the template file
and the names must correspond to the ones in the corresponding string in main.py
order does not matter

Open a terminal in the folder with the simulation files and run:

python ~/path/to/src/main.py

