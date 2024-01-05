# readMe

In this folder, the entire piece of Lennard-Jones potential code has been scripted. Entirely in Python and optimized in [Jupyter lab](https://jupyter.org/install). There are sections that have been divided into functions for better readability and further optimization. 

**lj_parameters.py** contains all the parameters and the initial variables that are essential to calculate a variety of parameters. It also contains those files in which the calculated results would be written for further processing. 

**lj_functions.py** contains all the functions that are required to generate particles, their motions, calculate temperature, energies, and a lot more.

**lj_mainloop.py** contains the primary loop that would allow the calculations and compilation. This script would calculate instantaneous temperature, pressure, kinetic, potential, and total energies. The loop would also show real-time progress and with adjustments to the internal parameters, this code can be customized to calculate even more parameters.   
