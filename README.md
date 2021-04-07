# OMEM
The Offline Marine Ecosystem Model is based on the Community Earth System Model (CESM) code. The OMEM's marine ecosystem is borrowed - except for some changes- from the Marine Biogeochemistry Library (MARBL). The offline advection-diffusion scheme is based on the Parallel Ocean Program (POP).

OMEM is fully written in Python. To look at the original MARBL fortran code visit https://github.com/marbl-ecosys/MARBL, and to obtain more information about POP go to https://github.com/ESCOMP/POP2-CESM. 


## Input data:

OMEM requires daily data for: zonal and meridional velocities, temperature and short wave radiation. The velocities and temperature are in 3D lat-lon-depth, and short wave radiation is 2D (lat-lon). 
If desired, 2D surfaces fluxes for nitrate, iron and dust can be input (also daily).
Interior fluxes of iron from deposition and hydrothermal vents can be added; these are only lat-lon-depth but do not vary on time.
Initial Conditions are borrowed from CESM ()
Open Boundary Conditions were borrowed from a global run of CESM2. But they can be adapted to any fields you may like. They need to be input as 3D arrays and the code calculates the north, south, east and west arrays.


## How to run:
First, the data needs to be masked in a propper manner. U and V determine the masking of the data. 
1- The input data (except U and V) cannot have any masked elements (it has to be interpolated over land), so the first part of the code multiplies the input fields by the correct mask.
2- Then vertical velocity (W) is calculated using the continuity equation.
3- Two latitudes and longitues are lost due to the processes of masking and W calculation. Also, the deepst most vertical level is removed, since the input initial conditions did not contain this level (this can be easily changed).

### Modules Files:

- utilies.py has several usefull tools to  inteprolate fro monthy to daily data, and contains a function to calculate latitudinal weights used on the x-grid spacing.
- continuity.py calculates W from U and V.
- run_saver.py saves the output at the desired frequency.

- process_inputs.py processes and masks the necesary input data

- run.py calls process_input.py if necesary to process the data, then calls continity.py to calculate W; it also calls main.py to run the model and run_saver.py to save the output.

- main.py contains the main code that evolves the equations on time, it calls ecosys.py and adv_diff.py

- ecosys.py contains the ecosystem equations
- ecosys_constants.py contains all the necesary constants used in the ecosys.py module

- adv_diff.py contains the offline advection-diffusion scheme.



## Model Outputs







