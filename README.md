# Offline Marine Ecosystem Model (OMEM)

## Contents:
[1. Model Description](#description)
2. How to Use
3. Example of results: Effect of a tropical cyclone on the Kuroshio region’s primary productivty

<a name="description"/>
## 1. Model Description


The Offline Marine Ecosystem Model (OMEM) consists on a physical and marine ecosystem components forced with ocean currents, temperature and radiation from a parent model. The physical component includes an advection-diffusion scheme forced by zonal and meridional velocities. 
The offline model is based on the Community Earth System Model (CESM) code. The OMEM's marine ecosystem is borrowed - except for some changes- from the Marine Biogeochemistry Library (MARBL) of CESM. The offline advection-diffusion scheme is based on the Parallel Ocean Program (POP). 
The offline model also features inputs for surface deposition of dust, iron, nitrate and ammonia and interior sources of hydrothermal and sedimentary iron.


OMEM is fully written in Python. To look at the original MARBL fortran code visit https://github.com/marbl-ecosys/MARBL, and to obtain more information about POP go to https://github.com/ESCOMP/POP2-CESM. 



## Input data:

OMEM requires daily data for: zonal and meridional velocities, temperature and short wave radiation. The velocities and temperature are in 3D lat-lon-depth, and short wave radiation is 2D (lat-lon). 
If desired, 2D surfaces fluxes for nitrate, iron and dust can be input (also daily).
Interior fluxes of iron from deposition and hydrothermal vents can be added; these are only lat-lon-depth but do not vary on time.
Initial Conditions are borrowed from CESM ()
Open Boundary Conditions were borrowed from a global run of CESM2. But they can be adapted to any fields you may like. They need to be input as 3D arrays and the code calculates the north, south, east and west arrays.

|Field | Long Name | Input File Units | Model Units (modified in code) |
|------|-----------|------------------|--------------------------------|
|U (uvelg) | Zonal velocity | cm/s | cm/s|
|V (vvelg)| Meridional velocity | cm/s| cm/s|
|W (wcontg) | Vertical velocity| cm/s  | cm/s (W is calculated from the continuity eq.)|
|hordiff | horizontal diffusivity (set to constant value)  |  |  2e^7 cm/s^2|
|KVMIX (kvmix) | Vertical mixing coefficient  | 0.25 cm/s^2|

|Interior Forcing | Long Name | Input File Units | Conversion Factor | Model Units |
|-----------------|-----------|------------------|-------------------|-------------|
|qswg | Shortwave radiation | W/m^2 | 1 |W/m^2|
|tempg | Tempreature | C| 1 |  C|
|feflux_sedg  | Iron flux from sediment deposition  | micromol Fe/m^2/day | 1.1574e^-6 | nmolFe/cm^2/s|
|feflux_ventg | Iron flux from vents | micromol Fe/m^2/day | 1.1574e^-6 | nmolFe/cm^2/s|


Model units obtained from https://marbl.readthedocs.io/en/latest/usr-guide/GCM-interface/GCM_requirements/forcing_fields.html; Input File does not conatin units, but the conversion factor and file units are explained at https://www.cesm.ucar.edu/models/cesm1.0/pop2/doc/users/POPecosys_main.html and https://github.com/ESCOMP/POP2-CESM/blob/master/bld/namelist_files/namelist_defaults_pop.xml 

| Surface Forcing | Long Name | Input File Units | Conversion factor | Model Units |
|-----------------|-----------|------------------|-------------------|-------------|
|feflux_solg | Iron dust flux| kg/m^2/s|  6.2668e^4 (3.5% iron per weight)  |  nmol/cm^2/s | 
|dustg  | Dust deposition |   |  kg N/m^2/s | 7.1429e^6 |nmol/cm^2/s | 
|NOxg  | Nitrogen deposition |  kg N/m^2/s  |  7.1429e^6 | nmol/cm^2/s  | 


## 2. How to Use 

First, the data needs to be masked in a propper manner. U and V determine the masking of the data. 
1- The input data (except U and V) cannot have any masked elements (it has to be interpolated over land), so the first part of the code multiplies the input fields by the correct mask.
2- Then vertical velocity (W) is calculated using the continuity equation.
3- Two latitudes and longitues are lost due to the processes of masking and W calculation. Also, the deepst most vertical level is removed, since the input initial conditions did not contain this level (this can be easily changed).

### Modules Files:

- utilies.py has several usefull tools to  inteprolate fro monthy to daily data, and contains a function to calculate latitudinal weights used on the x-grid spacing.
- continuity.py calculates W from U and V.
- saver.py saves the output at the desired frequency.

- process_inputs.py processes and masks the necesary input data

- run.py calls process_input.py if necesary to process the data, then calls continity.py to calculate W; it also calls main.py to run the model and run_saver.py to save the output.

- main.py contains the main code that evolves the equations on time, it calls ecosys.py and adv_diff.py

- ecosystem.py contains the ecosystem equations
- ecosys_constants.py contains all the necesary constants used in the ecosys.py module

- adv_diff.py contains the offline advection-diffusion scheme.

-density.py calculates the density at level k, and the density adiabatically displaced to level k+1, this is used if convective adjustment is activated


## Model Outputs


|Model Index | Component | Long Name | Input File Units |
|------------|-----------|-----------|------------------|
|0 | PO4 | Dissolved Inorganic Phosphate |  mmol/m^3 |
|1 | NO3 | Dissolved Inorganic Nitrate |  mmol/m^3 |
|2 | SiO3 | Dissolved Inorganic Silicate |  mmol/m^3 |
|3 | NH4 | Dissolved Inorganic Ammonia |  mmol/m^3 |
|4 | Fe | Dissolved Inorganic Iron |  mmol/m^3 |
|5 | Lig | Iron Binding Ligand |  mmol/m^3 |
|6 | DOC | Dissolved Organic Carbon |  mmol/m^3 |
|7 | DON | Dissolved Organic Nitrogen |  mmol/m^3|
|8 | DOP | Dissolved Organic Phosphorus |  mmol/m^3|
|9 | DOPr | Refractory DOP |  mmol/m^3|
|10 | DONr | Refractory DON |  mmol/m^3|
|11 | DOCr | Refractory DOC |  mmol/m^3|
|12 | zooC | Zooplankton Carbon |  mmol/m^3|
|13 | spC | Small Phytoplankton Carbon |  mmol/m^3 |
|14 | spP | Small Phytoplankton Phosphorus|  mmol/m^3 |
|15 | spChl | Small Phytoplankton Chlorophyll|  mmol/m^3 |
|16 | spFe | Small Phytoplankton Iron |  mmol/m^3 |
|17 | spCaCO_3 |Small Phytoplankton Calcium Carbonate |  mmol/m^3 |
|18 | diatC | Diatoms Carbon |  mmol/m^3 |
|19 | diatChl |Diatoms Chlorophyll |  mmol/m^3|
|20 | diatSi | Diatoms Silicate |  mmol/m^3|
|21 | diatFe | Diatoms Iron |  mmol/m^3|
|22 | diatP | Diatoms Phosphorus |  mmol/m^3|
|23 | diazC | Diazotrophs  Carbon |  mmol/m^3|
|24 | diazChl | Diazotrophs  Chlorophyll  |  mmol/m^3|
|25 | diazFe | Diazotrophs  Iron |  mmol/m^3|
|26 | diazP | Diazotrophs   Phosphorus |  mmol/m^3|
|27 | PIC soft | Particulate Inorganic Carbon (soft matter) |  mmol/m^3|
|28 | PIC hard | Particulate Inorganic Carbon (hard matter)  |  mmol/m^3|
|29 | PSiO2  soft| Particulate Inorganic Silicon  (soft matter)  |  mmol/m^3|
|30 |  PSiO2  hard | Particulate Inorganic Silicon  (hard matter) |  mmol/m^3|
|31 | PDUST soft | Sinking Dust (soft) |  mmol/m^3|
|32 | PDUST hard | Sinking Dust (hard)  |  mmol/m^3|
|33 | POC soft | Particulate Organic Carbon (soft) |  mmol/m^3|
|34 | POC hard | Particulate Organic Carbon (hard)  |  mmol/m^3|
|35 | PFe soft | Particulate Iron (soft) |  mmol/m^3|
|36 | PFe hard | Particulate Iron (hard)  |  mmol/m^3|
|37 | POP soft | Particulate Organic Phosphorus (soft)  |  mmol/m^3|
|38 | POP hard | Particulate Organic Phosphorus (hard)  |  mmol/m^3|
|39 | QA dust deficit | ? |  mmol/m^3|

Note that  mmol/m^3 = \mu mol/L (micromol/L) = nmol/cm^3, the input initial conditions are in  \mu mol/L. The units were obtained from the file marbl-ecosys-MARBL-e0d512d/src/default\_settings.yaml, the code was downloaded from https://zenodo.org/record/2541008#.YDROgs9KhJE



# How to use it:
Copy this repository into your computer. Add all the input files that you need in the same directory.

## 3. Example: Effect of a tropical cyclone on the Kuroshio region’s primary productivty

I created this model for my PhD thesis, which can be found at [https://scholarspace.manoa.hawaii.edu/items/b8ae1b00-a65a-425d-9b21-2acb495068c3/full]. The initial and open boundary conditions, as well as the surface and interior fluxes, are obtained from a historical climatology (1950-1960) of the low resolution (1 degree) run of the Community Earth System Model 2 (CESM2). One year of daily velocities, temperature and short wave radiation are borrowed from the parent model CESM1.2.2 (Small et al., 2014), which runs under present-day fixed CO2 concentration of 367 ppm. The oceanic currents featured by CESM 1.2.2 are further described by Chu et al. (2020). The online simulations of CESM2 and CESM1.2.2 were run on the Aleph supercomputer from the Institute for Basic Science (IBS) Center for Climate Physics (ICCP) in Busan, Korea. The offline model is run on a personal computer, with a 2.3 GHz 8-Core Intel Core i9 processor, and it takes 4.5 days to run 1 simulation year, with a time step of 400 s. 



