# SPARTA

## Commands

mpirun -np 12 spa_ < in.apollo
convert image_iso*ppm sim_iso.gif
convert image_front*ppm sim_front.gif

Binary STL to ASCII:
python2 stl_B2A.py ABE_sail.stl
STL to surface:
python2 stl2surf.py ASCII_ABE_sail.stl data.ABE_sail

## Conditions

The conditions in which the simulation has been made have been varied, as to gather drag values for different conditions, for interpolation later on.

First, the orbital altitudes have been varied between 95 km and 190 km.
This leads to the following parameters:

| Altitude [km] | Velocity [m/s] | Density [kg/m3] | Temperature [K] | Pressure [Pa] | Mixture [mol/mol]                                     |
|---------------|----------------|-----------------|-----------------|---------------|-------------------------------------------------------|
| 85            | 3494.17        | 7.1E-07         | 135             | 2.3E-02       | 90.5% CO2, 3.5% N2, 2.5% Ar, 1.5% CO, 1.5% O, 0.5% O2 |
| 115           | 3493.29        | 1.8E-08         | 115             | 3.7E-04       | 81% CO2, 4.5% N2, 3.5% Ar, 5% CO, 5.5% O, 0.5% O2     |
| 150           | 3483.82        | 1.6E-10         | 175             | 7.1E-06       | 42% CO2, 12.5% N2, 4.5% Ar, 15% CO, 25% O, 1% O2      |

The velocities have been taken from the orbital simulations made in [MCD/feasible_altitudes.py](../MCD/feasible_altitudes.py).
The mixture, temperatures, pressures, and densities have been obtained from the [online interface](http://www-mars.lmd.jussieu.fr/mcd_python) of the Mars Climate Database, taking time and position averages. values.

The script [comp_inputs.py](setup/comp_inputs.py) has then been used to compute the relevant inputs for the SPARTA simulation.
This leads to the inputs as in the table below, for the "Sail" type satellite, in a box of 0.5x0.5x0.7 meters.

| Altitude [km] | Velocity [m/s] | Density [#/m3] | Mixture fractions [-]                    | Minimum grid size (x, y, z) [m] | Timestep [s] | f_num [-] | Knudsen number |
|---------------|----------------|----------------|------------------------------------------|---------------------------------|--------------|-----------|----------------|
| 85            | 3494.17        | 1.003E+19      | 0.905, 0.035, 0.025, 0.015, 0.015, 0.005 | 2.403E+01, 1.717E+01, 1.717E+01 | 3.188E-06    | 7.673E+15 | 5.571E-01      |
| 115           | 3493.29        | 2.660E+17      | 0.81, 0.045, 0.035, 0.05, 0.055, 0.005   | 6.217E-01, 4.440E-01, 4.440E-01 | 1.233E-04    | 1.176E+19 | 2.154E+01      |
| 150           | 3483.82        | 2.983E+15      | 0.42, 0.125, 0.045, 0.15, 0.25, 0.01     | 5.994E-03, 4.281E-03, 4.281E-03 | 1.282E-02    | 1.472E+23 | 2.234E+03      |


85km: Drag = 1.26980e-01 N (use grid of 1.7E+01 1.2E+01 1.2E+01 to avoid SPARTA error)


## Results
