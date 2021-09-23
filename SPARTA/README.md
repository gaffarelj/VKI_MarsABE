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
| 95            | 3303.14        | 4.9E-07         | 125             | 8.5E-03       | 80% CO2, 9.5% N2, 4.5% Ar, 2.75% CO, 2.75% O, 0.5% O2 |
| 140           | 3312.32        | 1.2E-09         | 150             | 1.5E-05       | 67% CO2, 7% N2, 5% Ar, 7.5% CO, 12.5% O, 1% O2        |
| 190           | 3402.66        | 6.4E-12         | 155             | 3.5E-07       | 26% CO2, 13% N2, 2% Ar, 15% CO, 43% O, 1% O2          |

The velocities have been taken from the orbital simulations made in [MCD/feasible_altitudes.py](../MCD/feasible_altitudes.py).
The mixture, temperatures, pressures, and densities have been obtained from the [online interface](http://www-mars.lmd.jussieu.fr/mcd_python) of the Mars Climate Database, taking time and position averages. values.

The script [comp_inputs.py](setup/comp_inputs.py) has then been used to compute the relevant inputs for the SPARTA simulation.
This leads to the inputs as in the table below, for the "Sail" type satellite, in a box of 0.5x0.5x0.7 meters.

| Altitude [km] | Velocity [m/s] | Density [#/m3] | Mixture fractions [-]                    | Minimum grid size (x, y, z) [m] | Timestep [s] | f_num [-] | Knudsen number |
|---------------|----------------|----------------|------------------------------------------|---------------------------------|--------------|-----------|----------------|
| 95            | 3303.14        | 7.190E18       | 0.8, 0.095, 0.045, 0.0275, 0.0275, 0.005 | 1.734E01, 1.238E01, 1.238E01    | 4.676E-06    | 1.466E16  | 7.724E-01      |
| 140           | 3312.32        | 1.908E16       | 0.67, 0.07, 0.05, 0.075, 0.125, 0.01     | 4.222E-02, 3.016E-02, 3.016E-02 | 1.915E-03    | 2.694E21  | 3.172E02       |
| 190           | 3402.66        | 1.413E14       | 0.26, 0.13, 0.02, 0.15, 0.43, 0.01       | 2.259E-04, 1.614E-04, 1.614E-04 | 3.483E-01    | 1.302E26  | 2.927E04       |

## Results
