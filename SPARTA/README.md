# SPARTA

## Commands

mpirun -np 12 spa_ < in.apollo
convert image_iso*ppm sim_iso.gif
convert image_front*ppm sim_front.gif

Binary STL to ASCII:
python2 stl_B2A.py ABE_sail.stl
STL to surface:
python2 stl2surf.py ASCII_ABE_sail_2074triangles.stl data.ABE_sail

## Conditions

The conditions in which the simulation has been made have been varied, as to gather drag values for different conditions, for interpolation later on.

First, the orbital altitudes have been varied between 95 km and 190 km.
This leads to the following parameters:

| Altitude [km] | Velocity [m/s] | Density [kg/m3] | Mixture [mol/mol]                                     | Density [#/m3] |
|---------------|----------------|-----------------|-------------------------------------------------------|----------------|
| 95            | 3303.14        | 4.9E-07         | 80% CO2, 9.5% N2, 4.5% Ar, 2.75% CO, 2.75% O, 0.5% O2 |                |
| 140           | 3312.32        | 1.2E-09         | 67% CO2, 7% N2, 5% Ar, 7.5% CO, 12.5% O, 1% O2        |                |
| 190           | 3402.66        | 6.4E-12         | 26% CO2, 13% N2, 2% Ar, 15% CO, 43% O, 1% O2          |                |

The velocities have been taken from the orbital simulations made in [MCD/feasible_altitudes.py](../MCD/feasible_altitudes.py).
The mixture and densities have been obtained from the [online interface](http://www-mars.lmd.jussieu.fr/mcd_python) of the Mars Climate Database, taking time and position averages. values.

## Results

### Test
 * "Sail" satellite configuration
 * flow velocity of 3500 m/s
 * 300s, step of 0.0001s
 * 10% N, 10% O, 2% O2, 78% CO2
 * nrho of 1
 * 63M particles
 * 36M collisions (1.6B checked)
 * 15.4B particles move
 * Drag = 7.49525e-17 N