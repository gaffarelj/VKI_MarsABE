# SPARTA

The Stochastic PArallel Rarefied-gas Time-accurate Analyzer (SPARTA) software library has been used to get the drag coefficient of the satellite given different atmospheric conditions. SPARTA is a Direct Simulation Monte Carlo (DSMC) simulator, suited for rarefied flow simulations.

Its documentation can be accessed [here](https://sparta.github.io/), and how to compile and use it is explained in [this GitHub repository](https://github.com/sparta/sparta/blob/master/BUILD_CMAKE.md).

## Commands

Once SPARTA has been compiled, and the path to the `spa_` file added to path, the following commands can be used.

To run SPARTA on a given input script (for instance, satellite input script `in.sat`), the following command can be used:
```
spa_ < in.sat
```

SPARTA can be run on different CPU cores in parallel, if it was compiled with MPI.
Note that, doing so, the computer may run out of memory faster.
The following command can be used for this:
```
mpirun -np 12 spa_ < in.sat
```

In the input file, pictures can be generated at given time steps (with file extension `.ppm`).
With ImageMagick installed, the following command can be used to tranform them to an animated `.gif` file:
```
convert image*ppm movie.gif
```

Most CAD softwares such as CATIA generate a binary STL file. To convert it to an ASCII STL file, the following command can be used (with `stl_B2A.py` from the [tools folder](tools)):
```
python2 stl_B2A.py ABE_sail.stl
```

Finally, to convert an ASCII STL file to a surface that SPARTA understands, the following command can be used (with `stl2surf.py` from the [tools folder](tools)):
```
python2 stl2surf.py ASCII_ABE_sail.stl data.ABE_sail
```

## Conditions

The conditions in which the simulation has been made have been varied, as to gather different drag values, for interpolation later on.

First, the orbital altitudes have been varied between 85 km and 150 km.
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

## Input files
Three different SPARTA input files have been made in light of the different conditions.
These can be found in the [setup folder](setup), and are of the format `in.*`.

For the 85km altitude, the grid size has been manually replaced to `1.75E+01 1.25E+01 1.25E+01`. This is to avoid an error that occurs when the grid is too big compared to the geometry.

For the 115km altitude, the factor `f_num` has been lowered to `1.176E+16` to increase the number of simulated particles (and prevent it of being 0). Similarly, for 150km, `f_num` has been lowered to `1.472E+14`.

## Results

Running SPARTA for the different altitudes, the force in each direction has been saved in `.dat` files in the [results folder](results).

Then, [analyse_results.py](analyse_results.py) has been used to compute the drag force from all of the simulations.

At each altitude, the dynamic pressure (in Pa) has been computed as follows:

<img src="https://render.githubusercontent.com/render/math?math=q%20=%20\frac{1}{2}%20\cdot%20\rho%20\cdot%20V^2">

For the 'Sail' satellite type, the frontal area is computed as the central square plus the width of the solar panels multiplied by their thickness. This leads to a frontal area of:

<img src="https://render.githubusercontent.com/render/math?math=S%20=%20100%20\cdot%20100%20+%202%20\cdot%20191%20\cdot%202%20=%2010764%20mm^2%20=%200.010764%20m^2">

The drag coefficient has then been computed at each altitude by using the following equation:

<img src="https://render.githubusercontent.com/render/math?math=C_D%20=%20\frac{D}{q%20\cdot%20S}">

This leads to the drag coefficient of the table below:

| Altitude [km] | Velocity [m/s] | Density [kg/m3] | Dynamic pressure [Pa] | Drag [N]    | Drag coefficient [-] |
|---------------|----------------|-----------------|-----------------------|-------------|----------------------|
| 85            | 3494.17        | 7.1E-07         | 4.33E+00              | 1.27208E-01 | 2.7266               |
| 115           | 3493.29        | 1.8E-08         | 1.10E-01              | 2.74057E-03 | 2.3182               |
| 150           | 3483.82        | 1.6E-10         | 9.71E-04              | 2.68587E-05 | 2.5699               |
