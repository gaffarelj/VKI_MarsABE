# SPARTA

The Stochastic PArallel Rarefied-gas Time-accurate Analyzer (SPARTA) software library has been used to get the drag coefficient of the satellite given different atmospheric conditions. SPARTA is a Direct Simulation Monte Carlo (DSMC) simulator, suited for rarefied flow simulations.

## Installation
The SPARTA documentation can be accessed [here](https://sparta.github.io/).

Moreover, it can be installed and compiled by running the following commands:
```
git clone https://github.com/sparta/sparta
cd sparta
mkdir build
cd build
module load openmpi
cmake -LH ../cmake
make
```
The `spa_` compiled file will then be found in `build/src`.

This folder can be added to path using the following command:
```
echo "export PATH=\$PATH:/path/to/sparta/build/src" >> ~/.bash_profile
```

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
With ImageMagick installed, the following command can be used to transform them to an animated `.gif` file:
```
convert image*ppm movie.gif
```

Most CAD softwares such as CATIA generate a binary STL file. To convert it to an ASCII STL file, the following command can be used (with `stl_B2A.py` from the [tools folder](tools)):
```
python2 stl_B2A.py model.stl
```
Additionally, the `-rs` option can be added to rescale the binary STL from `mm` to `m` when translating it to ASCII:
```
python2 stl_B2A.py model.stl -rs
```

Finally, to convert an ASCII STL file to a surface that SPARTA understands, the following command can be used (with `stl2surf.py` from the [tools folder](tools)):
```
python2 stl2surf.py model_ascii.stl data.model
```

When all of the input files have been created, they can all be run by using the following command while in the `inputs` folder:
```
bash run_all.sh
```

The SPARTA data files describing the geometry can be converted to ParaView by using the following command (from the [paraview/surf](paraview/surf) folder):
```
pvpython ../../tools/surf2paraview.py ../../setup/data/data.CS_XXXX CS_XXXX
```

Similarly, the number of particles in the grid can be converted to ParaView by using the command below (from the [paraview/grid](paraview/grid) folder):
```
pvpython ../../tools/grid2paraview.py grid.CS_XXXX_YYkm CS_XXXX_YYkm -r ../../setup/results_sparta/CS_XXXX/npart_YYkm.*.gz
```
In the command above, grid.CS_XXXX_._YYkm is a file describing the grid geometry for satellite CS_XXXX at an altitude of YY km.

These two last commands can be run automatically by running `bash paraview_convert.sh` in the [paraview](paraview) folder.

## Satellite configurations

Different CubeSat configurations have been modeled, each with different solar panels configurations.

The naming scheme has been taken as CS ABCD:
 * CS: CubeSat
 * A: Solar panels above and below the main body
 * B: Extension of the solar panels in A (straight behind)
 * C: Solar panels on the main body
 * D: Extension of the solar panels in C (behind, at a 15deg angle from the centerline)

The following figure illustrates this naming scheme:

<img src="CS/naming.png" width="35%">

The number of solar panels (of size 300x100mm) has been counted for one side of the satellite. 

As an example, the [3U deployable solar array from EnduroSat](https://www.endurosat.com/cubesat-store/cubesat-solar-panels/3u-single-deployable-solar-array) could be used.
Including PCB, these have a mass of 0.27 kg and a thickness of 1.75 mm. At EOL, they have a power efficiency of at least 29%.

Finally, the area of the solar panels of each satellite has been computed across three different planes. This later simplifies solar power calculations.
These x, y, and z planes are defined as in the figure below:

<img src="CS/planes.png" width="35%">

| Satellite name | # of solar panels (x2) | Reference length [m] | Reference surface area [m2] | SA x-area [m2] | SA y-area [m2] | SA z-area [m2] | Illustration                           |
|----------------|------------------------|----------------------|-----------------------------|----------------|----------------|----------------|----------------------------------------|
| CS 0020        | 2                      | 0.3                  | 0.01                        | 0              | 0.042426       | 0.042426       | <img src="CS/CS_0020.png" width="85%"> |
| CS 1020        | 4                      | 0.341421             | 0.0104                      | 0              | 0.102426       | 0.042426       | <img src="CS/CS_1020.png" width="85%"> |
| CS 0021        | 4                      | 0.589778             | 0.041058285                 | 0.031058       | 0.083343       | 0.083343       | <img src="CS/CS_0021.png" width="85%"> |
| CS 2020        | 6                      | 0.541421             | 0.0108                      | 0              | 0.162426       | 0.042426       | <img src="CS/CS_2020.png" width="85%"> |
| CS 1021        | 6                      | 0.589778             | 0.041458285                 | 0.031058       | 0.143343       | 0.083343       | <img src="CS/CS_1021.png" width="85%"> |
| CS 3020        | 8                      | 0.741421             | 0.0112                      | 0              | 0.222426       | 0.042426       | <img src="CS/CS_3020.png" width="85%"> |
| CS 2021        | 8                      | 0.589778             | 0.041858285                 | 0.031058       | 0.203343       | 0.083343       | <img src="CS/CS_2021.png" width="85%"> |
| CS 2120        | 10                     | 0.6                  | 0.0108                      | 0              | 0.282426       | 0.042426       | <img src="CS/CS_2120.png" width="85%"> |
| CS 3021        | 10                     | 0.741421             | 0.042258285                 | 0.031058       | 0.263343       | 0.083343       | <img src="CS/CS_3021.png" width="85%"> |

## Conditions

The conditions in which the simulation has been made have been varied, as to gather different drag values, for interpolation later on.

First, the orbital altitudes have been varied between 85 km and 150 km.
This leads to the following parameters:

| Altitude [km] | Velocity [m/s] | Density [kg/m3]             | Dynamic pressure [Pa] | Temperature [K]     | Pressure [Pa]               | Mixture [mol/mol]                                     |
|---------------|----------------|-----------------------------|-----------------------|---------------------|-----------------------------|-------------------------------------------------------|
| 85            | 3510.90        | 1.27915E-06 (σ=5.54578E-07) | 7.88373E+00           | 140.264 (σ=10.3875) | 3.47635E-02 (σ=1.65223E-02) | 95.298% CO2, 1.909% N2, 1.960% Ar, 0.422% CO, 0.250% O, 0.162% O2 <br> (σ=1.199% CO2, 0.293% N2, 0.303% Ar, 0.376% CO, 0.232% O, 0.039% O2) |
| 115           | 3495.84        | 2.21588E-08 (σ=1.45383E-08) | 1.35400E-01           | 114.185 (σ=9.53312) | 4.84327E-04 (σ=3.18527E-04) | 87.489% CO2, 4.340% N2, 4.041% Ar, 2.341% CO, 1.409% O, 0.380% O2 <br> (σ=6.234% CO2, 1.968% N2, 1.308% Ar, 1.931% CO, 1.018% O, 0.180% O2) |
| 150           | 3475.51        | 1.94535e-10 (σ=8.33132e-11) | 1.17491E-03           | 172.404 (σ=23.5187) | 7.27188E-06 (σ=3.21932E-06) | 65.826% CO2, 11.923% N2, 4.910% Ar, 7.796% CO, 8.533% O, 1.012% O2  <br> (σ=13.260% CO2, 4.501% N2, 0.592% Ar, 3.837% CO, 4.956% O, 0.198% O2) |

These values have been gathered by running the `run_atmo_study` in [feasible_altitudes.py](../MCD/feasible_altitudes.py).
As can be seen in the code, they are atmospheric values that are averaged for a satellite that flies in a circular orbit inclined at 45deg for 2 Martian years.

## Input files
The script [comp_inputs.py](setup/comp_inputs.py) has then been used to generate the relevant input files for the SPARTA simulations.

These input files have been made for each satellite, for each of the satellite configurations.

These can be found in the [inputs](setup/inputs), and are of the format `in.*`.

It is worth noting that grid size has been capped to a maximum of (10, 10, 10), to avoid having the grid size too big compared to the geometry (satellite) size.

Additionally, the `f_num` parameter has been tuned to ensure that enough simulated particles are present. At h=85km, `f_num` has been increased by a factor of 100. It has been decreased by a factor of 1E3 at h=115km, and decreased by a factor of 1E7 for h=150km.

Note that all SPARTA simulations assume that the satellite accommodation coefficient is based on the adsorption of atomic oxygen (as in [this paper](https://doi.org/10.2514/1.49330)).

## Results

Running SPARTA for the different altitudes, the force in each direction has been saved in `.gz` files in the [results folder](results).

Then, [analyse_results.py](analyse_results.py) has been used to compute the drag force from all of the simulations.

At each altitude, the dynamic pressure (in Pa) has been computed as follows:

<img src="https://render.githubusercontent.com/render/math?math=q%20=%20\frac{1}{2}%20\cdot%20\rho%20\cdot%20V^2">

For each satellite configuration, the reference surface `S` has been taken as the frontal area of the CubeSat.
Values for these reference surface area can be found in the table of [this section](#satellite-configurations).

The drag coefficients have then been computed at each altitude by using the following equation:

<img src="https://render.githubusercontent.com/render/math?math=C_D%20=%20\frac{D}{q%20\cdot%20S}=100%20\cdot%20\frac{D}{q}">

Later on, when using these drag coefficients, the same reference surface must be used.

This leads to the drag coefficients of the table below, with the Knudsen numbers included as well:

| Satellite name | Altitude [km] | Knudsen number [-] | Drag [N]    | Drag coefficient [-] |
|----------------|---------------|--------------------|-------------|----------------------|
| CS 0020        | 85            | 1.842E-01          | 4.845E-01   | 2.91898              |
| CS 0020        | 115           | 3.555E-01          | 2.450E-01   | 2.92598              |
| CS 0020        | 150           | 1.491E+03          | 5.269E-05   | 2.86794              |
| CS 1020        | 85            | 9.370E-02          | 7.845E-01   | 4.54414              |
| CS 1020        | 115           | 1.808E-01          | 4.545E-01   | 5.21927              |
| CS 1020        | 150           | 7.583E+02          | 1.326E-04   | 6.93935              |
| CS 0021        | 85            | 1.619E-01          | 6.884E-01   | 1.01005              |
| CS 0021        | 115           | 3.124E-01          | 3.529E-01   | 1.02643              |
| CS 0021        | 150           | 1.310E+03          | 6.439E-05   | 0.85362              |
| CS 2020        | 85            | 9.370E-02          | 9.501E-01   | 5.29961              |
| CS 2020        | 115           | 1.808E-01          | 5.365E-01   | 5.93222              |
| CS 2020        | 150           | 7.583E+02          | 1.426E-04   | 7.18779              |
| CS 1021        | 85            | 1.021E-01          | 8.913E-01   | 1.29512              |
| CS 1021        | 115           | 1.970E-01          | 4.594E-01   | 1.32335              |
| CS 1021        | 150           | 8.260E+02          | 7.674E-05   | 1.00754              |
| CS 3020        | 85            | 9.370E-02          | 1.169E+00   | 6.28985              |
| CS 3020        | 115           | 1.808E-01          | 6.493E-01   | 6.92407              |
| CS 3020        | 150           | 7.583E+02          | 1.548E-04   | 7.52191              |
| CS 2021        | 85            | 9.210E-02          | 1.095E+00   | 1.57640              |
| CS 2021        | 115           | 1.778E-01          | 5.884E-01   | 1.67884              |
| CS 2021        | 150           | 7.453E+02          | 9.557E-05   | 1.24272              |
| CS 2120        | 85            | 7.453E-02          | 1.081E+00   | 6.02714              |
| CS 2120        | 115           | 1.438E-01          | 5.604E-01   | 6.19718              |
| CS 2120        | 150           | 6.032E+02          | 8.909E-05   | 4.49018              |
| CS 3021        | 85            | 7.453E-02          | 1.368E+00   | 1.95032              |
| CS 3021        | 115           | 1.438E-01          | 7.566E-01   | 2.13827              |
| CS 3021        | 150           | 6.032E+02          | 1.670E-04   | 2.15149              |