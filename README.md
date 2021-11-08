# Mars Very Low Orbit feasibility study
Project of the Short Training Program at the von Karman Institute.
Study of the feasibility of using air-breathing engines for satellite altitude maintenance in very low Mars orbit.

## Repo structure

### Global Reference Atmospheric Model for Mars (Mars GRAM 2010)
The [GRAM](GRAM) folder contains a single Python script, [call_GRAM.py](GRAM/call_GRAM.py), that can be used to query the atmospheric density of Mars as a function of altitude, latitude, and time. These density are extrapolated based on the data from Mars GRAM 2010, by NASA.

The variable `file_path` at the beginning of the Python file refers to the path of the `tpdmsy11.txt` data file that is distributed with Mars GRAM 2010. It can be requested to NASA by filling the form specified at [this page](https://software.nasa.gov/software/MFS-33158-1).

### Mars Climate Database (MCD)
The [MCD](MCD) folder contains the MCD Fortran interface that has been compiled in Python, as well as a class written to load this interface in parallel for each Martian month. This makes calls to the MCD thousands of times faster.

The 3Go of data that constitues the MCD have to be obtained by contacting its developer, as explained [on this page](http://www-mars.lmd.jussieu.fr/mars/access.html).
By default, the MCD data files should be located in `/mnt/c/MCD/data/`. However, this can be changed in the `self.dset` variable in [this module](MCD/parallel_mcd.py).

More details about this folder are given in its own [README](MCD/README.md).

### Stochastic PArallel Rarefied-gas Time-accurate Analyzer (SPARTA)
The [SPARTA](SPARTA) folder contains everything that has been used to compute the drag coefficient of the different satellite configurations at different altitudes.

It contains the [STL files](SPARTA/setup/STL) corresponding to each of the satellites configurations, with their solar panels deployed.

Its [tools folder](SPARTA/tools) contain two script used to convert the STL files to SPARTA surfaces.

The [comp_inputs.py](SPARTA/setup/comp_inputs.py) script can be run to automatically create the SPARTA input files used to configure the simulations.
This script also takes care of running the conversion from STLs to SPARTA surfaces.

[analyse_results.py](SPARTA/analyse_results.py) can be used to compute the drag from the SPARTA simulation result files.

Finally, the [ParaView folder](SPARTA/paraview) contains the configuration files to convert the raw data from SPARTA to ParaView.

Explanation on how to install SPARTA, as well as a much higher level of detail about its use, can be found in its own [README](SPARTA/README.md).

### Figures
The [figures folder](figures) simply contains all of the plots that have been generated to support the research. These are both in PDF format, or in HTML when they are interactive.

### Setup selection
The [setup selection](setup_selection) folder first contains the [integrators_propagators folder](setup_selection/integrators_propagators) with scripts made to select an integrator and propagator, and their respective settings.
A baseline is made, and an combination of integrator and propagator that results in a fast simulation but low error is selected.
This way, the propagation time of a satellite orbiting one year around Mars has been lowered by a factor of around 20s, at the cost of a reasonable deviation in propagated position.

Then, the [environments folder](setup_selection/environments) contains scripts to select the appropriate environment models, and test them.

Both sub-folders of the setup selection contains their own READMEs: one for the [environments](setup_selection/environments/README.md) and one for the [integratos/propagators](setup_selection/integrators_propagators/README.md).

### Thrust
The [thrust](thrust) folder contains two test scripts: [test_simplest.py](thrust/test_simplest.py) is used to test the implementation of thrust, and [test_hall_thrust.py](thrust/test_hall_thrust.py) is used to test the implementation of the electric thrust.

### Tools
The [tools](tools) folder contains the tools that are commonly used:

 * the [mission geometry](tools/mission_geometry.py) module contains a function to compute the shadowing fraction of Mars between the Satellite and the Sun, and a function to compute the power from the solar panels given the satellite orientation, its position, and the position of the Sun.

 * the [plot utilities](tools/plot_utilities.py) module contains functions to make, and potentially save, various type of plots.

 * the [std](tools/std.py) module contains a class that can be use to absorbd (and silence) all outputs from a specified set of code. This prevents C++ code from polluting the console (especially during optimisation runs).

 * the [time conversions](tools/time_conversions.py) module includes a function to convert Julian dates to the corresponding Martian sol number and corresponding solar longitude. This function can also convert Julian dates from Tudat (in seconds since J2000) to Julian dates understood by the MCD (in days since J2023).

### Utils
The [utils](utils) folder also contains tools, except they are now called utilities, and they are related to the simulation itself.

Most importantly, this folder contains all of the orbital propagation code, the satellite models, as well as the thrust models.
More information can be found on these utilities in [their README](utils/README.md).

## Python requirements

All of the Python modules required to run this repository can be installed by using the provided Conda environment. 

### Conda environment
To ease the installation, a [conda environment file](environment.yaml) has been created.
In the same folder as this file, the following command can then be used to install the conda environment with its required packages:
```
conda env create -f environment.yaml
```
This environment then contains most importantly TUDAT(Py) (TU Delft Astrodynamics Toolbox Python) and Pygmo (Python version of the Parallel Global Multiobjective Optimizer).

Before running any code, one must make sure that this environment is activated.
This can be done using:
```
conda activate tudat-pygmo-vki
```

If errors arise when running part of the code, it may be wise to force conda to use TUDAT(Py) version 0.5.22 and Pygmo version 2.16.1. This can be done by uncommenting the version numbers in the [conda environment file](environment.yaml).

### Install with cmake
Both Tudat(Py) and Pygmo can be installed by building their C++ libraries and their Python interface, using Pybind and CMake. This allows to edit their code, but requires a more tedious process to be followed, as explained below.

#### TU Delft Astrodynamics Toolbox
In addition, it is required to install `TudatPy`. This is the TU Delft Astrodynamics Toolbox used to run the astrodynamic simulations.
Please note that, in my case, the Windows Subsystem for Linux (v2) has been used.

First, tudat can be cloned into the folder you are in by doing the following:
```
git clone https://github.com/tudat-team/tudat-bundle
cd tudat-bundle
git submodule update --init --recursive
```

The conda environment can then be setup by using:
```
conda env create -f environment.yaml
conda activate tudat-bundle
```

Before compiling Tudat on your machine, [this line](https://github.com/tudat-team/tudatpy/blob/4169c827eaa16bf4b6cc9b8626d29f54c6724a76/tudatpy/kernel/expose_simulation/expose_environment_setup/expose_atmosphere_setup.cpp#L86) shall be changed to:
```
m.def("custom_constant_temperature_detailed",
```

The warnings on [these lines](https://github.com/tudat-team/tudat/blob/fa30c49dca7ee27630717efb8546802589a4c8b7/include/tudat/astro/propulsion/thrustGuidance.h#L185-L187) and [these lines](https://github.com/tudat-team/tudat/blob/fa30c49dca7ee27630717efb8546802589a4c8b7/src/astro/reference_frames/aerodynamicAngleCalculator.cpp#L383-L416) have also been commented out, to prevent them from polluting the console.

In `build.sh`, the vonfiguration and build steps should be replaced by the followings:
```
# configuration step
cmake -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
  -DCMAKE_CXX_STANDARD=14 \
  -DBoost_NO_BOOST_CMAKE=ON \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DTUDAT_BUILD_TESTS="${BUILD_TESTS}" \
  -DTUDAT_BUILD_WITH_NRLMSISE00=OFF \
  ..

# build step
cmake --build . -j8
```

Then, the library can be compiled using:
```
bash build.sh
```

Finally, the correct tudat installation folder can be added to the conda environment path by using the following:
```
echo "<tudat-bundle installation dir>/build/tudatpy" > ~/miniconda3/envs/tudat-bundle/lib/python3.8/site-packages/include_path.pth
```

#### Pygmo
The Pygmo optimisation toolbox from ESA shall be installed following the instructions [here](https://esa.github.io/pygmo2/install.html).

In essence, the following commands can be used when the `tudat-bundle` environment is active, to install the required Pygmo packages:
```
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install pygmo
```