# Mars Very Low Orbit feasibility study
Project of the Short Training Program at the von Karman Institute.
Study of the feasibility of using air-breathing engines for satellite altitude maintenance in very low Mars orbit.

## Repo structure

### Mars Climate Database
The [MCD](MCD) folder contains the Mars Climat Database Fortran interface that has been compiled in Python, as well as a class written to load this interface in parallel for each Martian month.
This makes calls to the MCD thousands of times faster.

### Setup selection
The [setup selection](setup_selection) folder contains scripts made to select an integrator and propagator, and their supporting functions.
A baseline is made, and an combination of integrator and propagator that results in a fast simulation but low error is selected.

This way, the propagation time of a satellite orbiting one year around Mars has been lowered from 105 seconds to 10 seconds, at the cost of a maximal error of 600 m.

*To do*: add scripts to select appropriate environmental accelerations, such as a more detailled gravitational model of Mars.

### Tools
The [tools](tools) folder contains the tools that are commonly used.

The [time conversions](tools/time_conversions.py) module includes a function to convert Julian dates to the corresponding Martian sol number and corresponding solar longitude.
This function can also convert Julian dates from Tudat (in seconds since J2000) to Julian dates understood by the MCD (in days since J2023).

The [plot utilities](tools/plot_utilities.py) module contains functions for plotting.

The [mission geometry](tools/mission_geometry.py) module contains a function to compute the shadowing fraction of Mars between the Satellite and the Sun, and a function to compute the power from the solar panels given the satellite orientation, its position, and the position of the Sun.

### SPARTA
The [SPARTA](SPARTA) folder contains everything that has been used to compute the drag coefficient of the different satellite configurations.

It contains the [STL files](SPARTA/setup/STL) corresponding to each of the satellites configurations, with their solar panels deployed.

Its [tools folder](SPARTA/tools) contain two script used to convert the STL files to SPARTA surfaces.

The [comp_inputs.py](SPARTA/setup/comp_inputs.py) script can be run to automatically create the SPARTA input files used to configure the simulations.
This script also takes care of running the conversion from STLs to SPARTA surfaces.

Finally, [analyse_results.py](SPARTA/analyse_results.py) can be used to compute the drag from the SPARTA simulation result files.

### Thrust
The [thrust](thrust) folder contains a first [simple thrust class](thrust/simple_thrust.py) that is used to test the implementation of thrust and of the solar panel power.
This test class has been run using the [test_simplest.py](thrust/test_simplest.py) script.

## Requirements
This section lists the required Python packages, libraries, and software that shall be installed.

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

### Mars Climate Database
First of, the MCD data files are required. By default, the code wants them to be in `/mnt/c/MCD/data/`.
This can however be changed in the Python files of the [MCD](MCD) folder.

<!--
### TU Delft Astrodynamics Toolbox
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

### Pygmo
The Pygmo optimisation toolbox from ESA shall be installed following the instructions [here](https://esa.github.io/pygmo2/install.html).

In essence, the following commands can be used when the `tudat-bundle` environment is active, to install the required Pygmo packages:
```
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install pygmo
```
!-->

### SPARTA
The SPARTA library is required to run the simulations to obtain the drag coefficient of the different satellites.

More explanation on how to install it is given in the [SPARTA folder](SPARTA).