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
The [tools](tools) folder contains the common tools that are commonly used.

The [time conversions](tools/time_conversions.py) module includes a function to convert Julian dates to the corresponding Martian sol number and corresponding solar longitude.
This function can also convert Julian dates from Tudat (in seconds since J2000) to Julian dates understood by the MCD (in days since J2023).

The [plot utilities](tools/plot_utilities.py) module contains functions for plotting.

## Requirements
First of, the MCD data files are required. By default, the code wants them to be in `/mnt/c/MCD/data/`.
This can however be changed in the Python files of the [MCD](MCD) folder.

The following Python libraries are required:
```
numpy
matplotlib
scipy
```

In addition, it is required to install `TudatPy`. This is the TU Delft Astrodynamics Toolbox used to run the astrodynamic simulations.
Installing this environment can be done by following the steps described [on this page](https://github.com/tudat-team/tudat-bundle#readme).
Please note that, in my case, the Windows Subsystem for Linux (v2) has been used.

Also, before compiling Tudat on your machine, [this line](https://github.com/tudat-team/tudatpy/blob/4169c827eaa16bf4b6cc9b8626d29f54c6724a76/tudatpy/kernel/expose_simulation/expose_environment_setup/expose_atmosphere_setup.cpp#L86) shall be changed to:
```
m.def("custom_constant_temperature_detailed",
```