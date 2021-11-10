# Mars Climate Database interface

This folder contains different files and modules that can be used to call the Mars Climate Database (MCD).
The full version (v5.3) of this database has been obtained from [the MCD website](http://www-mars.lmd.jussieu.fr/mars/access.html).
Note that the MCD data files are not included in this Github repository, because they are above the file limit set by Github, and because they should not be openely shared.

## Fortran to Python interface
The `F2PY` Python module has been used to compile the Fortran interface into a Python module, [fmcd.so](fmcd.so).
Please not that, on a different machine, this module may need to be recompiled. Please refer to the MCD documentation to do this.

An example of call to the MCD can be found in [example_direct_call.py](example_direct_call.py).

A code that times a series of MCD calls can be found in [time_mcd_call.py](time_mcd_call.py).
This shows that loading the file for a Martian month takes around 3.5 seconds, and any successive call to the MCD for the same month takes around 0.025 ms.
However, using a time in a different month results in different files being loaded, costing 3.5 seconds again.

## Parallel loading of the MCD modules
To get around the waiting time when calls to the MCD are made for distinct Martian months, the `fmcd.so` interface can be loaded not once but 13 times.
This way, the call to the MCD can be made using the module that always has the files loaded for the same Martian month.
Such Martian month is defined every 30 deg of solar longitude.

Unfortunately, this lead to the need to compile the [fmcd.so](fmcd.so) interface 13 times individually, once for each Martian month (plus an extra).
This is because, if the same module is loaded under different names in Python, the same object is actually loaded.
Also, renaming the [fmcd.so](fmcd.so) does not suffice, as the name needs to be defined before the Fortran interface is compiled.
This explains files [fmcd_1.so](fmcd_1.so) to [fmcd_13.so](fmcd_13.so).

A script that explores and times this duplication of the modules can be seen in [test_parallel_mcd.py](test_parallel_mcd.py).

Finally, a Python class can be obtained from [parallel_mcd.py](parallel_mcd.py).
This class can be used to load the modules in parallel and call the appropriate one depending on the Martian month.
An example script that uses this class can be accessed in [test_parallel_mcd.py](test_parallel_mcd.py).

### Density and wind functions
The [parallel_mcd.py](parallel_mcd.py) class contains two functions that are used during the propagation. 

The `density` method returns the atmospheric density in kg/m3 at the given inputs:
 * `h`: altitude above the Martian geoid, in meters
 * `lon`: Martian longitude, in radians
 * `lat`: Martian latitude, in radians
 * `time`
   * with `time_is_JD=True` and `JD_Tudat=True`: Julian date in seconds since J2000
   * with `time_is_JD=True` and `JD_Tudat=False`: Julian date in days since J2023
   * with `time_is_JD=False`: array with the Solar longitude in degree as the first elements and the time of the day in hour as the second element

The `wind` method returns the wind in 3D for the same inputs. The three directions are as follows:
 * `x`: North+
 * `y`: East+
 * `z`: Point to Mars centre, **positive down**

## Feasible altitudes investigation
The [feasible_altitudes.py](feasible_altitudes.py) script investigates two distinct things.

First, by setting `run_alt_study=True`, a plot of the initial satellite altitude vs total orbital lifetime is generated. This is done using a satellite of 5kg, with a Cd of 3, and a reference surface area of 0.015 m2.

Secondly, by setting `run_atmo_study=True`, the average velocity, density, temperature, pressure, and mixture ratio of the atmosphere over time and position are shown.
This is done by propagating a satellite around Mars for 1300 days, taking only Mars as a point mass as an acceleration, and in an orbit with a 45deg inclination and an altitude of 85km, 115km, or 150km. This way, sampling the atmospheric values at the satellite position every 250s should give a good estimate of the average atmospheric conditions at a given altitude.

Please note that both methods take quite intensive CPU time. In addition, the second method also takes considerable memory space.