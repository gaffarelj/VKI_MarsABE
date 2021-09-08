# Mars Climate Database interface

This folder contains different files and modules that can be used to call the Mars Climate Database (MCD).
The full version (v5.3) of this database has been obtained from [the MCD website](http://www-mars.lmd.jussieu.fr/mars/access.html).
Note that the MCD data files are not included in this Github repository, because they are above the file limit set by Github, and because they should not be openely shared.

## Fortran to Python interface
The `F2PY` Python module has been used to compile the Fortran interface into a Python module, [fmcd.so](fmcd.so).
Please not that, on a different machine, this module may need to be recompiled. Please refer to the MCD documentation to do this.
Also, the MCD Fortran to Python interface has been compiled using Linux subsystem for Windows, as compilation did not work on Windows directly.
Because of this, any file calling the interface shall be run from a Linux terminal.
The Linux subsystem for Windows can be activated following the steps described [here](https://docs.microsoft.com/en-us/windows/wsl/install-win10).
The Linux terminal can be opened by pressing the `Shift` key and right-clicking in a folder on Windows, then clicking on `Open Linux shell here`.

An example of call to the MCD can be found in [example_direct_call.py](example_direct_call.py).

A code that times a series of MCD calls can be found in [time_mcd_call.py](time_mcd_call.py).
This shows that loading the file for a Martian month takes around 3.5 seconds, and any sucsessive call to the MCD for the same month takes around 0.025 ms.
However, using a time in a different month results in different files being loaded, costing 3.5 seconds again.

## Parallel loading of the MCD modules
To get around the waiting time when calls to the MCD are made for distinct Martian months, the `fmcd.so` interface can be loaded not once but 13 times.
This way, the call to the MCD can be made using the module that always has the files loaded for the same Martian month.
Such Martian month is defined every 30 deg of solar longitude.

Unortunately, this lead to the need to compile the [fmcd.so](fmcd.so) interface 13 times individually, once for each Martian month (plus an extra).
This is because, if the same module is loaded under different names in Python, the same object is actually loaded.
Also, renaming the [fmcd.so](fmcd.so) does not suffise, as the name needs to be defined before the Fortran interface is compiled.
This explains files [fmcd_1.so](fmcd_1.so) to [fmcd_13.so](fmcd_13.so).

A script that explores and times this duplication of the modules can be seen in [test_parallel_mcd.py](test_parallel_mcd.py).

Finally, a Python class can be obtained from [parallel_mcd.py](parallel_mcd.py).
This class can be used to load the modules in parallel and call the appropriate one depending on the Martian month.
An example script that uses this class can be accessed in [test_parallel_mcd.py](test_parallel_mcd.py).