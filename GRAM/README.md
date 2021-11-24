# Mars GRAM

The [call_GRAM.py](call_GRAM.py) script can be used to query the atmospheric density of Mars as a function of altitude, latitude, and time.

These density are extrapolated based on the data from Mars GRAM 2010, by NASA. This extrapolation is made using a 3d-interpolation of an atmosphere lookup table.

The variable `file_path` at the beginning of the Python file refers to the path of the `tpdmsy11.txt` data file that is distributed with Mars GRAM 2010.
It can be requested to NASA by filling the form specified at [this page](https://software.nasa.gov/software/MFS-33158-1).
