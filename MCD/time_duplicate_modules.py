import numpy as np
import time
# Import the same Fortran interface compiled with different names in Python.
# This is not an elegant solution, but a very efficient one.
# This way, the required file can be loaded only once per module,
# and data from the MCD can be obtained by switching modules depending on the time
from fmcd_1 import call_mcd as call_mcd_1
from fmcd_2 import call_mcd as call_mcd_2
from fmcd_3 import call_mcd as call_mcd_3
from fmcd_4 import call_mcd as call_mcd_4
from fmcd_5 import call_mcd as call_mcd_5
from fmcd_6 import call_mcd as call_mcd_6
from fmcd_7 import call_mcd as call_mcd_7
from fmcd_8 import call_mcd as call_mcd_8
from fmcd_9 import call_mcd as call_mcd_9
from fmcd_10 import call_mcd as call_mcd_10
from fmcd_11 import call_mcd as call_mcd_11
from fmcd_12 import call_mcd as call_mcd_12
from fmcd_13 import call_mcd as call_mcd_13

zkey = 1                    # xz is thus the radial distance from the centre of Mars
xz = 3389.5e3+300e3         # [m], distance from the centre of Mars (300km altitude for tests)
xlon = 137.4                # [deg], East longitude
xlat = -4.6                 # [deg], Latitude
hireskey = 0                # use the lower resolution grid of 5.625x3.75 deg
localtime = 0               # mendatory if datekey=0
dset = "/mnt/c/MCD/data/"   # path to the MCD data
scena = 1                   # Climatology Scenario, solar EUV average conditions; use 2 for minimum solar conditions, and 3 for maximum
perturkey = 0               # do not additional perturbations
seedin = 1                  # (not) used for the random number generator
gwlength = 0                # (not) used for gravity waves length
# select which external variables to return
extvarkeys = np.zeros(100)
n_species = (56, 65)
extvarkeys[n_species[0]:n_species[1]+1] = 1
datekey = 1                 # 0 = Julian, 1 = Ls

# Load all of the MCD data in different modules
t0 = time.time()
call_mcd_list = []
# The following two arrays contain the starting time in Ls and Julian date at which the given module should be used
call_mcd_Ls = np.arange(0, 330.01, 30)
for i_module, xdate in enumerate(call_mcd_Ls):
    # Select the appropriate module, and put it in a dict with the starting Julian epoch as the key
    call_mcd = eval("call_mcd_%i" % (i_module+1))
    call_mcd_list.append(call_mcd)
    # Load the MCD data
    call_mcd(zkey,xz,xlon,xlat,hireskey,datekey,xdate,localtime,dset,scena,perturkey,seedin,gwlength,extvarkeys)
call_time = time.time() - t0
print("It took %.2f seconds to load the modules from the entire Martian year" % call_time) # takes around 40 seconds to load all modules

## Call the MCD for different times, and time the calls
xlon, xlat = 137.4, -4.6  # reset the latitude and longitude to fixed values
N = int(1e5)
t0 = time.time()
for i in range(N):
    xdate = np.random.uniform(call_mcd_Ls[0], call_mcd_Ls[-1])
    i_module = np.where(call_mcd_Ls - xdate <= 15)
    i_module = i_module[0][-1]
    call_mcd = call_mcd_list[i_module]
    pres,dens,temp,zonwind,merwind,meanvar,extvars,seedout,ier = call_mcd(zkey,xz,xlon,xlat,hireskey,datekey,xdate,localtime,dset,scena,perturkey,seedin,gwlength,extvarkeys)
call_time = time.time() - t0
print("Calling the MCD with different times took %.5f seconds." % (call_time))
print("This is an average of %.5f milliseconds per call" % (call_time*1e3/N)) # takes around 0.025 ms per call with pre-loaded modules