from fmcd import call_mcd
import numpy as np
import time

zkey = 1                    # xz is thus the radial distance from the centre of Mars
xz = 3389.5e3+300e3         # [m], distance from the centre of Mars (300km altitude for tests)
xlon = 137.4                # [deg], East longitude
xlat = -4.6                 # [deg], Latitude
hireskey = 0                # use the lower resolution grid of 5.625x3.75 deg
datekey = 0                 # xdate is thus the Julian date (2459580.5 is January 1st 2025, midnight; 2461041.5 is January 1st 2026, midnight)
xdate = 2460676.5           # Julian date
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

# Call the MCD, and time the call
t0 = time.time()
pres,dens,temp,zonwind,merwind,meanvar,extvars,seedout,ier = call_mcd(zkey,xz,xlon,xlat,hireskey,datekey,xdate,localtime,dset,scena,perturkey,seedin,gwlength,extvarkeys)
call_time = time.time() - t0
print("Calling the MCD took %.5f seconds." % call_time) # takes between 3.4 and 3.6 seconds

# Call the MCD for different geolocalisation and time the calls
N = int(1e5)
t0 = time.time()
for i in range(N):
    xlat = np.random.uniform(-90, 90)
    xlon = np.random.uniform(-180, 180)
    pres,dens,temp,zonwind,merwind,meanvar,extvars,seedout,ier = call_mcd(zkey,xz,xlon,xlat,hireskey,datekey,xdate,localtime,dset,scena,perturkey,seedin,gwlength,extvarkeys)
call_time = time.time() - t0
print("Calling the MCD with different geolocalisation took %.5f seconds." % (call_time))
print("This is an average of %.5f milliseconds per call" % (call_time*1e3/N)) # takes around 0.025 ms per call