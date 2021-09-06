from fmcd import call_mcd
import numpy as np
import time

zkey = 1                    # xz is thus the radial distance from the centre of Mars
xz = 3389.5e3+300e3         # [m], distance from the centre of Mars (300km altitude for tests)
xlon = 137.4                # [deg], East longitude
xlat = -4.6                 # [deg], Latitude
hireskey = 0                # use the lower resolution grid of 5.625x3.75 deg
datekey = 0                 # xdate is thus the Julian date (2459580.5 is January 1st 2025, midnight)
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

# Call the MCD
pres,dens,temp,zonwind,merwind,meanvar,extvars,seedout,ier = call_mcd(zkey,xz,xlon,xlat,hireskey,datekey,xdate,localtime,dset,scena,perturkey,seedin,gwlength,extvarkeys)

# Extract the volumetric ratio of each species in the atmosphere
species_name = ["CO2", "N2", "Ar", "CO", "O", "O2", "O3", "H", "H2"] # note: He also accessible if required later, at index 77
species_frac = []
for n in range(*n_species):
    species_frac.append(extvars[n])
species_dict = dict(zip(species_name, species_frac))

# Convert the volumetric ratio of each species to their density (TODO)

# Print the results
print("Density = %.5e [kg/m3]" % dens)
print("Wind = %.5f E / %.5f N [m/s]" % (zonwind, merwind))
print("Species [mol/mol] = %s" % species_dict)