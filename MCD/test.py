# IMPORTANT: on my machine, the MCD was compiled for Linux subsystem on windows. Thus, this file shall be run in Lunix.

from mcd import mcd
import numpy as np

req = mcd()
req.lat = -4.6 # latitude
req.lon = 137.4 # longitude
req.loct = 0 # local time
req.datekey = 0
req.xz = 1. # vertical coordinate
req.xdate = 2460676.5 # Julian date
req.zkey = 1
req.xz = 3389.5e3+300e3
req.perturkey = 0
req.scena = 0
req.update()

#for Ls in np.arange(0, 360, 30):
#    req.xdate = Ls
#    req.update()
#    req.printcoord()
#    density = req.dens # density [kg/m3]
#    print("Density = %.5f kg/m3" % density)

density = req.dens # density [kg/m3]
wind_zon = req.zonwind # zonal wind [m/s]
wind_mer = req.merwind # meridional wind [m/s]

species_name = ["CO2", "N2", "Ar", "CO", "O", "O2", "O3", "H", "H2"]
species_frac = []
for n in range(56, 65):
    species_frac.append(req.extvar[n])
species_dict = dict(zip(species_name, species_frac))

print("Density = %.5f kg/m3" % density)
print("Wind = %.5f°E/%.5f°N m/s" % (wind_zon, wind_mer))
print("Species [mol/mol] = %s" % species_dict)