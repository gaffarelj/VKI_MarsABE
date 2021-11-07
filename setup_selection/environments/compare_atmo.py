import numpy as np
import sys
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
from tools import plot_utilities as PU
from GRAM import call_GRAM
from MCD import parallel_mcd as pmcd

GRAM_atmo = call_GRAM.GRAM_atmo()
MCD_atmo = pmcd.parallel_mcd()

def expo_model(h):
    # Exponential parameters taken from http://link.springer.com/content/pdf/10.1007%2F978-3-540-73647-9_3.pdf
    density_scale_height, density_at_zero_altitude = 7295, 0.0525
    return density_at_zero_altitude * np.exp(-h/density_scale_height)

def densitiy_all_models(h=125e3, lat=0, lon=0, time=[60, 10]):
    density_expo = expo_model(h)
    density_GRAM = GRAM_atmo.get_density(h, lon, lat, time, time_is_JD=False)
    density_MCD = MCD_atmo.density(h, lat, lon, time, time_is_JD=False)
    return density_expo, density_GRAM, density_MCD

# Compare density vs altitude
altitudes = np.arange(85e3, 150.01e3, 100)
densities_expo, densities_GRAM, densities_MCD = [], [], []
for h in altitudes:
    all_dens = densitiy_all_models(h=h)
    densities_expo.append(all_dens[0]), densities_GRAM.append(all_dens[1]), densities_MCD.append(all_dens[2])
PU.plot_multiple([densities_expo, densities_GRAM, densities_MCD], [altitudes/1e3]*3, "Density [kg/m3]", "Altitude [km]", \
    "atmosphere_comparison/altitude", ["Exponential model", "Mars GRAM", "Mars Climate Database"], xlog=True)

# Compare density vs time of day
hrs = np.arange(0, 24.01, 0.1)
densities_expo, densities_GRAM, densities_MCD = [], [], []
for hr in hrs:
    all_dens = densitiy_all_models(time=[60, hr])
    densities_expo.append(all_dens[0]), densities_GRAM.append(all_dens[1]), densities_MCD.append(all_dens[2])
PU.plot_multiple([densities_expo, densities_GRAM, densities_MCD], [hrs]*3, "Density [kg/m3]", "Time of the day [hr]", \
    "atmosphere_comparison/dayime", ["Exponential model", "Mars GRAM", "Mars Climate Database"], xlog=True)

# Compare density vs latitude
lats = np.arange(-90, 90.01, 0.1)
densities_expo, densities_GRAM, densities_MCD = [], [], []
for lat in lats:
    all_dens = densitiy_all_models(lat=np.deg2rad(lat))
    densities_expo.append(all_dens[0]), densities_GRAM.append(all_dens[1]), densities_MCD.append(all_dens[2])
PU.plot_multiple([densities_expo, densities_GRAM, densities_MCD], [lats]*3, "Density [kg/m3]", "Latitude [deg]", \
    "atmosphere_comparison/latitude", ["Exponential model", "Mars GRAM", "Mars Climate Database"], xlog=True)

# Compare density vs longitude
lons = np.arange(-180, 180.01, 0.1)
densities_expo, densities_GRAM, densities_MCD = [], [], []
for lon in lons:
    all_dens = densitiy_all_models(lon=np.deg2rad(lon))
    densities_expo.append(all_dens[0]), densities_GRAM.append(all_dens[1]), densities_MCD.append(all_dens[2])
PU.plot_multiple([densities_expo, densities_GRAM, densities_MCD], [lons]*3, "Density [kg/m3]", "Latitude [deg]", \
    "atmosphere_comparison/longitude", ["Exponential model", "Mars GRAM", "Mars Climate Database"], xlog=True)

# Compare density vs solar longitude
Ls_s = np.arange(0, 360, 0.5)
densities_expo, densities_GRAM, densities_MCD = [], [], []
for Ls in Ls_s:
    all_dens = densitiy_all_models(time=[Ls, 10])
    densities_expo.append(all_dens[0]), densities_GRAM.append(all_dens[1]), densities_MCD.append(all_dens[2])
PU.plot_multiple([densities_expo, densities_GRAM, densities_MCD], [Ls_s]*3, "Density [kg/m3]", "Solar longitude [deg]", \
    "atmosphere_comparison/Ls", ["Exponential model", "Mars GRAM", "Mars Climate Database"], xlog=True)