import numpy as np
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 13, 'figure.figsize': (10.5, 7), 'savefig.format': 'pdf'})
import sys
sys.path.insert(0,"\\".join(sys.path[0].split("\\")[:-2])) # get back to uppermost level of the project
from MCD.parallel_mcd import parallel_mcd as PMCD
from tools import plot_utilities as PU

# Load the MCD interface
mcd = PMCD()

altitudes = np.arange(50, 351, 0.1)
legends = mcd.species_name.copy()
legends.insert(0, "Air")
all_dens = np.empty((len(legends), len(altitudes)))
all_frac = np.empty((len(legends), len(altitudes)))

for i, h in enumerate(altitudes):
    mcd.call(h=h*1e3)
    all_dens[0,i] = mcd.dens
    for j, dens in enumerate(mcd.species_dens):
        all_dens[j+1,i] = dens
        all_frac[j+1,i] = mcd.species_frac[j]

# Remove O3 (density too low)
legends.pop(7)
all_dens = np.delete(all_dens, 7, axis=0)
all_frac = np.delete(all_frac, [0, 7], axis=0)

PU.plot_multiple(all_dens, [altitudes]*len(legends), "Density [kg/m$^3$]", "Altitude [km]", "MCD/all_density_map", legends, xlog=True, colors=["black"], ylim=[altitudes[0], altitudes[-1]])

PU.plot_multiple(all_frac, [altitudes]*len(legends), "Fraction [mol/mol]", "Altitude [km]", "MCD/all_fractions_map", legends[1:], ylim=[altitudes[0], altitudes[-1]])