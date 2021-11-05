import numpy as np
import sys
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
from tools import plot_utilities as PU

# Define the altitudes and satellite names for which there is results to analyse
hs = [85, 115, 150]
sat_names = ["CS_0020", "CS_0021", "CS_1020", "CS_1021", "CS_2020", "CS_2021", "CS_2120", "CS_3020", "CS_3021"]

# Loop trough the satellite names and the altitudes
for s_name in sat_names:
    for i_h, h in enumerate(hs):
        # Get the sorted file list corresponding to the satellite name and altitude
        try:
            results_data = open(sys.path[0]+"/SPARTA/setup/results_sparta/%s/stats_%skm.dat" % (s_name, h)).readlines()
        # Print a warning and skip this config if no results file could be found
        except FileNotFoundError:
            print("Warning, it seems that the simulation for %s at %skm was not run yet." % (s_name, h))
            continue

        results = []
        reading_results = False
        for res_line in results_data:
            if reading_results:
                results_row = res_line.split()
                if results_row[0] == "Loop":
                    reading_results = False
                else:
                    results.append(results_row)

            if res_line.strip()[:8] == "Step CPU":
                reading_results = True

        results = np.asarray(results, dtype=float)
        times, fx, fy, fz, ppc = results[:,0], -results[:,-4], -results[:,-3], -results[:,-2], results[:,-1]

        # Plot the force in each direction
        PU.plot_single(times, fx, "Timestep number [-]", "$F_x$ [N]", "SPARTA/fx_%s_%skm" % (s_name, h))
        PU.plot_single(times, fy, "Timestep number [-]", "$F_y$ [N]", "SPARTA/fy_%s_%skm" % (s_name, h))
        PU.plot_single(times, fz, "Timestep number [-]", "$F_y$ [N]", "SPARTA/fz_%s_%skm" % (s_name, h))
        # Print the drag by averaging the force in the x-direction for the last few steps of the simulation
        print("Drag = %.5e N for %s at %.1fkm" % (np.mean(fx[-3:]), s_name, h))

        # Plot the number of particles over time
        PU.plot_single(times, ppc, "Timestep number [-]", "Mean number of particles per cell [-]", "SPARTA/npart_%s_%skm" % (s_name, h))