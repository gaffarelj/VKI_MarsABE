import numpy as np
import glob
from natsort import natsorted
import sys
sys.path.insert(0,"\\".join(sys.path[0].split("\\")[:-1]))
from tools import plot_utilities as PU

check_part_cells = True     # Set to True to check the number of particles in each cells

dt = 5                      # Measurements time step

# Define the altitudes and satellite names for which there is results to analyse
hs = [85, 115, 150]
sat_names = ["CS_0020", "CS_0021", "CS_1020", "CS_1021", "CS_2020", "CS_2021", "CS_2120", "CS_3020", "CS_3021"]
# Loop trough the satellite names and the altitudes
for s_name in sat_names:
    for h in hs:
        # Get the sorted file list corresponding to the satellite name and altitude
        file_list = natsorted(glob.glob("SPARTA/setup/results_sparta/%s/force_%skm.*.dat" % (s_name, h)))
        
        if check_part_cells:
            file_list_pc = natsorted(glob.glob("SPARTA/setup/results_sparta/%s/npart_%skm.*.dat" % (s_name, h)))

        # Prepare the result lists
        times, results = [], []

        # Go trough each result file
        for i, res_file in enumerate(file_list):
            times.append(i*dt)
            # Read the file
            data = np.loadtxt(res_file, skiprows=9)
            # Sum data from all surface
            forces = np.sum(data, axis=0)
            # Save forces
            results.append(forces)

        # Convert the results to a dict
        results = np.array(results)

        # Plot the force in each direction
        PU.plot_single(times, -results[:,0], "Time [s]", "$F_x$ [N]", "SPARTA/fx_%s_%skm" % (s_name, h))
        PU.plot_single(times, -results[:,1], "Time [s]", "$F_y$ [N]", "SPARTA/fy_%s_%skm" % (s_name, h))
        PU.plot_single(times, -results[:,2], "Time [s]", "$F_y$ [N]", "SPARTA/fz_%s_%skm" % (s_name, h))

        # Print the drag by averaging the force in the x-direction for the last 10% of the simulation
        print("Drag = %.5e N for %s at %.1fkm" % (-np.mean(results[-len(times)//10:,0]), s_name, h))