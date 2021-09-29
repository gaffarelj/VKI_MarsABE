import numpy as np
import glob
from natsort import natsorted
import sys
sys.path.insert(0,"\\".join(sys.path[0].split("\\")[:-1]))
from tools import plot_utilities as PU

# Define the altitudes and satellite names for which there is results to analyse
hs = [85, 115, 150]
sat_names = ["CS_0020", "CS_0021", "CS_1020", "CS_1021", "CS_2020", "CS_2021", "CS_2120", "CS_3020", "CS_3021"]
# Loop trough the satellite names and the altitudes
for s_name in sat_names:
    for h in hs:
        # Get the sorted file list corresponding to the satellite name and altitude
        file_list = natsorted(glob.glob("SPARTA/setup/results_sparta/%s_%skm.*.dat" % (s_name, h)))
        
        # Prepare the result lists
        times, results = [], []

        # Go trough each result file
        for res_file in file_list:
            # Read the file
            res_lines = open(res_file).readlines()

            # Prepare results list for this file
            results_list = []
            res_vars = [None]*8

            # Extract the current time step
            t_step = float(res_lines[1])
            times.append(t_step)

            # Go trough each line in the file
            for line in res_lines[9:]:
                # Get all values as float, and add them to the results list
                results_list.append([float(_res) for _res in line.split(" ")[:-1]])
            results_list = np.array(results_list)

            # Sum the results from each of the satellite surface
            for i in range(4):#8):
                res_vars[i] = sum(results_list[:,i])
            
            # Save the summed results
            results.append(res_vars)

        # Convert the results to a dict
        results = np.array(results)
        res_labels = ["tot_press", "fx", "fy", "fz"]#, "px", "py", "pz", "etot"]
        res_dict = dict()
        for i in range(4):#7):
            res_dict[res_labels[i]] = results[:,i]

        # Plot the force in each direction
        PU.plot_single(times, -res_dict["fx"], "Time [s]", "$F_x$ [N]", "SPARTA/fx_%s_%skm" % (s_name, h))
        PU.plot_single(times, -res_dict["fy"], "Time [s]", "$F_y$ [N]", "SPARTA/fy_%s_%skm" % (s_name, h))
        PU.plot_single(times, -res_dict["fz"], "Time [s]", "$F_y$ [N]", "SPARTA/fz_%s_%skm" % (s_name, h))

        # Print the drag by averaging the force in the x-direction for the last 10% of the simulation
        print("Drag = %.5e N for %s at %.1fkm" % (-np.mean(res_dict["fx"][-len(times)//10:]), s_name, h))