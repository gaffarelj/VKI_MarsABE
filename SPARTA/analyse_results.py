import numpy as np
import glob
from natsort import natsorted
import gzip
import sys
sys.path.insert(0,"\\".join(sys.path[0].split("\\")[:-1]))
from tools import plot_utilities as PU

check_part_cells = True     # Set to True to check the number of particles in each cells

meas_dt = [25, 25, 25]      # Epoch between measurements

# Define the altitudes and satellite names for which there is results to analyse
hs = [85, 115, 150]
sat_names = ["CS_0020", "CS_0021", "CS_1020", "CS_1021", "CS_2020", "CS_2021", "CS_2120", "CS_3020", "CS_3021"]
# Loop trough the satellite names and the altitudes
for s_name in sat_names:
    for i_h, h in enumerate(hs):
        dt = meas_dt[i_h]
        # Get the sorted file list corresponding to the satellite name and altitude
        file_list = natsorted(glob.glob("SPARTA/setup/results_sparta/%s/force_%skm.*.gz" % (s_name, h)))
        
        if check_part_cells:
            file_list_pc = natsorted(glob.glob("SPARTA/setup/results_sparta/%s/npart_%skm.*.gz" % (s_name, h)))

        # Prepare the result lists
        times, results, times_np, results_np_mean, results_np_std = [], [], [], [], []
        
        # Go trough each result file
        for i, res_file in enumerate(file_list):
            print("Reading force file %i/%i..." % (i+1, len(file_list)), end="\r")
            times.append(i*dt)
            # Read the file
            data = np.array([l.decode().strip().split(" ") for l in gzip.open(res_file, "rb").readlines()[9:]], dtype=float)
            # Sum data from all surface
            forces = np.sum(data[:,:3], axis=0)
            # Save forces
            results.append(forces)
        print()
        # Convert the results to a dict
        results = np.array(results)

        # Plot the force in each direction
        PU.plot_single(times, -results[:,0], "Timestep number [-]", "$F_x$ [N]", "SPARTA/fx_%s_%skm" % (s_name, h))
        PU.plot_single(times, -results[:,1], "Timestep number [-]", "$F_y$ [N]", "SPARTA/fy_%s_%skm" % (s_name, h))
        PU.plot_single(times, -results[:,2], "Timestep number [-]", "$F_y$ [N]", "SPARTA/fz_%s_%skm" % (s_name, h))

        # Print the drag by averaging the force in the x-direction for the last 10% of the simulation
        print("Drag = %.5e N for %s at %.1fkm" % (-np.mean(results[-len(times)//10:,0]), s_name, h))

        if check_part_cells:
            max_T = 0
            for i, res_file in enumerate(file_list_pc):
                print("Reading grid file %i/%i..." % (i+1, len(file_list_pc)), end="\r")
                times_np.append(i*dt)
                # Read the file
                data = np.array([l.decode().strip().split(" ") for l in gzip.open(res_file, "rb").readlines()[9:]], dtype=float)
                # Get the average number of particles per cell
                mean_npart = np.mean(data[:,1])
                # Get the std of the number of particles per cell
                std_npart = np.std(data[:,1])
                # Get the maximum temperature
                #max_T = max(max_T, max(data[:,2]))
                # Save number of particles data
                results_np_mean.append(mean_npart), results_np_std.append(std_npart)
            print()
            # Convert the results to a dict
            results_np_mean, results_np_std = np.array(results_np_mean), np.array(results_np_std)

            #print("Maximum temperature of %.5e K" % max_T)

            # Plot the number of particles over time
            PU.plot_multiple([times_np]*3, [results_np_mean, results_np_mean+results_np_std], \
                "Timestep number [-]", "Number of particles per cell [-]", "SPARTA/npart_%s_%skm" % (s_name, h), legends=["$\mu$", "$\mu$+$\sigma$"])