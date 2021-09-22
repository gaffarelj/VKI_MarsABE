import numpy as np
import glob
from natsort import natsorted
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 13, 'figure.figsize': (10.5, 7), 'savefig.format': 'pdf'})
import sys
sys.path.insert(0,"\\".join(sys.path[0].split("\\")[:-1]))
from tools import plot_utilities as PU

file_list = natsorted(glob.glob("SPARTA/results/dump.*.dat"))
    
times, results = [], []

for res_file in file_list:
    res_lines = open(res_file).readlines()
    results_list = []

    res_vars = [None]*8

    t_step = float(res_lines[1])
    times.append(t_step)

    for line in res_lines[9:]:
        results_list.append([float(_res) for _res in line.split(" ")[:-1]])

    results_list = np.array(results_list)

    for i in range(8):
        res_vars[i] = sum(results_list[:,i])
        
    results.append(res_vars)

results = np.array(results)
res_labels = ["tot_press", "fx", "fy", "fz", "px", "py", "pz", "etot"]

res_dict = dict()
for i in range(7):
    res_dict[res_labels[i]] = results[:,i]
PU.plot_single(times, -res_dict["fx"], "Time [s]", "$F_x$ [N]", "SPARTA/f_x")
PU.plot_single(times, -res_dict["fy"], "Time [s]", "$F_y$ [N]", "SPARTA/f_y")
PU.plot_single(times, -res_dict["fz"], "Time [s]", "$F_y$ [N]", "SPARTA/f_z")

print("Drag = %.5e N" % (-res_dict["fx"][-1]))