import numpy as np
import matplotlib
import pygmo
from natsort import natsorted
import glob
import time
import os
import sys
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
from optimisation import drag_comp_problem as DCp
from utils import sat_models as SM
from tools import plot_utilities as PU

## Select the thrust model
print("The following thrust models can be used:")
print(" 1: BHT-100 Hall thruster (on when power > 107 W)")
print(" 2: μNRIT 2.5 Radiofrequency ion thruster with Xenon tank (on when power > 13.1 W)")
print(" 3: μNRIT 2.5 Radiofrequency ion thruster with atmosphere-breathing inlet (on when power > 13.1 W and engine inlet mass flow > 1.456e-8 kg/s)")
thrust_model = None
while thrust_model not in [1, 2, 3]:
    try:
        thrust_model = int(input("Thrust model selection [1, 2, 3]: "))
    except ValueError:
        thrust_model = None
satellites = SM.satellites if thrust_model == 3 else SM.satellites_with_tank
plots_path = "optimisation/ABE/" if thrust_model == 3 else "optimisation/with_tank/"

# Setup the design variables range
min_h_p, max_h_p = 85e3, 150e3
min_h_a, max_h_a = 85e3, 500e3
min_e, max_e = 0, 0.1
min_i, max_i = 0, np.pi/2
min_omega, max_omega = 0, np.pi
min_Omega, max_Omega = 0, 2*np.pi
min_sat_i, max_sat_i = 0, len(satellites) - 1
design_var_range = (
    [min_h_p, min_h_a, min_i, min_omega, min_Omega, min_sat_i],
    [max_h_p, max_h_a, max_i, max_omega, max_Omega, max_sat_i]
)

# Setup the optimisation problem
fitness_weights = [1, 1, 1, 1]
fitness_names = ["Mean power", "Periapsis decay", "Mean altitude", "Mean Drag/Thrust"]
current_HT_problem = DCp.WT_problem(design_var_range, fitness_weights, thrust_model=thrust_model, verbose=False)
problem = pygmo.problem(current_HT_problem)

# Select whether to run the optimisation or load the latest result file
run_opti = (input("Run the optimisation ? ([y]/n) (if no, load latest saved results): ").lower().strip() in ["", "y"])
# Configure population size and number of generations
sizes = [60, 10, 10, 10]
algo_idx = 0
pop_size = sizes[algo_idx]
n_generations = 100
if run_opti:    # Run a new optimisation
    # Setup Pygmo
    seed = 12345
    algo_list = [
        pygmo.algorithm(pygmo.nsga2(seed=seed, cr=0.95, eta_c=10, m=0.001, eta_m=1)), # seems best
        pygmo.algorithm(pygmo.moead(seed=seed, neighbours=5)),
        pygmo.algorithm(pygmo.nspso(seed=seed)),
        pygmo.algorithm(pygmo.ihs(seed=seed)) # maybe, but only 1 pop improved by generation (-> increase gen number)
    ]

    print("Generating starting population (of size %i)..." % pop_size)
    pop = pygmo.population(problem, size=pop_size, seed=seed)
    algo = algo_list[algo_idx]

    opti_hist = []

    f_name = "WT_%i-%i-%s-%s" % (thrust_model, pop_size, time.strftime("%d%m%y_%H%M%S"), seed)
    f_path = sys.path[0] + "/optimisation/results/" + f_name

    # Run the optimisation
    t0 = time.time()
    for i in range(1,n_generations+1):
        print("Running generation %2d / %2d" % (i, n_generations))
        # Evolve the population
        pop = algo.evolve(pop)
        
        # Show time it took
        dt, t0 = time.time() - t0, time.time()
        print(" -> took %.1f seconds" % dt)

        # Add best results to historic
        f = np.array(pop.get_f())
        best_f = [min(f[:,i]) for i in range(len(fitness_weights))]
        best_f.append(min(np.mean(np.fabs(f), axis=1)))
        opti_hist.append(best_f)
        print(" -> best fitness is", best_f)

        # Save the results
        np.savez(f_path+"_%i"%i, inputs=DCp.FIT_INPUTS, results=np.array(DCp.FIT_RESULTS), opti_hist=np.array(opti_hist))
        # Remove results from previous generation
        if i > 1:
            os.remove(f_path+"_%s.npz"%(i-1))
        
    fit_results, fit_inputs, opti_hist = np.array(DCp.FIT_RESULTS), DCp.FIT_INPUTS, np.array(opti_hist)

else:       # Load the last saved results
    file_list = natsorted(glob.glob(sys.path[0]+"/optimisation/results/WT_%i-%i-*.npz" % (thrust_model, pop_size)))
    if len(file_list) == 0:
        raise FileNotFoundError("It appears that no optimisation for this given thrust model and optimiser settings were already run and saved.")
    last_results = np.load(file_list[-1])
    fit_inputs = last_results["inputs"]
    fit_results = last_results["results"]
    opti_hist = last_results["opti_hist"]
    s_names = fit_inputs[:,0].reshape((len(fit_inputs),1))
    fit_inputs = np.array(fit_inputs[:,1:], dtype=float)
    fit_inputs = np.concatenate([s_names, fit_inputs], axis=1, dtype=object)

# Find the optimum if the sum of fitnesses is used
idx_best = np.where(np.sum(fit_results[:,:len(fitness_weights)], axis=1) == np.min(np.sum(fit_results[:,:len(fitness_weights)], axis=1)))[0][0]
optimum_input, optimum_result = fit_inputs[idx_best], fit_results[idx_best]

# Print the results
print("Optimum (when summing the fitness): %s \n  with initial state: h_p=%.2f km, h_a=%.2f km, i=%.2f, omega=%.2f, Omega=%.2f" % \
    (satellites[optimum_input[0]], min(optimum_input[1:3])/1e3, max(optimum_input[1:3])/1e3, \
        np.rad2deg(optimum_input[3]), np.rad2deg(optimum_input[4]), np.rad2deg(optimum_input[5])))
print("Resulting optimum fitness:", optimum_result[:len(fitness_weights)]*np.array(fitness_weights))
print("Resulting characteristics: mean power=%.2f W, total decay=%.2f km, mean altitude=%3d km, mean T/D=%.2f" % \
    (optimum_result[-4], optimum_result[-3]/1e3, optimum_result[-2]/1e3, optimum_result[-1]))
print("Explored %i different possibilities." % len(fit_inputs))

# Plot the fitness progress over time
PU.plot_multiple([list(range(1, opti_hist.shape[0]+1))]*(len(fitness_weights)+1), opti_hist.T, "Generation number", "Best fitness", \
    plots_path+"history", legends=fitness_names+["Average fitness"], colors=["darkorange", "seagreen", "royalblue", "#202020"], \
    lstyle=["solid"]*len(fitness_weights)+["dashed"])

# Generate Pareto fronts
for zoomed in [False]:#, True]:
    if zoomed:
        # Zoom in on region where mean altitude is below 150km and where periapsis decay is within -5km to 5km
        idx_remove = np.where((fit_results[:,-3] > 5e3) | (fit_results[:,-3] < -5e3) | (fit_results[:,-2] > 150e3))
        prefix = "_zoomed"
    else:
        # Filter periapsis decays above 100km
        idx_remove = np.where(fit_results[:,-3] > 100e3)
        prefix = ""
    # Remove at selected indexes
    obj_power, obj_decay, obj_h, obj_DT = np.delete(fit_results[:,-4], idx_remove), np.delete(fit_results[:,-3], idx_remove), np.delete(fit_results[:,-2], idx_remove), np.delete(fit_results[:,-1], idx_remove)
    s_names = np.delete(fit_inputs[:,0], idx_remove)
    # Select color as a function of the satellite
    s_numbers = np.arange(0, len(satellites), 1)
    s_nn_map = dict(zip(list(satellites.keys()), s_numbers))
    ## Make the plots
    # Classic Pareto plots, 2 objectives
    PU.plot_single(obj_power, obj_decay/1e3, "Mean Power [W]", "Periapsis decay [km]", plots_path+"Pareto_Pd"+prefix, \
        scatter=True, add_front=True, front_sign=[-1,1])
    PU.plot_single(obj_power, obj_h/1e3, "Mean Power [W]", "Mean altitude [km]", plots_path+"Pareto_Ph"+prefix, \
        scatter=True, add_front=True, front_sign=[-1,1])
    PU.plot_single(obj_h/1e3, obj_decay/1e3, "Mean altitude [km]", "Periapsis decay [km]", plots_path+"Pareto_hd"+prefix, \
        scatter=True, add_front=True)
    PU.plot_single(obj_h/1e3, obj_DT, "Mean altitude [km]", "Mean Drag/Thrust [-]", plots_path+"Pareto_hT"+prefix, \
        scatter=True, add_front=True)
    PU.plot_single(obj_power, obj_DT, "Mean Power [W]", "Mean Drag/Thrust [-]", plots_path+"Pareto_PT"+prefix, \
        scatter=True, add_front=True, front_sign=[-1,1])
    PU.plot_single(obj_DT, obj_decay/1e3, "Mean Drag/Thrust [-]", "Periapsis decay [km]", plots_path+"Pareto_Td"+prefix, \
        scatter=True, add_front=True)
    # Plot decay vs mean altitude with power in the colormap
    power_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
                    'trunc({n},{a:.2f},{b:.2f})'.format(n="plasma", a=0.0, b=0.9),
                    matplotlib.pyplot.get_cmap("plasma")(np.linspace(0.0, 0.9, 10)))
    PU.plot_single(obj_h/1e3, obj_decay/1e3, "Mean altitude [km]", "Periapsis decay [km]", plots_path+"Pareto_hdP"+prefix, \
        scatter=True, add_front=True, z_data=obj_power, z_label="Mean power [W]", cmap=power_cmap)
    # Plot decay vs mean altitude with satellite name in the colormap
    s_name_cmap = matplotlib.colors.ListedColormap(['red', 'green', 'blue', 'yellow', 'orange'])
    bounds = [0, 1, 2, 3, 4, 5]
    norm =  matplotlib.colors.BoundaryNorm(bounds, s_name_cmap.N)
    PU.plot_single(obj_h/1e3, obj_decay/1e3, "Mean altitude [km]", "Periapsis decay [km]", plots_path+"Pareto_hdS"+prefix, \
        scatter=True, add_front=True, z_data=[s_nn_map[s_n] for s_n in s_names], z_label="Satellite", \
            cmap=s_name_cmap, cticks=[0.5, 1.5, 2.5, 3.5, 4.5], \
            clabels=list(satellites.keys()), NB=(norm, bounds))
    # Plot decay vs mean altitude with D/T in the colormap
    PU.plot_single(obj_h/1e3, obj_decay/1e3, "Mean altitude [km]", "Periapsis decay [km]", plots_path+"Pareto_hdT"+prefix, \
        scatter=True, add_front=True, z_data=np.clip(obj_DT, 0, 10), z_label="Mean Drag/Thrust [-]", cmap="rainbow")

idx_remove = np.where(fit_results[:,-3] >= 100e3)
# Make a Panda dataframe from the results
import pandas as pd
import plotly.express as px
s_pd = pd.Series(np.delete(fit_results[:,-3], idx_remove), name="Periapsis decay")
s_mh = pd.Series(np.delete(fit_results[:,-2], idx_remove), name="Mean altitude")
s_mTD = pd.Series(np.delete(fit_results[:,-1], idx_remove), name="Mean D/T")
s_mp = pd.Series(np.delete(fit_results[:,-4], idx_remove), name="Mean power")
s_sn = pd.Series(np.delete(fit_inputs[:,0], idx_remove), name="Satellite")
s_i_hp = pd.Series(np.delete(fit_inputs[:,1], idx_remove), name="Initial h_p")
s_i_ha = pd.Series(np.delete(fit_inputs[:,2], idx_remove), name="Initial h_a")
s_i_i = pd.Series(np.delete(fit_inputs[:,3], idx_remove), name="Initial i")
s_i_omega = pd.Series(np.delete(fit_inputs[:,4], idx_remove), name="Initial omega")
s_i_Omega = pd.Series(np.delete(fit_inputs[:,5], idx_remove), name="Initial Omega")
df = pd.concat([s_pd, s_mh, s_mTD, s_mp, s_sn, s_i_hp, s_i_ha, s_i_i, s_i_omega, s_i_Omega], axis=1)
# Make interactive plot
fig = px.scatter(df, x="Mean altitude", y="Periapsis decay", hover_name="Satellite", \
    hover_data=["Mean D/T", "Mean power", "Initial h_p", "Initial h_a", "Initial i", "Initial omega", "Initial Omega"])
fig.write_html(sys.path[0]+"/figures/"+plots_path+"interactive_pareto.html")
