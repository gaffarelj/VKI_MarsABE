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


def ask_choice(q, r, t):
    # Ask user the question `q`, the input has to be in range `r`, and be of the type `t`
    choice = None
    while choice is None:
        try:
            choice = t(input(q))
            if r == "yn":
                if choice.lower().strip() == "y":
                    return True
                elif choice.lower().strip() == "n":
                    return False
                else:
                    choice = None
            else:
                if choice < r[0] or choice > r[1]:
                    choice = None
        except ValueError:
            choice = None
    return choice

def get_saved_res(raise_error=True):
    # Return the last saved results for the current set of problem parameters
    file_list = natsorted(glob.glob("-".join(f_path.split("-")[:-1]) + "*.npz"))
    if len(file_list) == 0:
        error = "It appears that no optimisation was already run to at least 1 generation with the following parameters:\n"+\
            " * %i objectives\n * thrust model = %i\n * use of battery = %s\n * ionisation efficiency = %.2f\n * population size = %i" % \
                (4 if all_obj else 2, thrust_model, "V" if use_battery else "X", ionisation_efficiency, pop_size)
        if raise_error:
            raise FileNotFoundError(error)
        else:
            print(error)
            return None, None, None, None
    last_results = np.load(file_list[-1])
    fit_inputs = last_results["inputs"]
    fit_results = last_results["results"]
    opti_hist = last_results["opti_hist"]
    return fit_inputs, fit_results, opti_hist, file_list[-1]

if __name__ == "__main__":
    ## Select the thrust model
    all_obj = ask_choice("Use all objectives (y) or only the altitude and periapsis decay objectives (n)? (y/n): ", "yn", str)
    print("The following thrust models can be used:")
    print(" 2: ??NRIT 2.5 Radiofrequency ion thruster with Xenon tank (on when power > 13.1 W)")
    print(" 3: ??NRIT 2.5 Radiofrequency ion thruster with atmosphere-breathing inlet (on when power > 13.1 W and engine inlet mass flow > 1.456e-8 kg/s)")

    thrust_model = ask_choice("Thrust model selection [2, 3]: ", [2, 3], int)
    if thrust_model == 3:
        satellites = SM.satellites
        ionisation_efficiency = ask_choice("Ionisation efficiency (in ]0, 1[): ", [0, 1], float)
    else:
        satellites = SM.satellites_with_tank
        ionisation_efficiency = 1
    use_battery = ask_choice("Use the battery? (y/n): ", "yn", str)
    plots_path = "optimisation/DC_%s_%i-%s-%.2f_" % ("4O" if all_obj else "2O", thrust_model, "V" if use_battery else "X", ionisation_efficiency)

    # Setup the design variables range
    min_h_p, max_h_p = 85e3, 150e3
    min_h_a, max_h_a = 85e3, 500e3
    min_i, max_i = 0, np.pi/2
    min_omega, max_omega = 0, np.pi
    min_Omega, max_Omega = 0, 2*np.pi
    min_sat_i, max_sat_i = 0, len(satellites) - 1
    if all_obj:
        design_var_range = (
            [min_h_p, min_h_a, min_i, min_omega, min_Omega, min_sat_i],
            [max_h_p, max_h_a, max_i, max_omega, max_Omega, max_sat_i]
        )
    else:
        design_var_range = (
            [min_h_p, min_h_a, min_i, min_sat_i],
            [max_h_p, max_h_a, max_i, max_sat_i]
        )

    # Setup the optimisation problem
    fitness_weights = [1, 1, 1, 1] if all_obj else [1, 1]
    fitness_names = ["Mean power", "Periapsis decay", "Mean altitude", "Mean Drag/Thrust"] if all_obj else ["Mean altitude", "Periapsis decay"]

    # Select whether to run the optimisation or load the latest result file
    run_opti = ask_choice("Run the optimisation (y) or load latest results (n)? (y/n) : ", "yn", str)
    # Configure population size and number of generations
    pop_size = 60
    n_generations = 50
    seed = 12345
    # Setup file save path
    f_name = "DC_%s_%i-%s-%.2f_%i-%s_%s" % ("4O" if all_obj else "2O", thrust_model, "V" if use_battery else "X", ionisation_efficiency, \
        pop_size, seed, time.strftime("%d%m%y_%H%M%S"))
    f_path = sys.path[0] + "/optimisation/results/" + f_name
    if run_opti:        # Run a new optimisation
        current_problem = DCp.DC_problem(design_var_range, fitness_weights, thrust_model=thrust_model, \
            ionisation_efficiency=ionisation_efficiency, use_battery=use_battery, all_obj=all_obj, verbose=False)
        problem = pygmo.problem(current_problem)
        # Load the last saved generation with the same parameters
        fit_inputs, fit_results, opti_hist, last_g_file = get_saved_res(raise_error=False)
        pop = None
        if fit_inputs is not None:
            # If results from a previous run could be found...
            last_inputs, last_results = fit_inputs[-pop_size:], fit_results[-pop_size:]
            n_gen = fit_inputs.shape[0]//pop_size
            start_from_last = ask_choice("Start from the last saved population, after %i generations? (y/n) : " % (n_gen-1), "yn", str)
            if start_from_last:
                # Populate the Pygmo population with results from last run
                pop = pygmo.population(problem, size=0, seed=seed, b=pygmo.default_bfe())
                for i in range(pop_size):
                    _in = last_inputs[i]
                    if all_obj:
                        dv = [float(_in[1]), float(_in[2]), \
                            float(_in[3]), float(_in[4]), float(_in[5]), list(satellites.keys()).index(_in[0])]
                    else:
                        dv = [float(_in[1]), float(_in[2]), float(_in[3]), list(satellites.keys()).index(_in[0])]
                    _fit = last_results[i]
                    fits = [float(_fit[_i]) for _i in range(len(fitness_weights))]
                    pop.push_back(dv, fits)
                DCp.FIT_INPUTS = fit_inputs.tolist()
                DCp.FIT_RESULTS = fit_results.tolist()
                opti_hist = opti_hist.tolist()
        if pop is None:
            n_gen = 1
            last_g_file = None
            opti_hist = []
            print("Generating starting population (of size %i)..." % pop_size)
            pop = pygmo.population(problem, size=pop_size, seed=seed, b=pygmo.default_bfe())
        algo = pygmo.nsga2(seed=seed, cr=0.95, eta_c=10, m=0.005, eta_m=10)
        algo.set_bfe(pygmo.bfe())
        algo = pygmo.algorithm(algo)

        # Run the optimisation
        t0 = time.time()
        for i in range(n_gen,n_generations+1):
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
            np.savez(f_path+"-%i"%i, inputs=DCp.FIT_INPUTS, results=np.array(DCp.FIT_RESULTS), opti_hist=np.array(opti_hist))
            # Remove results from previous generation
            if last_g_file is not None:
                os.remove(last_g_file)
                last_g_file = None
            elif i > 1:
                os.remove(f_path+"-%s.npz"%(i-1))

    # Load the last saved results
    fit_inputs, fit_results, opti_hist, _ = get_saved_res()
    s_names = fit_inputs[:,0].reshape((len(fit_inputs),1))
    fit_inputs = np.array(fit_inputs[:,1:], dtype=float)
    fit_inputs = np.concatenate([s_names, fit_inputs], axis=1, dtype=object)
    print("Explored %i different possibilities." % len(fit_inputs))

    if all_obj:
        power_f, decay_f, h_f, D_T_f = fit_results[:,:4].T
        pareto_opt = PU.pareto_optimums([power_f, decay_f, h_f, D_T_f])
    else:
        h_f, decay_f, = fit_results[:,:2].T
        pareto_opt = PU.pareto_optimums([h_f, decay_f])
    print("There are %i Pareto optimum solutions across all %i objectives." % (np.sum(pareto_opt), 4 if all_obj else 2))
    print("Saving various Pareto plots...")

    plots_title = "%i objectives, with %s, %s battery" % (4 if all_obj else 2, "tank" if thrust_model == 2 else "ABE", "with" if use_battery else "without")
    if thrust_model == 3:
        plots_title += ", ionisation efficiency of %s%%" % (ionisation_efficiency*100)

    ## Plot the fitness progress over time
    PU.plot_multiple([list(range(1, opti_hist.shape[0]+1))]*(len(fitness_weights)+1), opti_hist.T, "Generation number", "Best fitness", \
        plots_path+"history", legends=fitness_names+["Average fitness"], colors=["darkorange", "seagreen", "royalblue", "#202020"], \
        lstyle=["solid"]*len(fitness_weights)+["dashed"], title=plots_title)

    ## Generate Pareto fronts
    if all_obj:
        # Filter periapsis decays above 100km
        idx_remove = np.where(fit_results[:,-3] > 100e3)
        # Remove at selected indexes
        obj_power, obj_decay, obj_h, obj_TD = np.delete(fit_results[:,-4], idx_remove), np.delete(fit_results[:,-3], idx_remove), np.delete(fit_results[:,-2], idx_remove), np.delete(fit_results[:,-1], idx_remove)
    else:
        idx_remove = np.where(fit_results[:,-1] > 100e3)
        obj_h, obj_decay = np.delete(fit_results[:,-2], idx_remove), np.delete(fit_results[:,-1], idx_remove)
    s_names = np.delete(fit_inputs[:,0], idx_remove)
    # Select color as a function of the satellite
    s_numbers = np.arange(0, len(satellites), 1)
    s_nn_map = dict(zip(list(satellites.keys()), s_numbers))
    ## Make the plots
    # Classic Pareto plots, 2 objectives
    PU.plot_single(obj_h/1e3, obj_decay/1e3, "Mean altitude [km]", "Periapsis decay [km]", plots_path+"Pareto_hd", \
        scatter=True, add_front=True, title=plots_title)
    # Plot decay vs mean altitude with satellite name in the colormap
    s_name_cmap = matplotlib.colors.ListedColormap(['red', 'green', 'blue', 'yellow', 'orange'])
    bounds = [0, 1, 2, 3, 4, 5]
    norm =  matplotlib.colors.BoundaryNorm(bounds, s_name_cmap.N)
    PU.plot_single(obj_h/1e3, obj_decay/1e3, "Mean altitude [km]", "Periapsis decay [km]", plots_path+"Pareto_hdS", \
        scatter=True, add_front=True, z_data=[s_nn_map[s_n] for s_n in s_names], z_label="Satellite", \
            cmap=s_name_cmap, cticks=[0.5, 1.5, 2.5, 3.5, 4.5], \
            clabels=list(satellites.keys()), NB=(norm, bounds), title=plots_title)
    if all_obj:
        PU.plot_single(obj_power, obj_decay/1e3, "Mean Power [W]", "Periapsis decay [km]", plots_path+"Pareto_Pd", \
            scatter=True, add_front=True, front_sign=[-1,1], title=plots_title)
        PU.plot_single(obj_power, obj_h/1e3, "Mean Power [W]", "Mean altitude [km]", plots_path+"Pareto_Ph", \
            scatter=True, add_front=True, front_sign=[-1,1], title=plots_title)
        PU.plot_single(obj_h/1e3, obj_TD, "Mean altitude [km]", "Mean Thrust/Drag [-]", plots_path+"Pareto_hT", \
            scatter=True, add_front=True, front_sign=[1,-1], title=plots_title)
        PU.plot_single(obj_power, obj_TD, "Mean Power [W]", "Mean Thrust/Drag [-]", plots_path+"Pareto_PT", \
            scatter=True, add_front=True, front_sign=[-1,-1], title=plots_title)
        PU.plot_single(obj_TD, obj_decay/1e3, "Mean Thrust/Drag [-]", "Periapsis decay [km]", plots_path+"Pareto_Td", \
            scatter=True, add_front=True, front_sign=[-1,1], title=plots_title)
        # Plot decay vs mean altitude with power in the colormap
        power_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
                        'trunc({n},{a:.2f},{b:.2f})'.format(n="plasma", a=0.0, b=0.9),
                        matplotlib.pyplot.get_cmap("plasma")(np.linspace(0.0, 0.9, 10)))
        PU.plot_single(obj_h/1e3, obj_decay/1e3, "Mean altitude [km]", "Periapsis decay [km]", plots_path+"Pareto_hdP", \
            scatter=True, add_front=True, z_data=obj_power, z_label="Mean power [W]", cmap=power_cmap, title=plots_title)
        # Plot decay vs mean altitude with D/T in the colormap
        PU.plot_single(obj_h/1e3, obj_decay/1e3, "Mean altitude [km]", "Periapsis decay [km]", plots_path+"Pareto_hdT", \
            scatter=True, add_front=True, z_data=np.clip(obj_TD, 0, 10), z_label="Mean Thrust/Drag [-]", cmap="rainbow", title=plots_title)

    idx_remove = np.where(fit_results[:,-3] >= 100e3)
    # Make a Panda dataframe from the results
    import pandas as pd
    import plotly.express as px
    s_pd = pd.Series(np.delete(fit_results[:,-3], idx_remove)/1e3, name="Periapsis decay [km]")
    s_mh = pd.Series(np.delete(fit_results[:,-2], idx_remove)/1e3, name="Mean altitude [km]")
    s_sn = pd.Series(np.delete(fit_inputs[:,0], idx_remove), name="Satellite")
    s_i_hp = pd.Series(np.delete(fit_inputs[:,1], idx_remove)/1e3, name="Initial h_p [km]")
    s_i_ha = pd.Series(np.delete(fit_inputs[:,2], idx_remove)/1e3, name="Initial h_a [km]")
    s_i_i = pd.Series(np.delete(fit_inputs[:,3], idx_remove)*180/np.pi, name="Initial i [deg]")
    s_pareto_opt = pd.Series(np.delete(pareto_opt, idx_remove), name="Pareto optimum")
    s_f_D = pd.Series(np.delete(decay_f, idx_remove), name="Periapsis decay fitness [-]")
    s_f_H = pd.Series(np.delete(h_f, idx_remove), name="Altitude fitness [-]")
    if all_obj:
        s_i_omega = pd.Series(np.delete(fit_inputs[:,4], idx_remove)*180/np.pi, name="Initial omega [deg]")
        s_i_Omega = pd.Series(np.delete(fit_inputs[:,5], idx_remove)*180/np.pi, name="Initial Omega [deg]")
        s_mTD = pd.Series(np.delete(fit_results[:,-1], idx_remove), name="Mean T/D [-]")
        s_mp = pd.Series(np.delete(fit_results[:,-4], idx_remove), name="Mean power [W]")
        s_f_P = pd.Series(np.delete(power_f, idx_remove), name="Power fitness [-]")
        s_f_TD = pd.Series(np.delete(D_T_f, idx_remove), name="T/D fitness [-]")
        df = pd.concat([s_pd, s_mh, s_mTD, s_mp, s_sn, s_i_hp, s_i_ha, s_i_i, s_i_omega, s_i_Omega, s_pareto_opt, s_f_P, s_f_D, s_f_H, s_f_TD], axis=1)
        h_data = {"Initial h_p [km]": ":.7f", "Initial h_a [km]": ":.7f", "Initial i [deg]": ":.7f", "Initial omega [deg]": ":.7f", "Initial Omega [deg]": ":.7f", \
            "Mean power [W]": ":.3f", "Periapsis decay [km]": ":.3f", "Mean altitude [km]": ":.3f", "Mean T/D [-]": ":.5f", "Power fitness [-]": ":.5f", \
            "Periapsis decay fitness [-]": ":.5f", "Altitude fitness [-]": ":.5f", "T/D fitness [-]": ":.5f"}
    else:
        df = pd.concat([s_pd, s_mh, s_sn, s_i_hp, s_i_ha, s_i_i, s_pareto_opt, s_f_D, s_f_H], axis=1)
        h_data = {"Initial h_p [km]": ":.7f", "Initial h_a [km]": ":.7f", "Initial i [deg]": ":.7f", "Periapsis decay [km]": ":.3f", \
            "Mean altitude [km]": ":.3f", "Periapsis decay fitness [-]": ":.5f", "Altitude fitness [-]": ":.5f"}
    # Make interactive plot
    fig = px.scatter(df, hover_data=h_data, x="Mean altitude [km]", y="Periapsis decay [km]", hover_name="Satellite", \
        title=plots_title, color="Pareto optimum", color_discrete_map={True: "#2e7d32", False: "#b71c1c"})
    fig.write_html(sys.path[0]+"/figures/"+plots_path+"interactive_pareto.html")
