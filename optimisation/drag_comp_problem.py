import numpy as np
import sys
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
from utils import sat_models as SM
from utils import propagation as P
from optimisation import comp_fitness as CF
import multiprocessing as MP

FIT_INPUTS, FIT_RESULTS = [], []

# Drag Compensation problem
class DC_problem:

    def __init__(self, design_var_range, fitness_weights, thrust_model=1, ionisation_efficiency=1, use_battery=True, verbose=False, all_obj=True):
        self.design_var_range = design_var_range
        self.fitness_weights = fitness_weights
        self.verbose = verbose
        self.thrust_model = thrust_model
        self.ionisation_eff = ionisation_efficiency
        self.use_battery = use_battery
        self.all_obj = all_obj

    def get_bounds(self):
        """
        *** Pygmo-related function ***
        Return the bounds of the design variables
        """
        return self.design_var_range

    def get_nobj(self):
        """
        *** Pygmo-related function ***
        Return the number of objectives
        """
        return len(self.fitness_weights)

    def get_nix(self):
        return 1

    def batch_fitness(self, dvs):
        """
        *** Pygmo-related function ***
        Return the batch fitness of the given 1D list of Design Variables. This is the cost function, that Pygmo will minimise.
        """
        global FIT_INPUTS, FIT_RESULTS
        inputs, fitnesses = [], []
        # Reshape the design variables and go trough them
        n_dv = 6 if self.all_obj else 4
        design_variables = np.reshape(dvs, (len(dvs)//n_dv, n_dv))
        sats = SM.satellites if self.thrust_model == 3 else SM.satellites_with_tank
        for dv in design_variables:
            # Extract from design variable
            if self.all_obj:
                h_p_0, h_a_0, i_0, omega_0, Omega_0, sat_index = dv
            else:
                h_p_0, h_a_0, i_0, sat_index = dv
                omega_0, Omega_0 = 0, 0

            # Select the satellite
            sat_name = list(sats.keys())[int(sat_index)]
            satellite = sats[sat_name]

            # Construct the input
            inputs.append([satellite, h_p_0, h_a_0, i_0, omega_0, Omega_0, self.thrust_model, self.ionisation_eff, self.use_battery])
            FIT_INPUTS.append([satellite.name, h_p_0, h_a_0, i_0, omega_0, Omega_0])
        
        # Get the fitness by running the orbital simulations in parallel (use half the number of processors available)
        with MP.get_context("spawn").Pool(processes=int(MP.cpu_count()-4)) as pool:
            outputs = pool.starmap(CF.comp_fitness, inputs)

        # Save the entire output and return the 1D list of fitnesses
        for output in outputs:
            FIT_RESULTS.append(list(output))
            if self.all_obj:
                power_f, decay_f, h_f, D_T_f, *_ = output
                fitnesses.append(power_f), fitnesses.append(decay_f), fitnesses.append(h_f), fitnesses.append(D_T_f)
            else:
                h_f, decay_f, *_ = output
                fitnesses.append(h_f), fitnesses.append(decay_f)
        
        return fitnesses

    def fitness(self, design_variables):
        """
        NOT IN USE ANYMORE
        *** Pygmo-related function ***
        Return the fitness of the given problem. This is the cost function, that Pygmo will minimise.
        """
        print(0, design_variables), input()
        # Extract the individual design variables
        h_p_0, h_a_0, i_0, omega_0, Omega_0, sat_index = design_variables

        # Select the satellite
        sats = SM.satellites if self.thrust_model == 3 else SM.satellites_with_tank
        sat_name = list(sats.keys())[int(sat_index)]
        satellite = sats[sat_name]

        # Compute the fitnesses, and the simulation performance parameters
        power_f, decay_f, h_f, D_T_f, mean_P, decay, mean_h, mean_T_D  = \
            CF.comp_fitness(satellite, h_p_0, h_a_0, i_0, omega_0, Omega_0, self.thrust_model, self.ionisation_eff, self.use_battery)
        
        if self.verbose:
            print("Satellite %s starts from h_p=%3d, h_a=%.2f, i=%2d, omega=%3d, Omega=%.3d" % \
                (sat_name, min(h_p_0, h_a_0)/1e3, max(h_p_0, h_a_0)/1e3, np.rad2deg(i_0), np.rad2deg(omega_0), np.rad2deg(Omega_0)))
            print(" -> mean power=%.2f W, total decay=%4d km, mean altitude=%3d km, mean T/D=%.2f" % \
                (mean_P, decay/1e3, mean_h/1e3, mean_T_D))

        # Assemble and return the cost
        cost = np.array(self.fitness_weights) * np.array([power_f, decay_f, h_f, D_T_f]) # Mean power, periapsis decay, mean altitude, mean Drag/Thrust
        if self.verbose:
            print(" -> cost is", cost)
        return cost

test_fitness_comp = False
if test_fitness_comp:
    from tools import plot_utilities as PU
    # Plot mean power vs power fitness
    mean_Ps = np.arange(0, 35, 0.1)
    power_fs = []
    for mean_P in mean_Ps:
        power_scale = 35
        power_f = ( (power_scale-mean_P) / power_scale )
        power_fs.append(power_f)
    PU.plot_single(mean_Ps, power_fs, "Mean power [W]", "Power fitness [-]", "optimisation/power_scale")
    # Plot decay vs decay fitness
    decay = np.arange(-150e3, 150e3, 10)
    PU.plot_single(decay/1e3, decay/100e3, "Decay [km]", "Decay fitness [-]", "optimisation/decay_scale")
    # Plot mean altitude vs altitude fitness
    mean_hs = np.arange(50e3, 650e3, 10)
    h_fs = []
    for mean_h in mean_hs:
        h_scale = [50e3, 500e3, 1e6]
        if mean_h <= h_scale[1]:
            h_f = (mean_h - h_scale[0]) / (h_scale[1] - h_scale[0]) * 0.75
        else:
            h_f = 0.75 + 0.25 * (mean_h - h_scale[1]) / (h_scale[2] - h_scale[1])
        h_fs.append(h_f)
    PU.plot_single(mean_hs/1e3, h_fs, "Mean altitude [km]", "Altitude fitness [-]", "optimisation/altitude_scale")
    # Plot mean T/D vs T/D fitness
    mean_T_Ds = np.arange(0, 50, 0.01)
    D_T_fs = 1 / (mean_T_Ds + 1)
    PU.plot_single(mean_T_Ds, D_T_fs, "Mean T/D [-]", "T/D fitness [-]", "optimisation/TD_scale")
