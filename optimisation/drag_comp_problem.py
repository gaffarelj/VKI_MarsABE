import numpy as np
import sys
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
from utils import sat_models as SM
from utils import propagation as P
from tudatpy.kernel import constants


FIT_INPUTS, FIT_HASHS, FIT_RESULTS = [], [], []

def comp_fitness(sat, h_p, h_a, i, omega, Omega, thrust_model):
    # Save the inputs in a list, and compute their hash
    fit_input = [sat.name, h_p, h_a, i, omega, Omega]
    fit_hash = hash(frozenset(fit_input))
    # Search if the fitness was already computed for these inputs
    if fit_hash in FIT_HASHS:
        idx = FIT_HASHS.index(fit_hash)
        # Double check that the inputs were the same, not just the hash
        if FIT_INPUTS[idx] == fit_input:
            # Return the cached fitness
            return FIT_RESULTS[idx]

    ## Setup the simulation
    # Create the orbital simulation instance, setup to simulate 10 days
    sim_days = 10
    OS = P.orbit_simulation(sat, "Mars", sim_days*constants.JULIAN_DAY, save_power=True)
    # Create the simulation bodies, and use the MCD
    OS.create_bodies(use_MCD=[False, False], use_GRAM=False)
    # Create the initial state of the satellite
    a = OS.R_cb + (h_a+h_p)/2
    e = 1 - (OS.R_cb + min(h_p, h_a)) / a       # Use min because h_p could actually be higher than h_a due to the way the problem is setup)
    OS.create_initial_state(a=a, e=e, i=i, omega=omega, Omega=Omega)
    # Load the accelerations from default config 1: Central body spherical harmonics of degree/order 4 and aerodynamics, Solar radiation
    OS.create_accelerations(default_config=1, thrust=thrust_model)
    # Create the integrator, termination settings, dependent variables, and propagator
    OS.create_integrator()
    OS.create_termination_settings()
    OS.create_dependent_variables(to_save=["h_p", "h", "D", "F_T"])
    prop_mass = (thrust_model != 3)
    OS.create_propagator(prop_mass=prop_mass)

    # Simulate the satellite in orbit
    times, _, _ = OS.simulate()

    # Extract the results from the simulation
    power_hist = list(OS.power_dict.values())
    h_p_s = OS.get_dep_var("h_p")
    altitudes = OS.get_dep_var("h")
    drags = OS.get_dep_var("D")
    drags_norm = np.fabs(np.linalg.norm(drags, axis=1))
    thrusts = OS.get_dep_var("F_T")
    thrusts_norm = np.fabs(np.linalg.norm(thrusts, axis=1))

    # Compute the simulation performance parameters
    mean_P, decay, mean_h, mean_T_D = np.mean(power_hist), h_p_s[0] - h_p_s[-1], np.mean(altitudes), np.mean(thrusts_norm)/np.mean(drags_norm)

    # Give a penalty in the decay if the simulation stopped earlier and the decay was not properly registered
    if times[-1] - times[0] < (sim_days-1)*constants.JULIAN_DAY and decay < 50e3:
        decay = 300e3

    ## Compute the fitness (=cost); scaling is used because, ideally, all cost values would be in the same range (0-1 for instance)
    # Max mean power; lots of power = smaller value = better (use maximum observed value as scale)
    power_scale = 35
    power_f = ( (power_scale-mean_P) / power_scale )
    if np.mean(power_hist) > power_scale:
        print("Warning, mean power of %.3f W was above scaling value of %.3f W" % (np.mean(power_hist), power_scale))
    # Min decay; if final altitude < 50km, fitness=1 (re-entered atmosphere); else: scale with maximum 100km
    decay_f = 1 if h_p_s[-1] <= 50e3 else decay/100e3
    # Min mean altitude
    h_scale = [50e3, 500e3, 1e6]
    if mean_h <= h_scale[1]:
        h_f = (mean_h - h_scale[0]) / (h_scale[1] - h_scale[0]) * 0.75
    else:
        h_f = 0.75 + 0.25 * (mean_h - h_scale[1]) / (h_scale[2] - h_scale[1])
    # Min Drag/Thrust ratio
    if np.mean(thrusts) == 0:
        D_T_f = 1
    else:
        D_T_f = 1 / mean_T_D
        D_T_f = 1 - (1/(D_T_f + 1)) # scale from [0,inf) to [0, 1]

    fit_result = power_f, decay_f, h_f, D_T_f, mean_P, decay, mean_h, mean_T_D

    # Save the inputs and results to the cache list
    FIT_INPUTS.append(fit_input)
    FIT_RESULTS.append(list(fit_result))
    FIT_HASHS.append(fit_hash)

    return fit_result

class WT_problem:

    def __init__(self, design_var_range, fitness_weights, thrust_model=1, verbose=False):
        self.design_var_range = design_var_range
        self.fitness_weights = fitness_weights
        self.verbose = verbose
        self.thrust_model = thrust_model

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

    def fitness(self, design_variables):
        """
        *** Pygmo-related function ***
        Return the fitness of the given problem. This is the cost function, that Pygmo will minimise.
        """
        # Extract the individual design variables
        h_p_0, h_a_0, i_0, omega_0, Omega_0, sat_index = design_variables

        # Select the satellite
        sats = SM.satellites if self.thrust_model == 3 else SM.satellites_with_tank
        sat_name = list(sats.keys())[int(sat_index)]
        satellite = sats[sat_name]

        # Compute the fitnesses, and the simulation performance parameters
        power_f, decay_f, h_f, D_T_f, mean_P, decay, mean_h, mean_T_D  = comp_fitness(satellite, h_p_0, h_a_0, i_0, omega_0, Omega_0, self.thrust_model)
        
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