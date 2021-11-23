import numpy as np
import sys
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
from utils import sat_models as SM
from utils import propagation as P
from tudatpy.kernel import constants


def comp_fitness(sat, h_p, h_a, i, omega, Omega, thrust_model, ionisation_eff, use_battery):
    ## Setup the simulation
    # Create the orbital simulation instance, setup to simulate 100 days
    sim_days = 1
    OS = P.orbit_simulation(sat, "Mars", sim_days*constants.JULIAN_DAY, save_power=True)
    # Create the simulation bodies, and use the MCD
    OS.create_bodies(use_MCD=[False, False], use_GRAM=False)
    # Create the initial state of the satellite
    a = OS.R_cb + (h_a+h_p)/2
    e = 1 - (OS.R_cb + min(h_p, h_a)) / a       # Use min because h_p could actually be higher than h_a due to the way the problem is setup)
    OS.create_initial_state(a=a, e=e, i=i, omega=omega, Omega=Omega)
    # Load the accelerations from default config 1: Central body spherical harmonics of degree/order 4 and aerodynamics, Solar radiation
    OS.create_accelerations(default_config=1, thrust=thrust_model, ionisation_eff=ionisation_eff, use_battery=use_battery)
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

    # Remove the orbital simulation variable to save memory
    del OS

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
    D_T_f = 1 / (mean_T_D + 1)

    return power_f, decay_f, h_f, D_T_f, mean_P, decay, mean_h, mean_T_D