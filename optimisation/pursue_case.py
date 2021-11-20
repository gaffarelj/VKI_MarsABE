from matplotlib.pyplot import sca
import numpy as np
import sys
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
from tudatpy.kernel import constants
from utils import propagation as P
from utils import sat_models as SM
from tools import plot_utilities as PU

use_tank = True
use_problem_module = False
satellites = SM.satellites_with_tank if use_tank else SM.satellites
thrust_model = 2 if use_tank else 3

# Propagation main parameters
s_name = "CS_2021"
h_p, h_a = 115127.202, 155619.38
i, omega, Omega = 1.1571467, 0.8484065, 0.4792579

if use_problem_module:
    from optimisation import drag_comp_problem as DCp
    power_f, decay_f, h_f, D_T_f, mean_P, decay, mean_h, mean_T_D = DCp.comp_fitness(satellites[s_name], h_p, h_a, i, omega, Omega, thrust_model)
    print("Mean power of %.2f W; decay of %.2f km; mean altitude of %.2f km; mean T/D of %.2f" % (mean_P, decay/1e3, mean_h/1e3, mean_T_D))
    print("fitness: power=%.2f; decay=%.2f; altitude=%.2f; T/D=%.2f" % (power_f, decay_f, h_f, D_T_f))
else:
    # Create the orbital simulation instance, setup to simulate 10 days
    OS = P.orbit_simulation(satellites[s_name], "Mars", 10*constants.JULIAN_DAY, save_power=True, verbose=True)
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
    OS.create_dependent_variables(to_save=["h_p", "h", "D", "F_T", "m", "lat", "lon"])
    OS.create_propagator(prop_mass=use_tank)

    # Simulate the satellite in orbit
    times, states, _ = OS.simulate()

    # Save results to CSV
    results = np.array([times, states[:,0], states[:,1], states[:,2]]).T
    np.savetxt(sys.path[0]+"/optimisation/results/orbit_%s.csv" % "tank" if use_tank else "ABE", results, delimiter=',', header="time, X, Y, Z")

    # Plot results
    PU.plot_multiple([times/3600]*2, [OS.get_dep_var("h"), OS.get_dep_var("h_p")], "Time [hr]", "Altitude [km]", "SHOW", legends=["Satellite", "Periapsis"])
    PU.plot_single(times/3600, np.linalg.norm(OS.get_dep_var("F_T"), axis=1), "Time [hr]", "Thrust [N]", "SHOW", scatter=True)

    h_p_s = OS.get_dep_var("h_p")
    print("Periapsis decay:", h_p_s[0] - h_p_s[-1])
    power_hist = list(OS.power_dict.values())
    print("Mean power:", np.mean(power_hist))