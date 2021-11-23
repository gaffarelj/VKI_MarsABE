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

use_tank = False
use_problem_module = False

## Propagation main parameters
# With this set of parameters, the satellite reenters after 150 days when using the ABE  (160 days with the battery, 80 days without anything)
# This config uses thrust only when the satellite starts dipping in the atmosphere, and manages to bring it back up
s_name, use_tank = "CS_2021", False
h_p, h_a = 137015.7038854, 208853.817262
i, omega, Omega = 1.42566047, 0.76761797586, 2.19917889074

satellites = SM.satellites_with_tank if use_tank else SM.satellites
thrust_model = 2 if use_tank else 3

if use_problem_module:
    from optimisation import drag_comp_problem as DCp
    power_f, decay_f, h_f, D_T_f, mean_P, decay, mean_h, mean_T_D = DCp.comp_fitness(satellites[s_name], h_p, h_a, i, omega, Omega, thrust_model)
    print("Mean power of %.2f W; decay of %.2f km; mean altitude of %.2f km; mean T/D of %.2f" % (mean_P, decay/1e3, mean_h/1e3, mean_T_D))
    print("fitness: power=%.2f; decay=%.2f; altitude=%.2f; T/D=%.2f" % (power_f, decay_f, h_f, D_T_f))
else:
    # Create the orbital simulation instance
    OS = P.orbit_simulation(satellites[s_name], "Mars", 650*constants.JULIAN_DAY, save_power=True, verbose=True)
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
    np.savetxt(sys.path[0]+"/optimisation/results/orbit.csv", results, delimiter=',', header="time, X, Y, Z")

    ## Plot results
    # Battery capacity
    battery_times = np.array(list(OS.battery_capacity.keys()))
    battery_times_hr = (battery_times-battery_times[0])/3600
    battery_capacity = list(OS.battery_capacity.values())
    PU.plot_single(battery_times_hr, battery_capacity, "Time [hr]", "Battery capacity [Whr]", "SHOW")

    # Altitude and thrust
    PU.plot_multiple([times/3600]*2, [np.array(OS.get_dep_var("h"))/1e3, np.array(OS.get_dep_var("h_p"))/1e3], "Time [hr]", "Altitude [km]", "SHOW", legends=["Satellite", "Periapsis"])
    PU.plot_single(times/3600, np.linalg.norm(OS.get_dep_var("F_T"), axis=1), "Time [hr]", "Thrust [N]", "SHOW", scatter=True)

    # Print infos
    h_p_s = OS.get_dep_var("h_p")
    print("Periapsis decay:", h_p_s[0] - h_p_s[-1])
    power_hist = list(OS.power_dict.values())
    print("Mean power:", np.mean(power_hist))