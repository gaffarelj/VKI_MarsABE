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

# Propagation main parameters
s_name = "CS_3021" # diff with CS 2120 ?
h_p, h_a = 140.93e3, 235.41e3
i, omega, Omega = np.deg2rad(81.90), np.rad2deg(173.18), np.deg2rad(166.92)
thrust_model = 2

# Create the orbital simulation instance, setup to simulate 5 days
OS = P.orbit_simulation(SM.satellites_with_tank[s_name], "Mars", 5*constants.JULIAN_DAY, save_power=True, verbose=True)
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
OS.create_propagator(prop_mass=True)

# Simulate the satellite in orbit
times, states, _ = OS.simulate()

# Save results to CSV
results = np.array([times, states[:,0], states[:,1], states[:,2]]).T
np.savetxt("orbit.csv", results, delimiter=',', header="time, X, Y, Z")

# Plot results
PU.plot_multiple([times/3600]*2, [OS.get_dep_var("h"), OS.get_dep_var("h_p")], "Time [hr]", "Altitude [km]", "SHOW", legends=["Satellite", "Periapsis"])
PU.plot_single(times/3600, np.linalg.norm(OS.get_dep_var("F_T"), axis=1), "Time [hr]", "Thrust [N]", "SHOW", scatter=True)

h_p_s = OS.get_dep_var("h_p")
print("Periapsis decay:", h_p_s[0] - h_p_s[-1])
power_hist = list(OS.power_dict.values())
print("Mean power:", np.mean(power_hist))