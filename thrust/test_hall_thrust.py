import sys
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
from tudatpy.kernel import constants
from tools import plot_utilities as PU
import numpy as np
from utils import sat_models as SM
from utils import propagation as P

# Select a satellite to propagate
satellite = SM.satellites["CS_3021"]
sim_time = 2*constants.JULIAN_DAY

# Create the orbital simulation around Mars, for two days
OS = P.orbit_simulation(satellite, "Mars", 2*constants.JULIAN_DAY, verbose=False, save_power=True)

OS.create_bodies()                                              # Create the simulation bodies
OS.create_initial_state(h_p=140e3)                                # Start in a circular orbit at an altitude of 140km above Mars
OS.create_termination_settings()                                # Create the settings to terminate the simulation
OS.create_dependent_variables(to_save=["F_T", "rho", "h", "m"]) # Setup the dependent variables to be saved
OS.create_accelerations(default_config=1, thrust=1)             # Create the acceleration models, with the hall thrust
OS.create_propagator(prop_mass=True)                            # Setup the translational and mass propagator
OS.create_integrator()                                          # Create the integrator

# Run the simulation
time, states, dep_vars = OS.simulate()
# Print some states and times (to check consistency when modifying the code)
print([states[0], states[100], states[-1]])
print([time[0], time[100], time[-1]])

# Extract results values
positions = np.linalg.norm(states[:,:3], axis=1)
thrust_acc = OS.get_dep_var("F_T")
density = OS.get_dep_var("rho")
altitude = OS.get_dep_var("h")
sat_mass_hist = OS.get_dep_var("m")

power_vals = OS.power_dict.values()
time_power = list(OS.power_dict.keys())
time_power = (np.array(time_power) - time_power[0])/3600

print("Maximum power available:", max(power_vals))
print("Start/end mass of satellite:", sat_mass_hist[0][0], sat_mass_hist[-1][0])

# Make plots
PU.plot_single(time_power, power_vals, "Time [hr]", "Power [W]", "thrust/power_ht", scatter=True)
PU.plot_single(time/3600, altitude/1e3, "Time [hr]", "Altitude [m]", "thrust/alt_ht")
PU.plot_single(states[:,0], states[:,1], "x [m]", "y [m]", "thrust/pos_ht", scatter=True, equal_ax=True)
PU.plot_multiple([time/3600]*3, thrust_acc.T, "Time [hr]", "Thrust acceleration [m/s$^2$]", "thrust/acc_ht", ["x-direction", "y-direction", "z-direction"])
PU.plot_dual(time/3600, np.linalg.norm(thrust_acc, axis=1), sat_mass_hist, \
    "Time [hr]", "Thrust acceleration [m/s$^2$]", "Satellite mass [kg]", "thrust/mass_ht")