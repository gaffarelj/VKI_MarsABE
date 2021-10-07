import sys
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
from tudatpy.kernel import constants
from tools import plot_utilities as PU
import numpy as np
from utils import sat_models as SM

from utils import propagation as P

satellite = SM.satellites["CS_3021"]
sim_time = 2*constants.JULIAN_DAY

OS = P.orbit_simulation(satellite, "Mars", 2*constants.JULIAN_DAY, verbose=False, save_power=True)

OS.create_bodies()
OS.create_initial_state(h=140e3)
OS.create_termination_settings()
OS.create_dependent_variables(to_save=["F_T", "rho", "h", "m"])
OS.create_accelerations(default_config=1, thrust=1)
OS.create_propagator(prop_mass=True)
OS.create_integrator()

bodies_to_propagate = [OS.sat.name]
central_bodies = [OS.central_body]
bodies = OS.bodies
initial_state = OS.initial_state
termination_settings = OS.termination_settings
dependent_variables_to_save = OS.dependent_variables_to_save
acceleration_settings = {OS.sat.name: OS.accelerations}
accelerations = OS.acceleration_models
propagator_settings = OS.propagator_settings
integrator_settings = OS.integrator_settings


# Run the simulation
time, states, dep_vars = OS.simulate()
print([states[0], states[100], states[-1]])
print([time[0], time[100], time[-1]])
print(dep_vars.shape)

# Extract values for plotting
positions = np.linalg.norm(states[:,:3], axis=1)
thrust_acc = OS.get_dep_var("F_T")
density = OS.get_dep_var("rho")
altitude = OS.get_dep_var("h")
sat_mass_hist = OS.get_dep_var("m")

power_vals = OS.power_dict.values()
time_power = list(OS.power_dict.keys())
time_power = (np.array(time_power) - time_power[0])/3600

print("Maximum power available:", max(power_vals))

print("Start/end mass of satellite:", sat_mass_hist[0], sat_mass_hist[-1])

PU.plot_single(time_power, power_vals, "Time [hr]", "Power [W]", "thrust/power_ht", scatter=True)
PU.plot_single(time/3600, altitude/1e3, "Time [hr]", "Altitude [m]", "thrust/alt_ht")
PU.plot_single(states[:,0], states[:,1], "x [m]", "y [m]", "thrust/pos_ht", scatter=True, equal_ax=True)
PU.plot_multiple([time/3600]*3, thrust_acc.T, "Time [hr]", "Thrust acceleration [m/s$^2$]", "thrust/acc_ht", ["x-direction", "y-direction", "z-direction"])
PU.plot_dual(time/3600, np.linalg.norm(thrust_acc, axis=1), sat_mass_hist, \
    "Time [hr]", "Thrust acceleration [m/s$^2$]", "Satellite mass [kg]", "thrust/mass_ht")