import sys
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
from tudatpy.kernel import constants
from tools import plot_utilities as PU
import numpy as np
from utils import sat_models as SM
from utils import propagation as P

# Select a satellite to propagate
satellite = SM.satellites["CS_1021"]
sim_time = 2*constants.JULIAN_DAY

# Create the orbital simulation around Mars, for two days
OS = P.orbit_simulation(satellite, "Mars", 2*constants.JULIAN_DAY, verbose=False, save_power=True)

OS.create_bodies()                                      # Create the simulation bodies
OS.create_initial_state(h_p=135e3)                        # Start in a circular orbit at an altitude of 140km above Mars
OS.create_termination_settings()                        # Create the settings to terminate the simulation
OS.create_dependent_variables(to_save=["F_T", "rho"])   # Setup the dependent variables to be saved
OS.create_accelerations(default_config=1, thrust=0)     # Create the acceleration models, with the basic thrust
OS.create_propagator()                                  # Setup the translational and mass propagator
OS.create_integrator()                                  # Create the integrator

# Run the simulation
time, states, dep_vars = OS.simulate()

# Extract values for plotting
positions = np.linalg.norm(states[:,:3], axis=1)
thrust_acc = OS.get_dep_var("F_T")
densities = OS.get_dep_var("rho")
power_vals = OS.power_dict.values()
time_power = list(OS.power_dict.keys())
time_power = (np.array(time_power) - time_power[0])/3600

# Make plot
PU.plot_single(time/3600, (positions-positions[0])/1e3, "Time [hr]", "$|r(t) - r_0|$ [km]", "thrust/test_pos")
PU.plot_dual(time/3600, np.linalg.norm(thrust_acc, axis=1), densities, "Time [hr]", "Thrust acceleration [m/s$^2$]", "Density [kg/m$^3$]", "thrust/test_acc_dens")
PU.plot_dual([time/3600, time_power/3600], np.linalg.norm(thrust_acc, axis=1), power_vals, "Time [hr]", \
    "Thrust acceleration [m/s$^2$]", "Available power [W]", "thrust/test_acc_power", diff_x=True)
PU.plot_multiple([time/3600]*3, thrust_acc.T, "Time [hr]", "Thrust acceleration [m/s$^2$]", "thrust/test_acc", ["x-direction", "y-direction", "z-direction"])