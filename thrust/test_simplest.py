import sys
sys.path.insert(0,"\\".join(sys.path[0].split("\\")[:-2]))
from tudatpy.kernel import constants
from tudatpy.kernel.simulation import propagation_setup
from tudatpy.kernel.interface import spice_interface
from setup_selection import setup_utils as SU
from tools import plot_utilities as PU
from tools import time_conversions as TC
import thrust.simple_thrust as T
import numpy as np

simulation_days = 2
simulation_start_epoch = TC.MCD_to_Tudat(2459942)
simulation_end_epoch = simulation_start_epoch + simulation_days*constants.JULIAN_DAY

# Define the environment and bodies
bodies, bodies_to_propagate, central_bodies = SU.create_bodies(use_MCD_atmo=True, sat_name="CS_1021")
# Define the initial state of the satellite
initial_state = SU.get_initial_state(bodies, altitude=140e3, eccentricity=0)
# Define the termination settings
termination_settings = SU.simulation_settings(simulation_end_epoch)
# Define the dependent variables to save
dependent_variables_to_save = [
    propagation_setup.dependent_variable.single_acceleration(
        propagation_setup.acceleration.thrust_acceleration_type, "Satellite", "Satellite"),
    propagation_setup.dependent_variable.density("Satellite", "Mars"),
    propagation_setup.dependent_variable.angle_of_attack("Satellite", "Mars")
]

acceleration_settings = {"Satellite":
    dict(
        Mars=
        [
            propagation_setup.acceleration.spherical_harmonic_gravity(4, 4),
            propagation_setup.acceleration.aerodynamic()
        ],
        Sun = [ propagation_setup.acceleration.cannonball_radiation_pressure() ],
        Satellite = [  T.thrust_settings(bodies, simulation_start_epoch, save_power=True) ]
    )}
accelerations = propagation_setup.create_acceleration_models(bodies, acceleration_settings, bodies_to_propagate, central_bodies)


# Define the propagator settings
integrator_settings = SU.get_integrator_settings_thrust(simulation_start_epoch)
propagator_settings = propagation_setup.propagator.translational(
    central_bodies,
    accelerations,
    bodies_to_propagate,
    initial_state,
    termination_settings,
    propagation_setup.propagator.encke,
    dependent_variables_to_save
)

# Run the simulation
time, states, dep_vars = SU.run_simulation(bodies, integrator_settings, propagator_settings, return_raw=True)

# Print the minimum and maximum time step
# dt = time[1:] - time[:-1]
# print(max(dt), min(dt))

# Extract values for plotting
positions = np.linalg.norm(states[:,:3], axis=1)
thrust_acc = dep_vars[:,:3]
densities = dep_vars[:,3]
aoa = dep_vars[:,4]
time_power = np.array(list(T.solar_irradiances.keys()))
time_power = time_power - time_power[0]

# Make plot
PU.plot_single(time/3600, (positions-positions[0])/1e3, "Time [hr]", "$|r(t) - r_0|$ [km]", "thrust/test_pos")
PU.plot_dual(time/3600, np.linalg.norm(thrust_acc, axis=1), densities, "Time [hr]", "Thrust acceleration [m/s$^2$]", "Density [kg/m$^3$]", "thrust/test_acc_dens")
PU.plot_dual([time/3600, time_power/3600], np.linalg.norm(thrust_acc, axis=1), T.solar_irradiances.values(), "Time [hr]", \
    "Thrust acceleration [m/s$^2$]", "Solar irradiance [W/m$^2$]", "thrust/test_acc_solar", diff_x=True)
PU.plot_dual([time/3600, time_power/3600], np.linalg.norm(thrust_acc, axis=1), T.power_dict.values(), "Time [hr]", \
    "Thrust acceleration [m/s$^2$]", "Available power [W]", "thrust/test_acc_power", diff_x=True)
PU.plot_multiple([time/3600]*3, thrust_acc.T, "Time [hr]", "Thrust acceleration [m/s$^2$]", "thrust/test_acc", ["x-direction", "y-direction", "z-direction"])