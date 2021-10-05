import sys
sys.path.insert(0,"\\".join(sys.path[0].split("\\")[:-2]))
from tudatpy.kernel import constants
from tudatpy.kernel.simulation import propagation_setup
from utils import old_setup_utils as SU
from tools import plot_utilities as PU
from tools import time_conversions as TC
import thrust.simple_thrust as T
import numpy as np
from utils import sat_models as SM

simulation_days = 2
simulation_start_epoch = TC.MCD_to_Tudat(2459942)
simulation_end_epoch = simulation_start_epoch + simulation_days*constants.JULIAN_DAY
satellite = SM.satellites["CS_3021"]

# Define the environment and bodies
bodies, bodies_to_propagate, central_bodies = SU.create_bodies(use_MCD_atmo=False, sat=satellite)
# Define the initial state of the satellite
initial_state = SU.get_initial_state(bodies, altitude=140e3, eccentricity=0, inclination=np.deg2rad(0))
# Define the termination settings
termination_settings = SU.simulation_settings(simulation_end_epoch)
# Define the dependent variables to save
dependent_variables_to_save = [
    propagation_setup.dependent_variable.single_acceleration(
        propagation_setup.acceleration.thrust_acceleration_type, "Satellite", "Satellite"),
    propagation_setup.dependent_variable.density("Satellite", "Mars"),
    propagation_setup.dependent_variable.altitude("Satellite", "Mars"),
    propagation_setup.dependent_variable.body_mass("Satellite")
]

acceleration_settings = {"Satellite":
    dict(
        Mars=
        [
            propagation_setup.acceleration.spherical_harmonic_gravity(4, 4),
            propagation_setup.acceleration.aerodynamic()
        ],
        Sun = [ propagation_setup.acceleration.cannonball_radiation_pressure() ],
        Satellite = [  T.thrust_settings(bodies, satellite, simulation_start_epoch, save_power=True, thrust_mod=1) ]
    )}
accelerations = propagation_setup.create_acceleration_models(bodies, acceleration_settings, bodies_to_propagate, central_bodies)


# Define the mass propagator settings
mass_rate_model = {"Satellite": [propagation_setup.mass.from_thrust()]}
mass_propagator_settings = propagation_setup.propagator.mass(
    bodies_to_propagate,
    mass_rate_model,
    np.array([satellite.wet_mass]),
    termination_settings)
# Define the translational propagator settings
integrator_settings = SU.get_integrator_settings_thrust(simulation_start_epoch)
translational_propagator_settings = propagation_setup.propagator.translational(
    central_bodies,
    accelerations,
    bodies_to_propagate,
    initial_state,
    termination_settings,
    propagation_setup.propagator.encke,
    dependent_variables_to_save
)
# Define the full propagator settings
propagator_settings = propagation_setup.propagator.multitype(
    [translational_propagator_settings, mass_propagator_settings], termination_settings, dependent_variables_to_save)
# The following 3 lines are needed to properly propagate the mass
translational_propagator_settings.recreate_state_derivative_models(bodies)
translational_propagator_settings.reset_and_recreate_acceleration_models(acceleration_settings, bodies)
propagator_settings.recreate_state_derivative_models(bodies)
# Run the simulation
time, states, dep_vars = SU.run_simulation(bodies, integrator_settings, propagator_settings, return_raw=True)

# Extract values for plotting
positions = np.linalg.norm(states[:,:3], axis=1)
thrust_acc = dep_vars[:,:3]
density = dep_vars[:,3]
altitude = dep_vars[:,4]
sat_mass_hist = dep_vars[:,5]

power_vals = T.power_dict.values()
time_power = list(T.power_dict.keys())
time_power = (np.array(time_power) - time_power[0])/3600

print("Maximum power available:", max(power_vals))

print("Start/end mass of satellite:", sat_mass_hist[0], sat_mass_hist[-1])

PU.plot_single(time_power, power_vals, "Time [hr]", "Power [W]", "thrust/power_ht", scatter=True)
PU.plot_single(time/3600, altitude/1e3, "Time [hr]", "Altitude [m]", "thrust/alt_ht")
PU.plot_single(states[:,0], states[:,1], "x [m]", "y [m]", "thrust/pos_ht", scatter=True, equal_ax=True)
PU.plot_multiple([time/3600]*3, thrust_acc.T, "Time [hr]", "Thrust acceleration [m/s$^2$]", "thrust/acc_ht", ["x-direction", "y-direction", "z-direction"])
PU.plot_dual(time/3600, np.linalg.norm(thrust_acc, axis=1), sat_mass_hist, \
    "Time [hr]", "Thrust acceleration [m/s$^2$]", "Satellite mass [kg]", "thrust/mass_ht")