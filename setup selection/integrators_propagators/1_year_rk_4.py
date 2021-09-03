# Import modules
import numpy as np
from tudatpy.kernel import constants
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel.simulation import environment_setup
from tudatpy.kernel.simulation import propagation_setup
from tudatpy.kernel.astro import conversion
from tudatpy.kernel.math import interpolators
from matplotlib import pyplot as plt
import time
import os
dir_path = os.path.dirname(os.path.realpath(__file__))

# Setup kernel
spice_interface.load_standard_kernels()

simulation_start_epoch = 0.0
simulation_end_epoch = 1*constants.JULIAN_YEAR

# Create bodies
bodies_to_create = ["Mars"]

global_frame_orientation = "ECLIPJ2000"

body_settings = environment_setup.get_default_body_settings(
    bodies_to_create,
    base_frame_orientation=global_frame_orientation
)
density_scale_height = 11.1e3
density_at_zero_altitude = 0.020
body_settings.get("Mars").atmosphere_settings = environment_setup.atmosphere.exponential(
    density_scale_height, density_at_zero_altitude)

bodies = environment_setup.create_system_of_bodies(body_settings)

# Add satellite body
bodies.create_empty_body("Satellite")
bodies.get_body("Satellite").set_constant_mass(200)

# Add aerodynamic settings
S_ref = 10
C_d, C_l, C_m = 1.2, 0.5, 0.1
aero_coefficient_settings = environment_setup.aerodynamic_coefficients.constant(S_ref, [C_d, C_l, C_m])
environment_setup.add_aerodynamic_coefficient_interface(bodies, "Satellite", aero_coefficient_settings )

# Setup environment
bodies_to_propagate = ["Satellite"]

central_bodies = ["Mars"]

acceleration_settings = {"Satellite":
    dict(
        Mars=
        [
            propagation_setup.acceleration.point_mass_gravity(),
            propagation_setup.acceleration.aerodynamic()
        ]
    )
}

acceleration_models = propagation_setup.create_acceleration_models(
    bodies,
    acceleration_settings,
    bodies_to_propagate,
    central_bodies
)

# Define initial state
mars_gravitational_parameter = bodies.get_body("Mars").gravitational_parameter
mars_radius = spice_interface.get_average_radius("Mars")

initial_state = conversion.keplerian_to_cartesian(
    gravitational_parameter = mars_gravitational_parameter,
    semi_major_axis = mars_radius + 300e3,
    eccentricity = 0.01,
    inclination = np.deg2rad(0),
    argument_of_periapsis = np.deg2rad(0),
    longitude_of_ascending_node = np.deg2rad(0),
    true_anomaly = np.deg2rad(0)
)

# Propagator and integrator settings
dependent_variables_to_save = [
    propagation_setup.dependent_variable.altitude("Satellite", "Mars"),
    propagation_setup.dependent_variable.density("Satellite", "Mars")
]

termination_altitude = propagation_setup.dependent_variable.altitude("Satellite", "Mars")
termination_setting_altitude = propagation_setup.propagator.dependent_variable_termination(
        dependent_variable_settings = termination_altitude,
        limit_value = 10e3, # stop at 10 km
        use_as_lower_limit = True)
termination_setting_time = propagation_setup.propagator.time_termination(simulation_end_epoch)
termination_settings_list = [termination_setting_altitude, termination_setting_time]
termination_settings = propagation_setup.propagator.hybrid_termination( 
    termination_settings_list, fulfill_single_condition = True )

propagator_settings = propagation_setup.propagator.translational(
    central_bodies,
    acceleration_models,
    bodies_to_propagate,
    initial_state,
    termination_settings,
    output_variables = dependent_variables_to_save
)

fixed_step_size = 10
integrator_settings = propagation_setup.integrator.runge_kutta_4(
    simulation_start_epoch,
    fixed_step_size
)

# Run the simulation
t0 = time.time()
dynamics_simulator = propagation_setup.SingleArcDynamicsSimulator(
    bodies, integrator_settings, propagator_settings
)

states = dynamics_simulator.state_history

dependent_variables = dynamics_simulator.dependent_variable_history
cpu_time = time.time() - t0

print("Simulation took %.2f seconds." % cpu_time)

# Compute results
time = [t / 3600 for t in dependent_variables.keys()]

dependent_variable_list = np.vstack( list( dependent_variables.values( ) ) )

altitudes = dependent_variable_list[:,0]
densities = dependent_variable_list[:,1]

plt.plot(np.array(time)/24, altitudes/1e3)
plt.grid(), plt.xlabel("Time [days]"), plt.ylabel("Altitude [km]")
plt.show()

np.savetxt(dir_path + "/rk_4_baseline.dat", np.array([time, altitudes]), fmt="%.5e")
# baseline took 110 seconds