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
        limit_value = 50e3, # stop at 50 km
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
    propagation_setup.propagator.unified_state_model_quaternions,
    dependent_variables_to_save
)

initial_time = 0 # seconds since J2000
initial_time_step = 20 # seconds
coefficient_set = propagation_setup.integrator.RKCoefficientSets.rkdp_87
minimum_step_size = 10 # seconds
maximum_step_size = 1800 # seconds
relative_error_tolerance = 1e-9 # -
absolute_error_tolerance = 1e-9 # -

integrator_settings = propagation_setup.integrator.runge_kutta_variable_step_size(
    initial_time,
    initial_time_step,
    coefficient_set,
    minimum_step_size,
    maximum_step_size,
    relative_error_tolerance,
    absolute_error_tolerance,
    save_frequency= 1,
    assess_termination_on_minor_steps = False,
    safety_factor = 0.95,
    maximum_factor_increase = 4.0,
    minimum_factor_increase = 0.05 )

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

rk_4_baseline = np.loadtxt(dir_path + "/rk_4_baseline.dat")
baseline_t = rk_4_baseline[0,:]
baseline_h = rk_4_baseline[1,:]

latest_time = min(baseline_t[-1], time[-1])
interp_times = np.arange(0, latest_time, 0.1)

interpolator_settings = interpolators.lagrange_interpolation(8, boundary_interpolation=interpolators.use_boundary_value)
first_interpolator = interpolators.create_one_dimensional_interpolator(dict(zip(time, altitudes)), interpolator_settings)
second_interpolator = interpolators.create_one_dimensional_interpolator(dict(zip(baseline_t, baseline_h)), interpolator_settings)
# Calculate the difference between the first and second model at specific epochs
model_difference = {t: second_interpolator.interpolate(t)- first_interpolator.interpolate(t) 
                    for t in interp_times}

diff_times = list(model_difference.keys())
diff_vals = list(model_difference.values())

plt.plot(np.array(diff_times[:-1])/24, np.array(diff_vals[:-1])/1e3)
plt.grid(), plt.xlabel("Time [days]"), plt.ylabel("Altitude [km]")
plt.show()
print("Final altitude of %.2f km, at a time of %.2f years" % (altitudes[-1]/1e3, time[-1]/365/24))
print("Max difference with rk_4 is of %.3f km, with a cpu time of %.2f seconds." % (max(diff_vals)/1e3, cpu_time))

# rk_4 from baseline (step of 10s) is around 110 seconds
# rkdp_87 from 1_year_integrators is around 10.5 seconds, 0.9 km difference

# rkdp_87 with Cowell is around 10.5 seconds, 0.9 km difference
# rkdp_87 with unified_state_model_quaternions from 1_year_propagator is around 10 seconds, 0.05 km difference

# same with 2.5e-1 tolerance: same result -> limited by max step size
# max step size from 300 to 1800: 1.7 seconds, 3.6 km deviation
# tolerance back to 2.5e-8: 3.5 seconds, 1.2 km difference. Difference especially big when perturbation becomes higher
# tolereance to 1e-9: 5 seconds, 0.6 km difference, jumps in altitude
# tolereance to 5e-9: 4.2 seconds, 0.9 km difference, jumps in altitude
# safety_factor from 0.75 to 0.85: 4.3 seconds, 0.9 km difference, much smaller jumps
# safety_factor to 0.95: 4.3 seconds, 0.9 km difference, no visible jumps
# tolereance to 1e-9: 5 seconds, 0.9 km difference