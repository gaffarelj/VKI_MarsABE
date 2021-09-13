# Import modules
import numpy as np
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel.simulation import environment_setup
from tudatpy.kernel.simulation import propagation_setup
from tudatpy.kernel.astro import conversion
from tudatpy.kernel.math import interpolators
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
import time
import sys
sys.path.insert(0,"\\".join(sys.path[0].split("\\")[:-2])) # get back to uppermost level of the project
spice_interface.load_standard_kernels()

def create_bodies(use_MCD_atmo=False, use_MCD_winds=False):
    # Create bodies
    bodies_to_create = ["Mars", "Sun"]
    global_frame_origin = 'Mars'

    global_frame_orientation = "ECLIPJ2000"

    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create,
        global_frame_origin,
        base_frame_orientation=global_frame_orientation
    )

    # Select the atmospheric model
    if not use_MCD_atmo:
        # Exponential parameters taken from http://link.springer.com/content/pdf/10.1007%2F978-3-540-73647-9_3.pdf
        density_scale_height = 7.295e3
        density_at_zero_altitude = 0.0525
        body_settings.get("Mars").atmosphere_settings = environment_setup.atmosphere.exponential(
            density_scale_height, density_at_zero_altitude)
        if use_MCD_winds:
            print("Warning: MCD winds cannot be used sperately from the MCD atmosphere.")
    else:
        from MCD.parallel_mcd import parallel_mcd as PMCD
        mcd = PMCD()
        body_settings.get("Mars").atmosphere_settings = environment_setup.atmosphere.custom_constant_temperature_detailed(
            mcd.density, constant_temperature=210, specific_gas_constant=192, ratio_of_specific_heats=1.3)
        # Values taken from https://meteor.geol.iastate.edu/classes/mt452/Class_Discussion/Mars-physical_and_orbital_statistics.pdf
        # TODO: CHECK THAT THE CONSTANT VALUES ABOVE DO NOT INFLUENCE DRAG RESULTS

        # If specified, add winds from the MCD. Only possible if the MCD atmospheric model is used
        if use_MCD_winds and use_MCD_atmo:
            body_settings.get("Mars").atmosphere_settings.wind_settings = environment_setup.atmosphere.custom_wind_model(mcd.wind)

    bodies = environment_setup.create_system_of_bodies(body_settings)

    # Add satellite body
    bodies.create_empty_body("Satellite")
    bodies.get_body("Satellite").set_constant_mass(200)

    # Add aerodynamic settings
    S_ref = 10
    C_d, C_l, C_m = 1.2, 0.5, 0.1
    aero_coefficient_settings = environment_setup.aerodynamic_coefficients.constant(S_ref, [C_d, C_l, C_m])
    environment_setup.add_aerodynamic_coefficient_interface(bodies, "Satellite", aero_coefficient_settings )

    # Add solar radiation settings
    reference_area_radiation = 4.0
    radiation_pressure_coefficient = 1.2
    occulting_bodies = ["Mars"]
    radiation_pressure_settings = environment_setup.radiation_pressure.cannonball(
    "Sun", reference_area_radiation, radiation_pressure_coefficient, occulting_bodies)
    environment_setup.add_radiation_pressure_interface(bodies, "Satellite", radiation_pressure_settings)

    bodies_to_propagate = ["Satellite"]
    central_bodies = ["Mars"]
    return bodies, bodies_to_propagate, central_bodies

def setup_environment(bodies, bodies_to_propagate, central_bodies):
    # Setup environment
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
    return acceleration_models

def get_initial_state(bodies, altitude=300e3, inclination=np.deg2rad(0)):
    # Define initial state
    mars_gravitational_parameter = bodies.get_body("Mars").gravitational_parameter
    mars_radius = spice_interface.get_average_radius("Mars")

    initial_state = conversion.keplerian_to_cartesian(
        gravitational_parameter = mars_gravitational_parameter,
        semi_major_axis = mars_radius + altitude,
        eccentricity = 0.01,
        inclination = inclination,
        argument_of_periapsis = np.deg2rad(0),
        longitude_of_ascending_node = np.deg2rad(45),
        true_anomaly = np.deg2rad(0)
    )
    return initial_state

def simulation_settings(end_time, end_altitude=50e3):
    # Propagator and integrator settings
    dependent_variables_to_save = [
        propagation_setup.dependent_variable.altitude("Satellite", "Mars"),
        propagation_setup.dependent_variable.density("Satellite", "Mars")
    ]

    termination_altitude = propagation_setup.dependent_variable.altitude("Satellite", "Mars")
    termination_setting_altitude = propagation_setup.propagator.dependent_variable_termination(
            dependent_variable_settings = termination_altitude,
            limit_value = end_altitude, # stop at 50 km
            use_as_lower_limit = True)
    termination_setting_time = propagation_setup.propagator.time_termination(end_time)
    termination_settings_list = [termination_setting_altitude, termination_setting_time]
    termination_settings = propagation_setup.propagator.hybrid_termination( 
        termination_settings_list, fulfill_single_condition = True )
    return termination_settings

def run_simulation(bodies, integrator_settings, propagator_settings, verbose=False, return_raw=False):
    print("Starting simulation...")
    t0 = time.time()
    dynamics_simulator = propagation_setup.SingleArcDynamicsSimulator(
        bodies, integrator_settings, propagator_settings, print_dependent_variable_data = verbose
    )

    states = dynamics_simulator.state_history
    states_elements = np.vstack(list(states.values()))

    dependent_variables = dynamics_simulator.dependent_variable_history
    cpu_time = time.time() - t0

    print("Simulation took %.2f seconds." % cpu_time)

    # Compute results
    #time_l = [t / 3600 for t in dependent_variables.keys()]
    time_l = np.array(list(dependent_variables.keys()))
    time_l -= time_l[0]

    dependent_variable_list = np.vstack( list( dependent_variables.values( ) ) )

    altitudes = dependent_variable_list[:,0]
    densities = dependent_variable_list[:,1]

    if return_raw:
        return np.array(time_l), np.array(states_elements), np.array(dependent_variable_list)

    return time_l, altitudes, densities, cpu_time

def compare_to_baseline(time_l, altitudes, baseline_f="rk_4_baseline", trunc_ends=False):
    rk_4_baseline = np.loadtxt("setup_selection/integrators_propagators/%s.dat" % baseline_f)
    baseline_t = rk_4_baseline[0,:]
    baseline_h = rk_4_baseline[1,:]

    earliest_time = baseline_t[0]
    latest_time = min(baseline_t[-1], time_l[-1])
    interp_times = np.arange(earliest_time, latest_time, 0.1)
    if trunc_ends:
        trunc_length = max(int(len(interp_times)*0.025), 3)
        interp_times = interp_times[trunc_length:-trunc_length]

    interpolator_settings = interpolators.lagrange_interpolation(8, boundary_interpolation=interpolators.use_boundary_value)
    first_interpolator = interpolators.create_one_dimensional_interpolator(dict(zip(time_l, altitudes)), interpolator_settings)
    second_interpolator = interpolators.create_one_dimensional_interpolator(dict(zip(baseline_t, baseline_h)), interpolator_settings)
    # Calculate the difference between the first and second model at specific epochs
    model_difference = {t: second_interpolator.interpolate(t)- first_interpolator.interpolate(t) 
                        for t in interp_times}

    diff_times = np.array(list(model_difference.keys()))-earliest_time
    diff_vals = list(model_difference.values())
    return diff_times, diff_vals

def plot_altitude(times, altitudes):
    plt.plot(np.array(time)/24, altitudes/1e3)
    plt.grid(), plt.xlabel("Time [days]"), plt.ylabel("Altitude [km]")
    plt.show()

def plot_difference(diff_times, diff_vals):
    plt.plot(np.array(diff_times[:-1])/24, np.array(diff_vals[:-1])/1e3)
    plt.grid(), plt.xlabel("Time [days]"), plt.ylabel("Altitude [km]")
    plt.show()

def get_integrator_settings(settings_index=10, verbose=False):
    """
    Return integrators settings corresponding to a given index.
    Input:
     * settings_index: int in range 0-17:
        0: RK4 intergator, step size of 10s (baseline)
        1: RK4 intergator, step size of 30s
        2: RK4 intergator, step size of 60s
        3: RKF45 intergator, step size of 1-300s, tolerance of 1E-9
        4: RKF45 intergator, step size of 1-300s, tolerance of 1E-8
        5: RKF45 intergator, step size of 1-300s, tolerance of 1E-6
        6: RKF45 intergator, step size of 1-500s, tolerance of 1E-9
        7: RKF56 intergator, step size of 10-300s, tolerance of 1E-9
        8: RKF78 intergator, step size of 10-300s, tolerance of 1E-9
        9: RKDP87 intergator, step size of 10-300s, tolerance of 1E-9
        10: RKDP87 intergator, step size of 10-300s, tolerance of 2.5E-8
        11: RKDP87 intergator, step size of 10-300s, tolerance of 1E-8
        12: RKDP87 intergator, step size of 10-300s, tolerance of 1E-7
        13: ABM intergator, step size of 10-300s, tolerance of 1E-9, order of 6-11
        14: ABM intergator, step size of 30-300s, tolerance of 1E-9, order of 6-11
        15: ABM intergator, step size of 10-500s, tolerance of 1E-9, order of 6-11
        16: BS intergator, step size of 10-500s, tolerance of 1E-9, max 5 steps
        17: BS intergator, step size of 10-500s, tolerance of 1E-9, max 4 steps
    """
    initial_time = 0
    minimum_step_size = 10
    maximum_step_size = 300
    tolerance = 1e-9
    initial_time_step = 20
    if settings_index in range(0, 3):
        step_sizes = [10, 30, 60]
        integrator_settings = propagation_setup.integrator.runge_kutta_4(
            initial_time,
            step_sizes[settings_index],
            save_frequency = 1,
            assess_termination_on_minor_steps = False
        )
        if verbose: print(settings_index, "RK4", step_sizes[settings_index])
    elif settings_index in range(3, 13):
        if settings_index in range(3, 7):
            coefficient_set = propagation_setup.integrator.RKCoefficientSets.rkf_45
            minimum_step_size = 1
            if settings_index == 4:
                tolerance = 1e-8
            elif settings_index == 5:
                tolerance = 1e-6
            elif settings_index == 6:
                maximum_step_size = 500
        elif settings_index == 7:
            coefficient_set = propagation_setup.integrator.RKCoefficientSets.rkf_56
        elif settings_index == 8:
            coefficient_set = propagation_setup.integrator.RKCoefficientSets.rkf_78
        elif settings_index in range(9, 13):
            coefficient_set = propagation_setup.integrator.RKCoefficientSets.rkdp_87
            if settings_index == 10:
                tolerance = 2.5e-8
            elif settings_index == 11:
                tolerance = 1e-8
            elif settings_index == 12:
                tolerance = 1e-7
        integrator_settings = propagation_setup.integrator.runge_kutta_variable_step_size(
	        initial_time,
	        initial_time_step,
	        coefficient_set,
	        minimum_step_size,
	        maximum_step_size,
	        tolerance,
	        tolerance,
            save_frequency= 1,
            assess_termination_on_minor_steps = False,
            safety_factor = 0.8,
            maximum_factor_increase = 4.0,
            minimum_factor_increase = 0.1 )
        if verbose: print(settings_index, coefficient_set, minimum_step_size, maximum_step_size, tolerance)
    elif settings_index in range(13, 16):
        minimum_order = 6
        maximum_order = 11
        if settings_index == 14:
            minimum_step_size = 30
            initial_time_step = 60
        elif settings_index == 15:
            maximum_step_size = 500
        integrator_settings = propagation_setup.integrator.adams_bashforth_moulton(
	        initial_time,
	        initial_time_step,
	        minimum_step_size,
	        maximum_step_size,
	        tolerance,
	        tolerance,
	        minimum_order,
	        maximum_order
        )
        if verbose: print(settings_index, "ABM", minimum_step_size, maximum_step_size, tolerance, minimum_order, maximum_order)
    elif settings_index in range(16, 18):
        extrapolation_sequence = propagation_setup.integrator.ExtrapolationMethodStepSequences.bulirsch_stoer_sequence
        maximum_number_of_steps = 5
        if settings_index == 17:
            maximum_number_of_steps = 4
        maximum_step_size = 500
        integrator_settings = propagation_setup.integrator.bulirsch_stoer(
	        initial_time,
	        initial_time_step,
	        extrapolation_sequence,
	        maximum_number_of_steps,
	        minimum_step_size,
	        maximum_step_size,
	        tolerance,
	        tolerance
        )
        if verbose: print(settings_index, extrapolation_sequence, minimum_step_size, maximum_step_size, tolerance, maximum_number_of_steps)

    return integrator_settings

def get_best_integrator(simulation_start_epoch, extra_accurate=False):
    tolerance = 1e-9 if extra_accurate else 7.5e-7
    # Setup the optimal integrator settings
    initial_time = simulation_start_epoch # seconds since J2000
    initial_time_step = 20 # seconds
    coefficient_set = propagation_setup.integrator.RKCoefficientSets.rkdp_87
    minimum_step_size = 1 if extra_accurate else 10 # seconds
    maximum_step_size = 1800 # seconds
    relative_error_tolerance = tolerance # -
    absolute_error_tolerance = tolerance # -
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
        safety_factor = 0.65,
        maximum_factor_increase = 2.0,
        minimum_factor_increase = 0.1 )
    return integrator_settings