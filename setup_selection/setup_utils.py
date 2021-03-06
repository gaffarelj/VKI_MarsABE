# Import modules
import numpy as np
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel.numerical_simulation import environment_setup
from tudatpy.kernel.numerical_simulation import propagation_setup
from tudatpy.kernel.numerical_simulation import SingleArcSimulator
from tudatpy.kernel.astro import element_conversion
from tudatpy.kernel.math import interpolators
import time
import sys
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
spice_interface.load_standard_kernels()

# Drag coefficient of different satellite at given altitudes
sat_drags = {
    "CS_0020": {85e3: 3.08541, 115e3: 2.46291, 150e3: 2.34536},
    "CS_0021": {85e3: 4.49286, 115e3: 4.35434, 150e3: 4.35666},
    "CS_1020": {85e3: 4.05514, 115e3: 2.73544, 150e3: 2.67014},
    "CS_1021": {85e3: 5.56548, 115e3: 4.57014, 150e3: 4.59003},
    "CS_2020": {85e3: 3.69607, 115e3: 2.97278, 150e3: 2.93880},
    "CS_2021": {85e3: 5.14644, 115e3: 4.86814, 150e3: 4.92966},
    "CS_2120": {85e3: 5.33483, 115e3: 3.29540, 150e3: 3.22048},
    "CS_3020": {85e3: 4.27504, 115e3: 3.28570, 150e3: 3.12841},
    "CS_3021": {85e3: 5.43094, 115e3: 5.158049, 150e3: 5.19109},
}

def create_bodies(use_MCD_atmo=False, use_MCD_winds=False, sat_name="", sat_mass=200):
    # Create bodies
    bodies_to_create = ["Mars", "Sun", "Jupiter"]
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
        body_settings.get("Mars").atmosphere_settings = environment_setup.atmosphere.custom_four_dimensional_constant_temperature(
            mcd.density, constant_temperature=210, specific_gas_constant=192, ratio_of_specific_heats=1.3)
        # Values taken from https://meteor.geol.iastate.edu/classes/mt452/Class_Discussion/Mars-physical_and_orbital_statistics.pdf

        # If specified, add winds from the MCD. Only possible if the MCD atmospheric model is used
        if use_MCD_winds and use_MCD_atmo:
            body_settings.get("Mars").atmosphere_settings.wind_settings = environment_setup.atmosphere.custom_wind_model(mcd.wind)

    bodies = environment_setup.create_system_of_bodies(body_settings)

    # Add satellite body
    bodies.create_empty_body("Satellite")
    bodies.get_body("Satellite").set_constant_mass(sat_mass)

    # Add aerodynamic settings
    S_ref = 0.01
    if sat_name != "":
        # Create a cubic spline interpolator, capped at the boundaries
        interpolator_settings = interpolators.cubic_spline_interpolation(boundary_interpolation=interpolators.use_boundary_value)
        # Define the drag coefficient values at given altitudes
        drag_values = sat_drags[sat_name]
        # Setup the drag interpolator
        drag_interpolator = interpolators.create_one_dimensional_interpolator(drag_values, interpolator_settings)
        def force_coefficients(_):
            # Get the altitude from the flight conditions
            h = bodies.get_body("Satellite").flight_conditions.altitude
            # Interpolate the drag coefficient given the altitude
            C_d = drag_interpolator.interpolate(h)
            return [C_d, 0, 0]
        aero_coefficient_settings = environment_setup.aerodynamic_coefficients.custom(force_coefficients, S_ref, [])
    else:
        C_d, C_l, C_m = 2.5, 0, 0
        aero_coefficient_settings = environment_setup.aerodynamic_coefficients.constant(S_ref, [C_d, C_l, C_m])

    environment_setup.add_aerodynamic_coefficient_interface(bodies, "Satellite", aero_coefficient_settings)

    # Add solar radiation settings
    reference_area_radiation = 0.125
    radiation_pressure_coefficient = 1.2
    occulting_bodies = ["Mars"]
    radiation_pressure_settings = environment_setup.radiation_pressure.cannonball(
    "Sun", reference_area_radiation, radiation_pressure_coefficient, occulting_bodies)
    environment_setup.add_radiation_pressure_interface(bodies, "Satellite", radiation_pressure_settings)

    bodies_to_propagate = ["Satellite"]
    central_bodies = ["Mars"]
    return bodies, bodies_to_propagate, central_bodies

def setup_environment(bodies, bodies_to_propagate, central_bodies, detail_level=0):
    # Detail level: 0 = PM, aero; 1 = SH D/O 4, aero; 2 = SG D/O 8, aero, canonnball radiation
    # Setup environment
    if detail_level == 0:
        acceleration_settings = {"Satellite":
            dict(
                Mars=
                [
                    propagation_setup.acceleration.point_mass_gravity(),
                    propagation_setup.acceleration.aerodynamic()
                ]
            )
        }
    elif detail_level == 1:
        acceleration_settings = {"Satellite":
            dict(
                Mars=
                [
                    propagation_setup.acceleration.spherical_harmonic_gravity(4, 4),
                    propagation_setup.acceleration.aerodynamic()
                ]
            )
        }
    elif detail_level == 2:
        acceleration_settings = {"Satellite":
            dict(
                Mars=
                [
                    propagation_setup.acceleration.spherical_harmonic_gravity(8, 8),
                    propagation_setup.acceleration.aerodynamic()
                ],
                Sun =
                [
                    propagation_setup.acceleration.cannonball_radiation_pressure(),
                    propagation_setup.acceleration.point_mass_gravity()
                ],
                Jupiter =
                [
                    propagation_setup.acceleration.point_mass_gravity()
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

def get_initial_state(bodies, altitude=300e3, inclination=np.deg2rad(0), eccentricity=0.01):
    # Define initial state
    mars_gravitational_parameter = bodies.get_body("Mars").gravitational_parameter
    mars_radius = spice_interface.get_average_radius("Mars")

    initial_state = element_conversion.keplerian_to_cartesian_elementwise(
        gravitational_parameter = mars_gravitational_parameter,
        semi_major_axis = mars_radius + altitude,
        eccentricity = eccentricity,
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
    dynamics_simulator = SingleArcSimulator(
        bodies, integrator_settings, propagator_settings, print_dependent_variable_data = verbose
    )

    states = dynamics_simulator.state_history
    states_elements = np.vstack(list(states.values()))

    dependent_variables = dynamics_simulator.dependent_variable_history
    cpu_time = time.time() - t0

    print("Simulation took %.2f seconds." % cpu_time)

    # Compute results
    #time_l = [t / 3600 for t in dependent_variables.keys()]
    time_l = np.array(list(states.keys()))
    time_l -= time_l[0]
    
    if len(dependent_variables) > 0:
        dependent_variable_list = np.vstack( list( dependent_variables.values( ) ) )
    else:
        dependent_variable_list = []

    if return_raw:
        return np.array(time_l), np.array(states_elements), np.array(dependent_variable_list)

    altitudes = dependent_variable_list[:,0]
    densities = dependent_variable_list[:,1]

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
    first_interpolator = interpolators.create_one_dimensional_scalar_interpolator(dict(zip(time_l, altitudes)), interpolator_settings)
    second_interpolator = interpolators.create_one_dimensional_scalar_interpolator(dict(zip(baseline_t, baseline_h)), interpolator_settings)
    # Calculate the difference between the first and second model at specific epochs
    model_difference = {t: second_interpolator.interpolate(t)- first_interpolator.interpolate(t) 
                        for t in interp_times}

    diff_times = np.array(list(model_difference.keys()))-earliest_time
    diff_vals = list(model_difference.values())
    return diff_times, diff_vals

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
    tolerance = 1e-12 if extra_accurate else 1e-8
    # Setup the optimal integrator settings
    initial_time = simulation_start_epoch # seconds since J2000
    initial_time_step = 150  # seconds
    coefficient_set = propagation_setup.integrator.RKCoefficientSets.rkdp_87 
    minimum_step_size = 0.05 # seconds
    maximum_step_size = 600  # seconds
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
        safety_factor = 0.85,
        maximum_factor_increase = 3.0,
        minimum_factor_increase = 0.25 )
    return integrator_settings

def get_integrator_settings_thrust(simulation_start_epoch):
    tolerance = 1e-7
    # Setup the optimal integrator settings
    initial_time = simulation_start_epoch # seconds since J2000
    initial_time_step = 150  # seconds
    coefficient_set = propagation_setup.integrator.RKCoefficientSets.rkdp_87 
    minimum_step_size = 1e-5 # seconds
    maximum_step_size = 600  # seconds
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
        safety_factor = 0.85,
        maximum_factor_increase = 3.0,
        minimum_factor_increase = 0.25 )
    return integrator_settings