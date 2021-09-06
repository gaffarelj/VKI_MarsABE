import sys
sys.path.insert(0,"\\".join(sys.path[0].split("\\")[:-2]))
from tudatpy.kernel import constants
from tudatpy.kernel.simulation import propagation_setup
from setup_selection import setup_utils as SU
import matplotlib.pyplot as plt
import numpy as np

simulation_start_epoch = 0.0
simulation_end_epoch = 1*constants.JULIAN_YEAR

# Define the environment and bodies
bodies, bodies_to_propagate, central_bodies = SU.create_bodies()
# Define the accelerations to be included
acceleration_models = SU.setup_environment(bodies, bodies_to_propagate, central_bodies)
# Define the initial state of the satellite
initial_state = SU.get_initial_state(bodies)
# Define the termination settings
termination_settings = SU.simulation_settings(simulation_end_epoch)
# Define the dependent variables to save
dependent_variables_to_save = [
    propagation_setup.dependent_variable.altitude("Satellite", "Mars"),
    propagation_setup.dependent_variable.density("Satellite", "Mars")
]

# Define the propagator settings
propagator_settings = propagation_setup.propagator.translational(
    central_bodies,
    acceleration_models,
    bodies_to_propagate,
    initial_state,
    termination_settings,
    output_variables = dependent_variables_to_save
)

initial_time = 0 # seconds since J2000
initial_time_step = 20 # seconds
coefficient_set = propagation_setup.integrator.RKCoefficientSets.rkdp_87
minimum_step_size = 10 # seconds
maximum_step_size = 300 # seconds
relative_error_tolerance = 2.5e-8 # -
absolute_error_tolerance = 2.5e-8 # -

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
    safety_factor = 0.75,
    maximum_factor_increase = 4.0,
    minimum_factor_increase = 0.05 )

#extrapolation_sequence = propagation_setup.integrator.ExtrapolationMethodStepSequences.bulirsch_stoer_sequence
#maximum_number_of_steps = 3
#integrator_settings = propagation_setup.integrator.bulirsch_stoer(
#	initial_time,
#	initial_time_step,
#	extrapolation_sequence,
#	maximum_number_of_steps,
#	minimum_step_size,
#	maximum_step_size,
#	relative_error_tolerance,
#	absolute_error_tolerance
#)

#minimum_order = 6
#maximum_order = 11
#integrator_settings = propagation_setup.integrator.adams_bashforth_moulton(
#	initial_time,
#	initial_time_step,
#	minimum_step_size,
#	maximum_step_size,
#	relative_error_tolerance,
#	absolute_error_tolerance,
#	minimum_order,
#	maximum_order
#)


#fixed_step_size = 60
#integrator_settings = propagation_setup.integrator.runge_kutta_4(
#    simulation_start_epoch,
#    fixed_step_size
#)

# Run the simulation
time, altitudes, densities, cpu_time = SU.run_simulation(bodies, integrator_settings, propagator_settings)

# Make plots
plt.plot(np.array(time)/24, altitudes/1e3)
plt.grid(), plt.xlabel("Time [days]"), plt.ylabel("Altitude [km]")
plt.show()

# Compute the difference with the baseline
diff_times, diff_vals = SU.compare_to_baseline(time, altitudes)

plt.plot(np.array(diff_times[:-1])/24, np.array(diff_vals[:-1])/1e3)
plt.grid(), plt.xlabel("Time [days]"), plt.ylabel("Altitude [km]")
plt.show()
print("Final altitude of %.2f km, at a time of %.2f years" % (altitudes[-1]/1e3, time[-1]/365/24))
print("Max difference with rk_4 is of %.3f km, with a cpu time of %.2f seconds." % (max(diff_vals)/1e3, cpu_time))

# rk_4 from baseline (step of 10s) is around 110 seconds
# rk_4 (step of 30s) is around 35 seconds, 11 km difference
# rk_4 (step of 60s) is around 17 seconds, 70 km difference

# rkf_45 (step of 1-300s, tolerance of 1e-9) is around 42 seconds, 5 km difference
# rkf_45 (step of 1-300s, tolerance of 1e-8) is around 27 seconds, 48 km difference
# rkf_45 (step of 1-300s, tolerance of 1e-6) is around 11 seconds, 130 km difference
# rkf_45 (step of 1-500s, tolerance of 1e-9) is around 40 seconds, 5 km difference

# rkf_56 (step of 10-300s, tolerance of 1e-9) is around 27 seconds, 32 km difference

# rkf_78 (step of 10-300s, tolerance of 1e-9) is around 16 seconds, 9 km difference

# rkdp_87 (step of 10-300s, tolerance of 1e-9) is around 15 seconds, 0.16 km difference
# rkdp_87 (step of 10-300s, tolerance of 2.5e-8) is around 10.5 seconds, 0.9 km difference
# rkdp_87 (step of 10-300s, tolerance of 1e-8) is around 12 seconds, 0.8 km difference
# rkdp_87 (step of 10-300s, tolerance of 1e-7) is around 10 seconds, 1 km difference

# abm (step of 10-300s, tolerance of 1e-9, order 6-11) is around 44 seconds, 0.04 km difference
# abm (step of 30-300s, tolerance of 1e-9, order 6-11) is around 20 seconds, 25 km difference
# abm (step of 10-500s, tolerance of 1e-9, order 6-11) is around 45 seconds, 0.05 km difference

# bs (step of 10-300s, tolerance of 1e-9, max 5 steps) is around 28 seconds, 0.1 km difference
# bs (step of 10-300s, tolerance of 1e-9, max 4 steps) is around 22 seconds, 0.6 km difference

# suggested: 
# rkdp_87 (step of 10-300s, tolerance of 2.5e-8) is around 10.5 seconds, 0.9 km difference