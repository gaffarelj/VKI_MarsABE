import sys
sys.path.insert(0,"\\".join(sys.path[0].split("\\")[:-2]))
from tudatpy.kernel import constants
from tudatpy.kernel.simulation import propagation_setup
from setup_selection import setup_utils as SU
from tools import plot_utilities as PU
from tools import time_conversions as TC
import matplotlib.pyplot as plt
import numpy as np

simulation_days = 20
simulation_start_epoch = TC.MCD_to_Tudat(2459942)
simulation_end_epoch = simulation_start_epoch + simulation_days*constants.JULIAN_DAY

# Define the environment and bodies
bodies, bodies_to_propagate, central_bodies = SU.create_bodies(use_MCD_atmo=True)
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
    propagation_setup.propagator.encke,
    dependent_variables_to_save
)
# Define the integrator settings
initial_time = simulation_start_epoch # seconds since J2000
initial_time_step = 20 # seconds
coefficient_set = propagation_setup.integrator.RKCoefficientSets.rkdp_87
minimum_step_size = 10 # seconds
maximum_step_size = 1800 # seconds
relative_error_tolerance = 7.5e-7 # -
absolute_error_tolerance = 7.5e-7 # -
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


# Run the simulation
time, altitudes, densities, cpu_time = SU.run_simulation(bodies, integrator_settings, propagator_settings)

# Make plot
PU.plot_dual(np.array(time)/3600, altitudes/1e3, densities, "Time [hr]", "Altitude [km]", "Density [kg/m$^3$]", "integ_prop/test_best_combo_%sday" % simulation_days)

# Compute the difference with the baseline
diff_times, diff_vals = SU.compare_to_baseline(time, altitudes, baseline_f="rk_4_baseline_MCD_%sday" % simulation_days, trunc_ends=True)

# Plot difference
PU.plot_single(np.array(diff_times)/3600, np.array(diff_vals), "Time [hr]", "Difference in altitude [m]", "integ_prop/diff_best_combo_%sday" % simulation_days)
#plt.plot(np.array(diff_times[:-1])/24, np.array(diff_vals[:-1])/1e3)
#plt.grid(), plt.xlabel("Time [days]"), plt.ylabel("Altitude [km]")
#plt.show()
print("Final altitude of %.2f km, at a time of %.2f days" % (altitudes[-1]/1e3, time[-1]/3600/24))
print("Max difference with rk_4 is of %.3f m, with a cpu time of %.2f seconds." % (max(diff_vals), cpu_time))

## Note: the change in settings has been done manually, starting from the integrator and propagator suggested in the previous steps.
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
# tolereance to 1e-9: 5 seconds, 0.6 km difference