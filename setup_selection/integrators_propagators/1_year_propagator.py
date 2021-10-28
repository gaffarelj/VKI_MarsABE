import sys
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
from tudatpy.kernel import constants
from tudatpy.kernel.numerical_simulation import propagation_setup
from setup_selection import setup_utils as SU
#import matplotlib.pyplot as plt
import numpy as np
from tools import time_conversions as TC

simulation_days = 20
simulation_start_epoch = TC.MCD_to_Tudat(2459942)
simulation_end_epoch = simulation_start_epoch + simulation_days*constants.JULIAN_DAY

# Define the environment and bodies
bodies, bodies_to_propagate, central_bodies = SU.create_bodies(use_MCD_atmo=True)
# Define the accelerations to be included
acceleration_models = SU.setup_environment(bodies, bodies_to_propagate, central_bodies, detail_level=1)
# Define the initial state of the satellite
initial_state = SU.get_initial_state(bodies, inclination=np.deg2rad(0.01))
# Define the termination settings
termination_settings = SU.simulation_settings(simulation_end_epoch)
# Define the dependent variables to save
dependent_variables_to_save = [
    propagation_setup.dependent_variable.altitude("Satellite", "Mars"),
    propagation_setup.dependent_variable.density("Satellite", "Mars")
]

integrator_settings = SU.get_best_integrator(simulation_start_epoch)

# Define the propagator settings
available_propagators = [propagation_setup.propagator.cowell,
						 propagation_setup.propagator.encke,
						 propagation_setup.propagator.gauss_keplerian,
						 propagation_setup.propagator.gauss_modified_equinoctial,
						 propagation_setup.propagator.unified_state_model_quaternions,
						 propagation_setup.propagator.unified_state_model_modified_rodrigues_parameters,
						 propagation_setup.propagator.unified_state_model_exponential_map]
for propagator in available_propagators:
    propagator_settings = propagation_setup.propagator.translational(
        central_bodies,
        acceleration_models,
        bodies_to_propagate,
        initial_state,
        termination_settings,
        propagator,
        dependent_variables_to_save
    )
    print("Propagator used:", propagator)

    # Run the simulation
    time, altitudes, densities, cpu_time = SU.run_simulation(bodies, integrator_settings, propagator_settings)

    # Compute the difference with the baseline
    print("Computing difference from baseline...")
    diff_times, diff_vals = SU.compare_to_baseline(time, altitudes, baseline_f="rk_4_baseline_envs_%sday" % simulation_days, trunc_ends=True)

    print("Final altitude of %.2f km, at a time of %.2f years" % (altitudes[-1]/1e3, (time[-1]-time[0])/365/24))
    print("Max difference with rk_4 is of %.3f m, with a cpu time of %.2f seconds." % (max(diff_vals), cpu_time))
    print()