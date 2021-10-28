import sys
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
from tudatpy.kernel import constants
from tudatpy.kernel.numerical_simulation import propagation_setup
from setup_selection import setup_utils as SU
#import matplotlib.pyplot as plt
#import numpy as np

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

for i_integrator in range(1, 18):
    integrator_settings = SU.get_integrator_settings(i_integrator, verbose=True)

    # Run the simulation
    time, altitudes, densities, cpu_time = SU.run_simulation(bodies, integrator_settings, propagator_settings)

    # Compute the difference with the baseline
    print("Computing difference from baseline...")
    diff_times, diff_vals = SU.compare_to_baseline(time, altitudes)

    print("Final altitude of %.2f km, at a time of %.2f years" % (altitudes[-1]/1e3, time[-1]/365/24))
    print("Max difference with rk_4 is of %.3f km, with a cpu time of %.2f seconds." % (max(diff_vals)/1e3, cpu_time))
    print()

# suggested: 
# rkdp_87 (step of 10-300s, tolerance of 2.5e-8) is around 10.5 seconds, 0.9 km difference (setting at index number 10)