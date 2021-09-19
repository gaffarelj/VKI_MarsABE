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
acceleration_models = SU.setup_environment(bodies, bodies_to_propagate, central_bodies, detail_level=1) # takes 100 s for 20 days, MCD, envs
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

#fixed_step_size = 10
fixed_step_size = 10
integrator_settings = propagation_setup.integrator.runge_kutta_4(
    simulation_start_epoch,
    fixed_step_size
)

# Run the simulation
time, altitudes, densities, cpu_time = SU.run_simulation(bodies, integrator_settings, propagator_settings)

# Make plot
PU.plot_dual(np.array(time)/3600, altitudes/1e3, densities, "Time [hr]", "Altitude [km]", "Density [kg/m$^3$]", "integ_prop/test_rk4_%sday" % simulation_days)

np.savetxt("setup_selection/integrators_propagators/rk_4_baseline_envs_%sday.dat" % simulation_days, \
   np.array([time, altitudes]), fmt="%.5e")
#np.savetxt("setup_selection/integrators_propagators/rk_4_baseline.dat", np.array([time, altitudes]), fmt="%.5e")
# baseline took 110 seconds (225 seconds including MCD atmosphere, for 50 days, step of 10 s, 7.7 seconds for 1 day, 89 seconds for 20 days)