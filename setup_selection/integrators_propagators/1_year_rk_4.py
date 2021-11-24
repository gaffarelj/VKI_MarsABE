import sys
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
from tudatpy.kernel import constants
from tudatpy.kernel.numerical_simulation import propagation_setup
from setup_selection import setup_utils as SU
from tools import plot_utilities as PU
from tools import time_conversions as TC
#import matplotlib.pyplot as plt
import numpy as np

simulation_days = 20
# simulation_start_epoch = TC.MCD_to_Tudat(2459942)
# simulation_end_epoch = simulation_start_epoch + simulation_days*constants.JULIAN_DAY

# # Define the environment and bodies
# bodies, bodies_to_propagate, central_bodies = SU.create_bodies(use_MCD_atmo=True)
# # Define the accelerations to be included
# acceleration_models = SU.setup_environment(bodies, bodies_to_propagate, central_bodies, detail_level=1) # takes 100 s for 20 days, MCD, envs
# # Define the initial state of the satellite
# initial_state = SU.get_initial_state(bodies, inclination=np.deg2rad(0.01))
# # Define the termination settings
# termination_settings = SU.simulation_settings(simulation_end_epoch)
# # Define the dependent variables to save
# dependent_variables_to_save = [
#     propagation_setup.dependent_variable.altitude("Satellite", "Mars"),
#     propagation_setup.dependent_variable.density("Satellite", "Mars")
# ]

# # Define the propagator settings
# propagator_settings = propagation_setup.propagator.translational(
#     central_bodies,
#     acceleration_models,
#     bodies_to_propagate,
#     initial_state,
#     termination_settings,
#     output_variables = dependent_variables_to_save
# )

# Setup simulation to use thrust
from utils import propagation as P
from utils import sat_models as SM
OS = P.orbit_simulation(SM.satellites["CS_1021"], "Mars", simulation_days*constants.JULIAN_DAY)

OS.create_bodies(use_MCD=[False, False], use_GRAM=False)
h_a, h_p = 160e3, 140e3
a = OS.R_cb + (h_a+h_p)/2
e = 1 - (OS.R_cb + min(h_p, h_a)) / a
OS.create_initial_state(a=a, e=e, i=0.25)
OS.create_accelerations(default_config=1, thrust=3, ionisation_eff=0.7, use_battery=True)
OS.create_integrator()
#fixed_step_size = 10
fixed_step_size = 10
OS.integrator_settings = propagation_setup.integrator.runge_kutta_4(
    OS.init_time,
    fixed_step_size
)
OS.create_termination_settings()
OS.create_dependent_variables(to_save=["h"])
OS.create_propagator()

# Run the simulation
#time, altitudes, densities, cpu_time = SU.run_simulation(bodies, integrator_settings, propagator_settings)
time, states, _ = OS.simulate()
altitudes = OS.get_dep_var("h")

# Make plot
#PU.plot_dual(np.array(time)/3600, altitudes/1e3, densities, "Time [hr]", "Altitude [km]", "Density [kg/m$^3$]", "integ_prop/test_rk4_%sday" % simulation_days)

np.savetxt("setup_selection/integrators_propagators/rk_4_baseline_thrust_%sday.dat" % simulation_days, \
   np.array([time, altitudes]), fmt="%.5e")
# np.savetxt("setup_selection/integrators_propagators/rk_4_baseline_envs_%sday.dat" % simulation_days, \
#    np.array([time, altitudes]), fmt="%.5e")
#np.savetxt("setup_selection/integrators_propagators/rk_4_baseline.dat", np.array([time, altitudes]), fmt="%.5e")
# baseline took 110 seconds (225 seconds including MCD atmosphere, for 50 days, step of 10 s, 7.7 seconds for 1 day, 89 seconds for 20 days)