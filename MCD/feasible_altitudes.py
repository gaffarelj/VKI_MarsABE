import sys
sys.path.insert(0,"\\".join(sys.path[0].split("\\")[:-2]))
from tudatpy.kernel import constants
from tudatpy.kernel.simulation import propagation_setup
from tudatpy.kernel.interface import spice_interface
from setup_selection import setup_utils as SU
from tools import plot_utilities as PU
from tools import time_conversions as TC
import thrust.simple_thrust as T
import matplotlib.pyplot as plt
import numpy as np

altitudes = np.arange(80, 220.1, 2.5)
#altitudes = [85, 115, 150]
orbit_time = []

for h in altitudes:

    simulation_days = 600
    simulation_start_epoch = TC.MCD_to_Tudat(2459942)
    simulation_end_epoch = simulation_start_epoch + simulation_days*constants.JULIAN_DAY

    # Define the environment and bodies
    bodies, bodies_to_propagate, central_bodies = SU.create_bodies(use_MCD_atmo=True)
    # Define the initial state of the satellite
    initial_state = SU.get_initial_state(bodies, altitude=h*1e3, inclination=np.deg2rad(0.01), eccentricity=0)
    # Define the termination settings
    termination_settings = SU.simulation_settings(simulation_end_epoch)
    # Define the dependent variables to save
    dependent_variables_to_save = [
        propagation_setup.dependent_variable.relative_speed("Satellite", "Mars"),
        propagation_setup.dependent_variable.density("Satellite", "Mars")
    ]
    integrator_settings = SU.get_best_integrator(simulation_start_epoch)

    accelerations = SU.setup_environment(bodies, bodies_to_propagate, central_bodies, detail_level=1)

    # Define the propagator settings
    propagator_settings = propagation_setup.propagator.translational(
        central_bodies,
        accelerations,
        bodies_to_propagate,
        initial_state,
        termination_settings,
        propagation_setup.propagator.gauss_modified_equinoctial,
        dependent_variables_to_save
    )

    # Run the simulation
    time, states, dep_vars = SU.run_simulation(bodies, integrator_settings, propagator_settings, return_raw=True)
    orbit_time.append(time[-1]/constants.JULIAN_DAY)
    print("Starting from altitude of %i km, stay in orbit %.1e days" % (h, time[-1]/constants.JULIAN_DAY))
    print("Orbital velocity of %.3f m/s, air density of %.5e kg/m3" % (np.mean(dep_vars[:10,0]), np.mean(dep_vars[:10,1])))

PU.plot_single(orbit_time, altitudes, "Time in orbit [days]", "Starting altitude [km]", "MCD/feasible_altitudes", xlog=True)