import sys
import numpy as np
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
from utils import propagation as P
from utils import sat_models as SM
from tudatpy.kernel import constants
from tools import plot_utilities as PU

run_alt_study = True                    # Run the study of the altitude vs orbital time
run_atmo_study = False                   # Run the study of the atmospheric properties at given altitudes

sat = SM.satellite("Orbiter", 5, 3, S_ref=0.015)     # Create satellite of 5kg, with Cd of 3, and reference surface area of 0.015 m2

if run_alt_study:
    altitudes = np.arange(50, 160.1, 2.5)
    orbit_time = []
    # Run for each altitude
    for h in altitudes:
        simulation_days = 600

        # Setup the simulation
        OS = P.orbit_simulation(sat, "Mars", simulation_days*constants.JULIAN_DAY)
        OS.create_bodies(use_MCD=[True, False])
        OS.create_initial_state(h=h*1e3, i=np.deg2rad(0.01))
        OS.create_termination_settings(min_altitude=25e3)
        OS.create_dependent_variables("rho")
        OS.create_integrator()
        OS.create_accelerations(default_config=1)
        OS.create_propagator()
        
        # Run the simulation
        time, states, dep_vars = OS.simulate()
        orbit_time.append(time[-1]/constants.JULIAN_DAY)
        print("Starting from altitude of %.1f km, stay in orbit %.1e days" % (h, time[-1]/constants.JULIAN_DAY))

    # Plot the initial altitude vs orbital time
    PU.plot_single(orbit_time, altitudes, "Time in orbit [days]", "Starting altitude [km]", "MCD/feasible_altitudes", xlog=True)

if run_atmo_study:
    altitudes = [85, 115, 150]
