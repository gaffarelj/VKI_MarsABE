import sys
import numpy as np
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
from utils import propagation as P
from utils import sat_models as SM
from tudatpy.kernel import constants
from tudatpy.kernel.numerical_simulation import propagation_setup
from tools import plot_utilities as PU

run_alt_study = False       # Run the study of the altitude vs orbital time
run_atmo_study = True       # Run the study of the atmospheric properties at given altitudes

if run_alt_study:
    # Create satellite of 5kg, with Cd of 3, and reference surface area of 0.015 m2
    sat = SM.satellite("Orbiter", 5, 3, S_ref=0.015)
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
    from MCD import parallel_mcd as pmcd
    altitudes = [85, 115, 150]
    # Create satellite of 5kg, with Cd of 0, and reference surface area of 0.015 m2 (no drag, but still measure the density)
    sat = SM.satellite("Orbiter", 5, 0, S_ref=0)
    # Setup the constant parts of the simulation
    OS = P.orbit_simulation(sat, "Mars", 2*650*constants.JULIAN_DAY)
    OS.create_bodies(use_MCD=[True, False], preload_MCD=False, save_MCD_vals=True)
    OS.create_termination_settings(min_altitude=25e3)
    OS.create_dependent_variables(["h", "V", "rho"])
    OS.create_integrator(tolerance=1e-6, dt=[0.1, 250, 1e5])
    # Use a fixed step integrator with dt = 250s
    OS.integrator_settings = propagation_setup.integrator.runge_kutta_4(
        OS.init_time,
        250,
        save_frequency = 1
    )
    # Only use Mars as a point mass as the acceleration (stay in the same orbit this way)
    accelerations = [P.env_acceleration("Mars", PM=True, aero=True)]
    OS.create_accelerations(env_accelerations=accelerations)
    for h in altitudes:
        # Setup the simulation with the given altitude
        OS.create_initial_state(h_p=h*1e3, i=np.deg2rad(45))
        OS.create_propagator()
        # Make sure the array to save the MCD values are empty
        pmcd.ALL_VALUES, pmcd.TIMES = [], []
        # Run the simulation
        time, states, dep_vars = OS.simulate()
        print("* Altitude of %i km" % h)
        MCD_vals = np.array(pmcd.ALL_VALUES)
        MCD_times = np.array(pmcd.TIMES)
        density, temperature, pressure, mixture = MCD_vals[:,0], MCD_vals[:,1], MCD_vals[:,2], MCD_vals[:,3:]
        # Make sure sum of mixtures is 1
        mixture = mixture[:,:6] / sum(np.mean(mixture[:,:6], axis=0))
        velocity = np.linalg.norm(states[:,3:], axis=1)
        altitude_hist = OS.get_dep_var("h")
        titles = ["Altitude", "Velocity", "Density", "Temperature", "Pressure"]
        units = ["m", "m/s", "kg/m3", "K", "Pa"]
        vals = [altitude_hist, velocity, density, temperature, pressure]
        # Print average results (with their standard deviations)
        for i in range(len(titles)):
            print("%s: mean of %.5e %s (std of %.5e %s)" % (titles[i], np.mean(vals[i]), units[i], np.std(vals[i]), units[i]))
        print("Mixture (average): CO2, N2, Ar, CO, O, O2:", np.array2string(np.mean(mixture, axis=0)*100, \
            separator="%, ", formatter={"float_kind": lambda x: "%.3f" % x}))
        print("Mixture (std):     CO2, N2, Ar, CO, O, O2:", np.array2string(np.std(mixture, axis=0)*100, \
            separator="%, ", formatter={"float_kind": lambda x: "%.3f" % x}))