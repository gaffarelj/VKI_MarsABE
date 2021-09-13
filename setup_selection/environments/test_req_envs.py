import sys
sys.path.insert(0,"\\".join(sys.path[0].split("\\")[:-2]))
from tudatpy.kernel import constants
from tudatpy.kernel.simulation import propagation_setup
from setup_selection import setup_utils as SU
from tools import plot_utilities as PU
from tools import time_conversions as TC
import matplotlib.pyplot as plt
import numpy as np

simulation_days = 1
simulation_start_epoch = TC.MCD_to_Tudat(2459942)
simulation_end_epoch = simulation_start_epoch + simulation_days*constants.JULIAN_DAY

# Define the environment and bodies
bodies, bodies_to_propagate, central_bodies = SU.create_bodies(use_MCD_atmo=False, use_MCD_winds=True)
# Define the initial state of the satellite
initial_state = SU.get_initial_state(bodies, altitude=450e3, inclination=np.deg2rad(90))
# Define the termination settings
termination_settings = SU.simulation_settings(simulation_end_epoch)
# Define the dependent variables to save
dependent_variables_to_save = [
    propagation_setup.dependent_variable.altitude("Satellite", "Mars"),
    propagation_setup.dependent_variable.density("Satellite", "Mars")
]
integrator_settings = SU.get_best_integrator(simulation_start_epoch, extra_accurate=True)

check_grav = False
check_rel = False
check_rad = True

# Define the accelerations to be included
if check_grav:
    # Start with a base acceleration: Mars PM and aero
    acceleration_settings = {"Satellite":
        dict(
            Mars=
            [
                propagation_setup.acceleration.point_mass_gravity(),
                propagation_setup.acceleration.aerodynamic()
            ]
        )}
    accelerations_list = [propagation_setup.create_acceleration_models(bodies, acceleration_settings, bodies_to_propagate, central_bodies)]
    # Mars SH D/O 2 and aero
    acceleration_settings = {"Satellite":
        dict(
            Mars=
            [
                propagation_setup.acceleration.spherical_harmonic_gravity(2, 2),
                propagation_setup.acceleration.aerodynamic()
            ]
        )}
    accelerations_list.append(propagation_setup.create_acceleration_models(bodies, acceleration_settings, bodies_to_propagate, central_bodies))

    # Mars SH D/O 4 and aero
    acceleration_settings = {"Satellite":
        dict(
            Mars=
            [
                propagation_setup.acceleration.spherical_harmonic_gravity(4, 4),
                propagation_setup.acceleration.aerodynamic()
            ]
        )}
    accelerations_list.append(propagation_setup.create_acceleration_models(bodies, acceleration_settings, bodies_to_propagate, central_bodies))

    # Mars SH D/O 8 and aero
    acceleration_settings = {"Satellite":
        dict(
            Mars=
            [
                propagation_setup.acceleration.spherical_harmonic_gravity(8, 8),
                propagation_setup.acceleration.aerodynamic()
            ]
        )}
    accelerations_list.append(propagation_setup.create_acceleration_models(bodies, acceleration_settings, bodies_to_propagate, central_bodies))

    legends = ["Mars PM+aero", "Mars SH D/O 2, aero", "Mars SH D/O 4, aero", "Mars SH D/O 8, aero"]
    fname = "comparison_grav"
elif check_rel:
    # Mars SH D/O 4 and aero
    acceleration_settings = {"Satellite":
        dict(
            Mars=
            [
                propagation_setup.acceleration.spherical_harmonic_gravity(4, 4),
                propagation_setup.acceleration.aerodynamic()
            ]
        )}
    accelerations_list = [propagation_setup.create_acceleration_models(bodies, acceleration_settings, bodies_to_propagate, central_bodies)]

    # Mars SH D/O 4 and aero + Schwarzschild relativistic correction
    acceleration_settings = {"Satellite":
        dict(
            Mars=
            [
                propagation_setup.acceleration.spherical_harmonic_gravity(4, 4),
                propagation_setup.acceleration.aerodynamic(),
                propagation_setup.acceleration.relativistic_correction(use_schwarzschild = True)
            ]
        )}
    accelerations_list.append(propagation_setup.create_acceleration_models(bodies, acceleration_settings, bodies_to_propagate, central_bodies))

    legends = ["Mars SH D/O 4, aero", "+ Schwarzschild relativistic correction"]
    fname = "comparison_relcorr"
elif check_rad:
    # Mars SH D/O 4 and aero
    acceleration_settings = {"Satellite":
        dict(
            Mars=
            [
                propagation_setup.acceleration.spherical_harmonic_gravity(4, 4),
                propagation_setup.acceleration.aerodynamic()
            ]
        )}
    accelerations_list = [propagation_setup.create_acceleration_models(bodies, acceleration_settings, bodies_to_propagate, central_bodies)]

    # Mars SH D/O 4 and aero + Cannonball radiation pressure
    acceleration_settings = {"Satellite":
        dict(
            Mars=
            [
                propagation_setup.acceleration.spherical_harmonic_gravity(4, 4),
                propagation_setup.acceleration.aerodynamic()
            ],
            Sun = [ propagation_setup.acceleration.cannonball_radiation_pressure() ] 
        )}
    accelerations_list.append(propagation_setup.create_acceleration_models(bodies, acceleration_settings, bodies_to_propagate, central_bodies))

    legends = ["Mars SH D/O 4, aero", "+ Cannonball radiation pressure"]
    fname = "comparison_rad"

pos_t = []
vel_t = []
times = []

for i, accelerations in enumerate(accelerations_list):
    # Define the propagator settings
    print("Acceleration:", legends[i])
    propagator_settings = propagation_setup.propagator.translational(
        central_bodies,
        accelerations,
        bodies_to_propagate,
        initial_state,
        termination_settings,
        propagation_setup.propagator.cowell,
        dependent_variables_to_save
    )

    # Run the simulation
    time, states, dep_vars = SU.run_simulation(bodies, integrator_settings, propagator_settings, return_raw=True)

    # Compute the positions and velocities
    positions = np.linalg.norm(states[:,:3], axis=1)
    velocities = np.linalg.norm(states[:,3:], axis=1)

    times.append(time)
    pos_t.append(np.fabs(positions-positions[0])/1e3)
    vel_t.append(np.fabs(velocities-velocities[0]))

# Make plot
PU.plot_multiple(times, pos_t, "Time [hr]", "$|r(t) - r_0|$ [km]", "envs/%s_pos" % fname, legends=legends, legend_loc=4, lstyle="--")
PU.plot_multiple(times, vel_t, "Time [hr]", "$|v(t) - v_0|$ [m/s]", "envs/%s_vel" % fname, legends=legends, legend_loc=4, lstyle="--")