import sys
sys.path.insert(0,"\\".join(sys.path[0].split("\\")[:-2]))
from tudatpy.kernel import constants
from tudatpy.kernel.simulation import propagation_setup
from tudatpy.kernel.interface import spice_interface
from setup_selection import setup_utils as SU
from tools import plot_utilities as PU
from tools import time_conversions as TC
from tools import mission_geometry as MG
import matplotlib.pyplot as plt
import numpy as np

simulation_days = 1
simulation_start_epoch = TC.MCD_to_Tudat(2459942)
simulation_end_epoch = simulation_start_epoch + simulation_days*constants.JULIAN_DAY

# Define the environment and bodies
bodies, bodies_to_propagate, central_bodies = SU.create_bodies(use_MCD_atmo=False, use_MCD_winds=True)
# Define the initial state of the satellite
initial_state = SU.get_initial_state(bodies)
# Define the termination settings
termination_settings = SU.simulation_settings(simulation_end_epoch)
# Define the dependent variables to save
dependent_variables_to_save = [
    propagation_setup.dependent_variable.relative_position("Satellite", "Mars"),
    propagation_setup.dependent_variable.relative_position("Sun", "Mars")
]
integrator_settings = SU.get_best_integrator(simulation_start_epoch)

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

    times.append(time/3600)
    pos_t.append(np.fabs(positions-positions[0])/1e3)
    vel_t.append(np.fabs(velocities-velocities[0]))

    # Check radiation pressure over time
    if check_rad:
        sat_pos = dep_vars[:,:3]
        sun_pos = dep_vars[:,3:]
        r_mars = [0, 0, 0]
        R_mars = spice_interface.get_average_radius("Mars")
        R_sun = spice_interface.get_average_radius("Sun")
        shadow_vals = []
        for i in range(len(sat_pos)):
            r_sat, r_sun = sat_pos[i], sun_pos[i]
            shadow = MG.shadow_function(r_sun, R_sun, r_mars, R_mars, r_sat)
            shadow_vals.append(shadow)
        # Make plot
        PU.plot_single(time/3600, shadow_vals, "Time [hr]", "Shadow value [-]", "envs/radiation_shadow")

# Make plot
PU.plot_multiple(times, pos_t, "Time [hr]", "$|r(t) - r_0|$ [km]", "envs/%s_pos" % fname, legends=legends, legend_loc=4, lstyle="--")
PU.plot_multiple(times, vel_t, "Time [hr]", "$|v(t) - v_0|$ [m/s]", "envs/%s_vel" % fname, legends=legends, legend_loc=4, lstyle="--")