import numpy as np
import sys
sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
from tools import time_conversions as TC
import thrust as T
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel import constants
from tudatpy.kernel.simulation import environment_setup
from tudatpy.kernel.simulation import propagation_setup
from tudatpy.kernel.astro import conversion
# Load the SPICE kernel
spice_interface.load_standard_kernels()


class env_acceleration:
    
    def __init__(self, body_name, PM=False, SH=False, SH_do=[2,2], aero=False, rad=False):
        """
        Environmental accelerations class, containing all of the acceleration settings for one body
        Inputs:
         * PM (bool): Point Mass gravitational acceleration should be included
         * SH (bool): Spherical Harmonics gravitational acceleration should be included
         * SH_do ([int]*2): degree and order up to which the Spherical Harmonics should be computed
         * aero (bool): Aerodynamic acceleration should be included
         * rad (bool): Cannonball radiation pressure should be included
        """
        self.body_name = body_name
        self.PM = PM
        self.SH = SH
        self.SH_do = SH_do
        self.aero = aero
        self.rad = rad
        self.a = []
        if PM and SH:
            print("Warning: both the PM and SH gravity should not be used for the same body.")
        if self.PM:
            self.a.append(propagation_setup.acceleration.point_mass_gravity())
        if self.SH:
            self.a.append(propagation_setup.acceleration.spherical_harmonic_gravity(*tuple(self.SH_do)))
        if self.aero:
            self.a.append(propagation_setup.acceleration.aerodynamic())
        if self.rad:
            self.a.append(propagation_setup.acceleration.cannonball_radiation_pressure())

class orbit_simulation:
    
    def __init__(self, sat, central_body, sim_time, init_time=TC.MCD_to_Tudat(2459942)):
        """
        Orbital simulation class, containing all the code required for setup, and simulation run.
        Inputs:
         * sat (utils.sat_models.satellite): satellite for which the orbits is simulation
         * central_body (str): body around which the satellite orbits
         * init_time (float): initial simulation time in seconds since J2000
        """
        self.sat = sat
        self.central_body = central_body
        self.init_time = init_time
        self.end_time = init_time + sim_time

    def create_bodies(self, additional_bodies=["Sun", "Jupiter"], use_MCD=[False, False]):
        """
        Create the simulation bodies.
        Inputs:
         * additional_bodies ([string]): bodies to create in addition to the central body (by default, the two most massives in the solar system)
         * use_MCD ([bool, bool]): first boolean to indicate wether the MCD atmosphere model should be implemented. If True, the second boolean is used to indicate whether winds should also be added
        """
        bodies_to_create = [self.central_body] + additional_bodies      # Bodies that will be created and used in the simulation
        global_frame_origin = self.central_body                         # Body at the centre of the simulation
        global_frame_orientation = "ECLIPJ2000"                         # Orientation of the reference frame
        # Setup TUDAT body settings
        body_settings = environment_setup.get_default_body_settings(
            bodies_to_create,
            global_frame_origin,
            base_frame_orientation=global_frame_orientation
        )

        if use_MCD[0]:
            # Use the atmospheric model from the Mars Climate Database
            from MCD.parallel_mcd import parallel_mcd as PMCD
            mcd = PMCD()
            body_settings.get("Mars").atmosphere_settings = environment_setup.atmosphere.custom_constant_temperature_detailed(
                mcd.density, constant_temperature=210, specific_gas_constant=192, ratio_of_specific_heats=1.3)
            # Values taken from https://meteor.geol.iastate.edu/classes/mt452/Class_Discussion/Mars-physical_and_orbital_statistics.pdf

            if use_MCD[1]:
                # Add winds from the MCD. Only possible if the MCD atmospheric model is used
                body_settings.get("Mars").atmosphere_settings.wind_settings = environment_setup.atmosphere.custom_wind_model(mcd.wind)
        else:
            if use_MCD[1]:
                print("Warning: the MCD winds can only be added if the MCD atmosphere is used as well")
            # Use an exponential atmosphere model
            # Exponential parameters taken from http://link.springer.com/content/pdf/10.1007%2F978-3-540-73647-9_3.pdf
            density_scale_height, density_at_zero_altitude = 7.295e3, 0.0525
            body_settings.get("Mars").atmosphere_settings = environment_setup.atmosphere.exponential(
                density_scale_height, density_at_zero_altitude)
    	
        # Create the system of bodies
        self.bodies = environment_setup.create_system_of_bodies(body_settings)
        
        self.bodies.create_empty_body(self.sat.name)                                # Add satellite body
        self.bodies.get_body(self.sat.name).set_constant_mass(self.sat.wet_mass)    # Add mass to the satellite

        # Define the aerodynamic coefficient settings
        def force_coefficients(_):
            # Get the altitude from the flight conditions
            h = self.bodies.get_body(self.sat.name).get_flight_conditions().current_altitude
            # Interpolate the drag coefficient given the altitude
            C_d = self.sat.get_cd(h)
            return [C_d, 0, 0]
        aero_coefficient_settings = environment_setup.aerodynamic_coefficients.custom(force_coefficients, self.sat.S_ref, [])
        # Add the aerodynamic settings to the satellite
        environment_setup.add_aerodynamic_coefficient_interface(self.bodies, self.sat.name, aero_coefficient_settings)

        # Define solar radiation settings
        reference_area_radiation, radiation_pressure_coefficient = 0.125, 1.2
        occulting_bodies = [self.central_body]
        radiation_pressure_settings = environment_setup.radiation_pressure.cannonball(
            "Sun", reference_area_radiation, radiation_pressure_coefficient, occulting_bodies)
        # Add the solar radiation settings to the satellite
        environment_setup.add_radiation_pressure_interface(self.bodies, self.sat.name, radiation_pressure_settings)

    def create_initial_state(self, h=150e3, e=0, i=0, omega=0, Omega=0, theta=0):
        """
        Define the initial state of the Satellite.
        Inputs:
         * h (float): altitude of the satellite above the central body, in meters
         * e (float): eccentricity of the orbit
         * i (float): inclination of the orbit in radians
         * omega (float): argument of periapsis of the orbit in radians
         * Omega (float): longitude of the ascending noce in radians
         * theta (float): true anomaly in radians (note: this defines where the satellite starts in its orbit, and can most of the time be ignored)
        """
        # Get the gravitational parameter and radius of the central body
        mu = self.bodies.get_body(self.central_body).gravitational_parameter
        R = spice_interface.get_average_radius(self.central_body)
        # Convert the Keplerian orbit to an initial cartesian state
        self.initial_state = conversion.keplerian_to_cartesian(
            gravitational_parameter = mu,
            semi_major_axis = R + h,
            eccentricity = e,
            inclination = i,
            argument_of_periapsis = omega,
            longitude_of_ascending_node = Omega,
            true_anomaly = theta
        )

    def create_accelerations(self, env_accelerations=[], default_config=None, thrust=None):
        """
        Create the accelerations for the simulated orbit.
        Inputs:
         * env_accelerations ([env_acceleration]*N): list of environmental accelerations to use (class env_acceleration)
         * default_congif (None or int): index of a default configuration to use
         * thrust (None or int): index of the thrust model to use from the thrust class
        """
        # Load a default environmental acceleration setup
        if default_config is not None:
            if len(env_accelerations) != 0:
                # If accelerations were also provided, inform that they will be overwritten
                print("Warning: the provided environmental accelerations have been overwritten")
            if default_config == 0:
                # Central body point mass and aerodynamics
                env_accelerations = [env_acceleration(self.central_body, PM=True, aero=True)]
            elif default_config == 1:
                # Central body spherical harmonics of degree/order 4 and aerodynamics
                env_accelerations = [env_acceleration(self.central_body, SH=True, SH_do=[4,4], aero=True)]
            elif default_config == 2:
                # Central body spherical harmonics of degree/order 8, aerodynamics, Solar PM and radiation, Jupiter PM
                env_accelerations = [
                    env_acceleration(self.central_body, SH=True, SH_do=[8,8], aero=True),
                    env_acceleration("Sun", PM=True, rad=True),
                    env_acceleration("Jupiter", PM=True),
                ]
        # Convert the provided accelerations to a dict (format for TUDAT)
        self.accelerations = dict()
        for env_a in env_accelerations:
            self.accelerations[env_a.body_name] = env_a.a
        # Add the thrust as an acceleration
        if thrust is not None:
            self.accelerations[self.sat.name] =  [T.thrust_settings(self, thrust)]
        # Create the acceleration models
        self.acceleration_models = propagation_setup.create_acceleration_models(
            self.bodies,
            {self.sat.name: self.accelerations},
            [self.sat.name],
            [self.central_body]
        )

    def create_integrator(self, tolerance=1e-7, dt=[1e-5, 150, 600], sf=[0.85, 0.25, 3.0]):
        """
        Create the RK-DP 87 variable step integrator for the orbital simulation.
        Inputs:
         * tolerance (float): integrator tolerance, this implicitely tunes the integration steps
         * dt ([float]*3): time step settings: 0=minimum, 1=initial, 2=maximum
         * sf ([float]*3): safety factor settings: 0=initial, 1=minimum increase factor, 2=maximum increase factor
        """
        # Setup the optimal integrator settings
        initial_time = self.init_time           # seconds since J2000
        initial_time_step = dt[1]                 # seconds
        coefficient_set = propagation_setup.integrator.RKCoefficientSets.rkdp_87 
        minimum_step_size = dt[0]                # seconds
        maximum_step_size = dt[2]                 # seconds
        relative_error_tolerance = tolerance    # -
        absolute_error_tolerance = tolerance    # -
        self.integrator_settings = propagation_setup.integrator.runge_kutta_variable_step_size(
            initial_time,
            initial_time_step,
            coefficient_set,
            minimum_step_size,
            maximum_step_size,
            relative_error_tolerance,
            absolute_error_tolerance,
            save_frequency= 1,
            assess_termination_on_minor_steps = False,
            safety_factor = sf[0],
            maximum_factor_increase = sf[1],
            minimum_factor_increase = sf[2] )

    def create_termination_settings(self, min_altitude=50e3, max_altitude=500e3, cpu_time=60*60):
        """
        Create the termination settings; conditions for when the simulation must stop.
        Inputs:
         * min_altitude (float): altitude below which propagation stops, in meters (set to 0 to remove)
         * max_altitude (float): altitude above which propagation stops, in meters (set to Inf to remove)
         * cpu_time (float): CPU time after which the propagation is stopped, in seconds (set to Inf to remove)
        """
        # Create a lower altitude termination setting (when sat is too deep in the atmosphere or crashed)
        termination_altitude = propagation_setup.dependent_variable.altitude(self.sat.name, self.central_body)
        termination_setting_min_altitude = propagation_setup.propagator.dependent_variable_termination(
                dependent_variable_settings = termination_altitude, limit_value = min_altitude, use_as_lower_limit = True)
        # Create an upper altitude termination setting (when sat is too high above the planet)
        termination_setting_max_altitude = propagation_setup.propagator.dependent_variable_termination(
                dependent_variable_settings = termination_altitude, limit_value = max_altitude, use_as_lower_limit = False)
        # Create a time termination setting (stop after x days in orbit)
        termination_setting_time = propagation_setup.propagator.time_termination(self.end_time)
        # Create a CPU time termination setting (stop after x minutes of CPU time)
        termination_setting_cpu_time = propagation_setup.propagator.cpu_time_termination(cpu_time)
        # Assemble the termination settings together (stop after one of the termination setting is reached)
        termination_settings_list = [termination_setting_min_altitude, termination_setting_max_altitude, termination_setting_time, termination_setting_cpu_time]
        self.termination_settings = propagation_setup.propagator.hybrid_termination( 
            termination_settings_list, fulfill_single_condition = True )

    def create_dependent_variables(self, to_save=[]):
        """
        Dependent variables to save during the propagation.
        Inputs:
         * to_save ([string]*N): string representing which dependent variable to save:
           * h:     altitude of the satellite
           * rho:   atmospheric density at the position of the satellite
           * V:     airspeed of the Satellite
           * m:     mass of the satellite
           * F_T:   acceleration due to the Thrust
           * D:     acceleration due to the atmosphere
           * C_D:   aerodynamic coefficients
           * r_cb:  relative position of the central body w.r.t. the sun
           * Kep:   Keplerian state of the satellite
        """
        # Start with no dependent variables to save
        self.dependent_variables_to_save = []
        # Dictionnary to save which dependent variables can be accessed at what index
        self.dep_var_loc = dict()
        idx = 0
        for dep_key in to_save:
            if dep_key == "h":
                # Altitude of the satellite
                d_v = propagation_setup.dependent_variable.altitude(self.sat.name, self.central_body)
                size = 1
            elif dep_key == "rho":
                # Atmospheric density at the position of the satellite
                d_v = propagation_setup.dependent_variable.density(self.sat.name, self.central_body)
                size = 1
            elif dep_key == "V":
                # Airspeed of the Satellite
                propagation_setup.dependent_variable.body_fixed_airspeed_velocity(self.sat.name, self.central_body)
                size = 3
            elif dep_key == "m":
                # Mass of the satellite
                d_v = propagation_setup.dependent_variable.body_mass(self.sat.name)
                size = 1
            elif dep_key == "F_T":
                # Acceleration due to the Thrust
                d_v = propagation_setup.dependent_variable.single_acceleration_norm(
                    propagation_setup.acceleration.thrust_acceleration_type, self.sat.name, self.sat.name)
                size = 3
            elif dep_key == "D":
                # Acceleration due to the atmosphere
                d_v = propagation_setup.dependent_variable.single_acceleration_norm(
                    propagation_setup.acceleration.aerodynamic_type, self.sat.name, self.central_body)
                size = 3
            elif dep_key == "C_D":
                # Aerodynamic coefficients
                propagation_setup.dependent_variable.aerodynamic_force_coefficients(self.sat.name)
                size = 3
            elif dep_key == "r_cb":
                # Relative position of the central body w.r.t. the sun
                propagation_setup.dependent_variable.relative_position(self.sat.name, self.central_body)
                size = 3
            elif dep_key == "Kep":
                # Keplerian state of the satellite
                propagation_setup.dependent_variable.keplerian_state(self.sat.name, self.central_body)
                size = 6
            else:
                # Show a warning message if the dependent variable is not known
                list_d_v = ["h", "rho", "V", "m", "F_T", "D", "C_D", "r_cb", "Kep"]
                print("Warning: The dependent variable '%s' is not in known list:" % dep_key, list_d_v)
                size = 0
            if size != 0:
                # Add the dependent variable to the list
                self.dependent_variables_to_save.append(d_v)
                # Save the index at which this dependent variable can be save
                self.dep_var_loc[dep_key] = [idx, size]
                # Increment the index by the size of the dependent variable
                idx += size

    def create_propagator(self, propagator=propagation_setup.propagator.encke, prop_mass=False):
        """
        Create the simulation propagator.
        Inputs:
         * propagator: type of propagator to use (Encke is advised; it propagates the difference with a Keplerian orbit)
         * prop_mass (bool): whether or not to propagate the satellite mass based on thrust or not
        """
        # Define the translational propagator settings
        translation_propagator = propagation_setup.propagator.translational(
            [self.central_body],
            self.acceleration_models,
            [self.sat.name],
            self.initial_state,
            self.termination_settings,
            propagator,
            self.dependent_variables_to_save
        )
        # Also propagate the mass based on the acceleration from the thrust
        if prop_mass:
            # Define the mass propagator settings
            mass_rate_model = {self.sat.name: [propagation_setup.mass.from_thrust()]}
            mass_propagator = propagation_setup.propagator.mass(
                [self.sat.name],
                mass_rate_model,
                np.array([self.sat.wet_mass]),
                self.termination_settings)
            # Define the full propagator settings
            self.propagator_settings = propagation_setup.propagator.multitype(
                [translation_propagator, mass_propagator], self.termination_settings, self.dependent_variables_to_save)
            # The following 3 lines are needed to properly propagate the mass
            translation_propagator.recreate_state_derivative_models(self.bodies)
            translation_propagator.reset_and_recreate_acceleration_models({self.sat.name: self.accelerations}, self.bodies)
            self.propagator_settings.recreate_state_derivative_models(self.bodies)
        # Only propagate the vehicle translationally
        else:
            self.propagator_settings = translation_propagator

    def simulate(self):
        pass

test = True
if test:
    from utils import sat_models as SM
    OS = orbit_simulation(SM.satellites["CS_1021"], "Mars", 2*constants.JULIAN_DAY)
    OS.create_bodies()
    OS.create_initial_state()
    OS.create_accelerations(default_config=1, thrust=1)
    OS.create_integrator()
    OS.create_termination_settings()
    OS.create_dependent_variables(to_save=["m", "V", "h", "rho", "Kep", "F_T"])
    OS.create_propagator(prop_mass=True)
