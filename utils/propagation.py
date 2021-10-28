import numpy as np
import time
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel import constants
from tudatpy.kernel.numerical_simulation import environment_setup
from tudatpy.kernel.numerical_simulation import propagation_setup
from tudatpy.kernel.numerical_simulation import SingleArcSimulator
from tudatpy.kernel.astro import element_conversion
import sys
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
from tools import time_conversions as TC
from utils import thrust as T
# Load the SPICE kernel
if spice_interface.get_total_count_of_kernels_loaded() == 0:
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
    
    def __init__(self, sat, central_body, sim_time, init_time=TC.MCD_to_Tudat(2459942), verbose=False, save_power=False):
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
        self.verbose = verbose
        self.save_power = save_power
        # Solar irradiance dict
        self.solar_irradiances = dict()
        # Power dict
        self.power_dict = dict()

    def create_bodies(self, additional_bodies=["Sun", "Jupiter"], use_MCD=[False, False], preload_MCD=False, save_MCD_vals=False):
        """
        Create the simulation bodies.
        Inputs:
         * additional_bodies ([string]): bodies to create in addition to the central body (by default, the two most massives in the solar system)
         * use_MCD ([bool, bool]): first boolean to indicate wether the MCD atmosphere model should be implemented. If True, the second boolean is used to indicate whether winds should also be added
         * preload_MCD (bool): if True, load all of the MCD files at once
        """
        bodies_to_create = [self.central_body] + additional_bodies      # Bodies that will be created and used in the simulation
        global_frame_origin = self.central_body                         # Body at the centre of the simulation
        # Setup TUDAT body settings
        body_settings = environment_setup.get_default_body_settings(
            bodies_to_create,
            global_frame_origin,
            base_frame_orientation="ECLIPJ2000"
        )

        if use_MCD[0] and self.central_body == "Mars":
            # Use the atmospheric model from the Mars Climate Database
            from MCD.parallel_mcd import parallel_mcd as PMCD
            mcd = PMCD(load_on_init=preload_MCD, save_all_vals=save_MCD_vals)
            body_settings.get(self.central_body).atmosphere_settings = environment_setup.atmosphere.custom_four_dimensional_constant_temperature(
                mcd.density, constant_temperature=210, specific_gas_constant=192, ratio_of_specific_heats=1.3)
            # Values taken from https://meteor.geol.iastate.edu/classes/mt452/Class_Discussion/Mars-physical_and_orbital_statistics.pdf

            if use_MCD[1]:
                # Add winds from the MCD. Only possible if the MCD atmospheric model is used
                body_settings.get(self.central_body).atmosphere_settings.wind_settings = environment_setup.atmosphere.custom_wind_model(mcd.wind)
        else:
            if use_MCD[0] and self.verbose:
                print("Warning: the MCD should only be used if the central body is Mars")
            if use_MCD[1] and self.verbose:
                print("Warning: the MCD winds can only be added if the MCD atmosphere is used as well")
            # Use an exponential atmosphere model
            if self.central_body == "Mars":
                # Exponential parameters taken from http://link.springer.com/content/pdf/10.1007%2F978-3-540-73647-9_3.pdf
                density_scale_height, density_at_zero_altitude = 7295, 0.0525
                body_settings.get(self.central_body).atmosphere_settings = environment_setup.atmosphere.exponential(
                    density_scale_height, density_at_zero_altitude) 
            else:
                print("Warning, no atmosphere model not setup for %s" % self.central_body)
    	
        # Create the system of bodies
        self.bodies = environment_setup.create_system_of_bodies(body_settings)
        
        self.bodies.create_empty_body(self.sat.name)                                # Add satellite body
        self.bodies.get_body(self.sat.name).set_constant_mass(self.sat.wet_mass)    # Add mass to the satellite

        # Define the aerodynamic coefficient settings
        def force_coefficients(_):
            # Get the altitude from the flight conditions
            h = self.bodies.get_body(self.sat.name).flight_conditions.altitude
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

        # Get the gravitational parameter and radius of the central body
        self.mu = self.bodies.get_body(self.central_body).gravitational_parameter
        self.R_cb = spice_interface.get_average_radius(self.central_body)

    def create_initial_state(self, a=None, h_p=150e3, e=0, i=0, omega=0, Omega=0, theta=0):
        """
        Define the initial state of the Satellite.
        Inputs:
         * a (float): semi-major axis of the orbit in meters
         * h (float): altitude of the satellite periapsis above the central body, in meters
         * e (float): eccentricity of the orbit [0-1]
         * i (float): inclination of the orbit in radians [0-pi/2]
         * omega (float): argument of periapsis of the orbit in radians [0-pi]
         * Omega (float): longitude of the ascending noce in radians [0-2pi]
         * theta (float): true anomaly in radians (note: this defines where the satellite starts in its orbit, and can most of the time be ignored) [0-2pi]
        """
        # Get the semi-major axis from the periapsis altitude and eccentricity (if it is not specified)
        if a is None:
            a = (self.R_cb + h_p) / (1 - e)
        # Convert the Keplerian orbit to an initial cartesian state
        self.initial_state = element_conversion.keplerian_to_cartesian_elementwise(
            gravitational_parameter = self.mu,
            semi_major_axis = a,
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
            if len(env_accelerations) != 0 and self.verbose:
                # If accelerations were also provided, inform that they will be overwritten
                print("Warning: the provided environmental accelerations have been overwritten")
            if default_config == 0:
                # Central body point mass and aerodynamics
                env_accelerations = [env_acceleration(self.central_body, PM=True, aero=True)]
            elif default_config == 1:
                # Central body spherical harmonics of degree/order 4 and aerodynamics, Solar radiation
                env_accelerations = [
                    env_acceleration(self.central_body, SH=True, SH_do=[4,4], aero=True),
                    env_acceleration("Sun", rad=True)
                    ]
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
        initial_time_step = dt[1]               # seconds
        coefficient_set = propagation_setup.integrator.RKCoefficientSets.rkdp_87 
        minimum_step_size = dt[0]               # seconds
        maximum_step_size = dt[2]               # seconds
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
            maximum_factor_increase = sf[2],
            minimum_factor_increase = sf[1] )

    def create_termination_settings(self, min_altitude=50e3, max_altitude=500e3, cpu_time=60*60):
        """
        Create the termination settings; conditions for when the simulation must stop.
        Inputs:
         * min_altitude (float): altitude below which propagation stops, in meters (set to 0 to remove)
         * max_altitude (float): periapsis altitude above which propagation stops, in meters (set to Inf to remove)
         * cpu_time (float): CPU time after which the propagation is stopped, in seconds (set to Inf to remove)
        """
        # Create a lower altitude termination setting (when sat is too deep in the atmosphere or crashed)
        termination_altitude = propagation_setup.dependent_variable.altitude(self.sat.name, self.central_body)
        termination_setting_min_altitude = propagation_setup.propagator.dependent_variable_termination(
                dependent_variable_settings = termination_altitude, limit_value = min_altitude, use_as_lower_limit = True)
        termination_periapsis = propagation_setup.dependent_variable.periapsis_altitude(self.sat.name, self.central_body)
        # Create an upper altitude termination setting (when sat periapsis is too high above the planet)
        termination_setting_max_altitude = propagation_setup.propagator.dependent_variable_termination(
                dependent_variable_settings = termination_periapsis, limit_value = max_altitude, use_as_lower_limit = False)
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
           * Kep:   Keplerian state of the satellite (a, e, i, omega, Omega, theta)
           * h_p:   altitude of the periapsis of the satellite
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
                d_v = propagation_setup.dependent_variable.body_fixed_airspeed_velocity(self.sat.name, self.central_body)
                size = 3
            elif dep_key == "m":
                # Mass of the satellite
                d_v = propagation_setup.dependent_variable.body_mass(self.sat.name)
                size = 1
            elif dep_key == "F_T":
                # Acceleration due to the Thrust
                d_v = propagation_setup.dependent_variable.single_acceleration(
                    propagation_setup.acceleration.thrust_acceleration_type, self.sat.name, self.sat.name)
                size = 3
            elif dep_key == "D":
                # Acceleration due to the atmosphere
                d_v = propagation_setup.dependent_variable.single_acceleration(
                    propagation_setup.acceleration.aerodynamic_type, self.sat.name, self.central_body)
                size = 3
            elif dep_key == "C_D":
                # Aerodynamic coefficients
                d_v = propagation_setup.dependent_variable.aerodynamic_force_coefficients(self.sat.name)
                size = 3
            elif dep_key == "r_cb":
                # Relative position of the central body w.r.t. the sun
                d_v = propagation_setup.dependent_variable.relative_position(self.sat.name, self.central_body)
                size = 3
            elif dep_key == "Kep":
                # Keplerian state of the satellite
                d_v = propagation_setup.dependent_variable.keplerian_state(self.sat.name, self.central_body)
                size = 6
            elif dep_key == "h_p":
                # Altitude of the periapsis of the satellite
                d_v = propagation_setup.dependent_variable.periapsis_altitude(self.sat.name, self.central_body)
                size = 1
            else:
                # Show a warning message if the dependent variable is not known
                list_d_v = ["h", "rho", "V", "m", "F_T", "D", "C_D", "r_cb", "Kep"]
                if self.verbose:
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
        translation_propagator_settings = propagation_setup.propagator.translational(
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
            mass_rate_settings = {self.sat.name: [propagation_setup.mass_rate.from_thrust()]}
            mass_rate_model = propagation_setup.create_mass_rate_models(self.bodies, mass_rate_settings, self.acceleration_models)
            mass_propagator_settings = propagation_setup.propagator.mass(
                [self.sat.name],
                mass_rate_model,
                [self.sat.wet_mass],
                self.termination_settings)
            # Define the full propagator settings
            self.propagator_settings = propagation_setup.propagator.multitype(
                [translation_propagator_settings, mass_propagator_settings], self.termination_settings, self.dependent_variables_to_save)
        # Only propagate the vehicle translationally
        else:
            self.propagator_settings = translation_propagator_settings

    def simulate(self):
        # If the verbose is set to True, time the simulation run
        if self.verbose:
            print("Starting simulation...")
            t0 = time.time()

        # Run the simulation
        dynamics_simulator = SingleArcSimulator(
            self.bodies, self.integrator_settings, self.propagator_settings, print_dependent_variable_data=self.verbose
        )

        # Extract the states and time
        self.states = np.vstack(list(dynamics_simulator.state_history.values()))
        self.sim_times = np.array(list(dynamics_simulator.state_history.keys()))
        self.sim_times -= self.sim_times[0]
        # Extract the dependent variables
        self.dep_vars = np.vstack(list(dynamics_simulator.dependent_variable_history.values()))

        # If the verbose is set to True, show simulation run time
        if self.verbose:
            cpu_time = time.time() - t0
            print("Simulation took %.2f seconds." % cpu_time)

        return self.sim_times, self.states, self.dep_vars

    def get_dep_var(self, key):
        # Get the dependent variable with the given key
        dep_var_idx, dep_var_size = self.dep_var_loc[key]
        dv = self.dep_vars[:,dep_var_idx:dep_var_idx+dep_var_size]
        if dv.shape[1] == 1:
            return dv[:,0]
        return dv

test = False
if test:
    # As a test, simulate a satellite for 2 days around Mars
    from utils import sat_models as SM
    # Create a satellite with a ballistic coefficient of 50kg/m2 (Cd=2.5, m=5kg, S_ref=0.04m2)
    orbiter = SM.satellite("Orbiter", 5, 2.5, 0.04)
    OS = orbit_simulation(orbiter, "Mars", 2*constants.JULIAN_DAY)
    # Do not use the MCD for atmospheric properties, nor winds
    OS.create_bodies(use_MCD=[False, False])
    # Start from an altitude of 180km
    OS.create_initial_state(h=140e3)
    # Load the accelerations from default config 1: Central body spherical harmonics of degree/order 4 and aerodynamics, Solar radiation
    OS.create_accelerations(default_config=1, thrust=0)
    OS.create_integrator()
    OS.create_termination_settings()
    OS.create_dependent_variables(to_save=["V", "h", "rho"])
    OS.create_propagator(prop_mass=False)
    OS.simulate()
    print("Simulation ended at %.2f hrs" % (OS.sim_times[-1]/3600))
    print("Final state:", OS.states[-1])
    altitudes = OS.get_dep_var("h")
    airspeed = OS.get_dep_var("V")
    density = OS.get_dep_var("rho")
    print(density[0])
    print("Final altitude of %.2f km" % (altitudes[-1]/1e3))
    print("Final airspeed of %.2f m/s" % np.linalg.norm(airspeed[-1]))