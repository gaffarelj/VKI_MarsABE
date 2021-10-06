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
    # Orbital simulation class, containing all the code required for setup, and simulation run
    def __init__(self, sat, central_body, sim_time, init_time=TC.MCD_to_Tudat(2459942)):
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
        accelerations = dict()
        for env_a in env_accelerations:
            accelerations[env_a.body_name] = env_a.a
        # Add the thrust as an acceleration
        if thrust is not None:
            accelerations[self.sat.name] =  [T.thrust_settings(self, thrust)]
        # Create the acceleration models
        self.acceleration_models = propagation_setup.create_acceleration_models(
            self.bodies,
            {self.sat.name: accelerations},
            [self.sat.name],
            [self.central_body]
        )



from utils import sat_models as SM

OS = orbit_simulation(SM.satellites["CS_1021"], "Mars", 2*constants.JULIAN_DAY)
OS.create_bodies()
OS.create_initial_state()
OS.create_accelerations(default_config=1, thrust=1)