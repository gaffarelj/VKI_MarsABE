import numpy as np
from numpy.lib.npyio import save
from tudatpy.kernel.simulation import propagation_setup
from tudatpy.kernel.interface import spice_interface
import sys
sys.path.insert(0,"\\".join(sys.path[0].split("\\")[:-2]))
from tools import mission_geometry as MG

# Astronomical Unit distance in [m]
AU = 149597890e3
# Solar irradiance dict
solar_irradiances = dict()
# Power dict
power_dict = dict()
# Satellite areas in x, y, z planes (see SPARTA README)
sat_areas = {
    "CS_0020": [0, 0.042426, 0.042426],
    "CS_1020": [0, 0.102426, 0.042426],
    "CS_0021": [0.031058, 0.083343, 0.083343],
    "CS_2020": [0, 0.162426, 0.042426],
    "CS_1021": [0.031058, 0.143343, 0.083343],
    "CS_3020": [0, 0.222426, 0.042426],
    "CS_2021": [0.031058, 0.203343, 0.083343],
    "CS_2120": [0, 0.282426, 0.042426],
    "CS_3021": [0.031058, 0.263343, 0.083343]
}

class thrust_model:

    def __init__(self, bodies, vehicle_name, init_time=0, Isp_base=800, dens_treshold=1e-13, \
        save_power=False, solar_constant=1366, power_treshold=10, sat_name="CS_0020", power_eff=0.29*0.93):
        self.bodies = bodies
        self.vehicle = bodies.get_body(vehicle_name)
        self.init_time = init_time
        self.Isp_base = Isp_base
        self.dens_treshold = dens_treshold
        self.central_body = self.bodies.get_body("Mars")
        self.sun = self.bodies.get_body("Sun")
        self.save_power = save_power
        self.solar_constant = solar_constant
        self.power_treshold = power_treshold
        self.sat_area = sat_areas[sat_name]
        self.power_eff = power_eff
    
    def magnitude(self, time):
        # Return a constant thrust magnitude
        return 0.001

    def specific_impulse(self, time):
        # Return a constant specific impulse
        return self.Isp_base

    def is_thrust_on(self, time):
        """
        Define whether the engine is on or not
        """
        # Engine is on if the air density is above a given treshold
        density_ok = self.vehicle.get_flight_conditions().current_density > self.dens_treshold
        # Engine is on if the solar irradiance is above a given treshold
        power_ok = self.power_available(time) > self.power_treshold
        return density_ok and power_ok

    def power_available(self, time):
        """
        Compute the maximum available power from the solar panels
        """
        # Compute the solar irradiance (including shadowing of Mars)
        self.solar_irradiance(time)
        # Extract the state of the satellite and of the Sun
        sat_state = self.vehicle.state
        sun_state = self.sun.state
        # Compute the power available from the solar panels
        self.power = MG.sat_power(sat_state, sun_state, self.irradiance, self.sat_area)
        # If specified, save the power
        if self.save_power:
            power_dict[time] = self.power
        return self.power

    def solar_irradiance(self, time):
        """
        Compute the solar irradiance at the satellite position
        """
        # Get the position and radius of the required bodies
        sun_pos = self.sun.state[:3]
        Mars_pos = self.central_body.state[:3]
        sat_pos = self.vehicle.state[:3]
        Mars_R = spice_interface.get_average_radius("Mars")
        Sun_R = spice_interface.get_average_radius("Sun")
        # Compute the shadow due to Mars between the sat and the Sun
        shadow = MG.shadow_function(sun_pos, Sun_R, Mars_pos, Mars_R, sat_pos)
        # Compute the distance between the sun and the satellite
        sun_dist = np.linalg.norm(sun_pos-sat_pos)
        # Compute the solar irradiance at the satellite position, based on the shadow number, and the solar constant on Earth
        self.irradiance = self.solar_constant / (sun_dist/AU)**2 * shadow
        # If specified, save the solar irradiance
        if self.save_power:
            solar_irradiances[time] = self.irradiance
        return self.irradiance

def thrust_settings(bodies, init_time=0, save_power=False):
    # Define the thrust guidance function
    thrust_guidance = thrust_model(bodies, "Satellite", init_time, save_power=save_power)
    # Define the thrust settings (direction and magnitude)
    thrust_direction_s = propagation_setup.acceleration.thrust_direction_from_state_guidance(
        central_body="Mars", is_colinear_with_velocity=True, direction_is_opposite_to_vector=False)
    thrust_magnitude_s = propagation_setup.acceleration.custom_thrust_magnitude(
        thrust_guidance.magnitude, thrust_guidance.specific_impulse, thrust_guidance.is_thrust_on)
    # Return the acceleration for Tudat
    return propagation_setup.acceleration.thrust_acceleration(thrust_direction_s, thrust_magnitude_s)