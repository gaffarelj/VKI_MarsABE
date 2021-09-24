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

class thrust_model:
    def __init__(self, bodies, vehicle_name, init_time=0, Isp_base=800, dens_treshold=1e-13, \
        save_solar=False, solar_constant=1360, power_treshold=100, sat_area=0.125, power_eff=0.35*0.93):
        global solar_irradiances
        self.bodies = bodies
        self.vehicle = bodies.get_body(vehicle_name)
        self.init_time = init_time
        self.Isp_base = Isp_base
        self.dens_treshold = dens_treshold
        self.central_body = self.bodies.get_body("Mars")
        self.sun = self.bodies.get_body("Sun")
        self.save_solar = save_solar
        self.solar_constant = solar_constant
        solar_irradiances = dict()
        self.power_treshold = power_treshold
        self.sat_area = sat_area
        self.power_eff = power_eff
    
    def magnitude(self, time):
        return 0.001

    def specific_impulse(self, time):
        return self.Isp_base

    def is_thrust_on(self, time):
        # Engine is on if the air density is above a given treshold
        density_ok = self.vehicle.get_flight_conditions().current_density > self.dens_treshold
        # Engine is on if the solar irradiance is above a given treshold
        # TODO: convert to a power treshold
        solar_ok = self.power_available(time) > self.power_treshold
        return density_ok and solar_ok

    def power_available(self, time):
        self.solar_irradiance(time)
        sat_state = self.vehicle.state
        sun_state = self.sun.state
        self.power = MG.sat_power(sat_state, sun_state, self.irradiance, self.sat_area)
        print(self.power), input()
        return self.power

    def solar_irradiance(self, time):
        global solar_irradiances
        # Get the position and radius of the required bodies
        sun_pos = self.sun.state[:3]
        Mars_pos = self.central_body.state[:3]
        sat_pos = self.vehicle.state[:3]
        Mars_R = spice_interface.get_average_radius("Mars")
        Sun_R = spice_interface.get_average_radius("Sun")
        # Compute the shadow due to Mars between the sat and the Sun
        shadow = MG.shadow_function(sun_pos, Sun_R, Mars_pos, Mars_R, sat_pos)
        # Compute the solar irradiance at the satellite position
        sun_dist = np.linalg.norm(sun_pos-sat_pos)
        irradiance = self.solar_constant / (sun_dist/AU)**2 * shadow
        if self.save_solar:
            solar_irradiances[time] = irradiance
        self.irradiance = irradiance
        return irradiance

def thrust_settings(bodies, init_time=0, save_solar=False):
    # Define the thrust guidance function
    thrust_guidance = thrust_model(bodies, "Satellite", init_time, save_solar=save_solar)
    # Define the thrust settings
    thrust_direction_s = propagation_setup.acceleration.thrust_direction_from_state_guidance(
        central_body="Mars", is_colinear_with_velocity=True, direction_is_opposite_to_vector=False)
    thrust_magnitude_s = propagation_setup.acceleration.custom_thrust_magnitude(
        thrust_guidance.magnitude, thrust_guidance.specific_impulse, thrust_guidance.is_thrust_on)
    # Return the acceleration for Tudat
    return propagation_setup.acceleration.thrust_acceleration(thrust_direction_s, thrust_magnitude_s)