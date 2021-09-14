import numpy as np
import tudatpy
from tudatpy.kernel.simulation import propagation_setup

class thrust_model:

    def __init__(self, vehicle, init_time=0, Isp_base=800):
        self.vehicle = vehicle
        self.init_time = init_time
        self.Isp_base = Isp_base
    
    def magnitude(self, time):
        return 2

    def specific_impulse(self, time):
        return self.Isp_base

    def is_thrust_on(self, time):
        return True

def thrust_settings(bodies, init_time=0):
    # Define the thrust guidance function
    thrust_guidance = thrust_model(bodies.get_body("Satellite"), init_time)
    # Define the thrust settings
    thrust_direction_s = propagation_setup.acceleration.thrust_direction_from_state_guidance(
        central_body="Mars", is_colinear_with_velocity=True, direction_is_opposite_to_vector=False)
    thrust_magnitude_s = propagation_setup.acceleration.custom_thrust_magnitude(
        thrust_guidance.magnitude, thrust_guidance.specific_impulse, thrust_guidance.is_thrust_on)
    # Return the acceleration for Tudat
    return propagation_setup.acceleration.thrust_acceleration(thrust_direction_s, thrust_magnitude_s)