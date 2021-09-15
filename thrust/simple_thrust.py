import numpy as np
import tudatpy
from tudatpy.kernel.simulation import propagation_setup
import sys
sys.path.insert(0,"\\".join(sys.path[0].split("\\")[:-2]))
from MCD import parallel_mcd as PMCD

class thrust_model:

    def __init__(self, vehicle, init_time=0, Isp_base=800, dens_treshold=1e-13):
        self.vehicle = vehicle
        self.init_time = init_time
        self.Isp_base = Isp_base
        self.dens_treshold = dens_treshold
    
    def magnitude(self, time):
        return 0.001

    def specific_impulse(self, time):
        return self.Isp_base

    def is_thrust_on(self, time):
        if self.is_MCD_loaded:
            # Engine is on if the air density is above a given treshold
            return self.vehicle.get_flight_conditions().current_density > self.dens_treshold
        return True
    
    def is_MCD_loaded(self):
        # Check if the MCD is loaded.
        return PMCD.LAST_RESULTS is not None or len(PMCD.MODULE_LIST) == 0

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