import sys
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
import numpy as np
from tudatpy.kernel.numerical_simulation import propagation_setup
from tudatpy.kernel.interface import spice_interface
from tools import mission_geometry as MG
from utils.thrust_models import BHT_100
from utils.thrust_models import muNRIT_25

# Astronomical Unit distance in [m]
AU = 149597890e3

# List to save values of battery capacity used, and time for which the battery was used
capacity_taken = []
dt_capacity = []

class thrust_model:

    def __init__(self, orbit_sim, thrust_mod=None, solar_constant=1366, I_sp=800, use_battery=True, ionisation_eff=1):
        """
        Satellite thrust model
        Inputs:
         * orbit_sim (utils.propagation.orbit_simulation): orbit_sim class of the satellite orbit
         * thrust_mod (int): thrust model to use
           * None: no thrust
           * 0: constant thrust of 1 mN when power above 10 W and at any density
           * 1: thrust based on the BHT-100 hall thrusters (with tank), interpolated from power, on when power above 107 W, at any density
           * 2: thrust based on the μNRIT2.5 radiofrequency ion thruster (with tank), interpolated from power, on when power above 13.1 W, at any density
           * 3: thrust based on the μNRIT2.5 radiofrequency ion thruster (without tank), interpolated from power and mass flow, on when power above 13.1 W or mass flow at engine above 1.45E-08 kg/s
        """
        # Save the relevant variables in the class
        self.os_sim = orbit_sim
        self.vehicle = self.os_sim.bodies.get_body(self.os_sim.sat.name)
        self.init_time = self.os_sim.init_time
        self.sat = self.os_sim.sat
        self.central_body = self.os_sim.bodies.get_body(self.os_sim.central_body)
        self.sun = self.os_sim.bodies.get_body("Sun")
        self.I_sp = I_sp
        self.solar_constant = solar_constant
        self.thrust_mod = thrust_mod
        self.ionisation_eff = ionisation_eff # Ionisation efficiency
        self.m_treshold = 0
        self.use_battery = use_battery
        self.last_thrust_call = None
        self.power_treshold = [0, np.inf]
        # Select the thrust model (and associated operating conditions)
        if thrust_mod == 0:
            self.power_treshold = [10, 100]
        elif thrust_mod == 1:
            self.power_treshold = [107, 158]
        elif thrust_mod == 2:
            self.power_treshold = [13.1, 34.4]
        elif thrust_mod == 3:
            self.power_treshold = [13.1, 34.4]
            self.m_treshold = 1.456e-8
    
    def magnitude(self, time):
        global capacity_taken
        # If there is no more propellant, return 0
        if self.sat.dry_mass is not None and self.vehicle.mass < self.sat.dry_mass:
            return 0
        # Return a constant thrust magnitude
        if self.thrust_mod == 0:
            return 0.001
        # Return thrust from BHT 100 thruster, based on power
        elif self.thrust_mod == 1:
            # If there is more than a certain power available, keep it for the other systems (thus use of min())
            self.thrust, self.m_flow, self.I_sp = BHT_100.from_power(min(self.power, self.power_treshold[1]))
            if self.thrust == 0:
                thrust, power_use = 0, 0
            else:
                thrust, power_use = self.thrust, self.power
        # Return thrust from μNRIT 2.5 thruster, based on power
        elif self.thrust_mod == 2:
            self.thrust, self.m_flow, self.I_sp = muNRIT_25.from_power(min(self.power, self.power_treshold[1]))
            if self.thrust == 0:
                thrust, power_use = 0, 0
            else:
                thrust, power_use = self.thrust, self.power
        elif self.thrust_mod == 3:
            self.thrust, power, self.m_flow, self.Isp = muNRIT_25.from_P_m(self.power, self.m_flow_t)
            if self.thrust == 0:
                thrust, power_use = 0, 0
            else:
                thrust, power_use = self.thrust * self.ionisation_eff, power
        else:
            raise NotImplementedError("The thrust model %i has not been implemented." % self.thrust_mod)
            
        if self.use_battery:
            # Save the extra power that can be put into the battery or the energy taken from the battery
            if self.power_is_from_battery:
                # Scale the battery capacity used by the actual power use
                capacity_taken[-1] = power_use/self.power * capacity_taken[-1]
                self.sat.power_to_battery = 0
            else:
                # Extra power from the solar panel is available power minus power used for thrust
                self.sat.power_to_battery = max(self.power - power_use, 0)
        if self.os_sim.save_thrust:
            self.os_sim.thrusts[time] = thrust
        return thrust

    def specific_impulse(self, time):
        # Return the specific impulse
        return self.I_sp

    def is_thrust_on(self, time):
        """
        Define whether the engine is on or not
        """
        # Save the time of the first call to thrust
        if self.last_thrust_call is None:
            self.last_thrust_call = time
        # Engine is on if the atmospheric mass flow at engine inlet is above given treshold
        density_fs = self.vehicle.flight_conditions.density
        velocity_fs = self.vehicle.flight_conditions.airspeed
        altitude_fs = self.vehicle.flight_conditions.altitude
        comp_ratio = self.sat.get_comp_ratio(altitude_fs)
        cd = self.sat.get_cd(altitude_fs)
        self.m_flow_t = comp_ratio * density_fs * velocity_fs * self.sat.S_t
        m_flow_ok = self.m_flow_t >= self.m_treshold
        # Engine is on if the solar irradiance is above a given treshold
        power_ok = self.power_available(time) > self.power_treshold[0]
        if not m_flow_ok and power_ok and self.use_battery:
            # Charge the battery
            self.sat.power_to_battery = self.power
        if not (m_flow_ok and power_ok) and self.os_sim.save_thrust:
            self.os_sim.thrusts[time] = 0
        self.os_sim.drags[time] = 0.5*density_fs*velocity_fs**2*cd*self.sat.S_ref
        # Do this last to make sure that the power is computed (inefficient, but gets data needed)
        if self.thrust_mod is None:
            self.os_sim.thrusts[time] = 0
            return False
        return m_flow_ok and power_ok

    def power_available(self, time):
        """
        Compute the maximum available power from the solar panels
        """
        global dt_capacity, capacity_taken
        # Compute the solar irradiance (including shadowing of Mars)
        self.solar_irradiance(time)
        # Extract the state of the satellite and of the Sun
        sat_state = self.vehicle.state
        sun_state = self.sun.state
        # Compute the power available from the solar panels
        self.power = MG.sat_power(sat_state, sun_state, self.irradiance, [self.sat.area_x, self.sat.area_y, self.sat.area_z]) * self.sat.SA_eff * self.sat.EPS_eff
        # If specified, save the (solar) power
        if self.os_sim.save_power:
            self.os_sim.power_dict[time] = self.power
        if self.use_battery:
            self.power_is_from_battery = False
            # Compute the power available from the battery if there is no power from the SA and the battery is charged enough
            if self.power == 0 and self.sat.battery_capacity > self.sat.battery_total_capacity * self.sat.keep_battery_frac:
                # Compute the time during which the battery would be used
                dt = time - self.last_thrust_call
                dt_hr = dt / 3600
                # Compute the maximum battery capacity that would be used
                capacity_to_use = self.power_treshold[1]*dt_hr
                # Use the battery if it has enough capacity
                if capacity_to_use > 0 and self.sat.battery_capacity > capacity_to_use:
                    dt_capacity.append(dt_hr), capacity_taken.append(capacity_to_use)
                    self.power_is_from_battery = True
                    self.power = self.power_treshold[1] # Set the available power to the maximum one that the thrust would use
        self.last_thrust_call = time
        return self.power

    def solar_irradiance(self, time):
        """
        Compute the solar irradiance at the satellite position
        """
        # Get the position and radius of the required bodies
        sun_pos = self.sun.state[:3]
        cb_pos = self.central_body.state[:3]
        sat_pos = self.vehicle.state[:3]
        cb_R = spice_interface.get_average_radius(self.os_sim.central_body)
        Sun_R = spice_interface.get_average_radius("Sun")
        # Compute the shadow due to Mars between the sat and the Sun
        shadow = MG.shadow_function(sun_pos, Sun_R, cb_pos, cb_R, sat_pos)
        # Compute the distance between the sun and the satellite
        sun_dist = np.linalg.norm(sun_pos-sat_pos)
        # Compute the solar irradiance at the satellite position, based on the shadow number, and the solar constant on Earth
        self.irradiance = self.solar_constant / (sun_dist/AU)**2 * shadow
        # If specified, save the solar irradiance
        if self.os_sim.save_power:
            self.os_sim.solar_irradiances[time] = self.irradiance
        return self.irradiance

def thrust_settings(propagation, thrust_mod, ionisation_eff=1, use_battery=True):
    # Define the thrust guidance function
    thrust_guidance = thrust_model(propagation, thrust_mod, ionisation_eff=ionisation_eff, use_battery=use_battery)
    # Define the thrust settings (direction and magnitude)
    thrust_direction_s = propagation_setup.thrust.thrust_direction_from_state_guidance(
        central_body=propagation.central_body, is_colinear_with_velocity=True, direction_is_opposite_to_vector=False)
    thrust_magnitude_s = propagation_setup.thrust.custom_thrust_magnitude(
        thrust_guidance.magnitude, thrust_guidance.specific_impulse, thrust_guidance.is_thrust_on)
    # Return the acceleration for Tudat
    return propagation_setup.acceleration.thrust_from_direction_and_magnitude(thrust_direction_s, thrust_magnitude_s)