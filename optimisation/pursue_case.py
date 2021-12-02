from matplotlib.pyplot import sca
import numpy as np
import sys
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
from tudatpy.kernel import constants
from utils import propagation as P
from utils import sat_models as SM
from tools import plot_utilities as PU

sim_days = 650

## Propagation main parameters
s_name, thrust_model, ionisation_eff, use_battery = "CS_2021", 3, 0.50, True
h_p, h_a = 143e3, 405e3
i, omega, Omega = np.deg2rad(25), 0, 0

satellites = SM.satellites_with_tank if thrust_model == 2 else SM.satellites

# Create the orbital simulation instance
OS = P.orbit_simulation(satellites[s_name], "Mars", sim_days*constants.JULIAN_DAY, save_power=True, save_thrust=True, verbose=True)
# Create the simulation bodies, and use the MCD
OS.create_bodies(use_MCD=[True, False], use_GRAM=False)
# Create the initial state of the satellite
a = OS.R_cb + (h_a+h_p)/2
e = 1 - (OS.R_cb + min(h_p, h_a)) / a       # Use min because h_p could actually be higher than h_a due to the way the problem is setup)
OS.create_initial_state(a=a, e=e, i=i, omega=omega, Omega=Omega)
# Load the accelerations from default config 1: Central body spherical harmonics of degree/order 4 and aerodynamics, Solar radiation
OS.create_accelerations(default_config=2, thrust=thrust_model, ionisation_eff=ionisation_eff, use_battery=use_battery)
# Create the integrator, termination settings, dependent variables, and propagator
OS.create_integrator()
OS.create_termination_settings()
OS.create_dependent_variables(to_save=["h_p", "h", "rho", "m"])
OS.create_propagator(prop_mass=(thrust_model == 2))

# Simulate the satellite in orbit
times, states, _ = OS.simulate()

# Save results to CSV
results = np.array([times, states[:,0], states[:,1], states[:,2]]).T
np.savetxt(sys.path[0]+"/optimisation/results/orbit.csv", results, delimiter=',', header="time, X, Y, Z")

## Plot results
# Battery capacity
battery_times = np.array(list(OS.battery_capacity.keys()))
battery_times_hr = (battery_times-battery_times[0])/3600
battery_capacity = list(OS.battery_capacity.values())
PU.plot_single(battery_times_hr, battery_capacity, "Time [hr]", "Battery capacity [Whr]", "optimisation/pursue_case/battery")

# Power
power_times = np.array(list(OS.power_dict.keys()))
if len(power_times) > 0:
    power_times_hr = (power_times-power_times[0])/3600
    power = list(OS.power_dict.values())
    PU.plot_single(power_times_hr, power, "Time [hr]", "Power [W]", "optimisation/pursue_case/power", scatter=False, markersize=5)

# Solar irradiance
irradiance_times = np.array(list(OS.solar_irradiances.keys()))
if len(irradiance_times) > 0:
    irradiance_times_hr = (irradiance_times-irradiance_times[0])/3600
    irradiances = list(OS.solar_irradiances.values())
    PU.plot_single(irradiance_times_hr, irradiances, "Time [hr]", "Solar irradiance [W/m2]", "optimisation/pursue_case/irradiance", scatter=False, markersize=5)

# Altitude
PU.plot_single(times/3600, np.array(OS.get_dep_var("h"))/1e3, "Time [hr]", "Altitude [km]", "optimisation/pursue_case/h")

# Satellite mass
PU.plot_single(times/3600, OS.get_dep_var("m"), "Time [hr]", "Satellite mass [kg]", "optimisation/pursue_case/mass")

# Density
PU.plot_single(times/3600, OS.get_dep_var("rho"), "Time [hr]", "Density [kg/m3]", "optimisation/pursue_case/rho")

# Drag
drag_times = np.array(list(OS.drags.keys()))
drag_times_hr = (drag_times-drag_times[0])/3600
drags = np.array(list(OS.drags.values()))
PU.plot_single(drag_times_hr, drags, "Time [hr]", "Drag [N]", "optimisation/pursue_case/drag")

# Thrust
thrust_times = np.array(list(OS.thrusts.keys()))
if len(thrust_times) > 0:
    thrust_times_hr = (thrust_times-thrust_times[0])/3600
    thrusts = np.array(list(OS.thrusts.values()))
PU.plot_single(thrust_times_hr, thrusts, "Time [hr]", "Thrust [N]", "optimisation/pursue_case/thrust", scatter=False, markersize=5)

# Print infos
h_m = np.mean(OS.get_dep_var("h"))
print("Mean altitude:", h_m)
h_p_s = OS.get_dep_var("h_p")
print("Periapsis decay:", h_p_s[0] - h_p_s[-1])
power_hist = list(OS.power_dict.values())
print("Mean power:", np.mean(power_hist))
print("Mean T/D:", np.mean(thrusts/drags))
