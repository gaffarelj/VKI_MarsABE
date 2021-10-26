import sys
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
import numpy as np

# Lists containing the characteristics of the muNRIT 2.5 Grid Ion thruster (Power [W], mass flow [sccm], Thrust [microN], Isp [s])
muNRIT25_T = [50, 100, 200, 300, 400, 500, 575]
muNRIT25_P = [13.1, 15.2, 19.2, 21.6, 27, 31.2, 34.4]
muNRIT25_I = [363, 713, 1236, 1735, 2194, 2609, 2861]
muNRIT25_m = [0.148, 0.159, 0.174, 0.186, 0.196, 0.206, 0.216]

# Convert the lists to SI units
p = 101325  # [Pa]
T = 273.15  # [K]
Z = 0.9924  # [-]
R = 8314    # [J/K/kmol]
M = 131.293 # [g/mol]
muNRIT25_T = np.array(muNRIT25_T)/1e6                               # Thrust [N]
muNRIT25_P = np.array(muNRIT25_P)                                   # Power [W]
muNRIT25_I = np.array(muNRIT25_I)                                   # Isp [s]
muNRIT25_m = np.array(muNRIT25_m) * ( 5/3*1e-8 * p*M / (Z*R*T) )    # Mass flow [kg/s]

# Thrust the can be reached for the given power
def from_power(power):
    if power < min(muNRIT25_P) or power > max(muNRIT25_P):
        raise ValueError("The power value (%s) should be in the following range: %s ; %s" % (power, min(muNRIT25_P), max(muNRIT25_P)))
    closest_idx = np.argsort(abs(muNRIT25_P - power))[:2]
    thrust = np.interp(power, muNRIT25_P[closest_idx], muNRIT25_T[closest_idx])
    m_flow = np.interp(power, muNRIT25_P[closest_idx], muNRIT25_m[closest_idx])
    Isp = np.interp(power, muNRIT25_P[closest_idx], muNRIT25_I[closest_idx])
    return thrust, m_flow, Isp

# Thrust that can be reached for the given mass flow
def from_m_flow(m_flow):
    if m_flow < min(muNRIT25_m) or m_flow > max(muNRIT25_m):
        raise ValueError("The power value (%s) should be in the following range: %s ; %s" % (m_flow, min(muNRIT25_m), max(muNRIT25_m)))
    closest_idx = np.argsort(abs(muNRIT25_m - m_flow))[:2]
    thrust = np.interp(m_flow, muNRIT25_m[closest_idx], muNRIT25_T[closest_idx])
    power = np.interp(m_flow, muNRIT25_m[closest_idx], muNRIT25_P[closest_idx])
    Isp = np.interp(m_flow, muNRIT25_m[closest_idx], muNRIT25_I[closest_idx])
    return thrust, power, Isp

# Thrust that can be reached for the given combination of mass flow and power
def from_P_m(power, m_flow):
    # Take the maximum possible mass flow and power if the value specified is above it
    power = min(power, max(muNRIT25_P))
    m_flow = min(m_flow, max(muNRIT25_m))
    # Try to get the thrust from the power
    try:
        thrust_P, m_flow_P, Isp_P = from_power(power)
    except ValueError:
        thrust_P, m_flow_P, Isp_P = 0, 0, 0
    # Try to get the thrust from the mass flow
    try:
        thrust_m, power_m, Isp_m = from_m_flow(m_flow)
    except ValueError:
        thrust_m, power_m, Isp_m = 0, 0, 0
    # Take the minimum thrust that could be obtained
    if thrust_m < thrust_P:
        thrust, power, m_flow, Isp = thrust_m, power_m, m_flow, Isp_m
    else:
        thrust, power, m_flow, Isp =  thrust_P, power, m_flow_P, Isp_P
    # If the thrust is 0, the power and mass flow should also be 0
    if thrust == 0:
        return 0, 0, 0, 0
    # Return the thrust, as well as the power, mass flow, and Isp
    return thrust, power, m_flow, Isp

# TODO: from_P_m is returning the minimum thrust that can be reached from Power or from mass flow.
# This leads to underperformance. A more complex model should be implemented later on