import sys
sys.path.insert(0,"\\".join(sys.path[0].split("\\")[:-2]))
import numpy as np
import tools.plot_utilities as PU

# Lists containing the characteristics of the BHT 100 Hall thruster (Power [W], mass flow [mg/s], Thrust [mN], Isp [s])
BHT_100_P = [109, 109, 107, 115, 113, 114, 116, 117, 111, 113, 114, 113, 113, 133, 136, 133, 156, 157, 158, 158]
BHT_100_m = [0.663, 0.663, 0.683, 0.688, 0.688, 0.693, 0.693, 0.693, 0.693, 0.693, 0.693, 0.693, 0.693, 0.663, 0.663, 0.663, 0.771, 0.801, 0.786, 0.663]
BHT_100_T = [6.5, 6.57, 6.59, 6.68, 6.82, 6.97, 6.95, 6.91, 7.25, 6.91, 6.43, 6.91, 6.96, 7.29, 7.48, 7.39, 9.11, 9.4, 9.19, 8.29]
BHT_100_I = [1000, 1010, 985, 990, 1011, 1027, 1023, 1017, 1067, 1017, 947, 1017, 1025, 1120, 1150, 1137, 1203, 1197, 1192, 1274]

# PU.plot_4d(BHT_100_P, BHT_100_m, BHT_100_I, BHT_100_T, ["Power [W]", "Mass flow [mg/s]", "Isp [s]", "Thrust [mN]"], "thrust/BHT_100")

# Convert the lists to SI units
BHT_100_P = np.array(BHT_100_P)
BHT_100_m = np.array(BHT_100_m) / 1e6
BHT_100_T = np.array(BHT_100_T) / 1e3
BHT_100_I = np.array(BHT_100_I)


def from_power(power):
    if power < min(BHT_100_P) or power > max(BHT_100_P):
        raise ValueError("The power value (%s) should be in the following range: %s ; %s" % (power, min(BHT_100_P), max(BHT_100_P)))
    closest_idx = np.argsort(abs(BHT_100_P - power))[:2]
    thrust = np.interp(power, BHT_100_P[closest_idx], BHT_100_T[closest_idx])
    m_flow = np.interp(power, BHT_100_P[closest_idx], BHT_100_m[closest_idx])
    Isp = np.interp(power, BHT_100_P[closest_idx], BHT_100_I[closest_idx])
    return thrust, m_flow, Isp