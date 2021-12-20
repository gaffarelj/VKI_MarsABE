import numpy as np
import sys
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
from tools import plot_utilities as PU
from GRAM import call_GRAM
from utils import sat_models as SM
from utils.thrust import muNRIT_25 as T
from scipy.interpolate import CubicSpline as CS

# Define satellite to be used
satellite = SM.satellites["CS_1021"]
S_ref, S_t = satellite.S_ref, satellite.S_t
# Define interpolator of velocity as a function of altitude
velocities = CS([85e3, 115e3, 150e3], [3510.9, 3495.84, 3474.51], extrapolate=False)
# Define altitudes
altitudes = np.arange(85e3, 150.01e3, 1)
GRAM_atm = call_GRAM.GRAM_atmo()
# Loop trough different ionisation efficiencies for the atmopshere
for ion_eff in [0.01, 0.1, 0.25, 0.5]:
    D_s, T_s = [], []
    for h in altitudes:
        # Get atmosphere density for given altitude
        rho = GRAM_atm.get_density(h, 0, 0, 0)
        # Get satellite properties at given altitude
        Cd = satellite.get_cd(h)
        comp_ratio = satellite.get_comp_ratio(h)
        # Get orbital velocity at given altitude
        V = velocities(h)
        # Compute drag and mass flow
        drag = 0.5*rho*V**2*Cd*S_ref
        m_flow = comp_ratio * rho * V * S_t
        # Get thrust from mass flow
        m_flow = min(m_flow, max(T.muNRIT25_m))
        if m_flow > min(T.muNRIT25_m):
            thrust, *_ = T.from_m_flow(m_flow)
            thrust *= ion_eff
        else:
            thrust = 0
        D_s.append(drag), T_s.append(thrust)
    D_s, T_s = np.array(D_s), np.array(T_s)

    # Plot and save thrust/drag ratio
    PU.plot_single(T_s/D_s, altitudes/1e3, "Thrust/Drag [-]", "Altitude [km]", "thrust_drag_%.2f"%ion_eff, \
        xline=1, title="Ionisation efficiency of %i%%" % (ion_eff*100), size=(4, 6))