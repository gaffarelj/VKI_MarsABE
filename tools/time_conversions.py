import numpy as np
from scipy.optimize import fsolve

# This value is the time in seconds since J2000 and corresponds to January 1st 2023 at midnight (MCD JD_ref)
MCD_0_TUDAT = 725330764.8000107
JULIAN_DAY = 86400

def JD_to_Ls(JD, JD_ref=2459940.032, JD_Tudat=False):
    """
    This functions converts a Julian Date to a Solar Longitude for Mars.
    Input:
     * JD (float): Julian Date to be converted
     * JD_ref (float): optional, Julian date at which Ls=0deg.
     * JD_Tudat (bool): optional, set to True if the date comes from Tudat
    Output:
     * Ls (float): Solar Longitude in [deg], between 0 and 360.
    This function has been mostly taken from MARS CLIMATE DATABASE v5.3
    DETAILED DESIGN DOCUMENT, p.33
    """
    # If the Julian Date is from Tudat (in seconds from J2000), convert it in the MCD format (days since J2023)
    if JD_Tudat:
        JD += MCD_0_TUDAT
        JD /= JULIAN_DAY
    # Compute the Martian sol number Ds
    Ds = (JD - JD_ref) * 86400/88775.245 % 668.6
    # Compute the mean anomaly
    M = 2*np.pi*(Ds-485.35)/668.6
    # Compute the eccentric anomaly
    e = 0.09340 # eccentricity
    E = fsolve(lambda E: E - e*np.sin(E) - M, np.pi)[0]
    # Compute the true anomaly
    theta = 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2))
    # Compute Ls
    Ls = (theta*180/np.pi + 250.99) % 360
    return Ls, Ds