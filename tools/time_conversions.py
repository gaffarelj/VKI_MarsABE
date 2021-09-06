import numpy as np
from scipy.optimize import fsolve

def JD_to_Ls(JD, JD_ref=2459940.032):
    # Compute the Martian sol number Ds
    DS = (JD - JD_ref) * 86400/88775.245 % 668.6
    # Compute the mean anomaly
    M = 2*np.pi*(DS-485.35)/668.6
    # Compute the eccentric anomaly
    e = 0.09340 # eccentricity
    E = fsolve(lambda E: E - e*np.sin(E) - M, np.pi)[0]
    # Compute the true anomaly
    theta = 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2))
    # Compute Ls
    Ls = (theta*180/np.pi + 250.99) % 360
    return Ls