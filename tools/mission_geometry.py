import numpy as np
from numpy.linalg import norm

def shadow_function(occulted_pos, occulted_r, occulting_pos, occulting_r, sat_pos):
    """
    Shadow function, to compute the satellite is in the sun light or not
    Inputs:
     * occulted_pos: [float]*3: position of the occulted body (e.g. the Sun), in Cartesian coordinates
     * occulted_r: float: radius of the occulted body
     * occulting_pos: [float]*3: position of the occulting body (e.g. the Moon), in Cartesian coordinates
     * occulting_r: float: radius of the occulting body
     * sat_pos: [float]*3: position of the satellite in Cartesian coordinates
    Output:
     * shadow number: float between 1 (satellite fully exposed to Sun), and 0 (satellite fully in shadow)
    Function inspired by Section 3.4 from   
    """
    # Make sure the inputs are numpy arrays
    r_sun, r_body, r = np.array(occulted_pos), np.array(occulting_pos), np.array(sat_pos)
    R_sun, R_body = occulted_r, occulting_r
    # Position of the Sun w.r.y. occulting body
    s_sun = r_sun - r_body
    # Position of the Satellite w.r.y. occulting body
    s = r - r_body
    # Apparent radius of occulted body
    a = np.arcsin(R_sun / norm(r_sun - r))
    # Apparent radius of occulting body
    b = np.arcsin(R_body / norm(s))
    # Apparent separation between the centers of both bodies
    c = np.arccos( (-s.T @ (r_sun - r)) / (norm(s) * norm(r_sun - r)) )
    # Check for partial occultation
    if np.fabs(a - b) < c and c < (a + b):
        x = (c**2 + a**2 - b**2) / (2*c)
        y = np.sqrt(a**2 - x**2)
        # Compute the occulted are
        A = a**2 * np.arccos(x/a) + b**2 * np.arccos((c-x)/b) - c*y
        # return the shadow fraction value
        return 1 - A / (np.pi*a**2)
    # Check for total occultation
    elif (c < (a - b) and a > b) or (c < (b - a) and a < b):
        return 0
    return 1

# Tests
tests = False
if tests:
    print("0.0 ?", end="\t")
    sun_d = 149598e6
    sun_r = 6.96e8
    earth_r = 6378.137e3
    print(shadow_function(-sun_d*np.array([1,0,0]), sun_r, np.zeros(3), earth_r, (earth_r+1e6)*np.array([1,0,0])))
    print("1.0 ?", end="\t")
    print(shadow_function(-sun_d*np.array([1,0,0]), sun_r, np.zeros(3), earth_r, (earth_r+1e6)*np.array([0,1,0])))
    print("0.4547 ?", end="\t")
    dir = np.array([0.018,1,0])
    dir = dir / norm(dir)
    print(shadow_function(-sun_d*np.array([1,0,0]), sun_r, np.zeros(3), earth_r, (earth_r+1e3)*dir))