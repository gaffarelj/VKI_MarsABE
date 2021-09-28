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

def sat_power(sat_state, sun_state, sun_irradiance, sat_area, test=False):
    """
    Satellite power function, to compute the total power available by the satellite.
    Inputs:
     * sat_vel: [float]*3: velocity of the satellite in Cartesian coordinates (to compute its orientation)
     * sun_pos: [float]*3: position of the Sun in Cartesian coordinates
     * sun_irradiance: float: solar irradiance on the satellite (in W/m2)
     * sat_area: [float]*3: effective area of the solar panels in the x, y, z planes (see SPARTA README)
    """
    # Do no bother computing the power if the irradiance is 0 W/m2
    if sun_irradiance == 0:
        return 0
    # Get velocity vector of the satellite and of the Sun w.r.t. Mars
    sat_vel = sat_state[3:]
    # Compute the heading of the sat w.r.t. Mars, from its velocity
    sat_heading = np.pi/2 if sat_vel[0] == 0 else np.arctan(sat_vel[1]/sat_vel[0])
    # Compute the heading of the Sun w.r.t. Mars from its position
    sun_heading = 0 if sun_state[1] == 0 else np.arctan(-sun_state[0]/sun_state[1])
    heading = sat_heading - sun_heading
    # Compute the angle of attack of the satellite
    AoA = np.arctan(sat_vel[2]/np.sqrt(sat_vel[0]**2 + sat_vel[1]**2))
    # Compute scaling factor in each plane
    power_scale = [np.fabs(np.sin(heading)), np.fabs(np.cos(heading) * np.cos(AoA)), np.fabs(np.sin(AoA))]
    # Compute the effective solar panel surfaces
    eff_surf = np.array(sat_area) * np.array(power_scale)
    # Compute the satellite power
    power = sum(eff_surf) * sun_irradiance

    if test:
        # Test the Solar power computation
        print("Sat heading of %.2f deg, sun heading of %.2f deg"% (np.rad2deg(sat_heading), np.rad2deg(sun_heading)))
        print("Heading of %.2f deg, angle of attack of %.2f deg" % (np.rad2deg(heading), np.rad2deg(AoA)))
        print("Power scaled by %.4f x / %.4f y / %.4f z" % tuple(power_scale))
        
        do_plot = (input("Do plot? y/[n]: ") == "y")
        if do_plot:
            # Imports for plotting
            import matplotlib
            matplotlib.use("pdf")
            import matplotlib.pyplot as plt
            plt.rcParams.update({'font.size': 13, 'figure.figsize': (7, 7), 'savefig.format': 'pdf'})
            import sys
            sys.path.insert(0,"\\".join(sys.path[0].split("\\")[:-1]))
            
            # Plot the position of the Sun and the Satellite
            fig, ax = plt.subplots()
            ax.scatter(sat_state[0], sat_state[1], label="Sat")
            ax.scatter(sun_state[0], sun_state[1], label="Sun")

            # Plot the satellite direction, and the light rays directions
            sat_vel_unit = sat_vel / np.linalg.norm(sat_vel)
            sun_pos_unit = -sun_state[:3] / np.linalg.norm(-sun_state[:3])
            f = 5e10
            ax.plot([sat_state[0], sat_state[0]+sat_vel_unit[0]*f], [sat_state[1], sat_state[1]+sat_vel_unit[1]*f], label="Sat direction")
            ax.plot([sun_state[0], sun_state[0]+sun_pos_unit[0]*f], [sun_state[1], sun_state[1]+sun_pos_unit[1]*f], label="Sunlight direction")

            # Save the figure
            fig.tight_layout(), ax.grid(), ax.legend()
            ax.set_title("Power scaled by %.4f x / %.4f y / %.4f z" % tuple(power_scale))
            ax.set_ylim([-2.5e11, 0.5e11]), ax.set_xlim([-1.5e11, 1.5e11])
            plt.savefig("figures/test_sun_power.pdf")
            plt.close()

    return power

# Tests
tests = False
if tests:
    # Test the shadow function
    sun_d = 149598e6
    sun_r = 6.96e8
    earth_r = 6378.137e3
    print("0.0 ?", end="\t")
    print(shadow_function(-sun_d*np.array([1,0,0]), sun_r, np.zeros(3), earth_r, (earth_r+1e6)*np.array([1,0,0])))
    print("1.0 ?", end="\t")
    print(shadow_function(-sun_d*np.array([1,0,0]), sun_r, np.zeros(3), earth_r, (earth_r+1e6)*np.array([0,1,0])))
    print("0.4547 ?", end="\t")
    dir = np.array([0.018,1,0])
    dir = dir / norm(dir)
    print(shadow_function(-sun_d*np.array([1,0,0]), sun_r, np.zeros(3), earth_r, (earth_r+1e3)*dir))