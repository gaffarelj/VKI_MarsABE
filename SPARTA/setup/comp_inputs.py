import numpy as np
import sys
sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
import os

# Define conditions at different orbital altitudes
hs = [85, 115, 150]
rhos = [7.1E-07, 1.8E-08, 1.6E-10]
ps = [2.3E-02, 3.7E-04, 7.1E-06]
Ts = [135, 115, 175]
Vs = [3494.17, 3493.29, 3483.82]
fracs = [
    np.array([0.905, 0.035, 0.025, 0.015, 0.015, 0.005]),
    np.array([0.809, 0.045, 0.035, 0.05, 0.055, 0.005]), # first one should be 0.81 but this causes an error in SPARTA
    np.array([0.42, 0.125, 0.045, 0.15, 0.25, 0.01])
]

save_to_input = True
sat_names = ["CS_0020", "CS_0021", "CS_1020", "CS_1021", "CS_2020", "CS_2021", "CS_2120", "CS_3020"]
L_s = [0.3, 0.589778, 0.341421, 0.589778, 0.541421, 0.589778, 0.6, 0.741421]
for j, s_name in enumerate(sat_names):
    print("\n\n* Satellite", s_name)
    # Loop trough conditions
    for i, h in enumerate(hs):
        print("\nWith conditions at altitude of %i km:" % h)
        print("Velocity is of %.2f m/s, and the species are mixed as follows:" % Vs[i])
        print(fracs[i])
        # Inputs
        rho = rhos[i]   # density [kg/m3]
        p = ps[i]       # pressure [Pa]
        T = Ts[i]       # temperature [K]
        u_s = Vs[i]     # free-stream velocity [m/s]
        L = L_s[j]      # reference length [m] (satellite width)
        h_box = 1.0     # box height [m]
        w_box = 1.0     # box width [m]
        l_box = 1.5     # box length [m]
        # Fraction of each species
        species_frac = fracs[i]
        if round(sum(species_frac), 5) != 1:
            print("Warning, the sum of the species fraction does not add up to 1.")

        # Constants
        # Species name, mass, diameter, frontal area
        species_names = ["CO2", "N2", "Ar", "CO", "O", "O2"]
        species_m = np.array([7.31E-26, 4.65E-26, 6.63E-26, 4.65E-26, 2.65E-26, 5.31E-26])
        species_d = np.array([33e-11, 364e-12, 340e-12, 376e-12, 7.4e-11, 346e-12])
        species_sigma = np.pi*species_d**2
        k_b = 1.38e-23  # Bolzmann constant

        # Compute values
        weighted_m = sum(species_frac * species_m)              # weighted mass of species
        weighted_sigma = sum(species_frac * species_sigma)      # weighted frontal area of species
        nrho = rho / weighted_m                                 # number density [#/m3]
        lambda_f = 1 / (np.sqrt(2) * weighted_sigma * nrho)     # mean free path [m]
        Kn = lambda_f / L                                       # Knudsen number [-]
        T_ps = weighted_m * u_s**2 / (7*k_b)                    # post-shock T (Cv=2.5*R) [K]
        cr_ps = np.sqrt(16*k_b*T_ps / (np.pi*weighted_m))       # post-shock average relative velocity [m/s]\
        nrho_ps = 7.67*nrho                                     # post-shock number density [#/m3]
        lambda_ps = 1 / (np.sqrt(2) * weighted_sigma * nrho_ps) # post-shock mean free path [m]
        nu_ps = weighted_sigma * nrho_ps * cr_ps                # post-shock collision frequency [Hz]
        tau_ps = 1 / nu_ps                                      # post-shock collision time [s]
        dt = tau_ps / 4                                         # time step [s]
        grid_f = lambda_f / 4                                   # grid dimension before shock [m]
        grid_ps = lambda_ps / 4                                 # post-shock grid dimension [m]
        n_real = (nrho + nrho_ps) / 2 * h_box * l_box * w_box   # real number of particles
        n_x = l_box / ((grid_f + grid_ps)/2)                    # number of grid segments along x
        n_y = w_box / ((grid_f + grid_ps)/2)                    # number of grid segments along y
        n_z = h_box / ((grid_f + grid_ps)/2)                    # number of grid segments along z
        n_cells = n_x * n_y * n_z                               # number of cells
        n_sim = 7 * n_cells                                     # number of simulated particles
        f_num = n_real / n_sim                                  # f_num for SPARTA

        # Print the results
        print("Minimum grid size of x=%.3e, y=%.3e, z=%.3e" % (n_x, n_y, n_z))
        print("timestep of %.3e s, nrho=%.3e, f_num=%.3e" % (dt, nrho, f_num))
        print("Knudsen number is %.3e" % Kn)

        # Save the results to an input input
        if save_to_input:
            # Convert STL only for the first altitude
            if h == hs[0]:
                print("Converting binary STL to SPARTA surface...")
                # Convert STL from binary to asci
                os.system("python2 \"%s/tools/stl_B2A.py\" \"%s/setup/STL/%s.stl\"" % (sys.path[0], sys.path[0], s_name))
                # Convert STL to data surface for SPARTA
                os.system("python2 \"%s/tools/stl2surf.py\" \"%s/setup/STL/%s_ASCII.stl\" \"%s/setup/data/data.%s\"" % (sys.path[0], sys.path[0], s_name, sys.path[0], s_name))
            print("Saving input to file...")
            # Setup the SPARTA inputs
            input_s = "# SPARTA input file for satellite %s, for an altitude of %.1fkm\n" % (s_name, h)
            input_s += "seed               12345\n"
            input_s += "dimension          3\n"
            input_s += "\n"
            input_s += "global              gridcut 1e-6 comm/sort yes surfmax 10000 splitmax 100\n"
            input_s += "\n"
            input_s += "boundary           o r r\n"
            input_s += "create_box         -%.4f %.4f -%.4f %.4f -%.4f %.4f\n" % (l_box/2, l_box/2, w_box/2, w_box/2, h_box/2, h_box/2)
            input_s += "\n"
            # prevent grid too big by using min()
            input_s += "create_grid         %.4e %.4e %.4e\n" % (min(n_x, 15), min(n_y, 15), min(n_z, 15))
            input_s += "\n"
            input_s += "balance_grid        rcb cell\n"
            input_s += "\n"
            f = 1e3 if h == 115 else 1e7 if h == 150 else 1 # increase number of simulated particles to avoid having 0
            input_s += "global              nrho %.4e fnum %.4e\n" % (nrho, f_num/f)
            input_s += "\n"
            input_s += "species             ../atmo.species CO2 N2 Ar CO O O2\n"
            for n, sp_n in enumerate(species_names):
                input_s += "mixture             atmo %s vstream -%.4f 0.0 0.0 frac %.4f\n" % (sp_n, u_s, species_frac[n])
            input_s += "\n"
            input_s += "read_surf           ../data/data.%s\n" % (s_name)
            input_s += "surf_collide        1 diffuse 293.15 0.2\n"
            input_s += "surf_modify         all collide 1\n"
            input_s += "\n"
            input_s += "collide             vss atmo ../atmo.vss\n"
            input_s += "\n"
            input_s += "fix                 in emit/face atmo xhi twopass\n"
            input_s += "\n"
            input_s += "timestep            %.4e\n" % (dt)
            input_s += "\n"
            input_s += "compute             2 surf all all press fx fy fz px py pz etot\n"
            input_s += "fix                 save ave/surf all 1 5 5 c_2[*] ave running\n"
            input_s += "dump                1 surf all 5 ../tmp_result/%s_%skm.*.dat f_save[*]\n" % (s_name, h)
            input_s += "\n"
            input_s += "stats               100\n"
            input_s += "stats_style         step cpu np nscoll nscheck nexit\n"
            input_s += "run                 500\n"
            
            # Write SPARTA inputs to input
            with open(sys.path[0] + "/setup/inputs/in.%s_%skm" % (s_name, h), "w") as input_f:
                input_f.write(input_s)

