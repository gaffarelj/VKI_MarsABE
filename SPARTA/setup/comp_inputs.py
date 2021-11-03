import numpy as np
import sys
sys.path.insert(0,"/".join(sys.path[0].split("/")[:-2]))
import os
import shutil
from tools import plot_utilities as PU

tot_epochs = [1000, 1000, 1000]    # Number of simulation epochs for each altitude (should be multiple of 5)

# Define conditions at different orbital altitudes
hs = [85]#, 115, 150]
rhos = [2.69329E-06, 1.37030E-06, 3.04190E-10]
ps = [7.36261E-02, 3.73367E-02, 1.06037E-05]
Ts = [143.093, 128.038, 166.587]
Vs = [3510.90, 3495.84, 3475.51]
fracs = [
    np.array([95.874, 1.781, 1.822, 0.249, 0.126, 0.148])/100,
    np.array([90.556, 3.295, 3.292, 1.605, 0.971, 0.281])/100,
    np.array([69.792, 10.613, 5.061, 6.932, 6.655, 0.947])/100
]

run_all_cmd = "#!/bin/sh\nmodule load openmpi\n"
paraview_surf = ""
paraview_grid = ""
sat_names = ["CS_0021"]#["CS_0020", "CS_0021", "CS_1020", "CS_1021", "CS_2020", "CS_2021", "CS_2120", "CS_3020", "CS_3021"]
L_sats = [0.6 if _n[-1] == "1" else 0.3 for _n in sat_names]
L_s = [0.3, 0.589778, 0.341421, 0.589778, 0.541421, 0.589778, 0.6, 0.741421, 0.741421]
for j, s_name in enumerate(sat_names):
    print("\n\n* Satellite", s_name)
    # Create folder for results of this satellite
    try:
        os.mkdir(sys.path[0]+"/SPARTA/setup/results_sparta/"+s_name+"/")
    except (FileExistsError, OSError):
        try:
            # Un/comment the two following lines to always remove the previous results when new input files are made
            shutil.rmtree(sys.path[0]+"/SPARTA/setup/results_sparta/"+s_name+"/")
            os.mkdir(sys.path[0]+"/SPARTA/setup/results_sparta/"+s_name+"/")
        except (PermissionError, OSError):
            print("Warning: could not delete folder", s_name)
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
        h_box = 0.7     # box height [m]
        w_box = 0.7     # box width [m]
        l_box = 1.10    # box length [m]
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
        weighted_m = sum(species_frac * species_m)                      # weighted mass of species
        weighted_sigma = sum(species_frac * species_sigma)              # weighted frontal area of species
        nrho = rho / weighted_m                                         # number density [#/m3]
        lambda_f = 1 / (np.sqrt(2) * weighted_sigma * nrho)             # mean free path [m]
        Kn = lambda_f / L                                               # Knudsen number [-]
        T_ps = weighted_m * u_s**2 / (7*k_b)                            # post-shock T (Cv=2.5*R) [K]
        cr_ps = np.sqrt(16*k_b*T_ps / (np.pi*weighted_m))               # post-shock average relative velocity [m/s]\
        nrho_ps = 23/3*nrho                                             # post-shock number density [#/m3]
        lambda_ps = 1 / (np.sqrt(2) * weighted_sigma * nrho_ps)         # post-shock mean free path [m]
        nu_ps = weighted_sigma * nrho_ps * cr_ps                        # post-shock collision frequency [Hz]
        tau_ps = 1 / nu_ps                                              # post-shock collision time [s]
        dt = tau_ps / 5                                                 # time step [s]
        grid_f_mfp = lambda_f / 5                                       # grid dimension before shock [m] (based on mean free path)
        grid_ps_mfp = lambda_ps / 5                                     # post-shock grid dimension [m] (based on mean free path)
        grid_f_vel = u_s*dt                                             # grid dimension before shock [m] (based on velocity)
        grid_ps_vel = cr_ps*dt                                          # post-shock grid dimension [m] (based on velocity)
        grid_f = min(grid_f_mfp, grid_f_vel, L/25)                      # Take minimum grid dimension (or L_ref/25, to avoid grid too small)
        grid_ps = min(grid_ps_mfp, grid_ps_vel, L/25)                   # Take minimum grid dimension (or L_ref/25, to avoid grid too small)
        n_real = (nrho + nrho_ps) / 2 * h_box * l_box * w_box           # real number of particles
        n_x = int(l_box / ((grid_f + grid_ps)/2))                       # number of grid segments along x
        n_y = int(w_box / ((grid_f + grid_ps)/2))                       # number of grid segments along y
        n_z = int(h_box / ((grid_f + grid_ps)/2))                       # number of grid segments along z
        n_cells = n_x * n_y * n_z                                       # number of cells
        n_sim = 25 * n_cells                                            # number of simulated particles (int factor results from analysis to have 10 ppc)
        f_num = n_real / n_sim                                          # f_num for SPARTA
        # Check that dt is small enough given dx and v
        dx = min(l_box/n_x, w_box/n_y, h_box/n_z)
        dt = min(dt, dx/u_s*0.75, dx/cr_ps*0.75)                        # Take smallest dt of all (factor of 0.75 to make sure to be below the limit imposed by velocity)
        # Make sure that the initial grid size (before refinement) is not too small
        n_x, n_y, n_z = min(n_x, 50*l_box), min(n_y, 50*w_box), min(n_z, 50*h_box)
        # Compute the accomodation coefficient based on the adsorption of atomic oxygen
        # https://doi.org/10.2514/1.49330
        K = 7.5E-17                     # model fitting parameter
        n_0 = nrho * species_frac[-2]   # number density of the atomic oxygen
        P = n_0 * T                     # atomic oxygen partial pressure
        alpha = K*P/(1+K*P)             # accomodation coefficient
        test_accomodation = False
        if test_accomodation and s_name == sat_names[0] and h == hs[0]:
            P_s = np.linspace(1.5e17, 9e18, 200)
            alpha_s = [K*_P/(1+K*_P) for _P in P_s]
            PU.plot_single(P_s, alpha_s, "$n_O \cdot T [k m^3]$", "accomodation $\\alpha$", "test_accomodation")

        # Print the results
        print("Minimum grid size of x=%.3e, y=%.3e, z=%.3e" % (n_x, n_y, n_z))
        print("timestep of %.3e s, nrho=%.3e, f_num=%.3e" % (dt, nrho, f_num))
        print("Knudsen number is %.3e" % Kn)

        ## Save the results to an input input
        # Convert STL only for the first altitude
        if h == hs[0]:
            # Write command to convert surface to ParaView
            paraview_surf += "pvpython ../../tools/surf2paraview.py ../../setup/data/data.%s %s \n" % (s_name, s_name)
            print("Converting binary STL to SPARTA surface...")
            # Convert STL from binary to asci
            os.system("python2 \"%s/SPARTA/tools/stl_B2A.py\" \"%s/SPARTA/setup/STL/%s.stl\" -rs" % (sys.path[0], sys.path[0], s_name))
            # Convert STL to data surface for SPARTA
            os.system("python2 \"%s/SPARTA/tools/stl2surf.py\" \"%s/SPARTA/setup/STL/%s_ASCII.stl\" \"%s/SPARTA/setup/data/data.%s\"" % (sys.path[0], sys.path[0], s_name, sys.path[0], s_name))
        print("Saving input to file...")

        # Setup the SPARTA inputs
        input_s =  "# SPARTA input file for satellite %s, for an altitude of %.1fkm\n" % (s_name, h)
        input_s += "print \"\"\nprint \"***** Running SPARTA simulation for %s, at h=%ikm *****\"\nprint \"\"\n" % (s_name, h)
        input_s += "seed                12345\n"
        input_s += "dimension           3\n"
        grid_def = "dimension           3\n"
        input_s += "\n"
        input_s += "global              gridcut 1e-3 comm/sort yes surfmax 10000 splitmax 100\n"
        input_s += "\n"
        input_s += "boundary            o r r\n"
        input_s += "create_box          -%.4f %.4f -%.4f %.4f -%.4f %.4f\n" % (l_box/2, l_box/2, w_box/2, w_box/2, h_box/2, h_box/2)
        grid_def+= "create_box          -%.4f %.4f -%.4f %.4f -%.4f %.4f\n" % (l_box/2, l_box/2, w_box/2, w_box/2, h_box/2, h_box/2)
        input_s += "\n"
        input_s += "create_grid         %i %i %i\n" % (np.ceil(n_x), np.ceil(n_y), np.ceil(n_z))
        input_s += "\n"
        input_s += "balance_grid        rcb cell\n"
        input_s += "\n"
        input_s += "global              nrho %.4e fnum %.4e\n" % (nrho, f_num)
        input_s += "\n"
        input_s += "species             ../atmo.species CO2 N2 Ar CO O O2\n"
        for n, sp_n in enumerate(species_names):
            input_s += "mixture             atmo %s vstream -%.4f 0.0 0.0 frac %.4f\n" % (sp_n, u_s, species_frac[n])
        input_s += "collide             vss atmo ../atmo.vss\n"
        input_s += "\n"
        input_s += "read_surf           ../data/data.%s trans %.4f 0 0\n" % (s_name, (0.075+L_sats[j]/2))
        input_s += "surf_collide        1 diffuse 293.15 %.4f\n" % (alpha)
        input_s += "surf_modify         all collide 1\n"
        input_s += "\n"
        input_s += "fix                 in emit/face atmo xhi twopass\n"
        input_s += "\n"
        input_s += "timestep            %.4e\n" % dt
        input_s += "\n"
        input_s += "compute             forces surf all all fx fy fz\n"
        input_s += "fix                 avg ave/surf all %i %i %i c_forces[*] ave running\n" % (tot_epochs[i]/1000, tot_epochs[i]/250, tot_epochs[i]/50)
        input_s += "compute             sum_force reduce sum f_avg[*]\n"
        input_s += "\n"
        input_s += "compute             grid_data grid all all n nrho massrho u temp\n"
        input_s += "fix                 grid_avg ave/grid all %i %i %i c_grid_data[*]\n" % (tot_epochs[i]/1000, tot_epochs[i]/250, tot_epochs[i]/50)
        input_s += "compute             avg_ppc reduce ave f_grid_avg[1]\n"
        input_s += "\n"
        input_s += "compute             knudsen lambda/grid f_grid_avg[2] f_grid_avg[5] CO2\n"
        input_s += "\n"
        input_s += "stats               %i\n" % (tot_epochs[i]/50)
        input_s += "stats_style         step cpu np nscoll nexit c_sum_force[*] c_avg_ppc\n"
        input_s += "\n"
        input_s += "dump                0v grid all %i ../results_sparta/%s/gridvals_%ikm_0.*.dat id f_grid_avg[*]\n" % (tot_epochs[i]/10, s_name, h)
        input_s += "dump                0K grid all %i ../results_sparta/%s/gridKn_%ikm_0.*.dat id c_knudsen[*]\n" % (tot_epochs[i]/10, s_name, h)
        input_s += "write_grid          ../results_sparta/%s/grid_%ikm_0.dat\n" % (s_name, h)
        grid_def = [grid_def]*3
        grid_def[0] += "read_grid           ../../setup/results_sparta/%s/grid_%skm_0.dat\n" % (s_name, h)
        input_s += "run                 %i\n" % (tot_epochs[i] * 1/2)
        input_s += "\n"
        input_s += "adapt_grid          all refine coarsen value c_knudsen 0.05 0.1 combine min thresh less more\n"
        input_s += "dump                1v grid all %i ../results_sparta/%s/gridvals_%ikm_1.*.dat id f_grid_avg[*]\n" % (tot_epochs[i]/10, s_name, h)
        input_s += "dump                1K grid all %i ../results_sparta/%s/gridKn_%ikm_1.*.dat id c_knudsen[*]\n" % (tot_epochs[i]/10, s_name, h)
        input_s += "write_grid          ../results_sparta/%s/grid_%ikm_1.dat\n" % (s_name, h)
        grid_def[1] += "read_grid           ../../setup/results_sparta/%s/grid_%skm_1.dat\n" % (s_name, h)
        input_s += "run                 %i\n" % (tot_epochs[i] * 3/10)
        input_s += "\n"
        input_s += "adapt_grid          all refine coarsen value c_knudsen 0.05 0.1 combine min thresh less more\n"
        input_s += "dump                2v grid all %i ../results_sparta/%s/gridvals_%ikm_2.*.dat id f_grid_avg[*]\n" % (tot_epochs[i]/10, s_name, h)
        input_s += "dump                2K grid all %i ../results_sparta/%s/gridKn_%ikm_2.*.dat id c_knudsen[*]\n" % (tot_epochs[i]/10, s_name, h)
        input_s += "write_grid          ../results_sparta/%s/grid_%ikm_2.dat\n" % (s_name, h)
        grid_def[2] += "read_grid           ../../setup/results_sparta/%s/grid_%skm2.dat\n" % (s_name, h)
        input_s += "run                 %i\n" % (tot_epochs[i] * 1/5)
        
        run_all_cmd += "mpirun -np 16 spa_ < in.%s_%skm | tee ../results_sparta/%s/stats_%ikm.dat\n" % (s_name, h, s_name, h)
        for i_r in range(3):
            paraview_grid += "\n"
            paraview_grid += "rm -rf %s_%skm_%i \n" % (s_name, h, i_r)
            paraview_grid += "rm -rf %s_%skm_%i.pvd \n" % (s_name, h, i_r)
            paraview_grid += "echo 'Converting result grid of %s at %skm (refinement %i) to ParaView...'\n" % (s_name, h, i_r)
            paraview_grid += "pvpython ../../tools/grid2paraview_original.py def/grid.%s_%skm_%i vals_%s_%skm_%i -r ../../setup/results_sparta/%s/gridvals_%skm_%i.*.dat \n" % \
                (s_name, h, i_r, s_name, h, i_r, s_name, h, i_r)
            paraview_grid += "pvpython ../../tools/grid2paraview_original.py def/grid.%s_%skm_%i Kn_%s_%skm_%i -r ../../setup/results_sparta/%s/gridKn_%skm_%i.*.dat \n" % \
                (s_name, h, i_r, s_name, h, i_r, s_name, h, i_r)
        
        # Write SPARTA inputs to input
        with open(sys.path[0] + "/SPARTA/setup/inputs/in.%s_%skm" % (s_name, h), "w") as input_f:
            input_f.write(input_s)

        # Write grid definition in ParaView folder
        for i_g in range(3):
            with open(sys.path[0] + "/SPARTA/paraview/grid/def/grid.%s_%skm_%i" % (s_name, h, i_g), "w") as input_f:
                input_f.write(grid_def[i])


# Write command to run all SPARTA input files
try:
    with open(sys.path[0] + "/SPARTA/setup/inputs/run_all.sh", "r+") as run_f:
        run_f.seek(0)
        run_f.write(run_all_cmd)
        run_f.truncate()
except FileNotFoundError:
    with open(sys.path[0] + "/SPARTA/setup/inputs/run_all.sh", "w") as run_f:
        run_f.write(run_all_cmd)
# Write command to create ParaView files
paraview_cmd = "#!/bin/sh\n"
paraview_cmd += "cd surf\n"
paraview_cmd += "rm -rf *\n"
paraview_cmd += paraview_surf
paraview_cmd += "cd ../grid\n"
paraview_cmd += paraview_grid
try:
    with open(sys.path[0] + "/SPARTA/paraview/paraview_convert.sh", "r+") as run_f:
        run_f.seek(0)
        run_f.write(paraview_cmd)
        run_f.truncate()
except FileNotFoundError:
    with open(sys.path[0] + "/SPARTA/paraview/paraview_convert.sh", "w") as run_f:
        run_f.write(paraview_cmd)