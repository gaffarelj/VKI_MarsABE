import numpy as np

# Define conditions at different orbital altitudes
hs = [95, 140, 190]
rhos = [4.9E-07, 1.2E-09, 6.4E-12]
ps = [8.5E-03, 1.5E-05, 3.5E-07]
Ts = [125, 150, 155]
Vs = [3303.14, 3312.32, 3402.66]
fracs = [
    np.array([0.8, 0.095, 0.045, 0.0275, 0.0275, 0.005]),
    np.array([0.67, 0.07, 0.05, 0.075, 0.125, 0.01]),
    np.array([0.26, 0.13, 0.02, 0.15, 0.43, 0.01])
]

# Loop trough conditions
for i, h in enumerate(hs):
    print("\nWith conditions at altitude of %i km:" % h)
    print("Velocity is of %.2f m/s, and the species are mixed as follows:" % Vs[i])
    print(fracs[i])
    # Inputs
    rho = rhos[i]   # density [kg/m3]
    p = ps[i]       # pressure [Pa]
    T = Ts[i]       # temperature [K]
    L = 0.37        # reference length [m] (satellite width)
    u_s = Vs[i]     # free-stream velocity [m/s]
    h_box = 0.5     # box height [m]
    w_box = 0.5     # box width [m]
    l_box = 0.7     # box length [m]
    # Fraction of each species
    species_frac = fracs[i]
    if sum(species_frac) != 1:
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
    print("(Use f_num=%.3e for better accuracy)" % (f_num*50))
    print("Knudsen number is %.3e" % Kn)