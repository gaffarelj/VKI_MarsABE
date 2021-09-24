import numpy as np

# Define conditions at different orbital altitudes
hs = [85, 115, 150]
rhos = [7.1E-07, 1.8E-08, 1.6E-10]
ps = [2.3E-02, 3.7E-04, 7.1E-06]
Ts = [135, 115, 175]
Vs = [3494.17, 3493.29, 3483.82]
fracs = [
    np.array([0.905, 0.035, 0.025, 0.015, 0.015, 0.005]),
    np.array([0.81, 0.045, 0.035, 0.05, 0.055, 0.005]),
    np.array([0.42, 0.125, 0.045, 0.15, 0.25, 0.01])
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