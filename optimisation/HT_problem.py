import numpy as np
import sys
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
from utils import sat_models as SM
from utils import propagation as P
from tudatpy.kernel import constants
from tudatpy.kernel.interface import spice_interface


class HT_problem:

    def __init__(self, design_var_range, fitness_weights, verbose=False):
        self.design_var_range = design_var_range
        self.fitness_weights = fitness_weights
        self.verbose = verbose

    def get_bounds(self):
        """
        *** Pygmo-related function ***
        Return the bounds of the design variables
        """
        return self.design_var_range

    def get_nobj(self):
        """
        *** Pygmo-related function ***
        Return the number of objectives
        """
        return len(fitness_weights)

    def fitness(self, design_variables):
        """
        *** Pygmo-related function ***
        Return the fitness of the given problem. This is the cost function, that Pygmo will minimise.
        """
        # Extract the individual design variables
        h_a_0, h_p_0, i_0, omega_0, Omega_0 = design_variables

        ## Setup the simulation
        # Select the satellite
        satellite = SM.satellites["CS_0020"]
        # Create the orbital simulation instance, setup to simulate 5 days
        OS = P.orbit_simulation(satellite, "Mars", 5*constants.JULIAN_DAY, save_power=True)
        # Create the simulation bodies, and use the MCD
        OS.create_bodies(use_MCD=[False, False])
        # Create the initial state of the satellite
        R_Mars = spice_interface.get_average_radius("Mars")
        a_0 = R_Mars + (h_a_0+h_p_0)/2
        e_0 = 1 - (R_Mars + min(h_p_0, h_a_0)) / a_0    # Use min because h_p could actually be higher than h_a due to the way the problem is setup)
        OS.create_initial_state(a=a_0, e=e_0, i=i_0, omega=omega_0, Omega=Omega_0)
        # Load the accelerations from default config 1: Central body spherical harmonics of degree/order 4 and aerodynamics, Solar radiation
        OS.create_accelerations(default_config=1, thrust=1)
        # Create the integrator, termination settings, dependent variables, and propagator
        OS.create_integrator()
        OS.create_termination_settings()
        OS.create_dependent_variables(to_save=["h_p", "h"])
        OS.create_propagator(prop_mass=False)
        
        # Simulate the satellite in orbit
        times, states, dep_vars = OS.simulate()

        # Extract the results from the simulation
        power_hist = list(OS.power_dict.values())
        h_p_s = OS.get_dep_var("h_p")
        decay = h_p_s[0] - h_p_s[-1]
        altitudes = OS.get_dep_var("h")

        if self.verbose:
            print("Start from h_p=%3d, h_a=%.2f, i=%2d, omega=%3d, Omega=%.3d -> mean of power=%.2f W, total decay=%4d km, mean altitude=%3d km" % \
                (h_p_0/1e3, h_a_0/1e3, np.rad2deg(i_0), np.rad2deg(omega_0), np.rad2deg(Omega_0), np.mean(power_hist), decay/1e3, np.mean(altitudes)/1e3))

        ## Compute the fitness (=cost); scaling is used because, ideally, all cost values would be in the same range (0-1 for instance)
        # Max mean power; lots of power = smaller value = better (use maximum observed value as scale)
        power_scale = 10
        power_f = ( (power_scale-np.mean(power_hist)) / power_scale )
        if np.mean(power_hist) > power_scale:
            print("Warning, mean power of %.3f W was above scaling value of %.3f W" % (np.mean(power_hist), power_scale))
        # Min decay; if final altitude < 50km, fitness=1 (re-entered atmosphere); else: scale with maximum 100km
        decay_f = 1 if h_p_s[-1] <= 50e3 else decay/100e3
        # Min mean altitude
        h_scale = [50e3, 500e3, 1e6]
        if np.mean(altitudes) <= h_scale[1]:
            h_f = (np.mean(altitudes) - h_scale[0]) / (h_scale[1] - h_scale[0]) * 0.75
        else:
            h_f = 0.75 + 0.25 * (np.mean(altitudes) - h_scale[1]) / (h_scale[2] - h_scale[1])
        # Assemble and return the cost
        cost = np.array(self.fitness_weights) * np.array([power_f, decay_f, h_f]) # Mean power, periapsis decay, mean altitude
        if self.verbose:
            print("Cost =", cost)
        return cost


test = True
if test:
    import pygmo
    # Setup the design variables range
    min_h_p, max_h_p = 85e3, 150e3
    min_h_a, max_h_a = 85e3, 500e3
    min_e, max_e = 0, 0.1
    min_i, max_i = 0, np.pi/2
    min_omega, max_omega = 0, np.pi
    min_Omega, max_Omega = 0, 2*np.pi
    test_a, test_b = 0, 1
    design_var_range = (
        [min_h_p, min_h_a, min_i, min_omega, min_Omega],
        [max_h_p, max_h_a, max_i, max_omega, max_Omega]
    )

    # Setup the optimisation problem
    fitness_weights = [3, 3, 1]
    current_HT_problem = HT_problem(design_var_range, fitness_weights, verbose=False)
    problem = pygmo.problem(current_HT_problem)

    # Setup Pygmo
    seed = 12345
    pop = pygmo.population(problem, size=12, seed=seed)
    algo = pygmo.algorithm(pygmo.nsga2(seed=seed, cr=0.95, eta_c=10, m=0.001, eta_m=2))

    # Prepare variables to save fitnesses and populations
    all_f, all_p = [list(f) for f in pop.get_f()], [list(p) for p in pop.get_x()]

    # Run the optimisation
    n_generations = 10
    for i in range(1,n_generations+1):
        print("Running generation %2d/%2d" % (i, n_generations))
        # Evolve the population
        pop = algo.evolve(pop)

        # Get the new fitnesses and population
        f = pop.get_f()
        p = pop.get_x()

        # Save everything
        for i_f, _f in enumerate(f):
            if list(_f) not in all_f:
                all_p.append(list(p[i_f]))
                all_f.append(list(_f))

    all_f, all_p = np.array(all_f), np.array(all_p)

    # Find the optimum if the sum of fitnesses is used
    idx_best = np.where(np.sum(all_f, axis=1) == np.min(np.sum(all_f, axis=1)))[0][0]
    optimum_f, optimum_p = all_f[idx_best], all_p[idx_best]

    print("Optimum initial state: h_p=%3d km, h_a=%3d km, i=%2d, omega=%3d, Omega=%.3d" % \
        (min(optimum_p[0:1])/1e3, max(optimum_p[0:1])/1e3, np.rad2deg(optimum_p[2]), np.rad2deg(optimum_p[3]), np.rad2deg(optimum_p[4])))
    print("Resulting optimum fitness:", optimum_f)
    print("Explored %i different possibilities." % len(all_f))