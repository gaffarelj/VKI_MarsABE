import numpy as np
import sys
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
from utils import sat_models as SM
from utils import propagation as P
from tudatpy.kernel import constants


class HT_problem:

    def __init__(self, design_var_range, fitness_weights):
        self.design_var_range = design_var_range
        self.fitness_weights = fitness_weights

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
        h_0, test = design_variables

        ## Setup the simulation
        # Select the satellite
        satellite = SM.satellites["CS_0020"]
        # Create the orbital simulation instance, setup to simulate 5 days
        OS = P.orbit_simulation(satellite, "Mars", 5*constants.JULIAN_DAY, save_power=True)
        # Create the simulation bodies, and use the MCD
        OS.create_bodies(use_MCD=[False, False])
        # Create the initial state of the satellite
        OS.create_initial_state(h=h_0)
        # Load the accelerations from default config 1: Central body spherical harmonics of degree/order 4 and aerodynamics, Solar radiation
        OS.create_accelerations(default_config=1, thrust=0)
        # Create the integrator, termination settings, dependent variables, and propagator
        OS.create_integrator()
        OS.create_termination_settings()
        OS.create_dependent_variables(to_save=["h"])
        OS.create_propagator(prop_mass=False)
        
        # Simulate the satellite in orbit
        times, states, dep_vars = OS.simulate()

        # Extract the results from the simulation
        power_hist = list(OS.power_dict.values())
        altitudes = dep_vars[:,0]
        decay = altitudes[0] - altitudes[-1]

        print("Start from h=%3d, sum of power=%6d W, total decay=%3d km" % (h_0/1e3, sum(power_hist), decay/1e3))

        ## Compute the fitness (=cost); scaling is used because, ideally, all cost values would be in the same range (0-1 for instance)
        # Total power; lots of power = smaller value = better (use maximum observed value of 25kW to scale)
        power_f = 25e3/sum(power_hist)
        # Total decay; if final altitude < 50km, fitness=1 (re-entered atmosphere); else: scale with maximum 100km
        if altitudes[-1] <= 50e3:
            decay_f = 1
        else:
            decay_f = decay/100e3

        cost = np.array([power_f, decay_f]) * np.array(self.fitness_weights)
        return cost


test = True
if test:
    import pygmo
    # Setup the design variables range
    min_h, max_h = 85e3, 150e3
    test_a, test_b = 0, 1
    design_var_range = (
        [min_h, test_a],
        [max_h, test_b]
    )

    # Setup the optimisation problem
    fitness_weights = [1, 0]
    current_HT_problem = HT_problem(design_var_range, fitness_weights)
    problem = pygmo.problem(current_HT_problem)

    # Setup Pygmo
    seed = 12345
    pop = pygmo.population(problem, size=8, seed=seed)
    algo = pygmo.algorithm(pygmo.nsga2(seed=seed, cr=0.95, eta_c=10, m=0.001, eta_m=2))

    # Prepare variables to save optimum
    optimum_f, optimum_p, optimum_sum_f = None, None, np.inf

    # Run the optimisation
    n_generations = 1
    for i in range(1,n_generations+1):
        print("Running generation %2d/%2d" % (i, n_generations))
        # Evolve the population
        pop = algo.evolve(pop)

        # Get the new fitnesses and population
        f = pop.get_f()
        p = pop.get_x()
        
        # If a new optimum is reached, save it
        opti_sum = np.min(np.sum(f, axis=0))
        if opti_sum < optimum_sum_f:
            idx_best = np.where(np.sum(f, axis=0) == opti_sum)[0][0]
            optimum_f, optimum_p = f[idx_best], p[idx_best]

    
    print("Optimum design variables: h_0=%3d km" % (optimum_p[0]/1e3))
    print("Resulting optimum fitness:", optimum_f)