import numpy as np
import pygmo
import time
import sys
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
from optimisation import HT_problem as HTp
from utils import sat_models as SM


# Setup the design variables range
min_h_p, max_h_p = 85e3, 150e3
min_h_a, max_h_a = 85e3, 500e3
min_e, max_e = 0, 0.1
min_i, max_i = 0, np.pi/2
min_omega, max_omega = 0, np.pi
min_Omega, max_Omega = 0, 2*np.pi
min_sat_i, max_sat_i = 0, len(SM.satellites) - 1
design_var_range = (
    [min_h_p, min_h_a, min_i, min_omega, min_Omega, min_sat_i],
    [max_h_p, max_h_a, max_i, max_omega, max_Omega, max_sat_i]
)

# Setup the optimisation problem
fitness_weights = [5, 10, 1]    # Mean power, periapsis decay, mean altitude
current_HT_problem = HTp.HT_problem(design_var_range, fitness_weights, verbose=False)
problem = pygmo.problem(current_HT_problem)

# Setup Pygmo
seed = 12345
pop = pygmo.population(problem, size=12, seed=seed)
algo = pygmo.algorithm(pygmo.nsga2(seed=seed, cr=0.95, eta_c=10, m=0.001, eta_m=2))

# Run the optimisation
n_generations = 5
t0 = time.time()
for i in range(1,n_generations+1):
    print("Running generation %2d/%2d" % (i, n_generations))
    # Evolve the population
    pop = algo.evolve(pop)
    
    # Show time it took
    dt, t0 = time.time() - t0, time.time()
    print(" -> took %.1f seconds" % dt)

fit_results, fit_inputs = np.array(HTp.FIT_RESULTS), np.array(HTp.FIT_INPUTS)

# Find the optimum if the sum of fitnesses is used
idx_best = np.where(np.sum(fit_results[:,:len(fitness_weights)], axis=1) == np.min(np.sum(fit_results[:,:len(fitness_weights)], axis=1)))[0][0]
optimum_input, optimum_result = fit_inputs[idx_best], fit_results[idx_best]

# Print the results
print(optimum_input)
print("Optimum: %s \n  with initial state: h_p=%3d km, h_a=%3d km, i=%2d, omega=%3d, Omega=%.3d" % \
    (optimum_input[0], min(optimum_input[1:3])/1e3, max(optimum_input[1:3])/1e3, \
        np.rad2deg(optimum_input[3]), np.rad2deg(optimum_input[4]), np.rad2deg(optimum_input[5])))
print("Resulting optimum fitness:", optimum_result[:3])
print("Resulting characteristics: mean power=%.2f W, total decay=%.2f km, mean altitude=%3d km, mean T/D=%.2f" % \
    (optimum_result[-4], optimum_result[-3]/1e3, optimum_result[-2]/1e3, optimum_result[-1]))
print("Explored %i different possibilities." % len(fit_inputs))