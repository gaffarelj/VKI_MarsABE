import numpy as np
import pygmo
import time
import sys
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
from optimisation.HT_problem import HT_problem
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
fitness_weights = [5, 10, 1]
current_HT_problem = HT_problem(design_var_range, fitness_weights, verbose=False)
problem = pygmo.problem(current_HT_problem)

# Setup Pygmo
seed = 12345
pop = pygmo.population(problem, size=12, seed=seed)
algo = pygmo.algorithm(pygmo.nsga2(seed=seed, cr=0.95, eta_c=10, m=0.001, eta_m=2))

# Prepare variables to save fitnesses and populations
all_f, all_p = [list(f) for f in pop.get_f()], [list(p) for p in pop.get_x()]

# Run the optimisation
n_generations = 5
t0 = time.time()
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
    
    # Show time it took
    dt, t0 = time.time() - t0, time.time()
    print(" -> took %.1f seconds" % dt)

all_f, all_p = np.array(all_f), np.array(all_p)

# Find the optimum if the sum of fitnesses is used
idx_best = np.where(np.sum(all_f, axis=1) == np.min(np.sum(all_f, axis=1)))[0][0]
optimum_f, optimum_p = all_f[idx_best], all_p[idx_best]

print("Optimum: %s \n  with initial state: h_p=%3d km, h_a=%3d km, i=%2d, omega=%3d, Omega=%.3d" % \
    (SM.satellites[list(SM.satellites.keys())[int(optimum_p[5])]], min(optimum_p[0:1])/1e3, max(optimum_p[0:1])/1e3, \
        np.rad2deg(optimum_p[2]), np.rad2deg(optimum_p[3]), np.rad2deg(optimum_p[4])))
print("Resulting optimum fitness:", optimum_f)
print("Explored %i different possibilities." % len(all_f))