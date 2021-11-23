# Optimisation

This folder contains the files used to run the optimisation aiming at finding the best combination of satellite configuration and orbit to compensate for drag in very low Martian orbit.

## Optimisation setup

The optimisation setup is mainly composed of the optimisation objectives, the design variables, and the constraints.

### Objectives

The optimisation that has been setup is a multiple-objectives optimisation. These objectives are the followings:

 * Periapsis decay: this should be as low as possible, and is measured by taking the difference between the periapsis at the beginning and at the end of the orbital simulation.
 * Mean altitude: this should also be as low as possible since we want a very low orbit, and it is measured by average all the geopotential altitudes during the orbital simulation.
 * Mean power: this should be maximised, because the power available to the satellite should be as high as possible (to power both the satellite subsystems and the thruster).
 * Mean Thrust/Drag ratio: this should also be maximised (or at least be above 1), because the thrust should be able to compensate the drag.

The scaling function used to compute the fitness from the objectives are described below. These functions are used to make sure that the fitness is closed to values between 0 and 1, and that a low fitness means a better objective.

The *periapsis decay*, computed in meters, is mapped to its fitness value using the following function:

<img src="https://latex.codecogs.com/png.latex?\bg_white%20f_%5Ctext%7Bperiapsis%20decay%7D%20%3D%20%5Cfrac%7B%5CDelta%20h_p%7D%7B10%5E5%7D">

This function thus maps decay from 0 to 100km to 0 to 1. A decay above 100km will thus have a fitness above 1, and a negative decay (a raise in periapsis) will have a negative fitness. These behaviors are accepted as the optimisation algorithm will make sense of these fitnesses.

Mapping a periapsis decay from -150km to 150km then results in the fitnesses of the plot that can be found [here](../figures/optimisation/decay_scale.pdf).

The *mean altitude*, also computed in meters, is mapped differently whether it is below or above 500km.
The following function is used to do so:

<img src="https://latex.codecogs.com/png.latex?\bg_white%20%5Cleft%5C%7B%5Cbegin%7Bmatrix%7D%20f_%7B%5Cmu_h%7D%20%3D%200.75%20%5Ccdot%20%5Cfrac%7B%5Cmu_h%20-%2050%20%5Ccdot%2010%5E3%7D%7B%28500-50%29%5Ccdot%2010%5E3%7D%20%26%20%5Ctext%7Bif%7D%5C%20%5Cmu_h%20%5Cleq%2050%20%5Ccdot%2010%5E3%20%5C%5C%20f_%7B%5Cmu_h%7D%20%3D%200.75%20%2B%200.25%20%5Ccdot%20%5Cfrac%7B%5Cmu_h%20-%20500%20%5Ccdot%2010%5E3%7D%7B%281000-500%29%5Ccdot%2010%5E3%7D%20%26%20%5Ctext%7Botherwise%7D%20%5Cend%7Bmatrix%7D%5Cright.">

This mapping function results in the fitnesses from the plot [here](../figures/optimisation/altitude_scale.pdf) for mean altitudes between 50km and 650km.

A *mean power* between 0 and 35W is scaled to a fitness between 1 and 0 by using the following function:

<img src="https://latex.codecogs.com/png.latex?\bg_white%20f_%7B%5Cmu_P%7D%20%3D%20%5Cfrac%7B35%20-%20%5Cmu_P%7D%7B35%7D">

With this, a high mean power gets a lower (better) fitness. The value of 35W is used because no satellite configuration at no orbit could reach over 35W.
This scaling results in the mean power vs power fitness plot that is [here](../figures/optimisation/power_scale.pdf).

Finally, the unitless *mean Thrust/Drag*, that can range from 0 to high values, is mapped between 0 to 1, using the following inverse function:

<img src="https://latex.codecogs.com/png.latex?\bg_white%20f_%7B%5Cmu_%7BT%2FD%7D%7D%20%3D%20%5Cfrac%7B1%7D%7B%5Cmu_%7BT%2FD%7D%20%2B%201%7D">

This results in the fitnesses mapping illustrated [here](../figures/optimisation/TD_scale.pdf), for Thrust/Drag ratios between 0 and 50.

### Design variables

The design variables that the optimiser can play around with to change the problem parameters are the following:
 * Periapsis initial geopotential altitude (in meters): this value ranges between 85km and 150km.
 * Apoapsis initial geopotential altitude (in meters): this value ranged between 85km and 500km.
 * Initial orbital inclination (in radians): this value ranges between 0 (equatorial orbit) and pi/2 (polar orbit).
 * Initial argument of periapsis of the orbit (in radians): this value ranges between 0 and pi.
 * Initial longitude of the ascending node of the orbit (in radians): this value ranges between 0 and 2pi.
 * Index of the satellite to use (integer): this value can take any integer value between 0 and the number of satellites defined in [this module](../utils/sat_models.py), minus 1.

Note that, because the range for the periapsis and apoapsis initial altitudes overlaps, the second one could end up being the periapsis and the first one the apoapsis.

### Constraints

The satellite properties such as its mass, volume, or solar panel size are hidden because pre-determined satellite configurations are used instead.

Also, while constraints on the maximum thermo-mechanical loads was first thought, it has been decided to keep this free, and analyse the results afterwards to check if the required thrust and the drag do not over-stress the structure of the satellite.

Furthermore, constrains have been used on the initial orbital elements such as the inclination, argument of periapsis, and longitude of the ascending node, so that the values taken fit the limits imposed by definition.

Finally, constrains have been set on the altitude. The minimum periapsis and apoapsis altitude is set to 85km, which resorted from the feasible altitudes study [here](../MCD/feasible_altitudes.py).
The maximum periapsis altitude is set to 150km, also as found from the feasible altitudes study. The maximum apoapsis is at 500km, allowing for a still relatively low apoapsis, while letting the satellite get out of the atmosphere for a good part of its orbit.

## Files

Multiple files are present in this folder, relating to the optimisation problem itself, running it, and analysing the results.

### Drag compensation problem

The file [drag_comp_problem.py](drag_comp_problem.py) defines the problem of drag compensation in very low Mars orbit.

This file first defines a function, `comp_fitness()`, that takes the design variables and the thrust model as inputs. A thrust model equal to `2` means that a fuel tank is used, and the mass propagated. A thrust model equal to `3` uses the same thruster but with the atmosphere-breathing inlet.
Then, the orbital simulation is run using [propagation.py](../utils/propagation.py), and the scaled fitness values are computed based on the results.

Then, the `DC_problem()` class contains the problem definition, understood by Pygmo.

In its initialisation (`DC_problem.__init__()`), the design variables are inputs, as well as the thrust model, and whether to print information during the runs or not (the verbose).

The `DC_problem.get_bounds()` function returns the bounds in which each design variable can vary.
`DC_problem.get_nobj()` returns the number of objectives for the optimisation. `DC_problem.get_nix()` returns the integer dimension of the problem. Set to 1, it means that one of the design variables (the index of the satellite configuration), needs to be an integer.

Finally, `DC_problem.fitness()` returns the fitness values for given design variables. It mostly calls the `comp_fitness()` that is defined at the top of the file, and that runs the orbital simulation and computes the fitnesses.

### Run problem

Running the drag compensation problem described above can be done using the [run_problem.py](run_problem.py) file.

In this file, the option of which thrust model to use is first offered. Then, the range of the design variables is setup, and the `DC_problem` defined.

Follows the selection of the optimiser, selected to be best as the Non-dominated Sorting Genetic Algorithm, with a population of 60 (the number of design variables by a factor 10).

The optimisation is then either run for up to 100 generations, and the results saved between each generation, or the latest saved results can be analysed.

During the results analysis, a plot is first made of the fitness history over the generation number.
Then, various Pareto fronts are plotted, comparing two objectives each time (with sometime an added objective of design variable as a colormap).

Finally, an interactive Pareto front is plotted, using the mean altitude and periapsis decay objectives, and saving the plot as an html file where each of the dot can be hovered to see both the design variables resulting in it, and the resulting objective scores.

### Pursue case

The [pursue_case.py](pursue_case.py) file can be used to simulate a single case. This can be used for instance to simulate one of the Pareto-optimum solutions found during the optimisation, for a longer number of days, using the MCD for the climate and winds, and more detailed environmental accelerations.

This file also allows for more detailed analysis of some design variables or of the states of the propagated satellite.

Finally, this file saves the positions of the propagated satellite over time in a `csv` file.

## Optimisation results

The results from the optimisation problem can be analysed in various ways.

### Data file

First, as stated above, the population is saved between each generation during the optimisation.

More concretely, a numpy `npz` file is created containing all of the distinct population member ever created, with their associated design variables and objective scores. This file can be accessed in the [results](results) folder, and the file is named as `DC_1-2-3_4-5_6-7`:
 1. Thrust model (typically `2` or `3`).
 2. Indication of the use of battery (`X` = no battery, `V` = battery).
 3. Ionisation efficiency of the atmosphere (float between 0 and 1, only used if the thrust model is `3`).
 4. Population size (typically `60`).
 5. Seed used for the optimisation (usually `12345`).
 6. Time at which the optimisation starts (in the following format: `DMYHMS`).
 7. Generation number that was saved (this is the last generation that was fully run before the process was stopped).

This is done between each generation rather than at the end so that the precious time taken to generated this increasingly better population is not wasted in case the script is stopped before completion.
To save space, avoid data redundancy, and avoid confusion, the data file from the previous generation is erased after the data from the newest generation is saved.

### Pareto fronts

All of the generated Pareto fronts are saved in the [figures/optimisation](../figures/optimisation) folder. The `ABE` subfolder is used to save the plots related to the optimisation made with thrust model `3`, using the atmosphere as propellant, and the `with_tank` subfolder is used when the thrust model `2` is used, with a Xenon propellant tank.

### ParaView

The orbital state of the satellite saved as a `csv` file can also be accessed in the [results](results) folder. 
These `csv` files contain the same `ABE` or `with_tank` naming convention as specified in the sub-section above.

The [mars_visu.pvsm](results/paraview/mars_visu.pvsm) can be loaded in ParaView to visualise the orbit taken by the satellite around Mars.