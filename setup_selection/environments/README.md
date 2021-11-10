# Environments selection

This folder contains scripts written to analyse and select an appropriate environment for the satellite orbiting Mars.

## Atmosphere analysis

The [atmo_summary.py](atmo_summary.py) file has been used to plot a vertical density profile of the atmosphere of Mars, including the contribution from each of the major molecular species.

This script generates a second script showing the same vertical profile of the molecular species as volumetric fractions.

## Mars Climate Database inclusion

The [test_MCD_atmosphere.py](test_MCD_atmosphere.py) script is a simulation test including the atmospheric densities and winds from the MCD.

## Required environments

The most important file in this folder is [test_req_envs.py](test_req_envs.py).
In this file, different environmental models are compared to each other, to see which ones are relevant to the simulation.

### Gravitational model
First, the gravitational model of Mars has been analysed.

It has been found that using a Spherical Harmonics (SH) model of Degree and Order (D/O) up to 4 leads to the most accurate results.

A deviation of more than 15% has been observed between Mars as a Point Mass (PM) and the SH up to D/O 2.
A slight deviation is also seen when using D/O 4.
However, it appears that including more D/O does not lead to significant deviations, and only increase the simulation time.

### Relativistic correction

Various accelerations can be added to account for the difference between Newtonian dynamics, and General Relativity.

The most important one, the Schwarzschild term, accounts for the gravity well that a mass produces.

However, including it in the simulation lead to no difference.

### Radiation pressure

The radiation pressure from the Sun has also been tested, using a cannonball model for the satellite.
No significant deviations were found.

### Other bodies contribution

The acceleration from the Sun and Jupiter as Point Masses on the satellite have been investigated.

They did not lead to significant deviations for a 20 days simulations. For simulations in the order of years, it is still recommended to include them, due to their small effect.

## Recommendations

In light of these findings, the following environment is recommended:
 * Mars gravitational acceleration from SH up to D/O 4.
 * Aerodynamic acceleration from the atmosphere.
 
If more precision is required, the following can be used instead:
 * Mars gravitational acceleration from SH up to D/O 8.
 * Aerodynamic acceleration from the atmosphere.
 * Radiation pressure from the Sun.
 * Point Mass gravity from the Sun and Jupiter.