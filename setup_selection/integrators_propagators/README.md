# Integrators and propagators selection
This folders contains different various scripts that have been used to select a robust yet efficient combination of integrator and propagator for the simulation of a Martian orbit during one year.

## Benchmark
First, [1_year_rk_4.py](1_year_rk_4.py) establishes a benchmark of the orbit during one year. It does so using a Cowell propagator, and a Runge Kutta 4 integrator with a fixed step of 10 seconds.
This simulation takes around 110 seconds to run on my CPU, for a simulated orbit of 1 year.
Once this benchmark simulation is done, the resulting altitudes as a function of time are saved in [rk_4_baseline.dat](rk_4_baseline.dat).

The same benchmark has then been run with the inclusion of densities from the Mars Climate Database. The propagation has then been run for 50 days, 20 days, and 1 day. The results have also been saved in `.dat` files.

## Integrator selection
Various integration scheme and settings have been tested, as to reduce the simulation time without introducing significant deviations w.r.t. the benchmark.
This has been done in [1_year_integrators.py](1_year_integrators.py).

The comparison of the different integrators can be seen in the following table.

| Integration scheme | Settings changed                                           | Simulation time [s] | Maximum difference in altitude [km] |
|--------------------|------------------------------------------------------------|---------------------|-------------------------------------|
| RK4                | Step of 10s                                                | 104.6                 | 0.0 (benchmark)                   |
| RK4                | Step of 30s                                                | 34.5                  | 11.0                              |
| RK4                | Step of 60s                                                | 17.0                  | 70.0                              |
| RKF45              | Step of 1-300s<br> Tolerances of 1E-9                      | 40.5                  | 5.0                               |
| RKF45              | Step of 1-300s<br> Tolerances of 1E-8                      | 26.2                  | 48.2                              |
| RKF45              | Step of 1-300s<br> Tolerances of 1E-6                      | 11.4                  | 129.7                             |
| RKF45              | Step of 1-500s<br> Tolerances of 1E-9                      | 40.7                  | 5.0                               |
| RKF56              | Step of 10-300s<br> Tolerances of 1E-9                     | 27.2                  | 31.6                              |
| RKF78              | Step of 10-300s<br> Tolerances of 1E-9                     | 17.9                  | 8.4                               |
| RKDP87             | Step of 10-300s<br> Tolerances of 1E-9                     | 15.0                  | 0.34                              |
| RKDP87             | Step of 10-300s<br> Tolerances of 2.5E-8                   | 10.5                  | 0.63                              |
| RKDP87             | Step of 10-300s<br> Tolerances of 1E-8                     | 11.9                  | 0.49                              |
| RKDP87             | Step of 10-300s<br> Tolerances of 1E-7                     | 10.0                  | 0.69                              |
| ABM                | Step of 10-300s<br> Tolerances of 1E-9 <br>Order of 6-11   | 44.2                  | 0.35                              |
| ABM                | Step of 30-300s<br> Tolerances of 1E-9 <br>Order of 6-11   | 20.3                  | 25.1                              |
| ABM                | Step of 10-500s<br> Tolerances of 1E-9 <br>Order of 6-11   | 43.4                  | 0.35                              |
| BS                 | Step of 10-500s<br> Tolerances of 1E-9 <br> Max 5 steps    | 19.2                  | 1.28                              |
| BS                 | Step of 10-500s<br> Tolerances of 1E-9 <br> Max 4 steps    | 21.6                  | 0.89                              |

In light of these results, the current recommendation is to use a RKDP87 integrator with a step size varying between 10s and 300s, with a tolerance of 2.5E-8.
This leads to a low propagation time of 10.5 seconds, with a deviation of only 630m compared to the benchmark.

## Propagator selection
All of the propagators available in the simulation framework have been tested, using the integrator recommended above.
This has been done in [1_year_propagator.py](1_year_propagator.py), first without including the density from the MCD, and only using Mars as a point mass and the aerodynamic drag as the accelerations.
This lead to the results of the table below.

| Propagator                                             | Simulation time [s] | Maximum difference in altitude [km] | Comment                                      |
|--------------------------------------------------------|---------------------|-------------------------------------|----------------------------------------------|
| Cowell                                                 | 10.5                | 0.65                                |                                              |
| Encke                                                  | 11.4                | 0.66                                | Singularity for eccectricity of 0.           |
| Gauss Keplerian                                        | 11.5                | 0.35                                | Singularity for inclination of 0 deg.        |
| Gauss Modified Equinoctial                             | 10.1                | 0.35                                | Singularity for inclination of 0 or 180 deg. |
| Unified State Model with Quaternions                   | 10.4                | 0.35                                |                                              |
| Unified State Model with Modified Rodrigues Parameters | 10.5                | 0.35                                |                                              |
| Unified State Model with Exponential Map               | 10.6                | 0.35                                |                                              |

Then, using the density from the MCD, keeping only Mars as a point mass and the aerodynamic drag as accelerations, and simulating a 20 days orbit, the following table has been made.
Because this model results in highly more varying densities than the exponential model, the steps taken by the variable step integrator are much smaller.


| Propagator                                             | Simulation time [s] | Maximum difference in altitude [m] |
|--------------------------------------------------------|---------------------|------------------------------------|
| Cowell                                                 | 10.92               | 6.883                              |
| Encke                                                  | 10.62               | 4.616                              |
| Gauss Keplerian                                        | 6.57                | 4.887                              |
| Gauss Modified Equinoctial                             | 4.54                | 4.968                              |
| Unified State Model with Quaternions                   | 4.83                | 4.911                              |
| Unified State Model with Modified Rodrigues Parameters | 4.79                | 4.923                              |
| Unified State Model with Exponential Map               | 4.71                | 4.928                              |

From the table above, the Gauss Modified Equinoctial propagator has been used to select the appropriate environment to be used. This was made in [this file](../environments/test_req_envs.py).

Lastly, the same process has been repeated, but adding the Martian spherical harmonics up to degree and order 4 to the gravitational acceleration, and adding a cannonball radiation pressure. This lead to the table below.

| Propagator                                             | Simulation time [s] | Maximum difference in altitude [m] |
|--------------------------------------------------------|---------------------|------------------------------------|
| Cowell                                                 | 10.31               | 7.835                              |
| Encke                                                  | 10.87               | 8.059                              |
| Gauss Keplerian                                        | 6.73                | 3.754                              |
| Gauss Modified Equinoctial                             | 4.76                | 5.182                              |
| Unified State Model with Quaternions                   | 9.46                | 0.753                              |
| Unified State Model with Modified Rodrigues Parameters | 9.36                | 0.787                              |
| Unified State Model with Exponential Map               | 9.22                | 0.781                              |

The current recommendation is then to use the Gauss Modified Equinoctial propagator. However, if CPU time is not an issue, the Unified State Model with Quaternions propagator leads to the most accurate results.

## Overall 
Finally, [1_year_best_combo.py](1_year_best_combo.py) explores finer tuning of the integrator, using the suggested propagator.
Note that the change in settings has been done manually, starting from the integrator and propagator suggested in the previous steps

In this process, caution has been payed so that the integrator is limited by its tolerance and not the step size.
Also, irregulars jumps in altitude have sometime been observed. The integrator has been tuned to avoid these.

The settings that resulted from this study are thus the followings:

* A Runge Kutta Dormant Prince 87 integrator with:
    * A relative and absolute tolerance of 5E-8
    * A step size range of 0.05s to 600s
    * An initial step size of 150s
    * A minimum and maximum factor increase of 0.15 and 3.5

    * A safety factor of 0.6
* A Gauss Modified Equinoctial propagator
    * Inclinations of 0 deg and 180 deg are not possible. They can be manually changed to 0.01 deg or 179.99 deg.

This leads to a simulation time of 5 seconds, with a deviation from the benchmark of a maximum of 3.239 meters in altitude after 20 days of simulation.

For a more accurate simulation, the same settings can be used but with a relative and absolute integrator tolerance of 1E-12 instead, leading to a simulation time of 10.704 seconds and a deviation of only 0.746 meters.

Finally, these settings may be tweaked again later on when accelerations such as drag or thrust are changed.