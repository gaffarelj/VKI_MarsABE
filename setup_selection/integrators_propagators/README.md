# Integrators and propagators selection
This folders contains different various scripts that have been used to select a robust yet efficient combination of integrator and propagator for the simulation of a Martian orbit during one year.

## Benchmark
First, [1_year_rk_4.py](1_year_rk_4.py) establishes a benchmark of the orbit during one year. It does so using a Cowell propagators, and a Runge Kutta 4 integrator with a fixed step of 10 seconds.
This simulation takes around 110 seconds to run on my CPU.
Once this benchmark simulation is done, the resulting altitudes as a function of time are saved in [rk_4_baseline.dat](rk_4_baseline.dat).

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
This leads to a low propagation time of 10.5, with a deviation of only 600m compared to the benchmark.

## Propagator selection
All of the propagators available in the simulation framework have been tested, using the integrator recommended above.
This has been done in [1_year_propagator.py](1_year_propagator.py), and lead to the results of the table below.

| Propagator                                             | Simulation time [s] | Maximum difference in altitude [km] | Comment *(to be added)* |
|--------------------------------------------------------|---------------------|-------------------------------------|-------------------------|
| Cowell                                                 | 10.5                | 0.6                                 |                         |
| Encke                                                  | 11.5                | 0.6                                 |                         |
| Gauss Keplerian                                        | 11.5                | 5                                   |                         |
| Gauss Modified Equinoctial                             | 10                  | 5                                   |                         |
| Unified State Model with Quaternions                   | 10                  | 0.05                                |                         |
| Unified State Model with Modified Rodrigues Parameters | 10.5                | 0.05                                |                         |
| Unified State Model with Exponential Map               | 10.5                | 0.05                                |                         |

From these, it is advised to use the Unified State Model propagator, most likely the one using Quaternions as the state.

## Overall 
Finally, [1_year_best_combo.py](1_year_best_combo.py) explores finer tuning of the integrator, using the suggested propagator.
This lead to the findings of the following table.

| Change in integrator settings                 | Simulation time [s] | Maximum difference in altitude [km] | Comment                                                             |
|-----------------------------------------------|---------------------|-------------------------------------|---------------------------------------------------------------------|
| Decrease the tolerance to 2.5E-1 (extreme).   | 10.5                | 0.6                                 | This shows that the integrator is limited by the maximum step size. |
| Increase the maximum step from 300s to 1800s. | 1.7                 | 3.6                                 |                                                                     |
| Set the tolerance back to 2.5E-8              | 3.5                 | 1.2                                 | The difference increases with increasing perturbations.             |
| Set the tolerance to 1E-9                     | 5                   | 0.6                                 | Big jump in simulated altitude are seen.                            |
| Set the tolerance to 5E-9                     | 4.2                 | 0.9                                 | Big jump in simulated altitude are seen.                            |
| Change the safety factor from 0.75 to 0.85    | 4.3                 | 0.9                                 | Size of the jump in altitude lowered.                               |
| Change the safety factor from 0.85 to 0.95    | 4.3                 | 0.9                                 | No more altitude jumps are present.                                 |
| Change the tolerance back to 1E-9             | 5                   | 0.6                                 |                                                                     |

The settings that resulted from this study are thus the followings:

* A Runge Kutta Dormant Prince 87 integrator with:
    * A step size range of 10s to 1800s
    * A relative and absolute tolerance of 1E-9
    * A safety factor of 0.95
* An Unified State Model propagator with Quaternions

This leads to a simulation time of 5 seconds, with a deviation from the benchmark of a maximum of 600 meters in altitude.