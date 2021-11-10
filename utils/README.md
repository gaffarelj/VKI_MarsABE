# Utilities
This folder contain the utilities that are directly related to the orbital simulation.

## Propagation
The [propagation.py](propagation.py) module offers an interface to simulate the orbit of a satellite, working on top of the Tudat(Py) interface.

### Orbit simulation
First, [propagation.py](propagation.py) contains a `orbit_simulation` class that can be used to setup and simulate a satellite orbit.

#### Initialization
The `orbit_simulation()` class must be initialized with the following inputs:
 * `sat`: object of the type `satellite`, defining the satellite of which the orbit is simulated (see [below](#satellite-models)).
 * `central_body`: name of the body around which the satellite orbits (for instance: `Earth`, `Mars`, or `Sun`).
 * `sim_time`: number of seconds for which the simulation is to be run.
 * `init_time`: time at which to start the simulation, in seconds since J2000.
 * `verbose` (optional): if False, do not let Tudat print to the terminal and hide non-critical warnings and information.
 * `save_power` (optional): by default, the power received by the satellite is not saved. Set to True to save it.

Two empty dictionaries are also declared during initialization to contain, if specified, the solar irradiance and the solar power.

#### Bodies
All bodies that will be relevant during the simulation must be created, as well as their relevant environmental model (such as their atmosphere).

This can be done using the `create_bodies()` method. By default, the body around which the satellite orbits will always be setup.

Also, during bodies setup, the frame orientation is set to be `ECLIPJ2000`, and the origin of the frame is the central body.

The following inputs can be used:
 * `additional_bodies`: list of bodies to setup in addition to the central body. By default, this list is set to always setup Jupiter and the Sun. **Warning**, do not forget to include all of the bodies of which you want to include effects such as the acceleration, or compute values. Otherwise, they will simply not exist.
 * `use_MCD`: list of two booleans. The first one is used to indicate whether the Mars Climate Database (MCD) atmosphere model should be implemented. With the first one True, the second boolean is used to indicate whether winds should also be added from the MCD.
 * `preload_MCD`: setting this to True will immediately trigger all of the MCD files to be loaded (which can take up to a minute). This way, they do not need to be loaded on the fly during the propagation.
 * `save_MCD_vals`: by setting this to True, all of the values returned from the MCD module will be saved.
 * `use_GRAM`: if True, the Mars GRAM 2010 atmospheric model is used to get density of the Martian atmosphere.

Do note that, if the Martian atmosphere is used and neither `use_MCD` nor `use_GRAM` are set to True, an exponential atmospheric model based on Mars GRAM is used.

No atmospheric model has been implemented for bodies other than Mars. Trying to access one will then raise a `NotImplementedError`.
The code implemented for the Martian atmosphere can be copied and adapted to another planet. Also, the documentation to setup an atmospheric model in Tudat(Py) can be accessed [here](https://tudatpy.readthedocs.io/en/latest/atmosphere.html).

In the `create_bodies()` method, the drag coefficients (that can be function of the altitude) of the satellite, if any, are implemented.

Finally, it is worth noting that, if the cannonball radiation pressure model is used, the reference surface area has been set to a default of `0.125 m2`, and the radiation pressure coefficient to a default of 1.2. This can easily be changed.

#### Initial state
The initial state (position and velocity) of the satellite defines its initial orbit.

This can be setup using the following inputs to the `create_initial_state()` function:
 * `a`: semi-major axis of the orbit in meters
 * `h_p`: geopotential altitude of the satellite periapsis in meters
 * `e`: eccentricity of the orbit in [0-1]
 * `i`: inclination of the orbit in radians in [0-pi/2]
 * `omega`: argument of periapsis of the orbit in radians in [0-pi]
 * `Omega`: longitude of the ascending node in radians in [0-2pi]
 * `theta`: true anomaly in radians in [0-2pi]. Note that this defines where the satellite starts in its orbit, and can most of the time be ignored

The function automatically converts these orbital elements to a cartesian initial state (`orbit_simulation.initial_state`).

Also, one of the `a` or `h_p` inputs must be specified. If the semi-major `a` axis is not specified but `h_p` is, `a` is computed from the eccentricity `e` and the periapsis altitude `h_p`.

#### Accelerations
All of the accelerations on the satellite are created in the `create_accelerations()` function.

This can be done in two ways: either by loading a default acceleration configuration, or by setting up the individual accelerations.

A default configuration can be loaded by setting `default_config` to 0, 1, or 2:

 0. Use the central body gravitational acceleration as a Point Mass and the aerodynamic acceleration due to the atmosphere of the central body.
 1. Use the central body gravitational acceleration as Spherical Harmonics up to Degree/Order 4, the aerodynamic acceleration due to the atmosphere, and the Solar radiation using the cannonball radiation pressure model.
 2. Use the central body gravitational acceleration as Spherical Harmonics up to Degree/Order 8, the aerodynamic acceleration due to the atmosphere of the central body, the Solar radiation using the cannonball radiation pressure model, and the gravitational acceleration as a Point Mass of the Sun and Jupiter.

The individual accelerations (in Tudat(Py) types) can be input to `env_accelerations` as a list. This can be made easier by using the `env_acceleration()` class provided below, as discussed [here](#environmental-accelerations-definition).

Finally, one of the thrust models can be added to the satellite by setting the `thrust` input to the index corresponding to the thrust model, as described [below](#thrust-models).

#### Integrator
Two ways have been thought to setup the integrator used during the orbital simulation.

Firstly, the `create_integrator()` method can be used to create a RK-DP 87 variable step integrator.
It can be tuned using the following inputs:
 * `tolerance`: integrator tolerance, this implicitly tunes the integration steps.
 * `dt`: list of time step settings (all in seconds):
   0. Minimum time step.
   1. Initial time step.
   2. Maximum time step.
 * `sf`: list of safety factor settings:
   0. Initial safety factor.
   1. Minimum increase factor.
   2. Maximum increase factor.

Secondly, if one wishes to use a different integrator, `orbit_simulation.integrator_settings` can directly be set to one of the Tudat(Py) integrator described on [this page](https://tudatpy.readthedocs.io/en/latest/integrator.html).

#### Termination settings
The termination settings consist in parameters that, when reached, trigger the end of the simulation.

This can be setup by setting `orbit_simulation.termination_settings` to termination settings as specified on [this page](https://tudatpy.readthedocs.io/en/latest/propagator.html).

Alternatively, the `create_termination_settings()` function can be used. It has the following input settings:
 * `min_altitude`: critical altitude under which the simulation must stop, defaulted to 50km.
 * `max_altitude`: critical altitude above which the simulation must stop, defaulted to 1000km (set to a ridiculously high value to ignore this limit).
 * `cpu_time`: CPU time after which to stop the simulation, in seconds, defaulted to 1hr.

#### Dependent variables
Numerous variables can be saved during the simulation, sometime crucial for analysis afterwards.

The `create_dependent_variables()` function can be used by simply setting the list of dependent variables to be saved, `to_save`, amongst the followings:
 * `h`:     altitude of the satellite in meters.
 * `rho`:   atmospheric density at the position of the satellite in kg/m3.
 * `V`:     airspeed of the Satellite in m/s (size of 3).
 * `m`:     mass of the satellite in kg.
 * `F_T`:   acceleration due to the Thrust in N (size of 3).
 * `D`:     acceleration due to the atmosphere in N (size of 3).
 * `C_D`:   aerodynamic coefficients (size of 3).
 * `r_cb`:  relative position of the central body w.r.t. the sun in m (size of 3).
 * `Kep`:   Keplerian state of the satellite (size 6: a, e, i, omega, Omega, theta).
 * `h_p`:   altitude of the periapsis of the satellite.
 * `lat`:   latitude of the satellite.
 * `lon`:   longitude of the satellite.

#### Propagator
Both the position in 3D and the mass of the satellite can be propagated.

The position is always propagated, while the mass is optional and can be enabled by setting `prop_mass=True` when calling `create_propagator()`.

By default, the propagator used to propagate the position is an Encke propagator, with its state being the difference between the satellite position and the one it would be on with a perfect un-disturbed Keplerian orbit. This propagator thus keeps the propagated state small, helping prevent deviations due to numerical error.

A different propagator can be used by using the `propagator` input of `create_propagator()` to another one from [this list](https://tudatpy.readthedocs.io/en/latest/propagator.html#tudatpy.numerical_simulation.propagation_setup.propagator.TranslationalPropagatorType).

#### Simulation
After the simulation setup is finished, it can finally be run by calling `simulate()`. This function returns three elements:
 * `sim_times`: the times in the simulation.
 * `states`: the satellite state history.
 * `dep_vars`: the dependent variables history.
These three elements are also saved as `orbit_simulation.sim_times`, `orbit_simulation.states`, and `orbit_simulation.dep_vars`.

#### Reading dependent variables
To access the dependent variables more easily, the `get_dep_var()` function can be used.

As its sole input, `key` should be one of the dependent variable names specified in [the dependent variables section](#dependent-variables).

### Environmental accelerations definition
Then, [propagation.py](propagation.py) contains a `env_acceleration` class that can be used to define environmental accelerations models.

This class is to be used as follows:
```
acceleration_Mars_on_sat = env_acceleration("Mars", SH=True, SH_do=[8,8], aero=True)
acceleration_Sun_on_sat = env_acceleration("Sun", PM=True, rad=True)
```
In this example, `acceleration_Mars_on_sat` and `acceleration_Sun_on_sat` will then be two lists, containing the following:
 * `acceleration_Mars_on_sat`:
   * the gravitational acceleration of Mars as Spherical Harmonics (`SH=True`), up to Degree 8 and Order 8 (`SH_do=[8,8]).
   * the aerodynamic acceleration caused by the atmosphere of Mars (`aero=True`). This uses the environment and satellite properties that must be setup beforehand.
 * `acceleration_Sun_on_sat`:
   * the gravitational acceleration of the Sun as a Point Mass (`PM=True`).
   * the radiation pressure of the Sun approximating the satellite as a cannonball (`rad=True`). This uses the satellite properties that must be setup beforehand.

All of the acceleration models in these lists are defined using Tudat(Py).

## Satellite models
All of the parameters related to the satellite itself have been compiled in the `satellite` class of [sat_models.py](sat_models.py) file.

These are its mass(es), drag coefficient(s), solar array, and motor properties.

This class can be initialized with the following inputs, of which only the three first have no default and must be specified:
 * `name`: name of the satellite.
 * `mass`: mass of the satellite.
 * `Cd`: drag coefficient of the satellite, either a single float (constant Cd), or a list (Cd at given altitudes).
 * `Cd_h`: list of altitudes (in meters) at which the Cd have been specified.
 * `prop_mass`: propellant mass in the satellite in kg, defaulted to 0.
 * `S_ref`: reference surface area in m2, defaulted to 0.01 m2.
 * `SA_areas`: solar panel areas in the x/y/z planes, in m2 (more explanation [here](https://github.com/gaffarelj/VKI_MarsABE/tree/main/SPARTA#satellite-configurations)).
 * `SA_frac`: fraction of the solar panel areas that actually collects solar radiation, defaulted to 0.7042.
 * `SA_eff`: solar panel efficiency, defaulted to 0.29.
 * `EPS_eff`: Electric Power System efficiency, defaulted to 0.89.
 * `S_t`: area of the throat at the end of the air-breathing inlet in m2.
 * `comp_ratio`: ratio between the free-stream density and the density at the end of the air-breathing inlet.

The [sat_models.py](sat_models.py) script also contains two dictionaries of pre-defined satellites: `satellites` and `satellites_with_tank`.

## Thrust models