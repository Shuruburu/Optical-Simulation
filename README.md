# Jan Sokal

# Simulation Guide: Optomechanical Coupling to the HBAR

This guide accompanies the code developed for my master's thesis. The simulations
address the optimization of the optomechanical coupling to a High-overtone Bulk
Acoustic Resonator (HBAR).

## Background

Using the transfer-matrix formalism, I investigated how the beam waist of a
Gaussian beam inside an optical cavity depends on the cavity geometry. All
simulations assume a fixed cavity length of **11.84 mm**. The following
parameters were varied:

- Radius of curvature of the mirror
- Position of the crystal inside the cavity
- [Other parameters, e.g. crystal thickness, refractive index. See Chapter X of the thesis]

For the full derivation and physical motivation, please refer to the thesis
(Chapter X).

## What this guide covers

1. What the code does, and how it is structured
2. Which parameters to change, and where, to obtain the desired results
3. [Optional: how to reproduce the figures in the thesis]



## Configuring the simulation: `tracks.py`

Most of the constants can be changed in the file `tracks.py`, including:

- the radius of curvature of the cavity mirrors
- the placement of the crystal
- the speed of sound
- [other parameters, e.g. crystal thickness, refractive index]

The first function in this file is `Track`. It is the base of the entire
simulation, so it is worth understanding before anything else.

`Track` returns a dictionary that describes the optical cavity, with one entry
for each element (mirrors, crystal, antenna, etc.). Storing everything in a
single dictionary gives the rest of the code one data structure to read from,
so the parameters are defined in one place and used consistently throughout.

To run the simulation with different parameters, edit the values in `Track`
and re-run main.py.

The important thing is the sequence of the element the program reads things in an order, meaning each element comes after each other. 
You need to make sure that the positions follow right after one another. An example, would be if vacum ends then the crytsal starts.

## Configuring the simulation : `Material_constants`

This section contains all of the material constants used in the simulation.
The elastic constants are stored in Voigt notation.

### Adding a new material

1. Add a new entry for your material, following the same structure and
   Voigt-notation format as the existing ones.
2. In the `coupling_approximation` function, change the material name to the
   name of your new entry.

The simulation will then use the new constants in its calculations.


## GRIN lens track

This track is used to investigate the coupling when the cavity contains a
different type of medium: a gradient-index (GRIN) lens.

For the definition of the GRIN lens and the physics behind it, see
Appendix A. g is the grin index lens and l is the length, 
usually you can pass the grin index and using function pitch to length get the length.


## Other functions
    Other function serve as helping functions, you need to read the comments to understand what they do but they are rather self explenatory

# Gaussian beam: file 

This file contains most of the algorithms used to calculate the beam waist
inside the cavity. The waist can be calculated in three ways. [Name the other
two here, or describe them in the sections below.]

### 1. Analytic method

The analytic method is implemented in the class `Analytic`. It works as follows:

1. Start from an initial beam vector.
2. Propagate it through multiple cavity round trips using the `round_trip`
   function. Two round trips are enough in principle; more are used as a check
   that the vector is self-consistent (i.e. it returns to itself after a
   round trip).
3. From the converged vector, extract the beam parameters.

To obtain the results:

1. Call `solve`.
2. Call `get_waist` to get the beam waist.
3. Call `get_radius_of_curvature` to get the radius of curvature of the beam.

### 2. Numerical method 

This method is contained entirely in the class `Beam_evolution`. It searches
for the self-consistent beam, i.e. the beam whose q parameter reproduces
itself after one round trip (see [Section 2 of the thesis] for the definition
of the q parameter).

The method works as follows:

1. Choose a starting beam waist and construct the corresponding q parameter.
2. Propagate the beam through one round trip of the cavity.
3. Compare the q parameter after the round trip with the initial one and
   compute the error between them.
4. Repeat steps 1 to 3 for a range of starting waists.
5. Select the waist with the lowest error. This is the beam waist of the
   cavity mode.

The accuracy of the result depends on the range and resolution of the
starting waists.  look into the function if u need to change it.


### Numerical method( matching the radi of curvature of the back mirror )
This method is contained entirely in the class `Beam_evolution`. It searches
for the self-consistent beam, i.e. the beam whose q parameter reproduces
itself after one round trip (see [Section 2 of the thesis] for the definition
of the q parameter).

Same as in the numerical above just a different criteria for the selection




### 3. Step-by-step propagation

The file also contains the function `propagate`, which works differently from
the two methods above. Instead of only finding the beam waist, it tracks how
the beam evolves along the whole cavity.

The function works as follows:

1. Create an array of q parameters at points that are equally spaced in
   optical distance.
2. Advance the beam by one small step at a time using the propagation
   function.
3. At each step, call `check_medium` to determine which medium the beam is
   currently in, so that the correct propagation is applied there.

The result is the beam parameters at every point along the cavity. These are
needed for plotting the beam profile and for the numerical approximation of the
optomechanical coupling.

## Optomechanical coupling

This function calculates the full optomechanical coupling to the HBAR. For the
theory behind the calculation, see Section 2.4 of the thesis.

The calculation assumes the optimal conditions described in the thesis
[add the section number if different from 2.4]. The material properties are
taken from the material constants (see [Material constants](#material-constants)).

This function is called by `run_signal`, so in most cases you do not need to
call it directly: run `run_signal` instead.

## `run_signal`: running the full simulation

`run_signal` is the main function for running the simulation. It:

1. Adjusts the position of the crystal inside the cavity.
2. Runs the full coupling approximation for each position, in parallel.
3. Plots the results.


## file scripts

This file contains ready-made scripts that make the software easier to use.
The most important ones are:

- **Cavity waist**: calculates the beam waist in the cavity using all of the
  methods described above and plots it. The plotting can be turned off if you
  only need the values.
- **Coupling approximation**: simulates the coupling strength for different
  mirrors.
- **Parameter sweep**: runs the simulation for different parameters of the
  crystal.
- **GRIN lens track**: calculates the distance required for specific GRIN
  lenses and saves the results to a JSON file.

  There are more functions but you just need to read the code in case again u can drop me a message!

In case of any question, you can text me on instagram(I may change my number on whatsapp so I give the insta) @yan_s_advantures
Enjoy!