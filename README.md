# Lattice Boltzmann Method
<!--toc:start-->
- [Lattice Boltzmann Method](#lattice-boltzmann-method)
  - [Introduction](#introduction)
  - [Overview](#overview)
  - [The Time Step Algorithm](#the-time-step-algorithm)
  - [Run the Simulation](#run-the-simulation)
    - [Build](#build)
    - [Running the Simulation](#running-the-simulation)
    - [Visualize the results](#visualize-the-results)
  - [Report](#report)
  - [Conclusions](#conclusions)
<!--toc:end-->

This hands on objective is to implement the Lattice Boltzmann Method and applying it to solve the 2D lid cavity problem.

The velocity set used is the common D2Q9 and as collision operator the Bhatnagar-Gross-Krook (BGK) was selected.

> The system therefore is the 2D lid cavity problem solved using LBGK algorithm on a D2Q9 velocity set.

## Introduction

The lattice Boltzmann equation originated from the description of the dynamics of a gas on a mesoscopic scale;
though it is possible to define a relation with the equations of fluid dynamics on the macroscale.
Therefore, from a solution of the Boltzmann equation for a given case it is usually possible to find a solution to the
Navier-Stokes equations for the same case[^1].

[^1]: Since the Boltzmann equation is more general, it has also solutions that do not correspond to Navier-Stokes solutions.

Even though the Boltzmann equation is analitically more difficult to solve than the NSE, its numerical scheme is
simpler and this is why it is frequently used in numerical problems of CFD simulations.

## Overview

The lattice Boltzmann method revolves around the _discrete-velocity distribution function_ $f_{i}(x,t)$, 
often called the particle _populations_.

This density function is used to describe how the fluid behaves on a mesoscopic scale,
through the modeling of collisions of particles and their movement between different _cells_.

The discretized version of the lattice Boltzmann equation in velocity and time is:

$$f_{i}(x + c_{i}\Delta{t},t + \Delta{t})=f_{i}(x,t) + \Omega_{i}(x,t)$$

in which:

- $c_{i}$ is the speed at which particles move to a neighbouring point;
- $\Omega_{i}$ is the collision operator.

As previously stated, the collision operator used in this hands-on is the Bhatangar-Gross-Krook (BGK) operator:

$$ \Omega_{i}(f) = -\frac{f_{i}-f^{eq}_{i}}{\tau} \Delta{t}$$

with $\tau$ the reaxation time and $\Delta{t}$ the time step.

The equilibrium is given by:

$$ f^{eq}_{i}(x,t) = 
w_{i} \rho (1 + \frac{u * c_{i}}{c_{s}^2}  
+ \frac{(u * c_{i})^2}{2c^4_{s}} - \frac{u*u}{2c^2_{s}}$$

in which:

- $w_{i}$ is the weight and $\rho$ the density
- $c_{i}$ represents the set of discrete particle velocities;
- $u$ represents the macroscopic fluid velocity;
- $c^2_{s} = (1/3)\frac{\Delta{x}^2}{\Delta{t}^2}$ determines the relation between the pressure and the density (in basic isothermal LBE).

> The latter of these factor need not to be explicitly calculated in the algorithm.

### The Time Step Algorithm

The LBM algorithm consists of a cyclic sequence of operations; each cycle
corresponds to one time step. The substeps are:

1. Perform the streaming;
2. Apply boundary conditions;
1. Compute the macroscopic _moments_ $\rho(x,t)$ and $u(x,t)$ from $f_{i}(x,t)$;
4. Compute the equilibrium $f^{eq}_{i}(x,t)$;
5. Compute the collision operation;
6. Write the macroscopic fields.
7. Increase time step and check convergence.

### Key Parameters

- Reynolds Number: $Re$
- Kinematic viscosity: $\nu =\frac{u_{lid}L_x}{Re}$ 
- Relaxation time: $\tau = 3\nu +0.5$

## Parallelization
!!!!ADD !!!

## Run the Simulation

### Build

To build the executable run:

`cmake -B build/ && make -C build/ [-j]`

### Running the Simulation

To run, enter in directory `build/` and run:

`./lbm-2-lbm <grid_num_cells_x> <grid_num_cells_y> <reynold_number> <lid_initial_velocity>`

### Visualize the results



## Report

!!!!!!!!!!!! ADD THIS SECTION !!!!!!!!!!!!!!!

## Conclusions

!!!!!!!!!!!! ADD THIS SECTION !!!!!!!!!!!!!!

---

[Alessandro Frisone](https://github.com/DatemiUn30L) |
[Corrado Sciancalepore](https://github.com/CorradoSciancalepore) |
[Elisa Antonioli](https://github.com/ElisaAntonioli) |
[Chiara Nonino](https://github.com/ChiaraNonino) |
[Gabriele Santandrea](https://github.com/Gab-San)
