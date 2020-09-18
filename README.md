# SPectral element based EARthquake cycle simulator (SPEAR) 

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**Insert Badges here**

Numerical simulation of long-term earthquake cycles on a two-dimensional vertical strike-slip fault with dynamic treatment of inertial effects. Written in [Julia](https://julialang.org).

## Current Status
Under development. Working on MPI based parallelization.

## Features
- Antiplane strain deformation on planar fault with SH waves
- Fully dynamic treatment of seismic events
- Rate-and-state-dependent friction with aging laws 
- Adaptive time-stepping to switch between interseismic and seismic events
- Customizable to include off-fault and on-fault heterogeneites
- Complex geometry of fault damage zones including gaussian and trapezoid fault damage zones
- Time-dependent healing of fault damage zones constrained by observations

## Dependencies
- [Julia version >= 1.1](https://julialang.org)
- [Algebraic Mulltigrid](https://github.com/JuliaLinearAlgebra/AlgebraicMultigrid.jl)
- [Iterative Solvers](https://github.com/JuliaMath/IterativeSolvers.jl)
- [FEMSparse](https://github.com/ahojukka5/FEMSparse.jl)

## Quickstart guide

## Documentation

## Screenshots


Please cite the following if using this code:

[Thakur, P., Huang, Y., & Kaneko, Y. Effects of low‐velocity fault damage zones on long‐term earthquake behaviors on mature strike‐slip faults. Journal of Geophysical Research: Solid Earth, e2020JB019587.](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020JB019587)

[Kaneko, Y., Ampuero, J. P., & Lapusta, N. (2011). Spectral‐element simulations of long‐term fault slip: Effect of low‐rigidity layers on earthquake‐cycle dynamics. Journal of Geophysical Research: Solid Earth, 116(B10).](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2011JB008395)
