# SPEAR: SPectral element based EARthquake cycle simulator 

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![Build Status](https://travis-ci.com/thehalfspace/Spear.svg?branch=master)](https://travis-ci.com/thehalfspace/Spear)

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
1. Install Julia >= 1.1 (preferably the latest version).
2. Run `install_dependencies.jl` script to install the specific packages.
3. Open `run.jl` and edit the output file name and the resolution of the problem.
4. Run `run.jl` from the terminal or IDE of your choice.
5. The output files are saved in `data/$(simulation_name)`.
6. You can look at the output file from `analyze_results.jl`.
7. For basic testing run `julia tests/basic_test_01.jl` to look at two plots saved in the `plots/$(simulation_name)` directory.

## Documentation
- Coming soon

## Screenshots
<img src="https://github.com/thehalfspace/Spear/blob/master/plots/example/cumulative_slip.png" alt="drawing" width="600"/>
<img src="https://github.com/thehalfspace/Spear/blob/master/plots/example/Vfmax01.png" alt="drawing" width="600"/>


## References
Please cite the following if using this code:

[Thakur, P., Huang, Y., & Kaneko, Y. Effects of low‐velocity fault damage zones on long‐term earthquake behaviors on mature strike‐slip faults. Journal of Geophysical Research: Solid Earth, e2020JB019587.](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020JB019587)

[Kaneko, Y., Ampuero, J. P., & Lapusta, N. (2011). Spectral‐element simulations of long‐term fault slip: Effect of low‐rigidity layers on earthquake‐cycle dynamics. Journal of Geophysical Research: Solid Earth, 116(B10).](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2011JB008395)
