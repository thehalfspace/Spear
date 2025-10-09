# SPEAR: SPectral Element EARthquake Simulator - Function Reference

## Overview

SPEAR is a spectral element method simulator for earthquake cycle dynamics on a 2D vertical strike-slip fault with fully dynamic treatment of seismic events. The code implements rate-and-state friction laws with adaptive time-stepping for both interseismic and coseismic periods.

## Setup and Entry Point

### `run.jl` - Main execution script
- Entry point for simulations
- Sets simulation resolution and output directory
- Calls `setParameters()` to initialize simulation parameters
- Executes `main()` function to run the simulation

### `par.jl` - Parameter setup
- **`setParameters(FZdepth, resolution)`** - Main setup function that initializes all simulation parameters
  - Sets domain dimensions (LX=48km depth, LY=30km width)
  - Configures mesh resolution and element sizes
  - Defines time parameters (150 years total simulation time)
  - Sets up material properties and fault parameters
  - Returns structured parameter sets for simulation

## Initial Conditions and Material Properties

### `src/initialConditions/defaultInitialConditions.jl`
- **`fricDepth(FltX)`** - Sets rate-state friction parameters (a, b) as functions of depth
- **`SeffDepth(FltX)`** - Defines effective normal stress profile with depth
- **`tauDepth(FltX)`** - Sets initial shear stress distribution along fault
- **`Int1D(P1, P2, val)`** - Linear interpolation utility function

### `src/MaterialProperties.jl`
- **`material_properties(NelX, NelY, NGLL, dxe, dye, ThickX, ThickY, wgll2, rho1, vs1, rho2, vs2)`** - Sets up rectangular fault damage zone with reduced rigidity
- **`rigid(x, y)`** - Defines material properties for trapezoidal damage zone geometry
- **`mat_trap(NelX, NelY, NGLL, iglob, M, dxe, dye, x, y, wgll2)`** - Assembles mass matrix for trapezoidal damage zones
- **`line(x, y)`** - Geometric function defining trapezoid damage zone boundaries

## Mesh Generation and Assembly

### `src/MeshBox.jl`
- **`MeshBox!(NGLL, Nel, NelX, NelY, FltNglob, dxe, dye)`** - Generates 2D spectral element mesh
  - Creates global node connectivity (iglob)
  - Generates global coordinates for all GLL nodes
  - Handles element-to-element connectivity for rectangular domain

### `src/GetGLL.jl`
- **`GetGLL(ngll)`** - Loads precomputed Gauss-Legendre-Lobatto (GLL) quadrature points, weights, and derivative matrices from tabulated data files

### `src/Assemble.jl`
- **`Massemble!(NGLL, NelX, NelY, dxe, dye, ThickX, ThickY, rho1, vs1, rho2, vs2, iglob, M, x, y, jac)`** - Assembles global mass matrix and computes stable time step

### `src/Kassemble.jl`
- **`stiffness_assembly(NGLL, NelX, NelY, dxe, dye, nglob, iglob, W)`** - Assembles global sparse stiffness matrix
- **`K_element(W, dxe, dye, NGLL, H, Nel)`** - Computes element-level stiffness matrices
- **`FEsparse(Nel, Ke, iglob)`** - Converts local element matrices to global sparse format

### `src/BoundaryMatrix.jl`
- **`BoundaryMatrix!(NGLL, NelX, NelY, rho1, vs1, rho2, vs2, dy_deta, dx_dxi, wgll, iglob, side)`** - Creates absorbing boundary condition matrices for domain edges ('L', 'R', 'T', 'B')

## Core Simulation Engine

### `src/main.jl`
- **`main(P)`** - Main simulation loop
  - Initializes kinematic fields (displacement, velocity, acceleration)
  - Implements adaptive solver switching between quasi-static (isolver=1) and dynamic (isolver=2) regimes
  - Handles fault boundary conditions with rate-state friction
  - Manages output at specified time intervals
  - Controls earthquake event detection and recording

## Fault Mechanics and Solver Functions

### `src/NRsearch.jl` - Newton-Raphson Search
- **`FBC!(IDstate, P, NFBC, FltNglob, psi1, Vf1, tau1, psi2, Vf2, tau2, psi, Vf, FltVfree, dt)`** - Fault boundary condition solver using two-step Newton-Raphson iterations
- **`NRsearch!(fo, Vo, cca, ccb, Seff, tau, tauo, psi, FltZ, FltVfree)`** - Core Newton-Raphson algorithm for fault slip velocity computation

### `src/otherFunctions.jl` - Rate-State Utilities
- **`slrFunc!(P, NFBC, FltNglob, psi, psi1, Vf, Vf1, IDstate, tau1, dt)`** - Computes slip rates for quasi-static regime
- **`IDS!(xLf, Vo, psi, dt, Vf, cnd, IDstate)`** - Integrates state variable evolution (aging law)
- **`IDS2!(xLf, Vo, psi, psi1, dt, Vf, Vf1, IDstate)`** - Second-order state variable integration
- **`exp1(x)` / `log1(x)`** - Optimized exponential/logarithm functions

### `src/dtevol.jl` - Time Step Control
- **`dtevol!(dt, dtmin, XiLf, FaultNglob, NFBC, Vf, isolver)`** - Adaptive time step computation based on CFL condition and fault velocities
  - Maximum time step: 50 days
  - Time step increase factor: 1.2
  - CFL constraint: dt < XiLf/|Vf|

## Linear Algebra and Solvers

### `src/PCG.jl` - Preconditioned Conjugate Gradient
- **`PCG!(P, Nel, diagKnew, dnew, F, iFlt, FltNI, H, Ht, iglob, nglob, W, a_elem, Conn)`** - Preconditioned conjugate gradient solver for large sparse systems
- **`element_computation!(P, iglob, F_local, H, Ht, W, Nel)`** - Multi-threaded element-level computations

## Utilities and Analysis

### `src/FindNearestNode.jl`
- **`FindNearestNode(xin, yin, X, Y)`** - Finds nearest mesh nodes to specified coordinates for output locations

### `src/damageEvol.jl`
- **`damage_indx!(ThickX, ThickY, dxe, dye, NGLL, NelX, NelY, iglob)`** - Identifies damage indices for fault zone evolution

### `src/dump.jl`
- **`stiff_element(NGLL, NelX, NelY, nglob, iglob, dxe, dye)`** - Test function for element stiffness validation (development use only)

## Key Parameters and Setup

### Material Properties (in `par.jl`):
- **Host rock**: ρ₁ = 2670 kg/m³, vs₁ = 3464 m/s
- **Damage zone**: ρ₂ = 2670 kg/m³, vs₂ = 3464 m/s (adjustable)
- **Domain**: 48 km (depth) × 30 km (width)

### Fault Parameters:
- **Plate loading**: Vpl = 1×10⁻⁹ m/s
- **Reference friction**: f₀ = 0.6
- **Reference velocity**: V₀ = 1×10⁻⁶ m/s
- **Characteristic slip**: Dc = 8 mm
- **Rate-state parameters**: a, b values vary with depth

### Time Control:
- **Total simulation**: 150 years
- **CFL number**: 0.6
- **Velocity threshold**: 0.001 m/s (earthquake detection)
- **State evolution**: IDstate = 2 (aging law)

## Output Files

The simulation generates several output files in `data/$(simulation_name)/`:
- `stress.out` - Shear stress evolution
- `sliprate.out` - Fault slip rate time series
- `slip.out` - Cumulative slip
- `event_time.out` - Earthquake start/end times and hypocenters
- `coseismic_slip.out` - Slip during individual events
- `params.out` - Simulation parameters

## Usage

1. Edit `par.jl` to set desired parameters
2. Edit `src/initialConditions/defaultInitialConditions.jl` for custom initial conditions
3. Set resolution and output directory in `run.jl`
4. Run: `julia run.jl`
5. Analyze results using provided scripts in the analysis directory