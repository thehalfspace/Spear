# SPEAR Postprocessing Scripts - Function Reference

## Overview

The `post/` directory contains Julia scripts for analyzing and visualizing earthquake simulation results from SPEAR. These scripts process the output data files to extract earthquake event details, compute seismic parameters, and generate publication-quality plots.

## Main Postprocessing Scripts

### `post/event_details.jl` - Event Analysis and Calculations

**Description**: Core analysis functions for extracting earthquake event characteristics from simulation output data.

#### Main Functions

**`get_index(seismic_stress, taubefore)`**
- **Purpose**: Finds the index of rupture initiation for each earthquake event
- **Input**:
  - `seismic_stress` - stress time series during seismic events
  - `taubefore` - stress before each event
- **Method**: Uses L2 norm to match stress patterns between event start and seismic time series
- **Returns**: Array of starting indices for each earthquake event

**`Coslip(S, Slip, SlipVel, Stress, time_)`**
- **Purpose**: Computes coseismic slip, stress drops, and event timing for all earthquake events
- **Input**:
  - `S` - simulation parameters structure
  - `Slip` - cumulative slip time series
  - `SlipVel` - slip velocity time series
  - `Stress` - shear stress time series
  - `time_` - simulation time vector
- **Parameters**:
  - Event threshold: `Vthres = 0.001` m/s (earthquake detection)
  - Event start: max slip rate > 1.01 × threshold
  - Event end: max slip rate < 0.99 × threshold
- **Returns**:
  - `delfafter` - coseismic slip for each event
  - `stressdrops` - stress drop for each event
  - `tStart, tEnd` - event start and end times
  - `vhypo, hypo` - hypocenter velocity and location

**`moment_magnitude_new(mu, FltX, delfafter, stressdrops, time_)`**
- **Purpose**: Calculates moment magnitude and rupture characteristics for each earthquake
- **Input**:
  - `mu` - shear modulus
  - `FltX` - fault depth coordinates
  - `delfafter` - coseismic slip from `Coslip()`
  - `stressdrops` - stress drops from `Coslip()`
  - `time_` - time vector
- **Method**:
  - Slip threshold: 1% of maximum slip per event
  - Assumes square rupture area (depth dimension = along-strike dimension)
  - Seismic moment: M₀ = μ × Area × Depth
  - Moment magnitude: Mw = (2/3)×log₁₀(M₀×10⁷) - 10.7
- **Returns**:
  - `Mw` - moment magnitude
  - `del_sigma` - average stress drop
  - `fault_slip` - average coseismic slip
  - `rupture_len` - rupture length along depth

### `post/plotting_script.jl` - Visualization Functions

**Description**: Comprehensive plotting functions for creating publication-quality figures from simulation results.

#### Plot Configuration

**`plot_params()`**
- **Purpose**: Sets default matplotlib parameters for consistent publication-style plots
- **Settings**:
  - Font: STIXGeneral (scientific publication standard)
  - Fontsize: 15-17 pt for labels, 13 pt for legends
  - Line width: 2.0, axis width: 1.5
  - Tick directions: inward
  - Auto layout enabled

#### Main Plotting Functions

**`slipPlot(delfafter2, rupture_len, FltX, Mw, tStart)`**
- **Purpose**: Creates horizontal bar plots of coseismic slip vs time at multiple depths
- **Features**:
  - 4-panel subplot showing slip at 60m, 4km, 6km, and 8km depths
  - Color-coded by moment magnitude (inferno_r colormap)
  - Filters events with Mw > 2.8
  - Time axis inverted (older events at top)
- **Output**: `coseismic_slip.png` (300 DPI)

**`eqCyclePlot(sliprate, FltX)`**
- **Purpose**: Creates 2D heatmap of slip rate evolution with depth and time
- **Features**:
  - Logarithmic color scale (1e-9 to 1e0 m/s)
  - Shows earthquake cycles from 16 km depth to surface
  - Interpolated color mapping using bicubic interpolation
  - Inferno colormap for slip rate visualization
- **Output**: `interpolated_sliprate.png` (300 DPI)

**`VfmaxPlot(Vfmax, t, yr2sec)`**
- **Purpose**: Plots maximum fault slip rate over simulation time
- **Features**:
  - Logarithmic y-axis for slip rate
  - Time in years on x-axis
  - Shows transition between interseismic and coseismic periods
- **Output**: `Vfmax01.png` (300 DPI)

**`Vfmaxcomp(Vfmax1, t1, Vfmax2, t2, yr2sec)`**
- **Purpose**: Comparison plot of maximum slip rates from two different simulations
- **Features**:
  - Overlay plots with labels ("Thakur", "Abdelmeguid")
  - Useful for validation and parameter studies
- **Output**: `Vfmax01.png` (300 DPI)

**`alphaaPlot(alphaa, t, yr2sec)`**
- **Purpose**: Plots shear modulus contrast evolution over time
- **Application**: For simulations with time-dependent material properties
- **Output**: `alpha_01.png` (300 DPI)

**`cumSlipPlot(delfsec, delfyr, FltX)`**
- **Purpose**: Shows cumulative slip profiles with depth
- **Features**:
  - Two datasets: annual slip (blue) and seismic slip (brown)
  - Depth range: 0-24 km
  - Inverted y-axis (depth increases downward)
- **Output**: `cumulative_slip.png` (300 DPI)

**`icsPlot(a_b, Seff, tauo, FltX)`**
- **Purpose**: Plots initial conditions - friction parameters and stress with depth
- **Features**:
  - Left axis: Normal and shear stress (MPa)
  - Right axis: Rate-state friction parameter (a-b)
  - Twin x-axes with different colors
  - Depth range: 0-80 km
- **Output**: `ics_02.png` (300 DPI)

### `post/rough_script.jl` - Development and Testing Functions

**Description**: Experimental scripts for testing new analysis approaches and detailed event visualization.

#### Main Functions

**`event_indx(tStart, tEnd, time_)`**
- **Purpose**: Finds array indices corresponding to earthquake start and end times
- **Input**: Event times and simulation time vector
- **Returns**: Start and end indices for each event
- **Use**: Facilitates detailed analysis of individual earthquake events

**`test1(S, O, evno)`**
- **Purpose**: Creates detailed slip rate evolution plot for a specific earthquake event
- **Features**:
  - High temporal resolution (0.1 second intervals)
  - Slip rate vs depth profile during single event
  - Multiple time snapshots overlaid
  - Depth range: 0-24 km
- **Input**:
  - `S` - simulation parameters
  - `O` - output data structure
  - `evno` - event number to analyze
- **Output**: `slipvel.pdf` (300 DPI)

## Data Processing Workflow

### Typical Analysis Sequence:

1. **Load simulation output files:**
   ```julia
   # From data/simulation_name/
   stress = readdlm("stress.out")
   sliprate = readdlm("sliprate.out")
   slip = readdlm("slip.out")
   delfsec = readdlm("delfsec.out")
   delfyr = readdlm("delfyr.out")
   ```

2. **Extract event characteristics:**
   ```julia
   include("post/event_details.jl")
   delfafter, stressdrops, tStart, tEnd, vhypo, hypo = Coslip(S, slip, sliprate, stress, time_)
   Mw, del_sigma, fault_slip, rupture_len = moment_magnitude_new(mu, FltX, delfafter, stressdrops, time_)
   ```

3. **Generate visualizations:**
   ```julia
   include("post/plotting_script.jl")
   path = "plots/simulation_name/"  # Set output directory
   plot_params()  # Set plot styling

   # Create plots
   slipPlot(delfafter, rupture_len, FltX, Mw, tStart)
   VfmaxPlot(Vfmax, time_, yr2sec)
   cumSlipPlot(delfsec, delfyr, FltX)
   icsPlot(a_b, Seff, tauo, FltX)
   ```

## Output Files and Analysis

### Key Metrics Computed:
- **Earthquake catalogs**: Mw, rupture length, stress drop, recurrence intervals
- **Slip distributions**: Coseismic vs interseismic slip partitioning
- **Rupture dynamics**: Hypocenter locations, rupture velocities
- **Stress evolution**: Pre- and post-event stress states

### Visualization Outputs:
- **Event catalogs**: Slip vs time plots colored by magnitude
- **Spatiotemporal evolution**: 2D heatmaps of slip rate
- **Time series**: Maximum slip rate evolution
- **Depth profiles**: Cumulative slip and initial conditions
- **Individual events**: Detailed rupture progression

## Dependencies

**Required Packages:**
- `PyPlot` - Python matplotlib interface
- `StatsBase` - Statistical functions
- `LaTeXStrings` - LaTeX formatting for labels
- `PyCall` - Python integration
- `LinearAlgebra` - Matrix operations

**Python Dependencies:**
- `matplotlib` - Plotting backend
- Scientific color maps (inferno, etc.)

## Usage Notes

1. **Set output path**: Define `path` variable before calling plotting functions
2. **Color normalization**: Magnitude-based coloring automatically scales to data range
3. **Resolution**: All plots save at 300 DPI for publication quality
4. **Filtering**: Event analysis typically filters by magnitude thresholds
5. **Memory**: Large datasets may require chunked processing for memory efficiency

This postprocessing suite provides comprehensive tools for earthquake cycle analysis, from basic event detection to detailed rupture characterization and publication-ready visualization.