## User Guide for running simulations 


### Example 1: Elastic halfspace with a narrow fault zone

1. Open `run.jl` and edit the output file name and the resolution of the problem.
2. Open `par.jl` and set up the parameters. Choose the fault zone properties (stiffness, density, dimensions) in the section `MEDIUM PROPERTIES`.
3. Run `run.jl` from the terminal or IDE of your choice.
4. The output files are saved in `data/$(simulation_name)`.
5. You can look at the output file from `analyze_results.jl`.


### Example 2: Same setup as above, but with time-dependent healing of fault zone stiffness during the interseismic period.
1. Open `run.jl` and edit the output file name and the resolution of the problem. Use the function `src/main_damage_healing.jl` instead of `src/main.jl` by uncommenting this line. 
2. Open `par.jl` and set up the parameters. Choose the fault zone properties (stiffness, density, dimensions) in the section `MEDIUM PROPERTIES`.
3. The healing parameters (healing rate and healing amount) are set in `src/main_damage_healing.jl` in the function healing (lines 16-25).
4. Run `run.jl` from the terminal or IDE of your choice.
5. The output files are saved in `data/$(simulation_name)`.
6. You can look at the output file from `analyze_results.jl`.
