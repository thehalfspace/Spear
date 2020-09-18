#################################
# Run the simulations from here
#################################

# 1. Go to par.jl and change as needed
# 2. Go to src/initialConditions/defaultInitialConditions and change as needed
# 3. Change the name of the simulation in this file
# 4. Run the simulation from terminal. (julia run.jl)
# 5. Plot results from the scripts folder

using Printf, LinearAlgebra, DelimitedFiles, SparseArrays,
    AlgebraicMultigrid, StaticArrays, IterativeSolvers, FEMSparse
using Base.Threads
#  BLAS.set_num_threads(1)

include("$(@__DIR__)/par.jl")	    #	Set Parameters

# Put the resolution for the simulation here: should be an integer
resolution = 2

# Output directory to save data
out_dir = "$(@__DIR__)/data/test_01/"
mkpath(out_dir)

P = setParameters(0e3,resolution)      # args = fault zone depth, resolution

include("$(@__DIR__)/src/dtevol.jl")
include("$(@__DIR__)/src/NRsearch.jl")
include("$(@__DIR__)/src/otherFunctions.jl")

include("$(@__DIR__)/src/main.jl")

simulation_time = @elapsed @time main(P)

println("\n")

@info("Simulation Complete!");
