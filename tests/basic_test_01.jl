#######################################
# Basic testing to visualize results
# #####################################

include("$(@__DIR__)/../analyze_results.jl")

VfmaxPlot(Vfmax, t, yr2sec);
cumSlipPlot(delfsec[1:4:end,:], delfyr[1:4:end, :], FltX);
