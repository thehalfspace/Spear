#######################################
# Basic testing to visualize results
# #####################################

include("$(@__DIR__)/../analyze_results.jl")

VfmaxPlot(Vfmax, t, yr2sec);
cumSlipPlot(delfsec, delfyr, FltX);
