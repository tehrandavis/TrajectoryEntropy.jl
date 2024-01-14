#module TrajectoryEntropy

using CSV, 
    DataFrames, 
    #Statistics, 
    #Optimization, 
    Plots, 
    Distances, 
    Unitful,
    NaNMath,
    Optim,
    #GLM,
    #EntropyHub,
    #NumericalIntegration,
    StatsBase,
    StatsPlots
    #Interpolations
    
    
include("coord_to_deg.jl");
include("uniquetol.jl");
include("bin_optimization.jl");
include("optimization_functions.jl");
include("entropy_decomposition.jl");
include("trajectory_entropy_analysis.jl");

#end
