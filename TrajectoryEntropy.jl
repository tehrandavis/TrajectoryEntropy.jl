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
    
    
include("src/coord_to_deg.jl");
include("src/uniquetol.jl");
include("src/bin_optimization.jl");
include("src/optimization_functions.jl");
include("src/entropy_decomposition.jl");
include("src/trajectory_entropy_analysis.jl");


#end
