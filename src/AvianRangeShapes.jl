module AvianRangeShapes
Pkg.add(url="https://github.com/JesSonne/AvianRangeShapes.jl.git")
using ArchGDAL,Rasters, DataFrames,Plots, SpreadingDye, NearestNeighbors, Images, JLD2, SkipNan
Pkg.instantiate()

include("functions.jl")


export null_models
end