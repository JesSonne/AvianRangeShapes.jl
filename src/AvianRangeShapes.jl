module AvianRangeShapes
using Pkg
using ArchGDAL,Rasters, DataFrames,Plots, SpreadingDye, NearestNeighbors, Images, JLD2, SkipNan
Pkg.instantiate()

include("functions.jl")


export null_models
end