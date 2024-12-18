module AvianRangeShapes
Pkg.add(url="https://github.com/mkborregaard/SpreadingDye.jl")

using ArchGDAL,Rasters, DataFrames,Plots, SpreadingDye, NearestNeighbors, Images, JLD2, SkipNan
include("functions.jl")

export AvianRangeShapes
end
