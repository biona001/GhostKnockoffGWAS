module GhostKnockoffGWAS

using Knockoffs
using JLD2
using HDF5
using CSV
using DataFrames
using LinearAlgebra
using Random
using Distributions
using DelimitedFiles
using RCall
R"library(ghostbasil)" # todo: write C++ wrapper to avoid R
R"library(liftOver)"

export ghostbasil, graphical_group_S

include("ghostbasil.jl")
include("utilities.jl")
include("solve_blocks.jl")

end # module GhostKnockoffGWAS
