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

export ghostbasil, 
    ghostbasil_parallel_by_block,
    graphical_group_S,
    # utilities
    zscore2pval,
    zscore,
    zscore2pval,
    pval2zscore,
    pval

include("ghostbasil.jl")
include("utilities.jl")
include("solve_blocks.jl")

end # module GhostKnockoffGWAS
