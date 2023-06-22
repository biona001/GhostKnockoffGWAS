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

function __init__()
    # import R packages (todo: should we write C++ wrapper?)
    R"library(ghostbasil)"
    R"library(liftOver)"
end

end # module GhostKnockoffGWAS
