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
using Optim

export ghostbasil, 
    ghostbasil_parallel,
    # utilities
    zscore2pval,
    zscore,
    zscore2pval,
    pval2zscore,
    pval

include("ghostbasil.jl")
include("ghostbasil_parallel.jl")
include("utilities.jl")

function __init__()
    # import R packages (todo: should we write C++ wrapper?)
    R"library(ghostbasil)"
    R"library(liftOver)"
end

end # module GhostKnockoffGWAS
