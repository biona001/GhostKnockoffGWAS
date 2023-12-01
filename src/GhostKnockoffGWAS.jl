module GhostKnockoffGWAS

using Knockoffs
using Ghostbasil
using JLD2
using HDF5
using CSV
using DataFrames
using LinearAlgebra
using Random
using Distributions
using DelimitedFiles
using Optim

export ghostbasil, 
    ghostbasil_parallel,
    # utilities
    zscore2pval,
    zscore,
    zscore2pval,
    pval2zscore,
    pval

# include("ghostbasil.jl")
include("ghostbasil_parallel.jl")
include("utilities.jl")

end # module GhostKnockoffGWAS
