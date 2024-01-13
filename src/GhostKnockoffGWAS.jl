module GhostKnockoffGWAS

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
using ArgParse
using Statistics

export ghostbasil, 
    ghostbasil_parallel,
    # utilities
    zscore2pval,
    zscore,
    zscore2pval,
    pval2zscore,
    pval,
    read_zscores

# include("ghostbasil.jl")
include("ghostbasil_parallel.jl")
include("utilities.jl")
include("app.jl")

end # module GhostKnockoffGWAS
