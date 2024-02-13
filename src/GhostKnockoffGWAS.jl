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
using StatsBase

export ghostknockoffgwas,
    # utilities
    read_zscores

include("ghostbasil_parallel.jl")
include("utilities.jl")
include("app.jl")

end # module GhostKnockoffGWAS
