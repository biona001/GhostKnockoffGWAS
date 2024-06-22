module GhostKnockoffGWAS

using Ghostbasil
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

# needed for solve_blocks()
using Knockoffs
using VCFTools
using VariantCallFormat
using ElasticArrays

export ghostknockoffgwas,
    read_zscores,
    solve_blocks

include("ghostbasil_parallel.jl")
include("utilities.jl")
include("app.jl")
include("make_hdf5.jl")

end # module GhostKnockoffGWAS
