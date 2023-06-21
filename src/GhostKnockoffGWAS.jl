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

export ghostbasil

include("ghostbasil.jl")
include("utilities.jl")

end # module GhostKnockoffGWAS
