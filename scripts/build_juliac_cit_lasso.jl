#!/usr/bin/env julia

using Pkg

function arg_value(args::Vector{String}, name::String, default::String)::String
    idx = findfirst(==(name), args)
    if isnothing(idx)
        return default
    end
    idx < length(args) || error("Missing value for $name")
    return args[idx + 1]
end

function repo_root()::String
    return normpath(@__DIR__, "..")
end

function setup_build_env!(env_dir::String)::Nothing
    Pkg.activate(env_dir)
    Pkg.add("JuliaC")
    Pkg.add(["CodecZlib", "Distributions", "HDF5", "Optim", "StatsBase"])
    Pkg.add(url="https://github.com/biona001/ghostbasil_jll.jl")
    Pkg.add(url="https://github.com/biona001/Ghostbasil.jl")
    Pkg.resolve()
    Pkg.instantiate()
    Pkg.precompile()
    return nothing
end

function build_cit_lasso(; target::String, trim::String, build_dir::String)::String
    root = repo_root()
    env_dir = mktempdir(; prefix="citlasso-juliac-env-")
    setup_build_env!(env_dir)

    bundle_dir = joinpath(build_dir, "cit-lasso-$target")
    entrypoint = joinpath(root, "src", "juliac_cit_lasso.jl")
    args = [
        "--output-exe", "cit-lasso",
        "--bundle", bundle_dir,
        "--trim=$trim",
        "--project", env_dir,
        entrypoint,
    ]
    trim == "no" || insert!(args, length(args), "--experimental")

    println("Building JuliaC AOT cit-lasso for $target with trim=$trim")
    run(`$(Base.julia_cmd()) --project=$env_dir -e 'using JuliaC; JuliaC.main(ARGS)' -- $args`)

    exe = joinpath(bundle_dir, "bin", Sys.iswindows() ? "cit-lasso.exe" : "cit-lasso")
    isfile(exe) || error("JuliaC did not produce $exe")
    run(`$exe --help`)
    return bundle_dir
end

function main(args::Vector{String})::Nothing
    target = arg_value(args, "--target", string(Sys.KERNEL, "-", Sys.MACHINE))
    trim = arg_value(args, "--trim", get(ENV, "CITLASSO_JULIAC_TRIM", "no"))
    build_dir = abspath(arg_value(args, "--build-dir", joinpath(repo_root(), "build")))
    mkpath(build_dir)

    bundle_dir = build_cit_lasso(target=target, trim=trim, build_dir=build_dir)
    println("Built $bundle_dir")
    run(`du -sh $bundle_dir`)
    return nothing
end

main(ARGS)
