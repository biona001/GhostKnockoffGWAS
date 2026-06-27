module CITLassoJuliaC

using CodecZlib
using Distributions
using Ghostbasil
using HDF5
using LinearAlgebra
using Optim
using Random
using Statistics
using StatsBase

include("utilities.jl")
include("ghostbasil_parallel.jl")

end

using .CITLassoJuliaC

struct JuliaCCITLassoOptions
    zfile::String
    ld_files::String
    n::Int
    genome_build::Int
    outfile::String
    outdir::String
    chr_col::String
    pos_col::String
    ref_col::String
    alt_col::String
    z_col::String
    seed::Int
    verbose::Bool
    random_shuffle::Bool
    skip_shrinkage_check::Bool
end

function print_usage()::Nothing
    println("usage: cit-lasso --zfile FILE --LD-files DIR --N N --genome-build 19|38 --out OUT [options]")
    println()
    println("required:")
    println("  --zfile FILE          Tab- or comma-separated Z-score file, optionally .gz")
    println("  --LD-files DIR        Directory containing preprocessed LD and knockoff files")
    println("  --N N                 Target-study sample size")
    println("  --genome-build 19|38  Human genome build used by the Z-score positions")
    println("  --out OUT             Output file prefix")
    println()
    println("options:")
    println("  --CHR COL --POS COL --REF COL --ALT COL --Z COL")
    println("  --seed INT")
    println("  --verbose true|false")
    println("  --random-shuffle true|false")
    println("  --skip-shrinkage-check true|false")
    return nothing
end

function parse_bool(s::String)::Bool
    v = lowercase(s)
    if v == "true" || v == "1" || v == "yes"
        return true
    elseif v == "false" || v == "0" || v == "no"
        return false
    end
    error("Expected a boolean value, got $s")
end

function next_arg(args::Vector{String}, i::Int, option::String)::String
    i < length(args) || error("Missing value for $option")
    return args[i + 1]
end

function parse_options(args::Vector{String})::JuliaCCITLassoOptions
    zfile = ""
    ld_files = ""
    n = 0
    genome_build = 0
    out = ""
    chr_col = ""
    pos_col = ""
    ref_col = ""
    alt_col = ""
    z_col = ""
    seed = 2023
    verbose = true
    random_shuffle = false
    skip_shrinkage_check = false

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--zfile"
            zfile = next_arg(args, i, arg)
            i += 2
        elseif arg == "--LD-files"
            ld_files = next_arg(args, i, arg)
            i += 2
        elseif arg == "--N"
            n = parse(Int, next_arg(args, i, arg))
            i += 2
        elseif arg == "--genome-build"
            genome_build = parse(Int, next_arg(args, i, arg))
            i += 2
        elseif arg == "--out"
            out = next_arg(args, i, arg)
            i += 2
        elseif arg == "--CHR"
            chr_col = next_arg(args, i, arg)
            i += 2
        elseif arg == "--POS"
            pos_col = next_arg(args, i, arg)
            i += 2
        elseif arg == "--REF"
            ref_col = next_arg(args, i, arg)
            i += 2
        elseif arg == "--ALT"
            alt_col = next_arg(args, i, arg)
            i += 2
        elseif arg == "--Z"
            z_col = next_arg(args, i, arg)
            i += 2
        elseif arg == "--seed"
            seed = parse(Int, next_arg(args, i, arg))
            i += 2
        elseif arg == "--verbose"
            verbose = parse_bool(next_arg(args, i, arg))
            i += 2
        elseif arg == "--random-shuffle"
            random_shuffle = parse_bool(next_arg(args, i, arg))
            i += 2
        elseif arg == "--skip-shrinkage-check"
            skip_shrinkage_check = parse_bool(next_arg(args, i, arg))
            i += 2
        elseif arg == "--help" || arg == "-h"
            print_usage()
            return JuliaCCITLassoOptions("", "", 0, 0, "", "", "", "", "", "", "", 2023, false, false, false)
        else
            error("Unknown option: $arg")
        end
    end

    isempty(zfile) && error("Missing required option: --zfile")
    isempty(ld_files) && error("Missing required option: --LD-files")
    n > 0 || error("Missing or invalid required option: --N")
    genome_build == 19 || genome_build == 38 || error("Missing or invalid required option: --genome-build")
    isempty(out) && error("Missing required option: --out")

    return JuliaCCITLassoOptions(
        zfile,
        ld_files,
        n,
        genome_build,
        basename(out),
        abspath(dirname(out)),
        chr_col,
        pos_col,
        ref_col,
        alt_col,
        z_col,
        seed,
        verbose,
        random_shuffle,
        skip_shrinkage_check,
    )
end

function run_cit_lasso(opts::JuliaCCITLassoOptions)::Cint
    isempty(opts.zfile) && return 0

    if opts.verbose
        println()
        println("Welcome to CIT-Lasso analysis!")
        println("zfile           = ", abspath(opts.zfile))
        println("LD_files        = ", abspath(opts.ld_files))
        println("N (sample size) = ", opts.n)
        println("hg_build        = ", opts.genome_build)
        println("outdir          = ", opts.outdir)
        println("outfile         = ", joinpath(opts.outdir, opts.outfile))
        isempty(opts.chr_col) || println("chr_col         = ", opts.chr_col)
        isempty(opts.pos_col) || println("pos_col         = ", opts.pos_col)
        isempty(opts.ref_col) || println("ref_col         = ", opts.ref_col)
        isempty(opts.alt_col) || println("alt_col         = ", opts.alt_col)
        isempty(opts.z_col) || println("z_col           = ", opts.z_col)
        println("seed            = ", opts.seed)
        println("verbose         = ", opts.verbose)
        println("random_shuffle  = ", opts.random_shuffle)
        println("skip_shrinkage_check = ", opts.skip_shrinkage_check)
        println()
    end

    z, chr, pos, effect_allele, non_effect_allele = CITLassoJuliaC.read_zscores(
        opts.zfile;
        chr_col=opts.chr_col,
        pos_col=opts.pos_col,
        ref_col=opts.ref_col,
        alt_col=opts.alt_col,
        z_col=opts.z_col,
    )

    CITLassoJuliaC.ghostknockoffgwas(
        opts.ld_files,
        z,
        chr,
        pos,
        effect_allele,
        non_effect_allele,
        opts.n,
        opts.genome_build,
        opts.outdir;
        outname=opts.outfile,
        seed=opts.seed,
        verbose=opts.verbose,
        skip_shrinkage_check=opts.skip_shrinkage_check,
        random_shuffle=opts.random_shuffle,
    )
    return 0
end

function @main(args::Vector{String})::Cint
    try
        return run_cit_lasso(parse_options(args))
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
end
