include(joinpath(@__DIR__, "juliac_gwas_module.jl"))
using .GhostKnockoffGWASJuliaC

struct JuliaCGWASOptions
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

function juliac_usage()::Nothing
    println("usage: GhostKnockoffGWAS --zfile FILE --LD-files DIR --N N --genome-build 19|38 --out OUT [options]")
    println("options: --CHR COL --POS COL --REF COL --ALT COL --Z COL --seed INT --verbose BOOL --random-shuffle BOOL --skip-shrinkage-check BOOL")
    return nothing
end

function parse_juliac_bool(s::String)::Bool
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

function parse_juliac_options(args::Vector{String})::JuliaCGWASOptions
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
            verbose = parse_juliac_bool(next_arg(args, i, arg))
            i += 2
        elseif arg == "--random-shuffle"
            random_shuffle = parse_juliac_bool(next_arg(args, i, arg))
            i += 2
        elseif arg == "--skip-shrinkage-check"
            skip_shrinkage_check = parse_juliac_bool(next_arg(args, i, arg))
            i += 2
        elseif arg == "--help" || arg == "-h"
            juliac_usage()
            return JuliaCGWASOptions("", "", 0, 0, "", "", "", "", "", "", "", 2023, false, false, false)
        else
            error("Unknown option: $arg")
        end
    end

    isempty(zfile) && error("Missing required option: --zfile")
    isempty(ld_files) && error("Missing required option: --LD-files")
    n > 0 || error("Missing or invalid required option: --N")
    genome_build == 19 || genome_build == 38 || error("Missing or invalid required option: --genome-build")
    isempty(out) && error("Missing required option: --out")

    return JuliaCGWASOptions(
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

function run_juliac_gkgwas(opts::JuliaCGWASOptions)::Cint
    isempty(opts.zfile) && return 0

    if opts.verbose
        println()
        println("Welcome to GhostKnockoffGWAS analysis!")
        println("zfile           = ", abspath(opts.zfile))
        println("LD_files        = ", abspath(opts.ld_files))
        println("N (sample size) = ", opts.n)
        println("hg_build        = ", opts.genome_build)
        println("outdir          = ", opts.outdir)
        println("outfile         = ", joinpath(opts.outdir, opts.outfile))
    end

    z, chr, pos, effect_allele, non_effect_allele = GhostKnockoffGWASJuliaC.read_zscores(
        opts.zfile;
        chr_col=opts.chr_col,
        pos_col=opts.pos_col,
        ref_col=opts.ref_col,
        alt_col=opts.alt_col,
        z_col=opts.z_col,
    )

    GhostKnockoffGWASJuliaC.ghostknockoffgwas(
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
    return run_juliac_gkgwas(parse_juliac_options(args))
end
