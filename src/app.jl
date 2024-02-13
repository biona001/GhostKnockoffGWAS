# main function that PackageCompiler will build against
function julia_main()::Cint
    try
        # read command line arguments
        zfile, LD_files, N, hg_build, outfile, outdir, chr_col, pos_col, 
            ref_col, alt_col, z_col, seed, verbose,
            random_shuffle, skip_shrinkage_check = parse_commandline(true)

        println("\n\nWelcome to GhostKnockoffGWAS analysis!")
        println("You have specified the following options:")
        println("zfile           = ", abspath(zfile))
        println("LD_files        = ", abspath(LD_files))
        println("N (sample size) = ", N)
        println("hg_build        = ", hg_build)
        println("outdir          = ", outdir)
        println("outfile         = ", joinpath(outdir, outfile))
        !isnothing(chr_col) && 
        println("chr_col         = ", chr_col)
        !isnothing(pos_col) && 
        println("pos_col         = ", pos_col)
        !isnothing(ref_col) && 
        println("ref_col         = ", ref_col)
        !isnothing(alt_col) && 
        println("alt_col         = ", alt_col)
        !isnothing(z_col) && 
        println("z_col           = ", z_col)
        println("seed            = ", seed)
        println("verbose         = ", verbose)
        println("random_shuffle  = ", random_shuffle)
        println("skip_shrinkage_check = ", skip_shrinkage_check)
        println("\n")

        # read Z scores
        t1 = @elapsed begin
            z, chr, pos, effect_allele, non_effect_allele = read_zscores(
                zfile, chr_col=chr_col, pos_col=pos_col, ref_col=ref_col, 
                alt_col=alt_col, z_col=z_col
            )
        end

        # run ghost knockoff analysis
        t2 = @elapsed begin
            ghostknockoffgwas(LD_files, z, chr, pos, effect_allele, 
                non_effect_allele, N, hg_build, outdir, outname=outfile, 
                seed=seed, verbose=verbose, 
                skip_shrinkage_check=skip_shrinkage_check,
                random_shuffle=random_shuffle)
        end

        if verbose
            println("Done! Result saved to $(joinpath(outdir, outfile)). ")
            println("Overall runtime = $(t1 + t2) seconds, with ")
            println("   $t1 seconds spent on reading the Z score file")
            println("   $t2 seconds spent on doing the analysis")
        end
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

function parse_commandline(parseargs::Bool)
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--zfile"
            help = "Tab or comma separated summary Z-score file, which can be " * 
                   ".gz compressed. By default, we assume the first row is a header line that " *
                   "contains at least CHR, POS, REF, ALT, and Z (other columns " * 
                   "will be ignored). Each row should be a SNP. CHR is the " *
                   "chromosome column and must be integer valued (e.g. chr22, " *
                   ", sex chromosomes, and missing values are NOT valid). POS is " * 
                   "the SNP basepair position (aligned to HG19 or HG38) and cannot be missing. REF " *
                   "the position of and ALT are the reference and alternate " * 
                   "alleles, which will be treated as the non-effective and effect " * 
                   "alleles, respectively, and also cannot be missing. Finally, Z " * 
                   "is the Z-score column. Missing Z scores can be specified as " * 
                   "NaN or as an empty cell."
            required = true
            arg_type = String
        "--LD-files"
            help = "Path to the directory storing pre-processed LD and knockoff files"
            required = true
            arg_type = String
        "--N"
            help = "Sample size for target (original) study"
            required = true
            arg_type = Int
        "--genome-build"
            help = "Specifies the human genome build for the target (original) " * 
                   "study. Must be 19 (hg19) or 38 (hg38)."
            arg_type = Int
            required = true
        "--out"
            help = "Output file prefix (without extensions)"
            required = true
            arg_type = String
        "--CHR"
            help = "The column in `zfile` that will be read as chromosome number " * 
                   "(if not specified, we read the column with header `CHR` in " * 
                   "`zfile`)"
            default = nothing
        "--POS"
            help = "The column in `zfile` that will be read as position (basepair)." * 
                   "(if not specified, we read the column with header `POS` in " * 
                   "`zfile`)"
            default = nothing
        "--REF"
            help = "The column in `zfile` that will be read as reference (non-" * 
                   "effective) allele. (if not specified, we read the column " * 
                   "with header `REF` in `zfile`)"
            default = nothing
        "--ALT"
            help = "The column in `zfile` that will be read as ALT (non-effective)" * 
                   "allele. (if not specified, we read the column with header " * 
                   "`ALT in `zfile`)"
            default = nothing
        "--Z"
            help = "The column in `zfile` that will be read as the Z-score" *
                   "(if not specified, we read the column with header `Z` " * 
                   "in `zfile`)"
            default = nothing
        "--seed"
            help = "Sets the random seed"
            arg_type = Int
            default = 2023
        "--verbose"
            help = "Whether to print intermediate messages"
            arg_type = Bool
            default = true
        "--random-shuffle"
            help = "Whether to randomly permute the order of Z-scores and " * 
                   "their knockoffs to adjust for potential ordering bias."
            arg_type = Bool
            default = true
        "--skip-shrinkage-check"
            help = "Whether to allow Knockoff analysis to proceed even with " * 
                   "large (>0.25) LD shrinkages. Only use this option if you " *
                   "know what you are doing. "
            arg_type = Bool
            default = false
    end

    # This is for code pre-compilation, enabling fast printing of "help statement".
    # It runs the parser once when the module is loaded
    # so that Julia can compile everything including functions that are defined in
    # other modules. Of course, we need to make sure that the arguments are valid or
    # the module will fail to load
    if !parseargs
        _useless = parse_args(
            ["--zfile","testdir","--LD-files","testdir2",
            "--N","1","--genome-build","38","--out","testdir3",
            "--CHR","1","--POS","2","--EA","3","--NEA","4","--Z","5",
            "--seed","2024","--verbose","true"], s
        )
        _useless = parse_args(["--help"], s)
        return nothing
    end

    parsed_args = parse_args(s)
    zfile = parsed_args["zfile"]
    LD_files = parsed_args["LD-files"]
    N = parsed_args["N"]
    hg_build = parsed_args["genome-build"]
    out = parsed_args["out"]
    outfile = basename(out)
    outdir = abspath(dirname(out))
    chr_col = parsed_args["CHR"]
    pos_col = parsed_args["POS"]
    ref_col = parsed_args["REF"]
    alt_col = parsed_args["ALT"]
    z_col = parsed_args["Z"]
    seed = parsed_args["seed"]
    verbose = parsed_args["verbose"]
    random_shuffle = parsed_args["random-shuffle"]
    skip_shrinkage_check = parsed_args["skip-shrinkage-check"]

    return zfile, LD_files, N, hg_build, outfile, outdir, 
        chr_col, pos_col, ref_col, alt_col, z_col, seed, verbose, 
        random_shuffle, skip_shrinkage_check
end
