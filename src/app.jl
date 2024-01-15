# main function that PackageCompiler will build against
function julia_main()::Cint
    try
        # read command line arguments
        zfile, knockoff_dir, N, hg_build, outfile, outdir, seed = 
            parse_commandline(true)

        println("\n\nWelcome to GhostKnockoffGWAS analysis!")
        println("You have specified the following options:")
        println("zfile           = ", abspath(zfile))
        println("knockoff_dir    = ", abspath(knockoff_dir))
        println("N (sample size) = ", N)
        println("hg_build        = ", hg_build)
        println("outdir          = ", outdir)
        println("outfile         = ", joinpath(outdir, outfile))
        println("seed            = ", seed)
        println("\n")

        # read Z scores
        t1 = @elapsed begin
            z, chr, pos, effect_allele, non_effect_allele = read_zscores(zfile)
        end

        # run ghost knockoff analysis
        t2 = @elapsed begin
            ghostbasil_parallel(knockoff_dir, z, chr, pos, effect_allele, 
                non_effect_allele, N, hg_build, outdir, outname=outfile, 
                seed=seed)
        end

        println("Done! Result saved to $(joinpath(outdir, outfile)). ")
        println("Overall runtime = $(t1 + t2) seconds, with ")
        println("   $t1 seconds spent on reading the Z score file")
        println("   $t2 seconds spent on doing the analysis")
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
                   ".gz compressed. The first row must be a header line that " *
                   "contains at least CHR, POS, REF, ALT, and Z (other columns " * 
                   "will be ignored). Each row should be a SNP. CHR is the " *
                   "chromosome column and must be integer valued (e.g. chr22, " *
                   ", sex chromosomes, and missing values are NOT valid). POS is " * 
                   "the SNP (aligned to HG19 or HG38) and cannot be missing. REF " *
                   "the position of and ALT are the reference and alternate " * 
                   "alleles, which will be treated as the non-effective and effect " * 
                   "alleles, respectively, and also cannot be missing. Finally, Z " * 
                   "is the Z-score column. Missing Z scores can be specified as " * 
                   "NaN or as an empty cell."
            required = true
            arg_type = String
        "--knockoff-dir"
            help = "Path to the directory storing pre-processed knockoff files"
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
        "--seed"
            help = "Sets the random seed"
            arg_type = Int
            default = 2023
    end

    # This is for code pre-compilation, enabling fast printing of "help statement".
    # It runs the parser once when the module is loaded
    # so that Julia can compile everything including functions that are defined in
    # other modules. Of course, we need to make sure that the arguments are valid or
    # the module will fail to load
    if !parseargs
        _useless = parse_args(
            ["--zfile","testdir","--knockoff-dir","testdir2",
            "--N","1","--genome-build","38","--out","testdir3","--seed","2024"], s
        )
        _useless = parse_args(["--help"], s)
        return nothing
    end

    parsed_args = parse_args(s)
    zfile = parsed_args["zfile"]
    knockoff_dir = parsed_args["knockoff-dir"]
    N = parsed_args["N"]
    hg_build = parsed_args["genome-build"]
    out = parsed_args["out"]
    outfile = basename(out)
    outdir = abspath(dirname(out))
    seed = parsed_args["seed"]

    return zfile, knockoff_dir, N, hg_build, outfile, outdir, seed
end
