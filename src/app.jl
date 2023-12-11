# main function that PackageCompiler will build against
function julia_main()::Cint
    try
        # read args
        zfile, knockoff_dir, Neffect, hg_build, outfile, outdir, m, seed = 
            parse_commandline()

        println("Running GhostKnockoffGWAS analysis with the following options:")
        println("zfile = $zfile")
        println("knockoff_dir = $knockoff_dir")
        println("N (effective sample size) = $Neffect")
        println("hg_build = $hg_build")
        println("outfile = $outfile")
        println("m (number of knockoffs for stability) = $m")
        println("seed = $seed")

        # run ghost knockoff analysis
        z, chr, pos, effect_allele, non_effect_allele = read_zscores(zfile)
        ghostbasil_parallel(knockoff_dir, z, chr, pos, effect_allele, 
            non_effect_allele, Neffect, hg_build, outdir, outname=outfile, 
            m=m, seed=seed)

        println("Done! Result saved to $(joinpath(outdir, outname))")
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "zfile"
            help = "Tab or comma separated summary Z-score file. The first row " *
                   "must be a header line that contains CHR/POS/REF/ALT fields, " * 
                   "as well as something that indicates Z-score information. " *
                   "CHR must be integer valued (e.g. chr21 is not correct). " * 
                   "If a column with `Z` is available, it will be used as Z " * 
                   "scores. Otherwise, we need both `pvalue` and `beta` to be " *
                   "available, or both `OR` (odds ratio) and `SE` (standard " * 
                   "error) to be available. In these cases, we will convert " * 
                   "them to Z scores."
            required = true
            arg_type = String
        "knockoff_dir"
            help = "Path to the directory storing pre-processed knockoff files"
            required = true
            arg_type = String
        "Neffect"
            help = "Effective sample size for original data"
            required = true
            arg_type = Int
        "out"
            help = "Output file prefix (without extensions)"
            required = true
            arg_type = String
        "--genome-build"
            help = "Specifies the human genome build for the summary file. " * 
                   "Must be 19 (hg19) or 38 (hg38)."
            arg_type = Int
            default = 38
        "--m"
            help = "Number of knockoffs to generate"
            arg_type = Int
            default = 5
        "--seed"
            help = "Sets the random seed"
            arg_type = Int
            default = 2023
    end
    parsed_args = parse_args(s)

    zfile = parsed_args["zfile"]
    knockoff_dir = parsed_args["knockoff_dir"]
    Neffect = parsed_args["Neffect"]
    hg_build = parsed_args["genome-build"]
    out = parsed_args["out"]
    outfile = basename(out)
    outdir = dirname(outfile)
    m = parsed_args["m"]
    seed = parsed_args["seed"]

    return zfile, knockoff_dir, Neffect, hg_build, outfile, outdir, m, seed
end
