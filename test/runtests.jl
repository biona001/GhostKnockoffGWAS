using Test
using GhostKnockoffGWAS
using Ghostbasil
using LinearAlgebra
using DelimitedFiles
using CSV
using DataFrames
using Downloads
using VCFTools
using HDF5
using Distributions
using Random
using StatsBase
using SnpArrays

"""
    vcf_to_plink(plink_exe::String, vcffile::String, plinkprefix::String)

Converts a VCF file to binary PLINK file using `plink_exe` (which should point 
to the PLINK 1.9 executable). 

# A1/A2 allele in PLINK vs REF/ALT in VCF?
According to PLINK doc (https://www.cog-genomics.org/plink/1.9/data), adding
`--keep-allele-order` will convert REF alleles in VCF files to A2 alleles. Thus,
given our convention that `x_{ij} = 1` means there are 1 ALT allele, the resulting
PLINK file will have all `0 <--> 2` flipped if we import the original VCF and  
converted PLINK file into memory using VCFTools and SnpArrays. Fortunately, the
resulting correlation matrix is identical, so this is not an issue for knockoff
optimization. 
"""
function vcf_to_plink(plink_exe::String, vcffile::String, plinkprefix::String)
    run(`$plink_exe --vcf $vcffile --double-id --keep-allele-order --real-ref-alleles --make-bed --out $plinkprefix`)
end

@testset "utilities" begin
    # convert between z score and p-value
    p = 0.00695184371278
    beta = -0.010602
    z = 2.6991423324213346
    @test GhostKnockoffGWAS.zscore(p, beta) ≈ 2.6991423324213346
    @test GhostKnockoffGWAS.pval(z) ≈ 0.00695184371278

    # find_matching_indices
    a = [1, 2, 3, 2, 6]
    b = [2, 1, 2, 4, 3]
    c = GhostKnockoffGWAS.find_matching_indices(a, b)
    @test c == [[2], [1, 3], [5], [1, 3], Int[]]
    @test typeof(c) <: Vector{Vector{Int64}}
end

@testset "LD shrinkage" begin
    Σ = [1.0 -0.06603158126487506 -0.0833880805108651 -0.07423250848805772 -0.08520692260952517 -0.3018926160412726 -0.07228886018105062 -0.05163600016697706 -0.06831703956672884 0.1119291270044286; 
    -0.06603158126487506 1.0 -0.025855482863549748 -0.023593031867989712 -0.024785186697126534 0.17970487929732645 -0.0217395069925557 -0.015763968744450328 -0.0248758148261822 -0.01685099072623645; 
    -0.0833880805108651 -0.025855482863549748 1.0 -0.02953703088587616 -0.036090226242879775 -0.13347090146200266 -0.027912450562371478 -0.018630494841470058 -0.02981323926900486 -0.02095834431260501; 
    -0.07423250848805772 -0.023593031867989712 -0.02953703088587616 1.0 -0.03235265862653292 -0.12511846555868517 -0.026570109459531203 -0.0178496231359085 -0.027896176152173473 -0.01715817220843155; 
    -0.08520692260952517 -0.024785186697126534 -0.036090226242879775 -0.03235265862653292 1.0 0.2459334205306294 -0.030729029543925132 -0.023387666784556588 -0.03616407173066363 -0.02240494855276279; 
    -0.3018926160412726 0.17970487929732645 -0.13347090146200266 -0.12511846555868517 0.2459334205306294 1.0 0.19677801401533077 -0.04917243296590054 0.21936244494341658 -0.0015389401383875156; 
    -0.07228886018105062 -0.0217395069925557 -0.027912450562371478 -0.026570109459531203 -0.030729029543925132 0.19677801401533077 1.0 -0.020079411765085316 -0.027283609452267554 -0.019604609586933476; 
    -0.05163600016697706 -0.015763968744450328 -0.018630494841470058 -0.0178496231359085 -0.023387666784556588 -0.04917243296590054 -0.020079411765085316 1.0 -0.020149157087303162 -0.013591443637752606; 
    -0.06831703956672884 -0.0248758148261822 -0.02981323926900486 -0.027896176152173473 -0.03616407173066363 0.21936244494341658 -0.027283609452267554 -0.020149157087303162 1.0 -0.015041050757637644; 
    0.1119291270044286 -0.01685099072623645 -0.02095834431260501 -0.01715817220843155 -0.02240494855276279 -0.0015389401383875156 -0.019604609586933476 -0.013591443637752606 -0.015041050757637644 1.0]
    z = [-0.0830535851088237, -2.25068631502435, 1.49776720532444, 1.65978366398296, 1.06042628327983, -1.41526767195536, -0.522515227438188, -0.140465835528732, 0.484580620845268, -0.981466593482507]
    γ = GhostKnockoffGWAS.find_optimal_shrinkage(Σ, z)
    γtrue = 6.149509305168035e-15
    @test 0 ≤ γ ≤ 1
    @test isapprox(γ, γtrue, atol=1e-12)
end

# commented out since more extensive test of this is in Ghostbasil.jl
# @testset "ghostbasil C++ solver" begin
#     # basic test to make sure block_group_ghostbasil does not crash
#     p = 10 # number of features
#     m = 5  # number of knockoffs
#     x = randn(p, p)
#     y = randn(p, p)
#     Ci = x'*x
#     Si = y'*y
#     r = randn(6p)
#     lambda_path = [0.1, 0.05, 0.01, 0.001]
#     beta_i = block_group_ghostbasil(Ci, Si, r, lambda_path, m=m)
#     @test length(beta_i) == (m+1)*p

#     # Test correctness with a specific region from Pan-UKBB. 
#     # The true answer is obtained from running the ghostbasil R package on
#     # the given data. Note these data are based on LD summary statistics, 
#     # which can be freely distributed
#     datadir = normpath(Ghostbasil.datadir())
#     C = readdlm(joinpath(datadir, "C.txt"))
#     S = readdlm(joinpath(datadir, "S.txt"))
#     r = readdlm(joinpath(datadir, "r.txt")) |> vec
#     lambda_path = readdlm(joinpath(datadir, "lambda_path.txt")) |> vec
#     beta_i = block_group_ghostbasil(C, S, r, lambda_path)
#     @test count(!iszero, beta_i) == 14
#     non_zero_idx = findall(!iszero, beta_i)
#     @test non_zero_idx == [215, 340, 341, 431, 559, 560, 763, 878,
#         1085, 1528, 1770, 2726, 3157, 3387]
#     @test all(beta_i[non_zero_idx] .≈ 
#         [-0.0005451576807477974, -0.0002947665358369119, -0.014810108231841996, 
#         -0.00043238505905553924, -0.006856782960094195, 0.0016850133196618926, 
#         -0.00019274743491070077, -0.0003841629682186662, -0.0005758897376517483, 
#         -0.0016095272344540995, 0.0020146776266541312, -0.001666736964454875, 
#         0.0004269118551173737, 0.00010746291189108347]
#     )
# end

@testset "read_zscores" begin
    # no missing data
    filepath = tempname()
    data = """
    CHR\tPOS\tREF\tALT\tZ
    1\t758351\tA\tG\t1.05551882016081
    1\t779885\tC\tT\t2.12197306477
    1\t779987\tA\tG\t1.95791489337
    1\t782105\tC\tA\t1.91243829548
    """
    CSV.write(filepath, CSV.File(IOBuffer(data)))
    z, chr, pos, effect_allele, non_effect_allele = read_zscores(filepath)
    @test all(z .≈ [1.05551882016081, 2.12197306477, 1.95791489337, 1.91243829548])
    @test all(chr .== 1)
    @test all(pos .== [758351, 779885, 779987, 782105])
    @test all(effect_allele .== ["G", "T", "G", "A"])
    @test all(non_effect_allele .== ["A", "C", "A", "C"])
    rm(filepath, force=true)

    # throws error with chrX
    data = """
    CHR\tPOS\tREF\tALT\tZ
    1\t758351\tA\tG\t1.05551882016081
    2\t779885\tC\tT\t2.12197306477
    3\t779987\tA\tG\t1.95791489337
    X\t782105\tC\tA\t1.91243829548
    """
    CSV.write(filepath, CSV.File(IOBuffer(data)))
    @test_throws ArgumentError read_zscores(filepath)

    # missing Z scores
    data = """
    CHR\tPOS\tREF\tALT\tZ
    1\t758351\tA\tG\t1.05551882016081
    1\t779885\tC\tT\tNaN
    1\t779987\tA\tG\t1.95791489337
    1\t782105\tC\tA\t
    """
    CSV.write(filepath, CSV.File(IOBuffer(data)))
    z, chr, pos, effect_allele, non_effect_allele = read_zscores(filepath)
    @test all(z .≈ [1.05551882016081, 1.95791489337])
    @test all(chr .== 1)
    @test all(pos .== [758351, 779987])
    @test all(effect_allele .== ["G", "G"])
    @test all(non_effect_allele .== ["A", "A"])

    # extra columns with possibly different names
    data = """
    rsid\tCHR\tPOS\tREF\tALT\tZ\tZreal
    rs123\t1\t758351\tA\tG\t1.05\t111
    rs234\t2\t779885\tC\tT\tNaN\t111
    rs345\t3\t779987\tA\tG\t1.95\t111
    rs456\t4\t782105\tC\tA\t\t111
    """
    CSV.write(filepath, CSV.File(IOBuffer(data)))
    z, chr, pos, effect_allele, non_effect_allele = read_zscores(filepath, z_col=7)
    @test all(z .≈ [111, 111, 111, 111])
    @test all(chr .== [1, 2, 3, 4])
    @test all(pos .== [758351, 779885, 779987, 782105])
    @test all(effect_allele .== ["G", "T", "G", "A"])
    @test all(non_effect_allele .== ["A", "C", "A", "C"])

    # repeated SNPs
    # data = """
    # CHR\tPOS\tREF\tALT\tZ
    # 1\t758351\tA\tG\t1.05551882016081
    # 1\t758351\tA\tG\t1.05551882016081
    # 1\t779885\tC\tT\t2.12197306477
    # 1\t779885\tT\tC\t2.12197306477
    # """
    # CSV.write(filepath, CSV.File(IOBuffer(data)))
    # read_zscores(filepath)

    # cleanup
    rm(filepath)
end

@testset "qvalues" begin
    m = 5
    kappa = [0, 0, 1, 5, 2, 0, 2]
    groups = [1, 1, 2, 3, 4, 1, 5]
    tau = [0.1, 0.21, 0.01, 0.05, 0.02, 0.04, 0.2]
    qvalues = GhostKnockoffGWAS.get_knockoff_qvalue(kappa, tau, m, groups=groups)

    # true answer compared with Zihuai's code
    @test all(qvalues .≈ [0.3, 0.2, 1.0, 1.0, 1.0, 0.33333333333333333, 1.0])
end

@testset "solve_blocks (within julia) with VCF input" begin
    # download a test VCF file if not exist
    # This VCF file 
        # - have ~1300 records
        # - some rsIDs are missing
        # - some GT fields are missing?
        # - chr field as integers
    vcffile = joinpath(dirname(pathof(GhostKnockoffGWAS)), "..", 
        "test/test.08Jun17.d8b.vcf.gz")
    isfile(vcffile) || Downloads.download(
        "http://faculty.washington.edu/browning/beagle/test.08Jun17.d8b.vcf.gz",
        vcffile)

    # get VCF summary
    tmpfile = tempname()
    gtstats(vcffile, tmpfile)
    df = CSV.read(tmpfile, DataFrame, 
        header=[:chr, :pos, :id, :ref, :alt, :qual, :filt, :info, :missings, 
        :missfreq, :nalt, :altfreq, :nminor, :maf, :hwe])
    rm(tmpfile, force=true)
    @test size(df, 1) == nrecords(vcffile) == 1356

    # read the region into memory and solve knockoff optimization problem
    chr = 22
    start_bp = 1
    end_bp = 999999999999
    outdir = joinpath(dirname(pathof(GhostKnockoffGWAS)), "..", "test/LD_files")
    isdir(outdir) || mkpath(outdir)
    hg_build = 19
    @time solve_blocks(vcffile, chr, start_bp, end_bp, outdir, hg_build, 
        min_maf = 0.01, min_hwe = 0.0)

    # test basic output structure
    @test isdir(outdir)
    @test isdir(joinpath(outdir, "chr22"))
    @test length(readdir(joinpath(outdir, "chr22"))) == 3

    # test Info file
    file = joinpath(outdir, "chr22", "Info_start$(start_bp)_end$(end_bp).csv")
    @test isfile(file)
    info = CSV.read(file, DataFrame)
    @test "pos_hg$hg_build" in names(info)
    @test all(x -> 0.01 ≤ x ≤ 0.99, info[!, "AF"])

    # test summary file
    file = joinpath(outdir, "chr22", "summary_start$(start_bp)_end$(end_bp).csv")
    @test isfile(file)
    summary = CSV.read(file, DataFrame, header=false)
    summary[findfirst(x -> x == "chr", summary[!, 1]), 2] == chr
    summary[findfirst(x -> x == "start_bp", summary[!, 1]), 2] == start_bp
    summary[findfirst(x -> x == "end_bp", summary[!, 1]), 2] == end_bp
    summary[findfirst(x -> x == "m", summary[!, 1]), 2] == 5
    summary[findfirst(x -> x == "p", summary[!, 1]), 2] == size(info, 1)

    # test h5 file
    file = joinpath(outdir, "chr22", "LD_start$(start_bp)_end$(end_bp).h5")
    @test isfile(file)
    h5reader = h5open(file, "r")
    D = read(h5reader, "D")
    S = read(h5reader, "S")
    Sigma = read(h5reader, "Sigma")
    @test eigmin(Symmetric(6/5*Sigma - D)) > 0
    @test size(D) == size(Sigma)
    @test size(D, 1) == summary[findfirst(x -> x == "p", summary[!, 1]), 2]
    @test size(S, 1) == summary[findfirst(x -> x == "nreps", summary[!, 1]), 2]
end

@testset "solve_blocks (within julia) with PLINK input" begin
    # download plink 1.9 software
    remote = "https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip"
    file = joinpath(dirname(pathof(GhostKnockoffGWAS)), "..", 
        "test/plink_linux_x86_64_20231211.zip")
    Downloads.download(remote, file)
    run(`unzip $file`)

    # convert VCF file to binary PLINK using PLINK 1.9
    plink_exe = joinpath(dirname(pathof(GhostKnockoffGWAS)), "..", "test/plink")
    vcffile = joinpath(dirname(pathof(GhostKnockoffGWAS)), "..", 
        "test/test.08Jun17.d8b.vcf.gz")
    plinkfile = joinpath(dirname(pathof(GhostKnockoffGWAS)), "..", 
        "test/test.08Jun17.d8b")
    vcf_to_plink(plink_exe, vcffile, plinkfile)

    # check converted numeric matrix is the same after 0 <--> 2 flipping
    # (see vcf_to_plink function for why there is a 0 <--> 2 flipping)
    Xtrue = convert_gt(Float64, vcffile, impute=true)
    Xtest = convert(Matrix{Float64}, SnpArray(plinkfile * ".bed"), impute=true)
    @test size(Xtrue) == size(Xtest)
    @test all(Xtrue .== 2 .- Xtest)

    # run solve_block on PLINK file
    chr = 22
    start_bp = 1
    end_bp = 999999999999
    outdir = joinpath(dirname(pathof(GhostKnockoffGWAS)), "..", "test/LD_files2")
    isdir(outdir) || mkpath(outdir)
    hg_build = 19
    @time solve_blocks(plinkfile * ".bed", chr, start_bp, end_bp, outdir,
        hg_build, min_maf = 0.01, min_hwe = 0.0)

    # check PLINK vs VCF output is the same
    VCF_outdir = joinpath(dirname(pathof(GhostKnockoffGWAS)), "..", "test/LD_files")
    VCF_h5 = joinpath(VCF_outdir, "chr22", "LD_start$(start_bp)_end$(end_bp).h5")
    VCFh5reader = h5open(VCF_h5, "r")
    PLINK_outdir = joinpath(dirname(pathof(GhostKnockoffGWAS)), "..", "test/LD_files2")
    PLINK_h5 = joinpath(PLINK_outdir, "chr22", "LD_start$(start_bp)_end$(end_bp).h5")
    PLINKh5reader = h5open(PLINK_h5, "r")
    @test all(read(PLINKh5reader, "Sigma") .≈ read(VCFh5reader, "Sigma"))
    @test all(read(PLINKh5reader, "groups") .== read(VCFh5reader, "groups"))

    close(VCFh5reader)
    close(PLINKh5reader)

    # seems like group_reps could be slightly different due to floating point precision
    # (we have floating point precision issues due to PLINK's X is 2 .- X in VCF)
    # @test all(read(PLINKh5reader, "group_reps") .== read(VCFh5reader, "group_reps"))
    # @test all(read(PLINKh5reader, "D") .≈ read(VCFh5reader, "D"))
    # @test all(read(PLINKh5reader, "S") .≈ read(VCFh5reader, "S"))
end

@testset "ghostbasil (within julia)" begin
    testdir = joinpath(dirname(pathof(GhostKnockoffGWAS)), "../test")
    k = 10
    mu = 0
    sigma = 5.0 # beta ~ N(mu, sigma)

    # import VCF and remove SNPs with MAF < 0.1
    vcffile = joinpath(testdir, "test.08Jun17.d8b.vcf.gz")
    X, sampleID, chr, pos, rsid, ref, alt = 
        convert_gt(Float64, vcffile, impute=true, center=false, scale=false, 
        save_snp_info=true)
    mafs = mean.(skipmissing.(eachcol(X))) ./ 2
    idx = findall(x -> 0.01 < x < 0.99, mafs)
    X = X[:, idx]
    chr, pos, rsid, ref, alt = chr[idx], pos[idx], vcat(rsid...)[idx], 
        ref[idx], vcat(alt...)[idx]
    n, p = size(X)

    # simulate phenotypes and normalize it
    Random.seed!(2024)
    beta = zeros(p)
    beta[1:k] .= rand(Normal(mu, sigma), k)
    shuffle!(beta)
    y = X * beta + randn(n)
    zscore!(y, mean(y), std(y))

    # compute z scores and save for later use
    z = X'*y ./ sqrt(n)
    zfile = joinpath(testdir, "zfile.txt")
    LD_files = joinpath(testdir, "LD_files")
    CSV.write(zfile, DataFrame("CHR"=>chr,"POS"=>pos,"REF"=>ref,"ALT"=>alt,"Z"=>z))

    # GhostKnockoffGWAS function
    LDfiles = joinpath(dirname(pathof(GhostKnockoffGWAS)), "..", "test/LD_files")
    chr = parse.(Int, chr)
    N = size(X, 1)
    hg_build = 19
    outdir = joinpath(dirname(pathof(GhostKnockoffGWAS)), "..", "test")
    ghostknockoffgwas(LDfiles, z, chr, pos, alt, ref, N, hg_build, outdir)
    @test isfile(joinpath(outdir, "result.txt"))
    @test isfile(joinpath(outdir, "result_summary.txt"))
    result = CSV.read(joinpath(outdir, "result.txt"), DataFrame)
    summary = CSV.read(joinpath(outdir, "result_summary.txt"), DataFrame, header=false)
    # @test summary[findfirst(x -> x == "target_fdr_0.1_num_selected", summary[!, 1]), 2] == 0
    # @test summary[findfirst(x -> x == "target_fdr_0.2_num_selected", summary[!, 1]), 2] == 8
    @test summary[findfirst(x -> x == "m", summary[!, 1]), 2] == 5
    @test summary[findfirst(x -> x == "nregions", summary[!, 1]), 2] == 1
    @test summary[findfirst(x -> x == "nsnps", summary[!, 1]), 2] == 402
    @test summary[findfirst(x -> x == "m", summary[!, 1]), 2] == 5
    @test summary[findfirst(x -> x == "lasso_lambda", summary[!, 1]), 2] ≈ 0.14776789851770086
    # @test summary[findfirst(x -> x == "mean_LD_shrinkage", summary[!, 1]), 2] ≈ 3.324562949461316e-16
end

@testset "Download, unpack, and run exe" begin
    wd = pwd()

    # download and unpack 
    cd(joinpath(dirname(pathof(GhostKnockoffGWAS)), "..", "test"))
    run(`wget https://github.com/biona001/GhostKnockoffGWAS/releases/download/v0.2.1/app_linux_x86.tar.gz`)
    run(`tar -xvzf app_linux_x86.tar.gz`)
    @test isdir("app_linux_x86")
    @test isfile("app_linux_x86/bin/GhostKnockoffGWAS")
    @test isfile("app_linux_x86/bin/solveblock")

    # help messages
    help1 = run(`./app_linux_x86/bin/solveblock -h`)
    help2 = run(`./app_linux_x86/bin/GhostKnockoffGWAS -h`)
    @test help1.exitcode == 0
    @test help2.exitcode == 0

    # solveblock executable
    exe = joinpath(dirname(pathof(GhostKnockoffGWAS)), "../test/app_linux_x86/bin/solveblock")
    vcffile = joinpath(dirname(pathof(GhostKnockoffGWAS)), "..", "test/test.08Jun17.d8b.vcf.gz")
    chr = 22
    start_bp = 1
    end_bp = 999999999999
    outdir = joinpath(dirname(pathof(GhostKnockoffGWAS)), "..", "test/LD_files")
    isdir(outdir) || mkpath(outdir)
    hg_build = 19
    sb = run(`$exe --vcffile $vcffile --chr $chr --start_bp $start_bp --end_bp $end_bp --outdir $outdir --genome-build $hg_build`)
    @test sb.exitcode == 0

    # GhostKnockoffGWAS executable
    exe = joinpath(dirname(pathof(GhostKnockoffGWAS)), "../test/app_linux_x86/bin/GhostKnockoffGWAS")
    LDfiles = joinpath(dirname(pathof(GhostKnockoffGWAS)), "..", "test/LD_files")
    outfile = "GK_out"
    zfile = "zfile.txt"
    gk = run(`$exe --zfile $zfile --LD-files $LDfiles --N 191 --genome-build 19 --out $outfile`)
    @test gk.exitcode == 0

    # cleanup
    rm("app_linux_x86.tar.gz", force=true)
    rm("app_linux_x86", force=true, recursive=true)
    rm("GK_out_summary.txt", force=true)
    rm("GK_out.txt", force=true)
    rm("LD_files", force=true, recursive=true)
    rm("LD_files2", force=true, recursive=true)
    rm("result_summary.txt", force=true)
    rm("result.txt", force=true)
    rm("test.08Jun17.d8b.vcf.gz", force=true)
    rm("zfile.txt", force=true)
    rm("plink_linux_x86_64_20231211.zip", force=true)
    rm("LICENSE", force=true)
    rm("plink", force=true)
    rm("prettify", force=true)
    rm("test.08Jun17.d8b.bed", force=true)
    rm("test.08Jun17.d8b.bim", force=true)
    rm("test.08Jun17.d8b.fam", force=true)
    rm("test.08Jun17.d8b.nosex", force=true)
    rm("test.08Jun17.d8b.log", force=true)
    rm("toy.map", force=true)
    rm("toy.ped", force=true)
    cd(wd)
end

@testset "solve_blocks on bad VCF files" begin
    # breaks when chr have non-Int values, e.g. "chr22"

    # breaks when POS is not sorted within chrs (but not between chrs)

    # breaks when >1 ALT allele

    # breaks when too little sample size

    # breaks when too little SNPs within a region

    # breaks when outdir doesn't exist

    # is able to skip to chr2
end
