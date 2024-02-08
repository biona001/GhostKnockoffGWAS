using Test
using GhostKnockoffGWAS
using Ghostbasil
using LinearAlgebra
using DelimitedFiles
using CSV
using DataFrames

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

@testset "ghostbasil C++ solver" begin
    # basic test to make sure block_group_ghostbasil does not crash
    p = 10 # number of features
    m = 5  # number of knockoffs
    x = randn(p, p)
    y = randn(p, p)
    Ci = x'*x
    Si = y'*y
    r = randn(6p)
    lambda_path = [0.1, 0.05, 0.01, 0.001]
    beta_i = block_group_ghostbasil(Ci, Si, r, lambda_path, m=m)
    @test length(beta_i) == (m+1)*p

    # Test correctness with a specific region from Pan-UKBB. 
    # The true answer is obtained from running the ghostbasil R package on
    # the given data. Note these data are based on LD summary statistics, 
    # which can be freely distributed
    datadir = normpath(Ghostbasil.datadir())
    C = readdlm(joinpath(datadir, "C.txt"))
    S = readdlm(joinpath(datadir, "S.txt"))
    r = readdlm(joinpath(datadir, "r.txt")) |> vec
    lambda_path = readdlm(joinpath(datadir, "lambda_path.txt")) |> vec
    beta_i = block_group_ghostbasil(C, S, r, lambda_path)
    @test count(!iszero, beta_i) == 14
    non_zero_idx = findall(!iszero, beta_i)
    @test non_zero_idx == [215, 340, 341, 431, 559, 560, 763, 878,
        1085, 1528, 1770, 2726, 3157, 3387]
    @test all(beta_i[non_zero_idx] .≈ 
        [-0.0005451576807477974, -0.0002947665358369119, -0.014810108231841996, 
        -0.00043238505905553924, -0.006856782960094195, 0.0016850133196618926, 
        -0.00019274743491070077, -0.0003841629682186662, -0.0005758897376517483, 
        -0.0016095272344540995, 0.0020146776266541312, -0.001666736964454875, 
        0.0004269118551173737, 0.00010746291189108347]
    )
end

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
