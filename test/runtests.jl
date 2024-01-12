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
    @test zscore(p, beta) ≈ 2.6991423324213346
    @test pval(z) ≈ 0.00695184371278

    # find_matching_indices
    a = [1, 2, 3, 2, 6]
    b = [2, 1, 2, 4, 3]
    c = GhostKnockoffGWAS.find_matching_indices(a, b)
    @test c == [[2], [1, 3], [5], [1, 3], Int[]]
    @test typeof(c) <: Vector{Vector{Int64}}
end

@testset "LD shrinkage" begin
    function mvn_logl_under_null_naive(Σ::AbstractMatrix, z::AbstractVector)
        return -0.5logdet(Symmetric(Σ)) - dot(z, inv(Symmetric(Σ)), z)
    end
    function mvn_logl_under_null(L::Cholesky, z::AbstractVector, u=zeros(length(z)))
        ldiv!(u, UpperTriangular(L.factors)', z)
        return -0.5logdet(L) - dot(u, u)
    end

    p = 1000
    z = randn(p)
    x = randn(p, p)
    storage=zeros(length(z))
    Sigma = Symmetric(x'*x)
    L = cholesky(Symmetric(Sigma));
    a = mvn_logl_under_null_naive(Sigma, z)
    b = mvn_logl_under_null(L, z)
    @test a ≈ b

    γ = GhostKnockoffGWAS.find_optimal_shrinkage(Sigma, z)
    @test 0 ≤ γ ≤ 1
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
