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
    c = find_matching_indices(a, b)
    @test c == [[2], [1, 3], [5], [1, 3], Int[]]
    @test typeof(c) <: Vector{Vector{Int64}}
end

@testset "LD shrinkage" begin
    function mvn_logl_under_null_naive(Σ::AbstractMatrix, z::AbstractVector)
        return -0.5logdet(Symmetric(Σ)) - dot(z, inv(Symmetric(Σ)), z)
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

    γ = find_optimal_shrinkage(Sigma, z)
    @test 0 ≤ γ ≤ 1
end