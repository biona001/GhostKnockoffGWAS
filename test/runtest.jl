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
