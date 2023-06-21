@testset "utilities" begin
    # convert between z score and p-value
    p = 0.00695184371278
    beta = -0.010602
    z = 2.6991423324213346
    @test zscore(p, beta) ≈ 2.6991423324213346
    @test pval(z) ≈ 0.00695184371278
end
