# convert p-values and effect sizes to Z scores
function pval2zscore(pvals::Vector{T}, beta::Vector{T}) where T
    length(pvals) == length(beta) || 
        error("pval2zscore: pvals and beta should have the same length")
    return zscore.(pvals, beta)
end
zscore(p::T, beta::T) where T = sign(beta) * quantile(Normal(), p/2)

# converts Z score to p-values
zscore2pval(z::Vector{T}) where T = pval.(z)
pval(z::T) where T = 2ccdf(Normal(), abs(z))

