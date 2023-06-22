# convert p-values and effect sizes to Z scores
function pval2zscore(pvals::AbstractVector{T}, beta::AbstractVector{T}) where T
    length(pvals) == length(beta) || 
        error("pval2zscore: pvals and beta should have the same length")
    return zscore.(pvals, beta)
end
zscore(p::T, beta::T) where T = sign(beta) * quantile(Normal(), p/2)

# converts Z score to p-values
zscore2pval(z::AbstractVector{T}) where T = pval.(z)
pval(z::T) where T = 2ccdf(Normal(), abs(z))

"""
    find_matching_indices(a, b)

Returns output vector `c` such that `c[i]` contains all indices of `b` that 
matches `a[i]`. Generalizes `indexin` in Base. 

# Example
```julia
a = [1, 2, 3, 2, 6]
b = [2, 1, 2, 4, 3]
c = find_matching_indices(a, b)
println(c)  # Output: [[2], [1, 3], [5], [1, 3], Int[]]
```
"""
function find_matching_indices(a::AbstractVector, b::AbstractVector)
    dict_b = Dict()
    for (i, val) in enumerate(b)
        if haskey(dict_b, val)
            push!(dict_b[val], i)
        else
            dict_b[val] = [i]
        end
    end
    c = Vector{Vector{Int}}(undef, length(a))
    for (i, val) in enumerate(a)
        if haskey(dict_b, val)
            c[i] = dict_b[val]
        else
            c[i] = []
        end
    end
    return c
end

