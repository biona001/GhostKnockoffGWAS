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

function find_optimal_shrinkage(Σ::AbstractMatrix, z::AbstractVector)
    opt = optimize(
        γ -> neg_mvn_logl_under_null(Σ, z, γ), 
        0, 0.25, Brent(), show_trace=true, abs_tol=0.000001,
        iterations = 20,
    )
    return opt.minimizer
end

# function neg_mvn_logl_under_null_naive(Σ::AbstractMatrix, z::AbstractVector)
#     return 0.5logdet(Symmetric(Σ)) + dot(z, inv(Symmetric(Σ)), z)
# end
# function neg_mvn_logl_under_null_naive(Σ::AbstractMatrix, z::AbstractVector, γ::Number)
#     return mvn_logl_under_null_naive((1-γ)*Σ + γ*I, z)
# end
# @time neg_mvn_logl_under_null_naive(Σ, zscore_tmp)

function neg_mvn_logl_under_null(Σ::AbstractMatrix, z::AbstractVector, γ::Number)
    return mvn_logl_under_null((1-γ)*Σ + γ*I, z)
end
function neg_mvn_logl_under_null(Σ::AbstractMatrix, z::AbstractVector)
    L = cholesky(Symmetric(Σ))
    u = zeros(length(z))
    ldiv!(u, UpperTriangular(L.factors)', z)
    return 0.5logdet(L) + dot(u, u)
end
# @time neg_mvn_logl_under_null(Σ, zscore_tmp)
# @time neg_mvn_logl_under_null(Σ, zscore_tmp, 1.0)
# @time neg_mvn_logl_under_null(Σ, zscore_tmp, 0.5)
# @time neg_mvn_logl_under_null(Σ, zscore_tmp, 0.1)
# @time neg_mvn_logl_under_null(Σ, zscore_tmp, 0.0)
# γ = find_optimal_shrinkage(Σ, zscore_tmp)


# function mvn_logl_under_null(L::Cholesky, z::AbstractVector, u=zeros(length(z)))
#     ldiv!(u, UpperTriangular(L.factors)', z)
#     return -0.5logdet(L) - dot(u, u)
# end
