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
        0, 1.0, Brent(), show_trace=false, 
        iterations = 50,
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
    return neg_mvn_logl_under_null((1-γ)*Σ + γ*I, z)
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

# counts number of Z scores that can be matched to LD panel
# ~400 seconds is running on typed SNPs only
function count_matchable_snps(
    knockoff_dir::String,  # Directory that stores knockoff results (i.e. output from part 1)
    z::Vector{Float64},    # Z scores
    chr::Vector{Int},      # chromosome of each Z score
    pos::Vector{Int},      # position of each Z score (specify hg build with hg_build)
    effect_allele::Vector{String},       # effect allele of Z score
    non_effect_allele::Vector{String},   # non-effect allele of Z score
    hg_build::Int,
    target_chrs=1:22,
    )
    nregions, nsnps, nknockoff_snps = 0, 0, 0
    for c in target_chrs
        files = readdir(joinpath(knockoff_dir, "chr$c"))
        chr_idx = findall(x -> x == c, chr)
        GWAS_pos = pos[chr_idx]
        GWAS_ea = effect_allele[chr_idx]
        GWAS_nea = non_effect_allele[chr_idx]
        zscores = z[chr_idx]
        for f in files
            endswith(f, ".h5") || continue
            fname = f[4:end-3]

            # read knockoff results
            Sigma_info = CSV.read(
                joinpath(knockoff_dir, "chr$(c)", "Info_$fname.csv"), DataFrame
            )
            nknockoff_snps += size(Sigma_info, 1)
            # map reference LD panel to GWAS Z-scores by position
            LD_pos = Sigma_info[!, "pos_hg$(hg_build)"]
            shared_snps = intersect(LD_pos, GWAS_pos)
            # delete SNPs if ref/alt don't match
            remove_idx = Int[]
            for (i, snp) in enumerate(shared_snps)
                GWAS_idx = findfirst(x -> x == snp, GWAS_pos)
                LD_idx = findfirst(x -> x == snp, LD_pos)
                ref_match_ea = Sigma_info[LD_idx, "ref"] == GWAS_ea[GWAS_idx]
                alt_match_nea = Sigma_info[LD_idx, "alt"] == GWAS_nea[GWAS_idx]
                ref_match_nea = Sigma_info[LD_idx, "ref"] == GWAS_nea[GWAS_idx]
                alt_match_ea = Sigma_info[LD_idx, "alt"] == GWAS_ea[GWAS_idx]
                if ref_match_ea && alt_match_nea 
                    continue
                elseif ref_match_nea && alt_match_ea
                    # push!(remove_idx, i)
                    zscores[GWAS_idx] *= -1
                else # SNP cannot get matched to LD panel
                    push!(remove_idx, i)
                end
            end
            deleteat!(shared_snps, unique!(remove_idx))

            # update counters
            nsnps += length(shared_snps)
            nregions += 1
        end
        println("count_matchable_snps processed chr $c, cumulative SNPs = $nsnps")
    end
    if nsnps / nknockoff_snps < 0.05
        error("Less than 5% of SNPs in the pre-computed knockoff LD panel " * 
              "can be successfully mapped to target Z scores. Please check " *
              "if the human genome build of the target study is $hg_build")
    end
    return nsnps
end

function read_zscores(filepath::String)
    # read chr, pos, ref/alt alleles, and Z scores
    chr, pos, effect_allele, non_effect_allele, z = try
        info = CSV.read(filepath, DataFrame)
        chr = info[!, "CHR"] |> Vector{Int}
        pos = info[!, "POS"] |> Vector{Int}
        effect_allele = info[!, "ALT"] |> Vector{String}
        non_effect_allele = info[!, "REF"] |> Vector{String}
        z = info[!, "Z"] |> Vector{Float64}
        chr, pos, effect_allele, non_effect_allele, z
    catch
        error(
            "Error reading Z score file $filepath. Does the file contain " *
            "CHR, POS, REF, ALT, and Z as headers?"
        )
    end

    # remove NaN/Inf
    idx = findall(x -> !isnan(x) && !isinf(x), z)

    return z[idx], chr[idx], pos[idx], effect_allele[idx], non_effect_allele[idx]
end
