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

# eq 24 of https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010299
function find_optimal_shrinkage(Σ::AbstractMatrix, z::AbstractVector)
    opt = optimize(
        γ -> neg_mvn_logl_under_null(Σ, z, γ), 
        0, 1.0, Brent(), show_trace=false, 
        iterations = 50,
    )
    return opt.minimizer
end
function neg_mvn_logl_under_null(Σ::AbstractMatrix, z::AbstractVector, γ::Number)
    return neg_mvn_logl_under_null((1-γ)*Σ + γ*I, z)
end
function neg_mvn_logl_under_null(Σ::AbstractMatrix, z::AbstractVector)
    L = cholesky(Symmetric(Σ))
    u = zeros(length(z))
    ldiv!(u, UpperTriangular(L.factors)', z) # non-allocating ldiv!(u, L.L, z)
    return 0.5logdet(L) + dot(u, u)
end

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
                if ref_match_nea && alt_match_ea
                    continue
                elseif ref_match_ea && alt_match_nea
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
    if nsnps / nknockoff_snps < 0.01
        error("Less than 1% of SNPs in the pre-computed knockoff LD panel " * 
              "can be successfully mapped to target Z scores. Please check " *
              "if the human genome build of the target study is $hg_build")
    end
    return nsnps
end

"""
    read_zscores(filepath::String)

Helper function to read a Z-score file. First row must be a header column with 
CHR, POS, REF, ALT, and Z. All other columns will be ignored. 

# todo: detect duplicate SNPs
"""
function read_zscores(filepath::String)
    # create CSV.File object (this reads the file lazily)
    csv_file = CSV.File(filepath)

    # check the 5 required header is present
    if !issubset([:CHR, :POS, :REF, :ALT, :Z], propertynames(csv_file))
        error("Error reading Z score file $filepath. Does the file contain " *
            "CHR, POS, REF, ALT, and Z as headers?")
    end

    # read chr, pos, ref/alt alleles, and Z scores
    chr = csv_file.CHR
    pos = csv_file.POS
    effect_allele = csv_file.ALT
    non_effect_allele = csv_file.REF
    z = csv_file.Z

    # chr must be integer valued
    if eltype(chr) <: AbstractString
        try
            chr = parse.(Int, chr)
        catch e
            println(
                "Error parsing CHR field of Z score file $filepath. " * 
                "Note that CHR field must be integer valued (e.g. chr22 " *
                "and sex chromosomes are NOT valid!)"
            )
            rethrow(e)
        end
    end

    # Z must be Float64 valued
    if eltype(z) <: AbstractString
        try
            z = parse.(Float64, chr)
        catch e
            println(
                "Error parsing z-scores of Z score file $filepath. " * 
                "Note that Z scores should be floating point valued. Missing " * 
                "Z scores can be specified as NaN or as an empty cell."
            )
            rethrow(e)
        end
    end

    # detect duplicate SNPs
    check_all_snps_are_unique(chr, pos, effect_allele, non_effect_allele)

    # find missing/NaN/Inf
    chr_idx = findall(x -> !ismissing(x) && !isnothing(x) && !isnan(x) && !isinf(x), chr)
    pos_idx = findall(x -> !ismissing(x) && !isnothing(x) && !isnan(x) && !isinf(x), pos)
    ref_idx = findall(x -> !ismissing(x) && !isnothing(x), non_effect_allele)
    alt_idx = findall(x -> !ismissing(x) && !isnothing(x), effect_allele)
    z_idx = findall(x -> !ismissing(x) && !isnothing(x) && !isnan(x) && !isinf(x), z)
    idx = intersect(chr_idx, pos_idx, ref_idx, alt_idx, z_idx)

    # filter Z/CHR/POS/REF/ALT for valid values
    z = z[idx] |> Vector{Float64}
    chr = chr[idx] |> Vector{Int}
    pos = pos[idx] |> Vector{Int}
    effect_allele = effect_allele[idx] |> Vector{String}
    non_effect_allele = non_effect_allele[idx] |> Vector{String}

    return z, chr, pos, effect_allele, non_effect_allele
end

function check_all_snps_are_unique(chr, pos, effect_allele, non_effect_allele)
    return nothing # todo
end
