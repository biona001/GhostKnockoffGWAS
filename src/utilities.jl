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
function neg_mvn_logl_under_null(evals, evecs, z::AbstractVector, γ::Number)
    z_scaled = evecs' * z
    evals_scaled = (1-γ) .* evals .+ γ
    return sum(log.(evals_scaled)) + dot(z_scaled, Diagonal(1 ./ evals_scaled), z_scaled)
end
function find_optimal_shrinkage(Σ::AbstractMatrix, z::AbstractVector)
    size(Σ, 1) == size(Σ, 2) == length(z) || 
        error("find_optimal_shrinkage: Dimension mismatch")
    evals, evecs = eigen(Symmetric(Σ))
    opt = optimize(
        γ -> neg_mvn_logl_under_null(evals, evecs, z, γ), 
        0, 1.0, Brent()
    )
    return opt.minimizer
end

# counts number of Z scores that can be matched to LD panel
# ~400 seconds is running on typed SNPs only
function count_matchable_snps(
    LD_files::String,  # Directory that stores LD/knockoff files (i.e. output from part 1)
    z::Vector{Float64},    # Z scores
    chr::Vector{Int},      # chromosome of each Z score
    pos::Vector{Int},      # position of each Z score (specify hg build with hg_build)
    effect_allele::Vector{String},       # effect allele of Z score
    non_effect_allele::Vector{String},   # non-effect allele of Z score
    hg_build::Int,
    target_chrs=sort!(unique(chr)),
    )
    nregions, nsnps, nknockoff_snps = 0, 0, 0
    for c in target_chrs
        files = readdir(joinpath(LD_files, "chr$c"))
        chr_idx = findall(x -> x == c, chr)
        GWAS_pos = pos[chr_idx]
        GWAS_ea = effect_allele[chr_idx]
        GWAS_nea = non_effect_allele[chr_idx]
        zscores = z[chr_idx]
        for f in files
            endswith(f, ".h5") || continue
            fname = f[4:end-3]

            # read LD/knockoff files
            Sigma_info = CSV.read(
                joinpath(LD_files, "chr$(c)", "Info_$fname.csv"), DataFrame
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
    return nsnps, nregions
end

function count_total_snps(
    LD_files::String, # Directory that stores LD/knockoff files
    target_chrs=1:22,
    ) 
    nsnps = 0
    for c in target_chrs
        files = readdir(joinpath(LD_files, "chr$c"))
        for f in files
            startswith(f, "Info") || continue
            nsnps += countlines(joinpath(LD_files, "chr$c", f)) - 1
        end
    end
    return nsnps
end

"""
    read_zscores(filepath::String)

Helper function to read a Z-score file at `filepath`. This function is mainly 
intended for Julia users running `GhostKnockoffGWAS` in the REPL.

# Input
+ `filepath`: Full file path to the Z-score file. First row must be a header
    column with `CHR`, `POS`, `REF`, `ALT`, and `Z`. If the file contains these 
    information but has a different header name for them, use the optional input
    arguments below. All other columns will be ignored.

# Optional inputs
+ `chr_col`: An integer, specifying which column in `filepath` should be read as
    CHR (by default we search for a header `CHR`)
+ `pos_col`: An integer, specifying which column in `filepath` should be read as 
    POS (by default we search for a header `POS`)
+ `ref_col`: An integer, specifying which column in `filepath` should be read as 
    REF (by default we search for a header `REF`)
+ `alt_col`: An integer, specifying which column in `filepath` should be read as
    ALT (by default we search for a header `ALT`)
+ `z_col`: An integer, specifying which column in `filepath` should be read as Z
    (by default we search for a header `Z`)

# Output
+ `z`: The Z scores stored in the `Z` column of `filepath`
+ `chr`: The chromosome number stored in `CHR` column of `filepath`. Only integer
    values are allowed.
+ `pos`: The position number stored in `POS` column of `filepath`.
+ `effect_allele`: The allele stored in `ALT` column of `filepath`.
+ `non_effect_allele`: The allele stored in `REF` column of `filepath`.
"""
function read_zscores(
    filepath::String;
    chr_col::Union{Nothing, Int}=nothing, 
    pos_col::Union{Nothing, Int}=nothing, 
    ref_col::Union{Nothing, Int}=nothing, 
    alt_col::Union{Nothing, Int}=nothing, 
    z_col::Union{Nothing, Int}=nothing
    )
    # create CSV.File object (this reads the file lazily)
    csv_file = CSV.File(filepath)

    # read chr, pos, ref/alt alleles, and Z scores
    chr = try
        isnothing(chr_col) ? csv_file.CHR : [csv_file[i][chr_col] for i in eachindex(csv_file)]
    catch
        error("Error reading CHR from $filepath.")
    end
    pos = try
        isnothing(pos_col) ? csv_file.POS : [csv_file[i][pos_col] for i in eachindex(csv_file)]
    catch
        error("Error reading POS from $filepath.")
    end
    non_effect_allele = try
        isnothing(ref_col) ? csv_file.REF : [csv_file[i][ref_col] for i in eachindex(csv_file)]
    catch
        error("Error reading REF from $filepath.")
    end
    effect_allele = try
        isnothing(alt_col) ? csv_file.ALT : [csv_file[i][alt_col] for i in eachindex(csv_file)]
    catch
        error("Error reading ALT from $filepath.")
    end
    z = try
        isnothing(z_col) ? csv_file.Z : [csv_file[i][z_col] for i in eachindex(csv_file)]
    catch
        error("Error reading Z scores from $filepath.")
    end

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
            # try to parse some commonly used values for missing, even though 
            # in the docs we don't allow it
            replace!(z, ""=>"NaN")
            replace!(z, "NA"=>"NaN")
            replace!(z, "Na"=>"NaN")
            replace!(z, "na"=>"NaN")
            replace!(z, "missing"=>"NaN")
            replace!(z, "Missing"=>"NaN")
            z = parse.(Float64, z)
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
    # check_all_snps_are_unique(chr, pos, effect_allele, non_effect_allele)

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

"""
    sample_mvn_efficient(C::AbstractMatrix{T}, D::AbstractMatrix{T}, m::Int)

Efficiently samples from `N(0, A)` where
```math
\\begin{aligned}
A &= \\begin{pmatrix}
    C & C-D & \\cdots & C-D\\\\
    C-D & C & \\cdots & C-D\\\\
    \\vdots & & \\ddots & \\vdots\\\\
    C-D & C-D & & C
\\end{pmatrix}
\\end{aligned}
```
Note there are `m` blocks per row/col

# Source
https://github.com/biona001/Knockoffs.jl/blob/master/src/ghost.jl#L60
"""
function sample_mvn_efficient(C::AbstractMatrix{T}, D::AbstractMatrix{T}, m::Int) where T
    p = size(C, 1)
    L = cholesky(Symmetric(C - (m-1)/m * D))
    e1 = randn(p)
    e2 = Vector{T}[]
    d = MvNormal(Symmetric(D))
    for i in 1:m
        push!(e2, rand(d))
    end
    e2_avg = 1/m * sum(e2)
    Zko = T[]
    for i in 1:m
        append!(Zko, L.L*e1 + e2[i] - e2_avg)
    end
    return Zko
end

"""
    ghost_knockoffs(Zscores, D, Σinv; [m=1])

Generate Ghost knockoffs given a list of z-scores (GWAS summary statistic). 

# Inputs
+ `Zscores`: List of z-score statistics
+ `D`: Matrix obtained from solving the knockoff problem satisfying 
    `(m+1)/m*Σ - D ⪰ 0`
+ `Σinv`: Inverse of the covariance matrix

# optional inputs
+ `m`: Number of knockoffs

# Reference
He, Z., Liu, L., Belloy, M. E., Le Guen, Y., Sossin, A., Liu, X., ... & Ionita-Laza, I. (2021). 
Summary statistics knockoff inference empowers identification of putative causal variants in 
genome-wide association studies. 

# Source
https://github.com/biona001/Knockoffs.jl/blob/master/src/ghost.jl#L32
"""
function ghost_knockoffs(Zscores::AbstractVector{T}, D::AbstractMatrix{T}, 
    Σinv::AbstractMatrix{T}; m::Int = 1) where T
    p = size(D, 1)
    length(Zscores) == size(Σinv, 1) == size(Σinv, 2) == p || 
        error("Dimension mismatch")
    DΣinv = D * Σinv
    C = 2D - DΣinv * D
    v = sample_mvn_efficient(C, D, m) # Jiaqi's trick
    DΣinv .*= -1
    DΣinv += I
    Pz = DΣinv * Zscores
    return repeat(Pz, m) + v
end

"""
    MK_statistics(T0::Vector, Tk::Vector{Vector}; filter_method)

Computes the multiple knockoff statistics kappa, tau, and W. 

# Inputs
+ `T0`: p-vector of importance score for original variables
+ `Tk`: Vector storing T1, ..., Tm, where Ti is importance scores for 
    the `i`th knockoff copy
+ `filter_method`: Either `Statistics.median` (default) or max (original 
    function used in 2019 Gimenez and Zou)

# output
+ `κ`: Index of the most significant feature (`κ[i] = 0` if original feature most 
    important, otherwise `κ[i] = k` if the `k`th knockoff is most important)
+ `τ`: `τ[i]` stores the most significant statistic among original and knockoff
    variables minus `filter_method()` applied to the remaining statistics. 
+ `W`: coefficient difference statistic `W[i] = abs(T0[i]) - abs(Tk[i])`

# Source
https://github.com/biona001/Knockoffs.jl/blob/master/src/threshold.jl#L96
"""
function MK_statistics(
    T0::Vector{T}, 
    Tk::Vector{Vector{T}};
    filter_method::Function = Statistics.median
    ) where T
    p, m = length(T0), length(Tk)
    all(p .== length.(Tk)) || error("Length of T0 should equal all vectors in Tk")
    κ = zeros(Int, p) # index of largest importance score
    τ = zeros(p)      # difference between largest importance score and median of remaining
    W = zeros(p)      # importance score of each feature
    storage = zeros(m + 1)
    for i in 1:p
        storage[1] = abs(T0[i])
        for k in 1:m
            if abs(Tk[k][i]) > abs(T0[i])
                κ[i] = k
            end
            storage[k+1] = abs(Tk[k][i])
        end
        W[i] = (storage[1] - filter_method(@view(storage[2:end]))) * (κ[i] == 0)
        sort!(storage, rev=true)
        τ[i] = storage[1] - filter_method(@view(storage[2:end]))
    end
    return κ, τ, W
end

"""
    mk_threshold(τ::Vector{T}, κ::Vector{Int}, m::Int, q::Number)

Chooses the multiple knockoff threshold `τ̂ > 0` by setting
τ̂ = min{ t > 0 : (1/m + 1/m * {#j: κ[j] ≥ 1 and W[j] ≥ t}) / {#j: κ[j] == 0 and W[j] ≥ τ̂} ≤ q }.

# Inputs
+ `τ`: τ[i] stores the feature importance score for the ith feature, i.e. the value
    T0 - median(T1,...,Tm). Note in Gimenez and Zou, the max function is used 
    instead of median
+ `κ`: κ[i] stores which of m knockoffs has largest importance score. When original 
    variable has largest score, κ[i] == 0.
+ `m`: Number of knockoffs per variable generated
+ `q`: target FDR (between 0 and 1)
+ `rej_bounds`: Number of values of top τ to consider (default = 10000)

# Reference: 
+ Equations 8 and 9 in supplement of "Identification of putative causal loci in 
    wholegenome sequencing data via knockoff statistics" by He et al. 
+ Algorithm 1 of "Improving the Stability of the Knockoff Procedure: Multiple 
    Simultaneous Knockoffs and Entropy Maximization" by Gimenez and Zou.

# Source
https://github.com/biona001/Knockoffs.jl/blob/master/src/threshold.jl#L55
"""
function mk_threshold(τ::Vector{T}, κ::Vector{Int}, m::Int, q::Number,
    rej_bounds::Int=10000
    ) where T <: AbstractFloat
    0 ≤ q ≤ 1 || error("Target FDR should be between 0 and 1 but got $q")
    length(τ) == length(κ) || error("Length of τ and κ should be the same")
    p = length(τ) # number of features
    τ̂ = typemax(T)
    offset = 1 / m
    for (i, t) in enumerate(sort(τ, rev=true))
        numer_counter, denom_counter = 0, 0
        for i in 1:p
            κ[i] ≥ 1 && τ[i] ≥ t && (numer_counter += 1)
            κ[i] == 0 && τ[i] ≥ t && (denom_counter += 1)
        end
        ratio = (offset + offset * numer_counter) / max(1, denom_counter)
        ratio ≤ q && 0 < t < τ̂ && (τ̂ = t)
        i > rej_bounds && break
    end
    return τ̂
end

"""
    get_knockoff_qvalue(κ, τ, m, [groups], [rej_bounds])

Computes the knockoff q-value for each variable. The knockoff q-value is the 
minimum target FDR for a given variable to be selected. For details, see eq 19 
of https://www.nature.com/articles/s41467-022-34932-z and replace the 
knockoff-filter by the within-group knockoff filter proposed in Alg2 of "A
Powerful and Precise Feature-level Filter using Group Knockoffs" by Gu and He (2024).

Note: Code is directly translated from Zihuai's R code here:
https://github.com/biona001/ghostknockoff-gwas-reproducibility/blob/main/he_et_al/GKL_RunAnalysis_All.R#L36
"""
function get_knockoff_qvalue(κ::AbstractVector, τ::AbstractVector, m::Int;
    groups::AbstractVector=collect(1:length(τ)), rej_bounds::Int=10000
    )
    b = sortperm(τ, rev=true)
    c_0 = κ[b] .== 0
    offset = 1 / m
    # calculate ratios for top rej_bounds tau values
    ratio = Float64[]
    temp_0 = 0
    for i in eachindex(b)
        temp_0 = temp_0 + c_0[i]
        temp_1 = i - temp_0
        G_factor = maximum(values(countmap(groups[b][1:i])))
        temp_ratio = (offset*G_factor+offset*temp_1) / max(1, temp_0)
        push!(ratio, temp_ratio)
        i > rej_bounds && break
    end
    # calculate q values for top rej_bounds values
    qvalues = ones(length(τ))
    if any(x -> x > 0, τ)
        index_bound = maximum(findall(τ[b] .> 0))
        for i in eachindex(b)
            temp_index = i:min(length(b), rej_bounds, index_bound)
            length(temp_index) == 0 && continue
            qvalues[b[i]] = minimum(ratio[temp_index])*c_0[i]+1-c_0[i]
            i > rej_bounds && break
        end
        qvalues[qvalues .> 1] .= 1
    end
    return qvalues
end
