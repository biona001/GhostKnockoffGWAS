"""
    ghostknockoffgwas(LD_files::String, z::Vector{Float64}, chr::Vector{Int}, 
        effect_allele::Vector{String}, non_effect_allele::Vector{String}, N::Int,
        hg_build::Int, outdir::String; [outname="result"], [seed=2023], 
        [target_chrs=sort!(unique(chr))], [A_scaling_factor = 0.01], 
        [kappa_lasso=0.6], [LD_shrinkage=false], 
        [target_fdrs=[0.01, 0.05, 0.1, 0.2]], [verbose=true], 
        [skip_shrinkage_check=false])

Runs the main `GhostKnockoffGWAS` pipeline on the Z scores in `z` using 
pre-computed knockoff data in `LD_files`. 

# Inputs
+ `LD_files`: Directory that stores pre-computed knockoff results
+ `z`: Vector of Z scores
+ `chr`: Chromosome of each Z score (cannot be X/Y/M chromosomes)
+ `pos`: Position of each Z score (specify hg build with `hg_build`)
+ `effect_allele`: Effect allele of Z score (i.e. ALT)
+ `non_effect_allele`: Non-effect allele of Z score (i.e. REF)
+ `N`: sample size of original study
+ `hg_build`: Human genome build (must be 19 or 38)
+ `outdir`: output directory, which must exist. We will output 2 files, one 
    containing the full analysis results, as well as a summary file. 

# Optional inputs
+ `seed`: Random seed for reproducibility (default = `2023`)
+ `target_chrs`: Target chromosomes to analyze. For example, one can specify
    `target_chrs = 22` to only analyze 1 chromosome, or `target_chrs = [1, 2]`
    to only analyze 2 chromosomes (default = `sort!(unique(chr))`).
+ `A_scaling_factor`: The scaling factor for `A = [X X̃]'*[X X̃]` for improving
    numerical stability. Scaling proceeds by adding `A_scaling_factor*I` to `A`
    (default = 0.01). 
+ `kappa_lasso`: A constant between 0 and 1 for tuning Lasso's lambda path. 
    Larger value forces earlier exit in Lasso lambda path, resulting in stronger 
    shrinkage. See the "lasso-min" method of "Controlled Variable Selection from
    Summary Statistics Only? A Solution via GhostKnockoffs and Penalized 
    Regression" by Chen et al (default `0.6`).
+ `LD_shrinkage`: Whether to perform shrinkage to LD and S matrices following
    method in SuSiE paper (i.e. eq 24 of "Fine-mapping from summary data with 
    the “Sum of Single Effects” model" by Zou et al). If `false`, we will still
    compute the shrinkage level, but it will not be used to adjust the LD
    matrices (default `false`). 
+ `target_fdrs`: Default target FDR levels (default = [0.01, 0.05, 0.1, 0.2])
+ `verbose`: Whether to print progress and informative intermediate results (
    default = `true`)
+ `skip_shrinkage_check`: Forces a result output even if there is a high
    estimated LD shrinkage by SuSiE's method (default = `false`)
+ `random_shuffle`: Whether to randomly permute the order of the original and 
    knockoff variables (default `false`). The main purpose of this option is 
    to take care of potential ordering bias of Lasso solvers. However, in our 
    simulations we never observed such biases, so we turn this off by default. 

# Output
By default we output 2 files into `outdir`
+ A knockoff statistics file where each SNP occupies a row and the columns include 
    various SNP attributes include rsid, AF, chr, pos, zscores...etc. The 
    columns `selected_fdr_FDR` indicates whether the variant was ultimately
    selected under the false discovery rate threshold of `FDR`.
+ A summary statistics file. The first dozens of rows print, for each false 
    discovery rate threshold `FDR`, the knockoff threshold `τ̂` and the number of
    groups that pass this threshold. The next couple of lines print some 
    parameters used in the knockoff analysis, as well as some timing results. 
"""
function ghostknockoffgwas(
    LD_files::String,
    z::Vector{Float64},
    chr::Vector{Int},
    pos::Vector{Int},
    effect_allele::Vector{String}, 
    non_effect_allele::Vector{String}, 
    N::Int, 
    hg_build::Int,
    outdir::String;
    outname::String="result",
    seed::Int = 2023,
    target_chrs=sort!(unique(chr)),
    A_scaling_factor = 0.01,
    kappa_lasso::Number=0.6,
    LD_shrinkage::Bool=false,
    target_fdrs = [0.01, 0.05, 0.1, 0.2],
    verbose::Bool=true,
    skip_shrinkage_check::Bool=false,
    random_shuffle::Bool = false
    )
    # check for errors
    any(isnan, z) && error("Z score contains NaN!")
    any(isinf, z) && error("Z score contains Inf!")
    isdir(outdir) || error("outdir $outdir does not exist!")
    N > 0 || error("Sample size N should be >0")
    hg_build == 19 || hg_build == 38 || error("Human genome build must be 19 or 38")
    length(z) == length(chr) == length(pos) == 
        length(effect_allele) == length(non_effect_allele) ||
        error("Length of z, chr, pos, effect_allele, and non_effect_allele should be the same")

    # number of simultaneous knockoffs (this is used in computation of `LD_files` and should NOT be changed)
    m = 5

    # find lambda value for lasso
    Random.seed!(seed)
    lambdamax = maximum(abs, z) / sqrt(N)
    lambdamin = 0.0001lambdamax
    lambda_path = exp.(range(log(lambdamin), log(lambdamax), length=100)) |> reverse!
    nsnps, tregions = count_matchable_snps(LD_files, z, chr, pos, effect_allele, 
        non_effect_allele, hg_build, target_chrs) # ~400 seconds on typed SNPs
    lambda = kappa_lasso * maximum(abs, randn((m+1)*nsnps)) / sqrt(N)
    lambda_path = vcat(lambda_path[findall(x -> x > lambda, lambda_path)], lambda)

    # intermediate vectors
    beta = Float64[]                       # full beta vector over all regions
    groups = String[]                      # group membership vector over all SNPs (original + knockoffs)
    Zscores = Float64[]                    # Z scores (original + knockoffs) for SNPs that can be matched to LD panel
    Zscores_store = Float64[]
    t1, t2, t3 = 0.0, 0.0, 0.0             # some timers
    t21, t22, t23, t24 = 0.0, 0.0, 0.0, 0.0# some timers
    start_t = time()
    γ_mean = 0.0                           # keeps track of LD shrinkage across regions
    headers = ["rsid", "AF", "chr", "ref", "alt", "pos_hg$(hg_build)"]
    df = DataFrame([String[], Float64[], Int[], String[], String[], Int[]], headers)

    # assemble knockoff results across regions
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

            # read knockoff results in current region
            t1 += @elapsed begin
                h5reader = HDF5.h5open(joinpath(LD_files, "chr$c", f), "r")
                Sigma_info = CSV.read(joinpath(LD_files, "chr$(c)", "Info_$fname.csv"), DataFrame)
                nknockoff_snps += size(Sigma_info, 1)
                # map reference LD panel to GWAS Z-scores by position
                LD_pos = Sigma_info[!, "pos_hg$(hg_build)"]
                rep_pos = LD_pos[read(h5reader, "group_reps")]
                shared_pos = intersect(LD_pos, GWAS_pos)
                # delete SNPs if ref/alt don't match
                pos2remove = Int[]
                for pos in shared_pos
                    GWAS_idx = findfirst(x -> x == pos, GWAS_pos)
                    LD_idx = findfirst(x -> x == pos, LD_pos)
                    ref_match_ea = Sigma_info[LD_idx, "ref"] == GWAS_ea[GWAS_idx]
                    alt_match_nea = Sigma_info[LD_idx, "alt"] == GWAS_nea[GWAS_idx]
                    ref_match_nea = Sigma_info[LD_idx, "ref"] == GWAS_nea[GWAS_idx]
                    alt_match_ea = Sigma_info[LD_idx, "alt"] == GWAS_ea[GWAS_idx]
                    if ref_match_nea && alt_match_ea
                        continue
                    elseif ref_match_ea && alt_match_nea
                        zscores[GWAS_idx] *= -1
                    else # SNP cannot get matched to LD panel
                        push!(pos2remove, pos)
                    end
                end
                setdiff!(shared_pos, unique!(pos2remove))
                setdiff!(rep_pos, pos2remove)
                # save matching snps info
                LD_keep_idx = indexin(shared_pos, LD_pos)
                GWAS_keep_idx = indexin(shared_pos, GWAS_pos)
                append!(df, @view(Sigma_info[LD_keep_idx, headers]))
                if length(LD_keep_idx) == 0 || length(GWAS_keep_idx) == 0 
                    nregions += 1
                    if verbose
                        println("region $nregions / $tregions (f = $f) matched 0 SNPs!")
                        flush(stdout)
                    end
                    continue
                end
            end

            # generate knockoffs for z scores
            t2 += @elapsed begin
                # use original Sigma and D (S matrix for rep+nonrep variables)
                Si = read(h5reader, "D")[LD_keep_idx, LD_keep_idx]
                Σi = read(h5reader, "Sigma")[LD_keep_idx, LD_keep_idx]
                zi = @view(zscores[GWAS_keep_idx])

                # shrinkage for consistency
                t21 += @elapsed begin
                    if length(zi) > 1000 # use reps to compute shrinkage
                        shared_rep_pos = intersect(shared_pos, rep_pos)
                        Σi_idx = filter!(!isnothing, indexin(shared_rep_pos, rep_pos))
                        zi_idx = filter!(!isnothing, indexin(shared_rep_pos, shared_pos))
                        Σi_rep = read(h5reader, "Sigma_reps")[Σi_idx, Σi_idx]
                        zi_rep = zi[zi_idx]
                        γ = find_optimal_shrinkage(Σi_rep, zi_rep)
                    else
                        γ = find_optimal_shrinkage(Σi, zi)
                    end
                    γ_mean += γ
                    if LD_shrinkage
                        Σi = (1 - γ)*Σi + γ*I
                        Si = (1 - γ)*Si + (m+1)/m*γ*I
                    end
                end

                # sample ghost knockoffs knockoffs
                Random.seed!(seed)
                t23 += @elapsed Σi_inv = inv(Symmetric(Σi))
                t24 += @elapsed Zko = ghost_knockoffs(zi, Si, Σi_inv, m=m)
            end

            # save mapped Zscores and some other variables
            empty!(Zscores_store)
            append!(Zscores_store, zi)         # append original Z scores 
            append!(Zscores_store, Zko)        # append knockoffs
            append!(Zscores, Zscores_store)
            current_groups = read(h5reader, "groups")[LD_keep_idx]
            grp = ["chr$(c)_$(fname)_group$(g)_0" for g in current_groups] # the last _0 implies this is original variable
            append!(groups, grp)
            for k in 1:m
                append!(groups, ["chr$(c)_$(fname)_group$(g)_$k" for g in current_groups])
            end
            length(Zscores_store) == (m+1) * length(current_groups) || 
                error("Number of Zscores should match groups")

            # randomly permute order of Z and Zko to avoid ordering bias
            if random_shuffle
                p = length(zi)
                perms = [collect(1:m+1) for _ in 1:p]
                for i in eachindex(zi)
                    shuffle!(perms[i])
                    @views permute!(Zscores_store[i:p:end], perms[i])
                end
            end

            # run GhostBasil
            t3 += @elapsed begin
                Ci = Σi - Si
                Si_scaled = Si + A_scaling_factor*I
                r = Zscores_store ./ sqrt(N)
                beta_i = block_group_ghostbasil(Ci, Si_scaled, r, lambda_path, 
                    m=m, delta_strong_size = 500, use_strong_rule=false)
            end

            # undo shuffling of Z and Zko
            if random_shuffle
                for i in eachindex(zi)
                    @views invpermute!(beta_i[i:p:end], perms[i])
                    @views invpermute!(Zscores_store[i:p:end], perms[i])
                end
            end

            # update counters
            append!(beta, beta_i)
            nsnps += length(shared_pos)
            nregions += 1
            if verbose
                println("region $nregions / $tregions (f = $f): chr $c, nz beta = " *
                        "$(count(!iszero, beta_i)), nsnps = $(length(shared_pos))" * 
                        ", shrinkage = $(round(γ, digits=4))")
                flush(stdout)
            end
        end
    end

    # some checks
    verbose && println("Matched $nsnps SNPs with Z-scores to the reference panel")
    length(Zscores) == length(groups) || 
        error("Number of Zscores should match groups")
    γ_mean /= nregions
    if γ_mean ∈ [0.1, 0.25]
        verbose && @warn(
            "Mean LD shrinkage is $γ_mean, which is a bit high. " * 
            "This suggests the LD panel supplied in \"$LD_files\" does not " *
            "well capture the correlation structure of the original study " * 
            "genotypes. Consider using a different LD panel."
        )
    elseif γ_mean > 0.25
        skip_shrinkage_check || error(
            "Mean LD shrinkage is $γ_mean, which is too high. " *
            "This suggests the LD panel supplied in \"$LD_files\" does not " *
            "correctly capture the correlation structure of the original study " * 
            "genotypes. Consider using a different LD panel. To bypass this " * 
            "error, use the option skip_shrinkage_check."
        )
    else
        verbose && println("Mean LD shrinkage = $γ_mean.")
    end

    # knockoff filter
    t4 = @elapsed begin
        original_idx = findall(x -> endswith(x, "_0"), groups)
        T0 = beta[original_idx]
        Tk = [beta[findall(x -> endswith(x, "_$k"), groups)] for k in 1:m]
        kappa, tau, W = MK_statistics(T0, Tk)
        qvalues = get_knockoff_qvalue(kappa, tau, m, groups=groups)

        # append result
        df[!, :group] = groups[original_idx]
        df[!, :zscores] = Zscores[original_idx]
        df[!, :lasso_beta] = beta[original_idx]
        df[!, :kappa] = kappa
        df[!, :tau] = tau
        df[!, :W] = W
        df[!, :qvals] = qvalues
        df[!, :pvals] = zscore2pval(df[!, :zscores])

        # pre-compute whether a variant is selected for easier interpretation
        for fdr in target_fdrs
            selected = zeros(Int, size(df, 1))
            selected[findall(x -> x ≤ fdr, qvalues)] .= 1
            df[!, "selected_fdr$fdr"] = selected
        end

        # sort by chr and pos
        sort!(df, [:chr, Symbol("pos_hg$(hg_build)")])

        # write to output
        CSV.write(joinpath(outdir, outname * ".txt"), df)
    end

    # save summary info
    open(joinpath(outdir, outname * "_summary.txt"), "w") do io
        for fdr in target_fdrs
            println(io, "target_fdr_$(fdr)_num_selected,", count(x -> x ≤ fdr, qvalues))
        end
        println(io, "m,$m")
        println(io, "nregions,$nregions")
        println(io, "nsnps,$nsnps")
        println(io, "lasso_lambda,$lambda")
        println(io, "mean_LD_shrinkage,$γ_mean")
        println(io, "import_time,$t1")
        println(io, "sample_knockoff_time,$t2")
        println(io, "ghostbasil_time,$t3")
        println(io, "knockoff_filter_time,$t4")
        println(io, "total_time,", time() - start_t)
        println(io, "sample_knockoff_time_t21,$t21")
        println(io, "sample_knockoff_time_t22,$t22")
        println(io, "sample_knockoff_time_t23,$t23")
        println(io, "sample_knockoff_time_t24,$t24")
    end

    return nothing
end
