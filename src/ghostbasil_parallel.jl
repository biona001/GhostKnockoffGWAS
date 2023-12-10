function ghostbasil_parallel(
    knockoff_dir::String,  # Directory that stores knockoff results (i.e. output from part 1)
    z::Vector{Float64},    # Z scores
    chr::Vector{Int},      # chromosome of each Z score
    pos::Vector{Int},      # position of each Z score (specify hg build with hg_build)
    effect_allele::Vector{String},       # effect allele of Z score
    non_effect_allele::Vector{String},   # non-effect allele of Z score
    N::Int,                # effective sample size
    hg_build::Int,         # human genome build, must be 19 or 38
    outdir::String;        # output directory (must exist)
    m::Int=5,
    outname::String="result",
    ncores = 1,
    seed::Int = 2023,
    target_chrs=1:22,
    A_scaling_factor = 0.01,
    kappa::Number=0.6,     # for tuning lambda
    save_intermdiate_result::Bool=false, # if true, will save beta, group, Zscore, and SNP summary stats, and not run knockoff filter
    LD_shrinkage::Bool=true, # if true, we will try to perform shrinkage to LD matrix following method in susie
    target_fdrs = [0.01, 0.05, 0.1, 0.15, 0.2],
    )
    # check for errors
    any(isnan, z) && error("Z score contains NaN!")
    any(isinf, z) && error("Z score contains Inf!")
    isdir(outdir) || error("outdir $outdir does not exist!")
    N > 0 || error("Effective sample size N should be >0")
    hg_build == 19 || hg_build == 38 || error("Human genome build must be 19 or 38")
    length(z) == length(chr) == length(pos) == 
        length(effect_allele) == length(non_effect_allele) ||
        error("Length of z, chr, pos, effect_allele, and non_effect_allele should be the same")

    # find lambda value for lasso
    Random.seed!(seed)
    lambdamax = maximum(abs, z) / sqrt(N)
    lambdamin = 0.0001lambdamax
    lambda_path = exp.(range(log(lambdamin), log(lambdamax), length=100)) |> reverse!
    nsnps = count_matchable_snps(knockoff_dir, z, chr, pos, effect_allele, 
        non_effect_allele, hg_build, target_chrs) # ~400 seconds on typed SNPs
    lambda = kappa * maximum(abs, randn((m+1)*nsnps)) / sqrt(N)
    lambda_path = lambda_path[findall(x -> x > lambda, lambda_path)]

    # intermediate vectors
    beta = Float64[]                       # full beta vector over all regions
    groups = String[]                      # group membership vector over all SNPs (original + knockoffs)
    groups_original = Int[]                # integer group membership vector for original SNPs
    Zscores = Float64[]                    # Z scores (original + knockoffs) for SNPs that can be matched to LD panel
    Zscores_ko_train = Float64[]           # needed for pseudo-validation in ghostbasil
    Zscores_store = Float64[]
    t1, t2, t3 = 0.0, 0.0, 0.0             # some timers
    t21, t22, t23, t24 = 0.0, 0.0, 0.0, 0.0# some timers
    start_t = time()
    γ_mean = 0.0                           # keeps track of LD shrinkage across regions
    df = DataFrame(rsid=String[], AF=Float64[], chr=Int[], 
        ref=String[], alt=String[], pos_hg19=Int[], pos_hg38=Int[])

    # assemble knockoff results across regions
    nregions, nsnps = 0, 0
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

            # read knockoff results in current region
            t1 += @elapsed begin
                result = JLD2.load(joinpath(knockoff_dir, "chr$c", f))
                Sigma_info = CSV.read(joinpath(knockoff_dir, "chr$(c)", "Info_$fname.csv"), DataFrame)
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
                        zscores[GWAS_idx] *= -1
                    else # SNP cannot get matched to LD panel
                        push!(remove_idx, i)
                    end
                end
                deleteat!(shared_snps, unique!(remove_idx))
                # save matching snps info
                LD_keep_idx = indexin(shared_snps, LD_pos)
                GWAS_keep_idx = indexin(shared_snps, GWAS_pos)
                append!(df, Sigma_info[LD_keep_idx, :])
                if length(LD_keep_idx) == 0 || length(GWAS_keep_idx) == 0 
                    continue
                end
            end

            # generate knockoffs for z scores
            t2 += @elapsed begin
                # use original Sigma and D (S matrix for rep+nonrep variables)
                Si = result["D"][LD_keep_idx, LD_keep_idx]
                Σi = result["Sigma"][LD_keep_idx, LD_keep_idx]
                zscore_tmp = @view(zscores[GWAS_keep_idx])
                t21 += @elapsed if LD_shrinkage
                    γ = find_optimal_shrinkage(Σi, zscore_tmp)
                    γ_mean += γ
                    Σi = (1 - γ)*Σi + γ*I
                    Si = (1 - γ)*Si + (m+1)/m*γ*I
                end
                # sample ghost knockoffs knockoffs
                Random.seed!(seed)
                t22 += @elapsed Zko_train = Knockoffs.sample_mvn_efficient(Σi, Si, m + 1)
                t23 += @elapsed Σi_inv = inv(Symmetric(Σi))
                t24 += @elapsed Zko = ghost_knockoffs(zscore_tmp, Si, Σi_inv, m=m)
            end

            # save mapped Zscores and some other variables
            empty!(Zscores_store)
            empty!(Zscores_ko_train)
            append!(Zscores_store, zscore_tmp) # append original Z scores 
            append!(Zscores_store, Zko)        # append knockoffs
            append!(Zscores, Zscores_store)
            append!(Zscores_ko_train, Zko_train)
            current_groups = result["groups"][LD_keep_idx]
            grp = ["chr$(c)_$(fname)_group$(g)_0" for g in current_groups] # the last _0 implies this is original variable
            append!(groups, grp)
            for k in 1:m
                append!(groups, ["chr$(c)_$(fname)_group$(g)_$k" for g in current_groups])
            end
            length(Zscores_store) == (m+1) * length(current_groups) == 
                length(Zscores_ko_train) || 
                error("Number of Zscores should match groups")

            # randomly permute order of Z and Zko to avoid ordering bias
            p = length(zscore_tmp)
            perms = [collect(1:m+1) for _ in 1:p]
            for i in eachindex(zscore_tmp)
                shuffle!(perms[i])
                @views permute!(Zscores_store[i:p:end], perms[i])
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
            for i in eachindex(zscore_tmp)
                @views invpermute!(beta_i[i:p:end], perms[i])
                @views invpermute!(Zscores_store[i:p:end], perms[i])
            end

            # update counters
            append!(beta, beta_i)
            nsnps += length(shared_snps)
            nregions += 1
            println("region $nregions: chr $c, nz beta = $(count(!iszero, beta_i)), nsnps = $(length(shared_snps))")
            flush(stdout)
        end
    end

    # some checks
    nregions == 1703 || @warn("Number of successfully solved region is $nregions, expected 1703")
    println("Matched $nsnps SNPs with Z-scores to the reference panel")
    length(Zscores) == length(groups) || 
        error("Number of Zscores should match groups")

    if save_intermdiate_result
        writedlm(joinpath(outdir, outname * "_groups.txt"), groups)
        writedlm(joinpath(outdir, outname * "_betas.txt"), beta)
        writedlm(joinpath(outdir, outname * "_Zscores.txt"), Zscores)
        CSV.write(joinpath(outdir, outname * "_stats.csv"), df)
    end

    # knockoff filter immediately
    t4 = @elapsed begin
        original_idx = findall(x -> endswith(x, "_0"), groups)
        T0 = beta[original_idx]
        Tk = [beta[findall(x -> endswith(x, "_$k"), groups)] for k in 1:m]
        T_group0 = Float64[]
        T_groupk = [Float64[] for k in 1:m]
        groups_original = groups[findall(x -> endswith(x, "_0"), groups)]
        unique_groups = unique(groups_original)
        for idx in find_matching_indices(unique_groups, groups_original)
            push!(T_group0, sum(abs, @view(T0[idx])))
            for k in 1:m
                push!(T_groupk[k], sum(abs, @view(Tk[k][idx])))
            end
        end
        kappa, tau, W = Knockoffs.MK_statistics(T_group0, T_groupk)
        
        # save analysis result
        df[!, :group] = groups[original_idx]
        df[!, :zscores] = Zscores[original_idx]
        df[!, :lasso_beta] = beta[original_idx]
        W_full, kappa_full, tau_full = Float64[], Float64[], Float64[]
        for idx in indexin(groups_original, unique_groups)
            push!(W_full, W[idx])
            push!(kappa_full, kappa[idx])
            push!(tau_full, tau[idx])
        end
        df[!, :W] = W_full
        df[!, :kappa] = kappa_full
        df[!, :tau] = tau_full
        df[!, :pvals] = zscore2pval(df[!, :zscores])
        for fdr in target_fdrs
            q = mk_threshold(tau, kappa, m, fdr)
            selected = zeros(Int, size(df, 1))
            selected[findall(x -> x ≥ q, W)] .= 1
            df[!, "selected_fdr$fdr"] = selected
        end
        CSV.write(joinpath(outdir, outname * ".txt"), df)
    end

    # save summary info
    open(joinpath(outdir, outname * "_summary.txt"), "w") do io
        # Q-value (i.e. threshold for knockoff filter)
        for fdr in target_fdrs
            q = mk_threshold(tau, kappa, m, fdr)
            println(io, "target_fdr_$(fdr),$q")
            println(io, "target_fdr_$(fdr)_num_selected,", count(x -> x ≥ q, W))
        end
        println(io, "m,$m")
        println(io, "nregions,$nregions")
        println(io, "nsnps,$nsnps")
        println(io, "lasso_lambda,$lambda")
        println(io, "mean_LD_shrinkage,$(γ_mean / nregions)")
        println(io, "import_time,$t1")
        println(io, "sample_knockoff_time,$t2")
        println(io, "ghostbasil_time,$t3")
        println(io, "knockoff_filter_time,$t4")
        println(io, "total_time,", time() - start_t)
        println(io, "sample_knockoff_time_t21=$t21")
        println(io, "sample_knockoff_time_t22=$t22")
        println(io, "sample_knockoff_time_t23=$t23")
        println(io, "sample_knockoff_time_t24=$t24")
    end

    return nothing
end
