function ghostknockoffgwas(
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
    kappa::Number=0.6,     # for tuning lambda, only used when pseudo_validate = false
    pseudo_validate::Bool = false, # if true, uses pseudo-validation, otherwise use zhaomeng's new technique
    LD_shrinkage::Bool=true, # if true, we will try to perform shrinkage to LD matrix following method in susie
    save_intermdiate_result::Bool=false, # if true, will save beta, group, Zscore, and SNP summary stats, and not run knockoff filter
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

    # matrices to be passed into ghostbasil (~25GB to store these in memory)
    Sigma = Matrix{Float64}[]              # overall covariance matrix for all chrom
    S = Matrix{Float64}[]                  # overall block diagonal matrix for all chrom

    # intermediate vectors
    groups  = String[]                     # group membership vector over all SNPs (original + knockoffs)
    groups_original = Int[]                # integer group membership vector for original SNPs
    Zscores = Float64[]                    # Z scores (original + knockoffs) for SNPs that can be matched to LD panel
    Zscores_ko_train = Float64[]           # needed for pseudo-validation in ghostbasil (Zscores_ko_train is a sample from N(0, A))
    γ_mean = 0.0                           # keeps track of LD shrinkage across regions
    permutations = Vector{Vector{Int64}}[] # needed to randomly shuffle Z and Zko in each block
    t1, t2, t3 = 0.0, 0.0, 0.0             # some timers
    start_t = time()
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

            t1 += @elapsed begin
                # read knockoff results
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
            # todo: use 2 stage procedure for efficiency
            t2 += @elapsed begin
                # use original Sigma and D for pseudo-validation
                D = result["D"][LD_keep_idx, LD_keep_idx]
                Σ = result["Sigma"][LD_keep_idx, LD_keep_idx]
                zscore_tmp = @view(zscores[GWAS_keep_idx])
                if LD_shrinkage
                    γ = find_optimal_shrinkage(Σ, zscore_tmp)
                    γ_mean += γ
                    Σ = (1 - γ)*Σ + γ*I
                    D = (1 - γ)*D + (m+1)/m*γ*I
                end
                # sample ghost knockoffs knockoffs
                Random.seed!(seed)
                Zko_train = Knockoffs.sample_mvn_efficient(Σ, D, m + 1)
                Σinv = inv(Symmetric(Σ))
                Zko = ghost_knockoffs(zscore_tmp, D, Σinv, m=m)
            end

            # save mapped Zscores and some other variables
            push!(Sigma, Σ)
            push!(S, D)
            append!(Zscores, zscore_tmp)
            append!(Zscores, Zko)
            append!(Zscores_ko_train, Zko_train)
            current_groups = result["groups"][LD_keep_idx]
            grp = ["chr$(c)_$(fname)_group$(g)_0" for g in current_groups] # the last _0 implies this is original variable
            append!(groups, grp)
            for k in 1:m
                append!(groups, ["chr$(c)_$(fname)_group$(g)_$k" for g in current_groups])
            end

            # randomly permute order of Z and Zko to avoid ordering bias
            p = length(zscore_tmp)
            perms = [collect(1:m+1) for _ in 1:p]
            for i in 1:p
                shuffle!(perms[i])
                Zi = @view(Zscores[(m+1)*nsnps+1:end])
                @views permute!(Zi[i:p:end], perms[i])
            end

            # update counters
            nsnps += length(shared_snps)
            nregions += 1
            push!(permutations, perms)
            println("region $nregions: nsnps = $nsnps, f = $f")
            flush(stdout)
        end
    end

    # some checks
    nregions == 1703 || @warn("Number of successfully solved region is $nregions, expected 1703")
    if nsnps == 0
        error("Matched $nsnps SNPs with Z-scores to the reference panel")
    else
        println("Matched $nsnps SNPs with Z-scores to the reference panel")
    end
    length(Zscores) == length(groups) || error("Number of Zscores should match groups")
    sum(size.(Sigma, 1)) == sum(size.(S, 1)) == size(df, 1) || 
        error("Dimension of Sigma, S, or df doesn't match")

    # run GhostBasil
    if pseudo_validate
        # note: lambda is chosen via pseudo-summary statistic approach
        # todo: wrap ghostbasil C++ code with Julia
        t3 = @elapsed begin
            ntrain = 4N / 5
            nvalid = N / 5
            # r, r train, and r validate
            r = Zscores ./ sqrt(N)
            rt = r + sqrt(nvalid / (N*ntrain)) .* Zscores_ko_train
            rv = 1/nvalid * (N*r - ntrain*rt)
            @rput Sigma S r m ncores rt rv A_scaling_factor
            R"""
            B <- c()
            for(i in 1:length(S)){
                Si <- as.matrix(S[[i]]) + A_scaling_factor*diag(1,nrow(S[[i]]))
                Si <- BlockMatrix(list(Si))
                Ci <- Sigma[[i]] - S[[i]]
                Bi <- BlockGroupGhostMatrix(Ci, Si, m+1)
                B  <- append(B, list(Bi))
            }
            A <- BlockBlockGroupGhostMatrix(B)
            result <- ghostbasil(A, rt, delta.strong.size = 500, 
                max.strong.size = nrow(A), n.threads=ncores, use.strong.rule=F)
            betas <- as.matrix(result$betas)
            lambdas <- result$lmdas
            
            # find best lambda
            Get.f<-function(beta,A,r){return(t(beta)%*%r/sqrt(A$quad_form(beta)))}
            f.lambda<-apply(betas,2,Get.f,A=A,r=rv)
            f.lambda[is.na(f.lambda)]<--Inf
            lambda<-lambdas[which.max(f.lambda)]
            
            # refit ghostbasil
            lambda.seq <- lambdas[lambdas > lambda]
            lambda.seq <- c(lambda.seq, lambda)
            fit.basil<-ghostbasil(A, r=r, user.lambdas=lambda.seq, n.threads=ncores, 
                delta.strong.size = 500, max.strong.size = nrow(A), 
                use.strong.rule=F)
            beta<-fit.basil$betas[,ncol(fit.basil$betas)]
            ll <- fit.basil$lmdas
            """
            @rget lambdas
            lambda = lambdas[end]
        end
    else
        t3 = @elapsed begin
            # create lambda sequence for ghostbasil lasso
            Random.seed!(seed)
            lambdamax = maximum(abs, z) / sqrt(N)
            lambdamin = 0.0001lambdamax
            lambda_path = exp.(range(log(lambdamin), log(lambdamax), length=100)) |> reverse! # default lambda_path in ghostbasil
            lambda = kappa * maximum(abs, randn((m+1)*nsnps)) / sqrt(N)
            lambda_path = lambda_path[findall(x -> x > lambda, lambda_path)]
            r = Zscores ./ sqrt(N)
            @rput Sigma S r m ncores A_scaling_factor lambda_path
            R"""
            B <- c()
            for(i in 1:length(S)){
                Si <- as.matrix(S[[i]]) + A_scaling_factor*diag(1,nrow(S[[i]]))
                Si <- BlockMatrix(list(Si))
                Ci <- Sigma[[i]] - S[[i]]
                Bi <- BlockGroupGhostMatrix(Ci, Si, m+1)
                B  <- append(B, list(Bi))
            }
            A <- BlockBlockGroupGhostMatrix(B)

            # fit ghostbasil
            fit.basil<-ghostbasil(A, r=r, user.lambdas=lambda_path, 
                delta.strong.size = 500, max.strong.size = nrow(A), 
                use.strong.rule=F, n.threads=ncores)
            beta <- fit.basil$betas[,ncol(fit.basil$betas)]
            ll <- fit.basil$lmdas
            """
        end
    end
    @rget beta

    # undo shuffling of Z and Zko
    counter = 0
    for perms in permutations
        p = length(perms)
        for i in eachindex(perms)
            beta_i = @view(beta[counter+1:counter+(m+1)*p])
            Zscores_i = @view(Zscores[counter+1:counter+(m+1)*p])
            @views invpermute!(beta_i[i:p:end], perms[i])
            @views invpermute!(Zscores_i[i:p:end], perms[i])
        end
        counter += (m+1)*p
    end

    # testing group structure
#     Sigma_full = BlockDiagonal(Sigma) |> Matrix
#     S_full = BlockDiagonal(S) |> Matrix
#     non_zero_idx = CartesianIndex[]
#     for g in unique(groups)
#         idx = findall(x -> x == g, groups)
#         for i in idx, j in idx
#             push!(non_zero_idx, CartesianIndex(i, j))
#         end
#     end
#     @test length(intersect(findall(!iszero, S_full), non_zero_idx)) == length(non_zero_idx)
#     [findall(!iszero, S_full) non_zero_idx]
    
    # testing A and Ci structure in ghostbasil
#     @rget Ci
#     @test all(Ci .≈ Sigma[end] - S[end])    

    if save_intermdiate_result
        # save group and beta information and apply knockoff filter later
        writedlm(joinpath(outdir, outname * "_groups.txt"), groups)
        writedlm(joinpath(outdir, outname * "_betas.txt"), beta)
        writedlm(joinpath(outdir, outname * "_Zscores.txt"), Zscores)
        CSV.write(joinpath(outdir, outname * "_stats.csv"), df)
    end

    # knockoff filter
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
        df[!, :pvals] = 2ccdf.(Normal(), abs.(df[!, :zscores])) # convert zscores to marginal p-values
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
        println(io, "significant_marginal_pvals,", count(x -> x < 5e-8, df[!, :pvals]))
    end

    return nothing
end
