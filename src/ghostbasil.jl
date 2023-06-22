function ghostbasil(
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
    ncores = Threads.nthreads(),
    target_chrs=1:22,
    A_scaling_factor = 0.01
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
    Zscores_ko_train = Float64[]           # needed for pseudo-validation in ghostbasil
    flip_sign_idx = Int[]
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
                empty!(flip_sign_idx)
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
                        push!(remove_idx, i)
                        # zscores[GWAS_idx] *= -1
                        # push!(flip_sign_idx, LD_idx)
                    else # SNP cannot get matched to LD panel
                        push!(remove_idx, i)
                    end
                end
                deleteat!(shared_snps, unique!(remove_idx))
                # also flip the signs of Σ/S/etc if ref/alt were opposite
                # note: sign in result["S"] not flipped (only on reps)
                for var in ["Sigma", "Sigma2", "D", "SigmaInv", "SigmaInv2"]
                    mat = result[var]
                    for idx in flip_sign_idx
                        mat[:, idx] .*= -1
                        mat[idx, :] .*= -1
                    end
                end
                # save matching snps info
                LD_keep_idx = indexin(shared_snps, LD_pos)
                GWAS_keep_idx = indexin(shared_snps, GWAS_pos)
                append!(df, Sigma_info[LD_keep_idx, :])
                if length(LD_keep_idx) == 0 || length(GWAS_keep_idx) == 0 
                    continue
                end
                # println("matched $(length(shared_snps)) snps")
            end

            # generate knockoffs for z scores
            # todo: use 2 stage procedure for efficiency
            t2 += @elapsed begin
                # use original Sigma and D for pseudo-validation
                D = result["D"][LD_keep_idx, LD_keep_idx]
                Σ = result["Sigma"][LD_keep_idx, LD_keep_idx]
                Zko_train = Knockoffs.sample_mvn_efficient(Σ, D, m + 1)
                # sample ghost knockoffs knockoffs
                Σinv = inv(Symmetric(Σ))
                Zko = ghost_knockoffs(zscores[GWAS_keep_idx], D, Σinv, m=m)
            end

            # save mapped Zscores and some other variables
            push!(Sigma, Σ)
            push!(S, D)
            append!(Zscores, zscores[GWAS_keep_idx])
            append!(Zscores, Zko)
            append!(Zscores_ko_train, Zko_train)
            current_groups = result["groups"][LD_keep_idx]
            grp = ["chr$(c)_$(fname)_group$(g)_0" for g in current_groups] # the last _0 implies this is original variable
            append!(groups, grp)
            for k in 1:m
                append!(groups, ["chr$(c)_$(fname)_group$(g)_$k" for g in current_groups])
            end

            # update counters
            nsnps += length(shared_snps)
            nregions += 1
            println("region $nregions: nsnps = $nsnps")
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
        result <- ghostbasil(A, rt, delta.strong.size = 500, max.strong.size = nrow(A), n.threads=ncores, use.strong.rule=F)
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
        fit.basil<-ghostbasil(A, r=r,user.lambdas=lambda.seq, delta.strong.size = 500, max.strong.size = nrow(A), use.strong.rule=F)
        beta<-fit.basil$betas[,ncol(fit.basil$betas)]
        ll <- fit.basil$lmdas
        """
        @rget beta ll lambdas
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

    # knockoff filter
    t4 = @elapsed begin
        # todo: abstract away this logic into Knockoffs.jl
        original_idx = findall(x -> endswith(x, "_0"), groups)
        T0 = beta[original_idx]
        Tk = [beta[findall(x -> endswith(x, "_$k"), groups)] for k in 1:m]
        T_group0 = Float64[]
        T_groupk = [Float64[] for k in 1:m]
        groups_original = groups[findall(x -> endswith(x, "_0"), groups)]
        unique_groups = unique(groups_original)
        for g in unique_groups
            idx = findall(x -> x == g, groups_original)
            push!(T_group0, sum(abs, @view(T0[idx])))# / length(idx))
            for k in 1:m
                push!(T_groupk[k], sum(abs, @view(Tk[k][idx])))# / length(idx))
            end
        end
        kappa, tau, W = Knockoffs.MK_statistics(T_group0, T_groupk)

        # save analysis result
        df[!, :group] = groups[original_idx]
        df[!, :zscores] = Zscores[original_idx]
        df[!, :lasso_beta] = beta[original_idx]
        W_full, kappa_full, tau_full = Float64[], Float64[], Float64[]
        for g in groups_original
            idx = findfirst(x -> x == g, unique_groups)
            push!(W_full, W[idx])
            push!(kappa_full, kappa[idx])
            push!(tau_full, tau[idx])
        end
        df[!, :W] = W_full
        df[!, :kappa] = kappa_full
        df[!, :tau] = tau_full
        df[!, :pvals] = 2ccdf.(Normal(), abs.(df[!, :zscores])) # convert zscores to marginal p-values
        CSV.write(joinpath(outdir, outname * ".txt"), df)
    end

    # save summary info
    open(joinpath(outdir, outname * "_summary.txt"), "w") do io
        # Q-value (i.e. threshold for knockoff filter)
        for fdr in [0.01, 0.05, 0.1, 0.15, 0.2]
            q = mk_threshold(tau, kappa, m, fdr)
            println(io, "target_fdr_$(fdr),$q")
            println(io, "target_fdr_$(fdr)_num_selected,", count(x -> x ≥ q, W))
        end
        println(io, "m,$m")
        println(io, "nregions,$nregions")
        println(io, "nsnps,$nsnps")
        println(io, "import_time,$t1")
        println(io, "sample_knockoff_time,$t2")
        println(io, "ghostbasil_time,$t3")
        println(io, "knockoff_filter_time,$t4")
        println(io, "total_time,", time() - start_t)
        println(io, "significant_marginal_pvals,", count(x -> x < 5e-8, df[!, :pvals]))
    end

    # also save full lambda sequence
    writedlm(joinpath(outdir, outname * "_lambdas.txt"), lambdas)

    q = mk_threshold(tau, kappa, m, 0.1)
    selected_idx = findall(x -> x ≥ q, W)
    selected_groups = unique_groups[selected_idx]
    selected_snps = findall(x -> x in selected_groups, groups_original)
    return selected_groups, selected_snps
end

function ghostbasil_parallel_by_block(
    knockoff_dir::String,  # Directory that stores knockoff results (i.e. output from part 1)
    z::Vector{Float64},    # Z scores
    chr::Vector{Int},      # chromosome of each Z score
    pos::Vector{Int},      # position of each Z score (tentatively assuming hg19 coordinates)
    ref::Vector{String},   # ref allele of Z score
    alt::Vector{String},   # alt allele of Z score
    N::Int,                # effective sample size
    hg_build::Int,         # human genome build, must be 19 or 38
    outdir::String;        # output directory (must exist)
    m::Int=5,
    outname::String="result",
    ncores = Threads.nthreads(),
    target_chrs=1:22,
    A_scaling_factor = 0.01
    )
    # check for errors
    any(isnan, z) && error("Z score contains NaN!")
    any(isinf, z) && error("Z score contains Inf!")
    isdir(outdir) || error("outdir $outdir does not exist!")
    N > 0 || error("Effective sample size N should be >0")
    hg_build == 19 || hg_build == 38 || error("Human genome build must be 19 or 38")
    length(z) == length(chr) == length(pos) == length(ref) == length(alt) ||
        error("Length of z, chr, pos, ref, and alt should be the same")

    # compute lambda sequence
    lambdamax = maximum(abs.(z)) ./ sqrt(N)
    lambdamin = 0.01lambdamax # default in ghostbasil
    lambda_path = range(lambdamin, lambdamax, length=100) |> collect |> reverse!
#     lambda_path = exp.(range(log(lambdamin), log(lambdamax), length=100)) |> reverse!

    # matrices to be passed into ghostbasil (~25GB to store these in memory)
    Sigma = Matrix{Float64}[]              # overall covariance matrix for all chrom
    S = Matrix{Float64}[]                  # overall block diagonal matrix for all chrom

    # intermediate vectors
    beta = Float64[]                       # full beta vector over all regions
    groups = String[]                      # group membership vector over all SNPs (original + knockoffs)
    groups_original = Int[]                # integer group membership vector for original SNPs
    Zscores = Float64[]                    # Z scores (original + knockoffs) for SNPs that can be matched to LD panel
    Zscores_ko_train = Float64[]           # needed for pseudo-validation in ghostbasil
    Zscores_store = Float64[]
    flip_sign_idx = Int[]
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
        GWAS_ref = ref[chr_idx]
        GWAS_alt = alt[chr_idx]
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
                empty!(flip_sign_idx)
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
                        push!(remove_idx, i)
                        # zscores[GWAS_idx] *= -1
                        # push!(flip_sign_idx, LD_idx)
                    else # SNP cannot get matched to LD panel
                        push!(remove_idx, i)
                    end
                end
                # also flip the signs of Σ/S/etc if ref/alt were opposite
                # note: sign in result["S"] not flipped (only on reps)
                for var in ["Sigma", "Sigma2", "D", "SigmaInv", "SigmaInv2"]
                    mat = result[var]
                    for idx in flip_sign_idx
                        mat[:, idx] .*= -1
                        mat[idx, :] .*= -1
                    end
                end
                # save matching snps info
                LD_keep_idx = indexin(shared_snps, LD_pos)
                GWAS_keep_idx = indexin(shared_snps, GWAS_pos)
                append!(df, Sigma_info[LD_keep_idx, :])
                if length(LD_keep_idx) == 0 || length(GWAS_keep_idx) == 0 
                    continue
                end
                println("matched $(length(shared_snps)) snps")
            end

            # generate knockoffs for z scores
            t2 += @elapsed begin
                # use original Sigma and D (S matrix for rep+nonrep variables) for pseudo-validation
                Si = result["D"][LD_keep_idx, LD_keep_idx]
                Σi = result["Sigma"][LD_keep_idx, LD_keep_idx]
                Zko_train = Knockoffs.sample_mvn_efficient(Σi, Si, m + 1)
                # sample ghost knockoffs knockoffs
                Σi_inv = inv(Symmetric(Σi))
                Zko = ghost_knockoffs(zscores[GWAS_keep_idx], Si, Σi_inv, m=m)
            end

            # save mapped Zscores and some other variables
            empty!(Zscores_store)
            empty!(Zscores_ko_train)
            append!(Zscores_store, zscores[GWAS_keep_idx])
            append!(Zscores_store, Zko)
            append!(Zscores, Zscores_store)
            append!(Zscores_ko_train, Zko_train)
            current_groups = result["groups"][LD_keep_idx]
            grp = ["chr$(c)_$(fname)_group$(g)_0" for g in current_groups] # the last _0 implies this is original variable
            append!(groups, grp)
            for k in 1:m
                append!(groups, ["chr$(c)_$(fname)_group$(g)_$k" for g in current_groups])
            end
            length(Zscores_store) == (m+1) * length(current_groups) || 
                error("Number of Zscores should match groups")

            # run GhostBasil
            t31 = @elapsed begin
                ntrain = 4N / 5
                nvalid = N / 5
                Ci = Σi - Si
                Si_scaled = Si + A_scaling_factor*I
                # r, r train, and r validate
                r = Zscores_store ./ sqrt(N)
                rt = r + sqrt(nvalid / (N*ntrain)) .* Zscores_ko_train
                rv = 1/nvalid * (N*r - ntrain*rt)
                @rput Ci Si Si_scaled r m ncores rt rv lambda_path
                R"""
                # make S into a block matrix with a single block
                Si_scaled <- BlockMatrix(list(Si_scaled))
                # Si_test <- Si$to_dense()
                # form Bi
                Bi <- BlockGroupGhostMatrix(Ci, Si_scaled, m+1)
                # Bi_test <- Bi$to_dense()
                # @rget Bi_test
                # p = length(shared_snps)
                # @test all(Bi_test[1:p, 1:p] .≈ Σi + 0.01I)
                # @test all(Bi_test[1:p, p+1:2p] .≈ Σi - Si)
                # @test all(Bi_test[p+1:2p, p+1:2p] .≈ Σi + 0.01I)
                # @test all(Bi_test[p+1:2p, 2p+1:3p] .≈ Σi - Si)
                # Basil max coordinate descents reached at lambda index: 2.
                # stochasic: sometimes throw max iteration warning, sometimes never converge
                result <- ghostbasil(Bi, rt, delta.strong.size=500, user.lambdas=lambda_path,
                    max.strong.size = nrow(Bi), n.threads=ncores, use.strong.rule=F)
                betas <- as.matrix(result$betas)
                """
            end

            t32 = @elapsed begin
                R"""
                # find best lambda
                Get.f<-function(beta,A,r){return(t(beta)%*%r/sqrt(A$quad_form(beta)))}
                f.lambda<-apply(betas,2,Get.f,A=Bi,r=rv)
                f.lambda[is.na(f.lambda)]<--Inf
                lambda<-lambda_path[which.max(f.lambda)]
                """
            end

            t33 = @elapsed begin
                R"""
                # refit ghostbasil
                lambda.seq <- lambda_path[lambda_path > lambda]
                lambda.seq <- c(lambda.seq, lambda)
                fit.basil<-ghostbasil(Bi, r=r,user.lambdas=lambda.seq, 
                    delta.strong.size=500, max.strong.size = nrow(Bi), use.strong.rule=F)
                beta_i <- fit.basil$betas[,ncol(fit.basil$betas)]
                R2 <- fit.basil$rsqs
                """
                @rget beta_i R2 betas
            end
            t3 += t31 + t32 + t33

            # update counters
            append!(beta, beta_i)
            nsnps += length(shared_snps)
            nregions += 1
#             nzcnt = 0
#             for β in eachcol(betas)
#                 count(!iszero, β) == 0 && (nzcnt += 1)
#             end
            println("region $nregions: nsnps = $nsnps, t1=$t1, t2=$t2, t3=$t3, t31 = $t31, t32 = $t32, t33 = $t33")
            flush(stdout)
        end
    end

    # some checks
    nregions == 1703 || @warn("Number of successfully solved region is $nregions, expected 1703")
    println("Matched $nsnps SNPs with Z-scores to the reference panel")
    length(Zscores) == length(groups) || 
        error("Number of Zscores should match groups")

    # knockoff filter
    t4 = @elapsed begin
        # todo: abstract away this logic into Knockoffs.jl
        original_idx = findall(x -> endswith(x, "_0"), groups)
        T0 = beta[original_idx]
        Tk = [beta[findall(x -> endswith(x, "_$k"), groups)] for k in 1:m]
        T_group0 = Float64[]
        T_groupk = [Float64[] for k in 1:m]
        groups_original = groups[findall(x -> endswith(x, "_0"), groups)]
        unique_groups = unique(groups_original)
        for g in unique_groups
            idx = findall(x -> x == g, groups_original)
            push!(T_group0, sum(abs, @view(T0[idx])))# / length(idx))
            for k in 1:m
                push!(T_groupk[k], sum(abs, @view(Tk[k][idx])))# / length(idx))
            end
        end
        kappa, tau, W = Knockoffs.MK_statistics(T_group0, T_groupk)
        
        # save analysis result
        df[!, :group] = groups[original_idx]
        df[!, :zscores] = Zscores[original_idx]
        df[!, :lasso_beta] = beta[original_idx]
        W_full, kappa_full, tau_full = Float64[], Float64[], Float64[]
        for g in groups_original
            idx = findfirst(x -> x == g, unique_groups)
            push!(W_full, W[idx])
            push!(kappa_full, kappa[idx])
            push!(tau_full, tau[idx])
        end
        df[!, :W] = W_full
        df[!, :kappa] = kappa_full
        df[!, :tau] = tau_full
        df[!, :pvals] = zscore2pval(df[!, :zscores])
        CSV.write(joinpath(outdir, outname * ".txt"), df)
    end

    # save summary info
    open(joinpath(outdir, outname * "_summary.txt"), "w") do io
        # Q-value (i.e. threshold for knockoff filter)
        for fdr in [0.01, 0.05, 0.1, 0.15, 0.2]
            q = mk_threshold(tau, kappa, m, fdr)
            println(io, "target_fdr_$(fdr),$q")
            println(io, "target_fdr_$(fdr)_num_selected,", count(x -> x ≥ q, W))
        end
        println(io, "m,$m")
        println(io, "nregions,$nregions")
        println(io, "nsnps,$nsnps")
        println(io, "import_time,$t1")
        println(io, "sample_knockoff_time,$t2")
        println(io, "ghostbasil_time,$t3")
        println(io, "knockoff_filter_time,$t4")
        println(io, "total_time,", time() - start_t)
    end

    # also save full lambda sequence
    writedlm(joinpath(outdir, outname * "_lambdas.txt"), lambdas)

#     q = mk_threshold(tau, kappa, m, 0.1)
#     selected_idx = findall(x -> x ≥ q, W)
#     selected_groups = unique_groups[selected_idx]
#     selected_snps = findall(x -> x in selected_groups, groups_original)
#     return selected_groups, selected_snps
    return nothing
end