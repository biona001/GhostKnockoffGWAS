function gwas(
    bed::SnpArray,
    y::AbstractVector;
    # parameters for this function
    window_size::Int = 500, 
    snps_to_keep::Union{AbstractVector{Int}, Nothing}=nothing,
    # group knockoff options
    m = 5, 
    tol = 0.0001, 
    min_maf=0.01,
    min_hwe=0.0,
    force_block_diag=true,
    opt_method::String = "maxent",
    linkage::String = "average",
    force_contiguous::Bool=false,
    group_cor_cutoff::Float64=0.5,
    group_rep_cutoff::Float64=0.5,
    # parameters for lasso
    kappa_lasso::Float64 = 0.6,
    A_scaling_factor::Float64 = 0.01, 
    verbose=true,
    )
    # test code
    # bed = SnpArray("/oak/stanford/groups/zihuai/paisa/GLIMPSE2_maf0.01/chr22.bed")
    # y = randn(size(bed, 1))
    # window_size = 500
    # kappa_lasso = 0.6
    # m = 5
    # import GhostKnockoffGWAS.estimate_sigma
    # import GhostKnockoffGWAS.get_knockoff_qvalue
    # group_cor_cutoff = 0.5
    # group_rep_cutoff = 0.5
    # force_contiguous = false
    # linkage = :average
    # tol = 0.0001
    # verbose = false
    # opt_method = :maxent
    # A_scaling_factor = 0.01

    n, p = size(bed)
    windows = floor(Int, p / window_size)

    # compute Z scores. Note this does not load bed file into memory. 
    xla = SnpLinAlg{Float64}(bed, center=true, scale=true, impute=true)
    z = xla' * y ./ sqrt(n)

    # preallocated vectors
    Xi = Matrix{Float64}(undef, n, window_size)
    z_zko = Vector{Float64}(undef, (m+1) * window_size)
    beta = Float64[]  # full beta vector over all regions
    groups = String[] # group membership vector over all SNPs (original + knockoffs)

    # find lambda value for lasso
    lambdamax = maximum(abs, z) / sqrt(n)
    lambdamin = 0.0001lambdamax
    lambda_path = exp.(range(log(lambdamin), log(lambdamax), length=100)) |> reverse!
    lambda = kappa_lasso * maximum(abs, randn((m+1)*p)) / sqrt(n)
    lambda_path = vcat(lambda_path[findall(x -> x > lambda, lambda_path)], lambda)

    @showprogress for window in 1:windows
        # snps in current window
        idx = (window - 1) * window_size + 1:window * window_size

        # data in this window
        SnpArrays.copyto!(Xi, @view(bed[:, idx]), center=true, scale=true, 
            impute=true)
        zi = @view(z[idx])
        copyto!(z_zko, zi)

        # estimate sigma
        Sigma = estimate_sigma(Xi)

        # define groups and choose representatives
        group = hc_partition_groups(Symmetric(Sigma), cutoff=group_cor_cutoff, 
            linkage=linkage, force_contiguous=force_contiguous)
        group_reps = choose_group_reps(Symmetric(Sigma), group, 
            threshold=group_rep_cutoff)

        # group knockoff optimization
        S, D, obj = solve_s_graphical_group(Symmetric(Sigma), group, 
            group_reps, opt_method, m=m, tol=tol, verbose=verbose)

        # sample ghost knockoffs
        Sigma_inv = inv(Symmetric(Sigma))
        Zko = ghost_knockoffs(zi, D, Sigma_inv, m=m)
        z_zko[window_size+1:end] .= Zko

        # lasso
        Sigma .-= D
        for i in axes(D, 1)
            D[i, i] += A_scaling_factor
        end
        z_zko ./= sqrt(n)
        beta_i = block_group_ghostbasil(Sigma, D, z_zko, lambda_path, 
            m=m, delta_strong_size = 500, use_strong_rule=false)

        # save results
        append!(beta, beta_i)
        for k in 0:m
            append!(groups, ["window$(window)_$k" for g in group])
        end
    end

    # last window
    idx = windows * window_size + 1:p
    Xi = convert(Matrix{Float64}, @view(bed[:, idx]), center=true, scale=true, 
        impute=true)
    z_zko = z[idx]
    Sigma = estimate_sigma(Xi)
    group = hc_partition_groups(Symmetric(Sigma), cutoff=group_cor_cutoff, 
        linkage=linkage, force_contiguous=force_contiguous)
    group_reps = choose_group_reps(Symmetric(Sigma), group, 
        threshold=group_rep_cutoff)
    S, D, obj = solve_s_graphical_group(Symmetric(Sigma), group, 
        group_reps, opt_method, m=m, tol=tol, verbose=verbose)
    Sigma_inv = inv(Symmetric(Sigma))
    Zko = ghost_knockoffs(z_zko, D, Sigma_inv, m=m)
    append!(z_zko, Zko)
    Sigma .-= D
    for i in axes(D, 1)
        D[i, i] += A_scaling_factor
    end
    z_zko ./= sqrt(n)
    beta_i = block_group_ghostbasil(Sigma, D, z_zko, lambda_path, 
        m=m, delta_strong_size = 500, use_strong_rule=false)
    append!(beta, beta_i)
    for k in 0:m
        append!(groups, ["window$(windows + 1)_$k" for g in group])
    end

    # knockoff filter
    original_idx = findall(x -> endswith(x, "_0"), groups)
    T0 = beta[original_idx]
    Tk = [beta[findall(x -> endswith(x, "_$k"), groups)] for k in 1:m]
    kappa, tau, W = MK_statistics(T0, Tk)
    qvalues = get_knockoff_qvalue(kappa, tau, m, groups=groups)

    return qvalues
end
