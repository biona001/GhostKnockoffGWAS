function estimate_sigma(X::AbstractMatrix; min_eigval=1e-5)
    # n, p = size(X)
    # if n > p
    #     Sigma = cor(X)
    # else
    #     error("p > n case not handled yet!")
    # end

    # shrinkage based covariance, handles high dimensional case
    Sigma = cov(LinearShrinkage(DiagonalUnequalVariance(), :lw), X)

    # ensure Sigma is PSD
    evals, evecs = eigen(Sigma)
    evals[findall(x -> x < min_eigval, evals)] .= min_eigval
    Sigma = evecs * Diagonal(evals) * evecs'
    Statistics.cov2cor!(Sigma, sqrt.(diag(Sigma))) # scale to correlation matrix

    return Sigma
end

"""
    get_block(vcffile::String, chr::Int, start_bp::Int, end_bp::Int;
        [min_maf::Float64=0.01], [min_hwe::Float64=0.0], 
        [snps_to_keep::Union{AbstractVector{Int}, Nothing}=nothing])

Imports genotype data of a VCF file from pos `start_bp` to `end_bp` as 
double precision numeric matrix. Each row is a sample and each column is a SNP.
All entries are {0, 1, 2} except for missing entries which is imputed with 
column mean.

# Inputs
+ `vcffile`: A VCF file storing individual level genotypes. Must end in `.vcf` 
    or `.vcf.gz`. The ALT field for each record must be unique, i.e. 
    multiallelic records must be split first. Missing genotypes will be imputed
    by column mean. 
+ `chr`: Target chromosome. This must match the `CHROM` field in the VCF file. 
+ `start_bp`: starting basepair (position)
+ `end_bp`: ending basepair (position)

# Optional inputs
+ `min_maf`: minimum minor allele frequency for a SNP to be imported (default `0.01`)
+ `min_hwe`: minimum HWE p-value for a SNP to be imported (default `0.0`)
+ `snps_to_keep`: Vector of SNP positions to import. If specified, only SNPs
    whose position is listed in `snps_to_keep` will be kept (default `nothing`)

# output
+ `X`: Double precision matrix containing genotypes between `start_bp` and 
    `end_bp`. All entries are {0, 1, 2} except for missing entries which is 
    imputed with column mean.
+ `df`: A `DataFrame` containing meta information on the columns of `X`, 
    including `rsid` (SNP id), `AF` (alt allele frequency), `chr`, `pos`, 
    `ref`, and `alt`
"""
function get_block(
    vcffile::String, 
    chr::Int, 
    start_bp::Int, 
    end_bp::Int;
    min_maf::Float64=0.01,
    min_hwe::Float64=0.0,
    snps_to_keep::Union{AbstractVector{Int}, Nothing}=nothing
    )
    reader = VCF.Reader(openvcf(vcffile, "r"))
    T = Float64
    n = VCFTools.nsamples(vcffile)
    nsnps = 0
    X = ElasticArray{T}(undef, n, 0)
    df = DataFrame("rsid"=>String[], "AF"=>T[], "chr"=>Int[], "pos"=> Int[], 
                   "ref"=> String[], "alt"=>String[])
    model = :additive # additive genetic model
    prev_pos = 0

    record = VCF.Record()
    while !eof(reader)
        read!(reader, record)

        # current SNP info
        chr_i = VCF.chrom(record)
        pos_i = VCF.pos(record)
        alt_i = VCF.alt(record)
        gtkey = VariantCallFormat.findgenokey(record, "GT")
        _, _, _, _, _, alt_freq, _, _, _, maf, hwepval = gtstats(record, nothing)

        # quick exit
        validate(record, chr_i, pos_i, prev_pos, alt_i, gtkey)
        chr_i = parse(Int, chr_i)
        chr_i < chr && continue
        chr_i > chr && break
        pos_i < start_bp && continue
        pos_i > end_bp && break
        (isnan(maf) || maf < min_maf) && continue
        hwepval < min_hwe && continue

        # add column to X
        nsnps += 1
        resize!(X, n, nsnps)
        for i in 1:n
            geno = record.genotype[i]
            if gtkey > lastindex(geno) || geno_ismissing(record, geno[gtkey])
                # Missing: fill with MAF
                X[i, end] = 2alt_freq
            else # not missing
                # "0" (REF) => 0x30, "1" (ALT) => 0x31
                a1 = record.data[geno[gtkey][1]] == 0x31
                a2 = record.data[geno[gtkey][3]] == 0x31
                X[i, end] = convert_gt(T, (a1, a2), model)
            end
        end

        # save record info
        prev_pos = pos_i
        push!(df, [join(VCF.id(record), ','), alt_freq, chr_i, pos_i, 
                  VCF.ref(record), alt_i[1]])
    end
    close(reader)

    if !isnothing(snps_to_keep) 
        idx = filter!(!isnothing, indexin(snps_to_keep, df[!, "pos"]))
        df = df[idx, :]
        X = X[:, idx]
    end

    return Matrix(X), df
end

function validate(record, chr_i, pos_i, prev_pos, alt_i, gtkey)
    if gtkey === nothing
        error("chr $chr_i at pos $pos_i has no GT field!")
    end
    if length(alt_i) > 1
        error("Detected multiallelic marker at chr $chr_i pos $pos_i. " * 
              "Please split multiallelic markers first." )
    end
    if pos_i < prev_pos
        error("Detected VCF file is not sorted starting at pos $pos_i")
    end
    if !isdigit(chr_i[1])
        error("Detected non-integer chromosome `$chr_i` at pos $pos_i, " * 
            "please rename all CHROM field in the VCF file into an integer. " * 
            "E.g. `chr1` should be renamed into `1`."
            )
    end
end

# function to rearrange the SNP orders to that resulting S/D actually block diagonal
function rearrange_snps!(groups, group_reps, Sigma, Sigma_info)
    perm = sortperm(groups)
    iperm = invperm(perm)
    groups .= @views groups[perm]
    group_reps .= @views iperm[group_reps]
    sort!(group_reps)
    @assert issorted(groups) && issorted(group_reps)
    Sigma .= @views Sigma[perm, perm]
    Sigma_info .= @views Sigma_info[perm, :]
    return nothing
end

"""
    solve_blocks(vcffile::String, chr::Int, start_bp::Int, end_bp::Int, 
        outdir::String, hg_build::Int; [m=5], [tol=0.0001], [min_maf=0.01], 
        [min_hwe=0.0], [force_block_diag=true], 
        [method::String = "maxent"], [linkage::String="average"],
        [force_contiguous::Bool=false], [group_cor_cutoff::Float64=0.5], 
        [group_rep_cutoff::Float64=0.5], [verbose=true])

Solves the group knockoff optimization problem on provided individual-level data
and outputs the result into `outdir`. All variants that reside on chromosome 
`chr` with position between `start_bp` and `end_bp` (inclusive) will be included. 

# Note on large VCF files
Currently reading/parsing a VCF file is a single-threaded operation (even if 
it is indexed). Thus, we *strongly recommend* one to split the input VCF 
file by chromosomes, and possibly into smaller chunks, before running this
function. 

# Inputs
+ `vcffile`: A VCF file storing individual level genotypes. Must end in `.vcf` 
    or `.vcf.gz`. The ALT field for each record must be unique, i.e. 
    multiallelic records must be split first. Missing genotypes will be imputed
    by column mean. 
+ `chr`: Target chromosome. This MUST be an integer and it must match the `CHROM`
    field in your VCF file. Thus, if your VCF file has CHROM field like `chr1`, 
    `CHR1`, or `CHROM1` etc, each record must be renamed into `1`. 
+ `start_bp`: starting basepair (position)
+ `end_bp`: ending basepair (position)
+ `outdir`: Directory that the output will be stored in (must exist)
+ `hg_build`: human genome build for the VCF file, must be 19 (hg19) or 38 (hg38)

# Optional inputs (for group knockoff optimization)
+ `snps_to_keep`: Vector of SNP positions to import. If specified, only SNPs
    whose position is listed in `snps_to_keep` will be kept (default `nothing`)
+ `tol`: Convergence tolerlance for coordinate descent algorithm (default `0.0001`)
+ `min_maf`: Minimum minor allele frequency for a variable to be considered (
    default `0.01`)
+ `min_hwe`: Cutoff for hardy-weinburg equilibrium p-values. Only SNPs with 
    p-value > `min_hwe` will be included (default `0.0`)
+ `force_block_diag`: Whether to re-order the columns/rows of the correlation 
    matrix and corresponding `S` matrix so that features in the same group 
    are contiguous (default `true`). This has no impact on the final results, 
    it is simply for computational performance. 
+ `method`: group knockoff optimization algorithm, choices include "maxent" 
    (default), "mvr", "sdp", or "equi". See sec 2 of https://arxiv.org/abs/2310.15069
+ `linkage`: *cluster linkage* function to use for hierarchically clustering groups.
    It defines how the distances between features are aggregated into the distances 
    between groups. Valid choices include:
    + `:average` (default): use the mean distance between any of the cluster members
    + `:single`: use the minimum distance between any of the cluster members
    + `:complete`: use the maximum distance between any of the members
    + `:ward`: the distance is the increase of the average squared distance of a
        point to its cluster centroid after merging the two clusters
    + `:ward_presquared`: same as `:ward`, but assumes that the distances in d 
        are already squared.
+ `force_contiguous`: whether to force groups to be contiguous (default `false`).
    Note if `force_contiguous=true`, `linkage` must be `:single`)
+ `group_cor_cutoff`: correlation cutoff value for defining groups (default 
    `0.5`). Value should be between 0 and 1, where larger values correspond to 
    larger groups. 
+ `group_rep_cutoff`: cutoff value for selecting group-representatives (default
    `0.5`). Value should be between 0 and 1, where larger values correspond to 
    more representatives per group. 
+ `verbose`: whether to print informative intermediate results (default `true`)

# output
Calling `solve_blocks` will create 3 files in the directory `outdir/chr`:
+ `XXX.h5`:  This contains data (Sigma, S, groups, ..., etc) for region XXX. It 
    contains the following:
    - `D`: A `p × p` (dense) matrix corresponding to the S matrix for both the
        representative and non-representative variables. Knockoff sampling should 
        use this matrix. 
    - `S`: Matrix obtained from solving the group knockoff optimization problem 
        on the representative (group-key) variables.
    - `Sigma`: The original `p × p` correlation matrix estimated from `vcffile`
    - `SigmaInv`: Inverse of `Sigma`
    - `Sigma_reps`: The correlation matrix for the representative variables. This
        is the matrix actually used to solve the group-knockoff optimization
    - `Sigma_reps_inv`: Inverse of `Sigma_reps`
    - `group_reps`: Indices `groups` that are used as representatives (i.e. 
        group-key variables)
    - `groups`: The group membership vector
+ `Info_XXX.csv`: This includes information for each variant (chr/pos/etc) present 
    in the corresponding `.h5` file.
+ `summary_XXX.csv`: Summary file for the knockoff optimization problem
"""
function solve_blocks(
    vcffile::String,
    chr::Int,
    start_bp::Int, 
    end_bp::Int, 
    outdir::String,
    hg_build::Int; # 19 or 38
    snps_to_keep::Union{AbstractVector{Int}, Nothing}=nothing,
    # group knockoff options
    tol=0.0001, 
    min_maf=0.01,
    min_hwe=0.0,
    force_block_diag=true,
    method::String = "maxent",
    linkage::String = "average",
    force_contiguous::Bool=false,
    group_cor_cutoff::Float64=0.5,
    group_rep_cutoff::Float64=0.5,
    verbose=true
    )
    isdir(outdir) || error("output directory $outdir does not exist!")
    method = Symbol(method)
    linkage = Symbol(linkage)

    # number of simultaneous knockoffs to generate. Note this number is FIXED 
    # because in `ghostbasil_parallel.jl`, `m=5` is hard-coded into high 
    # dimensional lasso regression
    m = 5 

    # import VCF data and estimate Sigma
    import_time = @elapsed begin
        X, data_info = get_block(vcffile, chr, start_bp, end_bp, 
            min_maf=min_maf, min_hwe=min_hwe, snps_to_keep=snps_to_keep)
        size(X, 2) > 1 || 
            error("Detected 1 or fewer SNP(s) between start_bp=$start_bp and end_bp=$end_bp, exiting.")
        size(X, 1) ≥ 10 || 
            error("Detected less than 10 samples, not recommended")

        Sigma = estimate_sigma(X)
        rename!(data_info, "pos" => "pos_hg$hg_build") # associate pos with hg_build
    end

    # define groups and representatives
    def_group_time = @elapsed begin
        groups = hc_partition_groups(Symmetric(Sigma), cutoff=group_cor_cutoff, 
            linkage=linkage, force_contiguous=force_contiguous)
        group_reps = choose_group_reps(Symmetric(Sigma), groups, 
            threshold=group_rep_cutoff)
        force_block_diag && rearrange_snps!(groups, group_reps, Sigma, data_info)
    end

    # solve group knockoff optimization problem
    solve_S_time = @elapsed begin
        S, D, obj = solve_s_graphical_group(Symmetric(Sigma), groups, 
            group_reps, method, m=m, tol=tol, verbose=verbose)

        # solve S using modified Sigma (enforcing conditional independence)
        # Sigma2 = Symmetric(cond_indep_corr(Sigma, groups, group_reps))
        # S2, D2, obj2 = solve_s_graphical_group(Sigma2, groups, group_reps, method,
        #     m=m, tol=tol, verbose=verbose)
    end

    # save main result in .h5 format and summary information in .csv
    dir = joinpath(outdir, "chr$chr")
    isdir(dir) || mkpath(dir)
    h5open(joinpath(dir, "LD_start$(start_bp)_end$(end_bp).h5"), "w") do file
        write(file, "S", S)
        write(file, "D", D)
        write(file, "Sigma", Sigma)
        write(file, "SigmaInv", inv(Sigma))
        write(file, "groups", groups)
        write(file, "group_reps", group_reps)
        write(file, "Sigma_reps", Sigma[group_reps, group_reps])
        write(file, "Sigma_reps_inv", inv(Sigma[group_reps, group_reps]))
    end
    CSV.write(joinpath(dir, "Info_start$(start_bp)_end$(end_bp).csv"), data_info)
    open(joinpath(dir, "summary_start$(start_bp)_end$(end_bp).csv"), "w") do io
        println(io, "chr,$chr")
        println(io, "start_bp,$start_bp")
        println(io, "end_bp,$end_bp")
        println(io, "m,$m")
        println(io, "p,", size(Sigma, 1))
        println(io, "nreps,", length(group_reps))
        println(io, "max_group_size,", countmap(groups) |> values |> collect |> maximum)
        println(io, "max_rep_group_size,", countmap(groups[group_reps]) |> values |> collect |> maximum)
        println(io, "import_time,$import_time")
        println(io, "def_group_time,$def_group_time")
        println(io, "solve_S_time,$solve_S_time")
    end

    if verbose
        println("Time to read VCF file: $import_time")
        println("Time to define groups: $def_group_time")
        println("Time to perform group knockoff optimization: $solve_S_time")
        println("Done!")
    end

    return nothing
end
