"""
    estimate_sigma(X::AbstractMatrix; [enforce_psd=true], [min_eigval=1e-5])
    estimate_sigma(X::AbstractMatrix, C::AbstractMatrix; [enforce_psd=true],
        [min_eigval=1e-5])

Estimate LD matrices from data `X`, accounting for covariates `C` if there are any. 
We adopt the method for Pan-UKB described in 
`https://pan-dev.ukbb.broadinstitute.org/docs/ld/index.html#ld-matrices`.
If `enforce_psd=true`, then the correlation matrix will be scaled so that the 
minimum eigenvalue is `min_eigval`.
"""
function estimate_sigma(X::AbstractMatrix, C::AbstractMatrix;
    enforce_psd::Bool=true, min_eigval::Float64 = 1e-5)
    # check for errors
    n = size(X, 1)
    n == size(C, 1) || error("Sample size in X and C should be the same")

    # pan-ukb routine
    Xc = StatsBase.zscore(X, mean(X, dims=1), std(X, dims=1))
    Mc = size(C, 2) > 1 ? I - C * inv(Symmetric(C' * C)) * C' : Diagonal(ones(n))
    Xadj = Mc * Xc
    Sigma = Xadj' * Xadj / n

    # numerical stability
    if enforce_psd
        evals, evecs = eigen(Sigma)
        evals[findall(x -> x < min_eigval, evals)] .= min_eigval
        Sigma = evecs * Diagonal(evals) * evecs'
    end

    # scale to correlation matrix
    StatsBase.cov2cor!(Sigma, sqrt.(diag(Sigma)))

    return Sigma
end
estimate_sigma(X; enforce_psd::Bool=true, min_eigval::Float64 = 1e-5) = 
    estimate_sigma(X, zeros(size(X, 1), 0); enforce_psd=enforce_psd, min_eigval=min_eigval)

"""
    get_block(file::String, chr::Int, start_bp::Int, end_bp::Int;
        [min_maf::Float64=0.01], [min_hwe::Float64=0.0], 
        [snps_to_keep::Union{AbstractVector{Int}, Nothing}=nothing])

Imports genotype data of a VCF or binary PLINK file from pos `start_bp` to `end_bp` 
as double precision numeric matrix. Each row is a sample and each column is a SNP.
All entries are {0, 1, 2} except for missing entries which is imputed with 
column mean.

# Inputs
+ `file`: A VCF file (ending in `.vcf` or `.vcf.gz`) or binary PLINK (ending in 
    `.bed`) file storing individual level genotypes. If VCF file is used, the 
    ALT field for each record must be unique, i.e. multiallelic records must be
    split first. Missing genotypes will be imputed by column mean. 
+ `chr`: Target chromosome. This must be an integer and, for VCF files, it must 
    match the `CHROM` field in the VCF file (i.e. if the VCF CHROM field is e.g.
    `chr1`, it must be renamed into `1`). 
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
    `ref`, and `alt`. For PLINK inputs, we count the number of A2 alleles, and 
    we assign `A1 = ref` and `A2 = alt`
"""
function get_block(
    file::String, 
    chr::Int, 
    start_bp::Int, 
    end_bp::Int;
    min_maf::Float64=0.01,
    min_hwe::Float64=0.0,
    snps_to_keep::Union{AbstractVector{Int}, Nothing}=nothing
    )
    is_vcf = endswith(file, ".vcf") || endswith(file, ".vcf.gz")
    is_plink = endswith(file, ".bed")
    is_vcf || is_plink || error(
        "Invalid input file. Genotype file should be in VCF (ends in .vcf " *
        "or .vcf.gz) or binary PLINK (ends in .bed) format."
        )
    if is_plink
        plinkprefix = file[1:end-4]
        isfile(plinkprefix * ".bim") && isfile(plinkprefix * ".fam") || 
            error("Detected PLINK input but .bim or .fam file not found")
    end

    if is_vcf
        return get_VCF_block(file, chr, start_bp, end_bp, min_maf=min_maf,
                min_hwe = min_hwe, snps_to_keep=snps_to_keep)
    else 
        return get_PLINK_block(file, chr, start_bp, end_bp, min_maf=min_maf,
                min_hwe = min_hwe, snps_to_keep=snps_to_keep)
    end
end

function get_VCF_block(
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
        push!(df, [get_record_id(record), alt_freq, chr_i, pos_i, 
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

"""
    get_AF(x::AbstractMatrix)

Compute alternate allele frequency on dosage matrix `x`, 
skipping missing values if there are any
"""
function get_AF(x::AbstractMatrix)
    return mean.(skipmissing.(eachcol(x))) ./ 2
end

function get_PLINK_block(
    bedfile::String, 
    chr::Int, 
    start_bp::Int, 
    end_bp::Int;
    min_maf::Float64=0.01,
    min_hwe::Float64=0.0,
    snps_to_keep::Union{AbstractVector{Int}, Nothing}=nothing
    )
    xdata = SnpData(bedfile[1:end-4])

    # import target region of X as numeric matrix, after MAF filtering
    idx = findall(x -> 
        x.chromosome == string(chr) && start_bp ≤ x.position ≤ end_bp, 
        eachrow(xdata.snp_info))
    mafs = maf(@view(xdata.snparray[:, idx]))
    idx2 = idx[findall(x -> x ≥ min_maf, mafs)]
    X = convert(Matrix{Float64}, @view(xdata.snparray[:, idx2]), impute=true)

    # SNP info
    df = xdata.snp_info[idx2, [1, 2, 4, 5, 6]]
    rename!(df, "chromosome"=>"chr", "snpid"=>"rsid", "position"=>"pos", 
        "allele1"=>"ref", "allele2"=>"alt")
    df[!, "AF"] = get_AF(X)
    df = df[:, [2, 6, 1, 3, 4, 5]] # reorder columns into rsid, AF, chr, pos, ref, alt

    # hwe filtering
    if min_hwe > 0
        n00s = [count(iszero, x) for x in eachcol(X)]
        n01s = [count(isone, x) for x in eachcol(X)]
        n11s = [count(xi -> xi == 2, x) for x in eachcol(X)]
        hwes = SnpArrays.hwe.(n00s, n01s, n11s)
        idx3  = findall(x -> x > min_hwe, hwes)
        X = X[:, idx3]
        df = df[idx3, :]
    end

    if !isnothing(snps_to_keep) 
        idx4 = filter!(!isnothing, indexin(snps_to_keep, df[!, "pos"]))
        df = df[idx4, :]
        X = X[:, idx4]
    end

    return Matrix(X), df
end

# return "." if id is missing
function get_record_id(record)
    return isempty(record.id) ? "." : join(VCF.id(record), ',')
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

# function to rearrange the SNP orders so that groups are contiguous
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
    get_covariates(covfile::String, genotype_file::String)

Read covariates from `covfile` and reorders the rows so that the rows come in 
the same order as sample IDs specified in (VCF or binary PLINK) `genotype_file`.

# Inputs
+ `covfile`: A comma- or tab-separated file containing sample covariates. The
    first row should be a header row. The first column should be sample IDs (not
    necessary to be in the sample order as genotype files and all other columns
    will be used as additional covariates. Note if genotypes are stored in binary
    PLINK format, then the sample ID column in the covariate file should be 
    FID_IID (that is, the first 2 columns of the .fam file merged by an 
    underscore).
+ `genotype_file`: A VCF or binary PLINK file storing individual level genotypes.
    Must end in `.vcf`, `.vcf.gz`, or `.bed`. 
"""
function get_covariates(covfile::String, genotype_file::String)
    if endswith(genotype_file, ".vcf") || endswith(genotype_file, ".vcf.gz")
        sampleIDs = sampleID(genotype_file)
    elseif endswith(genotype_file, ".bed")
        famfile = genotype_file[1:end-4] * ".fam"
        fam_df = CSV.read(famfile, DataFrame, header=false)
        sampleIDs = string.(fam_df[!, 1], "_", fam_df[!, 2])
    else
        error("Genotype file should be in VCF (ends in .vcf " *
            "or .vcf.gz) or binary PLINK (ends in .bed) format.")
    end

    # read covariate data and match sample IDs
    covdata = CSV.read(covfile, DataFrame)
    cov_sampleIDs = string.(covdata[!, 1])
    idx = indexin(sampleIDs, cov_sampleIDs)
    if length(idx) != length(sampleIDs)
        error("A covariate file was supplied but >=1 genotyped sample(s)" * 
            " does not have covariate data. Please check if the covariate" * 
            " file has the correct sample IDs.")
    end
    C = Matrix(covdata[idx, 2:end])
    return C::Matrix{Float64}
end

"""
    solve_blocks(file::String, chr::Int, start_bp::Int, end_bp::Int, 
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
it is indexed). Thus, we *strongly recommend* one convert to binary PLINK format.

# Inputs
+ `file`: A VCF or binary PLINK file storing individual level genotypes. Must
    end in `.vcf`, `.vcf.gz`, or `.bed`. If a VCF file is used, the ALT field for
    each record must be unique, i.e. multiallelic records must be split first. 
    Missing genotypes will be imputed by column mean. 
+ `covfile`: An optional comma- or tab-separated file containing sample covariates 
    (e.g. sex, age, PCs). This argument can be an empty string. The supplied 
    covariates will be used to improve LD estimation. The first column should be
    sample IDs (not necessary to be in the sample order as VCF or PLINK files)
    and all other columns will be used as additional covariates.
+ `chr`: Target chromosome. This MUST be an integer and it must match the `CHROM`
    field in your VCF/PLINK file. For example, if your VCF file has CHROM field
    like `chr1`, `CHR1`, or `CHROM1` etc, they must be renamed into `1`. 
+ `start_bp`: starting basepair (position)
+ `end_bp`: ending basepair (position)
+ `outdir`: Directory that the output will be stored in (must exist)
+ `hg_build`: human genome build for position of each SNP, must be 19 (hg19) or
    38 (hg38)

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
    file::String,
    covfile::String,
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

    # import data and estimate Sigma
    import_time = @elapsed begin
        X, data_info = get_block(file, chr, start_bp, end_bp, 
            min_maf=min_maf, min_hwe=min_hwe, snps_to_keep=snps_to_keep)
        size(X, 2) > 1 || 
            error("Detected 1 or fewer SNP(s) between start_bp=$start_bp and " * 
                  "end_bp=$end_bp in chr $chr, exiting.")
        size(X, 1) ≥ 10 || 
            error("Detected less than 10 samples, not recommended")

        # read covariates, if any
        C = covfile == "" ? zeros(size(X, 1), 0) : get_covariates(covfile, file)

        # estimate correlation matrix
        Sigma = estimate_sigma(X, C)
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
        println("Time to read file: $import_time")
        println("Time to define groups: $def_group_time")
        println("Time to perform group knockoff optimization: $solve_S_time")
        println("Done!")
    end

    return nothing
end
