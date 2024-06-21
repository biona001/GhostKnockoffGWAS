function estimate_sigma(X::AbstractMatrix)
    n, p = size(X)
    if n > p
        return cor(X)
    else
        error("p > n case not handled yet!")
    end
end

"""
    get_block(vcffile, )

Imports genotype data of a VCF file from pos `start_bp` to `end_bp` as 
double precision numeric matrix.

# Inputs
+ `vcffile`: A VCF file storing individual level genotypes. Must end in `.vcf` 
    or `.vcf.gz`. The ALT field for each record must be unique, i.e. 
    multiallelic records must be split first. Missing genotypes will be imputed.
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
+ `chrs`: Chromosome for each column of `X`
+ `poss`: Basepair position for each column of `X`
+ `ref`: Reference allele for each column of `X`
+ `alt`: Alt allele for each column of `X`
"""
function get_block(
    vcffile::String, 
    chr::String, 
    start_bp::Int, 
    end_bp::Int;
    min_maf::Float64=0.01,
    min_hwe::Float64=0.0,
    snps_to_keep::Union{AbstractVector{Int}, Nothing}=nothing
    )
    reader = VCF.Reader(openvcf(vcffile, "r"))
    T = Float64
    n = nsamples(vcffile)
    nsnps = 0
    X = ElasticArray{T}(undef, n, 0)
    df = DataFrame("rsid"=>String[], "AF"=>T[], "chr"=>String[], "pos"=> Int[], 
                   "ref"=> String[], "alt"=>String[])
    model = :additive # additive genetic model

    record = VCF.Record()
    while !eof(reader)
        read!(reader, record)

        # check chr
        chr_i = VCF.chrom(record)
        chr_i > chr && break

        # check position
        pos_i = VCF.pos(record)
        pos_i < start_bp && continue
        pos_i > end_bp && break

        # check SNP is usable
        alt_i = VCF.alt(record)
        validate(record, alt_i)
        _, _, _, _, _, alt_freq, _, _, _, maf, hwepval = gtstats(record, nothing)
        maf < min_maf && continue
        hwepval < min_hwe && continue
        gtkey = VariantCallFormat.findgenokey(record, "GT")

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
        push!(df, [join(VCF.id(record), ','), alt_freq, chr_i, pos_i, 
                  VCF.ref(record), alt_i[1]])
    end

    if !isnothing(snps_to_keep) 
        idx = filter!(!isnothing, indexin(snps_to_keep, df[!, "pos"]))
        df = df[idx, :]
        X = X[:, idx]
    end

    return Matrix(X), df
end
# using VCFTools, VariantCallFormat, ElasticArrays, DataFrames
# vcffile = "/oak/stanford/groups/zihuai/paisa/chr1.vcf.gz"
# chr = "1"
# min_maf = 0.0
# min_hwe = 0.0
# start_bp = 58814 # first 10 records
# end_bp = 801536
# @time X, df = get_block(vcffile, chr, start_bp, end_bp, min_maf=min_maf, min_hwe=min_hwe)
# 0.012035 seconds (46.78 k allocations: 13.440 MiB)

# chr = "1"
# start_bp = 249208612 # last 10 records
# end_bp = 249222527
# @time X, df = get_block(vcffile, chr, start_bp, end_bp, min_maf=min_maf, min_hwe=min_hwe)
# 44.146655 seconds (196.81 M allocations: 53.676 GiB, 2.20% gc time)

function validate(record, alt_i)
    if VariantCallFormat.findgenokey(record, "GT") === nothing
        error("chr $chr_i at pos $pos_i has no GT field!")
    end
    if length(alt_i) > 1
        error("Detected multiallelic marker at chr $chr_i pos $pos_i. " * 
              "Please split multiallelic markers first." )
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
    Sigma.data .= @views Sigma[perm, perm]
    Sigma_info .= @views Sigma_info[perm, :]
    return nothing
end

function solve_blocks(
    vcffile::String,
    chr::String,
    start_bp::Int, 
    end_bp::Int, 
    outdir::String,
    hg_build::Int; # 19 or 38
    snps_to_keep::Union{AbstractVector{Int}, Nothing}=nothing, # list of position (if provided, only SNPs in snps_to_keep will be included in Sigma)
    # group knockoff options
    m=5,
    tol=0.0001, 
    min_maf=0.01, # only SNPs with maf > min_maf will be included in Sigma
    min_hwe=0.0, # only SNPs with hwe > min_hwe will be included in Sigma
    force_block_diag=true, # whether to reorder row/col of S/Sigma/etc so S is returned as block diagonal 
    method::String = "maxent",
    group_def::String="hc",
    group_cor_cutoff::Float64=0.5,
    group_rep_cutoff::Float64=0.5,
    verbose=true
    )
    # using VCFTools, VariantCallFormat, ElasticArrays, DataFrames, Statistics, LinearAlgebra, Knockoffs
    # vcffile = "/oak/stanford/groups/zihuai/paisa/chr1.vcf.gz"
    # chr = "1"
    # min_maf = 0.01
    # min_hwe = 1e-8
    # start_bp = 58814 # first 500 records
    # end_bp = 2338569
    # snps_to_keep = nothing
    # hg_build = 19
    # group_def = "hc"
    # group_cor_cutoff = 0.5
    # group_rep_cutoff = 0.5
    # method = "maxent"
    # m = 5
    # tol = 0.0001
    # verbose = true

    isdir(outdir) || error("output directory $outdir does not exist!")
    group_def ∈ ["hc", "id"] || error("group_def should be \"hc\" or \"id\"")
    method = Symbol(method)

    # import VCF data and estimate Sigma
    import_time = @elapsed begin
        X, data_info = get_block(vcffile, chr, start_bp, end_bp, 
            min_maf=min_maf, min_hwe=min_hwe, snps_to_keep=snps_to_keep)
        Sigma = Symmetric(estimate_sigma(X))
        rename!(data_info, "pos" => "pos_hg$hg_build") # associate pos with hg_build
    end

    # define groups and representatives
    def_group_time = @elapsed begin
        groups = group_def == "hc" ? 
            hc_partition_groups(Sigma, cutoff=group_cor_cutoff, linkage=:average) : 
            id_partition_groups(Sigma, rss_target=group_cor_cutoff)
        group_reps = choose_group_reps(Sigma, groups, threshold=group_rep_cutoff)
    end

    # reorder SNPs so D2/S2 is really block diagonal
    if force_block_diag
        rearrange_snps!(groups, group_reps, Sigma, data_info)
    end

    # solve for S
    solve_S_time = @elapsed begin
        S, D, obj = solve_s_graphical_group(Sigma, groups, group_reps, method,
            m=m, tol=tol, verbose=verbose)

        # solve S using modified Sigma (enforcing conditional independence)
        # Sigma2 = Symmetric(cond_indep_corr(Sigma, groups, group_reps))
        # S2, D2, obj2 = solve_s_graphical_group(Sigma2, groups, group_reps, method,
        #     m=m, tol=tol, verbose=verbose)
    end

    # save main result in .h5 format and summary information in .csv
    dir = joinpath(outdir, "chr$chr")
    isdir(dir) || mkpath(dir)
    JLD2.save(
        joinpath(dir, "LD_start$(start_pos)_end$(end_pos).h5"), 
        Dict(
            "S" => S,
            "D" => D,
            "Sigma" => Sigma,
            "SigmaInv" => inv(Sigma),
            "groups" => groups,
            "group_reps" => group_reps,
            "Sigma_reps" => Sigma[group_reps, group_reps],
            "Sigma_reps_inv" => inv(Sigma[group_reps, group_reps])
        )
    )
    CSV.write(joinpath(dir, "Info_start$(start_pos)_end$(end_pos).csv"), data_info)
    open(joinpath(dir, "summary_start$(start_pos)_end$(end_pos)"), "w") do io
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

    return nothing
end
