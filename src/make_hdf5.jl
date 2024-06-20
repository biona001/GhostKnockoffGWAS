function estimate_sigma(X::AbstractMatrix)
    return cor(X)
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
