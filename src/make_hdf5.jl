function estimate_sigma(X::AbstractMatrix)
    return cor(X)
end

"""
    get_block(vcffile, )

Imports genotype data of a VCF file from pos `start_bp` to `end_bp` as 
double precision numeric matrix.

# Inputs
+ `vcffile`: A VCF file storing individual level genotypes. Must end in `.vcf` 
    or `.vcf.gz`. The CHROM field for all SNPs must equal to the input `chr` 
    for computational reasons. Finally, the ALT field for each record must be 
    unique, i.e. multiallelic records must be split first. 
+ `chr`: 

# Optional inputs
+ 

# output
+ ``
"""
function get_block(
    vcffile::String, 
    chr::Union{String, Int}, 
    start_bp::Union{String, Int}, 
    end_bp::Union{String, Int};
    min_maf::Float64=0.01,
    min_hwe::Float64=0.0
    )
    chr = string(chr)
    start_bp = Int(start_bp)
    end_bp = Int(end_bp)

    reader = VCF.Reader(openvcf(vcffile, "r"))
    T = Float64
    n = nsamples(vcffile)
    nsnps = 0
    X = ElasticArray{T}(undef, n, 0)
    chrs, poss, refs, alts = String[], Int[], String[], String[]
    model = :additive # additive genetic model

    for record in reader
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
        push!(chrs, chr_i)
        push!(poss, pos_i)
        push!(refs, VCF.ref(record))
        push!(alts, alt_i[1])
    end

    return Matrix(X), chrs, poss, refs, alts
end
# using VCFTools, VariantCallFormat, ElasticArrays
# vcffile = "/oak/stanford/groups/zihuai/paisa/chr1.vcf.gz"
# chr = 1 
# min_maf = 0.0
# min_hwe = 0.0
# start_bp = 58814 # first 10 records
# end_bp = 801536
# @time X, chr, pos, ref, alt = get_block(vcffile, chr, start_bp, end_bp, min_maf=min_maf, min_hwe=min_hwe)
# 0.014761 seconds (60.13 k allocations: 18.505 MiB, 9.10% gc time)

# chr = 1 
# start_bp = 249208612 # last 10 records
# end_bp = 249222527
# @time X, chr, pos, ref, alt = get_block(vcffile, chr, start_bp, end_bp, min_maf=min_maf, min_hwe=min_hwe)
# 54.847166 seconds (263.11 M allocations: 78.248 GiB, 2.47% gc time)


function validate(record, alt_i)
    if VariantCallFormat.findgenokey(record, "GT") === nothing
        error("chr $chr_i at pos $pos_i has no GT field!")
    end
    if length(alt_i) > 1
        error("Detected multiallelic marker at chr $chr_i pos $pos_i. " * 
              "Please split multiallelic markers first." )
    end
end
