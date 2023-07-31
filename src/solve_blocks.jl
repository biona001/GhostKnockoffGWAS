# put this block of code in a file `solve_blocks.jl`
# usage: julia solve_blocks.jl chr start_pos end_pos method outdir
# where `chr` = 1, 2, ..., 22
#       `start_pos` = integer
#       `end_pos` = inteter
#       `method` = maxent, mvr, or sdp
#       `outdir` = output directory (must exist)
# if we want to include typed SNPs only, we must provide a list of typed SNP's position
# into the argument "snps_to_keep". Not providing means we will use all SNPs.

# For Sherlock users, load needed modules
# ml julia/1.8.4 R/4.0.2 openssl/3.0.7 java/11.0.11 python/3.9.0

sleep(60rand()) # prevent too many jobs starting at the same time

using CSV
using DataFrames
using EasyLD
using Knockoffs
using LinearAlgebra
using StatsBase
using DelimitedFiles
using JLD2
using HDF5
using RCall
R"library(liftOver)"

function graphical_group_S(
    bm_file::String, # path to the `.bm` hail block matrix folder
    ht_file::String, # path to the `.ht` hail variant index folder
    chr::String, # 1, 2, ..., 22, or X/Y/M
    start_pos::Int, 
    end_pos::Int, 
    outdir::String;
    m=5,
    tol=0.0001, 
    min_maf=0.01, # only SNPs with maf > min_maf will be included in Sigma
    snps_to_keep=nothing, # list of position (if provided, only SNPs in snps_to_keep will be included in Sigma)
    force_block_diag=true, # whether to reorder row/col of S/Sigma/etc so S is returned as block diagonal 
    add_hg38_coordinates=true, # whether to augment Sigma_info with hg38 coordinates. SNPs that cannot be mapped to hg38 will be deleted
    liftOver_chain_dir::String = "",
    method::Symbol = :maxent,
    group_def::String="hc",
    rk::Number=Inf, # minimum rank of Σ before truncating the remaining eigvals to `min_maf`
    verbose=true
    )
    isdir(outdir) || error("output directory $outdir does not exist!")
    group_def ∈ ["hc", "id"] || error("group_def should be \"hc\" or \"id\"")
    add_hg38_coordinates && !isfile(liftOver_chain_dir) && 
        error("To add hg38 coordinates, a liftOver chain file is required")

    # import Sigma
    import_sigma_time = @elapsed begin
        bm = hail_block_matrix(bm_file, ht_file);
        Sigma, Sigma_info = get_block(bm, chr, start_pos, end_pos, 
            min_maf=min_maf, snps_to_keep=snps_to_keep, enforce_psd=true, rk=rk)
    end

    # append hg38 coordinates to Sigma_info and remove SNPs that can't be converted to hg38
    if add_hg38_coordinates
        Sigma, Sigma_info = augment_hg38(Sigma_info, Sigma, chr, liftOver_chain_dir)
    end

    # define groups and representatives
    def_group_time = @elapsed begin
        groups = group_def == "hc" ? 
            hc_partition_groups(Symmetric(Sigma), cutoff=0.5, linkage=:average) : 
            id_partition_groups(Symmetric(Sigma), rss_target=0.5)
        group_reps = choose_group_reps(Symmetric(Sigma), groups, threshold=0.5)
    end

    # reorder SNPs so D2/S2 is really block diagonal
    if force_block_diag
        rearrange_snps!(groups, group_reps, Sigma, Sigma_info)
    end

    # solve for S
    solve_S_time = @elapsed begin
        # solve S using original Sigma
        S, D, obj = solve_s_graphical_group(Symmetric(Sigma), groups, group_reps, method,
            m=m, tol=tol, verbose=verbose)

        # solve S using modified Sigma (enforcing conditional independence)
        # Sigma2 = cond_indep_corr(Sigma, groups, group_reps)
        # S2, D2, obj2 = solve_s_graphical_group(Symmetric(Sigma2), groups, group_reps, method,
        #     m=m, tol=tol, verbose=verbose)
    end

    # save result in .h5 format
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
    CSV.write(joinpath(dir, "Info_start$(start_pos)_end$(end_pos).csv"), Sigma_info)
    open(joinpath(dir, "summary_start$(start_pos)_end$(end_pos).csv"), "w") do io
        println(io, "chr,$chr")
        println(io, "start_pos,$start_pos")
        println(io, "end_pos,$end_pos")
        println(io, "m,$m")
        println(io, "p,", size(Sigma, 1))
        println(io, "nreps,", length(group_reps))
        println(io, "max_group_size,", countmap(groups) |> values |> collect |> maximum)
        println(io, "max_rep_group_size,", countmap(groups[group_reps]) |> values |> collect |> maximum)
        println(io, "import_sigma_time,$import_sigma_time")
        println(io, "def_group_time,$def_group_time")
        println(io, "solve_S_time,$solve_S_time")
    end

    return nothing
end

# use liftOver R package to add hg38 coordinates to Sigma_info
# SNPs that can't be converted, or if chr does not match after conversion, are deleted
function augment_hg38(Sigma_info, Sigma, chr, liftOver_chain_dir)
    pos = Sigma_info[!, "pos"]
    @rput chr pos liftOver_chain_dir
    R"""
    df<-cbind(data.frame(paste0('chr',chr)),pos,pos)
    colnames(df)<-c('chr','start','end')
    temp.Granges<-makeGRangesFromDataFrame(df)
    chain <- import.chain(liftOver_chain_dir)
    converted<-data.frame(liftOver(temp.Granges,chain))
    """
    @rget converted
    success_idx, pos_hg38 = falses(size(Sigma, 1)), Int[]
    for row in eachrow(converted)
        row["seqnames"] == "chr$chr" || continue # check chr match
        success_idx[row["group"]] = true
        push!(pos_hg38, row["start"])
    end
    Sigma_info_new = Sigma_info[success_idx, :]
    Sigma_new = Sigma[success_idx, success_idx]
    Sigma_info_new[!, :pos_hg19] = pos[success_idx]
    Sigma_info_new[!, :pos_hg38] = pos_hg38
    select!(Sigma_info_new, Not(:pos))
    return Sigma_new, Sigma_info_new
end

# function to rearrange the SNP orders to that resulting S2 and D2 actually block diagonal
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

# some directories
bm_file = "/scratch/users/bbchu/LD_matrices/UKBB.EUR.ldadj.bm"
ht_file = "/scratch/users/bbchu/LD_matrices/UKBB.EUR.ldadj.variant.ht"
liftOver_chain_dir = "/oak/stanford/groups/zihuai/GeneticsResources/LiftOver/hg19ToHg38.over.chain"

# inputs
chr = ARGS[1]
start_pos = parse(Int, ARGS[2])
end_pos = parse(Int, ARGS[3])
method = Symbol(ARGS[4]) # maxent, mvr, or sdp
group_def = ARGS[5] # "hc" or "id"
rk = parse(Int, ARGS[6]) #50, or 300
outdir = ARGS[7]

# testing one region
# chr = "1"
# start_pos = 100826405
# end_pos = 102041015
# m=5
# tol=0.0001 
# min_maf=0.01
# force_block_diag=true
# add_hg38_coordinates=true
# method = :maxent
# group_def = "hc"
# verbose=true
# outdir = "/oak/stanford/groups/zihuai/pan_ukb_group_knockoffs/maxent_hc"

# import typed SNP positions and ref/alt 
# snps_to_keep = nothing
bimfile = "/oak/stanford/groups/zihuai/UKB_data/genotyped_call/ukb_cal_chr$(chr)_v2.bim"
typed_df = CSV.read(bimfile, DataFrame, header=false)
snps_to_keep = typed_df[!, 4]

graphical_group_S(bm_file, ht_file, chr, start_pos, 
    end_pos, outdir, snps_to_keep=snps_to_keep, 
    method=method, group_def=group_def, rk=rk,
    liftOver_chain_dir=liftOver_chain_dir)