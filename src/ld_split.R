# Please checkout the main documentation for use of this script:
# https://biona001.github.io/GhostKnockoffGWAS/dev/man/solveblocks/#Determining-start_bp-and-end_bp

# usage: Rscript --vanilla ld_split.R arg1 arg2 arg3 arg4 arg5 arg6
# arg1 = chromosome number (must be an integer)
# arg2 = path to PLINK binary file (must end in `.bed` extension)
# arg3 = path to FBM file (without extensions. If this file doesn't exist, it
#        will be generated)
# arg4 = path to output file 
# arg5 = thr_r2, this is the thr_r2 used by snp_ldsplit. All correlation
#        smaller than thr_r2 are set to 0
# arg6 = max_r2, this is the max_r2 used by snp_ldsplit. This is the maximum 
#        acceptable correlation for SNPs in different blocks. 

# Description:
# This R code runs the snp_ldsplit function: https://privefl.github.io/bigsnpr/reference/snp_ldsplit.html
# It operates directly on PLINK binary files and produces roughly independent
# blocks stored as a plain text file. Note that it is generally recommended to 
# run this code with at least 16GB of RAM, possibly more are needed for dense
# array genotypes. 

# required libraries
library("bigsnpr")
library("dplyr")

# enable command line arguments
args = commandArgs(TRUE)
chr = as.numeric(args[1])
plinkfile = args[2]
fbmfile = args[3]
outfile = args[4]
thr_r2 = as.numeric(args[5])
max_r2 = as.numeric(args[6]) 

# import PLINK data as FBM (file backed matrix) format
rdsfile <- paste0(fbmfile, ".rds")
if (!file.exists(rdsfile)){snp_readBed2(plinkfile, backingfile = fbmfile)} 
x <- snp_attach(rdsfile)

# estimate correlation matrix
corr <- snp_cor(x$genotypes, infos.pos=x$map$physical.pos, ncores=1)

# compute LD regions
m <- ncol(corr)
max_sizes <- c(1000, 1500, 3000, 6000, 10000)
max_sizes <- max_sizes[max_sizes <= dim(corr)[1]]
splits <- snp_ldsplit(corr, thr_r2 = thr_r2, min_size = 500, max_size = max_sizes, max_r2 = max_r2)

# run snp_ldsplit with default objective
splits$cost2 <- sapply(splits$all_size, function(sizes) sum(sizes^2))
best_split <- splits %>%
    arrange(cost2 * sqrt(5 + cost)) %>%
    print() %>%
    slice(1) %>%
    print()
all_size <- best_split$all_size[[1]]
best_grp <- rep(seq_along(all_size), all_size)

# get position of LD split
unique_grp <- unique(best_grp)
start_pos <- integer(length(unique_grp))
end_pos <- integer(length(unique_grp))
for (i in seq_along(unique_grp)) {
  start_pos[i] <- min(which(best_grp == unique_grp[i]))
  end_pos[i] <- max(which(best_grp == unique_grp[i]))
}

# save result
pos <- x$map$physical.pos
result <- data.frame(
    chr = rep(chr, length(start_pos)),
    start = pos[start_pos], 
    stop = pos[end_pos]
)
write.table(result, outfile, row.names = FALSE, quote=FALSE, sep="\t")
