# usage: Rscript --vanilla manhattan.R arg1 arg2 arg3 arg4
# arg1 = main output file from GhostKnockoffGWAS
# arg2 = output directory
# arg3 = output filename (without extensions) to be used for both plots, 
#       e.g. phenotype name
# arg4 = target FDR in percentage

# This is R code modified from https://github.com/biona001/ghostknockoff-gwas-reproducibility/blob/main/he_et_al/GKL_Manhattan.R
# It operates directly on the output of GhostKnockoffGWAS, and
# produces 2 Manhattan plots using the R package CMplot.
# The first is a regular manhattan plot based on marginal p-values. 
# The second is a knockoff-equivalent manhattan plot based on W statistics.
# Independently significant SNPs are highlighted, where 2 discovered SNP is
# considered independent if they are at least 1Mb apart. 

# required libraries
library(data.table)
library(plyr)
library(dplyr)
library(CMplot)

# Enable command line arguments
args = commandArgs(TRUE)
input_file = args[1]
out_dir = args[2]
out_filename = args[3]
target_fdr = as.numeric(args[4])
if (!file.exists(input_file)){stop(paste0("input file ", input_file, " does not exist"))}
if (!dir.exists(out_dir)){stop(paste0("output directory ", out_dir, " does not exist"))}

# CMPlot output plots in current directory, so we cd to the output dir, make plots, then cd back
original_dir = getwd()
setwd(out_dir)

# for testing
# input_file = "/scratch/users/bbchu/GhostKnockoffGWAS/data/example_output.txt"
# out_dir = "/scratch/users/bbchu/GhostKnockoffGWAS/data"
# out_filename = "example_plot"
# target_fdr = 0.1

# input_file = "/scratch/users/bbchu/GhostKnockoffGWAS/67studies/results/BMI_Pulit_2018:WHR-Adjusted-For-BMI-Males/seed4001.permuteZ.txt"
# out_dir = "/scratch/users/bbchu/GhostKnockoffGWAS/67studies/results/BMI_Pulit_2018:WHR-Adjusted-For-BMI-Males/"
# out_filename = "manhattan"
# target_fdr = 0.1

# input_file = "/scratch/users/bbchu/GhostKnockoffGWAS/data/test_alzheimers_meta.txt"
# out_dir = "/scratch/users/bbchu/GhostKnockoffGWAS/data"
# out_filename = "AD_meta"
# target_fdr = 0.1

# input_file = "GK_out.txt"
# out_dir = "."
# out_filename = "chr1_plt"
# target_fdr = 0.1

#############################################################################
############ FIRST, CREATE MANHATTAN PLOT FOR MARGINAL ANALYSIS #############
############ In particular, this will generate an intermediate .csv file ####
############ needed for creating knockoff manhattan plots ###################
#############################################################################

main.text = 'Marginal association test' # plot title
memo.text = paste('MarginalAssociationTest_',out_filename, sep='') # full file name

## read result file
x1 <- fread(input_file,header=T)
x1_cols = colnames(x1)
setnames(x1, grep("^chr", x1_cols)[1], "CHR")
setnames(x1, grep("^pos", x1_cols)[1], "BP")
setnames(x1, grep("^ref", x1_cols)[1], "REF")
setnames(x1, grep("^alt", x1_cols)[1], "ALT")
x1 <- x1[, SNP:=paste(CHR,BP,REF,ALT,sep=":")]
x1 <- x1[which(!is.na(x1[,'CHR'])),]
setorderv(x1, c("CHR","BP"))
x1 <- as.data.frame(x1)
x1 <- x1[match(unique(x1$SNP),x1$SNP),]

## selected variants
x1_sug <- x1[x1$pvals<=5e-8,]
x1_sug <- x1_sug[!is.na(x1_sug$pvals),]

if(nrow(x1_sug)>0){
    x1_sug$IND <- NA
    x1_sug$TOP <- NA
    x1_sug$RAst <- NA
    x1_sug$RAen <- NA

    loc_w <- 1000000
    ra_w <- 500000

    # find independent loci
    CHs <- unique(x1_sug$CHR)
    inds <- 0
    for (ch in CHs) {
        BPs <- x1_sug[which(x1_sug$CHR==ch),"BP"]
        BP_first <- BPs[1]
        for (bp in BPs) {
            if (bp==BP_first) {
            inds <- inds + 1
            } else {
            bp_sw <- bp - BP_prv
            if (bp_sw > loc_w) {
                inds <- inds + 1
            }
            }
            x1_sug[which(x1_sug$CHR==ch & x1_sug$BP==bp),"IND"] <- inds
            BP_prv <- bp
        }
    }
    num_discoveries = tail(x1_sug$IND, n=1)

    # find top SNPs in independent loci
    for (i in x1_sug$IND) {
        x1_sug_t <- x1_sug[which(x1_sug$IND==i),c("TOP","SNP","pvals")]
        x1_sug_t[which.max(-log10(x1_sug_t$pvals)),"TOP"] <- 
            x1_sug_t[which.max(-log10(x1_sug_t$pvals)),"SNP"]
        x1_sug[which(x1_sug$IND==i),c("TOP","SNP","pvals")] <- x1_sug_t
    }

    # set range +/- extension for independent loci
    for (i in x1_sug$IND) {
        x1_sug_t <- x1_sug[which(x1_sug$IND==i),c("TOP","BP","RAst","RAen")]
        x1_sug_t[,"RAst"] <- x1_sug_t[1,"BP"] - ra_w
        x1_sug_t[,"RAen"] <- x1_sug_t[dim(x1_sug_t)[1],"BP"] + ra_w
        x1_sug[which(x1_sug$IND==i),c("TOP","BP","RAst","RAen")] <- x1_sug_t
    }

    # flag new loci
    new_loc_w <- 1000000
    x1_sug$NEW_hit <- "Y"
    AD_loci<-x1_sug
    AD_loci<-AD_loci[AD_loci[,'pvals']<=5e-8,]
    temp.new<-c()
    for(i in 1:nrow(x1_sug)){
        temp.chr<-x1_sug[i,'CHR']
        temp.pos<-x1_sug[i,'BP']
        temp.new<-c(temp.new,sum((temp.chr==AD_loci[,'CHR']) & (temp.pos>=AD_loci[,'BP']-new_loc_w) & (temp.pos<=AD_loci[,'BP']+new_loc_w), na.rm=T))
    }
    for(i in 1:nrow(x1_sug)){
        if(sum(temp.new[which(x1_sug[,'IND']==x1_sug[i,'IND'])])>0){x1_sug[i,"NEW_hit"] <- "N"}
    }

    # extract signal from all variants for respective independent loci (for plotting)
    x1_sig<-x1_sug[x1_sug[,'pvals']<=5e-8,]
    signal <- c()
    for (t in x1_sig[which(!is.na(x1_sig$TOP)),"TOP"]) {
        ra <- x1_sig[which(x1_sig$TOP==t),c("CHR","RAst","RAen","rsid","pvals","IND")]
        indt <- ra$IND
        xs <- x1[which(x1$CHR==ra$CHR & (x1$BP >= ra$RAst) & (x1$BP <= ra$RAen) ),
                    c("SNP","pvals","BP","CHR","rsid")]
        if (length(xs$SNP)>0) {
            xs <- xs[,c("SNP",'pvals')]
            if (ra$pvals<=5e-8) {xs$col <- "purple"}
            xs$text_col <- "red"; 
            if (any(x1_sig[which(x1_sig$IND==indt),"NEW_hit"] %in% "N")) {xs$text_col <- "red"}
            xs$pch <- 19
            xs$text <- NA; xs[,"text"] <- ra$rsid
            signal <- rbind(signal,xs)
        }
    }

    # keep signal with highlights
    signal_high <- signal[which(signal$pvals<=5e-8),]
    signal_top <- signal[which(signal$SNP %in% x1_sug[which(!is.na(x1_sug$TOP)),"TOP"]),]
    signal_topp <- signal_top
} else { # no discoveries
    signal_topp <- x1_sug
}

# thresholds
ths <- c(-log10(5e-8))

## CMplot
x1t <- x1[,c("SNP","CHR","BP",'pvals')]
x1t[which(x1t$pvals<=1e-50),"pvals"] <- 1e-50
x1t <- x1t[,c("SNP","CHR","BP","pvals")]
x1t[,4]<--log10(x1t[,4])
if (num_discoveries > 100) {
    CMplot(x1t, plot.type="m", LOG10=FALSE, col=c("grey30","grey60"), ylab="-log10(p)", # bin.range=c(0,500),
            chr.den.col=c("darkgreen", "yellow", "red"), main=main.text,
            threshold=ths, threshold.lty=c(2), threshold.lwd=c(1), threshold.col=c("black"),
            highlight=signal_topp$SNP, highlight.cex=1, 
            highlight.col=signal_topp$col, highlight.text.col=signal_topp$text_col,
            signal.col=c("cornflowerblue"),signal.cex=c(1),
            file="jpg",file.name=memo.text,dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)
} else {
    CMplot(x1t, plot.type="m", LOG10=FALSE, col=c("grey30","grey60"), ylab="-log10(p)", # bin.range=c(0,500),
            chr.den.col=c("darkgreen", "yellow", "red"), main=main.text,
            threshold=ths, threshold.lty=c(2), threshold.lwd=c(1), threshold.col=c("black"),
            highlight=signal_topp$SNP, highlight.cex=1, 
            highlight.col=signal_topp$col, highlight.text.col=signal_topp$text_col, highlight.text=signal_topp$text,
            signal.col=c("cornflowerblue"),signal.cex=c(1),
            file="jpg",file.name=memo.text,dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)
}



############################################################################
############ Next, CREATE MANHATTAN PLOT FOR KNOCKOFF ANALYSIS #############
############################################################################

main.text <- 'Knockoff conditional independent test' # plot title
memo.text = paste('GhostKnockoffGWAS_',out_filename, sep='') # full file name

## read result file
x1 <- fread(input_file,header=T)
x1_cols = colnames(x1)
setnames(x1, grep("^chr", x1_cols)[1], "CHR")
setnames(x1, grep("^pos", x1_cols)[1], "BP")
setnames(x1, grep("^ref", x1_cols)[1], "REF")
setnames(x1, grep("^alt", x1_cols)[1], "ALT")
x1 <- x1[, SNP:=paste(CHR,BP,REF,ALT,sep=":")]
x1 <- x1[which(!is.na(x1[,'CHR'])),]
setorderv(x1, c("CHR","BP"))
x1 <- as.data.frame(x1)
x1 <- x1[match(unique(x1$SNP),x1$SNP),]

#find unique variants
x1<-x1[match(unique(x1$SNP),x1$SNP),]

############ selected variants
x1_sug <- x1[x1[,'qvals']<=target_fdr,,drop=F]
x1_sug <- x1_sug[!is.na(x1_sug[,'qvals']),,drop=F]

if(nrow(x1_sug)>0){
    x1_sug$IND <- NA
    x1_sug$TOP <- NA
    x1_sug$RAst <- NA
    x1_sug$RAen <- NA

    loc_w <- 1000000
    ra_w <- 500000
        
    # find independent loci
    CHs <- unique(x1_sug$CHR)
    inds <- 0
    for (ch in CHs) {
        BPs <- x1_sug[which(x1_sug$CHR==ch),"BP"]
        BP_first <- BPs[1]
        for (bp in BPs) {
        if (bp==BP_first) {
            inds <- inds + 1
        } else {
            bp_sw <- bp - BP_prv
            if (bp_sw > loc_w) {
            inds <- inds + 1
            }
        }
        x1_sug[which(x1_sug$CHR==ch & x1_sug$BP==bp),"IND"] <- inds
        BP_prv <- bp
        }
    }
    num_discoveries = tail(x1_sug$IND, n=1)
        
    # find top SNPs in independent loci (first by group W, then by zscores if multiple W conincide)
    # for (i in x1_sug$IND) {
    #     x1_sug_t <- x1_sug[which(x1_sug$IND==i),c("TOP","SNP","group_W","zscores")]
    #     max_W = max(x1_sug_t[,'group_W'])
    #     max_W_idx = which(x1_sug_t[,'group_W'] == max_W)
    #     max_Z = max(abs(x1_sug_t[max_W_idx,'zscores']))
    #     x1_sug_t[which(abs(x1_sug_t[,'zscores']) == max_Z),"TOP"] <- 
    #         x1_sug_t[which(abs(x1_sug_t[,'zscores']) == max_Z),"SNP"]
    #     x1_sug[which(x1_sug$IND==i),c("TOP","SNP",'group_W',"zscores")] <- x1_sug_t
    # }

    # find top SNPs in independent loci
    for (i in x1_sug$IND) {
        x1_sug_t <- x1_sug[which(x1_sug$IND==i),c("TOP","SNP","W")]
        x1_sug_t[which.max(x1_sug_t[,"W"]),"TOP"] <- 
            x1_sug_t[which.max(x1_sug_t[,"W"]),"SNP"]
        x1_sug[which(x1_sug$IND==i),c("TOP","SNP","W")] <- x1_sug_t
    }
        
    # set range +/- extension for independent loci
    for (i in x1_sug$IND) {
        x1_sug_t <- x1_sug[which(x1_sug$IND==i),c("TOP","BP","RAst","RAen")]
        x1_sug_t[,"RAst"] <- x1_sug_t[1,"BP"] - ra_w
        x1_sug_t[,"RAen"] <- x1_sug_t[dim(x1_sug_t)[1],"BP"] + ra_w
        x1_sug[which(x1_sug$IND==i),c("TOP","BP","RAst","RAen")] <- x1_sug_t
    }
        
    #flag new loci
    new_loc_w <- 1000000
    x1_sug$NEW_hit <- "Y"
    AD_loci<-x1_sug
    temp.new<-c()
    for(i in 1:nrow(x1_sug)){
        temp.chr<-x1_sug[i,'CHR']
        temp.pos<-x1_sug[i,'BP']
        temp.new<-c(temp.new,sum((temp.chr==AD_loci[,'CHR']) & (temp.pos>=AD_loci[,'BP']-new_loc_w) & (temp.pos<=AD_loci[,'BP']+new_loc_w), na.rm=T))
    }
    for(i in 1:nrow(x1_sug)){
        if(sum(temp.new[which(x1_sug[,'IND']==x1_sug[i,'IND'])])>0){x1_sug[i,"NEW_hit"] <- "N"}
    }
        
    # extract signal from all variants for respective independent loci (for plotting)
    x1_sig<-x1_sug[x1_sug[,'qvals']<=target_fdr,]
    signal <- c()
    for (t in x1_sig[which(!is.na(x1_sig$TOP)),"TOP"]) {
        ra <- x1_sig[which(x1_sig$TOP==t),c("CHR","RAst","RAen","qvals","W","IND","rsid")]
        indt <- ra$IND
        xs <- x1[which(x1$CHR==ra$CHR & (x1$BP >= ra$RAst) & (x1$BP <= ra$RAen) ),
                c("SNP","qvals","W","BP","CHR","rsid")]

        if (length(xs$SNP)>0) {
            xs <- xs[,c("SNP","W","qvals")]
            if (ra[,"qvals"]<=target_fdr) {xs$col <- "purple"}
            xs$text_col <- "blue";
            if (any(x1_sig[which(x1_sig$IND==indt),"NEW_hit"] %in% "N")) {xs$text_col <- "red"}
            xs$pch <- 19
            xs$text <- NA; xs[,"text"] <- ra$rsid
            signal <- rbind(signal,xs)
        }
    }

    # keep signal with highlights
    signal_high <- signal[which(signal[,"qvals"]<=target_fdr),]
    signal_top <- signal[which(signal$SNP %in% x1_sug[which(!is.na(x1_sug$TOP)),"TOP"]),]
    signal_topp <- signal_top

    # remove any empty gene names in signal_topp if applicable
    if (any(which(signal_topp$text %in% ""))) {
        signal_topp <- signal_topp[-which(signal_topp$text %in% ""),]
    }

    # thresholds
    ths <- min(x1[x1[,"qvals"]<=target_fdr,"W"])
    ylim=c(0,ths*8)
} else { # no discoveries
    ths <- Inf
    ylim <- c(0, 1.2 * max(x1[,"W"]))
    signal_topp <- x1_sug
}

# CMplot
x1t <- x1[,c("rsid","SNP","CHR","BP","W","W","qvals")] 
x1t[which(x1t[,"W"]>ths*8),"W"]<- ths*8
if(ths!=Inf){ths<-min(ths,max(x1t[,"W"]))}
ths <- c(ths)
signal_topp<-signal_topp[signal_topp[,"qvals"]<=0.1,]

# for debugging
# print(paste0("ths = ", ths))
# x1t[x1t$W >= ths,]

x1t <- x1t[,c("SNP","CHR","BP","W")]
colnames(x1t)[4]<-'Test statistic'

# do not label SNPs (code runs too slow) if there are too many discoveries
if (num_discoveries > 100) {
    CMplot(x1t, plot.type="m", LOG10=FALSE, col=c("grey30","grey60"), ylab="Test statistic", ylim=ylim, # bin.range=c(0,500),
        chr.den.col=c("darkgreen", "yellow", "red"), main=main.text,
        threshold=ths, threshold.lty=c(2), threshold.lwd=c(1), threshold.col=c("red"),
        highlight=signal_topp$SNP, highlight.cex=1, 
        highlight.col=signal_topp$col, highlight.text.col=signal_topp$text_col, 
        signal.col=c("cornflowerblue"),signal.cex=c(1),
        file="jpg",file.name=memo.text,dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)
} else {
    CMplot(x1t, plot.type="m", LOG10=FALSE, col=c("grey30","grey60"), ylab="Test statistic", ylim=ylim, # bin.range=c(0,500),
        chr.den.col=c("darkgreen", "yellow", "red"), main=main.text,
        threshold=ths, threshold.lty=c(2), threshold.lwd=c(1), threshold.col=c("red"),
        highlight=signal_topp$SNP, highlight.cex=1, 
        highlight.col=signal_topp$col, highlight.text.col=signal_topp$text_col, highlight.text=signal_topp$text,
        signal.col=c("cornflowerblue"),signal.cex=c(1),
        file="jpg",file.name=memo.text,dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)
}

setwd(original_dir)
print('Done!')
