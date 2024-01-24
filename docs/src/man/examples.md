
# Detailed Example

This page collect examples of running the ghost knockoff pipeline. We will cover topics such as installation, examining input data, running the software, and interpreting the output. 

## Step 0: Download example data

If you would like to follow along with this tutorial, feel free to download [example_zfile.txt]() (4MB) which contains ~200k Z scores from chromosome 17. The first few rows of this file is shown below

```
$ head example_zfile.txt

CHR	POS	REF	ALT	Z
17	150509	T	TA	1.08773561923134
17	151035	T	C	0.703898767202681
17	151041	G	A	1.10771707088118
17	151872	T	C	-0.299877259561085
17	152087	C	T	-0.371627135786605
17	152104	G	A	-0.28387322965385
17	152248	G	A	0.901618600934489
17	152427	G	A	1.10987516000804
17	152771	A	G	0.708492545266136
```

In this example data, 
+ The `POS` field corresponds to hg38 positions. GhostKnockoffGWAS requires the position to be either hg19 or hg38.
+ The sample size used for generating this data is `506200`. Thus one should specify `--N 506200`.

!!! note 

    See [acceptable Z score format](https://biona001.github.io/GhostKnockoffGWAS/dev/man/getting_started/#Acceptable-Z-scores-file-format) for detailed requirements on this file.

## Step 1: Download Knockoff files and binary executable

Proceed to the [Downloads page]() and download (1) the software as well as (2) a pre-processed knockoff dataset suitable for your analysis. After unzipping, the executable will be located inside `bin/GhostKnockoffGWAS`. We recommend adding the folder containing the `GhostKnockoffGWAS` executable to `PATH` for easier access.

## Step 2: Prepare a valid Z score file

One needs a [valid Z score file](https://biona001.github.io/GhostKnockoffGWAS/dev/man/getting_started/#Acceptable-Z-scores-file-format) as input. Here is an example:

```
$ head example_zfile.txt

CHR	POS	REF	ALT	Z
17	150509	T	TA	1.08773561923134
17	151035	T	C	0.703898767202681
17	151041	G	A	1.10771707088118
17	151872	T	C	-0.299877259561085
17	152087	C	T	-0.371627135786605
17	152104	G	A	-0.28387322965385
17	152248	G	A	0.901618600934489
17	152427	G	A	1.10987516000804
17	152771	A	G	0.708492545266136
```
+ The first row is a header row which includes `CHR`, `POS`, `REF`, `ALT`, `Z`. Other columns will be ignored. 
+ Each row is a different SNP and each column is separated by a tab (i.e. `\t` character) or a comma

## Step 3: Running the analysis

To see a list of available arguments, execute `GhostKnockoffGWAS --help`. 

To run the example analysis, one can do

```shell
$ GhostKnockoffGWAS \
    --zfile example_zfile.txt \
    --knockoff-dir EUR \
    --N 506200 \
    --genome-build 38 \
    --out example_output
```

Here is the expected output:
```
Welcome to GhostKnockoffGWAS analysis!
You have specified the following options:
zfile           = /scratch/users/bbchu/GhostKnockoffGWAS/data/example_zfile.txt
knockoff_dir    = /scratch/users/bbchu/GhostKnockoffGWAS/data/EUR
N (sample size) = 506200
hg_build        = 38
outdir          = /scratch/users/bbchu/GhostKnockoffGWAS/data/
outfile         = /scratch/users/bbchu/GhostKnockoffGWAS/data/example_output
seed            = 2023
verbose         = true


count_matchable_snps processed chr 17, cumulative SNPs = 21136
region 1 / 47: chr 17, nz beta = 10, nsnps = 425, shrinkage = 0.0043
region 2 / 47: chr 17, nz beta = 18, nsnps = 319, shrinkage = 0.0326
region 3 / 47: chr 17, nz beta = 1, nsnps = 232, shrinkage = 0.0135
region 4 / 47: chr 17, nz beta = 10, nsnps = 303, shrinkage = 0.0
region 5 / 47: chr 17, nz beta = 7, nsnps = 401, shrinkage = 0.0331
region 6 / 47: chr 17, nz beta = 3, nsnps = 285, shrinkage = 0.0235
region 7 / 47: chr 17, nz beta = 14, nsnps = 453, shrinkage = 0.0855
region 8 / 47: chr 17, nz beta = 12, nsnps = 385, shrinkage = 0.0256
region 9 / 47: chr 17, nz beta = 15, nsnps = 584, shrinkage = 0.0077
region 10 / 47: chr 17, nz beta = 8, nsnps = 490, shrinkage = 0.0234
region 11 / 47: chr 17, nz beta = 17, nsnps = 320, shrinkage = 0.0215
region 12 / 47: chr 17, nz beta = 8, nsnps = 424, shrinkage = 0.0102
region 13 / 47: chr 17, nz beta = 12, nsnps = 492, shrinkage = 0.0151
region 14 / 47: chr 17, nz beta = 14, nsnps = 438, shrinkage = 0.0275
region 15 / 47: chr 17, nz beta = 12, nsnps = 522, shrinkage = 0.0042
region 16 / 47: chr 17, nz beta = 17, nsnps = 450, shrinkage = 0.0224
region 17 / 47: chr 17, nz beta = 10, nsnps = 386, shrinkage = 0.0038
region 18 / 47: chr 17, nz beta = 6, nsnps = 327, shrinkage = 0.0098
region 19 / 47: chr 17, nz beta = 12, nsnps = 333, shrinkage = 0.008
region 20 / 47: chr 17, nz beta = 15, nsnps = 257, shrinkage = 0.0101
region 21 / 47: chr 17, nz beta = 15, nsnps = 773, shrinkage = 0.3052
region 22 / 47: chr 17, nz beta = 21, nsnps = 364, shrinkage = 0.0269
region 23 / 47: chr 17, nz beta = 21, nsnps = 332, shrinkage = 0.0147
region 24 / 47: chr 17, nz beta = 31, nsnps = 656, shrinkage = 0.0211
region 25 / 47: chr 17, nz beta = 2, nsnps = 123, shrinkage = 0.0101
region 26 / 47: chr 17, nz beta = 3, nsnps = 177, shrinkage = 0.0109
region 27 / 47: chr 17, nz beta = 7, nsnps = 375, shrinkage = 0.0122
region 28 / 47: chr 17, nz beta = 13, nsnps = 472, shrinkage = 0.0437
region 29 / 47: chr 17, nz beta = 17, nsnps = 538, shrinkage = 0.0318
region 30 / 47: chr 17, nz beta = 12, nsnps = 555, shrinkage = 0.0184
region 31 / 47: chr 17, nz beta = 10, nsnps = 642, shrinkage = 0.0048
region 32 / 47: chr 17, nz beta = 4, nsnps = 286, shrinkage = 0.0123
region 33 / 47: chr 17, nz beta = 13, nsnps = 403, shrinkage = 0.0249
region 34 / 47: chr 17, nz beta = 5, nsnps = 278, shrinkage = 0.0082
region 35 / 47: chr 17, nz beta = 11, nsnps = 483, shrinkage = 0.0276
region 36 / 47: chr 17, nz beta = 35, nsnps = 800, shrinkage = 0.0181
region 37 / 47: chr 17, nz beta = 7, nsnps = 358, shrinkage = 0.0054
region 38 / 47: chr 17, nz beta = 32, nsnps = 1278, shrinkage = 0.004
region 39 / 47: chr 17, nz beta = 14, nsnps = 453, shrinkage = 0.0105
region 40 / 47: chr 17, nz beta = 10, nsnps = 341, shrinkage = 0.005
region 41 / 47: chr 17, nz beta = 14, nsnps = 731, shrinkage = 0.0104
region 42 / 47: chr 17, nz beta = 15, nsnps = 545, shrinkage = 0.0463
region 43 / 47: chr 17, nz beta = 7, nsnps = 575, shrinkage = 0.0093
region 44 / 47: chr 17, nz beta = 4, nsnps = 385, shrinkage = 0.0141
region 45 / 47: chr 17, nz beta = 9, nsnps = 354, shrinkage = 0.0083
region 46 / 47: chr 17, nz beta = 12, nsnps = 615, shrinkage = 0.0108
region 47 / 47: chr 17, nz beta = 7, nsnps = 418, shrinkage = 0.0251
Matched 21136 SNPs with Z-scores to the reference panel
Mean LD shrinkage = 0.02387362024098241.
Done! Result saved to /scratch/users/bbchu/GhostKnockoffGWAS/data/example_output. 
Overall runtime = 34.12649257 seconds, with 
   1.456621308 seconds spent on reading the Z score file
   32.669871262 seconds spent on doing the analysis
```

**Explanation for intermediate outputs**:

+ `GhostKnockoffGWAS` first prints the user-specified parameters in the analysis. Verify that they are correct.
+ Next we print the output of `count_matchable_snps`. It is essentially matching user supplied Z scores to the pre-computed knockoff data and counting how many SNPs can be matched. This information will be used to quantify the level shrinkage in Lasso regression. 
+ Then for each region, it will try to analyze the genome in quasi-independent regions, e.g. 
```
region 1 / 47: chr 17, nz beta = 10, nsnps = 425, shrinkage = 0.0043
region 2 / 47: chr 17, nz beta = 18, nsnps = 319, shrinkage = 0.0326
region 3 / 47: chr 17, nz beta = 1, nsnps = 232, shrinkage = 0.0135
...
```
Here there are 47 regions in chromosome 17. For each region it prints the number of non-zero beta estimated in that region, the number of Z-scores that are present in that region, and finally the level of shrinkage. The shrinkage level is a number between 0 and 1. It quantifies how well the correlation matrices used in the analysis approximates the LD structure for the original GWAS study under the null ($z = 0$), see [SuSiE paper](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010299) equation 24 for details.
+ Finally, the program concludes by printing the number of Z scores successfully matched, the output path, as well as a rough estimate of runtime. In this simple example, the analysis finished in roughly half a minute. 

## Step 4: Interpreting the result

If you are following along, `GhostKnockoffGWAS` should have produced 2 outputs
1. `example_output_summary.txt`
2. `example_output.txt`

### `example_output_summary.txt` 

This file contains broad summary of the analysis, as shown below

```
target_fdr_0.01,Inf
target_fdr_0.01_num_selected,0
target_fdr_0.02,Inf
target_fdr_0.02_num_selected,0
target_fdr_0.03,Inf
target_fdr_0.03_num_selected,0
target_fdr_0.04,Inf
target_fdr_0.04_num_selected,0
target_fdr_0.05,Inf
target_fdr_0.05_num_selected,0
target_fdr_0.06,Inf
target_fdr_0.06_num_selected,0
target_fdr_0.07,Inf
target_fdr_0.07_num_selected,0
target_fdr_0.08,Inf
target_fdr_0.08_num_selected,0
target_fdr_0.09,Inf
target_fdr_0.09_num_selected,0
target_fdr_0.1,0.0038835278195177795
target_fdr_0.1_num_selected,4
target_fdr_0.11,0.0037439278466356142
target_fdr_0.11_num_selected,6
target_fdr_0.12,0.0037439278466356142
target_fdr_0.12_num_selected,6
target_fdr_0.13,0.0037439278466356142
target_fdr_0.13_num_selected,6
target_fdr_0.14,0.0031012282968098537
target_fdr_0.14_num_selected,6
target_fdr_0.15,0.0031012282968098537
target_fdr_0.15_num_selected,6
target_fdr_0.16,0.0031012282968098537
target_fdr_0.16_num_selected,6
target_fdr_0.17,0.002768189835339677
target_fdr_0.17_num_selected,6
target_fdr_0.18,0.002768189835339677
target_fdr_0.18_num_selected,6
target_fdr_0.19,0.002768189835339677
target_fdr_0.19_num_selected,6
target_fdr_0.2,0.0022878815488677467
target_fdr_0.2_num_selected,7
m,5
nregions,47
nsnps,21136
lasso_lambda,0.0035656428122281454
mean_LD_shrinkage,0.02387362024098241
import_time,2.9585762090000007
sample_knockoff_time,8.064913196999997
ghostbasil_time,0.461628368
knockoff_filter_time,2.326016636
total_time,25.699404001235962
sample_knockoff_time_t21,2.9308446070000005
sample_knockoff_time_t22,1.398813931
sample_knockoff_time_t23,0.9110162740000001
sample_knockoff_time_t24,2.6240137139999993
```

+ The first ~40 rows indicate the number of unique discoveries given by `GhostKnockoffGWAS`, for different target FDR levels. For example, when target $\text{FDR} = 0.1$, there are 4 unique discoveries and the knockoff threshold is $\hat{\tau} = 0.00388$. According to the knockoff procedure, these discoveries are conditionally independent, although one can apply a post-processing step to further count the number of independent discoveries. 
+ The next few rows contain parameters used in the analysis, as well as timing results. 
+ *The most important parameter* corresponds to the value of `mean_LD_shrinkage`, here it is $0.023873$. As discussed above, this value quantifies how well the correlation matrices used in the analysis approximates the LD structure for the original GWAS study. A value close to 0 is good, while larger values indicate deviation. `GhostKnokcoffGWAS` automatically terminates when this value exceeds a certain threshold.

### `example_output.txt`

This is a comma-separated file that contains the full knockoff analysis output. The first 5 rows are shown:
```
$ head -5 example_output.txt
rsid,AF,chr,ref,alt,pos_hg19,pos_hg38,group,zscores,lasso_beta,W,kappa,tau,pvals,selected_fdr0.01,selected_fdr0.02,selected_fdr0.03,selected_fdr0.04,selected_fdr0.05,selected_fdr0.06,selected_fdr0.07,selected_fdr0.08,selected_fdr0.09,selected_fdr0.1,selected_fdr0.11,selected_fdr0.12,selected_fdr0.13,selected_fdr0.14,selected_fdr0.15,selected_fdr0.16,selected_fdr0.17,selected_fdr0.18,selected_fdr0.19,selected_fdr0.2
rs2294076,0.74445,17,A,G,6157,156366,chr17_start56_end1172398_group1_0,0.398036369392557,0.0,0.0,0.0,0.0,0.6906033770337676,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
rs2396789,0.90858,17,T,C,8547,158756,chr17_start56_end1172398_group1_0,0.405389455204145,0.0,0.0,0.0,0.0,0.6851912607868653,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
rs36068254,0.76068,17,A,G,10583,160792,chr17_start56_end1172398_group2_0,0.454143632070033,0.0,0.0,0.0,0.0,0.6497254507764232,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
rs6565705,0.7967,17,G,A,13905,164114,chr17_start56_end1172398_group1_0,0.654603172923188,0.0,0.0,0.0,0.0,0.5127232807057284,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
```

The first row is a header row. Each proceeding row corresponds to a SNP that was used in the analysis. 

+ `selected_fdr*` columns: these inform whether the variable is selected. Its values are 0 (indicating the SNP does not belong to a group that has been selected) or 1 (this SNP has been selected, along with those in the same group ).
+ `group` column: defines group membership. Note that in GhostKnockoffGWAS, false discovery rate (FDR) is guaranteed at the group level, that is, the expected number of falsely discovered groups is less than the target FDR level.
+ `AF` column: stands for alternate-allele-frequency. 
+ `lasso_beta`: This is the Lasso's estimated effect size for each SNP conditional on the knockoffs. 
+ `W, kappa, tau`: these are knockoff statistics computed from the analysis, please refer to our paper for more detail. 
+ Other columns such as `rsid,chr,ref`...etc should be self-explanatory.
