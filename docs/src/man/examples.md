
# Detailed Example

This page collect examples of running the ghost knockoff pipeline. We will cover topics such as installation, examining input data, running the software, and interpreting the output. 

## Step 1: Download pre-processed LD files and binary executable

Proceed to the [Downloads page](https://biona001.github.io/GhostKnockoffGWAS/dev/man/download) and download (1) the software as well as (2) a pre-processed knockoff dataset suitable for your analysis. 

After downloading the software and LD data, e.g. `app_linux_x86.tar.gz` and `EUR.zip`, unzip the files in linux command line via:
```shell
tar -xvzf app_linux_x86.tar.gz
unzip EUR.zip # decompresses to ~8.7GB
```
This should create 2 folders `app_linux_x86/` and `EUR/` in the current directory. The executable is located inside `app_linux_x86/bin/GhostKnockoffGWAS`. We recommend adding the folder containing the `GhostKnockoffGWAS` executable to `PATH` for easier access.

!!! warning

    Do NOT modify the contents in unzipped folders! 

## Step 2: Prepare a valid Z score file

One needs a [valid Z score file](https://biona001.github.io/GhostKnockoffGWAS/dev/man/getting_started/#Acceptable-Z-scores-file-format) as input. 

If you would like to follow along with this tutorial, feel free to download this test data [example_zfile.txt](https://github.com/biona001/GhostKnockoffGWAS/blob/main/data/example_zfile.txt) (17MB). The first few rows is
```
CHR	POS	REF	ALT	Z
7	27916	T	C	1.82946485242
7	30580	C	T	0.877343668618
7	30581	A	T	0.876791309991
7	31273	G	C	-0.567289962351949
7	31439	T	A	-0.907002943915131
7	31627	A	C	0.577058407641
7	32858	C	T	1.80586134742
7	33482	T	G	0.47877317796
7	34215	T	C	-0.711135940901
```
+ The first row is a header row which includes `CHR`, `POS`, `REF`, `ALT`, `Z`. Other columns will be ignored. 
+ Each row is a different SNP and each column is separated by a tab (i.e. `\t` character) or a comma

In this example

+ The `POS` field corresponds to hg38 positions. GhostKnockoffGWAS requires the position to be either hg19 or hg38.
+ The sample size used for generating this data is `506200`. Thus one should specify `--N 506200`.

## Step 3: Running the analysis

To see a list of available arguments, execute `GhostKnockoffGWAS --help`. 

To run the example analysis, run the following in the terminal

```shell
GhostKnockoffGWAS --zfile example_zfile.txt --LD-files EUR --N 506200 --genome-build 38 --out example_output
```

Here is the expected output:
```
Welcome to GhostKnockoffGWAS analysis!
You have specified the following options:
zfile           = /scratch/users/bbchu/GhostKnockoffGWAS/data/example_zfile.txt
LD_files        = /scratch/users/bbchu/GhostKnockoffGWAS/data/EUR
N (sample size) = 506200
hg_build        = 38
outdir          = /scratch/users/bbchu/GhostKnockoffGWAS/data/
outfile         = /scratch/users/bbchu/GhostKnockoffGWAS/data/example_output
seed            = 2023
verbose         = true
random_shuffle  = true
skip_shrinkage_check = false

count_matchable_snps processed chr 7, cumulative SNPs = 35855
region 1 / 99 (f = LD_start100196651_end101199252.h5): chr 7, nz beta = 9, nsnps = 306, shrinkage = 0.1909
region 2 / 99 (f = LD_start101199253_end103197509.h5): chr 7, nz beta = 11, nsnps = 332, shrinkage = 0.0346
region 3 / 99 (f = LD_start103197510_end104159524.h5): chr 7, nz beta = 12, nsnps = 215, shrinkage = 0.0458
region 4 / 99 (f = LD_start104159525_end105682904.h5): chr 7, nz beta = 10, nsnps = 358, shrinkage = 0.0012
region 5 / 99 (f = LD_start105682905_end107780177.h5): chr 7, nz beta = 18, nsnps = 532, shrinkage = 0.0034
...<some output truncated>

Matched 35855 SNPs with Z-scores to the reference panel
Mean LD shrinkage = 0.020501422972314207.
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
region 1 / 99 (f = LD_start100196651_end101199252.h5): chr 7, nz beta = 9, nsnps = 306, shrinkage = 0.1909
region 2 / 99 (f = LD_start101199253_end103197509.h5): chr 7, nz beta = 11, nsnps = 332, shrinkage = 0.0346
region 3 / 99 (f = LD_start103197510_end104159524.h5): chr 7, nz beta = 12, nsnps = 215, shrinkage = 0.0458
...
```
    Here there are 99 regions in chromosome 7. For each region it prints the number of non-zero beta estimated in that region, the number of Z-scores that are present in that region, and finally the level of shrinkage. The shrinkage level is a number between 0 and 1. It quantifies how well the correlation matrices used in the analysis approximates the LD structure for the original GWAS study under the null ($z = 0$), see [SuSiE paper](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010299) equation 24 for details.
+ Finally, the program concludes by printing the number of Z scores successfully matched, the output path, as well as a rough estimate of runtime. In this simple example, the analysis finished in roughly half a minute. 

## Step 4: Interpreting the result

If you are following along, `GhostKnockoffGWAS` should have produced 2 outputs
1. `example_output_summary.txt`
2. `example_output.txt`

### `example_output_summary.txt` 

This file contains broad summary of the analysis, as shown below

```
target_fdr_0.01_num_selected,0
target_fdr_0.05_num_selected,11
target_fdr_0.1_num_selected,15
target_fdr_0.2_num_selected,19
m,5
nregions,99
nsnps,35855
lasso_lambda,0.003807185801078654
mean_LD_shrinkage,0.020501422972314207
import_time,14.982829126999999
sample_knockoff_time,8.674102106999996
ghostbasil_time,0.706785777
knockoff_filter_time,5.011729142
total_time,29.969953060150146
sample_knockoff_time_t21,4.213006547999999
sample_knockoff_time_t22,0.9331997300000002
sample_knockoff_time_t23,1.103306673
sample_knockoff_time_t24,2.3132567919999993
```

+ The first 4 rows indicate the number of unique (conditionally-independent) discoveries according to `GhostKnockoffGWAS`, for different target FDR levels. For example, when target $\text{FDR} = 0.1$, there are 15 conditionally inependent discoveries, with each discovery potentially covering >1 SNP. According to the knockoff procedure, these discoveries are conditionally independent, although one can apply a post-processing step to further count the number of independent discoveries. We will see this action later in step 5. 
+ The next few rows contain parameters used in the analysis, as well as timing results. 

!!! tip
    
    One should always check the value of `mean_LD_shrinkage`, here it is $0.02050$. As discussed above, this value quantifies how well the correlation matrices used in the analysis approximates the LD structure for the original GWAS study. A value close to 0 is good, while larger values indicate deviation. `GhostKnokcoffGWAS` automatically terminates when this value exceeds a certain threshold.

### `example_output.txt`

This is a comma-separated file that contains the full knockoff analysis output. The first 5 rows are shown:
```
$ head -5 example_output.txt
rsid,AF,chr,ref,alt,pos_hg19,pos_hg38,group,zscores,lasso_beta,variant_kappa,variant_tau,variant_W,variant_q,pvals,group_W,group_kappa,group_tau,group_qvals,selected_fdr0.01,selected_fdr0.05,selected_fdr0.1,selected_fdr0.2
rs4535687,0.15927,7,G,C,41892,41892,chr7_start16161_end972751_group1_0,-1.17940334810126,0.0,0,0.0,0.0,1.0,0.23823760256835697,0.0,0.0,0.0,1.0,0,0,0,0
rs62429406,0.031058,7,T,G,43748,43748,chr7_start16161_end972751_group2_0,0.636126444862832,0.0,0,0.0,0.0,1.0,0.5246940103826294,0.0,0.0,0.0,1.0,0,0,0,0
rs117163387,0.034958,7,C,T,43961,43961,chr7_start16161_end972751_group3_0,-0.548757491205702,0.0,0,0.0,0.0,1.0,0.5831718861307663,0.0,0.0,0.0,1.0,0,0,0,0
rs4247525,0.040199,7,T,C,44167,44167,chr7_start16161_end972751_group4_0,0.463442453535633,0.0,0,0.0,0.0,1.0,0.6430472544316368,0.0,0.0,0.0,1.0,0,0,0,0
```

The first row is a header row. Each proceeding row corresponds to a SNP that was used in the analysis. 

+ `rsid,AF,chr,ref,alt,pos_hg19,pos_hg38` is the SNP ID, alternate allele frequency, reference allele, alternate allele, basepair position in HG19 coordinates, and basepair position in HG38 coordinates.
+ `group` column: defines group membership. Note that in GhostKnockoffGWAS, false discovery rate (FDR) is guaranteed at the group level, that is, the expected number of falsely discovered groups is less than the target FDR level.
+ `zscores`: This is the user-provided Z-scores.
+ `lasso_beta`: This is the Lasso's estimated effect size for each SNP conditional on the knockoffs. 
+ `variant_kappa,variant_tau,variant_W,variant_q,pvals,group_W,group_kappa,group_tau`: these are knockoff statistics computed from the analysis, please refer to our paper for more detail. 
+ `variant_q,group_qvals`: This is the knockoff q-values, which is the minimum target FDR for a given variable to be selected, i.e. for a target FDR level $\alpha$, all variants with `group_qvals` $\le \alpha$ is selected. `GhostKnockoffGWAS` performs selection on the group-level while variant-level qvalue is used for labeling significant SNPs in downstream manhattan plots. For details, see eq 19 of [this paper](https://www.nature.com/articles/s41467-022-34932-z)
+ `pvals`: This is the p-value obtained by back-transforming the input Z-scores
+ `selected_fdr*` columns: these inform whether the variable is selected. Its values are 0 (indicating the SNP does not belong to a group that has been selected) or 1 (this SNP has been selected, along with those in the same group ).

## Step 5: Generating Manhattan plots

We can generate Manhattan plots by running [this R script](https://github.com/biona001/GhostKnockoffGWAS/blob/main/src/manhattan.R) in the terminal (this requires the `R` packages `data.table`, `plyr`, `dplyr`, `CMplot`). Usage:

```R
$ Rscript --vanilla manhattan.R arg1 arg2 arg3 arg4
```
+ `arg1`: Main output file from GhostKnockoffGWAS
+ `arg2`: Where output Manhattan plots should be stored (a `.` indicates store in current directory)
+ `arg3`: Output filename (without extensions) to be used for both plots, e.g. phenotype name
+ `arg4`: Target FDR in percentage

For example, 

```R
$ Rscript --vanilla manhattan.R example_output.txt . example_plot 0.1
```

This produced the following plots

![knockoff_manhattan](../assets/Rect_Manhtn.GhostKnockoffGWAS_chr7.jpg)
![marginal_manhattan](../assets/Rect_Manhtn.MarginalAssociationTest_chr7.jpg)

### Explanation:

+ The knockoff plot displays the knockoff W values on the y-axis, one dot for each SNP. As we are plotting the group-level W statistics, all variants within the same group possess the same W value. How then do we decide which SNP to label? In this R script, the most significant SNP (as determined by the *individual-variant-level W statistics*) within a 1Mb region is labeled and colored with purple. A light blue dot is within 1Mb region of another more significant SNP. Careful readers may recall that in the summary file (shown in step 4), knockoffs discovered 15 conditionally independent *groups*, but here only 11 SNPs were labelled. This is because some discovered groups are too close to each other. In this example, although there are 15 independent discoveries according to the knockoff methodology, there are only 11 discoveries that are physically greater than 1Mb apart. Finally, for the variant rs9640386, although it has a weaker *group-level* W score compared to a nearby variant, its individual level W statistic is stronger, and therefore it is the labeled SNP. 
+ The marginal plot is a standard Manhattan plot with the y-axis plotting the negative logged p-values. Similar to the knockoff plot, all dots above the dotted line are marginally significant and colored with light blue, while the most signicant SNP within 1Mb region is colored with purple. 
+ The color bars beneath the x-axis displays chromosome density.
