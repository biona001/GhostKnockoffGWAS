# Customizing your own LD files

One can customize the `--LD-files` input starting with individual level data stored in VCF or binary PLINK format. This feature is supported by the `solveblock` executable located within `GhostKnockoffGWAS/bin/solveblock`.

## Command-line documentation of `solveblock` executable

Simple run

```
# VCF format
solveblock --file test.vcf.gz --chr 1 --start_bp 10583 --end_bp 1892607 --outdir ./test_LD_files --genome-build 19 

# PLINK format
solveblock --file test.bed --chr 1 --start_bp 10583 --end_bp 1892607 --outdir ./test_LD_files --genome-build 19 
```

## Required inputs

| Option name              | Argument        | Description   |
| :---                     |    :----:       |   :---        |
| `--file`        | String | A VCF or binary PLINK file storing individual level genotypes. Must end in `.vcf`, `.vcf.gz`, or `.bed`. If a VCF file is used, the ALT field for each record must be unique, i.e. multiallelic records must be split first. Missing genotypes will be imputed by column mean.  |
| `--chr`            | Int    | Target chromosome. This MUST be an integer and it must match the `CHROM` field in your VCF/PLINK file. For example, if your VCF file has CHROM field like `chr1`, `CHR1`, or `CHROM1` etc, they must be renamed into `1`.  |
| `--start_bp` | Int    | starting basepair (position) |
| `--end_bp` | Int    | ending basepair (position) |
| `--outdir`          | String | Directory that the output will be stored in (must exist) |
| `--genome-build` | Int | human genome build for position of each SNP, must be 19 (hg19) or 38 (hg38) |

## Optional inputs

| Option name              | Argument         | Description   |
| :---                    |    :----:         |   :---     |
| `--covfile`        | String | An optional comma- or tab-separated file containing sample covariates (e.g. sex, age, PCs). `.gz` compressed text file is acceptable. These will be used to improve LD estimation. The first row should be a header row. The first column should be sample IDs (not necessary to be in the sample order as genotype files) and all other columns will be used as additional covariates. Note if genotypes are stored in binary PLINK format, then the sample ID column in the covariate file should be FID_IID (that is, the first 2 columns of the .fam file merged by an underscore) (default `""`) |
| `--tol`        | Float64 | Convergence tolerlance for group knockoff coordinate descent optimization (default `0.0001`) |
| `--min_maf`    | Float64 | Minimum minor allele frequency for a variable to be considered (default `0.01`) |
| `--min_hwe`    | Int     | Cutoff for hardy-weinburg equilibrium p-values. Only SNPs with p-value >= `min_hwe` will be included (default `0.0`) |
| `--method`     | String  | group knockoff optimization algorithm, choices include `maxent` (defualt), `mvr`, `sdp`, or `equi`. See sec 2 of https://arxiv.org/abs/2310.15069 |
| `--linkage`    | String  | Linkage function to use for defining group membership. It defines how the distances between features are aggregated into the distances between groups. Valid choices include `average` (default), `single`, `complete`, `ward`, and `ward_presquared`. Note if `force_contiguous=true`, `linkage` must be `:single`|
| `--force_contiguous` | Bool     | whether to force groups to be contiguous (default `false`). Note if `force_contiguous=true`, `linkage` must be `:single`) |
| `--group_cor_cutoff`    | Float64    | correlation cutoff value for defining groups (default `0.5`). Value should be between 0 and 1, where larger values correspond to larger groups. |
| `--group_rep_cutoff` | Float64 | cutoff value for selecting group-representatives (default `0.5`). Value should be between 0 and 1, where larger values correspond to more representatives per group.  |
| `--verbose` | Bool | Whether to print intermediate messages (default `true`) |


## Output format

Calling `solveblock` will create 3 files in the directory `outdir/chr` (the `chr` directory will be created if it doesn't exist, but `outdir` must exist):
1. `XXX.h5`:  This contains the following data for region XXX
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
2. `Info_XXX.csv`: This includes information for each variant (chr/pos/etc) present 
    in the corresponding `.h5` file.
3. `summary_XXX.csv`: Summary file for the knockoff optimization problem

## Determining `start_bp` and `end_bp`

There are 2 options:

1. One can defined each start and end position by leveraging existing quasi-independent regions for your target sample. For example, we previously used the European blocks of [ldetect](https://bitbucket.org/nygcresearch/ldetect-data/src/master/). 
2. Given individual level data, one can compute approximately independent LD blocks directly, see [reference](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8696101/) and [R software](https://privefl.github.io/bigsnpr/reference/snp_ldsplit.html).

For option 2, we provide an [R script](https://github.com/biona001/GhostKnockoffGWAS/blob/main/src/ld_split.R) which can be ran in the terminal (this requires the `R` packages `bigsnpr` and `dplyr`). Usage:

```R
$ Rscript --vanilla ld_split.R arg1 arg2 arg3 arg4 arg5 arg6
```
+ `arg1` = chromosome number (must be an integer)
+ `arg2` = path to PLINK binary file (must end in `.bed` extension)
+ `arg3` = path to FBM file (without extensions. If this file doesn't exist, it will be generated)
+ `arg4` = path to output file 
+ `arg5` = `thr_r2`, this is the `thr_r2` used by snp_ldsplit. All correlation smaller than `thr_r2` are set to 0
+ `arg6` = `max_r2`, this is the `max_r2` used by snp_ldsplit. This is the maximum acceptable correlation for SNPs in different blocks. 

For example, 

```R
$ Rscript --vanilla ld_split.R 1 my_plink.bed my_plink_fbm regions.txt 0.01 0.3
```

## A note on run-time

Because VCF files are plain text files, it is inherently slow to read even if it is indexed. Thus, we recommend one to convert VCFs to binary PLINK format via [PLINK 1.9](https://www.cog-genomics.org/plink/):

```
$plink_exe --vcf $vcffile --double-id --keep-allele-order --real-ref-alleles --make-bed --out $plinkprefix
```
Note the `--keep-allele-order` is crucial to prevent PLINK from randomly converting the minor allele into A1. 
