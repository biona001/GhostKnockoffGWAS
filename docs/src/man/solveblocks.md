
# Customizing your own LD files

One can customize the `--LD-files` input starting with individual level data stored in VCF format. This feature is supported by the `solveblock` executable located within `GhostKnockoffGWAS/bin/solveblock`.

## Command-line documentation of `solveblock` executable

Simple run

```
solveblock --vcffile test.vcf.gz --chr 1 --start_bp 10583 --end_bp 1892607 --outdir ./test_LD_files --genome-build 19 
```

## Required inputs

| Option name              | Argument        | Description   |
| :---                     |    :----:       |   :---        |
| `--vcffile`        | String | A VCF file storing individual level genotypes. Must end in `.vcf` or `.vcf.gz`. The ALT field for each record must be unique, i.e. multiallelic records must be split first. Missing genotypes `GT` will be imputed by column mean. |
| `--chr`            | Int    | Target chromosome. This MUST be an integer and it must match the `CHROM` field in your VCF file. Thus, if your VCF file has CHROM field like `chr1`, `CHR1`, or `CHROM1` etc, each record must be renamed into `1`. |
| `--start_bp` | Int    | starting basepair (position) |
| `--end_bp` | Int    | ending basepair (position) |
| `--outdir`          | String | Directory that the output will be stored in (must exist) |
| `--genome-build` | Int | human genome build for position field of the VCF file. This must be 19 (hg19) or 38 (hg38) |

## Optional inputs

| Option name              | Argument         | Description   |
| :---                    |    :----:         |   :---     |
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

+ In our papers, we defined each start and end position by adapting the [quasi-independent regions of ldetect](https://bitbucket.org/nygcresearch/ldetect-data/src/master/). 
+ Given individual level data, one can compute approximately independent LD blocks directly, see [reference](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8696101/) and [R software](https://privefl.github.io/bigsnpr/reference/snp_ldsplit.html).

## A note on run-time

Because VCF files text files, it is inherently slow to read, even if it is indexed. Thus, we *strongly recommend* one to split the input VCF file by chromosomes, and possibly into smaller chunks, before running `solveblock`. For optimal performance, it is best to filter the VCF file down to records between `start_bp` and `end_bp` (e.g. with `bcftools`) before running `solveblock`. 
