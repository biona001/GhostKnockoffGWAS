
# Getting started with Ghost Knockoff GWAS analysis

This package is conducts knockoff-based inference to perform genome-wide conditional independent tests based on summary statistics (e.g. p-values). The methodology is described in the following papers

> He Z, Chu BB, Yang J, Gu J, Chen Z, Liu L, Morrison T, Bellow M, Qi X, Hejazi N, Mathur M, Le Guen Y, Tang H, Hastie T, Ionita-laza, I, Sabatti C, Candes C. "In silico identification of putative causal genetic variants", bioRxiv 2024. 

The main working assumption is that we do not have access to individual level genotype or phenotype data. Rather, for each SNP, we have its Z-scores with respect to some phenotype from a GWAS, and access to LD (linkage disequilibrium) data. The user is expected supply the Z-scores, while we supply the LD data in addition to some pre-computed knockoff data.

## Q: When should I use GhostKnockoffGWAS?

Answer: If you already conducted a GWAS, have an output file that includes Z scores (or equivalent) for each SNP, and there exist pre-processed LD files in [downloads page](https://biona001.github.io/GhostKnockoffGWAS/dev/man/download/) in which the listed population matches the ethnicities for your original GWAS study.

+ If your original study had little (e.g. <5) discoveries, then `GhostKnockoffGWAS` may not give better results. The methodology works better for more polygenic traits. 
+ If your study subjects are somewhat admixed, one can try using the most suitable LD files, and check how much deviation there are from the LD files by examining the `LD_shrinkage` parameter in the output of `GhostKnockoffGWAS`, see [this FAQ](https://biona001.github.io/GhostKnockoffGWAS/dev/man/FAQ/#Is-the-result-is-trustworthy?).
+ If instead you have individual level genotypes, you should run a GWAS using standard tools (e.g. PLINK, BOLT, GCTA, SAIGE, GEMMA, ...etc) before running GhostKnockoffGWAS. 

## Typical Workflow

Most users are expected to follow this workflow. Those familiar with the Julia programming language can use GhostKnockoffGWAS as a regular julia package, see [usage within Julia](https://biona001.github.io/GhostKnockoffGWAS/dev/man/julia).

1. Go to [Download Page](https://biona001.github.io/GhostKnockoffGWAS/dev/man/download) and download (1) the binary executable file and (2) the pre-processed LD files.
3. Unzip them
4. Prepare your input Z score file into accepted format, see [Acceptable Z-scores](https://biona001.github.io/GhostKnockoffGWAS/dev/man/getting_started/#Acceptable-Z-scores-file-format) below. 
5. Run the executable, see [running the executable](https://biona001.github.io/GhostKnockoffGWAS/dev/man/getting_started/#Running-the-executable)

For a detailed example, see [Detailed Example](https://biona001.github.io/GhostKnockoffGWAS/dev/man/examples/)

## Running the executable

See [Detailed Examples](https://biona001.github.io/GhostKnockoffGWAS/dev/man/examples) for a analysis example. To see a list of available arguments, execute `GhostKnockoffGWAS --help`. Its output is:

```shell
usage: <PROGRAM> --zfile ZFILE --LD-files LD-FILES --N N
                 --genome-build GENOME-BUILD --out OUT [--seed SEED]
                 [--verbose VERBOSE]
                 [--skip_shrinkage_check SKIP_SHRINKAGE_CHECK] [-h]

optional arguments:
  --zfile ZFILE         Tab or comma separated summary Z-score file,
                        which can be .gz compressed. The first row
                        must be a header line that contains at least
                        CHR, POS, REF, ALT, and Z (other columns will
                        be ignored). Each row should be a SNP. CHR is
                        the chromosome column and must be integer
                        valued (e.g. chr22, , sex chromosomes, and
                        missing values are NOT valid). POS is the SNP
                        (aligned to HG19 or HG38) and cannot be
                        missing. REF the position of and ALT are the
                        reference and alternate alleles, which will be
                        treated as the non-effective and effect
                        alleles, respectively, and also cannot be
                        missing. Finally, Z is the Z-score column.
                        Missing Z scores can be specified as NaN or as
                        an empty cell.
  --LD-files LD-FILES
                        Path to the directory storing pre-processed
                        knockoff files
  --N N                 Sample size for target (original) study (type:
                        Int64)
  --genome-build GENOME-BUILD
                        Specifies the human genome build for the
                        target (original) study. Must be 19 (hg19) or
                        38 (hg38). (type: Int64)
  --out OUT             Output file prefix (without extensions)
  --seed SEED           Sets the random seed (type: Int64, default:
                        2023)
  --verbose VERBOSE     Whether to print intermediate messages (type:
                        Bool, default: true)
  --skip_shrinkage_check SKIP_SHRINKAGE_CHECK
                        Whether to allow Knockoff analysis to proceed
                        even with large (>0.25) LD shrinkages. Only
                        use this option if you know what you are
                        doing.  (type: Bool, default: true)
  -h, --help            show this help message and exit
```
