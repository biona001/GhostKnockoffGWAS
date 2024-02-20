
# Getting started with Ghost Knockoff GWAS analysis

This package conducts knockoff-based inference to perform genome-wide conditional independent tests based on GWAS summary statistics (e.g. p-values). The methodology is described in the following paper

> He Z, Chu BB, Yang J, Gu J, Chen Z, Liu L, Morrison T, Bellow M, Qi X, Hejazi N, Mathur M, Le Guen Y, Tang H, Hastie T, Ionita-laza, I, Sabatti C, Candes C. "In silico identification of putative causal genetic variants", bioRxiv 2024. 

The main working assumption is that we do not have access to individual level genotype or phenotype data. Rather, for each SNP, we have its Z-scores with respect to some phenotype from a GWAS, and access to LD (linkage disequilibrium) data. The user is expected supply the Z-scores, while we supply the LD data in addition to some pre-computed knockoff data.

## Q: When should I use GhostKnockoffGWAS?

Answer: If you already conducted a GWAS, have an output file that includes Z scores (or equivalent) for each SNP, and there exist pre-processed LD files in [downloads page](https://biona001.github.io/GhostKnockoffGWAS/dev/man/download/) in which the listed population matches the ethnicities for your original GWAS study.

+ If your original study had little (e.g. <5) discoveries, then `GhostKnockoffGWAS` may not give better results. The methodology works better for more polygenic traits. 
+ If your study subjects are somewhat admixed, one can try using the most suitable LD files, and check how much deviation there are from the LD files by examining the `LD_shrinkage` parameter in the output of `GhostKnockoffGWAS`, see [this FAQ](https://biona001.github.io/GhostKnockoffGWAS/dev/man/FAQ/#Is-the-result-is-trustworthy?).
+ If instead you have individual level genotypes, you should run a GWAS using standard tools (e.g. PLINK, BOLT, GCTA, SAIGE, GEMMA, ...etc) before running GhostKnockoffGWAS. 

## Typical Workflow

Most users are expected to follow this workflow. Those familiar with the Julia programming language can use GhostKnockoffGWAS as a regular julia package, see [usage within Julia](https://biona001.github.io/GhostKnockoffGWAS/dev/man/julia).

1. Go to [Download Page](https://biona001.github.io/GhostKnockoffGWAS/dev/man/download) and download (1) the software and (2) the pre-processed LD files.
3. Unzip them
4. Prepare your input Z score file into accepted format, see [Acceptable Z-scores](https://biona001.github.io/GhostKnockoffGWAS/dev/man/getting_started/#Acceptable-Z-scores-file-format) below. 
5. Run the executable, see [running the executable](https://biona001.github.io/GhostKnockoffGWAS/dev/man/getting_started/#Running-the-executable)

For a detailed example, see [Detailed Example](https://biona001.github.io/GhostKnockoffGWAS/dev/man/examples/)
