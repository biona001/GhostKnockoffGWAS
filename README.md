# GhostKnockoffGWAS

| **Documentation** | **Build Status** | **Code Coverage**  |
|-------------------|------------------|--------------------|
| [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://biona001.github.io/GhostKnockoffGWAS/dev/)| [![build Actions Status](https://github.com/biona001/GhostKnockoffGWAS/actions/workflows/CI.yml/badge.svg)](https://github.com/biona001/GhostKnockoffGWAS/actions) [![CI (Julia nightly)](https://github.com/biona001/GhostKnockoffGWAS/actions/workflows/JuliaNightly.yml/badge.svg)](https://github.com/biona001/GhostKnockoffGWAS.jl/actions/workflows/JuliaNightly.yml) | [![codecov](https://codecov.io/gh/biona001/GhostKnockoffGWAS/branch/main/graph/badge.svg)](https://codecov.io/gh/biona001/GhostKnockoffGWAS) |

This is a package for analyzing summary statistics data from [*genome-wide association studies (GWAS)*](https://en.wikipedia.org/wiki/Genome-wide_association_study). The main working assumption is that we do not have access to individual level genotype or phenotype data. Rather, for each SNP, we have its Z-scores with respect to some phenotype from a GWAS, and access to LD (linkage disequilibrium) data. The user is expected supply the Z-scores, while we supply the LD data. `GhostKnockoffGWAS` integrates numerous recent advances in the knockoffs literature to control the false discovery rate of variable selection, while maximizing power, statility, and its ability to prioritize causal variants. 

## New users

To get started, please refer to the [documentation](https://biona001.github.io/GhostKnockoffGWAS/dev). 

## Advantages/disadvantages of GhostKnockoffGWAS

Compared to existing knockoff methods for GWAS, the main advantages of GhostKnockoffGWAS is (1) its ease of use and (2) its computational efficiency. The only user-provided input is marginal Z-scores. Computationally, running a knockoff-based GWAS pipeline took approximately 15 minutes on 650,000 SNPs. The main limitation of GhostKnockoffGWAS is that it relies on the availability of pre-processed LD files suitable for the user's target samples. 

## Bug fixes and user support

If you encounter a bug or need user support, please open a new issue on Github. Please provide as much detail as possible for bug reports, ideally a sequence of reproducible code that lead to the error.

PRs and feature requests are welcomed!

## Citation

If you use `GhostKnockoffGWAS` in your research, please cite the following references:

> He Z, Chu BB, Yang J, Gu J, Chen Z, Liu L, Morrison T, Bellow M, Qi X, Hejazi N, Mathur M, Le Guen Y, Tang H, Hastie T, Ionita-laza, I, Sabatti C, Candes C. "In silico identification of putative causal genetic variants", bioRxiv 2024.
