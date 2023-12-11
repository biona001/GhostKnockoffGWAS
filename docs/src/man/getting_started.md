
# Getting started with Ghost Knockoff GWAS analysis

This tutorial is for generating *Ghost Knockoffs* for analyzing summary statistics from a genome-wide association studies (GWAS). The methodology is described in the following papers

> He Z, Liu L, Belloy ME, Le Guen Y, Sossin A, Liu X, Qi X, Ma S, Gyawali PK, Wyss-Coray T, Tang H. GhostKnockoff inference empowers identification of putative causal variants in genome-wide association studies. Nature Communications. 2022 Nov 23;13(1):7209.

The main working assumption is that we do not have access to individual level genotype or phenotype data. Rather, for each SNP, we have

1. Z-scores $Z_j$ with respect to some phenotype from a GWAS, and
2. Access to LD (linkage disequilibrium) matrix

## Typical Workflow

Most users are expected to follow this workflow. For advanced users, see [Developer documentation](https://biona001.github.io/GhostKnockoffGWAS/dev/man/developer).

1. Download the [binary executable file]() (XXX GB)
2. Download the [pre-computed knockoff statistics](https://drive.google.com/file/d/1_ajlxFWE2MCSgBXDgDbeZh9Lq721WANA/view) (8.2GB)
3. Unzip both datasets
4. Prepare your input Z score file into accepted format, see [Acceptable Z-scores](https://biona001.github.io/GhostKnockoffGWAS/dev/man/getting_started/#Acceptable-Z-scores-format) below. 
5. Run the executable, see [running the executable](https://biona001.github.io/GhostKnockoffGWAS/dev/man/getting_started/#Running-the-executable)

## Running the executable

To see required inputs (and optional inputs), invoke
```shell
./GhostKnockoffGWAS --help
```

Example run:
```shell
./GhostKnockoffGWAS \
    ../../data/AD_Zscores_Meta_modified.txt \
    ../../data/EUR \
    506200 \
    ../../data/test_alzheimers_meta
```

## Acceptable Z-scores format

The Z score file should satisfy the following requirements:
1. It is a comma- or tab-separated file
2. The first row should be a header line
3. The header line should include `CHR`, `POS`, `REF`, `ALT`. The `ALT` allele will be treated as the effect allele and `REF` be treated as non-effect allele. The position of each variant must be from HG19 or HG38.
4. Z scores can be given in 3 ways
    1. If a column with the header `Z` exist, then it will be used as Z scores.
    2. If `pvalue` and `beta` (effect size) columns both exist, we will transform them into Z scores
    3. If `OR` (odds ratio) and `SE` (standard error) columns exist, they will be transformed into Z scores

