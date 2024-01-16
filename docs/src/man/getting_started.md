
# Getting started with Ghost Knockoff GWAS analysis

This tutorial is for generating *Ghost Knockoffs* for analyzing summary statistics from a genome-wide association studies (GWAS). The methodology is described in the following papers

> He Z, Liu L, Belloy ME, Le Guen Y, Sossin A, Liu X, Qi X, Ma S, Gyawali PK, Wyss-Coray T, Tang H. GhostKnockoff inference empowers identification of putative causal variants in genome-wide association studies. Nature Communications. 2022 Nov 23;13(1):7209.

The main working assumption is that we do not have access to individual level genotype or phenotype data. Rather, for each SNP, we have

1. Z-scores $Z_j$ with respect to some phenotype from a GWAS, and
2. Access to LD (linkage disequilibrium) matrix

## Typical Workflow

Most users are expected to follow this workflow. Those familiar with the Julia programming language can use GhostKnockoffGWAS as a regular julia package, see [usage within Julia](https://biona001.github.io/GhostKnockoffGWAS/dev/man/julia).

1. Download the [binary executable file]() (XXX GB)
2. Download the [pre-computed knockoff statistics](https://drive.google.com/file/d/1_ajlxFWE2MCSgBXDgDbeZh9Lq721WANA/view) (8.2GB)
3. Unzip both datasets
4. Prepare your input Z score file into accepted format, see [Acceptable Z-scores](https://biona001.github.io/GhostKnockoffGWAS/dev/man/getting_started/#Acceptable-Z-scores-format) below. 
5. Run the executable, see [running the executable](https://biona001.github.io/GhostKnockoffGWAS/dev/man/getting_started/#Running-the-executable)

## Running the executable

```shell
usage: <PROGRAM> --zfile ZFILE --knockoff-dir KNOCKOFF-DIR --N N
                 --genome-build GENOME-BUILD --out OUT [--seed SEED]
                 [-h]

  --zfile ZFILE         Tab or comma separated summary Z-score file,
                        which can be .gz compressed. The first row
                        must be a header line that contains at least
                        CHR, POS, REF, ALT, and Z (other columns will
                        be ignored). Each row should be a SNP. CHR is
                        the chromosome column and must be integer
                        valued (e.g. chr22 is NOT valid). POS is the
                        position of the SNP (aligned to HG19 or HG38).
                        REF and ALT are the reference and alternate
                        alleles, which will be treated as the
                        non-effective and effect alleles,
                        respectively. Finally Z is the Z-score column
  --knockoff-dir KNOCKOFF-DIR
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
  -h, --help            show this help message and exit
```

Example run:
```shell
./GhostKnockoffGWAS \
    --zfile ../../data/AD_Zscores_Meta_modified.txt \
    --knockoff-dir ../../data/EUR \
    --N 506200 \
    --genome-build 38 \
    --out ../../data/test_alzheimers_meta
```

## Acceptable Z-scores file format

The Z score file should satisfy the following requirements:
1. It is a comma- or tab-separated file (.gz compressed is acceptable)
2. The first row should be a header line, and every row after the first will be treated as a different SNP. 
3. The header line should include `CHR`, `POS`, `REF`, `ALT`, and `Z`. The `ALT` allele will be treated as the effect allele and `REF` be treated as non-effect allele. The POS (position) field of each variant must be from HG19 or HG38, which must be specified by the `--genome-build` argument. CHR/POS/REF/ALT fields cannot have missing values. Missing Z scores can be specified as `NaN` or as an empty cell.

If you have p-values, effect sizes, odds ratios...etc but not Z scores, you can convert them into Z score, for example by following the **Notes on computing Z-scores** of [this blog post](https://huwenboshi.github.io/data%20management/2017/11/23/tips-for-formatting-gwas-summary-stats.html). 
