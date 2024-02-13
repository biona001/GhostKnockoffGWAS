
# Acceptable Z-scores file format

The Z score file should satisfy the following requirements:
1. It is a comma- or tab-separated text file (.gz compressed is acceptable)
2. The first row should be a header line, and every row after the first will be treated as a different SNP. 
3. By default `GhostKnockoffGWAS` will search for column names `CHR`, `POS`, `REF`, `ALT`, and `Z`. Alternatively, you can specify which column should be used for each of these fields by providing the corresponding optional inputs, e.g. `--CHR 6` tells `GhostKnockoffGWAS` to use column 6 as `CHR`. The `ALT` allele will be treated as the effect allele and `REF` be treated as non-effect allele. The POS (position) field of each variant must be from HG19 or HG38, which must be specified by the `--genome-build` argument. 

Here is a minimal example with 10 Z scores

```
CHR	POS	REF	ALT	Z
17	150509	T	TA	1.08773561923134
17	151035	T	C	0.703898767202681
17	151041	G	A	NaN
17	151872	T	C	-0.299877259561085
17	152087	C	T	-0.371627135786605
17	152104	G	A	-0.28387322965385
17	152248	G	A	0.901618600934489
17	152427	G	A	1.10987516000804
17	152771	A	G	0.708492545266136
```


!!! tip

    Missing Z scores can be specified as `NaN` or as an empty cell. If you do not want a SNP to be considered in the analysis, you can change the its Z-score to NaN. CHR/POS/REF/ALT fields cannot have missing values.

## Requirements on the input Z-scores

The input Z-scores (or p-values) must be valid Z-scores, i.e. each has $N(0, 1)$ distribution under the null. [This paper](https://arxiv.org/abs/2310.04030) shows that input Z-scores can be calculated from various marginal association tests, such as 
+ generalized linear mixed effect model to account for sample relatedness
+ saddle point approximation for extreme case-control imbalance
+ meta-analysis that aggregates multiple studies.

## Obtaining Z scores from p-values, odds ratios, or equivalent

If you have p-values, effect sizes, odds ratios...etc but not Z scores, you can convert them into Z score, for example by following the *Notes on computing Z-scores* of [this blog post](https://huwenboshi.github.io/data%20management/2017/11/23/tips-for-formatting-gwas-summary-stats.html). 
