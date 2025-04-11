# Common questions and Answers

Here is a collection of common questions & answers. If you have a question not listed here, do not hesitate to open a new issue on Github. 

## How do I obtain Z-scores from p-values, effect sizes, odds-ratios...etc?

See the **Notes on computing Z-scores** section of [this blog post](https://huwenboshi.github.io/data%20management/2017/11/23/tips-for-formatting-gwas-summary-stats.html)


## Is the result is trustworthy?

Knockoff's FDR guarantees requires that the correlation matrices used in the analysis approximates the LD structure for the original GWAS study. Their consistency is measured by the `mean_LD_shrinkage` parameter in the summary output. This value lies between 0 and 1. Values close to 0 indicates good performance. Larger values (e.g. >0.1) indicates deviation. Very larger values (e.g. >0.25) will cause the program to hault and users should download a different set of precomputed knockoff data instead. See equation 24 of [the SuSiE paper](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010299) for details. 

## Expected run time? 

On roughly 0.6 million Z-scores, our software completed a GhostKnockoff analysis in roughly 15 minutes on a single 2.3GHz CPU. If your analysis is taking much longer than this, please see Q&A on [software unexpectedly slow](https://biona001.github.io/GhostKnockoffGWAS/dev/man/FAQ/#Software-unexpectedly-slow?).

## Memory requirement?

Our software requires ~9.1GB of RAM on an Alzheimer's Diseases anslysis with ~0.6 million SNPs. 

## Software unexpectedly slow?

Because the knockoff pipeline requires reading the pre-computed knockoff statistics sequentially into memory, both the downloaded software and data should be stored at a high speed (e.g. Lustre) file system. On most HPC clusters, one should store the data in the `SCRATCH` directory to run `GhostKnockoffGWAS`. 

To check whether I/O is the bottleneck, one can check the CPU usage while `GhostKnockoffGWAS` is running. For example, one can examine CPU usage via the `top` or `htop` command. CPU usage should almost always be at 99% or above.  

For undiagnosed performance issues, please file a new issue. 

## Sex chromosome support?

We currently do not support X/Y/M chromosome analysis.

## When will LD files for population X be available?

We plan to release more pre-processed LD files for download, once we tested and verified the methodology against suitable datasets. For now, if no suitable LD files exist in the downloads page, users must build their own LD files using the `solveblock` executable, see [Customizing LD files](https://biona001.github.io/GhostKnockoffGWAS/dev/man/solveblocks).

## Admixed samples?

If your study subjects are somewhat admixed, one can try `GhostKnockoffGWAS` with the most suitable LD files, then check how much deviation there are by examining the `LD_shrinkage` parameter in the output of `GhostKnockoffGWAS`, see [Is-the-result-is-trustworthy?](https://biona001.github.io/GhostKnockoffGWAS/dev/man/FAQ/#Is-the-result-is-trustworthy?).

If your study subjects are extremely admixed, then we advise you to perform simulations to check that FDR is controlled. 

## I want to build my own LD files. How do I determine `start_bp` and `end_bp`?

See [this paragraph](https://biona001.github.io/GhostKnockoffGWAS/dev/man/solveblocks/#Determining-start_bp-and-end_bp) of the documentation.

## Can your pipeline integrate with existing LD panels from gnomAD and/or Pan-UKB?

Yes in principle, but one needs to convert those LD panels into the format that is accepted by `GhostKnockoffGWAS`. This process is more involved than using `solveblock` on individual level data. The issue is connected to the fact that [gnomAD](https://gnomad.broadinstitute.org/downloads#v2-linkage-disequilibrium) and [Pan-UKB](https://pan-dev.ukbb.broadinstitute.org/docs/hail-format/index.html) distribute LD panels in Hail's [block-matrix format](https://hail.is/docs/0.2/linalg/hail.linalg.BlockMatrix.html#blockmatrix), a secret format implemented in python that (as of 2025) can only be read by the Hail package. This means that we need "glue code" to interface it with knockoff optimization solvers, and it is difficult to compile these glue code into standalone executables that users can call easily. 

In [developer's documentation](https://biona001.github.io/GhostKnockoffGWAS/dev/man/developer/#1.-Processing-downloaded-LD-panels), we briefly discussed how pre-existing LD panels from gnomAD and Pan-UKB can be [downloaded and imported](https://github.com/biona001/ghostknockoff-gwas-reproducibility/blob/main/chu_et_al/ghostknockoff-part0.ipynb) by [EasyLD.jl](https://github.com/biona001/EasyLD.jl) within Julia, which can then be used for [knockoff optimization](https://github.com/biona001/ghostknockoff-gwas-reproducibility/blob/main/chu_et_al/ghostknockoff-part1.ipynb). While the relevant steps were achieved and documented in these 2 notebooks, it is not a beginner friendly pipeline. We suggest interested users to reach out to us before attemping to carryout the process. 

## Another question?

Feel free to email us directly, or [file an issue on Github](https://github.com/biona001/GhostKnockoffGWAS/issues)
