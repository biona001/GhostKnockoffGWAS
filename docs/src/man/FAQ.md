
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

## When will non-European LD files be available?

We will release more pre-processed LD files for download, once we tested and verified the methodology against suitable datasets. Most likely, the first non-EUR release will be on African and East Asian populations. 

## Admixed samples?

If your study subjects are somewhat admixed, one can try `GhostKnockoffGWAS` with the most suitable LD files, then check how much deviation there are by examining the `LD_shrinkage` parameter in the output of `GhostKnockoffGWAS`, see [Is-the-result-is-trustworthy?](https://biona001.github.io/GhostKnockoffGWAS/dev/man/FAQ/#Is-the-result-is-trustworthy?).

If your study subjects are extremely admixed, then it is unlikely that GhostKnockoffGWAS will return good results. The main difficulty in enabling analysis for admixed cohorts lies in pre-computing good LD files for admixed subjects. Computing required quantities on the fly is too computationally intensive. 

## How do I specify my own groups?

We will add this feature in the near future
