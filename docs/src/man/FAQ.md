
# Common questions and Answers

Here is a collection of common questions & answers. If you have a question not listed here, do not hesitate to open a new issue on Github. 

## How do I obtain Z-scores from p-values, effect sizes, odds-ratios...etc?

See the **Notes on computing Z-scores** section of [this blog post](https://huwenboshi.github.io/data%20management/2017/11/23/tips-for-formatting-gwas-summary-stats.html)

## Expected run time? 

On roughly 0.6 million Z-scores, our software completed a GhostKnockoff analysis in roughly 15 minutes on a single 2.3GHz CPU. If your analysis is taking much longer than this, please see Q&A on [software unexpectedly slow](https://biona001.github.io/GhostKnockoffGWAS/dev/man/FAQ/#Software-unexpectedly-slow?).

## Software unexpectedly slow?

Because the knockoff pipeline requires reading the pre-computed knockoff statistics sequentially into memory, both the downloaded software and data should be stored at a high speed (e.g. Lustre) file system. On most HPC clusters, one should store the data in the `SCRATCH` directory to run `GhostKnockoffGWAS`. 

To check whether I/O is the bottleneck, one can check the CPU usage while `GhostKnockoffGWAS` is running. For example, one can examine CPU usage via the `top` or `htop` command. CPU usage should almost always be at 99% or above.  

For undiagnosed performance issues, please file a new issue. 

## Memory requirement?

Our software requires ~9.1GB of RAM on an Alzheimer's Diseases anslysis with ~0.6 million SNPs. 

## Sex chromosome support?

We currently do not support X/Y/M chromosome analysis.

## Non-European knockoff statistics?

We will release more pre-computed knockoff statistics for download, once we tested and verified the methodology against suitable datasets. Most likely, the first non-EUR release will be on African and East Asian populations. 

## How do I specify my own groups?

We will add this feature in the near future
