
# Developer documentation

This is for advanced users who wish to build customized knockoff analysis pipelines. Currently, customization is not easy, but it is possible in principle at 2 levels: 

1. Constructing knockoff statistics for custom LD panels, including
    + Specifying which LD panel to use
    + Defining quasi-independent regions and groups
    + Solving the knockoff (convex) optimization problem
    + Saving the result in a easy-to-read format, which will be read in step 2
2. Ghost Knockoff sampling and high dimensional Lasso regression
    + Read pre-computed knockoff statistics from step 1
    + Sample Ghost Knockoffs
    + Fit a pseudo-lasso problem
    + Applying the knockoff filter

A full example is provided in 3 separate jupyter notebooks: [part 0](https://github.com/biona001/ghostknockoff-gwas-reproducibility/blob/main/chu_et_al/ghostknockoff-part0.ipynb), [part 1](https://github.com/biona001/ghostknockoff-gwas-reproducibility/blob/main/chu_et_al/ghostknockoff-part1.ipynb), and [part 2](https://github.com/biona001/ghostknockoff-gwas-reproducibility/blob/main/chu_et_al/ghostknockoff-part2.ipynb). If you need assistance on any of these steps, feel free to reach out to us. 

## 1. Constructing knockoff statistics for custom LD panels

This section gives a high-level overview up to part 1 of the GhostKnockoff pipeline (computation of knockoff statistics). 

+ Processing of LD panels (including downloading and importing the data matrices) is carried out by [EasyLD.jl](https://github.com/biona001/EasyLD.jl). This package should make it easy to import a region of the LD matrix into memory in Julia. 
+ To partition the extremely large LD matrix into manageable pieces, we directly adopted the output of [ldetect](https://bitbucket.org/nygcresearch/ldetect-data/src/master/) for which `AFR` (african), `ASN` (east Asians), and `EUR` (european) results are already available (position coordinates are given in HG19). For the EUR panel, the autosomes are partitioned into 1703 "quasi-independent" regions, see Figure S2 of [this paper](https://arxiv.org/abs/2310.15069) for summaries. 
+ Knockoff optimization problem was carried out by [Knockoffs.jl](https://github.com/biona001/Knockoffs.jl). In particular, we defined groups via average-linkage hierarchical clustering, chose group-key variants within each group via Algorithm A1 in the paper with threshold value $c=0.5$, and employed the maximum-entropy group-knockoff solver.

For details, please see section 5.1 and 5.2 of [this paper](https://arxiv.org/pdf/2310.15069.pdf). Note that the precomputed knockoff statistics includes everything up to this point. 

## 2. Ghost Knockoff sampling and high dimensional Lasso regression

Over 1703 quasi-independent blocks, we have assembled
```math
\begin{aligned}
    \Sigma =
    \begin{bmatrix}
        \Sigma_1 & & \\
        & \ddots & \\
        & & \Sigma_{1703}
    \end{bmatrix}, \quad
    S = 
    \begin{bmatrix}
        S_1 & & \\
        & \ddots & \\
        & & S_{1703}
    \end{bmatrix}, \quad
    S_i = 
    \begin{bmatrix}
        S_{i,1} & & \\
        & \ddots & \\
        & & S_{i,G_i}
    \end{bmatrix}
\end{aligned}
```
where $\Sigma_i$ are LD matrices obtained from the Pan-UKBB panel and $S_i$ is the group-block-diagonal matrices obtained by solving the knockoff optimization problem. Given a Z-score vector $z$, we can compute $r = \frac{1}{\sqrt{n}} z$, and `ghostbasil` will solve the following optimization problem with $\lambda \ge 0, p_i \ge 0$, and $0 \le \alpha \le 1$.
```math
\begin{aligned}
    \min \frac{1}{2}\beta^t A \beta - \beta^tr + \lambda\sum_ip_i\left(\alpha|\beta_i| + \frac{1-\alpha}{2}\beta_i^2\right)
\end{aligned}
```
In `GhostKnockoffGWAS`, we set $\alpha = 1$ (i.e. a Lasso problem) and $p_i = 1$ for all $i$. $A = \frac{1}{n}[X,\tilde{X}]'[X,\tilde{X}]$ and $\beta$ contains the effect size for both original variables and their knockoffs. 

To solve this problem, we leverage the fact that Lasso's objective is seprable over the blocks: as long as we can find a lambda sequence to be used for all blocks, we can fit each block separately. Since the max lambda is only related to the marginal correlation between each feature and $y$, and the knockoffs are exchangeable to the original features, we can use the original genome-wide Z-scores to compute the lambda sequence. 

Thus, for each block $i \in \{1,...,1703\}$, we will call `ghostbasil(Bi, r)` where
```math
\begin{aligned}
    B_i &= \text{BlockGroupGhostMatrix}(C_i, S_i, m+1)\\
    C_i &= \Sigma_i - S_i
\end{aligned}  
```
Note that, since we use representative variant approach, $S_i$ is generally a dense matrix. To input a dense matrix, we use Jame's function `BlockGroupGhostMatrix` with a single block. 
```math
\begin{aligned}
    B_i = \text{BlockGroupGhostMatrix}(C_i, S_i, m+1) = 
    \begin{bmatrix}
        C_i+S_i & C_i & ... & C_i\\
        C_i & C_i+S_i & ... & \\
        \vdots & & \ddots & \vdots\\
        C_i & C_i & & C_i + S_i
    \end{bmatrix}
\end{aligned}
```
with the understanding that $B_i$ is the covariance matrix for $(Z, \tilde{Z}_1,...,\tilde{Z}_m)$
```math
\begin{aligned}
    B_i = 
    \begin{bmatrix}
        \Sigma_i & \Sigma_i-S_i & ... & \Sigma_i-S_i\\
        \Sigma_i-S_i & \Sigma_i & ... & \\
        \vdots & & \ddots & \vdots\\
        \Sigma_i-S_i & \Sigma_i-S_i & & \Sigma_i
    \end{bmatrix} = 
    \begin{bmatrix}
        C_i+S_i & C_i & ... & C_i\\
        C_i & C_i+S_i & ... & \\
        \vdots & & \ddots & \vdots\\
        C_i & C_i & & C_i + S_i
    \end{bmatrix}
\end{aligned}
```
where $C_i = \Sigma_i - S_i$. In Julia, this functionality is supported via the [Ghostbasil.jl](https://github.com/biona001/ghostbasil.jl) package. 

## Compiling the binaries

1. Make sure `gcc` is available. We recommend version 7.1, but [avoid using GCC 11+](https://github.com/JuliaPackaging/Yggdrasil/blob/0d38df8bc8ad10cff5fba1c19a5932a84286fcd2/CONTRIBUTING.md#compatibility-tips).
2. Make sure `GhostKnockoffGWAS` is installed within Julia. 
3. `dev` the package via
```julia
]dev GhostKnockoffGWAS
```
4. compile using [PackageCompiler.jl](https://github.com/JuliaLang/PackageCompiler.jl)
```julia
using PackageCompiler, GhostKnockoffGWAS
src = normpath(pathof(GhostKnockoffGWAS), "../..")
des = normpath(pathof(GhostKnockoffGWAS), "../../app_linux_x86")
precompile_script = normpath(pathof(GhostKnockoffGWAS), "../precompile.jl")
@time create_app(src, des, 
    include_lazy_artifacts=true, 
    force=true, 
    precompile_execution_file=precompile_script
)
```
The last step takes 1-2 hours. 