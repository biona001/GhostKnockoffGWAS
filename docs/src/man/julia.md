
# Usage within Julia

`GhostKnockoffGWAS` is a regular Julia package, which can be used directly within Julia for greater flexibility. To install it, execute the following in Julia:
```julia
using Pkg
Pkg.add(url="https://github.com/biona001/ghostbasil_jll.jl")
Pkg.add(url="https://github.com/biona001/Ghostbasil.jl")
Pkg.add(url="https://github.com/biona001/GhostKnockoffGWAS")
```

!!! warning

    `GhostKnockoffGWAS` uses upstream `HDF5.jl` directly. A forked `HDF5.jl` install is no longer required.

## Usage example

The following example performs summary-statistics GWAS under the GhostKnockoff framework.


```julia
using GhostKnockoffGWAS

# file paths and directories
LD_files = "/home/groups/sabatti/.julia/dev/GhostKnockoffGWAS/data/EUR"
zfile = "/home/groups/sabatti/.julia/dev/GhostKnockoffGWAS/data/AD_Zscores_Meta_modified.txt"
outdir = "/home/groups/sabatti/.julia/dev/GhostKnockoffGWAS/data"

# specify sample size and human genome build
N = 506200
hg_build = 38

# read Z-scores using built-in function read_zscores
z, chr, pos, effect_allele, non_effect_allele = GhostKnockoffGWAS.read_zscores(zfile)

# run analysis
@time ghostknockoffgwas(LD_files, z, chr, pos, effect_allele, 
    non_effect_allele, N, hg_build, outdir, outname="test_alzheimers_meta")
```

## Function API

```@docs
ghostknockoffgwas
read_zscores
solve_blocks
```
