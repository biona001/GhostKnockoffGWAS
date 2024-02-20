
# Usage within Julia

`GhostKnockoffGWAS` is a regular Julia package, which can be used directly within Julia for greater flexibility. To install it, execute the following in Julia
```julia
using Pkg
Pkg.add(url="https://github.com/biona001/ghostbasil_jll.jl")
Pkg.add(url="https://github.com/biona001/Ghostbasil.jl")
Pkg.add(url="https://github.com/biona001/GhostKnockoffGWAS")
```

!!! warning

    This package currently only works on Linux machines with Julia 1.8.x, 1.9.x, and 1.10.0. If you need it to work on a different Julia version, let us know by filing an issue on Github. 

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
```

## Compiling GhostKnockoffGWAS

1. `]add libcxxwrap_julia_jll` (note: as of Feb 2024, `libcxxwrap_julia_jll` must be v0.11.x)
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
