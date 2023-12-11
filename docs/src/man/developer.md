
# Developer documentation

This is for advanced users who wishes to build customized knockoff analysis pipelines. This is separated into 2 parts. If you would like assistance on any of these steps, feel free to reach out to us. 


## Constructing knockoff statistics for custom LD panels

+ Processing of LD panels (including downloading and importing the data matrices) is carried out by [EasyLD.jl](https://github.com/biona001/EasyLD.jl). This package should make it easy to import a region of the LD matrix into memory in Julia.
+ To partition the extremely large LD matrix into manageable pieces, we directly adopted the output of [ldetect](https://bitbucket.org/nygcresearch/ldetect-data/src/master/) for which `AFR` (african), `ASN` (east Asians), and `EUR` (european) results are already available (position coordinates are given in HG19). 
+ Knockoff optimization problem was carried out by [Knockoffs.jl](https://github.com/biona001/Knockoffs.jl).

Because pre-computed knockoff statistics are available for download, users do not have to manually install EasyLD.jl nor Knockoffs.jl to carry out this step.


## Using GhostKnockoffGWAS as a Julia package

`GhostKnockoffGWAS` is a regular Julia package, which can be used directly within Julia for greater flexibility. To install it, execute the following in Julia
```julia
using Pkg
Pkg.add(url="https://github.com/biona001/ghostbasil_jll.jl")
Pkg.add(url="https://github.com/biona001/Ghostbasil.jl")
Pkg.add(url="https://github.com/biona001/GhostKnockoffGWAS")
```

### Julia API

```@autodocs
Modules = [GhostKnockoffGWAS]
Order   = [:function, :type]
```

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
@time create_app(src, des, include_lazy_artifacts=true, force=true)
```
The last step takes 1-2 hours. 
