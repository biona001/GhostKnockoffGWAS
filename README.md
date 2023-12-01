# GhostKnockoffGWAS

This is a package to perform Ghost Knockoff analysis for summary statistics data. 

:warning: **Package is in early phase of development**. 

:warning: **This package currently only works on Julia v1.9.0 or v1.9.3** with Linux aarch64 or x86_64. Note Julia version **must** be v1.9.0 or v1.9.3. 

## Usage and documentation

Coming soon

## Installation

Download and install [Julia](https://julialang.org/downloads/) v1.9.0 or v1.9.3. Within Julia, copy and paste the following: 
```julia
using Pkg
Pkg.add(url="https://github.com/biona001/ghostbasil_jll.jl")
Pkg.add(url="https://github.com/biona001/Ghostbasil.jl")
Pkg.add(url="https://github.com/biona001/GhostKnockoffGWAS")
```

## Compiling the binaries

1. Make sure `GhostKnockoffGWAS` is installed
2. `dev` the package via
```julia
]dev GhostKnockoffGWAS
```
3. compile using [PackageCompiler.jl](https://github.com/JuliaLang/PackageCompiler.jl)
```julia
using PackageCompiler
src = normpath(pathof(GhostKnockoffGWAS), "../..")
des = normpath(pathof(GhostKnockoffGWAS), "../../app_linux_x86")
@time create_app(src, des, include_lazy_artifacts=true) # takes 1-2h
```
