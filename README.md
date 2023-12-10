# GhostKnockoffGWAS

This is a package to perform Ghost Knockoff analysis for summary statistics data. 

:warning: **Package is in early phase of development**. 

:warning: **This package currently only works on Julia v1.8.x or 1.9.x** with Linux aarch64 or x86_64. If you need it to work on a different Julia version, please file an issue. 

## Usage and documentation

Coming soon

## Installation

Download and install [Julia](https://julialang.org/downloads/). Within Julia, copy and paste the following: 
```julia
using Pkg
Pkg.add(url="https://github.com/biona001/ghostbasil_jll.jl")
Pkg.add(url="https://github.com/biona001/Ghostbasil.jl")
Pkg.add(url="https://github.com/biona001/GhostKnockoffGWAS")
```

## Compiling the binaries


1. `ml gcc/7.1`
2. Make sure `GhostKnockoffGWAS` is installed
3. `dev` the package via
```julia
]dev GhostKnockoffGWAS
```
4. compile using [PackageCompiler.jl](https://github.com/JuliaLang/PackageCompiler.jl)
```julia
using PackageCompiler, GhostKnockoffGWAS
src = normpath(pathof(GhostKnockoffGWAS), "../..")
des = normpath(pathof(GhostKnockoffGWAS), "../../app_linux_x86")
@time create_app(src, des, include_lazy_artifacts=true) # takes 1-2h
```
