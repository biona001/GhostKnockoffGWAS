# fast GhostKnockoffGwAS --help
using GhostKnockoffGWAS
GhostKnockoffGWAS.parse_commandline(false)

# needed to resolve https://github.com/biona001/GhostKnockoffGWAS/issues/7
using HDF5, JLD2
mktempdir() do dir
    file = joinpath(dir, "LD_test.h5")
    JLD2.save(file, Dict("S" => randn(2, 2), "groups" => rand(Int, 2)))
    result = h5open(file, "r")
    S = read(result, "S")
    groups = read(result, "groups")
end
