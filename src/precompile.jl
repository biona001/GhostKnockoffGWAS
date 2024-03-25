# fast GhostKnockoffGwAS --help
using GhostKnockoffGWAS
GhostKnockoffGWAS.parse_commandline(false)

# trying to resolve https://github.com/biona001/GhostKnockoffGWAS/issues/7
using HDF5
mktempdir() do dir
    file = joinpath(dir, "LD_test.h5")
    h5file = h5open(file, "w")
    h5file["S"] = randn(2, 2)
    h5file["groups"] = rand(Int, 2)
    close(h5file)
    result = h5open(file, "r")
    S = read(result, "S")
    groups = read(result, "groups")
    close(result)
end
