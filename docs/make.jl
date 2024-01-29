using GhostKnockoffGWAS
using Documenter

makedocs(
    doctest = false,
    checkdocs = :none,
    sitename = "GhostKnockoffGWAS",
    format = Documenter.HTML(size_threshold = nothing),
    modules = [GhostKnockoffGWAS],
    authors = "Benjamin Chu",
    clean = true,
    pages = [
        "Home" => "index.md",
        "Getting started" => "man/getting_started.md",
        "Downloads" => "man/download.md",
        "Acceptable Z-scores file" => "man/zfile.md",
        "Detailed Example" => "man/examples.md",
        # "Video tutorials" => "man/video.md",
        "FAQ" => "man/FAQ.md",
        "Usage within Julia" => "man/julia.md",
        "Developer documentation" => "man/developer.md",
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo   = "github.com/biona001/GhostKnockoffGWAS.git",
    target = "build",
    deps   = nothing,
    make   = nothing,
)
