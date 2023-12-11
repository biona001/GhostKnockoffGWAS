using GhostKnockoffGWAS
using Documenter

makedocs(
    sitename = "GhostKnockoffGWAS",
    format = Documenter.HTML(size_threshold = nothing),
    modules = [GhostKnockoffGWAS],
    authors = "Benjamin Chu",
    clean = true,
    pages = [
        "Home" => "index.md",
        "Getting started" => "man/getting_started.md",
        "Examples" => "man/examples.md",
        "FAQ" => "man/FAQ.md",
        "Developer documentation" => "man/developer.md",
    ]
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
