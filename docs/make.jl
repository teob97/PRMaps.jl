push!(LOAD_PATH,"src/")

using Documenter
using PRMaps

makedocs(
    sitename = "PRMaps",
    format = Documenter.HTML(),
    modules = [PRMaps],
    pages = [
        "Introduction" => "index.md",
        "Simple map making" => "simpleMap.md",
        "Polarized map" => "polarizedMap.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/teob97/PRMaps.jl.git",
)
