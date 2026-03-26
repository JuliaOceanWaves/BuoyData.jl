using Documenter
using DocumenterCitations
using BuoyData

DocMeta.setdocmeta!(BuoyData, :DocTestSetup, :(import BuoyData); recursive = true)
bib = CitationBibliography(joinpath(@__DIR__, "references.bib"))

makedocs(
    modules = [BuoyData],
    sitename = "BuoyData.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", "false") == "true",
        assets = String["assets/citations.css"]),
    pages = [
        "Home" => "index.md",
        "API" => "api.md"
    ],
    plugins = [bib]
)

deploydocs(
    repo = "github.com/JuliaOceanWaves/BuoyData.jl.git",
    devbranch = "main"
)
