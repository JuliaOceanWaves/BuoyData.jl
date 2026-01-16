using Documenter
import BuoyData

DocMeta.setdocmeta!(BuoyData, :DocTestSetup, :(import BuoyData); recursive = true)

makedocs(
    modules = [BuoyData],
    sitename = "BuoyData.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", "false") == "true"),
    pages = [
        "Home" => "index.md",
        "API" => "api.md"
    ]
)

deploydocs(
    repo = "github.com/JuliaOceanWaves/BuoyData.jl",
    devbranch = "main"
)
