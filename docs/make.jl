using Documenter
using NDBC

DocMeta.setdocmeta!(NDBC, :DocTestSetup, :(using NDBC); recursive = true)

makedocs(
    modules = [NDBC],
    sitename = "NDBC.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", "false") == "true"),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
    ],
)

deploydocs(
    repo = "github.com/JuliaOceanWaves/NDBC.jl",
    devbranch = "main",
)
