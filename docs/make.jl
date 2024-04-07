using Documenter, ParameterizedFunctions

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

include("pages.jl")

makedocs(sitename = "ParameterizedFunctions.jl",
    authors = "Chris Rackauckas",
    modules = [ParameterizedFunctions],
    clean = true, doctest = false, linkcheck = true,
    format = Documenter.HTML(assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/ParameterizedFunctions/stable/"),
    pages = pages)

deploydocs(repo = "github.com/SciML/ParameterizedFunctions.jl.git";
    push_preview = true)
