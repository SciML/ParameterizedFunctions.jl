using Documenter, ParameterizedFunctions

include("pages.jl")

makedocs(sitename = "ParameterizedFunctions.jl",
         authors = "Chris Rackauckas",
         modules = [ParameterizedFunctions],
         clean = true, doctest = false,
         format = Documenter.HTML(analytics = "UA-90474609-3",
                                  assets = ["assets/favicon.ico"],
                                  canonical = "https://docs.sciml.ai/ParameterizedFunctions/stable/"),
         pages = pages)

deploydocs(repo = "github.com/SciML/ParameterizedFunctions.jl.git";
           push_preview = true)
