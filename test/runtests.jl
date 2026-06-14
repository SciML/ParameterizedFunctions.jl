using SafeTestsets

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "QA"
    using Pkg
    Pkg.activate(joinpath(@__DIR__, "qa"))
    Pkg.develop(path = joinpath(@__DIR__, ".."))
    Pkg.instantiate()
    include(joinpath(@__DIR__, "qa", "qa.jl"))
end

if GROUP == "All" || GROUP == "Core"
    @safetestset "Core" begin
        include("./core.jl")
    end
end
