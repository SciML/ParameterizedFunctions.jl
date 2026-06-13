using ParameterizedFunctions, Aqua
using ExplicitImports
using Test
using ModelingToolkit
using ModelingToolkitBase
using Symbolics

@testset "Aqua" begin
    Aqua.find_persistent_tasks_deps(ParameterizedFunctions)
    Aqua.test_ambiguities(ParameterizedFunctions, recursive = false)
    Aqua.test_deps_compat(ParameterizedFunctions)
    Aqua.test_piracies(
        ParameterizedFunctions,
        treat_as_own = [ParameterizedFunctions.SciMLBase.AbstractParameterizedFunction]
    )
    Aqua.test_project_extras(ParameterizedFunctions)
    Aqua.test_stale_deps(ParameterizedFunctions)
    Aqua.test_unbound_args(ParameterizedFunctions)
    # MTK v11 has undefined exports (Variable, find_solvables!) that get re-exported
    # This is an upstream bug in ModelingToolkit, not ParameterizedFunctions
    Aqua.test_undefined_exports(ParameterizedFunctions; broken = true)
end

@testset "ExplicitImports" begin
    # Skip modules that are re-exported (names from @reexport using ModelingToolkit
    # are intentionally brought in implicitly for re-export purposes), and Base/Core
    # which are always available implicitly
    @test check_no_implicit_imports(
        ParameterizedFunctions;
        skip = (Base, Core, ModelingToolkit, ModelingToolkitBase, Symbolics)
    ) === nothing
    @test check_no_stale_explicit_imports(ParameterizedFunctions) === nothing
end
