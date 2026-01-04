using ExplicitImports
using ParameterizedFunctions
using Test
using ModelingToolkit
using ModelingToolkitBase
using Symbolics

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
