using SciMLTesting, ParameterizedFunctions, Test

# The ModelingToolkit interface ParameterizedFunctions deliberately reexports so that
# `using ParameterizedFunctions` on its own is enough to write the documented workflow
# -- define a model with `@ode_def`, hand it to `ODEProblem`, `solve` it, and inspect
# the generated `System` in the `sys` field. Owned and documented upstream; kept in
# sync with the reexport `export` block in src/ParameterizedFunctions.jl and the
# "Reexported ModelingToolkit interface" section of docs/src/api.md.
const REEXPORTS = (
    # Using the result of `@ode_def`.
    :ODEFunction, :ODEProblem, :solve,
    # The `System` in the generated function's `sys` field, and its accessors.
    :System, :equations, :independent_variable, :mtkcompile, :observed, :parameters,
    :unknowns,
    # Symbolic variables, for working with that system.
    Symbol("@parameters"), Symbol("@variables"), :Differential, :Num,
)

run_qa(ParameterizedFunctions; reexports_allow = REEXPORTS)

@testset "Reexport surface" begin
    # Every approved reexport must actually be reachable from `using
    # ParameterizedFunctions`, so the allow-list cannot drift into approving names the
    # package no longer provides. `isdefined(@__MODULE__, ...)` tests the property
    # directly: this file's `using ParameterizedFunctions` is what has to bring the
    # name into scope.
    @testset "$name" for name in REEXPORTS
        @test name in names(ParameterizedFunctions)
        @test isdefined(@__MODULE__, name)
    end
end
