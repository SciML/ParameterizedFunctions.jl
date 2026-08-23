# API

The public API of ParameterizedFunctions is the `@ode_def` family of macros,
documented on [The ode_def macro](@ref) page. `using ParameterizedFunctions` also
reexports a small part of the ModelingToolkit interface, described below.

## Reexported ModelingToolkit interface

`@ode_def` is a front end for [ModelingToolkit](https://docs.sciml.ai/ModelingToolkit/stable/):
it builds a `System` from the ODE DSL, keeps it in the generated function's `sys`
field, and returns a function whose documented use is to be handed to an
`ODEProblem` and solved. `using ParameterizedFunctions` brings the names that
workflow touches into scope, so they do not have to be imported separately:

```julia
using ParameterizedFunctions

lotka_volterra = @ode_def begin
    dx = a * x - b * x * y
    dy = -c * y + d * x * y
end a b c d

prob = ODEProblem(lotka_volterra, [1.0, 1.0], (0.0, 10.0), [1.5, 1.0, 3.0, 1.0])
equations(lotka_volterra.sys)
```

ParameterizedFunctions only reexports these names — they are owned and documented
upstream, at the links below.

  - Using the result of `@ode_def`, owned by
    [SciMLBase](https://docs.sciml.ai/SciMLBase/stable/) and
    [CommonSolve](https://github.com/SciML/CommonSolve.jl): `ODEProblem`,
    `ODEFunction`, `solve`
  - The generated system and its accessors, owned by
    [ModelingToolkitBase](https://docs.sciml.ai/ModelingToolkit/stable/) (the
    type of the `sys` field): `System`, `equations`, `unknowns`, `parameters`,
    `observed`, `independent_variable`, `mtkcompile`
  - Symbolic variables, for building expressions to use with that system, owned by
    [Symbolics](https://docs.sciml.ai/Symbolics/stable/) and ModelingToolkitBase:
    `@variables`, `@parameters`, `Differential`, `Num`

Anything else from ModelingToolkit, Symbolics or SciMLBase must be imported from
those packages directly. In particular, ParameterizedFunctions does not reexport
solver algorithms — `solve(prob, Tsit5())` still needs
[OrdinaryDiffEq](https://docs.sciml.ai/DiffEqDocs/stable/) (or
DifferentialEquations) — nor the acausal component-modelling surface of
ModelingToolkit (`@mtkmodel`, `@named`, connectors, and the rest), which the
`@ode_def` DSL does not use.
