# ParameterizedFunctions.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/ParameterizedFunctions/stable/)

[![codecov](https://codecov.io/gh/SciML/ParameterizedFunctions.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SciML/ParameterizedFunctions.jl)
[![Build Status](https://github.com/SciML/ParameterizedFunctions.jl/workflows/CI/badge.svg)](https://github.com/SciML/ParameterizedFunctions.jl/actions?query=workflow%3ACI)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

ParameterizedFunctions.jl is a component of the SciML ecosystem which allows
for easily defining parameterized ODE models in a simple syntax.

## Tutorials and Documentation

For information on using the package,
[see the stable documentation](https://docs.sciml.ai/ParameterizedFunctions/stable/). Use the
[in-development documentation](https://docs.sciml.ai/ParameterizedFunctions/dev/) for the version of
the documentation, which contains the unreleased features.

## Example

The following are valid ODE definitions.

```julia
using DifferentialEquations, ParameterizedFunctions

# Non-Stiff ODE

lotka_volterra = @ode_def begin
    d🐁 = α * 🐁 - β * 🐁 * 🐈
    d🐈 = -γ * 🐈 + δ * 🐁 * 🐈
end α β γ δ

p = [1.5, 1.0, 3.0, 1.0];
u0 = [1.0; 1.0];
prob = ODEProblem(lotka_volterra, u0, (0.0, 10.0), p)
sol = solve(prob, Tsit5(), reltol = 1e-6, abstol = 1e-6)

# Stiff ODE

rober = @ode_def begin
    dy₁ = -k₁ * y₁ + k₃ * y₂ * y₃
    dy₂ = k₁ * y₁ - k₂ * y₂^2 - k₃ * y₂ * y₃
    dy₃ = k₂ * y₂^2
end k₁ k₂ k₃

prob = ODEProblem(rober, [1.0, 0.0, 0.0], (0.0, 1e5), [0.04, 3e7, 1e4])
sol = solve(prob)
```
