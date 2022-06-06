# ParameterizedFunctions.jl

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Build Status](https://github.com/SciML/ParameterizedFunctions.jl/workflows/CI/badge.svg)](https://github.com/SciML/ParameterizedFunctions.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/SciML/ParameterizedFunctions.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SciML/ParameterizedFunctions.jl)
[![Coverage Status](https://coveralls.io/repos/github/SciML/ParameterizedFunctions.jl/badge.svg?branch=master)](https://coveralls.io/github/SciML/ParameterizedFunctions.jl?branch=master)

ParameterizedFunctions.jl is a component of the SciML ecosystem which allows
for easily defining parameterized ODE models in a simple syntax.

## Tutorials and Documentation

For information on using the package,
[see the stable documentation](https://parameterizedfunctions.sciml.ai/stable/). Use the
[in-development documentation](https://parameterizedfunctions.sciml.ai/dev/) for the version of
the documentation, which contains the unreleased features.

## Example

The following are valid ODE definitions.

```julia
f = @ode_def begin
  dy₁ = -k₁*y₁+k₃*y₂*y₃
  dy₂ =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
  dy₃ =  k₂*y₂^2
end k₁ k₂ k₃

f = @ode_def begin
  d🐁  = α*🐁  - β*🐁*🐈
  d🐈 = -γ*🐈 + δ*🐁*🐈
end α β γ δ
```