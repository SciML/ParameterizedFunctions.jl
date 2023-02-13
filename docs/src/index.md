# ParameterizedFunctions.jl: Simple High-Level ODE Definitions

ParameterizedFunctions.jl is a component of the SciML ecosystem which allows
for easily defining parameterized ODE models in a simple syntax.

## Installation

To install ParameterizedFunctions.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("ParameterizedFunctions")
```

## Example

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

## Contributing

  - Please refer to the
    [SciML ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://github.com/SciML/ColPrac/blob/master/README.md)
    for guidance on PRs, issues, and other matters relating to contributing to SciML.

  - There are a few community forums:
    
      + the #diffeq-bridged channel in the [Julia Slack](https://julialang.org/slack/)
      + [JuliaDiffEq](https://gitter.im/JuliaDiffEq/Lobby) on Gitter
      + on the [Julia Discourse forums](https://discourse.julialang.org)
      + see also [SciML Community page](https://sciml.ai/community/)

## Reproducibility

```@raw html
<details><summary>The documentation of this SciML package was built using these direct dependencies,</summary>
```

```@example
using Pkg # hide
Pkg.status() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>and using this machine and Julia version.</summary>
```

```@example
using InteractiveUtils # hide
versioninfo() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>
```

```@example
using Pkg # hide
Pkg.status(; mode = PKGMODE_MANIFEST) # hide
```

```@raw html
</details>
```

```@raw html
You can also download the 
<a href="
```

```@eval
using TOML
version = TOML.parse(read("../../Project.toml", String))["version"]
name = TOML.parse(read("../../Project.toml", String))["name"]
link = "https://github.com/SciML/" * name * ".jl/tree/gh-pages/v" * version *
       "/assets/Manifest.toml"
```

```@raw html
">manifest</a> file and the
<a href="
```

```@eval
using TOML
version = TOML.parse(read("../../Project.toml", String))["version"]
name = TOML.parse(read("../../Project.toml", String))["name"]
link = "https://github.com/SciML/" * name * ".jl/tree/gh-pages/v" * version *
       "/assets/Project.toml"
```

```@raw html
">project</a> file.
```
