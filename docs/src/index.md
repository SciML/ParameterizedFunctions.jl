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

### Latex

Latexify directly works on the generated expression. Example:

```julia
using ParameterizedFunctions, Latexify

lotka_volterra = @ode_def begin
    dx = a * x - b * x * y
    dy = -c * y + d * x * y
end a b c d

latexify(lotka_volterra.sys)
```

Generates:

```
L"\begin{align}
\frac{\mathrm{d} x\left( t \right)}{\mathrm{d}t} =& a x\left( t \right) - b x\left( t \right) y\left( t \right) \\
\frac{\mathrm{d} y\left( t \right)}{\mathrm{d}t} =&  - c y\left( t \right) + d x\left( t \right) y\left( t \right)
\end{align}
"`
```

## Contributing

  - Please refer to the
    [SciML ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://github.com/SciML/ColPrac/blob/master/README.md)
    for guidance on PRs, issues, and other matters relating to contributing to SciML.

  - See the [SciML Style Guide](https://github.com/SciML/SciMLStyle) for common coding practices and other style decisions.
  - There are a few community forums:
    
      + The #diffeq-bridged and #sciml-bridged channels in the
        [Julia Slack](https://julialang.org/slack/)
      + The #diffeq-bridged and #sciml-bridged channels in the
        [Julia Zulip](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
      + On the [Julia Discourse forums](https://discourse.julialang.org)
      + See also [SciML Community page](https://sciml.ai/community/)

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

```@eval
using TOML
using Markdown
version = TOML.parse(read("../../Project.toml", String))["version"]
name = TOML.parse(read("../../Project.toml", String))["name"]
link_manifest = "https://github.com/SciML/" * name * ".jl/tree/gh-pages/v" * version *
                "/assets/Manifest.toml"
link_project = "https://github.com/SciML/" * name * ".jl/tree/gh-pages/v" * version *
               "/assets/Project.toml"
Markdown.parse("""You can also download the
[manifest]($link_manifest)
file and the
[project]($link_project)
file.
""")
```
