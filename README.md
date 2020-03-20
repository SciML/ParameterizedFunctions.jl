# ParameterizedFunctions.jl

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Travis](https://travis-ci.org/JuliaDiffEq/ParameterizedFunctions.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/ParameterizedFunctions.jl) [![AppVeyor](https://ci.appveyor.com/api/projects/status/k6b7d86ddbas1ajk?svg=true)](https://ci.appveyor.com/project/ChrisRackauckas/parameterizedfunctions-jl)
[![codecov](https://codecov.io/gh/JuliaDiffEq/ParameterizedFunctions.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaDiffEq/ParameterizedFunctions.jl)


## Note that this library is deprecated in favor of [ModelingToolkit.jl](https://github.com/JuliaDiffEq/ModelingToolkit.jl). ParameterizedFunctions.jl will continue to work but is in maitanance mode.

ParameterizedFunctions.jl is a component of the JuliaDiffEq ecosystem which allows
for parameters to be explicitly present within functions. The interface which
ParameterizedFunctions describes allows for functionality which requires parameters,
such as parameter sensitivity analysis and parameter estimation, to be added to
the differential equation solvers of [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl).
While the interface itself is of importance to ecosystem developers,
ParameterizedFunctions.jl provides user-facing macros which make a
`ParameterizedFunction` easy to define, and automatically
include optimizations like explicit Jacobian functions and explicit inverse Jacobian
functions for the differential equation solvers to take advantage of. The result
is an easy to use API which allows for more functionality and more performance
optimizations than could traditionally be offered.

test

## The Basic Idea

`ParameterizedFunction` is a type which can be used in various JuliaDiffEq solvers where
the parameters must be accessible by the solver function. These use call overloading
generate a type which acts like a function `f(t,u,du)` but has access to many more
features. For example, a `ParameterizedFunction` can contain a function for the Jacobian
or Inverse Jacobian. If such functions exist, the solvers can use them to increase
the speed of computations. If they don't exist, the solvers will ignore them. Since
`ParameterizedFunction` is a subtype of `Function`, these can be used anywhere that
a function can be used, just with the extra functionality ignored.

## Basic Usage

### ODE Macros

A helper macro is provided to make it easier to define a `ParameterizedFunction`,
and it will symbolically compute a bunch of extra functions to make the differential
equation solvers run faster. For example, to define the previous `LotkaVolterra`,
you can use the following command:

```julia
f = @ode_def LotkaVolterra begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a b c d
```

or you can define it anonymously:

```julia
f = @ode_def begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a b c d
```

The macro also defines the Jacobian `f'`. This is defined as an in-place Jacobian `f(Val{:jac},t,u,J)`.
This is calculated using SymEngine.jl automatically, so it's no effort on your part.
The symbolic inverse of the Jacobian is also computed, and an in-place function
for this is available as well as `f(Val{:invjac},t,u,iJ)`. If the Jacobians cannot be
computed, a warning is thrown and only the function itself is usable. The functions
`jac_exists(f)` and `invjac_exists(f)` can be used to see whether the Jacobian
and the function for its inverse exist.

#### Extra Options

In most cases the `@ode_def` macro should be sufficient. This is because by default
the macro will simply calculate each function symbolically, and if it can't it
will simply throw a warning and move on. However, in extreme cases the symbolic
calculations may take a long time, in which case it is necessary to turn them
off. To do this, use the `ode_def_opts` function. The `@ode_def` macro simply defines the specifiable options:

```julia
opts = Dict{Symbol,Bool}(
      :build_tgrad => true,
      :build_jac => true,
      :build_expjac => false,
      :build_invjac => true,
      :build_invW => true,
      :build_invW_t => true,
      :build_hes => false,
      :build_invhes => false,
      :build_dpfuncs => true)
```

and calls the function `ode_def_opts(name::Symbol,opts,ex::Expr,params)`. Note that
params is an iterator holding expressions for the parameters.

In addition, one can also use their own function inside of the macro. For example:

```julia
f(x,y,d) = erf(x*y/d)
NJ = @ode_def FuncTest begin
  dx = a*x - b*x*y
  dy = -c*y + f(x,y,d)
end a b c d
```

will do fine. The symbolic derivatives will not work unless you define a derivative
for `f`.

#### Extra Macros

Instead of using `ode_def_opts` directly, one can use one of the following macros
to be more specific about what to not calculate. In increasing order of calculations:

```julia
@ode_def_bare
@ode_def
@ode_def_all
```

### Extra Functions

#### Jacobian Function

The Jacobian overload is provided by overloading in the following manner:

```julia
function (p::LotkaVolterra)(::Type{Val{:jac}},t,u,J)
  J[1,1] = p.a - p.b * u[2]
  J[1,2] = -(p.b) * u[1]
  J[2,1] = 1 * u[2]
  J[2,2] = -3 + u[1]
  nothing
end
```

#### Inverse Jacobian

The Inverse Jacobian overload is provided by overloading in the following manner:

```julia
function (p::LotkaVolterra)(::Type{Val{:invjac}},t,u,J)
  J[1,1] = (1 - (p.b * u[1] * u[2]) / ((p.a - p.b * u[2]) * (-3 + u[1] + (p.b * u[1] * u[2]) / (p.a - p.b * u[2])))) / (p.a - p.b * u[2])
  J[1,2] = (p.b * u[1]) / ((p.a - p.b * u[2]) * (-3 + u[1] + (p.b * u[1] * u[2]) / (p.a - p.b * u[2])))
  J[2,1] = -(u[2]) / ((p.a - p.b * u[2]) * (-3 + u[1] + (p.b * u[1] * u[2]) / (p.a - p.b * u[2])))
  J[2,2] = (-3 + u[1] + (p.b * u[1] * u[2]) / (p.a - p.b * u[2])) ^ -1
  nothing
end
```

#### Parameter Jacobian

For solvers which need parameters derivatives, specifying the functions can increase
performance. For our example, we allow the solvers to use the explicit derivatives
in the parameters by:

```julia
function (p::LotkaVolterra)(::Type{Val{:paramjac}},J,u,p,t)
    J[1, 1] = u[1] * 1
    J[1, 2] = -(u[1]) * u[2]
    J[1, 3] = 0 * 1
    J[1, 4] = 0 * 1
    J[2, 1] = 0 * 1
    J[2, 2] = 0 * 1
    J[2, 3] = -(u[2])
    J[2, 4] = u[1] * u[2]
    nothing
end
```
