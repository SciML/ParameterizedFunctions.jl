# ParameterizedFunctions.jl

[![Travis](https://travis-ci.org/JuliaDiffEq/ParameterizedFunctions.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/ParameterizedFunctions.jl) [![AppVeyor](https://ci.appveyor.com/api/projects/status/k6b7d86ddbas1ajk?svg=true)](https://ci.appveyor.com/project/ChrisRackauckas/parameterizedfunctions-jl)

`ParameterizedFunction` is a type which can be used in various JuliaDiffEq solvers where
the parameters must be accessible by the solver function. These use call overloading
generate a type which acts like a function `f(t,u,du)` but has access to many more
features. For example, a `ParameterizedFunction` can contain a function for the Jacobian
or Inverse Jacboian. If such functions exist, the solvers can use them to increase
the speed of computations. If they don't exist, the solvers will ignore them. Since
`ParameterizedFunction` is a subtype of `Function`, these can be used anywhere that
a function can be used, just with the extra functionality ignored.

## Basic Usage via Macros

### ODEs

A helper macro is provided to make it easier to define a `ParameterizedFunction`.
For example, to define the previous `LotkaVolterra`, you can use the following command:

```julia
f = @ode_def LotkaVolterra begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a=>1.5 b=>1 c=3 d=1
```

Note that the syntax  for parameters here is that `=>` will
put these inside the parameter type, while `=` will inline the number (i.e. replace
each instance of `c` with `3`). Inlining slightly decreases the function cost and
so is preferred in any case where you know that the parameter will always be constant.
This will silently create the `LotkaVolterra` type and thus `g=LotkaVolterra(a=1.0,b=2.0)`
will create a different function where `a=1.0` and `b=2.0`. However, at any time
the parameters of `f` can be changed by using `f.a =` or `f.b = `, or by using
symbols: `f[:a]=` etc.

The macro also defines the Jacobian `f'`. This is defined as an in-place Jacobian `f(t,u,J,:Jac)`.
This is calculated using SymEngine.jl automatically, so it's no effort on your part.
The symbolic inverse of the Jacobian is also computed, and an in-place function
for this is available as well as `f(t,u,iJ,:InvJac)`. If the Jacobians cannot be
computed, a warning is thrown and only the function itself is usable. The booleans
`f.jac_exists` and `f.invjac_exists` can be used to see whether the Jacobian
and the function for its inverse exist.

### Finite Element PDEs

Similar macros for finite element method definitions also exist. For the finite
element solvers, the definitions use `x[:,1]` instead of `x` and `x[:,2]` instead of `y`.
To more easily define systems of equations for finite element solvers, we can
use the `@fem_def` macro. The first argument is the function signature. This
is required in order to tell the solver linearity. Other than that, the macro
usage is similar to before. For example,

```julia
l = @fem_def (t,x,u) BirthDeath begin
  du = 1-x*α*u
  dv = 1-y*v
end α=0.5
```

defines a system of equations

```julia
l = (t,x,u)  -> [1-.5*x[:,1]*u[:,1]   1-x[:,2]*u[:,2]]
```

which is in the form for the FEM solver.

## The ParameterizedFunction Interface

The ParameterizedFunction Interface is as follows:

- ParameterizedFunction is a type which is a subtype of Function
- The type must hold the parameters.
- The type must also define 6 booleans to declare the existence of explicit Jacobians,
  Hessians, Inverse Jacobians, Inverse Hessians, explicit parameter functions,
  and parameter derivatives.
- The standard call `(p::TypeName)(t,u,du)` must be overloaded for the function
  calculation.

Solvers can interface with ParameterizedFunctions as follows:

```julia
f.a # or f[:a], accesses the parameter a
f.jac_exists # Checks for the existence of the explicit Jacobian function
f(t,u,du) # Call the function
f(t,u,2.0,du,:a) # Call the explicit parameter function with a=2.0
f(t,u,2.0,df,:a,:Deriv) # Call the explicit parameter derivative function with a=2.0
f(t,u,J,:Jac) # Call the explicit Jacobian function
f(t,u,iJ,:InvJac) # Call the explicit Inverse Jacobian function
f(t,u,H,:Hes) # Call the explicit Hessian function
f(t,u,iH,:InvHes) # Call the explicit Inverse Hessian function
```

It is requested that solvers should only use the explicit functions when they exist.

## Manually Defining `ParameterizedFunction`s

### Template

An example of explicitly defining a parameterized function is as follows. This serves
as a general template for doing so:

```julia
type  LotkaVolterra <: ParameterizedFunction
         a::Float64
         b::Float64
         jac_exists::Bool
         invjac_exists::Bool
         hes_exists::Bool
         invhes_exists::Bool
         pfuncs_exists::Bool
         pderiv_exists::Bool
end
f = LotkaVolterra(0.0,0.0,true,true,true)
(p::LotkaVolterra)(t,u,du) = begin
         du[1] = p.a * u[1] - p.b * u[1]*u[2]
         du[2] = -3 * u[2] + u[1]*u[2]
end
(p::LotkaVolterra)(t,u,du,sym::Symbol) = p(t,u,du,Val{sym})
(p::LotkaVolterra)(t,u,param,du,sym::Symbol) = p(t,u,param,du,Val{sym})
(p::LotkaVolterra)(t,u,param,du,sym::Symbol,sym2::Symbol) = p(t,u,param,du,Val{sym},Val{sym2})
# Add Functions
```

### Explanation

Let's go step by step to see what this template does. The first part defines a
type:

```julia
type  LotkaVolterra <: ParameterizedFunction
         a::Float64
         b::Float64
         jac_exists::Bool
         invjac_exists::Bool
         hes_exists::Bool
         invhes_exists::Bool
         pfuncs_exists::Bool
         pderiv_exists::Bool
end
```

The fields are the parameters for our function, and 6 booleans used by the solvers:

- Does the explicit Jacobian function exist?
- Does the explicit Inverse Jacobian function exist?
- Does the explicit Hessian function exist?
- Does the explicit Inverse Hessian function exist?
- Do the explicit parameter functions exist?
- Does the Parameter Derivative function exist?

Then we built the type:

```julia
f = LotkaVolterra(0.0,0.0,true,true,true,true)
```

We put in values for the parameters and told it that we will be defining each of
those functions. First we define the main overload. This is required even if none
of the other functions are provided. The function for the main overload is the
differential equation, so for the Lotka-Volterra equation:

```julia
(p::LotkaVolterra)(t,u,du) = begin
         du[1] = p.a * u[1] - p.b * u[1]*u[2]
         du[2] = -3 * u[2] + u[1]*u[2]
end
```

Note how we represented the parameters in the equation. If you did this and set
the booleans to false, the result is `f` is a `ParameterizedFunction`,
but `f(t,u,du)` acts like the function:

```julia
function f(t,u,du)
         du[1] = 0.0 * u[1] - 0.0 * u[1]*u[2]
         du[2] = -3 * u[2] + u[1]*u[2]
end
```

At anytime the function parameters can be accessed by the fields (`f.a`, `f.b`) or
through `getindex` on the field symbol (`f[:a]`) by pre-defined dispatches.

### Extra Functions

Now we have to add the dispatches for the extra functions. The code from the template:

```julia
(p::LotkaVolterra)(t,u,du,sym::Symbol) = p(t,u,du,Val{sym})
(p::LotkaVolterra)(t,u,param,du,sym::Symbol) = p(t,u,param,du,Val{sym})
(p::LotkaVolterra)(t,u,param,du,sym::Symbol,sym2::Symbol) = p(t,u,param,du,Val{sym},Val{sym2})
```

sets up the value dispatches by symbols, and so we only need to define the function
Val-type dispatches.

#### Jacobian Function

The Jacobian overload is provided by overloading in the following manner:

```julia
function (p::LotkaVolterra)(t,u,du,J,::Type{Val{:Jac}})
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
function (p::LotkaVolterra)(t,u,du,J,::Type{Val{:Jac}})
  J[1,1] = (1 - (p.b * u[1] * u[2]) / ((p.a - p.b * u[2]) * (-3 + u[1] + (p.b * u[1] * u[2]) / (p.a - p.b * u[2])))) / (p.a - p.b * u[2])
  J[1,2] = (p.b * u[1]) / ((p.a - p.b * u[2]) * (-3 + u[1] + (p.b * u[1] * u[2]) / (p.a - p.b * u[2])))
  J[2,1] = -(u[2]) / ((p.a - p.b * u[2]) * (-3 + u[1] + (p.b * u[1] * u[2]) / (p.a - p.b * u[2])))
  J[2,2] = (-3 + u[1] + (p.b * u[1] * u[2]) / (p.a - p.b * u[2])) ^ -1
  nothing
end
```

#### Hessian and Inverse Hessian

These are the same as the Jacobians, except with value types `:Hes` and `:InvHes`.

#### Explicit Parameter Functions

For solvers which need to autodifferentiate parameters (local sensitivity analysis),
explicit parameter functions are required. For our example, we do the following:

```julia
function (p::LotkaVolterra)(t,u,a,du,::Type{Val{:a}})
  du[1] = a * u[1] - p.b * u[1] * u[2]
  du[2] = -3 * u[2] + 1 * u[1] * u[2]
  nothing
end
function (p::LotkaVolterra)(t,u,b,du,::Type{Val{:b}})
  du[1] = p.a * u[1] - b * u[1] * u[2]
  du[2] = -3 * u[2] + 1 * u[1] * u[2]
  nothing
end
```

#### Explicit Parameter Derivatives

For solvers which need parameters derivatives, specifying the functions can increase
performance. For our example, we allow the solvers to use the explicit derivatives
in the parameters `a` and `b` by:

```julia
function (p::LotkaVolterra)(t,u,a,du,::Type{Val{:a}},::Type{Val{:Deriv}})
  du[1] = 1 * u[1]
  du[2] = 1 * 0
  nothing
end
function (p::LotkaVolterra)(t,u,b,du,::Type{Val{:b}},::Type{Val{:Deriv}})
  du[1] = -(u[1]) * u[2]
  du[2] = 1 * 0
  nothing
end
```
