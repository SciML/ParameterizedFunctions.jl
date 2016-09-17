# ParameterizedFunctions.jl

`ParameterizedFunction` is a type which can be used in various JuliaDiffEq solvers where
the parameters must be accessible by the solver function. These use call overloading
generate a type which acts like a function `f(t,u,du)` but has access to the model
parameters.

## Basic Usage

### ODEs

A helper macro is provided to make it easier to define a `ParameterizedFunction`.
For example, to define the previous `LotkaVoltera`, you can use the following command:

```julia
f = @ode_def LotkaVoltera begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a=>1.5 b=>1 c=3 d=1
```

Note that the syntax  for parameters here is that `=>` will
put these inside the parameter type, while `=` will inline the number (i.e. replace
each instance of `c` with `3`). Inlining slightly decreases the function cost and
so is preferred in any case where you know that the parameter will always be constant.
This will silently create the `LotkaVoltera` type and thus `g=LotkaVoltera(1.0,2.0)`
will create a different function where `a=1.0` and `b=2.0`. However, at any time
the parameters of `f` can be changed by using `f.a =` or `f.b = `, or by using
symbols: `f[:a]=` etc.

The macro also defines the Jacobian `f'`. This is defined as an in-place Jacobian `f'(t,u,J)`.
This is calculated using SymPy.jl automatically, so it's no effort on your part.

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

## Manually Defining `ParameterizedFunction`s

An example of explicitly defining a parameterized function is as follows:

```julia
type  LotkaVoltera <: ParameterizedFunction
         a::Float64
         b::Float64
end
f = LotkaVoltera(0.0,0.0)
(p::LotkaVoltera)(t,u,du) = begin
         du[1] = p.a * u[1] - p.b * u[1]*u[2]
         du[2] = -3 * u[2] + u[1]*u[2]
end
```

In this case, the result is `f` is a `ParameterizedFunction`, but `f(t,u,du)` acts
like the function:

```julia
function f(t,u,du)
         du[1] = 0.0 * u[1] - 0.0 * u[1]*u[2]
         du[2] = -3 * u[2] + u[1]*u[2]
end
```

At anytime the function parameters can be accessed by the fields (`f.a`, `f.b`) or
through `getindex` on the field symbol (`f[:a]`). To make an instance of the function
with different parameters, one just has to call the type constructor.
For example, one can define now define the function

```julia
function g(t,u,du)
         du[1] = 1.0 * u[1] - 2.0 * u[1]*u[2]
         du[2] = -3 * u[2] + u[1]*u[2]
end
```

by using the command `g = LotkaVoltera(1.0,2.0)`.

Note that the Jacobian overload is achieved by overloading
`Base.ctranspose`, and in the example corresponds to

```julia
function Base.ctranspose(p::LotkaVoltera) = (t,u,J) -> begin
  J[1,1] = p.a-p.b
  J[1,2] = -p.b
  J[2,1] = u[2]
  J[2,2] = -3 + u[1]
end
```
