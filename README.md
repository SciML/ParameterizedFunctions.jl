# ParameterizedFunctions.jl

## Basic Usage

`ParameterizedFunction` is a type which can be used in various JuliaDiffEq solvers where
the parameters must be accessible by the solver function. These use call overloading
generate a type which acts like a function `f(t,u,du)` but has access to the model
parameters. An example of explicitly defining a parameterized function is as follows:

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

## Helper Macro

A helper macro is provided to make it easier to define a `ParameterizedFunction`.
For example, to define the previous `LotkaVoltera`, you can use the following command:

```julia
f = @ode_def LotkaVoltera begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a=>1.5 b=>1 c=3 d=1
```

This will silently create the `LotkaVoltera` type and thus `g=LotkaVoltera(1.0,2.0)`
will work as before. Note that the syntax for parameters here is that `=>` will
put these inside the parameter type, while `=` will inline the number (i.e. replace
each instance of `c` with `3`). Inlining slightly decreases the function cost and
so is preferred in any case where you know that the parameter will always be constant.
