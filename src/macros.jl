"""
```
@ode_def name begin
    differential equation
end parameters :: ODEFunction
```

## Definition of the Domain-Specific Language (DSL)

A helper macro is provided to make it easier to define a `ParameterizedFunction`,
and it will symbolically compute a bunch of extra functions to make the differential
equation solvers run faster. For example, to define the previous `LotkaVolterra`,
you can use the following command:

```julia
f = @ode_def LotkaVolterra begin
    dx = a * x - b * x * y
    dy = -c * y + d * x * y
end a b c d
```

or you can define it anonymously:

```julia
f = @ode_def begin
    dx = a * x - b * x * y
    dy = -c * y + d * x * y
end a b c d
```

`@ode_def` uses ModelingToolkit.jl internally and returns an `ODEFunction` with the
extra definitions (Jacobian, parameter Jacobian, etc.) defined through the MTK
symbolic tools.
"""
macro ode_def(name, ex, params...)
    opts = Dict{Symbol, Bool}(:build_tgrad => true,
        :build_jac => true,
        :build_expjac => false,
        :build_invjac => false,
        :build_invW => false,
        :build_hes => false,
        :build_invhes => false,
        :build_dpfuncs => true)
    name isa Expr ? ode_def_opts(gensym(), opts, __module__, name, ex, params...) :
    ode_def_opts(name, opts, __module__, ex, params...)
end

"""
```
@ode_def_bare name begin
    differential equation
end parameters :: ODEFunction
```

Like `@ode_def` but the `opts` options are set so that no symbolic functions are generated.
See the `@ode_def` docstring for more details.
"""
macro ode_def_bare(name, ex, params...)
    opts = Dict{Symbol, Bool}(:build_tgrad => false,
        :build_jac => false,
        :build_expjac => false,
        :build_invjac => false,
        :build_invW => false,
        :build_hes => false,
        :build_invhes => false,
        :build_dpfuncs => false)
    name isa Expr ? ode_def_opts(gensym(), opts, __module__, name, ex, params...) :
    ode_def_opts(name, opts, __module__, ex, params...)
end

"""
```
@ode_def_all name begin
    differential equation
end parameters :: ODEFunction
```

Like `@ode_def` but the `opts` options are set so that all possible symbolic functions are generated.
See the `@ode_def` docstring for more details.
"""
macro ode_def_all(name, ex, params...)
    opts = Dict{Symbol, Bool}(:build_tgrad => true,
        :build_jac => true,
        :build_expjac => false,
        :build_invjac => false,
        :build_invW => true,
        :build_hes => false,
        :build_invhes => false,
        :build_dpfuncs => true)
    name isa Expr ? ode_def_opts(gensym(), opts, __module__, name, ex, params...) :
    ode_def_opts(name, opts, __module__, ex, params...)
end
