"""
    @ode_def [name] begin
        dstate = expression
        # additional derivative assignments
    end parameter_names...

Define a parameterized ordinary differential equation and return an
`SciMLBase.AbstractParameterizedFunction` suitable for an `ODEProblem`. The macro
generates in-place and out-of-place right-hand sides, a time gradient, and a
Jacobian through ModelingToolkit.

# Arguments

- `name`: optional name for the generated function type. Omit it to return an
  anonymous parameterized function.
- derivative assignments: a `begin` block whose statements have the form
  `dstate = expression`. The identifier after `d` becomes a state variable. The
  independent variable is `t`.
- `parameter_names...`: symbols used as parameters in the derivative expressions.
  Their numerical values are supplied by the problem's `p` argument.

# Keywords

This macro does not accept keyword arguments.

# Returns

An `AbstractParameterizedFunction` with generated `f`, `tgrad`, and `jac`
methods. Pass it as the `f` argument of an `ODEProblem`.

# Examples

```julia
lotka_volterra = @ode_def begin
    dx = a * x - b * x * y
    dy = -c * y + d * x * y
end a b c d
```
"""
macro ode_def(name, ex, params...)
    opts = Dict{Symbol, Bool}(
        :build_tgrad => true,
        :build_jac => true,
        :build_expjac => false,
        :build_invjac => false,
        :build_invW => false,
        :build_hes => false,
        :build_invhes => false,
        :build_dpfuncs => true
    )
    return name isa Expr ? ode_def_opts(gensym(), opts, __module__, name, ex, params...) :
        ode_def_opts(name, opts, __module__, ex, params...)
end

"""
    @ode_def_bare [name] begin
        dstate = expression
        # additional derivative assignments
    end parameter_names...

Define a parameterized ODE without generating symbolic derivative functions.

# Arguments

- `name`: optional name for the generated function type.
- derivative assignments: the same assignment-based ODE DSL accepted by
  `@ode_def`.
- `parameter_names...`: symbols used as parameters in the expressions.

# Keywords

This macro does not accept keyword arguments.

# Returns

An `AbstractParameterizedFunction` with a generated right-hand side and `nothing`
for generated derivative-function fields. Use it when symbolic Jacobian and time
gradient generation is unnecessary.

# Examples

```julia
linear = @ode_def_bare begin
    dx = a * x
end a
```
"""
macro ode_def_bare(name, ex, params...)
    opts = Dict{Symbol, Bool}(
        :build_tgrad => false,
        :build_jac => false,
        :build_expjac => false,
        :build_invjac => false,
        :build_invW => false,
        :build_hes => false,
        :build_invhes => false,
        :build_dpfuncs => false
    )
    return name isa Expr ? ode_def_opts(gensym(), opts, __module__, name, ex, params...) :
        ode_def_opts(name, opts, __module__, ex, params...)
end

"""
    @ode_def_all [name] begin
        dstate = expression
        # additional derivative assignments
    end parameter_names...

Define a parameterized ODE and request every derivative function supported by this
DSL, including factorized `W` functions for small systems.

# Arguments

- `name`: optional name for the generated function type.
- derivative assignments: the same assignment-based ODE DSL accepted by
  `@ode_def`.
- `parameter_names...`: symbols used as parameters in the expressions.

# Keywords

This macro does not accept keyword arguments.

# Returns

An `AbstractParameterizedFunction` with generated right-hand side, time-gradient,
Jacobian, and supported factorized-`W` functions.

# Examples

```julia
linear = @ode_def_all begin
    dx = a * x
end a
```
"""
macro ode_def_all(name, ex, params...)
    opts = Dict{Symbol, Bool}(
        :build_tgrad => true,
        :build_jac => true,
        :build_expjac => false,
        :build_invjac => false,
        :build_invW => true,
        :build_hes => false,
        :build_invhes => false,
        :build_dpfuncs => true
    )
    return name isa Expr ? ode_def_opts(gensym(), opts, __module__, name, ex, params...) :
        ode_def_opts(name, opts, __module__, ex, params...)
end
