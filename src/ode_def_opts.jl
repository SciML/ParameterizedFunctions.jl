findreplace(ex::Symbol, dict) = get(dict, ex, ex)
function findreplace(ex::Expr, dict)
    Expr(ex.head, map(x -> findreplace(x, dict), ex.args)...)
end
findreplace(ex, dict) = ex

"""
```julia
ode_def_opts(name::Symbol, opts::Dict{Symbol, Bool}, curmod, ex::Expr, params...;
    depvar = :t)
```

The core internal. Users should only interact with this through the `@ode_def_*` macros.

Options are self-explanatory by name mapping to `ODEFunction`:

  - build_tgrad
  - build_jac
  - build_expjac
  - build_invjac
  - build_invW
  - build_invW_t
  - build_hes
  - build_invhes
  - build_dpfuncs

`depvar` sets the symbol for the dependent variable.

Example:

```julia
opts = Dict{Symbol, Bool}(:build_tgrad => true,
    :build_jac => true,
    :build_expjac => false,
    :build_invjac => true,
    :build_invW => true,
    :build_invW_t => true,
    :build_hes => false,
    :build_invhes => false,
    :build_dpfuncs => true)
```
"""
function ode_def_opts(name::Symbol, opts::Dict{Symbol, Bool}, curmod, ex::Expr, params...;
        depvar = :t)
    # depvar is the dependent variable. Defaults to t
    # M is the mass matrix in RosW, must be a constant!

    origex = copy(ex) # Save the original expression

    if !(eltype(params) <: Symbol)
        error("The syntax for ParameterizedFunctions has changed. Simply list the parameters at the end, i.e. `a b c d`, instead of `a=5.0 b=>3.0 ...`. Parameters are defined in the problem type. See the documentation for more information.")
    end
    params = Symbol[params...]

    ## Build independent variable dictionary
    indvar_dict, syms = build_indvar_dict(ex, depvar)
    ####

    t = Symbolics.unwrap((@variables t)[1])
    vars = Symbolics.unwrap.([(@variables $x(t))[1] for x in syms])
    params = Symbolics.unwrap.([(@parameters $x)[1] for x in Symbol[params...]])

    vars_dict = Dict(x => Symbol(v) for (x, v) in zip(syms, vars))

    # replace x with x(t) if it's a var
    ex = findreplace(ex, vars_dict)

    # Build the Expressions

    # Run find replace to make the function expression
    symex = copy(ex) # Different expression for symbolic computations
    #ode_findreplace(ex,symex,indvar_dict,params)
    funcs = build_component_funcs(symex)
    mtk_ops = modelingtoolkitize_expr.(funcs, ([t; vars; params],), (curmod,))

    D = ModelingToolkit.Differential(t)

    mtk_diffeqs = [D(vars[i]) ~ mtk_ops[i] for i in 1:length(vars)]

    sys = ODESystem(mtk_diffeqs, t, vars, params, name = gensym(:Parameterized))
    ODEFunctionExpr(sys, tgrad = opts[:build_tgrad], jac = opts[:build_jac])
end
