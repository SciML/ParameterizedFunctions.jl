findreplace(ex::Symbol, dict) = get(dict, ex, ex)
function findreplace(ex::Expr, dict)
    return Expr(ex.head, map(x -> findreplace(x, dict), ex.args)...)
end
findreplace(ex, dict) = ex

"""
    ode_def_opts(name::Symbol, opts::Dict{Symbol, Bool}, curmod,
        ex::Expr, params...; depvar = :t)

Build the expression emitted by the `@ode_def` family of macros.

This is a developer API for packages that implement an ODE-definition macro. Most
users should call `@ode_def`, `@ode_def_bare`, or `@ode_def_all` instead. Callers
must pass the assignment-based DSL expression accepted by those macros; this
function returns an expression and does not evaluate the generated definition.

# Arguments

- `name`: name of the generated `AbstractParameterizedFunction` subtype.
- `opts`: derivative-generation options. `:build_tgrad`, `:build_jac`, and
  `:build_invW` control generated functions. The compatibility keys
  `:build_expjac`, `:build_invjac`, `:build_invW_t`, `:build_hes`,
  `:build_invhes`, and `:build_dpfuncs` must also be present and have `Bool`
  values, but do not currently change the generated expression.
- `curmod`: module in which functions named in `ex` are resolved.
- `ex`: `begin` expression containing `dstate = expression` assignments.
- `params...`: parameter symbols used in `ex`.

# Keywords

- `depvar = :t`: symbol used for the independent variable. It must not also name a
  dependent state.

# Returns

An expression that defines a concrete `AbstractParameterizedFunction` subtype,
generated callable methods, and a singleton instance.

# Examples

```julia
opts = Dict{Symbol, Bool}(
    :build_tgrad => true, :build_jac => true, :build_expjac => false,
    :build_invjac => false, :build_invW => false, :build_invW_t => false,
    :build_hes => false, :build_invhes => false, :build_dpfuncs => true,
)
```
"""
function ode_def_opts(
        name::Symbol, opts::Dict{Symbol, Bool}, curmod, ex::Expr, params...;
        depvar = :t
    )
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

    t = ModelingToolkitBase.t_nounits
    vars = SymbolicUtils.unwrap.([(@variables $x(t))[1] for x in syms])
    params = SymbolicUtils.unwrap.([(@parameters $x)[1] for x in Symbol[params...]])

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

    sys = System(mtk_diffeqs, t, vars, params, name = gensym(:Parameterized))

    f_ex_oop, f_ex_iip = ModelingToolkit.generate_rhs(sys)

    if opts[:build_tgrad]
        try
            tgrad_ex_oop, tgrad_ex_iip = ModelingToolkit.generate_tgrad(sys)
        catch
            @warn "tgrad construction failed"
            tgrad_ex_oop, tgrad_ex_iip = nothing, nothing
        end
    else
        tgrad_ex_oop, tgrad_ex_iip = nothing, nothing
    end

    if opts[:build_jac]
        try
            J_ex_oop, J_ex_iip = ModelingToolkit.generate_jacobian(sys)
        catch
            @warn "Jacobian construction failed"
            J_ex_oop, J_ex_iip = nothing, nothing
        end
    else
        J_ex_oop, J_ex_iip = nothing, nothing
    end

    if opts[:build_invW] && length(mtk_diffeqs) < 4
        try
            W_exs = ModelingToolkit.generate_factorized_W(sys, vars, params, false)
            W_ex_oop, W_ex_iip = W_exs[1]
            W_t_ex_oop, W_t_ex_iip = W_exs[2]
        catch
            @warn "W-expression construction failed"
            W_ex_oop, W_ex_iip = (nothing, nothing)
            W_t_ex_oop, W_t_ex_iip = (nothing, nothing)
        end
    else
        W_ex_oop, W_ex_iip = (nothing, nothing)
        W_t_ex_oop, W_t_ex_iip = (nothing, nothing)
    end

    fname = gensym(:ParameterizedDiffEqFunction)
    tname = gensym(:ParameterizedTGradFunction)
    jname = gensym(:ParameterizedJacobianFunction)
    Wname = gensym(:ParameterizedWFactFunction)
    W_tname = gensym(:ParameterizedW_tFactFunction)
    funcname = gensym(:ParameterizedODEFunction)

    if tgrad_ex_oop !== nothing
        full_tex = quote
            $tname($(tgrad_ex_oop.args[1].args...)) = $(tgrad_ex_oop.args[2])
            $tname($(tgrad_ex_iip.args[1].args...)) = $(tgrad_ex_iip.args[2])
        end
    else
        full_tex = quote
            $tname = nothing
        end
    end

    if J_ex_oop !== nothing
        full_jex = quote
            $jname($(J_ex_oop.args[1].args...)) = $(J_ex_oop.args[2])
            $jname($(J_ex_iip.args[1].args...)) = $(J_ex_iip.args[2])
        end
    else
        full_jex = quote
            $jname = nothing
        end
    end

    if W_ex_oop !== nothing
        full_wex = quote
            $Wname($(W_ex_oop.args[1].args...)) = $(W_ex_oop.args[2])
            $Wname($(W_ex_iip.args[1].args...)) = $(W_ex_iip.args[2])
            $W_tname($(W_t_ex_oop.args[1].args...)) = $(W_t_ex_oop.args[2])
            $W_tname($(W_t_ex_iip.args[1].args...)) = $(W_t_ex_iip.args[2])
        end
    else
        full_wex = quote
            $Wname = nothing
            $W_tname = nothing
        end
    end

    return quote
        struct $name{F, TG, TJ, TW, TWt, S} <:
            ParameterizedFunctions.SciMLBase.AbstractParameterizedFunction{true}
            f::F
            mass_matrix::ParameterizedFunctions.LinearAlgebra.UniformScaling{Bool}
            analytic::Nothing
            tgrad::TG
            jac::TJ
            jvp::Nothing
            vjp::Nothing
            jac_prototype::Nothing
            sparsity::Nothing
            Wfact::TW
            Wfact_t::TWt
            paramjac::Nothing
            syms::Vector{Symbol}
            indepvar::Symbol
            colorvec::Nothing
            sys::S
            initialization_data::Nothing
            nlprob_data::Nothing
        end

        (f::$name)(args...) = f.f(args...)

        function ParameterizedFunctions.SciMLBase.remake(func::$name; kwargs...)
            return func
        end

        $fname($(f_ex_oop.args[1].args...)) = $(f_ex_oop.args[2])
        $fname($(f_ex_iip.args[1].args...)) = $(f_ex_iip.args[2])
        $full_tex
        $full_jex
        $full_wex

        $name(
            $fname, ParameterizedFunctions.LinearAlgebra.I, nothing, $tname, $jname,
            nothing, nothing,
            nothing, nothing, $Wname, $W_tname, nothing, $syms, $(Meta.quot(depvar)),
            nothing, $sys, nothing, nothing
        )
    end |> esc
end
