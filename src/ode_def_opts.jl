findreplace(ex::Symbol, dict) = get(dict, ex, ex)
function findreplace(ex::Expr, dict)
    Expr(ex.head, map(x->findreplace(x, dict), ex.args)...)
end
findreplace(ex, dict) = ex

function ode_def_opts(name::Symbol,opts::Dict{Symbol,Bool},curmod,ex::Expr,params...;depvar=:t)
  # depvar is the dependent variable. Defaults to t
  # M is the mass matrix in RosW, must be a constant!

  origex = copy(ex) # Save the original expression

  if !(eltype(params) <: Symbol)
    error("The syntax for ParameterizedFunctions has changed. Simply list the parameters at the end, i.e. `a b c d`, instead of `a=5.0 b=>3.0 ...`. Parameters are defined in the problem type. See the documentation for more information.")
  end
  params = Symbol[params...]

  ## Build independent variable dictionary
  indvar_dict,syms = build_indvar_dict(ex,depvar)
  ####


  t = Sym{Number}(:t)
  vars = [Sym{FnType{Tuple, Number}}(x)(t) for x in syms]
  params = [Sym{Number}(x) for x in Symbol[params...]]

  vars_dict = Dict(x=>Symbol(v) for (x, v) in zip(syms, vars))

  # replace x with x(t) if it's a var
  ex = findreplace(ex, vars_dict)

  # Build the Expressions

  # Run find replace to make the function expression
  symex = copy(ex) # Different expression for symbolic computations
  #ode_findreplace(ex,symex,indvar_dict,params)
  funcs = build_component_funcs(symex)
  mtk_ops = modelingtoolkitize_expr.(funcs,([t;vars;params],),(curmod,))

  D = ModelingToolkit.Differential(t)

  mtk_diffeqs = [D(vars[i]) ~ mtk_ops[i] for i in 1:length(vars)]

  sys = ODESystem(mtk_diffeqs,t,vars,params)

  f_ex_oop,f_ex_iip = ModelingToolkit.generate_function(sys,vars,params)

  if opts[:build_tgrad]
    try
      tgrad_ex_oop,tgrad_ex_iip = ModelingToolkit.generate_tgrad(sys,vars,params)
    catch
      @warn "tgrad construction failed"
      tgrad_ex_oop,tgrad_ex_iip = nothing,nothing
    end
  else
    tgrad_ex_oop,tgrad_ex_iip = nothing,nothing
  end

  if opts[:build_jac]
    try
      J_ex_oop,J_ex_iip = ModelingToolkit.generate_jacobian(sys,vars,params)
    catch
      @warn "Jacobian construction failed"
      J_ex_oop,J_ex_iip = nothing,nothing
    end
  else
    J_ex_oop,J_ex_iip = nothing,nothing
  end

  if opts[:build_invW] && length(mtk_diffeqs) < 4
    try
      W_exs = ModelingToolkit.generate_factorized_W(sys,vars,params,false)
      W_ex_oop,W_ex_iip = W_exs[1]
      W_t_ex_oop,W_t_ex_iip = W_exs[2]
    catch
      @warn "W-expression construction failed"
      W_ex_oop,W_ex_iip = (nothing,nothing)
      W_t_ex_oop,W_t_ex_iip = (nothing,nothing)
    end
  else
    W_ex_oop,W_ex_iip = (nothing,nothing)
    W_t_ex_oop,W_t_ex_iip = (nothing,nothing)
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

  quote
    struct $name{F,TG,TJ,TW,TWt,S} <: ParameterizedFunctions.DiffEqBase.AbstractParameterizedFunction{true}
      f::F
      mass_matrix::ParameterizedFunctions.LinearAlgebra.UniformScaling
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
    end

    (f::$name)(args...) = f.f(args...)

    $fname($(f_ex_oop.args[1].args...)) = $(f_ex_oop.args[2])
    $fname($(f_ex_iip.args[1].args...)) = $(f_ex_iip.args[2])
    $full_tex
    $full_jex
    $full_wex

    $name($fname,ParameterizedFunctions.LinearAlgebra.I,nothing,$tname,$jname,nothing,nothing,
          nothing,nothing,$Wname,$W_tname,nothing,$syms,$depvar,nothing,$sys)
  end |> esc
end
