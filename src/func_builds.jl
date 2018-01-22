function build_jac_func(symjac,indvar_dict,params;params_from_function=true)
  Jex = :()
  for i in 1:size(symjac,1)
    for j in 1:size(symjac,2)
      ex = parse(string(symjac[i,j]))
      if typeof(ex) <: Expr
        ode_findreplace(ex,copy(ex),indvar_dict,params,params_from_function=params_from_function)
      else
        ex = ode_symbol_findreplace(ex,indvar_dict,params,params_from_function=params_from_function)
      end
      push!(Jex.args,:(internal_var___J[$i,$j] = $ex))
    end
  end
  Jex.head = :block
  push!(Jex.args,nothing)
  Jex
end

function build_tgrad_func(symtgrad,indvar_dict,params;params_from_function=false)
  tgradex = :()
  for i in 1:length(symtgrad)
    ex = parse(string(symtgrad[i]))
    if typeof(ex) <: Expr
      ode_findreplace(ex,copy(ex),indvar_dict,params,params_from_function=params_from_function)
    else
      ex = ode_symbol_findreplace(ex,indvar_dict,params,params_from_function=params_from_function)
    end
    push!(tgradex.args,:(internal_var___grad[$i] = $ex))
  end
  tgradex.head = :block
  push!(tgradex.args,nothing)
  tgradex
end

function build_p_funcs(paramfuncs,indvar_dict,params)
  pfuncs = Vector{Expr}(length(params))
  params_type = typeof(params)
  for i in 1:length(params)
    pfunc = :()
    params_drop_cur = copy(params)
    deleteat!(params_drop_cur,i)
    for j in 1:length(paramfuncs[1])
      ex = paramfuncs[i][j]
      if typeof(ex) <: Expr
        ode_findreplace(ex,copy(ex),indvar_dict,params_drop_cur)
      else
        ex = ode_symbol_findreplace(ex,indvar_dict,params_drop_cur)
      end
      push!(pfunc.args,:(internal_var___du[$j] = $ex))
    end
    pfunc.head = :block
    push!(pfunc.args,nothing)
    pfuncs[i] = pfunc
  end
  pfuncs
end

function build_component_funcs(symex)
  funcs = Vector{Expr}(0) # Get all of the functions for symbolic computation
  for (i,arg) in enumerate(symex.args)
    if i%2 == 0
      ex = arg.args[2]
      if (typeof(ex) <: Symbol) || (typeof(ex) <: Number)
        push!(funcs,:($ex*1))
      else # It's an expression, just push
        fix_ex = copy(ex)
        flip_mult!(fix_ex)
        push!(funcs,fix_ex)
      end
    end
  end
  funcs
end
