function build_jac_func(symjac,indvar_dict,param_dict,inline_dict;params_from_function=true)
  Jex = :()
  for i in 1:size(symjac,1)
    for j in 1:size(symjac,2)
      ex = parse(string(symjac[i,j]))
      if typeof(ex) <: Expr
        ode_findreplace(ex,copy(ex),indvar_dict,param_dict,inline_dict,params_from_function=params_from_function)
      else
        ex = ode_symbol_findreplace(ex,indvar_dict,param_dict,inline_dict,params_from_function=params_from_function)
      end
      push!(Jex.args,:(J[$i,$j] = $ex))
    end
  end
  Jex.head = :block
  push!(Jex.args,nothing)
  Jex
end

function build_tgrad_func(symtgrad,indvar_dict,param_dict,inline_dict;params_from_function=false)
  tgradex = :()
  for i in 1:length(symtgrad)
    ex = parse(string(symtgrad[i]))
    if typeof(ex) <: Expr
      ode_findreplace(ex,copy(ex),indvar_dict,param_dict,inline_dict,params_from_function=params_from_function)
    else
      ex = ode_symbol_findreplace(ex,indvar_dict,param_dict,inline_dict,params_from_function=params_from_function)
    end
    push!(tgradex.args,:(grad[$i] = $ex))
  end
  tgradex.head = :block
  push!(tgradex.args,nothing)
  tgradex
end

function build_p_funcs(paramfuncs,indvar_dict,param_dict,inline_dict)
  params = param_dict.keys
  pfuncs = Vector{Expr}(length(params))
  param_dict_type = typeof(param_dict)
  for i in 1:length(params)
    pfunc = :()
    param = params[i]
    param_dict_drop_cur = deepcopy(param_dict)
    delete!(param_dict_drop_cur,param)
    for j in 1:length(paramfuncs[1])
      ex = paramfuncs[i][j]
      if typeof(ex) <: Expr
        ode_findreplace(ex,copy(ex),indvar_dict,param_dict_drop_cur,inline_dict)
      else
        ex = ode_symbol_findreplace(ex,indvar_dict,param_dict_drop_cur,inline_dict)
      end
      push!(pfunc.args,:(du[$j] = $ex))
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
        push!(funcs,:(1*$ex))
      else # It's an expression, just push
        push!(funcs,arg.args[2])
      end
    end
  end
  funcs
end
