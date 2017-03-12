macro fem_def(sig,name,ex,params...)

  origex = ex
  ## Build Symbol dictionary
  indvar_dict,syms = build_indvar_dict(ex)

  param_dict, inline_dict = build_paramdicts(params)
  # Run find replace
  fem_findreplace(ex,indvar_dict,syms,param_dict,inline_dict)
  fex = ex
  funcs = Vector{Expr}(0) # Get all of the functions
  for (i,arg) in enumerate(ex.args)
    if i%2 == 0
      push!(funcs,arg.args[2])
    end
  end
  if length(syms)==1
    ex = funcs[1]
  else
    ex = Expr(:hcat,funcs...)
  end

  exprs = Vector{Expr}(0)

  # Build the type
  typeex,constructorex = maketype(name,param_dict,origex,funcs,syms,fex)

  push!(exprs,typeex)
  push!(exprs,constructorex)

  if typeof(sig) == Symbol
    overloadex = :(((p::$name))($(sig)) = $ex) |> esc
  else
    overloadex = :(((p::$name))($(sig.args...)) = $ex) |> esc
  end
  # Overload the Call

  push!(exprs,overloadex)
  # Return the type from the default consturctor
  def_const_ex = :(($name)()) |> esc
  push!(exprs,def_const_ex)
  expr_arr_to_block(exprs)
end

function fem_findreplace(ex,indvar_dict,syms,param_dict,inline_dict)
  for (i,arg) in enumerate(ex.args)
    if isa(arg,Expr)
      fem_findreplace(arg,indvar_dict,syms,param_dict,inline_dict)
    elseif isa(arg,Symbol)
      if haskey(indvar_dict,arg)
        ex.args[i] = :(u[:,$(indvar_dict[arg])])
      elseif haskey(inline_dict,arg)
        ex.args[i] = :($(inline_dict[arg])) # Inline if in inline_dict
      elseif haskey(param_dict,arg)
        ex.args[i] = :(p.$arg) # replace with p.arg
      elseif haskey(FEM_SYMBOL_DICT,arg)
        ex.args[i] = FEM_SYMBOL_DICT[arg]
      end
    end
  end
end
