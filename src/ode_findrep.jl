function ode_findreplace(ex,symex,indvar_dict,param_dict,inline_dict;params_from_function=true,
                         vectorized_form=false,vectorized_returns=:standard)
  for (i,arg) in enumerate(ex.args)
    if isa(arg,Expr)
      ode_findreplace(arg,symex.args[i],indvar_dict,param_dict,inline_dict;
                      params_from_function=params_from_function,vectorized_form=vectorized_form,
                      vectorized_returns=vectorized_returns)
    elseif isa(arg,Symbol)
      s = string(arg)
      if haskey(indvar_dict,arg)
        if vectorized_form
          ex.args[i] = :(u[:,$(indvar_dict[arg])]) # replace with u[:,i]
        else
          ex.args[i] = :(u[$(indvar_dict[arg])]) # replace with u[i]
        end
      elseif haskey(inline_dict,arg)
        ex.args[i] = :($(inline_dict[arg])) # inline from inline_dict
        symex.args[i] = :($(inline_dict[arg])) # also do in symbolic
      elseif haskey(param_dict,arg)
        if params_from_function
          ex.args[i] = :(p.$arg) # replace with p.arg
        else
          idx = findfirst(param_dict.keys .== arg)
          ex.args[i] = :(params[$idx]) # replace with params[arg]
        end
        symex.args[i] = arg # keep arg
      elseif length(string(arg))>1 && haskey(indvar_dict,Symbol(s[nextind(s, 1):end])) && Symbol(s[1])==:d
        tmp = Symbol(s[nextind(s, 1):end]) # Remove the first letter, the d
        if vectorized_returns == :slice
          ex.args[i] = :(du[:,$(indvar_dict[tmp])])
          symex.args[i] = :(du[:,$(indvar_dict[tmp])]) #also do in symbolic
        elseif vectorized_returns == :vals
          ex.args[i] = Symbol("du$(indvar_dict[tmp])")
          symex.args[i] = Symbol("du$(indvar_dict[tmp])") #also do in symbolic
        else # vectorized_returns == :standard
          ex.args[i] = :(du[$(indvar_dict[tmp])])
          symex.args[i] = :(du[$(indvar_dict[tmp])]) #also do in symbolic
        end
      end
    end
  end
end

function bad_derivative(ex)
  for (i,arg) in enumerate(ex.args)
    if isa(arg,Expr)
      bad_derivative(arg)
    elseif arg == :Derivative
      warn("Undefined derivative found. If you are using a non-elementary function, you must define the derivative in order to calculate Jacobians et. al. Please refer to the documentation.")
      error("Failed")
    end
  end
end



function ode_symbol_findreplace(ex,indvar_dict,param_dict,inline_dict;params_from_function=true)
  if haskey(indvar_dict,ex)
    ex = :(u[$(indvar_dict[ex])]) # replace with u[i]
  elseif haskey(param_dict,ex)
    if params_from_function
      ex = :(p.$ex) # replace with u[i]
    else
      idx = findfirst(param_dict.keys .== ex)
      ex = :(params[$idx])
    end
  end
  :($ex*1) # Add the 1 to make it an expression not a Symbol
end


function flip_mult!(ex)
  for (i,arg) in enumerate(ex.args)
    if isa(arg,Expr)
      if arg.args[1] == :(*) && length(arg.args) >= 3 && (isa(arg.args[2], Number) || (isa(arg.args[2], Expr) && arg.args[2].args[1] == :-))
        arg.args[3],arg.args[2] = arg.args[2], arg.args[3]
      else
        flip_mult!(arg)
      end
    end
  end
end
