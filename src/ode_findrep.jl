function ode_findreplace(ex, symex, indvar_dict, params; params_from_function = true,
                         vectorized_form = false, vectorized_returns = :standard)
    for (i, arg) in enumerate(ex.args)
        if isa(arg, Expr)
            ode_findreplace(arg, symex.args[i], indvar_dict, params;
                            params_from_function = params_from_function,
                            vectorized_form = vectorized_form,
                            vectorized_returns = vectorized_returns)
        elseif isa(arg, Symbol)
            s = string(arg)
            if haskey(indvar_dict, arg)
                if vectorized_form
                    ex.args[i] = :(internal_var___u[:, $(indvar_dict[arg])]) # replace with internal_var___u[:,i]
                else
                    ex.args[i] = :(internal_var___u[$(indvar_dict[arg])]) # replace with internal_var___u[i]
                end
            elseif arg in params
                idx = findfirst(params .== arg)
                ex.args[i] = :(internal_var___p[$idx]) # replace with params[arg]
                symex.args[i] = arg # keep arg
            elseif length(string(arg)) > 1 &&
                   haskey(indvar_dict, Symbol(s[nextind(s, 1):end])) && Symbol(s[1]) == :d
                tmp = Symbol(s[nextind(s, 1):end]) # Remove the first letter, the d
                if vectorized_returns == :slice
                    ex.args[i] = :(internal_var___du[:, $(indvar_dict[tmp])])
                    symex.args[i] = :(internal_var___du[:, $(indvar_dict[tmp])]) #also do in symbolic
                elseif vectorized_returns == :vals
                    ex.args[i] = Symbol("internal_var___du$(indvar_dict[tmp])")
                    symex.args[i] = Symbol("internal_var___du$(indvar_dict[tmp])") #also do in symbolic
                else # vectorized_returns == :standard
                    ex.args[i] = :(internal_var___du[$(indvar_dict[tmp])])
                    symex.args[i] = :(internal_var___du[$(indvar_dict[tmp])]) #also do in symbolic
                end
            end
        end
    end
end

function bad_derivative(ex)
    for (i, arg) in enumerate(ex.args)
        if isa(arg, Expr)
            bad_derivative(arg)
        elseif arg == :Derivative
            warn("Undefined derivative found. If you are using a non-elementary function, you must define the derivative in order to calculate Jacobians et. al. Please refer to the documentation.")
            error("Failed")
        end
    end
end

function ode_symbol_findreplace(ex, indvar_dict, params; params_from_function = true)
    if haskey(indvar_dict, ex)
        ex = :(internal_var___u[$(indvar_dict[ex])]) # replace with internal_var___u[i]
    elseif ex in params
        idx = findfirst(params .== ex)
        ex = :(internal_var___p[$idx])
    end
    :($ex * 1) # Add the 1 to make it an expression not a Symbol
end

function flip_mult!(ex)
    for (i, arg) in enumerate(ex.args)
        if isa(arg, Expr)
            if arg.args[1] == :(*) && length(arg.args) >= 3 &&
               (isa(arg.args[2], Number) ||
                (isa(arg.args[2], Expr) && arg.args[2].args[1] == :-))
                arg.args[3], arg.args[2] = arg.args[2], arg.args[3]
            else
                flip_mult!(arg)
            end
        end
    end
end
