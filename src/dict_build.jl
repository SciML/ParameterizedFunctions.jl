function build_indvar_dict(ex, depvar)
    indvar_dict = OrderedDict{Symbol, Int}()
    cur_sym = 0
    for i in 2:2:length(ex.args) #Every odd line is line number
        arg = ex.args[i].args[1] #Get the first thing, should be dsomething
        firstarg = Symbol(first(string(arg))) # Check for d
        if firstarg == :d
            nodarg = Symbol(join(Base.Iterators.drop(string(arg), 1)))
            if nodarg == depvar
                warn("$depvar is fixed as the independent variable but is also used as a dependent variable. Results my be incorrect.")
            end
            if !haskey(indvar_dict, nodarg)
                cur_sym += 1
                indvar_dict[nodarg] = cur_sym
            else
                error("The derivative term for $nodarg is repeated. This is not allowed.")
            end
        end
    end
    syms = indvar_dict.keys
    indvar_dict, syms
end

function build_param_list(params)
    param_list = Vector{Symbol}()
    for i in 1:length(params)
        if params[i] isa Symbol
            push!(param_list, params[i])
        elseif params[i].head == :call || params[i].head == :(=)
            warn("p=>val and p=val are deprecated. Simply list the parameters. See the DifferentialEquations.jl documentation for more information on the syntax change.")
        end
    end
    param_dict
end
