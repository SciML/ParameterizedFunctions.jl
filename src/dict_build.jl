function build_indvar_dict(ex,depvar)
  indvar_dict = OrderedDict{Symbol,Int}()
  cur_sym = 0
  for i in 2:2:length(ex.args) #Every odd line is line number
    arg = ex.args[i].args[1] #Get the first thing, should be dsomething
    firstarg = Symbol(first(string(arg))) # Check for d
    if firstarg == :d
      nodarg = Symbol(join(Base.Iterators.drop(string(arg), 1)))
      if nodarg == depvar
        warn("$depvar is fixed as the independent variable but is also used as a dependent variable. Results my be incorrect.")
      end
      if !haskey(indvar_dict,nodarg)
        cur_sym += 1
        indvar_dict[nodarg] = cur_sym
      else
        error("The derivative term for $nodarg is repeated. This is not allowed.")
      end
    end
  end
  syms = indvar_dict.keys
  indvar_dict,syms
end

function build_paramdicts(params)
  param_dict = OrderedDict{Symbol,Any}(); inline_dict = OrderedDict{Symbol,Any}()
  for i in 1:length(params)
    if VERSION < v"0.6.0-dev.2613" && params[i].head == :(=>)
      param_dict[params[i].args[1]] = params[i].args[2] # works for k=3, or k=>3
    elseif params[i].head == :call
      param_dict[params[i].args[2]] = params[i].args[3] # works for k=3, or k=>3
    elseif params[i].head == :(=)
      inline_dict[params[i].args[1]] = params[i].args[2] # works for k=3, or k=>3
    end
  end
  param_dict,inline_dict
end
