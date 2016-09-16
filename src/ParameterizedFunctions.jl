module ParameterizedFunctions
import Base.getindex

### Basic Functionality

abstract ParameterizedFunction <: Function
getindex{s}(p::ParameterizedFunction,::Val{s}) = getfield(p,s) ## Val for type-stability

### Macros

macro ode_def(name,ex,params...)
  ## Build independent variable dictionary
  dict = Dict{Symbol,Int}()
  for i in 2:2:length(ex.args) #Every odd line is line number
    arg = ex.args[i].args[1] #Get the first thing, should be dsomething
    nodarg = Symbol(string(arg)[2:end]) #Take off the d
    if !haskey(dict,nodarg)
      s = string(arg)
      dict[Symbol(string(arg)[2:end])] = i/2 # and label it the next int if not seen before
    end
  end
  syms = keys(dict)

  pdict = Dict{Symbol,Any}(); idict = Dict{Symbol,Any}()
  ## Build parameter and inline dictionaries
  for i in 1:length(params)
    if params[i].head == :(=>)
      pdict[params[i].args[1]] = params[i].args[2] # works for k=3, or k=>3
    elseif params[i].head == :(=)
      idict[params[i].args[1]] = params[i].args[2] # works for k=3, or k=>3
    end
  end
  # Run find replace to make the function expression
  ode_findreplace(ex,dict,syms,pdict,idict)
  push!(ex.args,nothing) # Make the return void
  # Build the type
  f = maketype(name,pdict)
  # Overload the Call
  overloadex = :(((p::$name))(t,u,du) = $ex)
  @eval $overloadex
  # Export the type
  exportex = :(export $name)
  @eval $exportex
  return f
end

function ode_findreplace(ex,dict,syms,pdict,idict)
  for (i,arg) in enumerate(ex.args)
    if isa(arg,Expr)
      ode_findreplace(arg,dict,syms,pdict,idict)
    elseif isa(arg,Symbol)
      s = string(arg)
      if haskey(dict,arg)
        ex.args[i] = :(u[$(dict[arg])]) # replace with u[i]
      elseif haskey(idict,arg)
        ex.args[i] = :($(idict[arg])) # inline from idict
      elseif haskey(pdict,arg)
        ex.args[i] = :(p.$arg) # replace with p.arg
      elseif length(string(arg))>1 && haskey(dict,Symbol(s[nextind(s, 1):end])) && Symbol(s[1])==:d
        tmp = Symbol(s[nextind(s, 1):end]) # Remove the first letter, the d
        ex.args[i] = :(du[$(dict[tmp])])
      end
    end
  end
end

function maketype(name, pdict)
    @eval type $name <: ParameterizedFunction
        $((:($x::$(typeof(t))) for (x, t) in pdict)...)
    end
    eval(name)(values(pdict)...)
end


export ParameterizedFunction, @ode_def
end # module
