module ParameterizedFunctions
using SymEngine
import Base: getindex,ctranspose

### Basic Functionality

abstract ParameterizedFunction <: Function
getindex{s}(p::ParameterizedFunction,::Val{s}) = getfield(p,s) ## Val for type-stability

### Macros

macro ode_def(name,ex,params...)
  ## Build independent variable dictionary
  dict = Dict{Symbol,Int}()
  syms = Vector{Symbol}(0)
  for i in 2:2:length(ex.args) #Every odd line is line number
    arg = ex.args[i].args[1] #Get the first thing, should be dsomething
    nodarg = Symbol(string(arg)[2:end]) #Take off the d
    if !haskey(dict,nodarg)
      s = string(arg)
      dict[Symbol(string(arg)[2:end])] = i/2 # and label it the next int if not seen before
      push!(syms,Symbol(string(arg)[2:end]))
    end
  end

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
  symex = copy(ex) # Different expression for symbolic computations
  ode_findreplace(ex,symex,dict,pdict,idict)
  push!(ex.args,nothing) # Make the return void
  # Build the type
  f = maketype(name,pdict)
  # Overload the Call
  overloadex = :(((p::$name))(t,u,du) = $ex)
  @eval $overloadex
  # Export the type
  exportex = :(export $name)
  @eval $exportex

  # Now do the Jacobian. Get the component functions
  funcs = Vector{Expr}(0) # Get all of the functions for symbolic computation
  for (i,arg) in enumerate(symex.args)
    if i%2 == 0
      ex = arg.args[2]
      if typeof(ex) <: Symbol
        push!(funcs,:(1*$ex))
      else # It's an expression, just push
        push!(funcs,arg.args[2])
      end
    end
  end
  # Declare the symbols

  symstr = symarr_to_sympy(syms)
  symdefineex = Expr(:(=),parse("("*symstr*")"),SymEngine.symbols(symstr))
  symtup = parse("("*symstr*")")
  @eval $symdefineex
  symtup = @eval $symtup # symtup is the tuple of symPy symbols
    # Build the Jacobian Matrix of SymPy Expressions
  numsyms = length(symtup)
  symjac = Matrix(numsyms,numsyms)
  println("here")
  for i in eachindex(funcs)
    funcex = funcs[i]
    println(funcs[i])
    symfunc = @eval $funcex
    for j in eachindex(symtup)
      symjac[i,j] = diff(symfunc,symtup[j])
    end
  end
  println("here2")
  # Build the Julia function
  Jex = :()
  for i in 1:numsyms
    for j in 1:numsyms
      ex = parse(string(symjac[i,j]))
      if typeof(ex) <: Expr
        ode_findreplace(ex,ex,dict,pdict,idict)
      else
        ex = ode_symbol_findreplace(ex,dict,pdict,idict)
      end
      push!(Jex.args,:(J[$i,$j] = $ex))
    end
  end
  Jex.head = :block
  push!(Jex.args,nothing)
  Jex = :(jac = (t,u,J)->$Jex)
  @eval Base.ctranspose(p::$name) = $Jex
  return f
end

function ode_findreplace(ex,symex,dict,pdict,idict)
  for (i,arg) in enumerate(ex.args)
    if isa(arg,Expr)
      ode_findreplace(arg,symex.args[i],dict,pdict,idict)
    elseif isa(arg,Symbol)
      s = string(arg)
      if haskey(dict,arg)
        ex.args[i] = :(u[$(dict[arg])]) # replace with u[i]
      elseif haskey(idict,arg)
        ex.args[i] = :($(idict[arg])) # inline from idict
        symex.args[i] = :($(idict[arg])) # also do in symbolic
      elseif haskey(pdict,arg)
        ex.args[i] = :(p.$arg) # replace with p.arg
        symex.args[i] = :($(pdict[arg])) # also do in symbolic
      elseif length(string(arg))>1 && haskey(dict,Symbol(s[nextind(s, 1):end])) && Symbol(s[1])==:d
        tmp = Symbol(s[nextind(s, 1):end]) # Remove the first letter, the d
        ex.args[i] = :(du[$(dict[tmp])])
        symex.args[i] = :(du[$(dict[tmp])]) #also do in symbolic
      end
    end
  end
end

function ode_symbol_findreplace(ex,dict,pdict,idict)
  if haskey(dict,ex)
    ex = :(u[$(dict[ex])]) # replace with u[i]
  end
  :(1*$ex) # Add the 1 to make it an expression not a Symbol
end

function maketype(name, pdict)
    @eval type $name <: ParameterizedFunction
        $((:($x::$(typeof(t))) for (x, t) in pdict)...)
    end
    eval(name)(values(pdict)...)
end

macro fem_def(sig,name,ex,params...)
  ## Build Symbol dictionary
  dict = Dict{Symbol,Int}()
  for (i,arg) in enumerate(ex.args)
    if i%2 == 0
      dict[Symbol(string(arg.args[1])[2:end])] = i/2 # Change du->u, Fix i counting
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
  # Run find replace
  fem_findreplace(ex,dict,syms,pdict,idict)
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
  # Build the type
  f = maketype(name,pdict)
  # Overload the Call
  newsig = :($(sig.args...))
  overloadex = :(((p::$name))($(sig.args...)) = $ex)
  @eval $overloadex
  # Export the type
  exportex = :(export $name)
  @eval $exportex
  return f
end

function fem_findreplace(ex,dict,syms,pdict,idict)
  for (i,arg) in enumerate(ex.args)
    if isa(arg,Expr)
      fem_findreplace(arg,dict,syms,pdict,idict)
    elseif isa(arg,Symbol)
      if haskey(dict,arg)
        ex.args[i] = :(u[:,$(dict[arg])])
      elseif haskey(idict,arg)
        ex.args[i] = :($(idict[arg])) # Inline if in idict
      elseif haskey(pdict,arg)
        ex.args[i] = :(p.$arg) # replace with p.arg
      elseif haskey(FEM_SYMBOL_DICT,arg)
        ex.args[i] = FEM_SYMBOL_DICT[arg]
      end
    end
  end
end


### Utility Functions

"""
symarr_to_sympy(symarr::Vector{Symbol})

Converts a Vector{Symbol} into a string for sympy parsing

Symbol[:x,:y] --> "x,y"
"""
function symarr_to_sympy(symarr::Vector{Symbol})
  str = ""
  for sym in symarr
    str = str*string(sym)*","
  end
  str[1:end-1]
end

const FEM_SYMBOL_DICT = Dict{Symbol,Expr}(:x=>:(x[:,1]),:y=>:(x[:,2]),:z=>:(x[:,3]))

export ParameterizedFunction, @ode_def, @fem_def
end # module
