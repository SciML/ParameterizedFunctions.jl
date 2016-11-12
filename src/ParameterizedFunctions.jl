module ParameterizedFunctions
using SymEngine, DataStructures, DiffEqBase
import Base: getindex

### Macros

jac_exists(f) = !isempty(methods(f,Tuple{Type{Val{:Jac}}, Any,Any,Any}))
invjac_exists(f) = !isempty(methods(f,Tuple{Type{Val{:InvJac}}, Any,Any,Any}))
hes_exists(f) = !isempty(methods(f,Tuple{Type{Val{:Hes}}, Any,Any,Any}))
invhes_exists(f) = !isempty(methods(f,Tuple{Type{Val{:InvHes}}, Any,Any,Any}))
paramjac_exists(f) = !isempty(methods(f,Tuple{Type{Val{:param_Jac}}, Any,Any,Any,Any}))
pfunc_exists(f,p::Symbol) = !isempty(methods(f,Tuple{Type{Val{p}}, Any,Any,Any}))
pderiv_exists(f,p::Symbol) = !isempty(methods(f,Tuple{Type{Val{p}},Val{:Deriv},Any,Any,Any}))

pfunc_exists(f) = !isempty(methods(f,Tuple{Type{Val{f.params[1]}}, Any,Any,Any,Any}))
pderiv_exists(f) = !isempty(methods(f,Tuple{Type{Val{f.params[1]}},Val{:Deriv},Any,Any,Any}))

function ode_def_opts(name::Symbol,opts::Dict{Symbol,Bool},ex::Expr,params...)
  origex = copy(ex) # Save the original expression

  ## Build independent variable dictionary
  indvar_dict,syms = build_indvar_dict(ex)
  ## Build parameter and inline dictionaries
  param_dict, inline_dict = build_paramdicts(params)

  ####
  # Build the Expressions

  # Run find replace to make the function expression
  symex = copy(ex) # Different expression for symbolic computations
  ode_findreplace(ex,symex,indvar_dict,param_dict,inline_dict)
  push!(ex.args,nothing) # Make the return void
  fex = ex # Save this expression as the expression for the call

  # Parameter-Explicit Functions
  pex = copy(origex) # Build it from the original expression
  # Parameter find/replace
  ode_findreplace(pex,copy(ex),indvar_dict,param_dict,inline_dict;params_from_function=false)

  ######
  # Build the Functions

  # Get the component functions
  funcs = build_component_funcs(symex)

  # Declare the SymEngine symbols
  symtup,paramtup = symbolize(syms,param_dict.keys)
  numsyms = length(indvar_dict)
  numparams = length(param_dict)
  # Jacobian Calculation
  symjac = Matrix{SymEngine.Basic}(0,0)
  invjac = Matrix{SymEngine.Basic}(0,0)
  symhes = Matrix{SymEngine.Basic}(0,0)
  invhes = Matrix{SymEngine.Basic}(0,0)
  param_symjac = Matrix{SymEngine.Basic}(0,0)
  Jex = :(error("Jacobian Does Not Exist"))
  jac_exists = false
  invJex = :(error("Inverse Jacobian Does Not Exist"))
  invjac_exists = false
  Hex = :(error("Hessian Does Not Exist"))
  hes_exists = false
  invHex = :(error("Inverse Hessian Does Not Exist"))
  invhes_exists = false
  param_Jex = :(error("Parameter Jacobian Does Not Exist"))
  param_jac_exists = false

  if opts[:build_Jac]
    try #Jacobians and Hessian
      # Build the Jacobian Matrix of SymEngine Expressions
      symjac = Matrix{SymEngine.Basic}(numsyms,numsyms)
      for i in eachindex(funcs)
        funcex = funcs[i]
        symfunc = @eval $funcex
        for j in eachindex(symtup)
          symjac[i,j] = diff(SymEngine.Basic(symfunc),symtup[j])
        end
      end

      # Build the Julia function
      Jex = build_jac_func(symjac,indvar_dict,param_dict,inline_dict)
      jac_exists = true

      if opts[:build_InvJac]
        try # Jacobian Inverse
          invjac = inv(symjac)
          invJex = build_jac_func(invjac,indvar_dict,param_dict,inline_dict)
          invjac_exists = true
        catch err
          warn("Jacobian could not invert")
        end
      end
      if opts[:build_Hes]
        try # Hessian
          symhes = Matrix{SymEngine.Basic}(numsyms,numsyms)
          for i in eachindex(funcs), j in eachindex(symtup)
            symhes[i,j] = diff(symjac[i,j],symtup[j])
          end
          # Build the Julia function
          Hex = build_jac_func(symhes,indvar_dict,param_dict,inline_dict)
          hes_exists = true
          if opts[:build_InvHes]
            try # Hessian Inverse
              invhes = inv(symhes)
              invHex = build_jac_func(invhes,indvar_dict,param_dict,inline_dict)
              invhes_exists = true
            catch err
              warn("Hessian could not invert")
            end
          end
        end
      end
    catch err
      warn("Failed to build the Jacoboian. This means the Hessian is not built as well.")
    end
  end

  # Parameter Functions
  paramfuncs = Vector{Vector{Expr}}(numparams)
  for i in eachindex(paramtup)
    tmp_pfunc = Vector{Expr}(length(funcs))
    for j in eachindex(funcs)
      tmp_pfunc[j] = copy(funcs[j])
    end
    paramfuncs[i] = tmp_pfunc
  end
  pfuncs = build_p_funcs(paramfuncs,paramtup,indvar_dict,param_dict,inline_dict)

  d_pfuncs = Vector{Expr}(0)
  param_symjac = Matrix{SymEngine.Basic}(numsyms,numparams)
  pderiv_exists = false
  if opts[:build_dpfuncs]
    try # Parameter Gradients

      d_paramfuncs  = Vector{Vector{Expr}}(numparams)
      for i in eachindex(paramtup)
        tmp_dpfunc = Vector{Expr}(length(funcs))
        for j in eachindex(funcs)
          funcex = funcs[j]
          symfunc = @eval $funcex
          d_curr = diff(SymEngine.Basic(symfunc),paramtup[i])
          param_symjac[j,i] = d_curr
          symfunc_str = parse(string(d_curr))
          if typeof(symfunc_str) <: Number
            tmp_dpfunc[j] = :(1*$symfunc_str)
          elseif typeof(symfunc_str) <: Symbol
            tmp_dpfunc[j] = :(1*$symfunc_str)
          else
            tmp_dpfunc[j] = symfunc_str
          end
        end
        d_paramfuncs[i] = tmp_dpfunc
      end
      d_pfuncs = build_p_funcs(d_paramfuncs,paramtup,indvar_dict,param_dict,inline_dict)
      pderiv_exists = true

      # Now build the parameter Jacobian
      param_symjac_ex = Matrix{Expr}(numsyms,numparams)
      for i in 1:numparams
        param_symjac_ex[:,i] = d_paramfuncs[i]
      end

      param_Jex = build_jac_func(param_symjac_ex,indvar_dict,param_dict,inline_dict,params_from_function=false)
      param_jac_exists = true
    catch err
      warn("Failed to build the parameter derivatives.")
    end
  end
  # Build the type
  f = maketype(name,param_dict,origex,funcs,syms,fex,pex=pex,
               symjac=symjac,Jex=Jex,invjac=invjac,
               invJex=invJex,symhes=symhes,invhes=invhes,Hex=Hex,
               invHex=invHex,params=param_dict.keys,
               pfuncs=pfuncs,d_pfuncs=d_pfuncs,
               param_symjac=param_symjac,param_Jex=param_Jex)
  # Overload the Call
  overloadex = :(((p::$name))(t::Number,u,du) = $fex)
  @eval $overloadex
  # Value Dispatches for the Parameters
  for i in 1:length(paramtup)
    param = Symbol(paramtup[i])
    param_func = pfuncs[i]
    param_valtype = Val{param}
    overloadex = :(((p::$name))(::Type{$param_valtype},t,u,$param,du) = $param_func)
    @eval $overloadex
  end

  # Build the Function
  overloadex = :(((p::$name))(t::Number,u,du,params) = $pex)
  @eval $overloadex

  # Value Dispatches for the Parameter Derivatives
  if pderiv_exists
    for i in 1:length(paramtup)
      param = Symbol(paramtup[i])
      param_func = d_pfuncs[i]
      param_valtype = Val{param}
      overloadex = :(((p::$name))(::Type{$param_valtype},::Type{Val{:Deriv}},t,u,$param,du) = $param_func)
      @eval $overloadex
    end
  end

  # Add the Jacobian
  if jac_exists
    overloadex = :(((p::$name))(::Type{Val{:Jac}},t,u,J) = $Jex)
    @eval $overloadex
  end
  # Add the Inverse Jacobian
  if invjac_exists
    overloadex = :(((p::$name))(::Type{Val{:InvJac}},t,u,J) = $invJex)
    @eval $overloadex
  end
  # Add the Hessian
  if hes_exists
    overloadex = :(((p::$name))(::Type{Val{:Hes}},t,u,J) = $Hex)
    @eval $overloadex
  end
  # Add the Inverse Hessian
  if invhes_exists
    overloadex = :(((p::$name))(::Type{Val{:InvHes}},t,u,J) = $invHex)
    @eval $overloadex
  end
  # Add Parameter Jacobian
  if param_jac_exists
    overloadex = :(((p::$name))(::Type{Val{:param_Jac}},t,u,J,params) = $param_Jex)
    @eval $overloadex
  end

  return f
end

function build_component_funcs(symex)
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
  funcs
end

function symbolize(syms,param_dict_keys)
  symstr = symarr_to_symengine(syms)
  paramstr = symarr_to_symengine(param_dict_keys)
  full_symstr = symarr_to_symengine([syms;param_dict_keys])
  symdefineex = Expr(:(=),parse("("*full_symstr*")"),SymEngine.symbols(full_symstr))
  symtup = parse("("*symstr*")")
  @eval $symdefineex
  symtup = @eval $symtup # symtup is the tuple of SymEngine symbols for independent variables
  if length(param_dict_keys) == 1
    paramtup = parse("("*paramstr*",)")
  else
    paramtup = parse("("*paramstr*")")
  end
  paramtup = @eval $paramtup
  symtup,paramtup
end

function ode_findreplace(ex,symex,indvar_dict,param_dict,inline_dict;params_from_function=true)
  for (i,arg) in enumerate(ex.args)
    if isa(arg,Expr)
      ode_findreplace(arg,symex.args[i],indvar_dict,param_dict,inline_dict;params_from_function=params_from_function)
    elseif isa(arg,Symbol)
      s = string(arg)
      if haskey(indvar_dict,arg)
        ex.args[i] = :(u[$(indvar_dict[arg])]) # replace with u[i]
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
        ex.args[i] = :(du[$(indvar_dict[tmp])])
        symex.args[i] = :(du[$(indvar_dict[tmp])]) #also do in symbolic
      end
    end
  end
end

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

function build_p_funcs(paramfuncs,paramtup,indvar_dict,param_dict,inline_dict)
  pfuncs = Vector{Expr}(length(paramtup))
  param_dict_type = typeof(param_dict)
  for i in 1:length(paramtup)
    pfunc = :()
    param = Symbol(paramtup[i])
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
  :(1*$ex) # Add the 1 to make it an expression not a Symbol
end

function maketype(name,param_dict,origex,funcs,syms,fex;
                  pex=:(),
                  symjac=Matrix{SymEngine.Basic}(0,0),
                  Jex=:(),invjac=Matrix{SymEngine.Basic}(0,0),
                  invJex=:(),
                  symhes = Matrix{SymEngine.Basic}(0,0),
                  invhes = Matrix{SymEngine.Basic}(0,0),
                  Hex = :(),
                  invHex = :(),
                  params = Symbol[],
                  pfuncs=Vector{Expr}(0),
                  d_pfuncs = Vector{Expr}(0),
                  param_symjac=Matrix{SymEngine.Basic}(0,0),
                  param_Jex=:())

    @eval type $name <: ParameterizedFunction
        origex::Expr
        funcs::Vector{Expr}
        pfuncs::Vector{Expr}
        d_pfuncs::Vector{Expr}
        syms::Vector{Symbol}
        symjac::Matrix{SymEngine.Basic}
        invjac::Matrix{SymEngine.Basic}
        symhes::Matrix{SymEngine.Basic}
        invhes::Matrix{SymEngine.Basic}
        param_symjac::Matrix{SymEngine.Basic}
        Jex::Expr
        param_Jex::Expr
        invJex::Expr
        Hex::Expr
        invHex::Expr
        fex::Expr
        pex::Expr
        params::Vector{Symbol}
        $((:($x::$(typeof(t))) for (x, t) in param_dict)...)
    end

    # Export the type
    exportex = :(export $name)
    @eval $exportex

    # Make the default constructor
    new_ex = Meta.quot(origex)
    Jex_ex = Meta.quot(Jex)
    invJex_ex = Meta.quot(invJex)
    Hex_ex = Meta.quot(Hex)
    invHex_ex = Meta.quot(invHex)
    fex_ex = Meta.quot(fex)
    pex_ex = Meta.quot(pex)
    param_Jex_ex = Meta.quot(param_Jex)
    constructorex = :($(name)(;$(Expr(:kw,:origex,new_ex)),
                  $(Expr(:kw,:funcs,funcs)),
                  $(Expr(:kw,:pfuncs,pfuncs)),
                  $(Expr(:kw,:d_pfuncs,d_pfuncs)),
                  $(Expr(:kw,:syms,syms)),
                  $(Expr(:kw,:symjac,symjac)),
                  $(Expr(:kw,:invjac,invjac)),
                  $(Expr(:kw,:symhes,symhes)),
                  $(Expr(:kw,:invhes,invhes)),
                  $(Expr(:kw,:param_symjac,param_symjac)),
                  $(Expr(:kw,:Jex,Jex_ex)),
                  $(Expr(:kw,:param_Jex,param_Jex_ex)),
                  $(Expr(:kw,:invJex,invJex_ex)),
                  $(Expr(:kw,:Hex,Hex_ex)),
                  $(Expr(:kw,:invHex,invHex_ex)),
                  $(Expr(:kw,:fex,fex_ex)),
                  $(Expr(:kw,:pex,pex_ex)),
                  $(Expr(:kw,:params,params)),
                  $((Expr(:kw,x,t) for (x, t) in param_dict)...)) =
                  $(name)(origex,funcs,pfuncs,d_pfuncs,syms,
                  symjac,invjac,symhes,invhes,param_symjac,
                  Jex,param_Jex,invJex,Hex,invHex,fex,pex,
                  params,
                  $(((x for x in keys(param_dict))...))))
    eval(constructorex)

    # Make the type instance using the default constructor
    eval(name)()
end

function build_indvar_dict(ex)
  indvar_dict = OrderedDict{Symbol,Int}()
  for i in 2:2:length(ex.args) #Every odd line is line number
    arg = ex.args[i].args[1] #Get the first thing, should be dsomething
    nodarg = Symbol(string(arg)[2:end]) #Take off the d
    if !haskey(indvar_dict,nodarg)
      s = string(arg)
      indvar_dict[Symbol(string(arg)[2:end])] = i/2 # and label it the next int if not seen before
    end
  end
  syms = indvar_dict.keys
  indvar_dict,syms
end

function build_paramdicts(params)
  param_dict = OrderedDict{Symbol,Any}(); inline_dict = OrderedDict{Symbol,Any}()
  for i in 1:length(params)
    if params[i].head == :(=>)
      param_dict[params[i].args[1]] = params[i].args[2] # works for k=3, or k=>3
    elseif params[i].head == :(=)
      inline_dict[params[i].args[1]] = params[i].args[2] # works for k=3, or k=>3
    end
  end
  param_dict,inline_dict
end

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
  # Build the type
  f = maketype(name,param_dict,origex,funcs,syms,fex)
  # Overload the Call
  newsig = :($(sig.args...))
  overloadex = :(((p::$name))($(sig.args...)) = $ex)
  @eval $overloadex
  return f
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


### Utility Functions

"""
symarr_to_symengine(symarr::Vector{Symbol})

Converts a Vector{Symbol} into a string for SymEngine parsing

Symbol[:x,:y] --> "x,y"
"""
function symarr_to_symengine(symarr::Vector{Symbol})
  str = ""
  for sym in symarr
    str = str*string(sym)*","
  end
  str[1:end-1]
end

const FEM_SYMBOL_DICT = Dict{Symbol,Expr}(:x=>:(x[:,1]),:y=>:(x[:,2]),:z=>:(x[:,3]))

macro ode_def(name,ex,params...)
  name_ex = Meta.quot(name)
  ex_ex = Meta.quot(ex)
  params = Meta.quot(params)
  quote
    opts = Dict{Symbol,Bool}(
    :build_Jac => true,
    :build_InvJac => true,
    :build_Hes => true,
    :build_InvHes => true,
    :build_dpfuncs => true)
    ode_def_opts($(esc(name_ex)),opts,$(esc(ex_ex)),$(esc(params))...)
  end
end

macro ode_def_bare(name,ex,params...)
  name_ex = Meta.quot(name)
  ex_ex = Meta.quot(ex)
  params = Meta.quot(params)
  quote
    opts = Dict{Symbol,Bool}(
    :build_Jac => false,
    :build_InvJac => false,
    :build_Hes => false,
    :build_InvHes => false,
    :build_dpfuncs => false)
    ode_def_opts($(esc(name_ex)),opts,$(esc(ex_ex)),$(esc(params))...)
  end
end

macro ode_def_nohes(name,ex,params...)
  name_ex = Meta.quot(name)
  ex_ex = Meta.quot(ex)
  params = Meta.quot(params)
  quote
    opts = Dict{Symbol,Bool}(
    :build_Jac => true,
    :build_InvJac => true,
    :build_Hes => false,
    :build_InvHes => false,
    :build_dpfuncs => true)
    ode_def_opts($(esc(name_ex)),opts,$(esc(ex_ex)),$(esc(params))...)
  end
end

macro ode_def_noinvhes(name,ex,params...)
  name_ex = Meta.quot(name)
  ex_ex = Meta.quot(ex)
  params = Meta.quot(params)
  quote
    opts = Dict{Symbol,Bool}(
    :build_Jac => true,
    :build_InvJac => true,
    :build_Hes => true,
    :build_InvHes => false,
    :build_dpfuncs => true)
    ode_def_opts($(esc(name_ex)),opts,$(esc(ex_ex)),$(esc(params))...)
  end
end

macro ode_def_noinvjac(name,ex,params...)
  name_ex = Meta.quot(name)
  ex_ex = Meta.quot(ex)
  params = Meta.quot(params)
  quote
    opts = Dict{Symbol,Bool}(
    :build_Jac => true,
    :build_InvJac => false,
    :build_Hes => false,
    :build_InvHes => false,
    :build_dpfuncs => true)
    ode_def_opts($(esc(name_ex)),opts,$(esc(ex_ex)),$(esc(params))...)
  end
end

export ParameterizedFunction, @ode_def, @fem_def, ode_def_opts,
       @ode_def_bare, @ode_def_nohes, @ode_def_noinvjac, @ode_def_noinvhes

export jac_exists, invjac_exists, hes_exists, invhes_exists,
      paramjac_exists, pfunc_exists, pderiv_exists

end # module




##### Extra

# Jacobian Factorization

#=
local fsymjac_L
local fsymjac_U
local Jex_L
local Jex_U
try
  # Factorize the Jacobian
  fsymjac = lufact(symjac)
  fsymjac_L = fsymjac[:L]
  fsymjac_U = fsymjac[:U]
  Jex_L,Jex_U = build_fjac_func(fsymjac_L,fsymjac_U,indvar_dict,param_dict,inline_dict)
catch
  fsymjac_L = Matrix{SymEngine.Basic}(0,0)
  fsymjac_U = Matrix{SymEngine.Basic}(0,0)
end

"""
Builds the LU-factorized Jacobian functions
"""
function build_fjac_func(symjac_L,symjac_U,indvar_dict,param_dict,inline_dict)
  # Lower Triangle
  Jex = :()
  for i in 1:size(symjac_L,1)
    for j in 1:i
      ex = parse(string(symjac_L[i,j]))
      if typeof(ex) <: Expr
        ode_findreplace(ex,ex,indvar_dict,param_dict,inline_dict)
      else
        ex = ode_symbol_findreplace(ex,indvar_dict,param_dict,inline_dict)
      end
      push!(Jex.args,:(J[$i,$j] = $ex))
    end
  end
  Jex.head = :block
  push!(Jex.args,nothing)
  Jex_L = :(jac = (t,u,J)->$Jex)

  # Upper Triangle
  Jex = :()
  for j in 1:size(symjac_U,2)
    for i in 1:j
      ex = parse(string(symjac_U[i,j]))
      if typeof(ex) <: Expr
        ode_findreplace(ex,ex,indvar_dict,param_dict,inline_dict)
      else
        ex = ode_symbol_findreplace(ex,indvar_dict,param_dict,inline_dict)
      end
      push!(Jex.args,:(J[$i,$j] = $ex))
    end
  end
  Jex.head = :block
  push!(Jex.args,nothing)
  Jex_U = :(jac = (t,u,J)->$Jex)

  Jex_L,Jex_U
end

=#
