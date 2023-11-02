__precompile__()

"""
$(DocStringExtensions.README)
"""
module ParameterizedFunctions

using DocStringExtensions
using DataStructures, DiffEqBase, Latexify
using Reexport
@reexport using ModelingToolkit
using ModelingToolkit: Sym, FnType, tosymbol

import LinearAlgebra

import Base: getindex

include("ode_def_opts.jl")
include("utils.jl")
include("dict_build.jl")
include("macros.jl")
include("latexify.jl")

export @ode_def, ode_def_opts, @ode_def_bare, @ode_def_all
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
      if ex isa Expr
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
      if ex isa Expr
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
