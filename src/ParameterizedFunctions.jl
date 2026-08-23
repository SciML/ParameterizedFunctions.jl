"""
$(DocStringExtensions.README)
"""
module ParameterizedFunctions

using DocStringExtensions: DocStringExtensions
using DataStructures: DataStructures, OrderedDict
using DiffEqBase: DiffEqBase
using ModelingToolkit: ModelingToolkit, System, tosymbol
using ModelingToolkitBase: ModelingToolkitBase, @parameters
using PrecompileTools: @compile_workload, @setup_workload
using Symbolics: Symbolics, @variables
using SymbolicUtils: SymbolicUtils, BasicSymbolic

import LinearAlgebra
import SciMLBase

# ---------------------------------------------------------------------------
# Reexported interface (see the second `export` block below).
#
# `@ode_def` is a front end for ModelingToolkit: it builds a `System`, keeps it in
# the generated function's `sys` field, and returns something whose documented use
# is `ODEProblem(f, u0, tspan, p); solve(prob, alg)`. Those names reached the user
# through `@reexport using ModelingToolkit`, so `using ParameterizedFunctions` was
# enough to write the documented workflow and to inspect `f.sys`.
#
# The blanket reexport (520 names) is gone; this is the subset that normal
# documented use needs. Each name is imported from the module that owns it and
# stays documented there.
# ---------------------------------------------------------------------------
using CommonSolve: solve
using SciMLBase: ODEFunction, ODEProblem
using ModelingToolkitBase: equations, independent_variable, mtkcompile, observed,
    parameters, unknowns
using Symbolics: Differential, Num

include("ode_def_opts.jl")
include("utils.jl")
include("dict_build.jl")
include("macros.jl")

export @ode_def, ode_def_opts, @ode_def_bare, @ode_def_all

# Reexported ModelingToolkit interface; approved via `reexports_allow` in
# test/qa/qa.jl and documented in docs/src/api.md.
export ODEFunction, ODEProblem, solve
export System, equations, independent_variable, mtkcompile, observed, parameters,
    unknowns
export @parameters, @variables, Differential, Num

@setup_workload begin
    @compile_workload begin
        f = @ode_def begin
            dx = a * x - b * x * y
            dy = -c * y + d * x * y
        end a b c d
        u = [2.0, 3.0]
        p = [1.5, 1.0, 3.0, 1.0]
        f(u, p, 1.0)
        f.jac(u, p, 1.0)
        f.tgrad(u, p, 1.0)
    end
end

end # module
