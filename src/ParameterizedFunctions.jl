"""
$(DocStringExtensions.README)
"""
module ParameterizedFunctions

    using DocStringExtensions: DocStringExtensions
    using DataStructures: DataStructures, OrderedDict
    using DiffEqBase: DiffEqBase
    using Latexify: Latexify, @latexrecipe, latexify
    using Reexport: @reexport
    @reexport using ModelingToolkit
    using ModelingToolkit: ModelingToolkit, System, tosymbol
    using ModelingToolkitBase: @parameters
    using Symbolics: Symbolics, @variables
    using SymbolicUtils: SymbolicUtils, BasicSymbolic

    import LinearAlgebra
    import SciMLBase

    include("ode_def_opts.jl")
    include("utils.jl")
    include("dict_build.jl")
    include("macros.jl")
    include("latexify.jl")

    export @ode_def, ode_def_opts, @ode_def_bare, @ode_def_all
end # module
