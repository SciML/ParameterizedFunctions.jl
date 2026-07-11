using SciMLTesting, ParameterizedFunctions, Test

dependency_reexports(pkg) = Tuple(
    name for name in public_api_names(pkg)
        if isdefined(pkg, name) && parentmodule(getfield(pkg, name)) !== pkg
)

const DEPENDENCY_REEXPORTS = dependency_reexports(ParameterizedFunctions)

run_qa(
    ParameterizedFunctions;
    explicit_imports = true,
    api_docs_kwargs = (;
        rendered = true,
        ignore = DEPENDENCY_REEXPORTS,
        rendered_ignore = DEPENDENCY_REEXPORTS,
    ),
    aqua_kwargs = (;
        piracies = (;
            treat_as_own = [ParameterizedFunctions.SciMLBase.AbstractParameterizedFunction],
        ),
    ),
    # MTK v11 has undefined exports (Variable, find_solvables!) that get re-exported.
    # This is an upstream bug in ModelingToolkit, not ParameterizedFunctions.
    aqua_broken = (:undefined_exports,),
)
