using SciMLTesting, ParameterizedFunctions, Test

run_qa(
    ParameterizedFunctions;
    explicit_imports = true,
    aqua_kwargs = (;
        piracies = (;
            treat_as_own = [ParameterizedFunctions.SciMLBase.AbstractParameterizedFunction],
        ),
    ),
    # MTK v11 has undefined exports (Variable, find_solvables!) that get re-exported.
    # This is an upstream bug in ModelingToolkit, not ParameterizedFunctions.
    aqua_broken = (:undefined_exports,),
    ei_kwargs = (;
        # AbstractParameterizedFunction is owned by SciMLBase and re-exported by
        # DiffEqBase; the @latexrecipe dispatch accesses it via DiffEqBase (a direct
        # dep). Both names are non-public in their respective modules.
        all_qualified_accesses_via_owners = (; ignore = (:AbstractParameterizedFunction,)),
        all_qualified_accesses_are_public = (;
            # AbstractParameterizedFunction: non-public in SciMLBase/DiffEqBase
            ignore = (:AbstractParameterizedFunction,),
        ),
    ),
)
