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
)
