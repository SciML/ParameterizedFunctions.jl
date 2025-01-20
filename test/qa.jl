using ParameterizedFunctions, Aqua
@testset "Aqua" begin
    Aqua.find_persistent_tasks_deps(ParameterizedFunctions)
    Aqua.test_ambiguities(ParameterizedFunctions, recursive = false)
    Aqua.test_deps_compat(ParameterizedFunctions)
    Aqua.test_piracies(ParameterizedFunctions,
        treat_as_own = [ParameterizedFunctions.SciMLBase.AbstractParameterizedFunction])
    Aqua.test_project_extras(ParameterizedFunctions)
    Aqua.test_stale_deps(ParameterizedFunctions)
    Aqua.test_unbound_args(ParameterizedFunctions)
    Aqua.test_undefined_exports(ParameterizedFunctions)
end
