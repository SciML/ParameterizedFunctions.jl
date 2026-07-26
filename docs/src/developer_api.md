# Developer API

`ode_def_opts` is the low-level implementation used by the public `@ode_def`,
`@ode_def_bare`, and `@ode_def_all` macros. It is exported for compatibility,
but package users should use those macros. It is intended for SciML package
developers that need to define a macro with a different set of generated
derivative functions; its input expression must follow the assignment-based
ODE DSL accepted by `@ode_def`.

```@docs
ParameterizedFunctions.ode_def_opts
```
