# The ode_def macro

```@docs
ParameterizedFunctions.@ode_def
ParameterizedFunctions.@ode_def_bare
ParameterizedFunctions.@ode_def_all
```

```@example ode_def
using ParameterizedFunctions

linear = @ode_def begin
    dx = a * x
end a

linear([2.0], [3.0], 0.0)
```
