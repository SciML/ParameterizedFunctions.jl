function param_values(f::AbstractParameterizedFunction)
  [getfield(f,s) for s in f.params]
end
