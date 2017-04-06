function param_values(f::AbstractParameterizedFunction)
  [getfield(f,s) for s in f.params]
end

num_params(f::AbstractParameterizedFunction) = length(f.params)
