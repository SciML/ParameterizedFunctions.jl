@latexrecipe function f(func::DiffEqBase.AbstractParameterizedFunction)
    return latexify(func.sys)
end
