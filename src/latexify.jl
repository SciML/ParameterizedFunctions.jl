@latexrecipe function f(func::SciMLBase.AbstractParameterizedFunction)
    return latexify(func.sys)
end
