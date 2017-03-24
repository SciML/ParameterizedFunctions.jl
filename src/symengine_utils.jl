
function symbolize(syms,param_dict_keys)
  symstr = symarr_to_symengine(syms)
  paramstr = symarr_to_symengine(param_dict_keys)
  full_symstr = symarr_to_symengine([syms;param_dict_keys])
  symdefineex = Expr(:(=),parse("("*full_symstr*")"),SymEngine.symbols(full_symstr))
  symtup = parse("("*symstr*")")
  @eval $symdefineex
  symtup = @eval $symtup # symtup is the tuple of SymEngine symbols for independent variables
  if length(param_dict_keys) == 1
    paramtup = parse("("*paramstr*",)")
  else
    paramtup = parse("("*paramstr*")")
  end
  paramtup = @eval $paramtup
  symtup,paramtup
end

### Utility Functions

"""
symarr_to_symengine(symarr::Vector{Symbol})

Converts a Vector{Symbol} into a string for SymEngine parsing

Symbol[:x,:y] --> "x,y"
"""
function symarr_to_symengine(symarr::Vector{Symbol})
  str = ""
  for sym in symarr
    str = str*string(sym)*","
  end
  str[1:end-1]
end
