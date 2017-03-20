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
