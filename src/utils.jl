### Utility Functions

function build_component_funcs(symex)
  funcs = Vector{Expr}(undef,0) # Get all of the functions for symbolic computation
  for (i,arg) in enumerate(symex.args)
    if i%2 == 0
      ex = arg.args[2]
      if (typeof(ex) <: Symbol) || (typeof(ex) <: Number)
        push!(funcs,:($ex*1))
      else # It's an expression, just push
        fix_ex = copy(ex)
        flip_mult!(fix_ex)
        push!(funcs,fix_ex)
      end
    end
  end
  funcs
end

function modelingtoolkitize_expr(ex::Expr,vars,curmod)
    names = [x.op.name for x in vars]
    ex.head === :if && (ex = Expr(:call, ifelse, ex.args...))
    ex.head === :call || throw(ArgumentError("internal representation does not support non-call Expr"))
    op = ex.args[1] ∈ names ? vars[findfirst(x->ex.args[1] == x.op.name,vars)] : getproperty(curmod,ex.args[1]) # HACK
    args = Expression[modelingtoolkitize_expr(x,vars,curmod) for x in ex.args[2:end]]
    return Operation(op, args)
end

function modelingtoolkitize_expr(ex::Variable,vars,curmod)
    convert(Operation,ex)
end

function modelingtoolkitize_expr(ex::Symbol,vars,curmod)
    names = [x.op.name for x in vars]
    op = ex ∈ names ? vars[findfirst(x->ex == x.op.name,vars)] : getproperty(curmod,ex) # HACK
    convert(Expression,op)
end

function modelingtoolkitize_expr(ex::Number,vars,curmod)
    convert(Expression,ex)
end
