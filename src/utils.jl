### Utility Functions

function build_component_funcs(symex)
    funcs = Vector{Expr}(undef, 0) # Get all of the functions for symbolic computation
    for (i, arg) in enumerate(symex.args)
        if i % 2 == 0
            ex = arg.args[2]
            if (typeof(ex) <: Symbol) || (typeof(ex) <: Number)
                push!(funcs, :($ex * 1))
            else # It's an expression, just push
                fix_ex = copy(ex)
                flip_mult!(fix_ex)
                push!(funcs, fix_ex)
            end
        end
    end
    return funcs
end

function modelingtoolkitize_expr(ex::Expr, vars, curmod)
    names = [tosymbol(x) for x in vars]
    ex.head === :if && (ex = Expr(:call, ifelse, ex.args...))
    ex.head === :call ||
        throw(ArgumentError("internal representation does not support non-call Expr"))
    op = ex.args[1] ∈ names ? vars[findfirst(x -> ex.args[1] == tosymbol(x), vars)] :
         getproperty(curmod, ex.args[1]) # HACK
    return op((modelingtoolkitize_expr(x, vars, curmod) for x in ex.args[2:end])...)
end

function modelingtoolkitize_expr(ex::Sym, vars, curmod)
    return ex
end

function modelingtoolkitize_expr(ex::Symbol, vars, curmod)
    names = tosymbol.(vars)
    idx = findfirst(x -> ex == x, names)
    op = idx !== nothing ? vars[idx] : getproperty(curmod, ex)
    return op
end

function modelingtoolkitize_expr(ex::Number, vars, curmod)
    return ex
end
