function maketype(name,param_dict,origex,funcs,syms,fex;
                  pex=:(), vector_ex = :(), vector_ex_return = :(),
                  tgradex = :(),
                  Jex=:(),
                  expJex=:(),
                  invJex=:(),invWex=:(),invWex_t=:(),
                  Hex = :(),
                  invHex = :(),
                  params = Symbol[],
                  pfuncs=Vector{Expr}(0),
                  d_pfuncs = Vector{Expr}(0),
                  param_Jex=:())

    typeex = :(mutable struct $name <: DiffEqBase.AbstractParameterizedFunction{true}
        origex::Expr
        funcs::Vector{Expr}
        pfuncs::Vector{Expr}
        d_pfuncs::Vector{Expr}
        syms::Vector{Symbol}
        tgradex::Expr
        Jex::Expr
        expJex::Expr
        param_Jex::Expr
        invJex::Expr
        invWex::Expr
        invWex_t::Expr
        Hex::Expr
        invHex::Expr
        fex::Expr
        pex::Expr
        vector_ex::Expr
        vector_ex_return::Expr
        params::Vector{Symbol}
    end)

    # Make the default constructor
    new_ex = Meta.quot(origex)
    tgradex_ex = Meta.quot(tgradex)
    Jex_ex = Meta.quot(Jex)
    expJex_ex = Meta.quot(expJex)
    invJex_ex = Meta.quot(invJex)
    invWex_ex = Meta.quot(invWex)
    invWex_t_ex = Meta.quot(invWex_t)
    Hex_ex = Meta.quot(Hex)
    invHex_ex = Meta.quot(invHex)
    fex_ex = Meta.quot(fex)
    pex_ex = Meta.quot(pex)
    vector_ex_ex = Meta.quot(vector_ex)
    vector_ex_return_ex = Meta.quot(vector_ex_return)
    param_Jex_ex = Meta.quot(param_Jex)
    constructorex = :($(name)(;$(Expr(:kw,:origex,new_ex)),
                  $(Expr(:kw,:funcs,funcs)),
                  $(Expr(:kw,:pfuncs,pfuncs)),
                  $(Expr(:kw,:d_pfuncs,d_pfuncs)),
                  $(Expr(:kw,:syms,syms)),
                  $(Expr(:kw,:tgradex,tgradex_ex)),
                  $(Expr(:kw,:Jex,Jex_ex)),
                  $(Expr(:kw,:expJex,expJex_ex)),
                  $(Expr(:kw,:param_Jex,param_Jex_ex)),
                  $(Expr(:kw,:invJex,invJex_ex)),
                  $(Expr(:kw,:invWex,invWex_ex)),
                  $(Expr(:kw,:invWex_t,invWex_t_ex)),
                  $(Expr(:kw,:Hex,Hex_ex)),
                  $(Expr(:kw,:invHex,invHex_ex)),
                  $(Expr(:kw,:fex,fex_ex)),
                  $(Expr(:kw,:pex,pex_ex)),
                  $(Expr(:kw,:vector_ex,vector_ex_ex)),
                  $(Expr(:kw,:vector_ex_return,vector_ex_return_ex)),
                  $(Expr(:kw,:params,params))) =
                  $(name)(origex,funcs,pfuncs,d_pfuncs,syms,
                  tgradex,Jex,expJex,param_Jex,
                  invJex,invWex,invWex_t,
                  Hex,invHex,fex,pex,vector_ex,vector_ex_return,params)) |> esc

    # Make the type instance using the default constructor
    typeex,constructorex
end
