function maketype(name,param_dict,origex,funcs,syms,fex;
                  pex=:(), vector_ex = :(), vector_ex_return = :(),
                  symtgrad=Vector{SymEngine.Basic}(0),
                  tgradex = :(),
                  symjac=Matrix{SymEngine.Basic}(0,0),
                  Jex=:(),
                  expjac=Matrix{SymEngine.Basic}(0,0),
                  expJex=:(),
                  invjac=Matrix{SymEngine.Basic}(0,0),
                  invJex=:(),invWex=:(),invWex_t=:(),
                  symfuncs = Vector{SymEngine.Basic}(0),
                  symhes   = Matrix{SymEngine.Basic}(0,0),
                  invhes   = Matrix{SymEngine.Basic}(0,0),
                  syminvW  = Matrix{SymEngine.Basic}(0,0),
                  syminvW_t= Matrix{SymEngine.Basic}(0,0),
                  Hex = :(),
                  invHex = :(),
                  params = Symbol[],
                  pfuncs=Vector{Expr}(0),
                  d_pfuncs = Vector{Expr}(0),
                  param_symjac=Matrix{SymEngine.Basic}(0,0),
                  param_Jex=:())

    typeex = :(mutable struct $name <: DiffEqBase.AbstractParameterizedFunction{true}
        origex::Expr
        funcs::Vector{Expr}
        symfuncs::Vector{SymEngine.Basic}
        pfuncs::Vector{Expr}
        d_pfuncs::Vector{Expr}
        syms::Vector{Symbol}
        symtgrad::Vector{SymEngine.Basic}
        symjac::Matrix{SymEngine.Basic}
        expjac::Matrix{SymEngine.Basic}
        invjac::Matrix{SymEngine.Basic}
        syminvW::Matrix{SymEngine.Basic}
        syminvW_t::Matrix{SymEngine.Basic}
        symhes::Matrix{SymEngine.Basic}
        invhes::Matrix{SymEngine.Basic}
        param_symjac::Matrix{SymEngine.Basic}
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
                  $(Expr(:kw,:symfuncs,symfuncs)),
                  $(Expr(:kw,:pfuncs,pfuncs)),
                  $(Expr(:kw,:d_pfuncs,d_pfuncs)),
                  $(Expr(:kw,:syms,syms)),
                  $(Expr(:kw,:symtgrad,symtgrad)),
                  $(Expr(:kw,:symjac,symjac)),
                  $(Expr(:kw,:expjac,expjac)),
                  $(Expr(:kw,:invjac,invjac)),
                  $(Expr(:kw,:syminvW,syminvW)),
                  $(Expr(:kw,:syminvW_t,syminvW_t)),
                  $(Expr(:kw,:symhes,symhes)),
                  $(Expr(:kw,:invhes,invhes)),
                  $(Expr(:kw,:param_symjac,param_symjac)),
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
                  $(name)(origex,funcs,symfuncs,pfuncs,d_pfuncs,syms,
                  symtgrad,symjac,expjac,invjac,syminvW,syminvW_t,symhes,invhes,
                  param_symjac,tgradex,Jex,expJex,param_Jex,invJex,invWex,invWex_t,
                  Hex,invHex,fex,pex,vector_ex,vector_ex_return,params)) |> esc

    # Make the type instance using the default constructor
    typeex,constructorex
end
