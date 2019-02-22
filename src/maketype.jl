function maketype(name,param_dict,origex,funcs,syms,fex;
                  pex=:(), vector_ex = :(), vector_ex_return = :(),
                  tgradex = :(),
                  Jex=:(),
                  expJex=:(),
                  invJex=:(),invWex=:(),invWex_t=:(),
                  Hex = :(),
                  invHex = :(),
                  params = Symbol[],
                  pfuncs=Vector{Expr}(undef,0),
                  d_pfuncs = Vector{Expr}(undef,0),
                  symjac = Matrix{Expr}(undef,0,0),
                  symtgrad = Vector{Expr}(undef,0),
                  param_Jex=:(),
                  f_expr=:(),
                  tgrad_expr=:(),
                  jac_expr=:(),
                  invjac_expr=:(),
                  invW_expr=:(),
                  invW_t_expr=:(),
                  param_jac_expr=:())

    typeex = :(mutable struct $name{F,J,T,W,Wt,PJ,TT1,TT2} <: DiffEqBase.AbstractParameterizedFunction{true}
        f::F
        analytic::Nothing
        jac::J
        tgrad::T
        invW::W
        invW_t::Wt
        paramjac::PJ
        mass_matrix::LinearAlgebra.UniformScaling{Bool}
        jac_prototype::Nothing
        origex::Expr
        funcs::Vector{Expr}
        pfuncs::Vector{Expr}
        d_pfuncs::Vector{Expr}
        syms::Vector{Symbol}
        symjac::Matrix{TT1} # https://github.com/symengine/SymEngine.jl/issues/122
        symtgrad::Vector{TT2}
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
        nf::Ref{Int} # number of function evals
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

    constructorex = :($(name)() =
                  $(name)($f_expr,nothing,
                  $jac_expr,$tgrad_expr,$invW_expr,$invW_t_expr,$param_jac_expr,
                  $(LinearAlgebra.I),nothing,$new_ex,$funcs,$pfuncs,$d_pfuncs,
                  $syms,$symjac,$symtgrad,
                  $tgradex_ex,$Jex_ex,$expJex_ex,$param_Jex_ex,
                  $invJex_ex,$invWex_ex,$invWex_t_ex,
                  $Hex_ex,$invHex_ex,$fex_ex,$pex_ex,$vector_ex_ex,
                  $vector_ex_return_ex,$params,Ref(0))) |> esc

    callex = :(((f::$name))(args...) = (f.nf[] += 1; f.f(args...))) |> esc
    callex2 = :(((f::$name))(u,p,t::Number) = (f.nf[] += 1;du=similar(u);f.f(du,u,p,t);du)) |> esc

    # Make the type instance using the default constructor
    typeex,constructorex,callex,callex2
end
