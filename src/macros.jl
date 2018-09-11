macro ode_def(name,ex,params...)
    opts = Dict{Symbol,Bool}(
    :build_tgrad => true,
    :build_jac => true,
    :build_expjac => false,
    :build_invjac => true,
    :build_invW => true,
    :build_hes => false,
    :build_invhes => false,
    :build_dpfuncs => true)
    name isa Expr ? ode_def_opts(gensym(),opts,name,params...) :
        ode_def_opts(name,opts,ex,params...)
end

macro ode_def_bare(name,ex,params...)
    opts = Dict{Symbol,Bool}(
    :build_tgrad => false,
    :build_jac => false,
    :build_expjac => false,
    :build_invjac => false,
    :build_invW => false,
    :build_hes => false,
    :build_invhes => false,
    :build_dpfuncs => false)
    name isa Expr ? ode_def_opts(gensym(),opts,name,params...) :
        ode_def_opts(name,opts,ex,params...)
end

macro ode_def_nohes(name,ex,params...)
  opts = Dict{Symbol,Bool}(
  :build_tgrad => true,
  :build_jac => true,
  :build_expjac => false,
  :build_invjac => true,
  :build_invW => true,
  :build_hes => false,
  :build_invhes => false,
  :build_dpfuncs => true)
  name isa Expr ? ode_def_opts(gensym(),opts,name,params...) :
    ode_def_opts(name,opts,ex,params...)
end

macro ode_def_noinvhes(name,ex,params...)
    opts = Dict{Symbol,Bool}(
    :build_tgrad => true,
    :build_jac => true,
    :build_expjac => false,
    :build_invjac => true,
    :build_invW => true,
    :build_hes => false,
    :build_invhes => false,
    :build_dpfuncs => true)
    name isa Expr ? ode_def_opts(gensym(),opts,name,params...) :
        ode_def_opts(name,opts,ex,params...)
end

macro ode_def_noinvjac(name,ex,params...)
    opts = Dict{Symbol,Bool}(
    :build_tgrad => true,
    :build_jac => true,
    :build_expjac => false,
    :build_invjac => false,
    :build_invW => false,
    :build_hes => false,
    :build_invhes => false,
    :build_dpfuncs => true)
    name isa Expr ? ode_def_opts(gensym(),opts,name,params...) :
        ode_def_opts(name,opts,ex,params...)
end
