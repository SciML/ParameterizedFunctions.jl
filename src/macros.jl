macro ode_def(name,ex,params...)
  name_ex = Meta.quot(name)
  ex_ex = Meta.quot(ex)
  params = Meta.quot(params)
  quote
    opts = Dict{Symbol,Bool}(
    :build_tgrad => true,
    :build_jac => true,
    :build_expjac => false,
    :build_invjac => true,
    :build_invW => true,
    :build_hes => true,
    :build_invhes => true,
    :build_dpfuncs => true)
    ode_def_opts($(esc(name_ex)),opts,$(esc(ex_ex)),$(esc(params))...)
  end
end

macro ode_def_mm(name,ex,M,params...)
  name_ex = Meta.quot(name)
  ex_ex = Meta.quot(ex)
  params = Meta.quot(params)
  quote
    opts = Dict{Symbol,Bool}(
    :build_tgrad => true,
    :build_jac => true,
    :build_expjac => false,
    :build_invjac => true,
    :build_invW => true,
    :build_hes => true,
    :build_invhes => true,
    :build_dpfuncs => true)
    ode_def_opts($(esc(name_ex)),opts,$(esc(ex_ex)),$(esc(params))...;M=$(esc(M)))
  end
end

macro ode_def_bare(name,ex,params...)
  name_ex = Meta.quot(name)
  ex_ex = Meta.quot(ex)
  params = Meta.quot(params)
  quote
    opts = Dict{Symbol,Bool}(
    :build_tgrad => false,
    :build_jac => false,
    :build_expjac => false,
    :build_invjac => false,
    :build_invW => false,
    :build_hes => false,
    :build_invhes => false,
    :build_dpfuncs => false)
    ode_def_opts($(esc(name_ex)),opts,$(esc(ex_ex)),$(esc(params))...)
  end
end

macro ode_def_nohes(name,ex,params...)
  name_ex = Meta.quot(name)
  ex_ex = Meta.quot(ex)
  params = Meta.quot(params)
  quote
    opts = Dict{Symbol,Bool}(
    :build_tgrad => true,
    :build_jac => true,
    :build_expjac => false,
    :build_invjac => true,
    :build_invW => true,
    :build_hes => false,
    :build_invhes => false,
    :build_dpfuncs => true)
    ode_def_opts($(esc(name_ex)),opts,$(esc(ex_ex)),$(esc(params))...)
  end
end

macro ode_def_nohes_mm(name,ex,M,params...)
  name_ex = Meta.quot(name)
  ex_ex = Meta.quot(ex)
  params = Meta.quot(params)
  quote
    opts = Dict{Symbol,Bool}(
    :build_tgrad => true,
    :build_jac => true,
    :build_expjac => false,
    :build_invjac => true,
    :build_invW => true,
    :build_hes => false,
    :build_invhes => false,
    :build_dpfuncs => true)
    ode_def_opts($(esc(name_ex)),opts,$(esc(ex_ex)),$(esc(params))...;M=$(esc(M)))
  end
end

macro ode_def_noinvhes(name,ex,params...)
  name_ex = Meta.quot(name)
  ex_ex = Meta.quot(ex)
  params = Meta.quot(params)
  quote
    opts = Dict{Symbol,Bool}(
    :build_tgrad => true,
    :build_jac => true,
    :build_expjac => false,
    :build_invjac => true,
    :build_invW => true,
    :build_hes => true,
    :build_invhes => false,
    :build_dpfuncs => true)
    ode_def_opts($(esc(name_ex)),opts,$(esc(ex_ex)),$(esc(params))...)
  end
end

macro ode_def_noinvhes_mm(name,M,ex,params...)
  name_ex = Meta.quot(name)
  ex_ex = Meta.quot(ex)
  params = Meta.quot(params)
  quote
    opts = Dict{Symbol,Bool}(
    :build_tgrad => true,
    :build_jac => true,
    :build_expjac => false,
    :build_invjac => true,
    :build_invW => true,
    :build_hes => true,
    :build_invhes => false,
    :build_dpfuncs => true)
    ode_def_opts($(esc(name_ex)),opts,$(esc(ex_ex)),$(esc(params))...;M=$(esc(M)))
  end
end

macro ode_def_noinvjac(name,ex,params...)
  name_ex = Meta.quot(name)
  ex_ex = Meta.quot(ex)
  params = Meta.quot(params)
  quote
    opts = Dict{Symbol,Bool}(
    :build_tgrad => true,
    :build_jac => true,
    :build_expjac => false,
    :build_invjac => false,
    :build_invW => false,
    :build_hes => false,
    :build_invhes => false,
    :build_dpfuncs => true)
    ode_def_opts($(esc(name_ex)),opts,$(esc(ex_ex)),$(esc(params))...)
  end
end
