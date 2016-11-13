macro ode_def(name,ex,params...)
  name_ex = Meta.quot(name)
  ex_ex = Meta.quot(ex)
  params = Meta.quot(params)
  quote
    opts = Dict{Symbol,Bool}(
    :build_tgrad => true,
    :build_Jac => true,
    :build_InvJac => true,
    :build_InvW => true,
    :build_Hes => true,
    :build_InvHes => true,
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
    :build_Jac => true,
    :build_InvJac => true,
    :build_InvW => true,
    :build_Hes => true,
    :build_InvHes => true,
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
    :build_Jac => false,
    :build_InvJac => false,
    :build_InvW => false,
    :build_Hes => false,
    :build_InvHes => false,
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
    :build_Jac => true,
    :build_InvJac => true,
    :build_InvW => true,
    :build_Hes => false,
    :build_InvHes => false,
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
    :build_Jac => true,
    :build_InvJac => true,
    :build_InvW => true,
    :build_Hes => false,
    :build_InvHes => false,
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
    :build_Jac => true,
    :build_InvJac => true,
    :build_InvW => true,
    :build_Hes => true,
    :build_InvHes => false,
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
    :build_Jac => true,
    :build_InvJac => true,
    :build_InvW => true,
    :build_Hes => true,
    :build_InvHes => false,
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
    :build_Jac => true,
    :build_InvJac => false,
    :build_InvW => false,
    :build_Hes => false,
    :build_InvHes => false,
    :build_dpfuncs => true)
    ode_def_opts($(esc(name_ex)),opts,$(esc(ex_ex)),$(esc(params))...)
  end
end
