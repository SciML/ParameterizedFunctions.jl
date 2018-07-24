function ode_def_opts(name::Symbol,opts::Dict{Symbol,Bool},ex::Expr,params...;_M=I,depvar=:t)
  # depvar is the dependent variable. Defaults to t
  # M is the mass matrix in RosW, must be a constant!
  origex = copy(ex) # Save the original expression

  if !(eltype(params) <: Symbol)
    error("The syntax for ParameterizedFunctions has changed. Simply list the parameters at the end, i.e. `a b c d`, instead of `a=5.0 b=>3.0 ...`. Parameters are defined in the problem type. See the documentation for more information.")
  end
  params = Symbol[params...]

  ## Build independent variable dictionary
  indvar_dict,syms = build_indvar_dict(ex,depvar)
  ####
  # Build the Expressions

  # Run find replace to make the function expression
  symex = copy(ex) # Different expression for symbolic computations
  ode_findreplace(ex,symex,indvar_dict,params)
  push!(ex.args,nothing) # Make the return void
  fex = ex # Save this expression as the expression for the call

  # Parameter-Explicit Functions
  pex = copy(origex) # Build it from the original expression
  # Parameter find/replace
  ode_findreplace(pex,copy(ex),indvar_dict,params;params_from_function=false)

  # Vectorized Functions
  vector_ex = copy(origex) # Build it from the original expression
  ode_findreplace(vector_ex,copy(origex),indvar_dict,params;
                  vectorized_form=true,vectorized_returns=:slice)
  vector_ex_return = copy(origex) # Build it from the original expression
  ode_findreplace(vector_ex_return,copy(origex),indvar_dict,params;
                  vectorized_form=true,vectorized_returns=:vals)
  dus = [Symbol("internal_var___du$i") for i in 1:length(keys(indvar_dict))] # TODO: vectorized forms need @. to work
  push!(vector_ex_return.args,:(hcat($(dus...)))) # Make the return void

  ######
  # Build the Functions

  # Get the component functions
  funcs = build_component_funcs(symex)

  numsyms = length(indvar_dict)
  numparams = length(params)
  M = Matrix(1*_M,numsyms,numsyms)
  # Parameter Functions
  paramfuncs = Vector{Vector{Expr}}(undef, numparams)
  for i in 1:numparams
    tmp_pfunc = Vector{Expr}(undef,length(funcs))
    for j in eachindex(funcs)
      tmp_pfunc[j] = copy(funcs[j])
    end
    paramfuncs[i] = tmp_pfunc
  end
  pfuncs = build_p_funcs(paramfuncs,indvar_dict,params)

  # Symbolic Setup
  symfuncs = Vector{SymEngine.Basic}(undef, 0)
  symtgrad = Vector{SymEngine.Basic}(undef, 0)
  symjac   = Matrix{SymEngine.Basic}(undef, 0,0)
  expjac   = Matrix{SymEngine.Basic}(undef, 0,0)
  invjac   = Matrix{SymEngine.Basic}(undef, 0,0)
  symhes   = Matrix{SymEngine.Basic}(undef, 0,0)
  invhes   = Matrix{SymEngine.Basic}(undef, 0,0)
  syminvW  = Matrix{SymEngine.Basic}(undef, 0,0)
  syminvW_t= Matrix{SymEngine.Basic}(undef, 0,0)
  param_symjac = Matrix{SymEngine.Basic}(undef, 0,0)
  tgradex = :(error("t-gradient Does Not Exist"))
  tgrad_exists = false
  Jex = :(error("Jacobian Does Not Exist"))
  jac_exists = false
  expJex = :(error("Exponential Jacobian Does Not Exist"))
  expjac_exists = false
  invJex = :(error("Inverse Jacobian Does Not Exist"))
  invjac_exists = false
  invWex = :(error("Inverse Rosenbrock-W Does Not Exist"))
  invW_exists = false
  invWex_t = :(error("Inverse Rosenbrock-W Transformed Does Not Exist"))
  invW__t_exists = false
  Hex = :(error("Hessian Does Not Exist"))
  hes_exists = false
  invHex = :(error("Inverse Hessian Does Not Exist"))
  invhes_exists = false
  param_Jex = :(error("Parameter Jacobian Does Not Exist"))
  param_jac_exists = false

  d_pfuncs = Vector{Expr}(undef, 0)
  param_symjac = Matrix{SymEngine.Basic}(undef,numsyms,numparams)
  pderiv_exists = false

  if opts[:build_tgrad] || opts[:build_jac] || opts[:build_dpfuncs]
    try #do symbolic calculations

      # Set Internal γ, used as a symbol for letting users pass an extra scalar
      γ = symbols("internal_γ")

      # Build the symbolic functions

      symfuncs = [SymEngine.Basic(f) for f in funcs]

      if opts[:build_tgrad]
        try
          symtgrad = [diff(f,depvar) for f in symfuncs]
          tgrad_exists = true
          tgradex = build_tgrad_func(symtgrad,indvar_dict,params)
        catch err
          @warn("Time Derivative Gradient could not be built")
        end
      end

      if opts[:build_jac]
        try #Jacobians and Hessian
          # Build the Jacobian Matrix of SymEngine Expressions
          symjac = Matrix{SymEngine.Basic}(undef,numsyms,numsyms)
          for i in eachindex(funcs)
            for j in eachindex(syms)
              symjac[i,j] = diff(symfuncs[i],syms[j])
            end
          end

          # Build the Julia function
          Jex = build_jac_func(symjac,indvar_dict,params)
          bad_derivative(Jex)
          jac_exists = true

          if opts[:build_expjac]
            try
              expjac = exp(γ*symjac) # This does not work, which is why disabled
              expJex = build_jac_func(expjac,indvar_dict,params)
              bad_derivative(expJex)
              expjac_exists = true
            catch
              @warn("Jacobian could not exponentiate")
            end
          end

          if opts[:build_invjac]
            try # Jacobian Inverse
              L = Base.convert(SymEngine.CDenseMatrix,symjac)
              invjac = inv(L)
              invJex = build_jac_func(invjac,indvar_dict,params)
              bad_derivative(invJex)
              invjac_exists = true
            catch err
              @warn("Jacobian could not invert")
            end
          end
          if opts[:build_invW]
            try # Rosenbrock-W Inverse
              L = Base.convert(SymEngine.CDenseMatrix,M - γ*symjac)
              syminvW = inv(L)
              L = Base.convert(SymEngine.CDenseMatrix,M/γ - symjac)
              syminvW_t = inv(L)
              invWex = build_jac_func(syminvW,indvar_dict,params)
              bad_derivative(invWex)
              invW_exists = true
              invWex_t = build_jac_func(syminvW_t,indvar_dict,params)
              bad_derivative(invWex_t)
              invW_t_exists = true
            catch err
              @warn("Rosenbrock-W could not invert")
            end
          end
          if opts[:build_hes]
            try # Hessian
              symhes = Matrix{SymEngine.Basic}(numsyms,numsyms)
              for i in eachindex(funcs), j in eachindex(syms)
                symhes[i,j] = diff(symjac[i,j],syms[j])
              end
              # Build the Julia function
              Hex = build_jac_func(symhes,indvar_dict,params)
              bad_derivative(Hex)
              hes_exists = true
              if opts[:build_invhes]
                try # Hessian Inverse
                  L = Base.convert(SymEngine.CDenseMatrix,symhes)
                  invhes = inv(L)
                  invHex = build_jac_func(invhes,indvar_dict,params)
                  bad_derivative(invHex)
                  invhes_exists = true
                catch err
                  @warn("Hessian could not invert")
                end
              end
            catch err
              @warn("Hessians failed to build.")
            end
          end
        catch err
          @warn("Failed to build the Jacobian. This means the Hessian is not built as well.")
        end
      end # End Jacobian tree

      if opts[:build_dpfuncs]
        try # Parameter Gradients
          d_paramfuncs  = Vector{Vector{Expr}}(undef,numparams)
          for i in eachindex(params)
            tmp_dpfunc = Vector{Expr}(undef,length(funcs))
            for j in eachindex(funcs)
              funcex = funcs[j]
              d_curr = diff(symfuncs[j],params[i])
              param_symjac[j,i] = d_curr
              symfunc_str = Meta.parse(string(d_curr))
              if typeof(symfunc_str) <: Number
                tmp_dpfunc[j] = :($symfunc_str*1)
              elseif typeof(symfunc_str) <: Symbol
                tmp_dpfunc[j] = :($symfunc_str*1)
              else
                tmp_dpfunc[j] = symfunc_str
              end
            end
            d_paramfuncs[i] = tmp_dpfunc
          end
          d_pfuncs = build_p_funcs(d_paramfuncs,indvar_dict,params)
          pderiv_exists = true

          # Now build the parameter Jacobian
          param_symjac_ex = Matrix{Expr}(undef,numsyms,numparams)
          for i in 1:numparams
            param_symjac_ex[:,i] = d_paramfuncs[i]
          end

          param_Jex = build_jac_func(param_symjac_ex,indvar_dict,params,params_from_function=false)
          param_jac_exists = true
        catch err
          @warn("Failed to build the parameter derivatives.")
        end
      end
    catch err
      @warn("Symbolic calculations could not initiate. Likely there's a function which is not differentiable by SymEngine.")
    end
  end

  # Build the Function
  f_expr = :((internal_var___du,internal_var___u,internal_var___p,t::Number) -> $pex)

  # Add the mass matrix
  mm_expr = quote $M end

  # Add the t gradient
  if tgrad_exists
    tgrad_expr = :((internal_var___grad,internal_var___u,internal_var___p,t) -> $tgradex)
  else
    tgrad_expr = :(nothing)
  end

  # Add the Jacobian
  if jac_exists
    jac_expr = :((internal_var___J,internal_var___u,internal_var___p,t) -> $Jex)
  else
    jac_expr = :(nothing)
  end

  # Add the Inverse Jacobian
  if invjac_exists
    invjac_expr = :((internal_var___J,internal_var___u,internal_var___p,t) -> $invJex)
  else
    invjac_expr = :(nothing)
  end

  # Add the Inverse Rosenbrock-W
  if invW_exists
    invW_expr = :((internal_var___J,internal_var___u,internal_var___p,internal_γ,t) -> $invWex)
  else
    invW_expr = :(nothing)
  end

  # Add the Inverse Rosenbrock-W Transformed
  if invW_exists
    invW_t_expr = :((internal_var___J,internal_var___u,internal_var___p,internal_γ,t) -> $invWex_t)
  else
    invW_t_expr = :(nothing)
  end

  # Add Parameter Jacobian
  if param_jac_exists
    param_jac_expr = :((internal_var___J,internal_var___u,internal_var___p,t) -> $param_Jex)
  else
    param_jac_expr = :(nothing)
  end

  # Build the type
  exprs = Vector{Expr}(undef, 0)

  typeex,constructorex,callex,callex2 = maketype(name,params,origex,
               funcs,syms,fex,
               pex=pex,
               vector_ex = vector_ex,vector_ex_return = vector_ex_return,
               tgradex=tgradex,expJex=expJex,Jex=Jex,
               invWex=invWex,invWex_t=invWex_t,
               invJex=invJex,Hex=Hex,
               invHex=invHex,params=params,
               pfuncs=pfuncs,d_pfuncs=d_pfuncs,
               param_Jex=param_Jex,
               f_expr=f_expr,tgrad_expr=tgrad_expr,
               jac_expr=jac_expr,
               invjac_expr=invjac_expr,
               invW_expr=invW_expr,
               invW_t_expr=invW_t_expr,
               param_jac_expr=param_jac_expr,
               mm_expr=mm_expr,
               symjac = convert.(Expr,symjac),
               symtgrad = convert.(Expr,symtgrad))

  push!(exprs,typeex)
  push!(exprs,constructorex)
  push!(exprs,callex)
  push!(exprs,callex2)

  # Return the type from the default consturctor
  def_const_ex = :(($name)()) |> esc
  push!(exprs,def_const_ex)
  expr_arr_to_block(exprs)
end
