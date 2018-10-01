using ParameterizedFunctions, DiffEqBase
using Test, LinearAlgebra, InteractiveUtils

using SpecialFunctions

### ODE Macros

println("Build some examples")
f_t = @ode_def SymCheck begin # Checks for error due to symbol on 1
  dx = x
  dy = -c*y + d*x*y*t^2
end a b c d

@test DiffEqBase.__has_syms(f_t)

f_t2 = @ode_def_noinvjac SymCheck2 begin # Checks for error due to symbol on 1
  dx = 1
  dy = -c*y + d*x*y*t^2
end a b c d

f_t3 = @ode_def_noinvjac ExprCheck begin # Checks for error due to symbol on 1
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a b c d # Change to π after unicode fix

f = @ode_def LotkaVolterra begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a b c d

f_2 = @ode_def_nohes LotkaVolterra3 begin
  dx = a*x - b^2*x*y
  dy = -c*y + d*x*y
end a b c d

println("Test Values")
t = 1.0
u = [2.0,3.0]
p = [1.5,1,3,1]
du = zeros(2)
grad = similar(du)
J = zeros(2,2)
iJ= zeros(2,2)
iW= zeros(2,2)
f(du,u,p,t)
@test du == [-3.0,-3.0]
@test du == f(u,p,t)
f_t(du,u,p,t)
@test du == [2.0,-3.0]
f_t2(du,u,p,t)
@test du == [1.0,-3.0]

println("Test t-gradient")
f.tgrad(grad,u,p,t)
@test grad == zeros(2)
f_t.tgrad(grad,u,p,t)
@test grad == [0.0;12.0]

println("Test Jacobians")
f.jac(J,u,p,t)
@test J  == [-1.5 -2.0
             3.0 -1.0]
#@test f.jac(u, p, t) == [-1.5 -2.0; 3.0 -1.0]

println("Test Inv Rosenbrock-W")
f.invW(iW,u,p,2.0,t)
@test minimum(iW - inv(I - 2*J) .< 1e-10)

f.invW_t(iW,u,p,2.0,t)
@test minimum(iW - inv(I/2 - J) .< 1e-10)

println("Parameter Jacobians")
pJ = zeros(Float64,2,4)
f.paramjac(pJ,u,[2.0;2.5;3.0;1.0],t)
@test pJ == [2.0 -6.0  0.0 0.0
             0.0  0.0 -3.0 6.0]

@code_llvm DiffEqBase.has_jac(f)

println("Test booleans")
@test DiffEqBase.has_jac(f) == true
@test DiffEqBase.has_paramjac(f) == true

@code_llvm DiffEqBase.has_paramjac(f)

println("Test difficult differentiable")
NJ = @ode_def_nohes DiffDiff begin
  dx = a*x - b*x*y
  dy = -c*y + erf(x*y/d)
end a b c d
NJ(du,u,[1.5,1,3,4],t)
@test du == [-3.0;-3*3.0 + erf(2.0*3.0/4)]
@test du == NJ(u, [1.5,1,3,4], t)
# NJ(Val{:jac},t,u,J) # Currently gives E not defined, will be fixed by the next SymEgine

test_fail(x,y,d) = erf(x*y/d)
println("Test non-differentiable")
NJ = @ode_def NoJacTest begin
  dx = a*x - μ*x*y
  dy = -c*y + test_fail(x,y,d)
end a μ c d
NJ(du,u,[1.5,1,3,4],t)
@test du == [-3.0;-3*3.0 + erf(2.0*3.0/4)]
@test du == NJ(u,[1.5,1,3,4],t)
@test_throws MethodError NJ.jac(iJ,u,p,t)
# NJ(Val{:jac},t,u,J) # Currently gives E not defined, will be fixed by the next SymEgine

println("Test anonymous definition")

f_t_noname = @ode_def begin # Checks for error due to symbol on 1
  dx = x
  dy = -c*y + d*x*y*t^2
end a b c d

@test DiffEqBase.__has_syms(f_t_noname)

f = @ode_def begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a b c d
@test_nowarn f([0.1,0.2], [1,2], [1,2,3,4], 1)

# Test failures of derivatives should not have #undef
sir_ode = @ode_def SIRModel begin
  dS = -b*S*I
  dI =  b*S*I - g*I
  dR =  g*I
end b g


println("Make the problems in the problem library build")

using DiffEqProblemLibrary
using DiffEqProblemLibrary.ODEProblemLibrary: importodeproblems
using DiffEqProblemLibrary.SDEProblemLibrary: importsdeproblems

importodeproblems()
importsdeproblems()
