using ParameterizedFunctions
using Base.Test

### ODE Macros

println("Build some examples")
f = @ode_def_nohes SymCheck begin # Checks for error due to symbol on 1
  dx = x
  dy = -c*y + d*x*y
end a=>1.5 b=>1 c=3 d=1

f = @ode_def LotkaVolterra begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a=>1.5 b=>1 c=>3 d=1

f_2 = @ode_def_nohes LotkaVolterra3 begin
  dx = a*x - b^2*x*y
  dy = -c*y + d*x*y
end a=>1.5 b=>1 c=>3 d=1

type  LotkaVolterra2 <: ParameterizedFunction
         a::Float64
         b::Int64
end
(p::LotkaVolterra2)(t,u,du) = begin
         du[1] = p.a * u[1] - p.b * u[1]*u[2]
         du[2] = -3 * u[2] + u[1]*u[2]
end

println("Test Values")
t = 1.0
u = [2.0,3.0]
du = zeros(2)
J = zeros(2,2)
iJ= zeros(2,2)
f(t,u,du)
@test du == [-3.0,-3.0]

println("Test Explicit Parameter Functions")
f(t,u,2.0,du,:a)
@test du == [-2.0,-3.0]
f(t,u,2.0,du,:a,:Deriv)
@test du == [2.0,0.0]
f(t,u,du,[2.0;2.5;3.0])
@test du == [-11.0;-3.0]

println("Test Jacobians")
f(t,u,J,:Jac)
f(t,u,iJ,:InvJac)
@test J  == [-1.5 -2.0
             3.0 -1.0]
@test minimum(iJ - inv(J) .< 1e-10)
pJ = Matrix{Float64}(2,3)
f(t,u,pJ,[2.0;2.5;3.0],:param_Jac)
@test pJ == [2.0 -6.0 0
             0 0 -3.0]

println("Test Hessians")
H = J
iH = iJ
f(t,u,H,:Hes)
@test_throws ErrorException f(t,u,iH,:InvHes)
@test J  == [0.0 0.0
             0.0 0.0]

println("Test using new parameters")
g = LotkaVolterra(a=1.0,b=2.0)
@test g.b == 2.0
@test g.a == 1.0
@test g.a * u[1] - g.b * u[1]*u[2] == -10.0
g(t,u,du)
@test du == [-10.0,-3.0]
h = LotkaVolterra2(1.0,2.0)
h(t,u,du)
@test du == [-10.0,-3.0]

println("Test booleans")
@test f.jac_exists == true
@test f.invjac_exists == true
@test f.hes_exists == true
@test f.invhes_exists == false
@test f.pderiv_exists == true
@test f.pfuncs_exists == true

println("Test non-differentiable")
NJ = @ode_def NoJacTest begin
  dx = a*x - b*x*y
  dy = -c*y + erf(x*y/d)
end a=>1.5 b=>1 c=3 d=4
NJ(t,u,du)
@test du == [-3.0;-3*3.0 + erf(2.0*3.0/4)]
@test NJ.jac_exists == false

### FEM Macros

println("Test FEM")
f = @fem_def (t,x) TestType begin
  du = exp.(-t-5*(1-2x+2x.^2 - 2y +2y.^2)).*(-161 + a*(x - x.^2 + y - y.^2))
end a=400
g = (t,x) -> exp.(-t-5*(1-2x[:,1]+2x[:,1].^2 - 2x[:,2] +2x[:,2].^2)).*(-161 + 400*(x[:,1] - x[:,1].^2 + x[:,2] - x[:,2].^2))
x = rand(4,2)

h = (t,x,u)  -> [1-.5*2.0*u[:,1]   -1-u[:,2]]

l = @fem_def (t,x,u) TestType2 begin
  du = 1-α*β*u
  dv = -1-v
end α=>0.5 β=2.0

@test f(1.0,x) == g(1.0,x)
@test h(1.0,x,x) == l(1.0,x,x)
