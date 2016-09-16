using ParameterizedFunctions
using Base.Test

f = @ode_def LotkaVoltera begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a=>1.5 b=>1 c=3 d=1


type  LotkaVoltera2 <: ParameterizedFunction
         a::Float64
         b::Int64
end
(p::LotkaVoltera2)(t,u,du) = begin
         du[1] = p.a * u[1] - p.b * u[1]*u[2]
         du[2] = -3 * u[2] + u[1]*u[2]
end

t = 1.0
u = [2.0,3.0]
du = zeros(2)
f(t,u,du)
@test du == [-3.0,-3.0]
g = LotkaVoltera(1.0,2.0)
g(t,u,du)
@test du == [-10.0,-3.0]
h = LotkaVoltera2(1.0,2.0)
h(t,u,du)
@test du == [-10.0,-3.0]
