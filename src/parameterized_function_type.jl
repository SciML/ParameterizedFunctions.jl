type ParameterizedFunction{isinplace,F,P} <: AbstractParameterizedFunction{isinplace}
  f::F
  p::P
end

function ParameterizedFunction(f,p)
  isinplace = numparameters(f)>=4
  ParameterizedFunction{isinplace,typeof(f),typeof(p)}(f,p)
end

(pf::ParameterizedFunction{true,F,P}){F,P}(t,u,du) = pf.f(t,u,pf.p,du)
(pf::ParameterizedFunction{true,F,P}){F,P}(t,u,params,du) = pf.f(t,u,params,du)
(pf::ParameterizedFunction{false,F,P}){F,P}(t,u) = pf.f(t,u,pf.p)
(pf::ParameterizedFunction{false,F,P}){F,P}(t,u,params) = pf.f(t,u,params)
