type ParameterizedFunction{isinplace,F,P} <: AbstractParameterizedFunction{isinplace}
  f::F
  p::P
end

Base.@pure function ParameterizedFunction(f,p)
  isinplace = numparameters(f)>=4
  ParameterizedFunction{isinplace,typeof(f),typeof(p)}(f,p)
end

(pf::ParameterizedFunction{true,F,P}){F,P}(t,u,du) = pf.f(t,u,pf.p,du)
(pf::ParameterizedFunction{true,F,P}){F,P}(t,u,params,du) = pf.f(t,u,params,du)
(pf::ParameterizedFunction{false,F,P}){F,P}(t,u) = pf.f(t,u,pf.p)
(pf::ParameterizedFunction{false,F,P}){F,P}(t,u,params) = pf.f(t,u,params)

type DAEParameterizedFunction{isinplace,F,P} <: AbstractParameterizedFunction{isinplace}
  f::F
  p::P
end

Base.@pure function DAEParameterizedFunction(f,p)
  isinplace = numparameters(f)>=5
  DAEParameterizedFunction{isinplace,typeof(f),typeof(p)}(f,p)
end

(pf::DAEParameterizedFunction{true,F,P}){F,P}(t,u,du,out) = pf.f(t,u,pf.p,du,out)
(pf::DAEParameterizedFunction{true,F,P}){F,P}(t,u,params,du,out) = pf.f(t,u,params,du,out)
(pf::DAEParameterizedFunction{false,F,P}){F,P}(t,u,du) = pf.f(t,u,pf.p,du)
(pf::DAEParameterizedFunction{false,F,P}){F,P}(t,u,params,du) = pf.f(t,u,params,du)

type DDEParameterizedFunction{isinplace,F,P} <: AbstractParameterizedFunction{isinplace}
  f::F
  p::P
end

Base.@pure function DDEParameterizedFunction(f,p)
  isinplace = numparameters(f)>=5
  DDEParameterizedFunction{isinplace,typeof(f),typeof(p)}(f,p)
end

(pf::DDEParameterizedFunction{true,F,P}){F,P}(t,u,h,du) = pf.f(t,u,h,pf.p,du)
(pf::DDEParameterizedFunction{true,F,P}){F,P}(t,u,h,params,du) = pf.f(t,u,h,params,du)
(pf::DDEParameterizedFunction{false,F,P}){F,P}(t,u,h) = pf.f(t,u,h,pf.p)
(pf::DDEParameterizedFunction{false,F,P}){F,P}(t,u,h,params) = pf.f(t,u,h,params)
