function [t,y]=sim_ct1(a,f,u,x0)
% function [t,y]=sim_ct1(a,f,u,x0)
%
% simulator for polynomial ODE models
%
% INPUTS:
%  a  - n-by-n msspoly of [x;u],x=msspoly('x',[n 1]),u=msspoly('u',[m 1]))
%  f  - n-by-1 msspoly of [x;u] 
%  u  - real-valued piecewise polynomal of dimension m
%  x0 - n-by-1 real vector 
%
% OUTPUT:
%   y  - n-by-length(u.breaks) real
%
% DESCRIPTION: 
%   yo(:,i)=x(to(i)) are samples of numerical (ode45) solution of the ODE
%      a(x(t),u(t))x'+f(x(t),u(t))=0 with x(t0)=x0 

if ~isa(a,'msspoly'), error('input 1 not a msspoly'); end
[n,na]=size(a);
if na~=n, error('input 1 not a square matrix'); end
if ~isa(f,'msspoly'), error('input 2 not a msspoly'); end
[mf,nf]=size(f);
if mf~=n, error('inputs 1,2 not compatible'); end
if nf~=1, error('input 2 not a column'); end
if ~sim_isppvec(u), error('input 3 not a pp column'); end
t0=u.breaks(1);
t1=u.breaks(end);
m=u.dim;
if ~isa(x0,'double'), error('input 4 not a double'); end
if ~isreal(x0), error('input 4 not real'); end
if ~isequal([n 1],size(x0)), error('inputs 1,4 not compatible'); end
xu=[msspoly('x',[n 1]);msspoly('u',[m 1])];
if ~isfunction(a,xu), error('invalid a'); end
if ~isfunction(f,xu), error('invalid f'); end 
[t,y]=ode45(@(T,Y)ctsimi(f,a,xu,[Y;ppval(u,T)]),[t0 t1],x0);

function v=ctsimi(f,a,xu,z)
v=-double(subs(a,xu,z))\double(subs(f,xu,z));