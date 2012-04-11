function y=nlid_ionsim(f,Mu,My,u,y0,r)
% function y=nlid_ionsim(f,Mu,My,u,y0,r)
%
% calculate response y=[y(1);...;y(T)]  (y(t)'s are scalars)
% to input u=[u(1);...;u{T}] of the polynomial model 
%  f(v[t],...,v[t-m],w[t],...,w[t-m])=0, u=Mu*[w;1], y=My*[v;1]
% with initial conditions  y0=[y(1);y(2);...;y(m)]
%
% INPUTS:
%   f  - msspoly of z=[v;w], where v=msspoly('v',m+1), w=msspoly('w',m+1),
%   Mu - 1-by-2 double with Mu(1)>0
%   My - 1-by-2 double with My(1)>0
%   u  - T-by-1 double, T>m
%   y0 - m-by-1 double 
%   r  - positive integer (number of Newton method iterations per step)

if nargin<6, error('6 inputs required'); end
if ~isa(y0,'double'), error('input 5 not a "double"'); end
[m,ny]=size(y0);
if ny~=1, error('input 5 is not a column'); end
if m<1, error('input 5 is empty'); end
if ~isa(u,'double'), error('input 4 not a "double"'); end
[T,nu]=size(u);
if nu~=1, error('input 4 is not a column'); end
if T<m, error('input 4 is too short'); end
if ~isa(My,'double'), error('input 3 not a "double"'); end
if ~isequal(size(My),[1 2]), error('input 3 is not 1-by-2'); end
if My(1)<=0, error('first element of input 3 is not positive'); end
if ~isa(Mu,'double'), error('input 2 not a "double"'); end
if ~isequal(size(Mu),[1 2]), error('input 2 is not 1-by-2'); end
if Mu(1)<=0, error('first element of input 2 is not positive'); end
if ~isa(r,'double'), error('input 6 not a "double"'); end
if ~isscalar(r), error('input 6 is not a scalar'); end
if r~=max(1,round(r)), error('input 6 not a positive integer'); end
v=msspoly('v',m+1);          % variables admissible in f
w=msspoly('w',m+1);
z=[v(2:m+1);w];              
if ~isfunction(f,[v;w]), error('illegal variables in the 1st input'); end
vv=zeros(T,1);               % to be the normalized output
v0=(y0-My(2))/My(1);
v0=v0./max(1,abs(v0));
vv(1:m)=v0;
ww=(u-Mu(2))/Mu(1);          % normalized input
ww=ww./max(1,abs(ww));
fprintf('\nnlid_ionsim.')
for t=m+1:T,
    xt=[vv(t-1:-1:t-m);ww(t:-1:t-m)];
    F=subs(f,z,xt);
    vv(t)=nnewton(F,v(1),r,vv(t-1),diff(F,v(1)));
    if mod(t,100)==0, fprintf('.'); end
end
fprintf('done.\n')
y=My(1)*vv+My(2);
