function sim_findzero_test(f,t,k)
% function sim_findzero_test(f,t,k)
%
% testing sim_findzero.m 
% on polynomial f (default f=[1 3 4 1]) 
% initializing at t (default 0)
% using k (default k=5) iterations

if nargin<1, f=[1 3 4 1]; end
if nargin<2, t=0; end
if nargin<3, k=5; end

n=length(f)-1;
g=f(1:n).*(n:-1:1);
p=@(x)[polyval(f,x) polyval(g,x)];
y=sim_findzero(p,t,k);
fprintf(' root: %f,  residual: %e\n',y,polyval(f,y))