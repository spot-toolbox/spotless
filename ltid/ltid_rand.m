function G=ltid_rand(n,k,m)
% function G=ltid_rand(n,k,m)
%
% G is random stable DT system of order 2n with m inputs and k outputs

if nargin<1, n=5; end
if nargin<2, k=3; end
if nargin<3, m=2; end
r=0.95*rand(1,n).^(1/3);
t=pi*rand(1,n);
a=r.*cos(t);
b=r.*sin(t);
A=[diag(a) diag(b);-diag(b) diag(a)];
B=randn(2*n,m);
C=randn(k,2*n);
D=randn(k,m);
G=ss(A,B,C,D,-1);