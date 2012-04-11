function ltid_tpmin_test(n,N)
% function ltid_tpmin_test(n,N)
%
% testing ltid_tpmin.m on random trigonometric polynomials of degree n

if nargin<1, n=20; end
if nargin<2, N=500; end
a=randn(n,1);
[y,t]=ltid_tpmin(a);
tt=linspace (0,pi,N);
yy=cos(repmat(tt',1,n).*repmat(0:n-1,N,1))*a;
close(gcf);
plot(tt,yy,tt,y,[t t],[min(yy) max(yy)])
grid