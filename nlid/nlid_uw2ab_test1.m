function nlid_uw2ab_test1(n,k,d,r)
% function nlid_uw2ab_test1(n,k,d,r)
%
% tests nlid_uw2ab.m on n random samples of f(t)=1/(1+r*sum(cos(t).^2))
% where t is k dimensional

if nargin<1, n=200; end
if nargin<2, k=2; end
if nargin<3, d=2; end
if nargin<4, r=1; end

x=rand(n,k);
u=1./(1+r*sqrt(sum(x.^2,2)));
z=msspoly('z',[k 1]);
[aa,bb,dd,L]=nlid_uw2ab(u,acos(x),sum(z)^d,1);
close(gcf);plot(u,nlid_abdx2u(aa,bb,dd,x),'.');grid
