function nlid_siso_test2(d,n)
% function nlid_siso_test2(d,n)
%
% Use nlid_siso.m to fit polynomial equations
%  f(y(t),y(t-1),u(t))=0  (deg(f)=2*d+1, default d=1)
% to i/o data of system
%  y(t)=y(t-1)+0.1*tanh(u(t)-y(t-1)) (t=2,...,n, default n=100)

h=0.1;
T=10;
if nargin<1, d=1; end
if nargin<2, n=100; end
y=msspoly('y',2);
u=msspoly('u',1);
q=1+(y'*y+u^2)^d;
tt=linspace(0,T,n)';
uu=sin(tt).^5;
yy=zeros(n,1);
for t=2:n, yy(t)=yy(t-1)+h*uu(t)^2*exp(-yy(t-1)); end
[vv,ww,emsg]=nlid_siso_c2m({[uu yy]},1,0);
if ~isempty(emsg), error(emsg); end
[f,e]=nlid_siso(ww,vv,q);
%x=decomp(f)
yh=nlid_siso_sim(f,q,uu,yy(1),5);
%close(gcf); plot(tt,yy,tt,yh); grid
close(gcf); plot(1:n,yy,1:n,yh); grid
fprintf('\n nlid_siso_sim error: bound=%f, true=%f\n', ...
    e,sqrt(mean((yy-yh).^2)))