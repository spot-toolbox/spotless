function nlid_siso_test1(q,a,T,n)
% function nlid_siso_test1(q,a,T,n)
%
% Test f=nlid_siso(U,Y,q) on U=linspace(-1,1,n)', Y=max(0,1-T*abs(U))
% 
% INPUTS:
%   q  -  msspoly in u=msspoly('u',1) and y=msspoly('y'),
%         default q=(1+u+u^2)^3
%   a  -  msspoly in z=[y;u], default a=msspoly(1)
%   T  -  positive real number, default T=3
%   n  -  positiv integer, default n=300

y=msspoly('y',1);
u=msspoly('u',1);
if nargin<1, q=1+u^4; end
if nargin<2, a=1+u^2; end
if nargin<3, T=3; end
if nargin<4, n=100; end

U=linspace(-1,1,2*n+1)';
Y=max(0,1-T*abs(U));
%Y=1./(1+U.^2);
[f,e]=nlid_siso(U,Y,q,a);

Yh=nlid_siso_sim(f,q,U,zeros(0,1),3);
[x,p]=decomp(f);
fprintf('\n nlid_siso error: bound=%f, true=%f  (%d terms)\n', ...
    e,sqrt(mean((Y-Yh).^2)),size(p,1))
close(gcf);plot(U,Y,U,Yh);grid