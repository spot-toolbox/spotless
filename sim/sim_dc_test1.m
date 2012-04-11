function sim_dc_test1(m,d,r)
% function sim_dc_test1(m,d,r)
%
% test sim_dc.m with m samples, degrees d (1-by-2), noise level r

if nargin<1, m=1000; end
if nargin<2, d=[4 4]; end
if nargin<3, r=0; end

t=linspace(0,2*pi,m)';
x=[cos(3*t) sin(2*t)]+r*randn(m,2);
%x=randn(2*m,2);
[h,u,f]=sim_dc(x,d);
h1=max(abs(f(u)));
fprintf('  h=%f, h1=%f\n',h,h1)
N=300;
t=linspace(0,pi,N);
w=[vec(repmat(t',1,N)) vec(repmat(t,N,1))];
w=w(abs(f(w))<h1,:);
close(gcf)
plot(u(:,1),u(:,2),'.',w(:,1),w(:,2),'.')
