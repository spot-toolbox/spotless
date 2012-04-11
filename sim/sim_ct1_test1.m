function sim_ct1_test1(n,m,N,cse)
% function sim_ct1_test1(n,m,N,cse)
%
% test sim_ct1.m on a system with dim(x)=n, dim(u)=m, 
% with N breakpoints and random generator initial state cse

if nargin<4, cse=mod(round(100*sum(clock)),1000); end
if nargin<3, N=10; end
if nargin<2, m=3; end
if nargin<1, n=2; end
fprintf(' sim_ct1_test1: n=%d, case=%d\n',n,cse)
rand('state',cse)
t=cumsum(rand(1,N));
ut=rand(m,N)-0.5;
us=spline(t,ut);
x0=rand(n,1);

x=msspoly('x',[n 1]); 
u=msspoly('u',[m 1]);
a=x*u'*u*x'+eye(n);
f=a*x*sum(u);
%

[to,yo]=sim_ct1(a,f,us,x0);

v=sim_ppint(spline(t,sum(ut,1)));
yy=x0*exp(-ppval(v,to'));

for i=1:n,
    close(gcf);plot(to,yo(:,i),to,yy(i,:),'.');grid;pause
end