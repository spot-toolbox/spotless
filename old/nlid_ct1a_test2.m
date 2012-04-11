function nlid_ct1a_test2(N,d,cse)
% function nlid_ct1a_test2(N,d,cse)
%
% testing nlid_ct1a.m on data samples representing the response of
%   sqrt(1+y^4)*dy/dt=u^2-y*u^2-y^3

if nargin<1, N=10; end 
if nargin<2, d=[2 1]; end
if nargin<3, cse=mod(round(100*sum(clock)),1000); end
fprintf(' nlid_ct1a_test: N=%d,  d=[%d,%d],  case=%d\n',N,d,cse)
rand('state',cse)

uu=spline(0:N,rand(1,N+1));
func=@(t,y)(((1-y)*ppval(uu,t).^2-y^3)/sqrt(1+y^4));
[t,y]=ode45(func,[0 N],0);
y=y';
t=t';
u=ppval(uu,t);
v=((1-y).*u.^2-y.^3)./sqrt(1+y.^4);
[a,F,g]=nlid_ct1a([v;y;u],d);
as=double(msubs(a,msspoly('x',[1 1]),y));
fs=double(msubs(F,[msspoly('x',[1 1]);msspoly('u',[1 1])],[y;u]));
vs=-fs./as;
er=sqrt(sum((v-vs).^2)/length(v));
fprintf('  min(a)=%f,  opt. error = %f,  std error = %f\n',min(as),g,er)
%close(gcf);plot(v,vs,'.');grid;

[ts,ys]=sim_ct1(a,F,uu,0);
close(gcf);plot(t,y,'.',ts,ys);grid
