function nlid_siso_sim_test1(n,T)
% function nlid_siso_sim_test1(n,T)
%
% testing nlid_siso_sim.m by simulating the response y of system
% y''+y/(1+y'^2)=u, y(0)=1, y'(0)=-2 to input u(t)=sin(t) 
% using n samples on interval [0,T] 

if nargin<1, n=100; end
if nargin<2, T=10; end
tt=linspace(0,T,n);
uu=sin(tt)';
h=uu(2)-uu(1); h2=h^2;
yy=zeros(n,1); yy(1)=1; yy(2)=yy(1)-2*h;
for t=3:n,
    yy(t)=(2*yy(t-1)-yy(t-2)+h2*uu(t))/(1+h2/(1+(yy(t-1)-yy(t-2))^2/h2));
end
y=msspoly('y',3);
u=msspoly('u',1);
p=1+(y(2)-y(3))^2/h2;
q=p+h2;
f=q*y(1)-(2*y(2)-y(3)+h2*u)*p;
yh=nlid_siso_sim(f,q,uu,yy(1:2),1);
close(gcf); plot(tt,yy,tt,yh); grid
fprintf('\n nlid_siso_sim error: %f\n',var(yy-yh))