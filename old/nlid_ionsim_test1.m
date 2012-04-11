function nlid_ionsim_test1(n,T)
% function nlid_ionsim_test1(n,T)
%
% testing nlid_ionsim.m by simulating the response y of system
% y''+y/(1+y'^2)=u, y(0)=1, y'(0)=-2 to input u(t)=sin(t) 
% using n samples on interval [0,T] 

if nargin<1, n=100; end
if nargin<2, T=10; end
tt=linspace(0,T,n);
u=sin(tt)';
h=u(2)-u(1); h2=h^2;
y=zeros(n,1); y(1)=1; y(2)=y(1)-2*h;
for t=3:n,
    y(t)=(2*y(t-1)-y(t-2)+h2*u(t))/(1+h2/(1+(y(t-1)-y(t-2))^2/h2));
end
a=max(abs(y));       % normalization coefficients
b=max(abs(u));
My=[a 0];
Mu=[b 0];
v=msspoly('v',3);
w=msspoly('w',1);
f=(1+a^2*(v(2)-v(3))^2/h2)*(a*(v(1)-2*v(2)+v(3))-b*h2*w)+a*h2*v(1);
yh=nlid_ionsim(f,Mu,My,u,y(1:2),1);
close(gcf); plot(tt,y,tt,yh); grid
fprintf('\n nlid_ionsim error: %f\n',var(y-yh))