function nlid_siso_old_test1(pf,pr,ph,n)
% function nlid_siso_old_test1(pf,pr,ph,n)
%
% testing nlid_siso_old.m on y(t)=y(t-1)/(1+0.1*u(t)^2)+u(t-1)

if nargin<4, n=100; end
tt=linspace(0,10*pi,n)';
uu=sin(tt);
yy=zeros(n,1);
for t=2:n, yy(t)=yy(t-1)/(1.1+0.3*uu(t)^2)+0.1*uu(t-1); end
Y=[yy(2:n) yy(1:n-1)];
U=[uu(2:n) uu(1:n-1)];
y=msspoly('y',2);
u=msspoly('u',2);
if nargin<1, pf=1+y(2)+u(1)+u(2); end
if nargin<2, pr=msspoly(1); end
if nargin<3, ph=msspoly(1); end

[f,r,h]=nlid_siso_old(U,Y,pf,pr,ph);

z=[y(2);u];
f0=subs(f,y(1),-1);
f1=subs(f,y(1),1);
yh=zeros(n,1);
for t=2:n,
    F0=double(subs(f0,z,[yh(t-1);uu(t);uu(t-1)]));
    F1=double(subs(f1,z,[yh(t-1);uu(t);uu(t-1)]));
    yh(t)=-(F0+F1)/(F1-F0);
end
close(gcf);plot(tt,yy,tt,yh); grid

