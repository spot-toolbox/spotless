function nlid_fl_old_test1(D,T,n)
% function nlid_fl_old_test1(D,T,n)
%
% test nlid_fl.m by approximating the relation y=sat(u)=u/max(1,T*abs(u))
% (default T=5) for n samples (default n=300) on the interval {u}=[-1,1]  
% using implicit polynomial equality f(v,u)=0, where D is 1-by-3 
% matrix of power settings (default D=[1 1 5])

if nargin<1, D=[1 0;0 5]; end
if nargin<2, T=5; end
if nargin<3, n=300; end

uu=linspace(-1,1,n)';
%vv=uu./max(1,T*abs(uu))+0.00*randn(n,1);
%vv=(1-uu+uu.^2)./(1+uu+uu.^2);
vv=max(0,1-T*abs(uu));
[f,r,U0,U1,Q]=nlid_fl_old([vv uu],D);
v=msspoly('v',[1 1]);
u=msspoly('u',[1 1]);
vh=zeros(n,1);
ee=zeros(n,1);
k=1+2*(deg(f,v)-1);
g=subs(f,u,uu(1));
y=newton(g,v,10,0);
for i=1:n,
    g=subs(f,u,uu(i));
    y=newton(g,v,k,y);
    vh(i)=y;
    ee(i)=double(subs(g,v,y));
end
[x,p]=decomp(f);
p
fprintf('\n true error: %f  (%d terms approximation)\n', ...
    sum((vh-vv).^2),size(p,1))
close(gcf);
subplot(2,1,1);plot(uu,vv,uu,vh);grid
subplot(2,1,2);plot(uu,ee);grid