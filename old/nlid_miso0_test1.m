function nlid_miso0_test1(p,T,n)
% function nlid_miso0_test1(p,T,n)
%
% test nlid_miso0.m for approximating relation v=sat(w)=w/max(1,T*abs(w))
% for n samples (default n=100) on the interval {w}=[-1,1] (default T=5) 
% using implicit polynomial equality f(v,w)=0, where p is an mp-by-2 
% matrix of power settings (default p=[0 2])

if nargin<1, p=[1 5]; end
if nargin<2, T=5; end
if nargin<3, n=300; end

ww=linspace(-1,1,n)';
vv=ww./max(1,T*abs(ww))+0.00*randn(n,1);
%vv=ww./(1+ww.^2);
f=nlid_miso0([vv ww],p);
v=msspoly('v');
w=msspoly('w',1);

ww=linspace(-2,2,n)';
vv=ww./max(1,T*abs(ww))+0.00*randn(n,1);
vh=zeros(n,1);
ee=zeros(n,1);
if max(p(:,1))==0,
    k=1;
else
    k=3;
end
g=subs(f,w,ww(1));
y=newton(g,v,10,0);
for i=1:n,
    g=subs(f,w,ww(i));
    y=newton(g,v,k,y);
    vh(i)=y;
    ee(i)=double(subs(g,v,y));
end
close(gcf);
subplot(2,1,1);plot(ww,vv,ww,vh);grid
subplot(2,1,2);plot(ww,ee);grid