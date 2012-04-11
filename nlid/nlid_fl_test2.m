function e=nlid_fl_test2(db,da,n)
% function e=nlid_fl_test2(db,da,n)
%
% test nlid_fl.m by approximating the relation y=|u|
% for n samples (default n=300) on the interval {u}=[-1,1]  
% by rational function   f(u)=b(u)/a(u)
% using implicit polynomial equality f(v,u)=0, where pf,pr are msspoly
% in z=[v;u], v=msspoly('v',[1 1]), u=msspoly('u',[1 1]) and have same 
% meaning as in nlid_fl.m (default pf=(1+v)*(1+u)^3, pr=(1+u)^3), fl is
% the non-boundedness flag (as in nlid_fl.m)

v=msspoly('v',[1 1]);
u=msspoly('u',[1 1]);
%if nargin<1, pf=(1+v)^2*(1+u)^10; end
%if nargin<2, pr=(1+v+u)^4; end
if nargin<1, pf=(1+v)*(1+u)^4; end
if nargin<2, pr=(1+u)^8; end
if nargin<3, fl=1; end
if nargin<4, T=3; end
if nargin<5, n=300; end

uu=linspace(-1,1,n)';
%vv=uu./max(1,T*abs(uu))+0.00*randn(n,1);
%vv=(1-uu+uu.^2)./(1+uu+uu.^2);
vv=max(0,1-T*abs(uu));
[f,e,r]=nlid_fl([vv uu],pf,pr,1,fl);

if nargout<1,
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
    fprintf('\n nlid_fl error: bound=%f, true=%f  (%d terms)\n', ...
        e,sqrt(mean((vh-vv).^2)),size(p,1))
    close(gcf);
    subplot(2,1,1);plot(uu,vv,uu,vh);grid
    subplot(2,1,2);plot(uu,ee);grid
end