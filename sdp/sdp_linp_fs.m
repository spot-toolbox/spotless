function [t,s]=sdp_linp_fs(a)
% function [t,s]=sdp_linp_fs(a)
%
% assuming a>0, maximize t (0<t<2) such that a*[1;-t;-s]>0 ,1] (0<s<2)

a1=a(a(:,3)>0,:);
a2=a(a(:,3)<0,:);
ss=(a(:,3)==0)&(a(:,2)>0);
tmax=min([2,min(a(ss,1)./a(ss,2))]);
t=0;
s=0;
while (tmax-t>0.01)||(t==0),
    tt=0.5*(tmax+t);
    s0=max([0 max((a2(:,1)-a2(:,2)*tt)./a2(:,3))]);
    s1=min([2 min((a1(:,1)-a1(:,2)*tt)./a1(:,3))]);
    if s1<s0,
        tmax=tt;
    else
        t=tt;
        s=0.5*(s0+s1);
    end
end