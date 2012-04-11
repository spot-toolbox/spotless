function H=ltid_th2H(t,h,w,r)
% function H=ltid_th2H(t,h,w,r)
%
% transforming the output of ltid_vwqt2ph.m to stable transfer matrix,
% dropping terms below absolute tolerance r

if nargin<2, error('2 inputs required'); end
if nargin<4,
    rr=zeros(size(t));
else
    w=real(w(:));
    minw=min(w);
    maxw=max(w);
    if (minw<=0)||(maxw>=pi), error('input 3 out of range'); end
    z=exp(1i*minw); 
    rr=abs((z^2+1-2*z*cos(t))/(z^2-1));
    z=exp(1i*maxw); 
    rr=r*min(rr,abs((z^2+1-2*z*cos(t))/(z^2-1)));
end
[n,N]=size(h);
k=round(sqrt(N));
if k^2~=N, error('input 2: number of rows is not a full square'); end
z=tf('z');
H=ss([],[],[],zeros(k,k),-1);
for i=1:length(t),
    Hi=reshape(h(i,:),k,k);
    [V,D]=eig(Hi);
    d=diag(D);
    Vr=V(:,d>rr(i))*diag(sqrt(d(d>rr(i))));
    kr=size(Vr,2);
    if t(i)==0,
        H=H+Vr*(ss((z+1)/(z-1))*Vr');
    else
        H=H+Vr*(ss((z^2-1)/(z^2-2*cos(t(i))*z+1))*Vr');
    end
end