function [a,b]=mint_v2i(x)
% function [a,b]=mint_v2i(x)
%
% for a double integer m-by-n matrix x such that [1 x] is left invertible,
% a (n-by-r) and b (1-by-r) are double integer matrices such that
% v*a<=b iff v is in the convex hull X of the rows of x

if nargin<1, error('one input required'); end
if ~mint_isint(x), error('input is not a double integer'); end
x=mss_unique(x);       % remove repeated rows
[m,n]=size(x);    
if rank([ones(m,1) x])<=n, error('singular input'); end
if n==1,                         % 1-dimensional case
    a=[1 -1];
    b=[min(x) max(x)];
    return; 
end
xx=sum(x,1);
e0=[-1;zeros(n,1)];
K=convhulln(x);                  % the faces of X
r=size(K,1);
a=zeros(n,r);
b=zeros(1,r);
ic=0;
for i=1:r,
    M=[m xx;ones(n,1) x(K(i,:),:)];
    if abs(det(M))>0.5,
        ic=ic+1;
        h=round(abs(det(M))*(M\e0));
        a(:,ic)=h(2:n+1);
        b(ic)=-h(1);
    end
end
a=a(:,1:ic);
b=b(1:ic);