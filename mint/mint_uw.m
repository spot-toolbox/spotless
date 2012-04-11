function [u,w]=mint_uw(x)
% function [u,w]=mint_uw(x)
%
% for a double integer m-by-n matrix x of rank k find 
% 0/1 matrices u (n-by-k lower triangular) and w (k-by-m upper triangular) 
% such that u'*u=I, ww'=I, and w*x*u is invertible

if nargin<1, error('input required'); end
if ~mint_isint(x), error('input is not a double integer'); end
[m,n]=size(x);         
k=rank(x);
kk=1:k;
u=zeros(n,k);
w=zeros(k,m);
rr=0;
i=0;
while rr<k,
    i=i+1;
    r=rank(x(:,1:i));
    if r>rr,
        rr=r;
        u(i,r)=1;
    end
end
rr=0;
i=0;
while rr<k,
    i=i+1;
    r=rank(x(1:i,:));
    if r>rr,
        rr=r;
        w(r,i)=1;
    end
end