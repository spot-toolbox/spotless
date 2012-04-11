function v=sim_cconv(h,u,d)
% function v=sim_cconv(h,u,d)
%
% circular convolution with column vector h, and delay by d,
% applied column-wise to u
%
% INPUTS:
%   u  -  n-by-k double
%   h  -  m-by-1 double (m<=n)
%   
% OUTPUT:
%   v  -  n-by-1 double
%
% v=[e(d+1:n);e(1:d)] where e=w(m:m+n-1), w=conv(h,[u(n+2-m:n);u;u(1:m)])

if nargin<3, d=0; end
[n,k]=size(u);
[m,r]=size(h);
if r~=1, error('first argument not a column'); end
v=zeros(n,k);
for i=1:k,
    w=conv(h,[u(n+2-m:n,i);u(:,i);u(1:m,i)]);
    e=w(m:m+n-1);
    v(:,i)=[e(d+1:n);e(1:d)];
end