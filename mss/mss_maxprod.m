function x=mss_maxprod(M)
% function x=mss_maxprod(M)
%
% For m-by-n real non-negative matrix M, 
% x is the positive real vector maximizing prod(x) subject to Mx<1. 
% Zero columns of M are ignored, and the
% corresponding entries of x are set to Inf.
%
% Uses SeDuMi

error('mss_maxprod is not finished yet')

if nargin<1, error('1 input required'); end
if ~isa(M,'double'), error('1st input not a "double"'); end
if any(M(:)<0), error('1st input has negative elements'); end
[m,n]=size(M);
z=all(M==0,1);
if all(z),               % all zero elements
    x=repmat(Inf,n,1); 
    return
end
if n==1,
    x=1/max(M);
    return
end
if any(z),               % some zero columns
    x=repmat(Inf,n,1);
    x(~z)=mss_maxprod(M(:,~z));
    return
end
k=ceil(log2(n));
N=2^k;
if n~=N,                 % dimension of x is not a power of 2
    x=mss_maxprod([M zeros(m,N-n);zeros(N-n,n) eye(N-n)]);
    x=x(1:n);
    return
end
% With n=2^k, use SeDuMi with m non-negative variables y=1-M*x, 
% and n-1 rotated Lorentz cones v1{i,j}*v2{i,j}>v3{i,j}^2, where
% i=1,...,k, j=1,...,2^(k-i), v1{1,j}=x(2j-1), v2{i,j}=x(2j),
% v1{i+1,j}=v3{i,2j-1}, v2{i+1,j}=v3{i,2j}, v3{k,1}->max
% the overall vector xx of decision variables is
% [y;v{1,1};v{1,2};...;v{1,n/2};v{2,1};...;v{2,n/4};...;v{k,1}]
% where v{i,j} is the vector three scalar decision variables starting at
% (m+3*(2^k+...+2^(k-i+1))+3*j)-th position in xx. 
% Total: m+3*(n-1) variables, m+(n-2) equality constraints
A0=[speye(m) M]; 
B0=ones(m,1);                % y+Mx=1
ii=zeros(3*(n-2),1); 
jj=ii;
cc=ii;
ra=0;                        % so far, ra positions in ii,jj,cc are filled
rx=m;                        % current shift in xx
for i=1:k-1,
    ki2=2^(k-i-1);             % number of cones at level i+1
    ii(ra+(1:2*ki2))=rx+3*(1:2*ki2)';  % third

    



