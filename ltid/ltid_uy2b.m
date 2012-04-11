function b=ltid_uy2b(u,y,o)
% function b=ltid_uy2b(u,y,o)
%
% INPUTS:
%   u  -  n-by-m real
%   y  -  n-by-k real
%
% OUTPUTS:
%   b  -  m-by-k real
%
% Approximates y by  u*b in the uniform metric

if nargin<2, error('2 inputs required'); end
if nargin<3, o='L'; end

[n,m]=size(u);
[n1,k]=size(y);
if n~=n1, error('incompatible inputs'); end

pr=mssprog;
b=msspoly('b',k*m);              % coefficients of b
pr.free=b;
b=reshape(b,m,k);
x=msspoly('x',n*(k+1));
x=reshape(x,k+1,n);
pr.lor=x;
pr.eq=x(1,1:n-1)-x(1,2:n);
pr.eq=u*b-y-x(2:k+1,:)';
pr.sedumi=x(1,1);

b=pr({b});