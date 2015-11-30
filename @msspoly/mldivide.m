function q=mldivide(a,p)
% function q=mldivide(a,p)
%
% division by polynomials is not permitted!

if nargin<2, error('2 inputs required'); end
if ~isa(a,'double'), 
    error(['left division of msspoly by ' class(a) ' is not allowed']); end
b=inv(a);
if ~all(isfinite(b(:))), error('division by non-invertible matrix'); end
q=b*p;
