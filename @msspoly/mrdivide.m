function q=mrdivide(p,a)
% function q=mrdivide(p,a)
%
% division by polynomials is not permitted!

if nargin<2, error('2 inputs required'); end
if ~isa(a,'double'), 
    error(['division of msspoly by ' class(a) ' is not allowed']); end
b=inv(a);
if ~all(isfinite(b(:))), error('division by non-invertible matrix'); end
q=p*inv(a);
