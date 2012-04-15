function q=mpmonomials(varargin)
% function q=mpmonomials(x1,p1,x2,p2,...)
%
% q is the column of all monomials of free mss polynomials
% x1,x2,... of degrees specified in p1,p2,... respectively

ni=nargin;
mi=round(ni)/2;
if ni~=2*mi, error('number of inputs not even'); end
if mi<1, error('2 inputs required'); end
v=varargin;
q=mss_mpmonomials(v);
