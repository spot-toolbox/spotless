function b=match(x,y)
% function b=match(x,y)
%
% for free mtpoly x,y, b(i)=k if y(i)=x(k), b(i)=0 otherwise
[bx,xn]=isfree(x);
if ~bx, error('input 1 is not free'); end
if size(xn,2)~=1, error('input 1 must be a column'); end
if nargin<2, error('2 inputs required'); end
if ~isa(y,'mtpoly'), error('input 2 is not an "mtpoly"'); end
[by,yn]=isfree(y);
if ~by, error('input 2 is not free'); end
b=mss_match(xn,yn);
