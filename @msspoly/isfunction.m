function f=isfunction(p,x)
% function f=isfunction(p,x)
%
% true if p is a function of free mss polynomial x

if nargin<2, error('2 inputs required'); end
if ~isa(x,'msspoly'), error('2nd input not an mss polynomial'); end
[b,xn]=isfree(x);
if ~b, error('2nd input not free'); end
k=(size(p.s,2)-3)/2;
ee=mss_match([xn(:);0],p.s(:,3:2+k));
f=all(ee(:)>0);