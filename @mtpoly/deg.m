function d=deg(p,x)
% function d=deg(p,x)
%
% degree of p with respect to free mtpoly variable x

k=round((size(p.s,2)-3)/2);
if nargin<2,
    d=max(sum(p.s(:,3+k:2+2*k),2));
    return
end
if ~isa(x,'mtpoly'), error('2nd argument must be mtpoly'); end
[f,xn]=isfree(x);
if ~f, error('2nd argument must be a free mtpoly'); end
d=max(sum((mss_match(xn(:),p.s(:,3:2+k))>0).*p.s(:,3+k:2+2*k),2));
