function d=deg(p,x)
% function d=deg(p,x)
%
% degree of p with respect to free msspoly variable x

k=round((size(p.s,2)-3)/2);
if nargin<2,
    d=max(sum(p.s(:,3+k:2+2*k),2));
    return
end
if ~isa(x,'msspoly'), error('2nd argument must be msspoly'); end
[f,xn]=isfree(x);
if ~f, error('2nd argument must be a free msspoly'); end
d=max(sum((mss_match(xn(:),p.s(:,3:2+k))>0).*p.s(:,3+k:2+2*k),2));