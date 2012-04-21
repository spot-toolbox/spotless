function [o1,o2]=linear(p,x)
%
%   q=linear(p,x)
%
%   p -- A k-by-1 msspoly.
%   x -- A n-by-1 free msspoly.
%
%  Throws an error if p is not affine in x, otherwise returns
%  q -- A k-by-(n+1) msspoly such that q*[1;x] == p.
%
%
%  [l,q]=linear(p,x)
%
%  Same as above except l = 1 if p is affine in x and l = 0, q = [] if
%  p is not affine in x (no errors thrown).

if isa(p,'double'), p = msspoly(p); end

if ~isa(p,'msspoly') || (~isempty(p) && size(p,2) ~= 1)
    error('First argument must be a k-by-1 msspoly');
end
k = size(p,1);

[f,xn] = msspoly.isfreemsspoly(x);

if ~f || size(x,2) ~= 1
    error('Second argument must be a n-by-1 free msspoly');
end
n = size(x,1);

if isempty(p) || isempty(x) || deg(p,x) == 0
    l = 1;
    q = [ p zeros(k,n) ];
elseif deg(p,x) == 1
    l = 1;
    q = [ subs(p,x,0*x) diff(p,x) ];
else
    l = 0; q = [];
end

if nargout == 1
    if ~l
        error('First argument not affine in second argument.');
    else
        o1 = q;
    end
else
    o1 = l; o2 = q;
end


end



