function ik=match(x,y)
%
%
%   ik = match(x,y)
%
%   x -- free msspoly.
%   y -- n-by-m free msspoly.
%
%   Returns:
%   b -- n-by-m non-negative integers.
%
%   b(i) = k if y(i) = x(k) for some k, 0 otherwise.
%

[f,xn] = msspoly.isfreemsspoly(msspoly(x));

if ~f, error('First argument must be free msspoly.'); end
if nargin < 2, error('Two arguments required.'); end

[f,yn] = msspoly.isfreemsspoly(y);
if ~f, error('Second argument must be free msspoly.'); end

ik=msspoly.match_list(xn,yn);

end
