function f=isfunction(p,x)
% function f=isfunction(p,x)
%
% p -- msspoly
% x -- free msspoly
%
% f -- 1 if p depends on x and x alone.
%      0 otherwise.


if nargin < 2, error('Two inputs required.'); end

[f,xn] = msspoly.isfreemsspoly(x);

if ~f
    error('Second argument must be a free msspoly.');
end

dep = msspoly.match_list([xn(:);0],p.var);

f = full(all(dep(:) > 0));

end
