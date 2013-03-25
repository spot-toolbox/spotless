function d=deg(p,x)
if isempty(p.pow), d =0; return; end
if nargin < 2
    d = max(sum(abs(p.pow),2));
else
    if ~msspoly.isfreemsspoly(x), error('2nd argument must be free msspoly'); end
    [~,xn] = isfree(x);
    d=max(sum(abs((msspoly.match_list(xn(:),p.var)>0).*p.pow),2));
end
