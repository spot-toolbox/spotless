function q=linear(p,x)
% function q=linear(p,x)
%
% q='?' if column p is not linear with respect to free mss vector x, 
% otherwise p(y,x)=q(y)*[1;x]

[mp,np]=size(p);
if np~=1, error('1st argument must be a column vector'); end
if ~isa(x,'mtpoly'), error('2nd argument must be an mss polynomial'); end
[f,xn]=isfree(x);
nx=length(xn);
if ~f, error('2nd argument must be a free mss polynomial'); end
if size(xn,2)~=1, error('input 2 must be a column'); end
if deg(p)<1, q=[p zeros(mp,nx)]; return; end
[ms,ns]=size(p.s);
k=round((ns-3)/2);
vs=p.s(:,3:2+k);
ds=p.s(:,3+k:2+2*k);
ee=mss_match(xn,vs);
if any(ds(ee>0)>1), q='?'; return; end
s=[p.s(:,1) 1+max(ee,[],2) vs.*(ee==0) ds p.s(:,ns)];
q=mtpoly(mp,nx+1,s);
