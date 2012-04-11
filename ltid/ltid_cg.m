function [x,i]=ltid_cg(a,g,B,h)
% function [x,i]=ltid_cg(a,g,B,h)
%
% conjugate gradient for minimizing |A*x-B| with relative accuracy h 
% where A is given by a=@(x)A*x,  g=@(p)A'*p
%
% INPUTS:
%   a  -  function handle: a=@(x)(A*x) for n-by-1 double x
%   g  -  function handle: g=@(x)(A'*x) for m-by-1 double p
%   B  -  n-by-1 double
%   h  -  positive real: tolerance (relative to |B|, default h=1e-8)
%
% OUTPUTS:
%   x  -  n-by-1: |Ax-B| is close to its minimum
%   i  -  number of iteration steps performed

if nargin<3, error('3 inputs required'); end
if nargin<4, h=1e-8; end
if ~isa(a,'function_handle'), error('input 1 not a function handle'); end
if ~isa(g,'function_handle'), error('input 2 not a function handle'); end
if ~isa(B,'double'), error('input 3 not a double'); end
[n,nb]=size(B); if nb~=1, error('input 3 not a column'); end
h=double(h(1)); if h<=0, error('input 4 not positive'); end

% run iterations for x[k], e[k]=B-A*x[k], d[k]=A'*e[k] according to
%    x[k+1]=x[k]-(x[k]-x[k-1])*y[k]+d[k]*z[k],         x[-1]=x[0]=0;
%    e[k+1]=e[k]-(e[k]-e[k-1])*y[k]-A*d[k]*z[k],       e[-1]=e[0]=B
% where y[k],z[k] at k=1,2,3,... are selected to minimize norm(e[k])

i=0;                                 % step number
e=B;                                 % e[i]=B at i=0
d=g(B);                              % d[i]=A'*B at i=0
m=size(d,1);                         % dimension of x
N=min(n,m);                          % estimate of the rank of A
x=zeros(m,1);                        % x[0]=0
en=norm(B);
hB=h*en;                             % absolute tolerance
v=a(d);                              % v=A*A'*B
vn=norm(v);
if vn==0, return; end                % x=0 if AA'b=0
z=(norm(d)/vn)^2;                    % first correction coefficient
x=d*z;
eo=B;
eon=en;
e=eo-v*z;
en=norm(e);
if en>eon,  x=zeros(m,1); return; end   % no progress at all
if en>eon-hB,  return; end              % slow progress
xo=zeros(m,1);
for i=1:N,
    d=g(e);                            % d=A'*e[i]
    v=a(d);                            % v=A*A'*e[i]
    u=e-eo;                            % u=e[i]-e[i-1]
    ue=u'*e;
    ve=v'*e;
    uv=u'*v;
    u2=norm(u)^2;
    v2=norm(v)^2;
    D=u2*v2-abs(uv)^2;
    if D==0, return; end
    y=(ue*v2-ve*uv)/D;
    z=(ve*u2-ue*uv')/D;
    xx=x-(x-xo)*y+d*z;
    ee=e-u*y-v*z;
    een=norm(ee);
    %fprintf('y=%1.3f, |e|=%2.6f, |Ax-B+e|=%e\n',y,en,norm(a(x)-B+e))
    if een>en, return; end
    if een>en-hB, x=xx; return; end
    xo=x;
    x=xx;
    eo=e;
    e=ee;
    en=een;
end