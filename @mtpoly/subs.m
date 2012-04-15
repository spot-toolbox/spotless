function q=subs(p,a,b)
% function q=subs(p,a,b)
%
% p is the result of substituting the entries of the simple mss polynomial 
% column b in place of the corresponding entries of free mss polynomial
% column a within ms polynomial p

% AM 09.01.09

%fprintf('\n****SUBS****\n')
if nargin<3, error('three inputs required'); end
p=mtpoly(p);
a=mtpoly(a);
b=mtpoly(b);
[f,xn]=isfree(a);
if ~f, error('2nd argument must be a free mtpoly column vector'); end
if size(xn,2)~=1, error('input 2 must be a column'); end
na=length(xn);
if na<1, q=p; return; end
[f,x]=issimple(b);
if ~f, error('3rd argument must be a simple mtpoly column vector'); end
if size(x,1)~=size(xn,1), error('unequal number of variables'); end
[ms,ns]=size(p.s);
k=round((ns-3)/2);
vs=p.s(:,3:2+k);
ds=p.s(:,3+k:2+2*k);
cs=p.s(:,ns);
Cs=ones(size(vs));
ee=mss_match(xn,vs);
eee=(ee>0);
eeee=ee(eee);
Cs(eee)=x(eeee,2);
vs(eee)=x(eeee,1);
cs=cs.*prod(Cs.^ds,2);
q=mtpoly(p.m,p.n,[p.s(:,1:2) vs ds cs]);


%for i=1:na,
%    e=(vs==xn(i));
%    Vs(e)=x(i,1);
%    if x(i,1)==0,
%        Cs(e)=x(i,2);
%    end
%end
%q=mtpoly(p.m,p.n,[p.s(:,1:2) Vs p.s(:,3+k:2+2*k) ...
%    p.s(:,ns).*prod(Cs.^(p.s(:,3+k:2+2*k)),2)]);



