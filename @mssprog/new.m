function [q,v]=new(h,x,t,s)
% function [q,v]=new(h,x,t,s)
%
% INPUTS:
%   h  -  mss program (class(h)='mssprog')
%   x  -  free msspoly or a positive integer
%   t  -  {'free','pos','lor','rlor','psd'}, default 'free' 
%   s  -  a single character (default '@'), not used when x is msspoly         
%
% OUTPUTS:
%   q  -  updated mss program
%   v  -  msspoly
%
% Registers msspoly v as SeDuMi type t with mss program h (q after update)
% v is the result of converting free msspoly w according to the type t
% (actual conversion required when t='psd', in which case the condition
% length(w)=m*(m+1)/2 must be satisfied for some positive integer m
% w=x when x is a free msspoly, w=msspoly(x,s) otherwise

if nargin<2, error('2 inputs required'); end
p=struct(h);
if nargin<3, t='free'; end
if ~ischar(t), error('3rd argument must be a string'); end
k=strmatch(t,{'free','pos','lor','rlor','psd','hpsd'},'exact');
if isempty(k), error(['variable type "' t '" not supported']); end
if nargin<4, s='@'; end

if isa(x,'double') && prod(size(x)) == 2 && k == 5 & strcmp(s,'@')
    n=x(2)*round(x(1)*(x(1)+1)/2);
    v=reshape(msspoly(s,[n,p.n+1]),round(x(1)*(x(1)+1)/2),x(2));
    p.n=n+p.n;
    q = mssprog(p);
    q = new(q,v,'psd');
    return;
end

if isa(x,'double'),
    N=1;
    x=x(1);
    if x~=max(1,round(x)), error('dimension not a positive integer'); end
    m=x; 
        
    if k==5, 
        n=round(m*(m+1)/2); 
    elseif k==6
        n=round(2*m*(2*m+1)/2); 
    else
        n=m;
    end
    v=msspoly(s,[n,p.n+1]);
    if strcmp(s,'@'), p.n=n+p.n; end
elseif isa(x,'msspoly'),
    if ~isfree(x), error('2nd argument must be a free mss polynomial'); end
    [n,N]=size(x);
    if k==5
        m=round((sqrt(1+8*n)-1)/2);
        if m*(m+1)~=2*n, error('wrong dimension for the psd type'); end
    else
        m=n;
    end
    v=x;
else
    error('2nd argument must be "double" or "msspoly"');
end
NN=n*N;
if (k==5|k==6)&&(m==1), k=2; end
u=[p.v;v(:)];
if ~isfree(u), error('registered variable detected'); end
p.v=u;
p.o=[p.o sum(p.t==k)+(1:NN)];
p.t=[p.t repmat(k,1,NN)];
switch k,
  case 3,
    p.q=[p.q repmat(m,1,N)];
  case 4,
    p.r=[p.r repmat(m,1,N)];
  case 5,
    p.s=[p.s repmat(m,1,N)];
    v=mss_v2s(v(1:n));
end
q=mssprog(p);
