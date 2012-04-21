function q=monomials(x,p)
%
% q=monomials(x,p)
%
%  x -- free msspoly.
%  p -- array of non-negative integers.
%
%  q -- All monomials of total degree equal to p(i) for some i.
%

if nargin < 2, error('Two arguments required.'); end
[f,xn]=isfree(x);

if ~f, error('First argument must be free msspoly'); end
xn = xn(:);
nx = length(xn);

if ~msspoly.isIntGE(p,0), error(['Second argument must be array of non-negative ' ...
                        'integers']); 
end

if isempty(p)
    q = msspoly();
else
    ps = unique(p(:));
    np=length(ps);
    
    dd=mss_asd(nx,p);
    mq=size(dd,1);
    
    q = msspoly([ mq  1],...
               [ (1:mq)' ones(mq,1) ],...
               repmat(xn',mq,1),...
               dd, ones(mq,1));
end
end
