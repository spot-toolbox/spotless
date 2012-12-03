function [q,dd]=monomials(x,p)
%
% [q,exp]=monomials(x,p)
%
%  x -- n-by-1 free msspoly.
%  p -- array of non-negative integers.
%
%  q   -- N-by-1 msspoly: all monomials of total deg. == p(i).
%  exp -- N-by-n non-negative integer: vector degree of each monomial.
%

if nargin < 2, error('Two arguments required.'); end
[f,xn]=isfree(x);

if ~f, error('First argument must be free msspoly'); end

xn = xn(:);
nx = length(xn);

if ~spot_isIntGE(p,0), error(['Second argument must be array of non-negative ' ...
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
