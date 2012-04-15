%%% Partial Decomposition of a Polynomial
% [R,p] = pdecomp(q,x)
%
%  Arguments:
%    q  -- An n-by-m mtpoly.
%    x  -- A free mtpoly.
%  Output:
%    R -- P-by-n*m mtpoly, deg(R,x) = 0.
%    p -- P associated power of x.
%
%  Decomposes q into polynomial coefficients of the powers of x.
%
function [R,p] = pdecomp(q,x)
%     if ~all(size(q) == [1 1])
%         error('q must be 1-by-1');
%     else
    if ~isfree(x)
        error('x must be free');
    elseif deg(q) == 0
        p = 0;
        R = q;
        return;
    end
        
    q = reshape(q,prod(size(q)),1);
    xid = x.s;
    xid = xid(3);
    st = q.s;
    n = (size(st,2)-3)/2;
    ind = st(:,1:2);
    ids = st(:,2+(1:n));
    pp  = st(:,2+n+(1:n));
    coeff = st(:,2*(1+n)+1);
    
    msk = ids == xid;
    
    ps = max(pp.*msk,[],2);
    [p,I,J] = unique(ps);
    
    pp = pp.*(ids ~= xid);
    ind(:,2) = J;
    
    p = p';
    R = mtpoly(size(q,1),length(p),[ind ids pp coeff]);
end
