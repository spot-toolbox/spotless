function [R,p] = pdecomp(q,x)
% [R,p] = pdecomp(q,x)
%
%  Arguments:
%    q  -- An n-by-m msspoly.
%    x  -- A 1-by-1 free msspoly.
%  Output:
%    R -- n*m-by-P msspoly, deg(R,x) = 0.
%    p -- P associated power of x.
%
%  Decomposes q into polynomial coefficients of the powers of x.
%
    
    [f,xn] = isfree(x);
    
    if ~f || ~spot_hasSize(x,[1 1])
        error('Second argument must be free 1-by-1 msspoly.'); 
    end
    
    msk = q.var == xn;
    J = sub2ind(q.dim,q.sub(:,1),q.sub(:,2));
    P = max(max(q.pow.*(msk),[],2),0) + min(min(q.pow.*(msk),[],2),0);
    
    [p,~,I] = unique(P);

    if isempty(p) 
        p = 0;
        R = q(:);
    else
        pow = q.pow.*(~msk);
        var = q.var;
        
        R = msspoly([prod(q.dim) length(p)],...
                    [J I],...
                    var,pow,q.coeff);
    end
end
    
    


