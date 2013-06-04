function pr = withSOSTrigPoly(pr,p,cs,basis)
%
%  prout = withSOSTrig(pr,p,c,s,basis)
%
%  pr    -- spotsosprog
%  p     -- 1-by-1 msspoly
%  cs    -- k-by-2 free msspoly
%  basis -- l-by-1 msspoly, monomials.
    
    if ~spot_hasSize(p,[1 1]),
        error('Second argument must be scalar msspoly.');
    end
    
    if ~pr.isRealPolyLinearInDec(p)
        error(['Coefficients must be real, indeterminates ' ...
               'non-trigonometric, and expression must ' ...
               'be linear in decision variables.']);
    end
    
    if ~isfree(cs)
        error('Third argument must be free msspoly.');
    end

    pr.sosTrigExpr{pr.numTrigSOS+1} = {p,cs,basis};
end