function [pr,p,coeff,basis] = newSOSTrigPoly(pr,basis,cs,basis2,n)
    if nargin < 4, n = 1; end
    
    [pr,p,coeff,basis] = pr.newFreeTrigPoly(cs,basis,n);
    for i = 1:n
        pr = pr.withSOSTrigPoly(p(i),cs,basis2);
    end
end