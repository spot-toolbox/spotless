function [pr,poly,coeff,basis] = newFreeTrigPoly(pr,cs,basis,n)
    
    if ~isfree(cs)
        error('cs  must be free msspoly.');
    end
    
    if ~isempty(cs) && size(cs,2) ~= 2
        error('cs must have two columns.');
    end
    
    if ~isempty(cs)
        basis = spot_trigPolySpan(basis,cs);
    end
    
    [pr,poly,coeff] = pr.newFreePoly(basis,n);
    
end