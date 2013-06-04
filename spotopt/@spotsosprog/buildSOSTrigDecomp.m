function [pr,Q,phi,y,basis] = buildSOSTrigDecompPrimal(pr,expr,cs,basis)
    if nargin < 3, cs = []; end
    if nargin < 4, basis = []; end
    if ~spot_hasSize(expr,[1 1])
        error('buildSOSTrigDecomp expects a scalar polynomial.');
    end

    decvar = pr.variables;
    
    if ~isempty(cs)
        c = cs(:,1);
        s = cs(:,2);
    
        % Multipliers with dummy variables.
        [~,m,coeff] = pr.newFreeTrigPoly(cs,basis,size(cs,1));

        phi = spotsosprog.buildGramBasis(expr+m.'*(c.^2+s.^2-1),...
                                         [decvar;coeff(:)]);
        phi = spot_trigPolySpan(phi,cs);
    else
        phi = spotsosprog.buildGramBasis(expr,decvar);
    end
    keyboard
    [pr,Q] = pr.newPSD(length(phi));
    
    decvar = [decvar ; mss_s2v(Q)];
    sosCnst = spot_trigPolySpan(expr-phi'*Q*phi,cs);

    A = diff(sosCnst,decvar);
    b = subs(sosCnst,decvar,0*decvar);
    [var,pow,Coeff] = decomp([b A].');
    
    [pr,y] = pr.withEqs(Coeff'*[1;decvar]);
    basis = recomp(var,pow,eye(size(pow,1)));
end
