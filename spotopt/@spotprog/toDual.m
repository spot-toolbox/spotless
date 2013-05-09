function [dl,dobj] = toDual(pr,pobj)
%
%  [dl,dobj] = toDual(pr,pobj)
%
%  (pr,dl) are the primal and dual optimization problems given below:
%
% Primal:
%
%    min (c1,f)+(c2,x), 
%    x in K1,  
%    A1*f+A2*x = b, 
%    g - F1*f-F2*x = r in K2.
%
% Dual:
%  
%    max. (y,b) - (q,g),
%    q in K2*,
%    c2-A2'y+F2'*q = z in K1*,
%    c1 = A1'*y-F1'*q
%
    
%%
%%
%%  NEED TO IMPLEMENT ADJOINT CORRECTLY
%%
%%
    
    
    f = pr.freeVar;
    x = pr.varToDec(pr.coneVar,pr.coneToVar);
    nf = length(f);
    [c,~] = spot_decomp_linear(pobj,[f;x]);
    
    c1 = c(1:nf);   c2 = c(nf+1:end);
    A1 = pr.A(:,1:nf); A2 = pr.A(:,nf+1:end);
    F1 = pr.F(:,1:nf); F2 = pr.F(:,nf+1:end);
    
    P1diag = spotprog.coneInnerProduct(pr.K1);
    P1inv = spdiags(1./P1diag,0,length(P1diag),length(P1diag));
    P2diag = spotprog.coneInnerProduct(pr.K2);
    P2 = spdiags(P2diag,0,length(P2diag),length(P2diag));
    dl = spotprog(pr.name);
    dl.K1 = pr.K2;
    dl.K2 = pr.K1;
    dl.b = -c1';
    dl.A = -[ A1' -F1'*P2 ];
    dl.g = P1inv*(c2');
    dl.F = P1inv*[ A2' F2'*P2 ];
    
    % Now, move around the indices so that sol.dualEval on the old variables will
    % give the same answer.
    
    dl.freeVar   = pr.dualVar;
    dl.coneVar   = pr.coneCstrVar;
    dl.coneToVar = pr.coneToCstrVar;
    dl.dualVar   = pr.freeVar;
    dl.coneCstrVar   = pr.coneVar;
    dl.coneToCstrVar = pr.coneToVar;
    
    y = pr.dualVar; % Associated with equations.
    q = spotprog.varToDec(pr.coneCstrVar,pr.coneToCstrVar);

    dobj = msspoly(0);
    if ~isempty(pr.b)
        dobj = dobj + pr.b'*y;
    end
    if ~isempty(pr.g)
        dobj = dobj - pr.g'*P2*q;
        %
        % g in K2
        %
    end
    
end