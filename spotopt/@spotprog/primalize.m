function [pr] = toPrimalWithFree(pr)
    if pr.isPrimalWithFree()
        return;
    end
    g = pr.g;
    F = pr.F;
    coneCstrVar = pr.coneCstrVar;
    dualCstrVar = pr.dualCstrVar;
    coneToCstrVar = pr.coneToCstrVar;
    K2 = pr.K2;
            
    pr.g = [];
    pr.F = [];
    pr.coneCstrVar = [];
    pr.dualCstrVar = [];
    pr.coneToCstrVar = [];
    pr.K2 = struct('l',0,'q',[],'r',[],'s',[]);
            
    % Next, we need to move over the coneCstrVar.
    coneVar = pr.coneVar;
    coneToVar = pr.coneToVar;
    
    % Migrate over the linear constraints.
    K1len = [ pr.K1.l
              sum(pr.K1.q)
              sum(pr.K1.r)
              sum(spotprog.psdDimToNo(pr.K1.s)) ];
    K2len = [ K2.l
              sum(K2.q)
              sum(K2.r)
              sum(spotprog.psdDimToNo(K2.s)) ];
    
    % Insert gaps into the old mapping.
    insertPt = 0;
    coneToCstrVar = coneToCstrVar;
    shiftPt = 0;
    for i = 1:length(K1len)
        % Add gap in original variables.
        insertPt = insertPt + K1len(i);
        toShift = coneToVar > insertPt;
        coneToVar(toShift) = coneToVar(toShift)+K2len(i);
        pr.A = [pr.A(:,1:insertPt+pr.numFree) sparse(size(pr.A,1),K2len(i)) pr.A(:,insertPt+pr.numFree+1:end)];
        insertPt = insertPt + K2len(i);
        
        % Move other variables over.
        toShift = coneToCstrVar > shiftPt;
        coneToCstrVar(toShift) = coneToCstrVar(toShift) + K1len(i);
        shiftPt = shiftPt + K1len(i) + K2len(i);
    end
    
    pr.coneVar = [ coneVar ; coneCstrVar ];
    pr.coneToVar = [ coneToVar ; coneToCstrVar ];
    
    % Next we need to square away the equations.
    [~,Ivar] = sort(coneToVar);
    [~,Icstr] = sort(coneToCstrVar);
    [Anew,bnew] = spot_decomp_linear(g-F*[pr.freeVar;coneVar(Ivar)]-coneCstrVar(Icstr),[pr.freeVar;pr.coneVar]);
    
    [~,I] = sort(pr.coneToVar);
    
    nf = pr.numFree;
    
    pr.A = [ pr.A ; Anew(:,[(1:nf)';nf+I])];
    pr.b = [ pr.b ; bnew];
    pr.dualVar = [ pr.dualVar ; dualCstrVar ];
    
    pr.K1.l = pr.K1.l + K2.l;
    pr.K1.q = [ pr.K1.q  K2.q];
    pr.K1.r = [ pr.K1.r  K2.r];
    pr.K1.s = [ pr.K1.s  K2.s];
end