function [x,y,z,info] = spot_gurobi_dsos(A,b,c,K,options)
    
    if isfield(options,'solveroptions')
       options = options.solveroptions;
    else
       options = struct();
    end
    
%     %  First, construct structure with counts of SeDuMi
%     %  variables.
%     K = struct();
%     K.f = pr.freeNum;
%     K.l = pr.posNum;
%     K.q = pr.lorDim;
%     K.r = pr.rlorDim;
%     K.s = pr.psdDim;
%     % K.dd = pr.ddDim;    
    
    v = [ pr.freeVariables
        pr.posVariables
        pr.lorVariables
        pr.rlorVariables];
    % pr.ddVariables];
            
    % Create LP inequality matrices (Aineq*x <= bineq)
    Aineq =  sparse([],[],[],length(v),length(v),0);
    Aineq(K.f+1:K.f+K.l,K.f+1:K.f+K.l) = -speye(K.l);
    bineq =  sparse([],[],[],length(v),1,0);    
    
    model.obj = full(c);
    model.A = [Aineq;A];
    model.rhs = full([bineq;b]);
    model.lb = -Inf*ones(size(Aineq,2),1);
    model.sense = [repmat('<',1,size(Aineq,1)), repmat('=',1,size(A,1))];
    
    sprintf('\nDone setting up problem. Running LP now...');
    
    output = gurobi(model,options);
    if strcmp(output.status,'OPTIMAL')
        info.primalInfeasible = 0;
        info.dualInfeasible = 0;
        x = output.x;
    else
        info.primalInfeasible = 0;
        info.dualInfeasible = 0;
        x = NaN*ones(length(varNo));
    end
    
    y = []; % I need to FIX THESE
    z = [];
   