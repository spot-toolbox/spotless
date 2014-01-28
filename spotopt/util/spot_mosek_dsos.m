function [x,y,z,info] = spot_mosek_dsos(A,b,c,K,options)
    
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
    
    if isempty(K.l)
        K.l = 0;
    end
    
    if isempty(K.q)
        K.q = 0;
    end
    
    if isempty(K.r)
        K.r = 0;
    end
    
    if isempty(K.s)
        K.s = 0;
    end
    
    if isempty(K.f)
        K.f = 0;
    end
    
    
    lv = K.l + K.q + K.r + K.s + K.f;
%     v = [ pr.freeVariables
%         pr.posVariables
%         pr.lorVariables
%         pr.rlorVariables];
    % pr.ddVariables];
            
    % Create LP inequality matrices (Aineq*x <= bineq)
    Aineq = sparse(K.f+(1:K.l)',K.f+(1:K.l)',-ones(K.l,1),lv,lv);
    %Aineq =  sparse([],[],[],lv,lv,0);
    %Aineq(K.f+1:K.f+K.l,K.f+1:K.f+K.l) = -speye(K.l);
    bineq =  sparse([],[],[],lv,1,0);
    
    Amsk = [Aineq;A];
    buc = [bineq;b];
    
    blc = [-inf*ones(size(bineq));b];
    
    sprintf('\nDone setting up problem. Running LP now...');
    
    % Run mosek problem
    start = spot_now();
    % save msk_problem_data.mat c Amsk blc buc 
    res = msklpopt(c,Amsk,blc,buc,[],[],options);
    [info.ctime,info.wtime] = spot_etime(spot_now(),start);
    
    
    switch res.sol.itr.prosta
      case 'PRIMAL_AND_DUAL_FEASIBLE',
        info.primalInfeasible = 0;
        info.dualInfeasible   = 0;
      case 'PRIMAL_INFEASIBLE',
        info.primalInfeasible = 1;
        info.dualInfeasible = 0;
      case 'DUAL_INFEASIBLE',
        info.dualInfeasible = 1;
        info.primalInfeasible = 0;
      case 'PRIMAL_AND_DUAL_INFEASIBLE',
        info.primalInfeasible = 0;
        info.dualInfeasible = 1;
      case 'UNKNOWN',
        info.primalInfeasible = 0;
        info.dualInfeasible = 0;
    end
    
    if ~info.primalInfeasible,
         x = res.sol.itr.xx;
    else
        x = [];
    end
    
    y = []; % I need to FIX THESE
    z = [];
   