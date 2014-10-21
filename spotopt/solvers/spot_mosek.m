function [x,y,z,info] = spot_mosek(A,b,c,K,options)


    if ~isfield(K,'f') 
        K.f = 0;
    end

    if ~isfield(K,'l') 
        K.l = 0;
    end

    if ~isfield(K,'q') || all(K.q == 0)
        K.q = [];
    end

    if ~isfield(K,'s') || all(K.s == 0)
       K.s = [];
    else
       K.s = K.s(K.s ~= 0); 
    end

    if ~isfield(K,'r') || all(K.r == 0)
        K.r = [];
    end


    if nargin < 5, options = spot_sdp_default_options(); end
    
    
    [n,nf,nl,nq,nr,ns] = spot_sdp_cone_dim(K);
    
    %    if ns > 0, error('SPOT: semidefinite constraints not supported for MOSEK.'); end
    

    % First add the non-semidefinite cones.
    evalc('[r,res] = mosekopt(''symbcon'');');
    
    nn = nf+nl+nq+nr;
    prob = struct();
    prob.c = full(c(1:nn)');
    prob.a = A(:,1:nn);
    prob.blc = full(b');
    prob.buc = full(b');
    
    if nr > 0
        prob.cones.type = [];
        prob.cones.subptr = [];
        prob.cones.sub = [];
    end

    if ~isempty(K.q)
        prob.cones.type = repmat(res.symbcon.MSK_CT_QUAD,1,length(K.q));
        prob.cones.subptr = [ 1+[0 cumsum(K.q(1:end-1))] ];
        prob.cones.sub    = [ nf+nl+(1:nq) ];
    end
    if ~isempty(K.r)
        prob.cones.type = [ prob.cones.type repmat(res.symbcon.MSK_CT_RQUAD,1,length(K.r))];
        prob.cones.subptr = [ prob.cones.subptr 1+nq+[0 cumsum(K.r(1:end-1))] ];
        prob.cones.sub    = [ prob.cones.sub nf+nl+nq+(1:nr) ];
    end
    
    function [j,k,l] = jklOf(s,jj)
        j = zeros(size(jj,1),size(jj,2));
        k = zeros(size(jj,1),size(jj,2));
        l = zeros(size(jj,1),size(jj,2));
        
        mtxOff = [0 cumsum(s.^2)];
        for i = 1:length(jj)
            % Which matrix?
            I = find(jj(i)<=mtxOff,1,'first');
            if isempty(I)
                error('Dimension error for indexing.');
            end
            j(i) = I;
            n = s(j(i)-1);
            [k(i),l(i)] = ind2sub([n,n],jj(i)-mtxOff(j(i)-1));
            if k(i) < l(i),
                tmp = l(i);
                l(i) = k(i);
                k(i) = tmp;
            end
        end
        j = j - 1;
    end

    if ~isempty(K.s)
        prob.bardim = K.s;

        barc = c(nn+1:end);
        bara = A(:,nn+1:end);
        
        [jj,~,val] = find(barc(:));
        [j,k,l] = jklOf(K.s,jj');
        
        val(k ~= l) = val(k ~= l)/2;
        
        prob.barc.subj = j;
        prob.barc.subk = k;
        prob.barc.subl = l;
        prob.barc.val  = val';
        
        [jj,ii,val] = find(bara');
        [j,k,l] = jklOf(K.s,jj');
        val(k ~= l) = val(k ~= l)/2;
        
        prob.bara.subi = ii';
        prob.bara.subj = j;
        prob.bara.subk = k;
        prob.bara.subl = l;
        prob.bara.val  = val';
        
    end
    
    if nl > 0
        prob.blx    = [ -Inf*ones(1,nf) ...
                        zeros(1,nl) ...
                        -Inf*ones(1,nq+nr)];
    end
    
    if options.verbose
        cmd = 'minimize info';
    else
        cmd = 'minimize echo(0)';
    end
    
    if isfield(options.solver_options,'mosek')
        param = options.solver_options.mosek;
    else
        param = struct();
    end
    
    start = spot_now();
    if options.verbose
        [r,res] = mosekopt(cmd,prob);
    else
        [info.console,r,res] = evalc('mosekopt(cmd, prob);');
    end
    [info.ctime,info.wtime] = spot_etime(spot_now(),start);

    if ~isfield(res, 'sol')
        status = spotsolstatus.STATUS_SOLVER_ERROR;
    else
        switch res.sol.itr.prosta
          case 'PRIMAL_AND_DUAL_FEASIBLE',
            status = spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE;
          case 'PRIMAL_INFEASIBLE',
            status = spotsolstatus.STATUS_PRIMAL_INFEASIBLE;
          case 'DUAL_INFEASIBLE',
            status = spotsolstatus.STATUS_DUAL_INFEASIBLE;
          case 'PRIMAL_AND_DUAL_INFEASIBLE',
            status = spotsolstatus.STATUS_PRIMAL_AND_DUAL_INFEASIBLE;
          case 'UNKNOWN',
            status = spotsolstatus.STATUS_NUMERICAL_PROBLEMS;
        end
    end

    if spotprogsol.statusIsPrimalFeasible(status)
        x = res.sol.itr.xx;
        if ~isempty(K.s)
            barx = zeros(ns,1);
            moff = 0;
            soff = 0;
            
            map = containers.Map('KeyType','uint32','ValueType','any');
            ss = unique(K.s);
            for i = 1:length(ss)
                n = ss(i);
                V = reshape((1:n^2),n,n);
                [row,col,~] = find(tril(V));
                
                offD = row~=col;
                                
                I = sub2ind([n n],...
                            [ row ; col(offD) ],...
                            [ col ; row(offD) ]);

                map(ss(i)) = sparse(I,[1:length(row) find(offD')]',ones(length(I),1));
            end

            soff = 0;
            moff = 0;
            for i = 1:length(K.s)
                n = K.s(i);
                d = spot_psdDimToNo(n);

                barx(soff + (1:n^2)) = map(n)*res.sol.itr.barx(moff+(1:d));
                soff = soff + n^2;
                moff = moff + d;
            end
            x = [ x ; barx ];

        end
    else
        x = []; % NaN*ones(n,1);
    end
    
    if spotprogsol.statusIsDualFeasible(status)
        y = res.sol.itr.y;
        z = c(:)-A'*y;
    else
        y = []; % NaN*ones();
        z = []; % NaN*ones();
    end
    
    info.solverName = 'mosek';
    if status ~= spotsolstatus.STATUS_SOLVER_ERROR
        info.solverInfo = res.sol;
    end
    info.status = status;

end
