function [x,y,z,info] = spot_mosek_sdsos(A,b,c,K,options)

if  isfield(options,'solveroptions')
    options = options.solveroptions;
else
    options = struct();
end

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


lv = K.l + sum(K.q) + K.r + K.s + K.f; % keyboard;

% Create LP inequality matrices (Aineq*x <= bineq)
Aineq = sparse(K.f+(1:K.l)',K.f+(1:K.l)',-ones(K.l,1),lv,lv);
%Aineq =  sparse([],[],[],lv,lv,0);
%Aineq(K.f+1:K.f+K.l,K.f+1:K.f+K.l) = -speye(K.l);
bineq =  sparse([],[],[],lv,1,0);

Amsk = [Aineq;A];
buc = [bineq;b];
    
blc = [-inf*ones(size(bineq));b];


[r,res] = mosekopt('symbcon');

% Setup non-conic part of problem
prob.c = c;
prob.a = sparse(Amsk);
prob.blc = blc;
prob.buc = buc;

prob.blx = [];
prob.bux = [];

% Setup conic part

prob.cones.sub = zeros(1,3*length(K.q));
prob.cones.subptr = zeros(1,length(K.q));
for k = 1:(length(K.q))
    prob.cones.sub((3*(k-1)+1):(3*(k-1)+3)) = (K.f + K.l) + [3*k-2 3*k-1 3*k];
    prob.cones.subptr(k) = 3*(k-1)+1;
end

prob.cones.type = repmat(res.symbcon.MSK_CT_QUAD, 1, length(K.q));

sprintf('\nDone setting up problem. Running SOCP now...');% keyboard;

% Run mosek problem
start = spot_now();
[r,res]=mosekopt('minimize',prob,options);
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
   


