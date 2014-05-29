function [x,y,z,info] = spot_gurobi(A,b,c,K,options)

if  isfield(options,'solveroptions')
    options = options.solveroptions;
else
    options = struct();
end

if isempty(K.l)
    K.l = 0;
end

if ~isempty(K.s)
  error('Gurobi does not support semidefinite constraints')
end

if ~isempty(K.r)
  error('spot_gurobi does not support semidefinite constraints (to be implemented)')
end

if isempty(K.f)
    K.f = 0;
end


lv = K.l + sum(K.q) + sum(K.r) + sum(K.s.^2) + K.f; % keyboard;

% Create LP inequality matrices (Aineq*x <= bineq)
Aineq = sparse(K.f+(1:K.l)',K.f+(1:K.l)',-ones(K.l,1),lv,lv);
%Aineq =  sparse([],[],[],lv,lv,0);
%Aineq(K.f+1:K.f+K.l,K.f+1:K.f+K.l) = -speye(K.l);
bineq =  sparse([],[],[],lv,1,0);

model.obj = full(c);
model.A = [Aineq;A];
model.rhs = full([bineq;b]);
model.lb = -Inf*ones(size(Aineq,2),1);
model.sense = [repmat('<',1,size(Aineq,1)), repmat('=',1,size(A,1))];

for k = 1:(length(K.q))
    model.cones(k).index = (K.f + K.l) + [3*k-2 3*k-1 3*k];
end

% model.cones = pr.cones;

sprintf('\nDone setting up problem. Running Gurobi now...');% keyboard;

tic
output = gurobi(model,options);
info.runtime = toc;

switch output.status
  case {'LOADED', 'INTERRUPTED', 'IN_PROGRESS'}
    status = spotsolstatus.STATUS_UNSOLVED;
  case 'OPTIMAL',
    status = spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE;
  case 'INFEASIBLE'
    status = spotsolstatus.STATUS_PRIMAL_INFEASIBLE;
  case 'UNBOUNDED'
    status = spotsolstatus.STATUS_DUAL_INFEASIBLE;
  case {'CUTOFF', 'NUMERIC', 'SUBOPTIMAL'}
    status = spotsolstatus.STATUS_NUMERICAL_PROBLEMS;
  case {'INF_OR_UBND', 'ITERATION_LIMIT', 'NODE_LIMIT', 'TIME_LIMIT', 'SOLUTION_LIMIT'}
    status = spotsolstatus.STATUS_SOLVER_ERROR;
  otherwise
    status = spotsolstatus.STATUS_SOLVER_ERROR;
end
info.status = status;

if isfield(output,'runtime')
  info.runtime = output.runtime;
else
  %sometimes Gurobi doesn't output a runtime
  info.runtime = NaN;
end

if isfield(output,'x')
  x = output.x;
else
  x = [];
end

y = []; % I need to FIX THESE
z = [];



