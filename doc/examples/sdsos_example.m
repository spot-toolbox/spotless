%% 
% This example demonstrates how to setup DSOS and SDSOS constraints. The 
% specific problem considered in this example is the following:
% Given a (quartic) homogeneous polynomial in N variables, we want to find
% a lower bound gamma on the unit sphere (i.e., p(x) >= gamma, for all x with ||x|| = 1). 
% We'd like to find the largest such gamma. Then the DSOS/SDSOS/SOS constraint
% we can write down is:
% p(x) - gamma*(x'x)^2 is DSOS (or SDSOS/SOS)
% So when (x'x) = 1, p(x) - gamma >= 0, and thus p(x) >= gamma.
% This DSOS/SDSOS problem is considered in the paper:
% Amir Ali Ahmadi and Anirudha Majumdar, "DSOS and SDSOS: More Tractable
% Alternatives to Sums of Squares and Semidefinite Programming", In
% Preparation.
%%


% Number of variables (indeterminates) in p
N = 5;
x = msspoly('x',N); % x of dimension N

randn('state',0);

vx = monomials(x,4:4); % Degree 4 homogeneous
% Generate random polynomial
cp = randn(1,length(vx));
p = cp*vx;

%% DSOS program

% Initialize program
prog = spotsosprog;
prog = prog.withIndeterminate(x);

% New free variable gamma
[prog,gamma] = prog.newFree(1);

% DSOS constraint
prog = prog.withDSOS((p - gamma*(x'*x)^2)); % Only line that changes between DSOS,SDSOS,SOS programs

% MOSEK options
options = spot_sdp_default_options();
options.solveroptions.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT'; % Use just the interior point algorithm to clean up
options.solveroptions.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER'; % Don't use basis identification (it's slow)

% Solve program
sol = prog.minimize(-gamma, @spot_mosek, options);

% Optimal value
opt_dsos = double(sol.eval(gamma));

disp(['Optimal value (DSOS): ' num2str(opt_dsos)])

%% SDSOS program (Same as DSOS program, except the line prog = prog.withDSOS(...) is now prog = prog.withSDSOS(...)

% Initialize program
prog = spotsosprog;
prog = prog.withIndeterminate(x);

% New free variable gamma
[prog,gamma] = prog.newFree(1);

% SDSOS constraint
prog = prog.withSDSOS((p - gamma*(x'*x)^2)); % Only line that changes between DSOS,SDSOS,SOS programs

% MOSEK options
options = spot_sdp_default_options();

% Solve program
sol = prog.minimize(-gamma, @spot_mosek, options);

% Optimal value
opt_sdsos = double(sol.eval(gamma));

disp(['Optimal value (SDSOS): ' num2str(opt_sdsos)])


%% SOS program (Same as DSOS/SDSOS program, except the line prog = prog.withDSOS(...) is now prog = prog.withSOS(...)


% Initialize program
prog = spotsosprog;
prog = prog.withIndeterminate(x);

% New free variable gamma
[prog,gamma] = prog.newFree(1);

% SOS constraint
prog = prog.withSOS((p - gamma*(x'*x)^2)); % Only line that changes between DSOS,SDSOS,SOS programs

% MOSEK options
options = spot_sdp_default_options();

% Solve program
sol = prog.minimize(-gamma, @spot_mosek, options);

% Optimal value
opt_sos = double(sol.eval(gamma));

disp(['Optimal value (SOS): ' num2str(opt_sos)])













