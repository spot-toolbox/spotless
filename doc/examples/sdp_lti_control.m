% Problem data matrices.
A = [2 1 0 ; 0 2 1 ; 0 0 2];
B = [ 0 0 1]';
n = size(B,1);
m = size(B,2);

prog = spotprog;
[prog,S] = prog.newPSD(n);
[prog,L] = prog.newFree(n*m);
L = reshape(L,m,n);

rho = 0.1;
prog = prog.withPSD([ (1-rho)*S (A*S+B*L) ; (A*S+B*L)' S]);

sol = prog.minimize(0);

S = double(sol.eval(S));
L = double(sol.eval(L));
K = L/S;