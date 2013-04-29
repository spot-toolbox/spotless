% Problem Data:  A random symmetric matrix.
n = 10;
A = randn(n,n);
A = A + A';

prog = spotprog;
[prog,P] = prog.newPSD(n);
[prog,obj] = prog.newFree(1);

prog = prog.withLor([ obj ; P(:)-A(:)]);

sol = prog.minimize(obj);

Popt = double(sol.eval(P));