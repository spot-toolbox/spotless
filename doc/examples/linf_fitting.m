% Problem Data:  absolute value evaluated at 101 points.
d = 5;
N = 50;
tt = linspace(-1,1,2*N+1);
yy = abs(tt);


% Fit a degree d polynomial.
prog = spotprog;

t = msspoly('t');
basis = monomials(t,0:d);
[prog,coeff] = prog.newFree(length(basis));
f = coeff'*basis;

err = yy - msubs(f,t,tt);

[prog,obj] = prog.newFree(1);
prog = prog.withPos(obj - err');
prog = prog.withPos(obj + err');

sol = prog.minimize(obj);

fopt = sol.eval(f);