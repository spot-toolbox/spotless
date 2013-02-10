pr = SPOTSQLProg;

xx = linspace(-1,1,100)';
yy = sin(pi*xx);%+ 0.1*randn(size(xx));

x = msspoly('x');
mn = monomials(x,0:5);
[pr,coeff] = pr.newFree(length(mn));

yhat = mn'*coeff;

err = msubs(yhat,x,xx')' - yy;

[pr,t] = pr.newFree(1);

[pr1,u] = pr.newPos(length(err));
[pr1,l] = pr1.newPos(length(err));
pr1 = pr1.withEqs(u - (t - err));
pr1 = pr1.withEqs(l - (t + err));

pr2 = pr.withPos(t - err);
pr2 = pr2.withPos(t + err);

sol1 = pr1.optimizeStandardDual(-t);
sol2 = pr2.optimizeStandardPrimal(t);

[sol1.eval(t) sol2.eval(t)]