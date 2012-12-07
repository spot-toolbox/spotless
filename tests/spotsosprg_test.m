clear classes

pr = spotsosprg;

x = msspoly('x',2);
A = randn(2,2);
p = x'*A'*A*x;

[pr,token] = pr.withSOS(p);

sol = pr.optimize(0);

[Q,phi] = sol.getSOSDecomposition(token);