clear classes

pr = spotsosprg;





n = 3;
K = 10;

x = msspoly('x',n);
p = msspoly(zeros(K,1));
v = zeros(K,1);

[pr,gamma] = pr.newFree(1);

for i = 1:K
    A = randn(n+1,n+1);
    Q = A'*A;
    p(i) = [1;x]'*Q*[1;x] - gamma;
    v(i) = Q(1,1) - Q(1,2:end)*inv(Q(2:end,2:end))*Q(2:end,1);
end

[pr,tokens] = pr.withSOS(p);

sol = pr.optimize(-gamma);

%[Q,phi] = sol.getSOSDecomposition(tokens);