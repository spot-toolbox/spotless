prg = SPOTSQLProg;

[prg,t] = prg.newFree(1);

n=20;
N = 100;
e = zeros(N,1);
H = msspoly(zeros(SPOTSQLProg.psdDimToNo(n+1),N));
for I = 1:N
    v = ones(n,1);
    A = randn(n,n);
    A = A'*A + 0.1*eye(n);
    H(:,I) = mss_s2v([ t v' ; v A]);
    %     prg = prg.withPSD();
    e(I) = v'*inv(A)*v;
end
prg = prg.withBlkPSD(H);
sol = prg.optimizeStandardPrimal(t,struct('fid',1));
sol = prg.optimizeStandardDual(-t,struct('fid',1));
[max(e) sol.eval(t)]