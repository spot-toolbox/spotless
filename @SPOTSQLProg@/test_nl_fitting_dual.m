
prog = SPOTSQLProg;

[prog,Evec] = prog.newFree(SPOTSQLProg.psdDimToNo(n));
E = mss_v2s(Evec);

[prog,Svec] = prog.newFree(SPOTSQLProg.psdDimToNo(n));
S = mss_v2s(Svec);

In = eye(n);
prog = prog.withPSD([ E In ; In S]);
prog = prog.withPos(100 - trace(S));

fmonom = [monomials(x,1);monomials(u,1)];
[prog,fcoeff] = prog.newFree(n*length(fmonom));
f = reshape(fcoeff,n,length(fmonom))*fmonom;
F = diff(f,x);

[prog,r] = prog.newFree(N);

err = E*v - f;
z = zeros(n,1);

H = [ 0   z'   err'
      z   E-Q  F'
      err F    E ];

vars = [ v ; x ; u];
Data = [ V X U ].';

vecH = mss_s2v(H);
Hdata = msubs(vecH,vars,Data);
Hdata(1,:) = r';

prog = prog.withBlkPSD(Hdata);


sol_dual   = prog.optimizeStandardDual(-sum(r),struct('solver','sedumi',...
                                                  'solver_options',struct('fid',1)));
