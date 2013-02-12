prog = spotsqlprg;

In = eye(n);

[prog,ES] = prog.newPSD(2*n);
E = ES(1:n,1:n);
S = ES(n+(1:n),n+(1:n));
prog = prog.withEqs(In - ES(1:n,n+(1:n)));
prog = prog.withPos(100 - trace(S));

fmonom = [monomials(x,1);monomials(u,1)];
[prog,fcoeff] = prog.newFree(n*length(fmonom));
f = reshape(fcoeff,n,length(fmonom))*fmonom;
F = diff(f,x);

err = E*v - f;
z = zeros(n,1);

H = [ 0   z'   err'
      z   E-Q  F'
      err F    E ];

vars = [ v ; x ; u];
Data = [ V X U ].';

vecH = mss_s2v(H);
Hdata = msubs(vecH,vars,Data);

[prog,r] = prog.newFree(N);
Hdata(1,:) = r.';

[prog,HH] = prog.withBlkPSD(Hdata);

sol_primal = prog.optimize(sum(r),struct('fid',1));