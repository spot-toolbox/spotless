%
%
% Compare primal and dual standard form solutions for a fitting problem:
%
%  [ r_t  0        e_t'
%    0    E-C'C  F_t'
%    e_t  F_T      E   ]
%
%  e_t = Ex(t+1) - f(x(t),u(t))
%  F_t = df/dx(x(t),u(t))
%

m = 1;
n = 3;

A = randn(3,3);
A = A*0.9/max(abs(eig(A)));

B = randn(n,1)/5;

C = randn(1,n)/5;
Q = C'*C;

R = dlyap(A',Q);

N = 1e5;
us = 0.5*randn(1,N);
xs = zeros(n,N+1);
for i = 1:N
    xs(:,i+1) = A*xs(:,i)  + B*us(:,i);
end


U = us(:,1:N).';
X = xs(:,1:N).';
V = xs(:,2:N+1).';


x = msspoly('x',n);
v = msspoly('v',n);
u = msspoly('u',m);

prog = SPOTSQLProg;

[prog,Evec] = prog.newFree(SPOTSQLProg.psdDimToNo(n));
E = mss_v2s(Evec);

prog = prog.withPSD(E);

% [prog,Svec] = prog.newFree(SPOTSQLProg.psdDimToNo(n));
% S = mss_v2s(Svec);
% In = eye(n);
% prog = prog.withPSD([ E In ; In S]);
% prog = prog.withPos(1 - trace(S));

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


sol_primal = prog.optimizeStandardPrimal(sum(r),struct('fid',1));
sol_dual   = prog.optimizeStandardDual(-sum(r),struct('fid',1));

