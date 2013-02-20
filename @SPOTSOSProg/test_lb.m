% x = msspoly('x',2);
%  p = (4-2.1*x(1)^2+x(1)^4/3)*x(1)^2 + x(1)*x(2) + (-4+4*x(2)^2)*x(2)^2;

% x = msspoly('x');
% p = x^2 - 1;
n=10;
x = msspoly('x',n);

A = rand(n,n);

pmonom = monomials(x,0:2);
p = subs(rand(1,length(pmonom))*pmonom,x,(1+x).^2);



prg = SPOTSOSProg;
[prg,s] = prg.newFree(1);
prg = prg.withSOS(p + s);
psol = prg.minimizePrimalForm(s,struct('solver','sedumi','solver_options',struct('fid',1)));
dsol = prg.minimizeDualForm(s,struct('solver','sedumi','solver_options',struct('fid',1)));
double([psol.eval(s) dsol.eval(s)])


sdpvar x1 x2 x3 x4 x5 x6 x7 x8 x9 x10;
sdpvar r;

F = sos(p + r);

solvesos(F,r,sdpsettings('sos.model',1),r)

x = msspoly('x',4); x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4);
prg = SPOTSOSProg;
[prg,s] = prg.newFree(1);
prg = prg.withSOS(p + s);
psol = prg.minimizePrimalForm(s,struct('solver','sedumi','solver_options',struct('fid',1)));
