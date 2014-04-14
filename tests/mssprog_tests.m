% A test of mssprog.  Very course.
er = 'Test failed.';

n = 4;
I = eye(n);
A = randn(n,n);

prog = mssprog;
[prog,t] = new(prog,1,'pos');
Q = mss_v2s(msspoly('psd',nchoosek(2*n+1,2)));
prog.psd = mss_s2v(Q);
%[prog,Q] = new(prog,2*n,'psd');
prog.eq = mss_s2v(Q - [ t*I A ; A' t*I ]); % These are things to avoid!
[prog,info] = sedumi(prog,t,0,struct());

if abs(svds(A,1) - double(prog(t))) > 1e-6, error(er); end

%% Simple test w/ LPs
prog = mssprog;
x  = msspoly('x');
mn = monomials(x,0:10);
[prog,coeff] = new(prog,length(mn),'free');
f = coeff'*mn;

K = 100;
tt = linspace(-1,1,K);
yy = double(tt < 0);

err = yy - msubs(f,x,tt);

[prog,q] = new(prog,K+1,'lor');
t = q(1);
prog.eq = q(2:K+1)' - err;

prog = sedumi(prog,t,0,struct());

fopt = prog(f);
plot(tt,msubs(fopt,x,tt),tt,yy)