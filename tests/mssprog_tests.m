% A test of mssprog.  Very course.
er = 'Test failed.';

n = 4;
I = eye(n);
A = randn(n,n);

prog = mssprog;
[prog,t] = new(prog,1,'pos');
[prog,Q] = new(prog,2*n,'psd');
prog.eq = mss_s2v(Q - [ t*I A ; A' t*I ]); % These are things to avoid!
[prog,info] = sedumi(prog,t,0,struct());

if abs(svds(A,1) - double(prog(t))) > 1e-6, error(er); end