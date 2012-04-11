function ltid_cg_test1(m,n)

if nargin<1, m=310; end
if nargin<2, n=300; end

A=sprandn(m,n,0.1);
B=randn(m,1);
a=@(x)(A*x);
g=@(u)(A'*u);

tic;[x,i]=ltid_cg(a,g,B,1e-9);t=toc;e=norm(A*x-B);
tic;x1=A\B;t1=toc;e1=norm(A*x1-B);

fprintf('\n  cg: |Ax-B|=%f in %f sec (%d steps)  %f',e,t,i)
fprintf('\n  vs: |Ax-B|=%f in %f sec\n',e1,t1)
