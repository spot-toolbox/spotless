function ltid_sfls_test1(m,n)

if nargin<1, m=310; end
if nargin<2, n=300; end

A=randn(m,n);
B=randn(m,1);
a=@(x)(A'*(A*x));
b=A'*B;
h=norm(B);

tic;[x,i]=ltid_sfls(a,b,h,1e8,500);t=toc;e=norm(A*x-B);
tic;x1=A\B;t1=toc;e1=norm(A*x1-B);

fprintf('\n fls: |Ax-B|=%f in %f sec (%d steps)  %f',e,t,i)
fprintf('\n  vs: |Ax-B|=%f in %f sec\n',e1,t1)
