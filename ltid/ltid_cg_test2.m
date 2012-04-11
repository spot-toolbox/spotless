function ltid_cg_test2(m,n)
% function ltid_cg_test2(m,n)
%
% testing ltid_cg.m on the task of minimizing |conv(x,q)-B|
% for random m-by-1 q and random (n+m)-by-1 B

if nargin<1, m=2000; end
if nargin<2, n=3000; end

q=randn(m,1);
p=q(m:-1:1);
B=randn(n+m,1);
h=norm(B);

%A=toeplitz([q;zeros(n,1)],[q(1) zeros(1,n)]);
%W=sparse(A);
%a=@(x)full(W'*(W*x));
%b=A'*B;
%tic;[x,i]=ltid_sfls(a,b,h,1e8,500);t=toc;e=norm(conv(q,x)-B);
%fprintf('\n fls: |Ax-B|=%f in %f sec (%d steps)  %f',e,t,i)

s.type='()';
s.subs={m:m+n};

a=@(x)conv(q,x);
g=@(u)subsref(conv(p,u),s);

%a=@(x)subsref(conv(p,conv(q,x)),s);
%b=subsref(conv(p,B),s);

%tic;[x,i]=ltid_sfls(a,b,h,1e6,500);t=toc;e=norm(conv(q,x)-B);
%fprintf('\n fls: |Ax-B|=%f in %f sec (%d steps)  %f',e,t,i)

tic;[x,i]=ltid_cg(a,g,B);t=toc;e=norm(conv(q,x)-B);
fprintf('\n cg: |Ax-B|=%f in %f sec (%d steps)',e,t,i)

A=toeplitz([q;zeros(n,1)],[q(1) zeros(1,n)]);
tic;x1=A\B;t1=toc;e1=norm(A*x1-B);
fprintf('\n vs: |Ax-B|=%f in %f sec\n',e1,t1)