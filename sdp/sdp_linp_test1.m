function sdp_linp_test1(m,n,cs)
% function sdp_linp_test1(m,n,cs)
%
% test sdp_linp.m and its clones
% on the feasibility of the primal/dual linear programs
%
% B+A*x>0, C'*x->min, y>0, C=A'*y, B'*y->min

if nargin<1, m=100000; end
if nargin<2, n=100; end
if nargin<3, 
    cs=round(10000*rand); 
    fprintf('random seed: %d\n',cs)
end
RandStream.setDefaultStream(RandStream('mt19937ar','seed',cs));

A=sprandn(m,n,0.01);
B=rand(m,1)-A*randn(n,1);
C=A'*rand(m,1);

tic;[x,y,kf,ko] = sdp_linp1(A,B,C);t1=toc;
fprintf('sdp_linp1: %f sec, Cx=%f, chk: %1.1e>0, %1.1e~0, kf=%d, ko=%d\n\n\n', ...
    t1,C'*x,min([A*x+B;y]),norm(C-A'*y),kf,ko)

K.l=m; pars.fid=0;
tic; [y,x]=sedumi(-A',-C,B,K,pars);t1=toc; x=full(x);y=full(y);
fprintf('sedumi:    %f sec, Cx=%f, chk: %1.1e>0, %1.1e~0\n', ...
    t1,C'*x,min([A*x+B;y]),norm(C-A'*y))
