function speedtest(n,m,k,p)
%  solving [Q F;F' 0]*x=y, for x ((n+m)-by-k), given n-by-n Q=Q'>0,
%  n-by-m F, and y ((n+m)-by-k).

if nargin<1, n=1000; end
if nargin<2, m=50; end
if nargin<3, k=3; end
if nargin<4, p=0.01; end
%A=randn(n);F=[eye(m);eye(m);zeros(n-2*m,m)];
A=sprandn(2*n,n,p); F=[speye(m);speye(m);sparse(n-2*m,m)];
y=randn(n+m,k);
Q=A'*A+0.01*speye(n);
clear A
fprintf('nnz(Q)=%d/%d\n\n',nnz(Q),n^2)

tic;
v=(Q+F*F')\[F y(1:n,:)];
w=(F'*v(:,1:m))\(F'*v(:,m+1:m+k)-y(n+1:n+m,:));
x1=v(:,m+1:m+k)-v(:,1:m)*w;
x=[x1;F'*x1+w];
toc
fprintf('|err|=%e\n\n',norm([Q F;F' zeros(m)]*x-y))

tic;
x=[Q F;F' zeros(m)]\y;
toc
fprintf('|err|=%e\n',norm([Q F;F' zeros(m)]*x-y))


