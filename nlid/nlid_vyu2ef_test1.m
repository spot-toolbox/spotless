function nlid_vyu2ef_test1(m,k,N)
% function nlid_vyu2ef_test1(m,k,N)
%
% testing nlid_vyu2ef.m on system  y[t+1]=sin(y[t]+u[t])
% with approximation of y-degree 2m+1, u-degree k, on N-samples data

if nargin<1, m=1; end
if nargin<2, k=3; end
if nargin<3, N=300; end

fprintf('   Preparind test data ...')
Y=2*rand(N,1)-1;                       
U0=2*rand(N,1)-1;
V=sin(Y+U0)+0.01*randn(N,1);
%V=0.5*Y+U0;
U=repmat(U0,1,k).^repmat(1:k,N,1);

fprintf('\n   Running nlid_vyu2ef.m ...\n')
[E,F]=nlid_vyu2ef(V,Y,U,m);            % system id

fprintf('\n   Checking the result ...')
fyu=sum([U ones(N,1)].*(F*(repmat(Y,1,2*m+2).^repmat(2*m+1:-1:0,N,1)).').',2);
Vh=zeros(N,1);
for i=1:N, Vh(i)=sim_findzero([E -fyu(i)]); end
fprintf('\n')
close(gcf); plot(V,Vh,'.'); grid