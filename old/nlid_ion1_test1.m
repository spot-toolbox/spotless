function nlid_ion1_test1(k,G,M,n)
% function nlid_ion1_test1(k,G,M,n)
%
% test nlid_ion1.m on to find a match of order <=k (default k=1)
% to the i/o data of system G (default G(z)=1/(z-0.5)) with 
% normalized white noise input and measurement noise of amplitude M
% (default M=0), using n samples (default n=100)

if nargin<1, k=1; end
if nargin<2, z=tf('z'); G=1/(z); end
if nargin<3, M=0; end
if nargin<4, n=100; end

u=randn(n,1);
y=lsim(G,u)+M*randn(n,1);
uu=toeplitz(u(k+1:n),u(k+1:-1:1));
yy=toeplitz(y(k+1:n),y(k+1:-1:k));

v=msspoly('v',2);
w=msspoly('w',k+1);
pf=1+sum(v)+sum(w);
ph=msspoly(1);
pr=msspoly(1);
th=1;
[Mu,My,f,r,h,Q,Df,Dr,Dh]=nlid_ion1(uu,yy,pf,pr,ph,th);
f
yh=nlid_ionsim(f,Mu,My,u,y(1),1);
close(gcf);
plot(y,yh,'.',[min(y) max(y)],[min(y) max(y)]); grid