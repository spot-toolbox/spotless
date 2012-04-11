function ltid_uy2b_test(m,k,n)

if nargin<1, m=3; end
if nargin<2, k=4; end
if nargin<3, n=200; end


b=randn(m,k);
u=randn(n,m);
y=u*b+0.01*(2*rand(n,k)-1);
bs=ltid_uy2b(u,y);
ys=u*bs;
er=sqrt(max(sum((ys-y).^2,2)));
fprintf(' error: %f\n',er)
for i=1:k, close(gcf);plot(ys(:,i),y(:,i));grid;pause;end