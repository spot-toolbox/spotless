function sim_colb_test(d,r,n)

if nargin<1, d=5; end
if nargin<2, r=0.1; end
if nargin<3, n=200; end

t=linspace(0,2*pi,n)';
NN=(1:2*n)';
u=[cos(t);min(t,pi-t)];
y=u+0.1*(randn(2*n,1)-0.5);
v=0.5*(sim_colb(y+r*NN.^2,d)-sim_colb(-y+r*NN.^2,d));
close(gcf);
subplot(2,1,1); plot(NN,y,NN,v); grid
subplot(2,1,2); plot(NN,u,NN,v); grid
