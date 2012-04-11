function sdp_tau2_test(n)

if nargin<1, n=100; end

p=[randn(n,2) rand(n,1)];
tau=sdp_tau2(p);
N=1000;
t=linspace(0,tau,N); 
y=min(p*[t.^2;t;ones(1,N)],[],1);
close(gcf);plot(t,y);grid