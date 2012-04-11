function dsgn_dtsf2_test(T,m,N)

if nargin<1, T=pi/4; end
if nargin<2, m=4; end
if nargin<3, N=300; end

[h,r]=dsgn_dtsf2(T,m,N);
t=linspace(0,pi,5000)';
H=abs(exp((1i*t)*((1-2*m):(2*m-1)))*h);
r0=max(H(t>pi-T));

dr=-20*log10(r);
dr0=-20*log10(r0);

fprintf('\n %2.1fdB~%2.1fdB for %d taps\n',dr,dr0,m)

close(gcf)
subplot(2,1,1); bar((1-2*m):(2*m-1),h); grid
subplot(2,1,2); semilogy(t/pi,H); grid