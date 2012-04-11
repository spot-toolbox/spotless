function dsgn_dtcf2_test(T,m,N)

if nargin<1, T=pi/4; end
if nargin<2, m=4; end
if nargin<3, N=300; end

[h,r]=dsgn_dtcf2(T,m,N);
t=linspace(0,pi,5000)';
H=abs(exp((1i*t)*((-2*m):(2*m)))*h);
r0=max(H((t<T)|(t>pi-T)));

dr=-20*log10(r);
dr0=-20*log10(r0);

fprintf('\n %2.1fdB~%2.1fdB for %d taps\n',dr,dr0,m)

close(gcf)
subplot(2,1,1); bar((-2*m:2*m),h); grid
subplot(2,1,2); semilogy(t/pi,H); grid