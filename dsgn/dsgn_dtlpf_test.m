function dsgn_dtlpf_test(m,N,T)


if nargin<1, m=200; end
if nargin<2, N=2000; end
if nargin<3, 
    T=0;
    h=dsgn_dtlpf(m,N);
else
    h=dsgn_dtlpf(m,N,T);
end
t=linspace(0,pi,5000)';
H=exp((1i*t)*(-m:m))*h;
r0=abs(H(t<pi-T)-1i*t(t<pi-T));

fprintf('\n L2: %2.1fdB, LInf: %2.1fdB\n', ...
    20*log10(norm(r0)/sqrt(N)),20*log10(max(r0)/sqrt(N)))

close(gcf)
subplot(2,1,1); bar((-m:m),h); grid
subplot(2,1,2); plot(t/pi,imag(H),t/pi,t); grid