function ltid_mim_test(n,k)

if nargin<1, n=50; end
if nargin<2, k=3; end
z=tf('z');
w=linspace(0.5,2.6,n)';
c0=2.64+rand*0.4;
F0=randn(3,k);
H0=[(1+z)/(1-z) (1-z)/(1+z) (z-1/z)/(z+1/z-2*cos(c0))]*F0;
v0=squeeze(freqresp(H0,w)).'+(0.01*1i)*randn(n,k);
[F,c]=ltid_mim(w,v0);
H=[(1+z)/(1-z) (1-z)/(1+z) (z-1/z)/(z+1/z-2*cos(c))]*F;
vh=squeeze(freqresp(H,w)).';
close(gcf);plot(w,imag(v0),'.',w,imag(vh));grid