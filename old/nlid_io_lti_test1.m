function nlid_io_lti_test1(m,k,G,M,n)
% function nlid_io_lti_test1(m,k,G,M,n)
%
% test nlid_io_lti.m on to find a match of order <m (default 1)
% to the i/o data of system G (default G(z)=1/(z-0.5)) with 
% normalized white noise input and measurement noise of amplitude M
% (default M=0), using n samples (default n=100)

if nargin<1, m=0; end
if nargin<2, k=1; end
if nargin<3, z=tf('z'); G=1/z; end
if nargin<4, M=0; end
if nargin<5, n=100; end

u=randn(n,1);
y=lsim(G,u)+M*randn(n,1);
Gh=nlid_io_lti({[u y]},m,k);
yh=lsim(Gh,u);
w=linspace(0,pi,200);
g=squeeze(freqresp(G,w));
gh=squeeze(freqresp(Gh,w));
close(gcf);
subplot(3,1,1);plot(w,real(g),w,real(gh));grid
subplot(3,1,2);plot(w,imag(g),w,imag(gh));grid
subplot(3,1,3);plot(y,yh,'.',[min(y) max(y)],[min(y) max(y)]); grid