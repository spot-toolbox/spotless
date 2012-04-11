function nlid_siso_lti_test1(m,k,G,M,n,c0)
% function nlid_siso_lti_test1(m,k,G,M,n,c0)
%
% test nlid_sisolti.m on to find a good LTI match 
%   y(t)=v(t)+c,  v(t)+a1*v(t-1)+...+am*v(t-m)=b0*u(t)+...+bk*u(t-k),
% to the i/o data of system G (default G(z)=1/z) with 
% normalized white noise input and measurement noise of amplitude M
% (default M=0), using n samples (default n=100), and offset c0 (default 7)

if nargin<1, m=1; end
if nargin<2, k=10; end
if nargin<3, z=tf('z'); G=1/(z-0.8)^2; end
if nargin<4, M=0; end
if nargin<5, n=100; end
if nargin<6, c0=7; end

u=randn(n,1);
y=lsim(G,u)+M*randn(n,1)+c0;
[Gh,c,e]=nlid_siso_lti({[u y]},m,k);
yh=lsim(Gh,u)+c;
fprintf('\n nlid_sisolti error: bound=%f, actual=%f\n',e,sqrt(mean((yh-y).^2)))
w=linspace(0,pi,200);
g=squeeze(freqresp(G,w));
gh=squeeze(freqresp(Gh,w));
close(gcf);
subplot(3,1,1);plot(w,real(g),w,real(gh));grid
subplot(3,1,2);plot(w,imag(g),w,imag(gh));grid
subplot(3,1,3);plot(y,yh,'.',[min(y) max(y)],[min(y) max(y)]); grid