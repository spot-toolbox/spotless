function ltid_vw2abc_test(n,m,N)
% function ltid_vw2abc_test(n,m,N)
%
% testing ltid_vw2abc.m on random row G of order 2*n, m inputs, N samples

if nargin<1, n=5; end
if nargin<2, m=2; end
if nargin<3, N=100; end
G=ltid_rand(n,1,m);
w=linspace(0,pi,N)';
v=squeeze(freqresp(G,w)).';
v=v+0.01*(randn(size(v))+1i*randn(size(v)));
[aa,bb,cc,L]=ltid_vw2abc(v,w,2*n,2,0.1);
g=ltid_wabc2v(w,aa,bb,cc);
er=sqrt(max(sum(abs(g-v).^2,2)));
amin=ltid_tpmin(aa);
fprintf('\n error: %f>%f, passivity check %f>0\n',er,L,amin)
ltid_chk_vwg(v,w,g);