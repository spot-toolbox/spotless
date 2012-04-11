function ltid_vwq2p_test(n,m,N)
% function ltid_vwq2p_test(n,m,N)
%
% testing ltid_vwq2p.m on random row G of order 2*n, m inputs, N samples

if nargin<1, n=5; end
if nargin<2, m=2; end
if nargin<3, N=100; end
G=ltid_rand(n,1,m);
A=ssdata(G);
q=poly(A);
w=linspace(0,pi,N)';
v=squeeze(freqresp(G,w)).';
v=v+0.1*randn(size(v));

p=ltid_vwq2p(v,w,q);
g=ltid_pqw2v(p,q,w);
er=sqrt(max(sum(abs(v-g).^2,2)));
fprintf(' Least Squares error: %f\n',er)
ltid_chk_vwg(v,w,g);

p=ltid_vwq2p(v,w,q,2);
g=ltid_pqw2v(p,q,w);
er=sqrt(max(sum(abs(v-g).^2,2)));
fprintf(' SeDuMi error: %f\n',er)
ltid_chk_vwg(v,w,g);
