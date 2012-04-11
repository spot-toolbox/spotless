function ltid_vwa2b_test(n,m,N)
% function ltid_vwa2b_test(n,m,N)
%
% testing ltid_vwa2b.m on random row G of order 2*n, m inputs, N samples

if nargin<1, n=5; end
if nargin<2, m=2; end
if nargin<3, N=100; end
G=ltid_rand(n,1,m);
[p,q]=tfdata(G);
aa=ltid_q2a(q{1,1});
w=linspace(0,pi,N)';
v=squeeze(freqresp(G,w)).';
v=v+0.01*randn(size(v));
[bb,L]=ltid_vwa2b(v,w,aa,1);
g=ltid_abw2v(aa,bb,w);
er=sqrt(max(sum((g-real(v)).^2,2)));
amin=ltid_tpmin(aa);
fprintf('\n error: %f>%f, check %f>0\n',er,L,amin)
ltid_chk_vwg(v,w,g,1);