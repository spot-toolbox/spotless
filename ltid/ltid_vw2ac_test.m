function ltid_vw2ac_test(n,m,N)
% function ltid_vw2ac_test(n,m,N)
%
% testing ltid_vw2ac.m on random row G of order 2*n, m inputs, N samples

if nargin<1, n=5; end
if nargin<2, m=2; end
if nargin<3, N=100; end
G=ltid_rand(n,1,m);
w=linspace(0.5,1.5,N)';
z=exp(1i*w);
p=w+0.1;
v=squeeze(freqresp(G,w)).';
v=v+(0.01i)*randn(size(v));

[aa,cc,L]=ltid_vw2ac(v,w,2*n,1,1e-1,p);
g=imag(ltid_wabc2v(w,aa,[],cc));
er=max(sqrt(sum((g-imag(v)).^2,2))./p);
amin=ltid_tpmin(aa);
fprintf('\n Algorithm 1 error: %f < L < %f, check %f>0\n',L,er,amin)
ltid_chk_vwg(imag(v),w,g,1);

[aa,cc,L]=ltid_vw2ac(v,w,2*n,2,1e-1,p);
g=imag(ltid_wabc2v(w,aa,[],cc));
er=max(sqrt(sum((g-imag(v)).^2,2))./p);
amin=ltid_tpmin(aa);
fprintf('\n Algorithm 2 error: %f < L < %f, check %f>0\n',L,er,amin)
ltid_chk_vwg(imag(v),w,g,1);
