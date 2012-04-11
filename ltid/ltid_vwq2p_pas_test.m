function ltid_vwq2p_pas_test(n,k,N)
% function ltid_vwq2p_pas_test(n,k,N)
%
% testing ltid_vwq2p_pas.m on a randomly generated k-by-k
% passive symmetric transfer matrix G of order 2*n, with N samples

if nargin<1, n=3; end
if nargin<2, k=2; end
if nargin<3, N=100; end

STREAM = RandStream.getDefaultStream;
reset(STREAM);

G=ltid_rand_pas(n,k);
G=G+G.';

H=ltid_reshape(G,1,k^2);
w=linspace(0,2,N)';
v=squeeze(freqresp(H,w)).';
v=v+0.01*(randn(size(v))+1i*randn(size(v)));

m=2*n;
aa=ltid_vw2abc_psd(v,w,m,1);
q=ltid_a2q(aa);
t=ltid_qw2t(q,w,5);
p=ltid_vwq2p_pas(v,w,q,t,2);

t=[w;linspace(2.001,pi,N)'];
zt=exp(t*(1i*(m:-1:0)));
zw=exp(w*(1i*(m:-1:0)));
gt=(zt*p)./repmat(zt*q',1,k^2);
g=(zw*p)./repmat(zw*q',1,k^2);
ix=mss_s2v(reshape(1:k^2,k,k));
er=sqrt(max(sum(abs(v(:,ix)-g(:,ix)).^2,2)));
r=zeros(2*N,1);
for i=1:2*N, r(i)=min(eig(real(reshape(gt(i,:),k,k)))); end
fprintf('\n error: %f,  passivity check %f>0\n', ...
        er,min(r))

close(gcf);plot(t,r); grid;pause
ltid_chk_vwg(v(:,ix),w,g(:,ix),0);

