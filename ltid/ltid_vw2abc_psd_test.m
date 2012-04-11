function ltid_vw2abc_psd_test(n,k,N)
% function ltid_vw2abc_psd_test(n,k,N)
%
% testing ltid_vw2abc_psd.m on a randomly generated k-by-k
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

[aa,bb,cc,L]=ltid_vw2abc_psd(v,w,2*n,2);
amin=ltid_tpmin(aa);
t=[w;linspace(2.001,pi,N)'];
cst=cos(t*(0:2*n));
snw=sin(w*(1:2*n));
csw=cos(w*(0:2*n));
bt=cst*bb;
g=(csw*bb+1i*(snw*cc))./repmat(csw*aa,1,k^2);
ix=mss_s2v(reshape(1:k^2,k,k));
er=sqrt(max(sum(abs(v(:,ix)-g(:,ix)).^2,2)));
r=zeros(2*N,1);
for i=1:2*N, r(i)=min(eig(reshape(bt(i,:),k,k))); end
fprintf('\n error: %f>%f, check min(a)=%f>0,  min(eig(b))=%f>0\n', ...
        er,L,amin,min(r))

close(gcf);plot(t,r); grid;pause
ltid_chk_vwg(v(:,ix),w,g(:,ix),0);

