function ltid_ab2pq_test(m,k)

if nargin<1, m=5; end
if nargin<2, k=7; end
aa=randn(m+1,1);
aa(1)=aa(1)-ltid_tpmin(aa)+0.1;
bb=randn(m+1,k);
[p,q]=ltid_ab2pq(aa,bb);
t=linspace(0,pi,50)';
z=exp(1i*t*(0:m));
cs=real(z);
vpq=real((z*p.')./repmat((z*q.'),1,k));
vab=(cs*bb)./repmat(cs*aa,1,k);
fprintf('\n ltid_ab2pq error: %f out of %f\n',norm(vpq-vab),norm(vab))