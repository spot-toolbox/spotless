function ltid_a2q_test(m)

if nargin<1, m=5; end
q=[1;randn(m,1)];
z=roots(q);
q=rand*poly(z/(0.01+max(abs(z))))';
qi=q(m+1:-1:1);
a=conv(q,qi);
a=[a(m+1);2*a(m+2:2*m+1)];
qq=ltid_a2q(a)';

fprintf('\n ltid_bq2G relative error: %f\n',norm(qq-q)/norm(q))