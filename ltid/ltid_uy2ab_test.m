function as=ltid_uy2ab_test(m,k,n)

if nargin<1, m=3; end
if nargin<2, k=4; end
if nargin<3, n=200; end
if nargin<4, st0=1; end

STREAM = RandStream.getDefaultStream;
reset(STREAM);


a=[0;randn(m-1,1)];
a(1)=-ltid_tpmin(a);
a=a/a(1);
b=randn(k,1);
u=randn(n,k);
y=randn(n,m-1);
y=[u*b-y*a(2:m)+0.0*randn(n,1) y];

[as,bs]=ltid_uy2ab(u,y);

if nargout<1,
    fprintf(' relative error(L2):   %f\n', ...
        (norm(a-as)+norm(b-bs))/(norm(a)+norm(b)))
    [as,bs]=ltid_uy2ab(u,y,'LInf');
    fprintf(' relative error(LInf): %f\n', ...
        (norm(a-as)+norm(b-bs))/(norm(a)+norm(b)))
end