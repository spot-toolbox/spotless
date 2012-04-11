function nlid_tpls_test(n,d)
% function nlid_tpls_test(n,d)
%
% test nlid_tpls.m on n samples of y=min(|x|,1-|x|)

if nargin<1, n=10000; end
if nargin<2, d=[4 4 4]; end

x=randn(n,length(d));
y=sqrt(sum(x.*x,2));
y=min(y,1-y);
v=nlid_tpls(x,y,d);
fprintf(' relative errors: L2=%f, LInf=%f\n', ...
    norm(y-v)/norm(y),max(abs(y-v))/max(abs(y)))
close(gcf);plot(y,v,'.');grid