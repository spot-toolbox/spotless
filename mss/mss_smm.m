function y=mss_smm(x)
% function y=mss_smm(x)
%
% symmetric mismatch for square matrices
[n,m]=size(x);
if n~=m, error('argument must be a square matrix'); end
y=mss_s2v(x(1:n-1,2:n));
x=x';
y=y-mss_s2v(x(1:n-1,2:n));
