function [y,s]=mss_unique(x)
% function [y,s]=mss_unique(x)
%
% y=s*x is the matrix of unique rows in x, s is a sparse 0/1 matrix

[m,n]=size(x);
if (m<2)||(n==0), y=x;s=eye(m); return; end
[y,jj]=sortrows(x);
e=[1==1;any(y(1:m-1,:)~=y(2:m,:),2)];
ii=cumsum(e);
y=y(e,:);
k=size(y,1);
s=sparse(jj,ii,ones(m,1),m,k);