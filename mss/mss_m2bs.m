function y=mss_m2bs(x)
% function y=mss_m2bs(x)
%
% re-arrange the rows of x into block symmetric matrix 
%     [ x(1,:) x(2,:) x(4,:)    ]
%     [ x(2,:) x(3,:) x(5,:)    ]
% y = [ x(4,:) x(5,:) x(6,:)    ]
%     [                     ... ] 

% AM  24.02.10

[n,k]=size(x);
m=round((sqrt(1+8*n)-1)/2);
if m*(m+1)~=2*n, error('unable to convert: wrong dimension'); end
ii=repmat(1:m,m,1);
jj=repmat((1:m)',1,m);
a=min(ii,jj);
b=max(ii,jj);
h=a+round(b.*(b-1)/2);
y=x(kron(h,ones(1,k))+repmat(n*(0:k-1),m,m));
