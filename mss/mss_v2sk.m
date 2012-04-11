function y=mss_v2sk(x)
% function y=mss_v2sk(x)
%
% re-arrange vector x into skew symmetric matrix 
%     [   0     x(1)   x(2)   x(4)    ]
%     [ -x(1)     0    x(3)   x(5)    ]
% y = [ -x(2)  -x(3)     0    x(6)    ]
%     [ -x(4)  -x(5)  -x(6)     0     ]
%     [                           ... ]
%     

% AM 27.06.09

x=x(:);
n=size(x,1);
x=[x;0];
m=round((sqrt(1+8*n)-1)/2);
if m*(m+1)~=2*n, error('unable to convert: wrong dimension'); end
ii=repmat(1:m,m,1);
jj=repmat((1:m)',1,m);
a=min(ii,jj); 
b=max(ii,jj); 
h=a+round(b.*(b-1)/2);
h(ii>jj)=n+1;
y=x(h);
y=[zeros(1,m+1);y zeros(m,1)];
y=y-y.';
