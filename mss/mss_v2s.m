function y=mss_v2s(x,z)
% function y=mss_v2s(x,z)
%
% with nargin==1, re-arrange the rows of x into block symmetric matrix 
%     [ x(1) x(2) x(4)    ]
%     [ x(2) x(3) x(5)    ]
% y = [ x(4) x(5) x(6)    ]
%     [               ... ] 
% with nargin>1, use the scheme 
%                                 [  4  5  6  7  ]
%          [ 2 3 4 ]              [  5  1  2  8  ]
% (1:6) -> [ 3 1 5 ],   (1:10) -> [  6  2  3  9  ],   etc.
%          [ 4 5 6 ]              [  7  8  9  10 ]
%

% AM 10.01.09

x=x(:);
n=size(x,1);
m=round((sqrt(1+8*n)-1)/2);
if m*(m+1)~=2*n, error('unable to convert: wrong dimension'); end

ii=repmat(1:m,m,1);         % column number in an m-by-m matrix
jj=repmat((1:m)',1,m);      % row number
a=min(ii,jj);
b=max(ii,jj);
if nargin<2,
    h=a+round(b.*(b-1)/2);
else
    e=2*min(a,m-b+1)-1;
    c=m-e;
    d=b+a-e;
    h=d+round(c.*(c-1)/2);
end
    
y=x(h);
