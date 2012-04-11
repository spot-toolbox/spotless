function q=uminus(p)
%

% AM 09.01.09

n=size(p.s,2);
q=msspoly(p.m,p.n,[p.s(:,1:n-1) -p.s(:,n)]);