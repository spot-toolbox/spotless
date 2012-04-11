function q=mpower(p,n)
% function q=mpower(p,n)
% 
% p must be square, n converted to a non-negative integer

% AM 09.01.09

if p.n~=p.m, error('non-square matrix'); end
if ~isa(n,'double'), error('2nd input not "double"'); end
if ~isequal(size(n),[1 1]), error('2nd input not 1-by-1'); end
if n~=max(0,round(double(n))), error('2nd input not integer>=0'); end
q=eye(p.n);
if n==0, q=msspoly(speye(p.n)); return; end
if n==1, q=p; return; end
if n==2, q=p*p; return; end
if n==3, q=p*p*p; return; end
k=floor(log2(n))+1;
h{1}=p; for i=2:k, h{i}=h{i-1}*h{i-1}; end
q=p;
n=n-1;
while n>=1,
    r=floor(log2(n))+1;
    q=q*h{r};
    n=n-2^(r-1);
end

