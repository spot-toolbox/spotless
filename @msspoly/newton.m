function x=newton(p,v,k,x0,q)
% function x=newton(p,v,k,x0,q)
%
% for mss polynomial p in free msspoly variable v,
% performs k (default k=1) iterations x -> x-inv(q(x))*p(x),
% (default q=diff(p)) starting at x=x0 (default x0=0)

if nargin<2, error('2 inputs required'); end
[mp,np]=size(p);
if np~=1, error('1st input must be a column vector'); end
if ~isa(v,'msspoly'), error('2nd input not msspoly'); end
if ~isfree(v), error('2nd input must be free'); end
if size(v,2)~=1, error('input 2 must be a column'); end
if mp~=length(v), error('inputs 1,2: incompatible dimensions'); end
if ~isfunction(p,v), error('1st input not a function of 2nd'); end
if nargin<3, k=1; end; 
if ~isa(k,'double'), error('3rd input not a "double"'); end
if ~isscalar(k), error('3rd input not a scalar'); end
if k~=max(1,round(k)), error('3rd input not integer>0'); end
if nargin<4, x0=zeros(mp,1); end
if ~isa(x0,'double'), error('4th input not a "double"'); end
if ~isequal(size(x0),[mp,1]), error('inputs 1,4: incompatible dimensions'); end
if nargin<5, q=diff(p,v); end
if ~isa(q,'msspoly'), error('5th input not msspoly'); end
if ~isequal(size(q),[mp,mp]), error('inputs 1,5: incompatible dimensions'); end
if ~isfunction(q,v), error('5th input not a function of 2nd'); end

x=x0;
for i=1:k,
    x=x-inv(double(subs(q,v,x)))*double(subs(p,v,x));
end