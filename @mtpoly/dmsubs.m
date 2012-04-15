function q=msubs(p,x,v)
% function q=msubs(p,x,v)
%
% INPUTS:
%   p  -  m-by-1 mtpoly 
%   x  -  k-by-1 free mtpoly 
%   v  -  k-by-n real double 
%
% OUTPUT:
%   q  -  m-by-n double
%
% DESCRIPTION: q(:,i) is the result of substituting v(:,i) for x in p
% p must be an mtpoly in x alone.

% mmt 3.24.11
if nargin<3, error('three inputs required'); end
if ~isa(x,'mtpoly'), error('input 2 not a mtpoly'); end
[f,xn]=isfree(x);
if ~f, error('input 2 is not free'); end
if ~isa(v,'double'), error('input 3 not a double'); end
if ~isreal(v), error('input 3 not real'); end
[k,n]=size(v);
if ~isequal(size(xn),[k 1]), error('inputs 2,3 not compatible'); end
[m,np]=size(p);
if np~=1, error('input 1 is not a column'); end

[xd,pd,Md] = decomp(p);

N = size(v,2);
po = size(pd,1);
n = size(pd,2);

% First, test that xd is a subset of xn.
% Sort out indicies.
[~,xdn] = isfree(xd);
perm = mss_match(xn,xdn);

perm = perm(perm ~= 0);
% Second, generate matrix of evaluated monomials.
if length(perm) ~= length(xd)
    error('p must only be a function of x');
end

pd = pd';
pd = pd(:);
pow = repmat(v(perm,:),po,1).^repmat(pd(:),1,N); 

sz = [po N];
subs = repmat(1:prod(sz),n,1);

mo = accumarray(subs(:),pow(:),[prod(sz) 1],@prod);
q = Md*reshape(mo,sz);




