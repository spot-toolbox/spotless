function q=recomp(x,p,M)
% function q=recomp(x,p,M)
%
% inverse to @msspoly/decomp.m
 
if nargin<2, error('2 inputs required'); end
[b,xn]=isfree(x);
if ~b, error('input 1 is not free'); end
if size(xn,2)~=1, error('input 2 must be a column'); end
if ~mint_isint(p), error('input 2 not integer'); end
if any(p(:)<0), error('input 2 has negative entries'); end
nx=size(x,1);
[mp,np]=size(p);
if np~=nx, error('inputs 1,2 incompatible'); end
q=msspoly(mp,1,[(1:mp)' ones(mp,1) (p>0)*diag(xn) p ones(mp,1)]);
if nargin>2,
    if ~isa(M,'double'), error('input 3 not a "double"'); end
    if size(M,2)~=mp, error('inputs 2,3 incompatible'); end
    q=M*q;
end
    