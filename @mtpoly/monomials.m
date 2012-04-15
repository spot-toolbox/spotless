function q=monomials(x,p)
% function q=monomials(x,p)
%
% x is a free mss polynomial, 
% q is the column of all monomials of x of degrees specified in p

if nargin<2, error('2 inputs required'); end
[f,xn]=isfree(x);
nx=length(xn);
if ~f, error('1st input not a free mtpoly'); end
if size(xn,2)~=1, error('input 1 must be a column'); end
if ~isa(p,'double'), error('2nd input not a double'); end
if ~isequal(p,max(0,round(p))), error('2nd input not integer>=0'); end
p1=sort(p(:));
np=length(p1);
if any([-1;p1(1:np-1)]==p1), error('2nd input has repeated entries'); end
dd=mss_asd(nx,p);
mq=size(dd,1);
q=mtpoly(mq,1,[(1:mq)' ones(mq,1) repmat(xn',mq,1) dd ones(mq,1)]);
