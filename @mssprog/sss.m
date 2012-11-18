function [pr,U,Q]=sss(pr0,q)
% function [pr,U,Q]=sss(pr0,q)
% 
% INPUTS:
%   pr0 -  mssprog
%   q   -  msspoly (square sized)
%
% OUTPUTS:
%   pr  -  updated pr0
%   U   -  a column of msspoly monomials
%   Q   -  a semidefinite decision variable 
%
% registers new semidefinite decision variable Q, 
% as well as the equalities v'*q*v=U'*Q*U, to pr0 (updated to pr),
% where v=msspoly('#',size(q,1)) 

if nargin<2, error('2 inputs required'); end
if ~isa(q,'msspoly'), error('input 2 not an "msspoly"'); end
if isempty(q), error('input 2 is empty'); end
pr=pr0;
[m,n]=size(q);
if m~=n, error('input 2 not a square matrix'); end
v=msspoly('#',m-1);
[pr,U0,Q0]=sos(pr,[1;v]'*q*[1;v]);
U=U0{1};
Q=Q0{1};