function [pr,U,Q]=sos(pr0,q)
% function [pr,U,Q]=sos(pr0,q)
% 
% INPUTS:
%   pr0 -  mssprog
%   q   -  msspoly
%
% OUTPUTS:
%   pr  -  updated pr0
%   U   -  a cell array of columns of msspoly monomials
%   Q   -  a cell array of semidefinite decision variables 
%
% registers new semidefinite decision variables Q{i}, as well as
% the equalities q(i)=U{i}'*Q{i}*U{i} with pr (updated to pr1)

if nargin<2, error('2 inputs required'); end
if ~isa(q,'msspoly'), error('input 2 not an "msspoly"'); end
if isempty(q), error('input 2 is empty'); end
pr=pr0;
xx=variables(pr);
for i=1:length(q),
    [x,p]=decomp(q(i));
    b=match(xx,x);
    if any(any(p(:,b>0)>1)),
        error(['inequality no. ' num2str(i) ' is not linear in decision parameters']);
    end
    yy=x(b==0);
    pp=p(:,b==0);
    pp=mint_ch(pp,2);
    U{i}=recomp(yy,pp);
    [pr,Q{i}]=new(pr,size(U{i},1),'psd');
    pr=eq(pr,q(i)-U{i}'*Q{i}*U{i});
end
    
