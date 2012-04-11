function [pr,U,Q,D]=bsos(pr0,q)
% function [pr,U,Q,D]=bsos(pr0,q)
% 
% INPUTS:
%   pr0 -  mssprog
%   q   -  msspoly
%
% OUTPUTS:
%   pr  -  updated pr0
%   U   -  a cell array of columns of msspoly monomials
%   Q   -  a cell array of semidefinite decision variables 
%   D   -  a cell array of non-negative vector decision variables
%
% registers with pr=pr0 new semidefinite decision variables Q{i}, 
% positive vector variables D{i}, and equalities 
% q(i)=sum(D{i})+U{i}'*(Q{i}-diag(D{i}))*U{i} 

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
    pp0=p(:,b==0);
    pp=mint_ch([zeros(1,size(pp0,2));pp0;pp0+1],2);
    U{i}=recomp(yy,[pp]);
    ni=size(U{i},1);
    [pr,QQ]=new(pr,ni,'psd');
    Q{i}=QQ;
    [pr,DD]=new(pr,ni,'pos');
    D{i}=DD;
    pr=eq(pr,q(i)-sum(DD)-U{i}'*(QQ-diag(DD))*U{i});
end
    
