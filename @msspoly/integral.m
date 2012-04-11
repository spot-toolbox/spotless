function q=integral(p,y)
% function q=integral(p,y)
% 
% INPUTS:
%   p  -  msspoly
%   y  -  free scalar msspoly
%
% OUTPUT:
%   q  -  msspoly
%
% DESCRIPTION:
%   q(y,x)=int_0^y p(y,x)dy
%
% 18.09.09 by ameg@mit.edu

if ~isa(y,'msspoly'), error('input 2 not a msspoly'); end
[b,xn]=isfree(y);
if ~b, error('input 2 not free'); end
if size(xn,1)~=1, error('input 2 not scalar'); end

[ms,ns]=size(p.s);
k=round((ns-3)/2);
b=(p.s(:,3:2+k)==xn);       % entries with id(y)
bb=any(b,2);                % rows with id(y)
bi=~bb;                     % rows without id(y)
be=zeros(ms,1);             % column with xn in matching places
be(bi)=xn;
sd=p.s(:,3+k:2+2*k);        % degrees matrix
sc=p.s(:,ns);               % coefficients column
is=find(b);                 % linear indexes for id(y) matches
isr=mod(is-1,ms)+1;         % rows of id matches
sd(is)=sd(is)+1;            % increase matching degrees by 1
sc(isr)=sc(isr)./sd(is);    % divide matching coefficients by degree
q=msspoly(p.m,p.n,[p.s(:,1:2+k) be sd bi sc]);