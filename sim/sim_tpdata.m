function [D,R]=sim_tpdata(d)
% function [D,R]=sim_tpdata(t,d)
%
% INPUT:
%   d  -  1-by-n positive integer
% OUTPUTS:
%   D  -  k-by-n integer
%   R  -  k-by-1 real
%
% D lists all integer rows h such that -d-1<h<d+1, so that the functions
% fp: Rn->R defined by fp(t)=p'*cos(D*t+R) for n-by-1 real t,
% span all trigonometric polynomials of degree bounded by d, as p
% ranges over the set of all k-by-1 real columns
%
% NOTE:  y=cos(s*D'+repmat(R',size(s,1),1))*p' produces the column of
%        the values fp takes on the rows of real m-by-n matrix s

d=max(1,round(d(:)))';
n=size(d,2);
dd=2*d+1;
k=prod(dd);
D=zeros(k,n);
z=(0:(k-1))';
M=0;
for i=1:n,
    M=d(n+1-i)+dd(n+1-i)*M;
    D(:,i)=mod(z,dd(i));
    z=(z-D(:,i))/dd(i);
    D(:,i)=D(:,i)-d(i);
end
R=zeros(k,1);
R((0:(k-1))'>M)=pi/2;