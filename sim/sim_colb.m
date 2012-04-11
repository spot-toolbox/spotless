function v=sim_colb(y,d)
% function v=sim_colb(y,d)
%
% convex lower bound of data vector y using depth d 
%
% INPUTS:
%   y  -  n-by-k real
%   d  -  positive integer 
% OUTPUT:
%   v  -  convex lower bound of depth d

v=y;
n=size(y,1);
for i=1:d,
    for j=1:d,
        s=i+j;
        v(i+1:n-j,:)=min(v(i+1:n-j,:),(j/s)*v(1:n-s,:)+(i/s)*v(s+1:n,:));
    end
end