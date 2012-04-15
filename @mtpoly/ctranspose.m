function q=ctranspose(p)
% function q=ctranspose(p)
%
% AM 09.01.09

q=mtpoly(p.n,p.m,[p.s(:,2) p.s(:,1) p.s(:,3:size(p.s,2))]);
