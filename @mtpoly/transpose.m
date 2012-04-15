function q=transpose(p)
% function q=transpose(p)
%
% AM 09.01.09

q=mtpoly(p.n,p.m,[p.s(:,2) p.s(:,1) p.s(:,3:size(p.s,2))]);
