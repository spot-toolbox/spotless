function x=double(p)
% function x=double(p)
%
% x='?' unless p is constant

% AM 09.01.09

if size(p.s,2)==3,
    x=accumarray([p.s(:,1),p.s(:,2)],p.s(:,3),[p.m,p.n]);
else
    x='?';
end