function [b,xn]=isfree(p)
% [b,xn]=isfree(p)
%
% Returns b=1 if p is an msspoly whose entries are each unique
%           variables to the power 1 with coefficient 1.
%
%         if b = 1, xn is the matrix of variable ids in the order
%         they appear in p(:).
%
xn = [];


if ~isa(p,'msspoly'), 
    b = 0;
    return;
end

if isempty(p), b = 1;
elseif size(p.pow,1) == 0, b = 0;      % empty matrix not free.
elseif size(p.pow,2) ~= 1, b = 0; % more than one power.
elseif any(p.pow ~= 1), b = 0;    % powers other than one.
elseif all(p.dim == [1 1])          % scalar case
    b = 1;
    xn = p.var(1);
else
    if size(p.pow,1) ~= prod(p.dim) % some zero variables.
        b = 0;
    elseif any(p.var == 0)          % some constant terms.
        b = 0;
    elseif length(p.var) ~= length(unique(p.var))
        b = 0;
    else
        b = 1;
        xn = zeros(p.dim);
        xn(sub2ind(p.dim,p.sub(:,1),p.sub(:,2))) = p.var;
    end
end
