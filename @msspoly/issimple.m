function [s,x]=issimple(p)
%
%
%  [s,x] = issimple(p)
%
%  s -- 1 when p is a k-by-1 msspoly of constants and
%       variables. 0 otherwise.
%
%  x -- k-by-2 double array
%       x(:,1) are the variable ID numbers (0 for constant).
%       x(:,2) is the vector of values (1 for variables).
%
%
x = [];
if ~isa(p,'msspoly') || isempty(p) || size(p,2) ~= 1
    s = 0;
else
    ind = sub2ind(p.dim,p.sub(:,1),p.sub(:,2));
    if length(unique(ind)) ~= length(ind)
        s = 0;
    else
        if isempty(p.var), msk = logical(zeros(size(p.coeff))); 
        else
            msk = p.var ~= 0; 
        end
        
        if size(p.var,2) > 1 | ...
                ~all(p.pow == 1) |...
                ~all(p.coeff(msk) == 1)
            s = 0;
        else
            x = zeros(size(p,1),2);
            x(ind(msk),1) = p.var(msk);
            x(ind(~msk),2) = p.coeff(~msk);
            s = 1;
        end
    end
end
