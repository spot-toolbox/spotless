function [x,y]=size(p,n)
if nargin == 1
    if nargout == 2
        x = p.dim(1);
        y = p.dim(2);
    else
        x = p.dim;
    end
else
    x = p.dim(n);
end

end
