function [x,err]=double(p)
% x = double(p);
%
% p -- an msspoly
%
% Throws an error if p is not a constant matrix, o.w. x is 
% the corresponding matrix.
%
% [x,err] = double(p);
%
% Instead of throwing an error the flag err == 1 if p is not
% a constant matrix.
if size(p.var,2) == 0 || size(p.var,1) == 0
    x=accumarray(p.sub,p.coeff,p.dim);
    err = 0;
else
    if nargout == 2
        err = 1;
        x = [];
    else
        error('Cannot cast non-constant msspoly to double.');
    end
end
