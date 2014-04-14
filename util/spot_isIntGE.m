function flg = spot_isIntGE(var,bnd)
    if nargin < 2, bnd = -Inf; end
    if ~isa(var,'double')
        flg = 0;
    elseif ~isequal(fix(var),var) | any(var < bnd)
        flg = 0;
    else
        flg = 1;
    end
end