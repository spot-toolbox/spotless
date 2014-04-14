function [d,err] = spot_psdNoToDim(n)
    d=round((sqrt(1+8*n)-1)/2);
    if spot_psdDimToNo(d) ~= n
        d = NaN;
        err = 1;
    else
        err = 0;
    end
end