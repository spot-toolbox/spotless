function msk = istrig(q)
%
%  b = istrig(q)
%
%  q -- n-by-m free msspoly.
%  b -- n-by-m array, 1 if q(i,j) is a trigonometric variable, 0 o.w.
    if isempty(q),
        msk = [];
        return;
    end
    [x,p,C] = decomp(q);
    [flag,xid] = isfree(x);

    if ~flag
        error('Argument must be free');
    end
    
    msk = logical(double(subs(q,x,double(msspoly.isTrigId(xid)))));
end