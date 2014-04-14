function flg = spot_hasSize(v,sz)
    if length(size(v)) ~= length(sz)
        flg = 0;
    else
        flg = isequal(size(v),sz);
    end
end
