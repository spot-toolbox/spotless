function flg = hasSize(v,sz)
    if length(size(v)) ~= length(sz)
        flg = 0;
    else
        flg = all(size(v) == sz);
    end
end
