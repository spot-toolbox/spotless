function q = fliplr(p)
    sz = size(p);
    if sz(2) < 2 || sz(1) < 1
        q = p;
    else
        [i,j,s] = find(p);
        j = sz(2) - j + 1;
        I = sub2ind(size(p),i,j);
        q = msspoly(zeros(sz));
        q=assign(q,s,I);
    end
end