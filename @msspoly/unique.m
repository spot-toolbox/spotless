function q = unique(p)
    [x,pow,Coeff] = decomp(p);
    Coeff = unique(Coeff,'rows');
    q = recomp(x,pow,Coeff);
end