function q = real(p)
    q = (p + conj(p))/2;
    %    q = msspoly(p.dim,p.sub,p.var,p.pow,real(p.coeff));
end