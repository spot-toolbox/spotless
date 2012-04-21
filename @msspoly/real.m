function q = real(p)
    q = msspoly(p.dim,p.sub,p.var,p.pow,real(p.coeff));
end