function q = imag(p)
    q = msspoly(p.dim,p.sub,p.var,p.pow,imag(p.coeff));
end