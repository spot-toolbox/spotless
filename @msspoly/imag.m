function q = imag(p)
    q = (p - conj(p))/(2*j);
    %    q = msspoly(p.dim,p.sub,p.var,p.pow,imag(p.coeff));
end