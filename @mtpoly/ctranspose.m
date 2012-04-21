function q=ctranspose(p)
q = msspoly(fliplr(p.dim),p.sub(:,[2 1]),p.var,p.pow,conj(p.coeff));
end
