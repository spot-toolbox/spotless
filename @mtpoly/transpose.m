function q=transpose(p)
q = msspoly(fliplr(p.dim),p.sub(:,[2 1]),p.var,p.pow,p.coeff);
end
