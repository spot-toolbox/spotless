function q=uminus(p)
q = msspoly(p.dim,p.sub,p.var,p.pow,-p.coeff);
end
