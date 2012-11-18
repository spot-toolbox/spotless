function q=conj(p)
msk = msspoly.isTrigId(p.var);
pow = p.pow;
pow(msk) = -1*pow(msk);
q = msspoly(p.dim,p.sub,p.var,pow,conj(p.coeff));
end
