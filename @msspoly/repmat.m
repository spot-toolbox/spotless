function q=repmat(p,rep,m)

if nargin == 3
    rep = [ rep m ];
end

if ~spot_hasSize(rep,[1 2])
    error('Repetition for exactly two dimensions must be given.');
elseif ~spot_isIntGE(rep,1)
    error('Arguments must be positive integers.');
end

sz  = p.dim;
szq = rep.*sz;

if size(p.coeff,1) == 0
    q = msspoly.zeros(szq);
else
    E = size(p.coeff,1);
    vs = repmat(p.var,  prod(rep),1);
    ps = repmat(p.pow,  prod(rep),1);
    cs = repmat(p.coeff,prod(rep),1);
    ss = repmat(p.sub, prod(rep),1);
    
    ss(:,1) = ss(:,1) + kron((1:rep(1))'-1,repmat(p.dim(1),E*rep(2),1));
    ss(:,2) = ss(:,2) + repmat(kron((0:rep(2)-1)',...
                                    repmat(p.dim(2),E,1)),rep(1),1);

    q = msspoly(szq,ss,vs,ps,cs);
end

end


