function q=trace(p)
if p.dim(1) ~= p.dim(2)
    error('Matrix must be square.');
end

msk = p.sub(:,1) == p.sub(:,2);

if isempty(msk)
    q = msspoly(0);
else
    q = msspoly([1 1],repmat([1 1],nnz(msk),1),...
               p.var(msk,:),p.pow(msk,:),...
               p.coeff(msk,:));
end
end
