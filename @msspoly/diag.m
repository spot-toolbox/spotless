function q=diag(p)
M = max(p.dim);
m = min(p.dim);

if m == 1 % Construct a diagonal matrix from this.
    sub = p.sub;
    if p.dim(1) > p.dim(2),
        sub(:,2) = sub(:,1);
    else
        sub(:,1) = sub(:,2);
    end
    q = msspoly([M M],...
               sub,p.var,p.pow,p.coeff);
else % Extract the diagonal from this.
    msk = p.sub(:,1) == p.sub(:,2);
    q = msspoly([m 1],...
               [ p.sub(msk,1) ones(nnz(msk),1) ],...
               p.var(msk,:),p.pow(msk,:),...
               p.coeff(msk,:));
end

end
