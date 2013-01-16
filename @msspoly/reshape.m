function q=reshape(p,n,m)
if nargin == 2
    if ~spot_hasSize(n,[1 2])
        error('Single argument must be 1-by-2.');
    end
    if ~spot_isIntGE(n,1)
        error('Size argument must be positive integers.');
    end
    m = n(2);
    n = n(1);
end

N = prod(p.dim);

if isempty(n) || isempty(m)
    if isempty(n) && isempty(m)
        error('Only size argument can be empty.');
    elseif isempty(n) && spot_isIntGE(N/m,1)
        n = N/m;
    elseif isempty(m) && spot_isIntGE(N/n,1)
        m = N/n;
    else
        error('Number of elements not divisible by given size.');
    end
elseif N ~= n*m
    error('Number of elements must not change');
end

[qi,qj] = ind2sub([n m],sub2ind(p.dim,p.sub(:,1),p.sub(:,2)));

q = msspoly([n m],[qi qj],p.var,p.pow,p.coeff);
end
