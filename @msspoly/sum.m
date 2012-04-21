function q=sum(p,d)

if isempty(p), q = p;
else
    if nargin<2, 
        if size(p,1)==1,
            d=2;
        else
            d=1;
        end
    end
    dim = p.dim;
    dim(d) = 1;
    sub = p.sub;
    sub(:,d) = 1;
    
    q = msspoly(dim,sub,p.var,p.pow,p.coeff);
end

