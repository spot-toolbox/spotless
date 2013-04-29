function [pr,obj,zero,tol] = lorTest3(pr)
%
%  Test:  Multiple columns, posed in dual form.
    
    n= 10;
    k = 30;
    data = randn(n,k);
    [pr,l] = pr.newFree(k);
    pr = pr.withLor([l';data]);
    
    zero = sum(l) - sum(sqrt(sum(data.^2)));
    tol  = 1e-8*sum(sqrt(sum(data.^2)));
    obj = sum(l);
end