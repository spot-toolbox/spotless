function [pr,l,zero,tol] = lorTest1(pr)
%
%  Test:  A single column, posed in dual form.
    
    data = randn(10,1);
    [pr,l] = pr.newFree(1);
    pr = pr.withLor([l;data]);
    
    zero = l - norm(data);
    tol  = 1e-8*norm(data);
end