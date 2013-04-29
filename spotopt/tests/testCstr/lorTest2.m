function [pr,obj,zero,tol] = lorTest2(pr)
%
%  Test:  A single column, posed in primal form.
    
    data = randn(10,1);
    [pr,l] = pr.newLor(length(data)+1);
    pr = pr.withEqs(l(1+(1:length(data))) - data);
    
    obj = l(1);
    zero = l(1) - norm(data);
    tol  = 1e-8*norm(data);
end