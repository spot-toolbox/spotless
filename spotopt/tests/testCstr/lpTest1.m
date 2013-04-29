function [pr,obj,zero,tol] = lpTest1(pr)
%
%  Evalutate the nuclear norm of a matrix, posed
%  in dual form.
%
    n = 10;
    m = 4;

    data = randn(n,m);
    
    [pr,t] = pr.newFree(1);
    [pr,s] = pr.newFree(n*m);
    
    s = reshape(s,n,m);
    
    pr = pr.withPos(s - data);
    pr = pr.withPos(s + data);
    pr = pr.withPos(t - sum(s));
    
    obj = t;
    zero = obj - max(sum(abs(data)));
    tol  = 1e-7*max(sum(abs(data)));
end