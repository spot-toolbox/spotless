function [pr,obj,zero,tol] = lpTest1(pr)
%
%  Evalutate the nuclear norm of a matrix, posed
%  in dual form.
%
    n = 10;
    m = 4;

    data = randn(n,m);
    
    
    [pr,s] = pr.newPos(n*m);
    [pr,u] = pr.newPos(n*m);
    [pr,l] = pr.newPos(n*m);
    [pr,v] = pr.newPos(m);
    [pr,t] = pr.newPos(1);
    
    u = reshape(u,n,m);
    l = reshape(l,n,m);
    s = reshape(s,n,m);
    
    pr = pr.withEqs((s-data) - u);
    pr = pr.withEqs((s+data) - l);
    pr = pr.withEqs((t - sum(s))' - v);
    
    obj = t;
    zero = obj - max(sum(abs(data)));
    tol  = 1e-7*max(sum(abs(data)));
end