function [pr,obj,zero,tol] = psdTest2(pr)
%
%  Evalutate the nuclear norm of a matrix, posed
%  in dual form.
%
    n = 10;
    m = 4;

    A = randn(n,m);
    
    [pr,P] = pr.newPSD(n+m);
    
    W1 = P(1:n,1:n);
    W2 = P(n+(1:m),n+(1:m));
    
    pr = pr.withEqs(P(1:n,n+(1:m)) - A);
    
    obj = (trace(W1)+trace(W2))/2;
    zero = obj - sum(svd(A));
    tol  = 1e-7*sum(svd(A));
end