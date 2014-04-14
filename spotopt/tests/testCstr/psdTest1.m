function [pr,obj,zero,tol] = psdTest1(pr)
%
%  Evalutate the nuclear norm of a matrix, posed
%  in dual form.
%
    n = 10;
    m = 4;

    A = randn(n,m);
    
    [pr,w1] = pr.newFree(nchoosek(n+1,2));
    W1 = mss_v2s(w1);
    [pr,w2] = pr.newFree(nchoosek(m+1,2));
    W2 = mss_v2s(w2);
    
    pr = pr.withPSD([ W1 A ; A' W2]);
    
    obj = (trace(W1)+trace(W2))/2;
    zero = obj - sum(svd(A));
    tol  = 1e-7*sum(svd(A));
end