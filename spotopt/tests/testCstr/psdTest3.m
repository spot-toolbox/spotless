function [pr,obj,zero,tol] = psdTest3(pr)
%
%  Evalutate the nuclear norm of a matrix, posed
%  in dual form.
%
    [pr,f] = pr.newFree(1);
    [pr,z] = pr.withPSD([1 f ; f 1]);
    
    obj = f;
    
    zero = z(1,2) - (-1);
    tol  = 5e-8;
end