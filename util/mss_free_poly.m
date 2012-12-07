function [prog,poly,coeff] = mss_free_poly(prog,basis,n)
    if nargin < 3, n = 1; end
    [prog,coeff] = new(prog,length(basis)*n,'free');
    poly = reshape(coeff,n,length(basis))*basis;
end