function p = test_randomPoly(n,dim,dens,int)
    if nargin < 3
        dens = 4;
    end
    if nargin < 4
        int = 1;
    end
    x = msspoly('x',n);
    [~,xn] = isfree(x);
    % Pick random points in the matrix.
    K = floor(prod(dim)*dens);
    sub = ceil(repmat(dim,K,1).*rand(K,2));
    
    var = repmat(xn',K,1);
    pow = poissrnd(1,K,n);
    if int
        coeff = poissrnd(2,K,1);
    else
        coeff = randn(K,1);
    end
    
    p = msspoly(dim,sub,var,pow,coeff);
end