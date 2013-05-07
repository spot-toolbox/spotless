function [x,y] = spot_eq_solve(A,b,tol)
% 
% [x,y] = spot_eq_feasible(A,b)
% [x,y] = spot_eq_feasible(A,b,tol)
%
%  A    -- m-by-n matrix, m > 0, n > 0 or error.
%  b    -- m-by-1 matrix.
%  
%  x    -- Particular solution, Ax=b if feasible, empty o.w.
%  y    -- Empty if feasible, o.w. y'*b < 0 and y'*A = 0.
%  tol  -- norm(A*x-b) < tol*(1+sum(abs(b)))
%          indicates a feasible solution was found.
%
    
    if nargin < 3, tol = 1e-6; end
    
    if size(A,1) ~= size(b,1)
        error('Equations sizes do not match.');
    end
    
    if size(A,2) == 0
        error('spot_eq_solve: empty arguments not supported.');
    end
    
    if size(A,1) == 0
        x = zeros(size(A,2),1);
        y = [];
        return;
    end
    
    
    %-- begin code from qr help
    if issparse(A), 
        R = qr(A); 
    else 
        R = triu(qr(A)); 
    end

    s = warning('error','MATLAB:rankDeficientMatrix');
    warning('error','MATLAB:nearlySingularMatrix');
    try
        x = R\(R'\(A'*b));
        r = b - A*x;
        e = R\(R'\(A'*r));
        x = x + e;
    catch 
        x = pinv(A)*b;
    end
    warning(s);
    %-- end code from qr help
    
    
    if norm(A*x-b) < 1e-10*(1+sum(abs(b)))
        y    = [];
    else
        q = [ -1 ; zeros(size(A,2),1) ];
        Abar = [ b A ]';
        y = Abar\q;
        x = [];
    end
end