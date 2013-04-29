function [feas,E,F,g,U,V,w,Ad,bd,cd,Kd,dd] = spot_sdp_remove_redundant_eqs(A,b,c,K,d)
    if nargin < 5, d = 0 ; end
    
    n = size(A,2);
    m = size(A,1);
    E = speye(n); F = speye(n); g = sparse(n,1);
    cd = c; dd = d; Kd = K;
    U = speye(m);
    
    if isempty(A)
        Ad = [];
        bd = [];
        feas = 0;
        return;
    end
    
    [x,y] = spot_eq_solve(A,b);
    
    if isempty(x),
        feas = 1;
        w = y;
        V = [];
        Ad = A;
        bd = b;
    else
        feas = 0;
    
        if issparse(A)
            R = qr(A');
        else
            [~,R] = qr(A');
        end
        [I,J] = find(R);
        keep = accumarray(I,J,[max(I) 1],@min,0);
        tossCnt = size(A,1) - length(keep);
        
        w = sparse(size(A,1),1);
        V = speye(size(A,1));
        
        V = V(:,keep);
        Ad = A(keep,:);
        bd = b(keep);
    end
    
end