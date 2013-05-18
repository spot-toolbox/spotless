function errors = spot_sdp_dimacs(A,b,c,K,x,y,z)
    if nargin < 5, x = []; end
    if nargin < 6, y = []; end
    errors = NaN*ones(6,1);
    
    function lm = lambda_min(x,K)
        lm = Inf;
        
        off = K.f;
        
        if K.l > 0
            lm = min(x(1:K.l));
        end

        off = off + K.l;
        
        for i = 1:length(K.q)
            n = K.q(i);
            if n == 0, continue; end
            lm = min(x(off+1)^2-sum(x(off+(2:n)).^2),lm);
            off = off + K.q(i);
        end
        
        for i = 1:length(K.r)
            n = K.r(i);
            if n == 0, continue; end
            lm = min(2*x(off+1)*x(off+2)-sum(x(off+(3:n)).^2),lm);
            off = off + K.r(i);
        end
        
        for i = 1:length(K.s)
            n = K.s(i);
            if n == 0, continue; end            
            lm = min(min(eig(full(reshape(x(off+(1:n^2)),n,n)))),lm);
            off = off + n^2;
        end
    end
    
    if ~isempty(x)
        denom = 1+sum(abs(b));
        errors(1) = norm(A*x-b)/denom;
        errors(2) = max(0,-lambda_min(x,K)/denom);
    end
    
    if ~isempty(y)
        if nargin < 7,
            z = c-A'*y;
        end
        denom = 1+sum(abs(c));
        errors(3) = norm(z-(c-A'*y))/denom;
        errors(4) = max(0,-lambda_min(z,K)/denom);

        if ~isempty(x) 
            denom=(1 + abs(c'*x) + abs(b'*y));
            errors(5) = (c'*x - b'*y)/denom;
            
            if 0 == errors(4) && 0 == errors(2)
                errors(6) = x'*z/denom;
            end
        end
    end
end