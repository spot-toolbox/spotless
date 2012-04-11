% sim_p2f -- Convert an msspoly to a function of doubles in a
%        specified order.
%
% f = sim_p2f(p,x)
%  
% Arguments:
%     p -- A polynomial in x.
%     x -- A free n-by-1 msspoly of unique degree 1 monomials.
%
function f = sim_p2f(p,x)
% TODO Fix
    [xf,pf,Mf]=decomp(msspoly(p));    
    
    if isempty(pf)
        f = @(x) double(p);
        return;
    end
    
    sf = xf.s;
    s = x.s;
    
    [idf,If]= sort(sf(:,3));
    [idx,Ibx]= sort(s(:,3));
    % Test that xf \subset x
    i = 1; j = 1; %id = []; 
    Ib = [];

    while i <= length(idf)
        if j > length(idx)
            error('x must contain all arguments of p');
        end

        if idf(i) == idx(j)
            %            id = [id;idx(j)];
            Ib = [Ib;Ibx(j)];
            i = i+1;
        end
        j = j+1;
    end

    function r = preal(xreal)
        r = reshape(Mf*prod(repmat(xreal(Ib)',size(pf,1),1).^pf,2),size(p));
    end
    f = @preal;
end