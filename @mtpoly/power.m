function q=power(p,n)
% q=mpower(p,n)
%
% p -- k-by-m msspoly matrix.
% n -- non-negative integer.
% 
% Returns q k-by-m with q(i,j) = p(i,j)^n
%
    if msspoly.hasSize(n,[1 1])
        q = p.iter_binary(n,@times);
    elseif msspoly.hasSize(n,size(p))
        p0 = p(:);
        q = msspoly(ones(size(p0)));
        for i = 1:max(n(:))
            ni = n(:) >= i;
            q(ni) = q(ni).*p(ni);
        end
    else
        error('Dimensions do not match.');
    end



end
