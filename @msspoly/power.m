function q=power(p,n)
% q=mpower(p,n)
%
% p -- k-by-m msspoly matrix.
% n -- non-negative integer.
% 
% Returns q k-by-m with q(i,j) = p(i,j)^n
%
    if spot_hasSize(n,[1 1])
        q = p.iter_binary(n,@times);
    elseif spot_hasSize(p,[1 1])
        q = repmat(p,size(n,1),size(n,2)).^n;
    elseif spot_hasSize(n,size(p))
        q = msspoly(ones(size(p)));
        
        for i = 1:prod(size(p))
            q = assign(q,power(indexinto(p,i),n(i)),i);
        end
    else
        error('Dimensions do not match.');
    end



end
