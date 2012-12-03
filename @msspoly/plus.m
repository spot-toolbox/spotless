function p=plus(p1,p2)
p1=msspoly(p1);
p2=msspoly(p2);

if ~spot_hasSize(p1,size(p2))
    if spot_hasSize(p1,[1 1])
        p1 = repmat(p1,size(p2));
    elseif spot_hasSize(p2,[1 1])
        p2 = repmat(p2,size(p1));
    else
        error('incompatible dimensions');
    end
end

[var1,var2]=msspoly.padZeros(p1.var,p2.var);
[pow1,pow2]=msspoly.padZeros(p1.pow,p2.pow);

p = msspoly(p1.dim,...
           [p1.sub;p2.sub],...
           [var1;var2],...
           [pow1;pow2],...
           [p1.coeff;p2.coeff]); % 
end
