function y=isequal(a,b)

if ~isa(a,'msspoly') | ~isa(b,'msspoly') 
    y = 0;
elseif ~msspoly.hasSize(a,size(b)) 
    y = 0;
elseif isempty(a.coeff) & isempty(b.coeff)
    y = 1;
elseif ~msspoly.hasSize(a.coeff,size(b.coeff))
    y = 0;
elseif isempty(a.var)
    if ~isempty(b.var), y = 0;
    else
        y =  (size(a.coeff,1) == size(b.coeff,1)) & ...
             all(all([a.coeff] ==...
                     [b.coeff]));
    end
elseif ~msspoly.hasSize(a.var,size(b.var))
    y = 0;
else
    y =  (size(a.coeff,1) == size(b.coeff,1)) & ...
         all(all([a.sub a.var a.pow a.coeff] ==...
                 [b.sub b.var b.pow b.coeff]));
end

end
