function q=mss_mpmonomials(v)

ni=length(v);
mi=round(ni)/2;

q = [];
for i = 1:mi
    var = v{2*i-1}; exp = v{2*i};
    if ~isempty(var)
        h = monomials(var,exp);
        if isempty(q)
            q = h;
        else
            q = q*h.';
            q = q(:);
        end
    end
end

end