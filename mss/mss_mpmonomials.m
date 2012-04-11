function q=mss_mpmonomials(v)

ni=length(v);
mi=round(ni)/2;

q=monomials(v{1},v{2});
for i=2:mi,
    h=monomials(v{2*i-1},v{2*i});
    q=q*h';
    q=q(:);
end