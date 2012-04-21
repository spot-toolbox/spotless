function U=mono_root(p)
% function U=mono_root(p)
%
% U is the column vector of all monomials dominated by sqrt(p)

[x,d]=decomp(p);
if isempty(x),
    U=msspoly(1);
else
    d=mint_ch(mint_down(d),2);
    U=recomp(x,d);
end
