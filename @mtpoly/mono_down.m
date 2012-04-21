function U=mono_down(p)
% function U=mono_down(p)
%
% U is the column vector of all monomials dominated by p

[x,d]=decomp(p);
if isempty(x),
    U=msspoly(1);
else
    d=mint_ch(mint_down(d));
    U=recomp(x,d);
end
