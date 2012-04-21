function U=mono(p)
% function U=mono(p)
%
% U is the column vector of all monomials from p

[x,d]=decomp(p);
if isempty(x),
    U=msspoly(1);
else
    U=recomp(x,d);
end
