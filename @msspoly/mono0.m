function U=mono0(p)
% function U=mono0(p)
%
% U is the column vectors of all monomials from p, excluding 1 (if present)

[x,d]=decomp(p);
d=d(any(d~=0,2),:);
U=recomp(x,d);
