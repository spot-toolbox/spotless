function p=sim_p2pp(q)
% function p=sim_p2pp(q)
%
% convert polynomial to piecewise polynomial

[m,n]=size(q);

p.form='pp';
p.breaks=[0 1];
p.coefs=q;
p.pieces=1;
p.order=size(q,2);
p.dim=size(q,1);
