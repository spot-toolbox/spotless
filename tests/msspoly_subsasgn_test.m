nq = 3;
q=msspoly('q',nq);
s=msspoly('s',nq);
c=msspoly('c',nq);
qt=TrigPoly(q,s,c);

J = zeros(3,3) * qt(1)

pos_row_indices = 1:3;

v_or_qdot_indices = zeros(0,1);

Jpos = zeros(3,0);

J(pos_row_indices, v_or_qdot_indices) = Jpos
