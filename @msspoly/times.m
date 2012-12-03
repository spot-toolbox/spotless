function q=times(p1,p2)

p1=msspoly(p1);
p2=msspoly(p2);

sz = size(p2);
if ~spot_hasSize(p1,sz)
    error('Incorrect dimensions.');
end
  
L1 = sub2ind(sz,p1.sub(:,1),p1.sub(:,2));
L2 = sub2ind(sz,p2.sub(:,1),p2.sub(:,2));

% Build all pairs of indices into sparse vector which
% correspond to same entry of the matrix.
ik = msspoly.relate(L1,L2);

if isempty(ik)
    q = msspoly.zeros(sz);
else
    i = ik(:,1);
    k = ik(:,2);
    q = msspoly(sz,...
               p1.sub(i,:),...
               [p1.var(i,:) p2.var(k,:)],...
               [p1.pow(i,:) p2.pow(k,:)],...
               p1.coeff(i,:).*p2.coeff(k,:));
end

end
