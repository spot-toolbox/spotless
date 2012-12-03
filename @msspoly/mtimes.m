function q=mtimes(p1,p2)
p1=msspoly(p1);
p2=msspoly(p2);

if p1.dim(2) ~= p2.dim(1)
    if spot_hasSize(p1,[1 1])
        p1 = diag(repmat(p1,p2.dim(1),1)); %p1 scalar
    elseif spot_hasSize(p2,[1 1])
        p2 = diag(repmat(p2,p1.dim(2),1)); %p2 scalar
    else
        error('Incompatible dimensions.');
    end
end

sz1 = p1.dim;
sz2 = p2.dim;

% Find all indices where p1(.,i)p2(k,.) ~= 0
ik = msspoly.relate(p1.sub(:,2),p2.sub(:,1));

if isempty(ik)
    q = msspoly.zeros(p1.dim(1),p2.dim(2));
else
    i = ik(:,1);
    k = ik(:,2);
    
    q = msspoly([p1.dim(1) p2.dim(2)],...
               [p1.sub(i,1) p2.sub(k,2)],...
               [p1.var(i,:) p2.var(k,:)],...
               [p1.pow(i,:) p2.pow(k,:)],...
               p1.coeff(i,:).*p2.coeff(k,:));
    
end
