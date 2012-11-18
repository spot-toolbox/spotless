function q=indexinto(p,varargin)
% Grab linear indices corresponding to this subsref.
    pp = 1:prod(p.dim);

    Q=subsref(reshape(pp,p.dim),struct('type','()','subs',{varargin}));

    qq=Q(:);
    % subscripts for the values which come form this matrix.
    [qi,qj] = ind2sub(size(Q),1:length(qq));
    
    pind = sub2ind(p.dim,p.sub(:,1),p.sub(:,2));
    
    if isempty(pind) | isempty(Q)
        q = msspoly.zeros(size(Q));
    else
        
        rel = msspoly.relate(qq,pind); % all pairs (i,k) s.t. qq(i) = pind(k)
        
        q = msspoly(size(Q),...
                   [reshape(qi(rel(:,1)),[],1) reshape(qj(rel(:,1)),[],1)],...
                   p.var(rel(:,2),:),...
                   p.pow(rel(:,2),:),...
                   p.coeff(rel(:,2),:));
    end
end
