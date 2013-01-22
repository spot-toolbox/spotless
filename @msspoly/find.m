function [I,J,S] = find(p)
%  Ind = find(p)
%  [I,J] = find(p)
%  [I,J,V] = find(p)
%
%  p -- n-by-m msspoly.  
%
%  S   -- Indicies of non-zero entries.
%  I,J -- Subscripts of non-zero entries.
%  V   -- corresponding non-zero entries.
%
    
    ind = sub2ind(size(p),p.sub(:,1),p.sub(:,2));
    
    if nargout < 2
        I = unique(ind);
    elseif nargout == 2
        IJ = unique(p.sub,'rows');
        I = IJ(:,1); J = IJ(:,2);
    else
        [~,A,B] = unique(ind);
        nnz = length(A);
        S = msspoly([nnz 1],[B ones(size(p.var,1),1)],p.var,p.pow,p.coeff);
        IJ = p.sub(A,:);
        I = IJ(:,1);
        J = IJ(:,2);
    end    
end