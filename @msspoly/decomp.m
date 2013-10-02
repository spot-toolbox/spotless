function [x,p,M,sz]=decomp(q,v)
%
%  [x,p,M,sz]=decomp(q)
%
%  q -- An n-by-k msspoly.
%
%  Returns:
%  x  -- an v-by-1 free msspoly including all variables appearing in q.
%  p  -- an m-by-v array of non-negative integers.
%  M  -- an (nk)-by-m sparse array of doubles.
%  sz -- size(q)
%
%  Satisfying:
%
%  (*)    q(i) = sum_k M(i,k) prod_j x(j)^p(k,j)
%
%  [x,p,M]=decomp(q,y)
%  
%  y -- l-by-1 free msspoly.
%
%  Similar output, satisfying (*) except:
%  
%   x -- now includes no elements of y.
%   M -- now an msspoly that can depend on y.
%
%


    function [vary,powy] = subvp(num,var,pow)
    nn = msspoly.match_list(num,var);
    [i,j,s] = find(nn);
    vary = sparse(i,j,num(s),size(var,1),size(var,2));
    powy = sparse(i,j,pow(sub2ind(size(pow),i,j)));
    end

sz = size(q);

if length(q) == 0
    x = msspoly();
    p = [];
    M = [];
elseif nargin < 2
    % First flatten the matrix.
    sub = q.sub;
    ind = sub2ind(q.dim,sub(:,1),sub(:,2));
    nk = prod(q.dim);
    
    % Next, build the list of unique monomials.
    varpow = [q.var q.pow];
    [~,J,K] = unique(varpow,'rows');
    uvar = q.var(J,:);
    upow = q.pow(J,:);
    
    % K maps row number to corresponding monomial.
    m = length(J); % number of unique monomials.    
    M = sparse(ind,K,q.coeff,nk,m); 
    
    % Next construct the full list of variable id numbers.
    xn = unique(q.var(q.var(:)~=0));  % All the ID numbers.
    xn = reshape(xn,[],1);
    v = length(xn);                   % Count.
    
    varn = msspoly.match_list(xn,uvar);     % Replace entries of var with
                                            % indices into xn.
    
    vind = find(varn ~= 0); % index of location of variables.
    [i,~] = ind2sub(size(varn),vind); % find row number
                                      % (corresponds to unique monomial).
    p = sparse(i,varn(vind),upow(vind),m,v);
    x  = msspoly(size(xn),[(1:v)' ones(v,1)],xn,ones(v,1),ones(v,1)); 
    
    o1 = x;
    o2 = p;
    o3 = M;
else
    [x,p,M] = decomp(q);
    
    mtch = match(x,v);
    v = indexinto(v,mtch~=0);
    mtch = mtch(mtch ~= 0);
    
    if ~isempty(mtch)
        d = recomp(v,p(:,mtch),speye(size(p,1)));
        M = M.*repmat(d',size(M,1),1);
        
        b = 1:size(p,2);
        b(mtch) = [];
        
        x = indexinto(x,b');
        [p,~,J] = unique(p(:,b),'rows');
        M = M*sparse(1:size(M,2),J,ones(size(M,2),1),size(M,2),size(p,1));
    end
end


x = msspoly(x);
end

