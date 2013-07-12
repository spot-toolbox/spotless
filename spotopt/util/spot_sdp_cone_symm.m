function xprime = spot_sdp_cone_symm(x,K);
%  Utility for handling SDP formats.
%  
%  xprime = spot_sdp_cone_symm(x,K);
    [n,nf,nl,nq,nr,ns] = spot_sdp_cone_dim(K);
    
    if size(x,1) ~= n, 
        error('x must have appropriate number of rows');
    end
    
    if issparse(x),
        xprime = sparse(size(x,1),size(x,2));
    else
        xprime = zeros(size(x));
    end
    
    nid = nf+nl+nq+nr;
    xprime(1:nid,:) = x(1:nid,:);
    
    if ns == 0, return; end
    
    sMap = containers.Map('KeyType','uint32','ValueType','any');
    sz = unique(K.s);
    for i = 1:length(sz)
        n = sz(i);
        T = reshape(1:n^2,n,n)';
        I = speye(n^2);
        sMap(n) = (I + sparse(T(:),(1:n^2),ones(n^2,1)))/2;
    end
    
    maps = cell(length(K.s),1);
    for i = 1:length(K.s)
        maps{i} = sMap(K.s(i));
    end
    T = blkdiag(maps{:});
    xprime(nid+1:end,:) = T*x(nid+1:end,:);
end