function [U2V] = spot_sdp_cone_utproj(K)
%    [V2U,U2V] = spot_sdp_cone_ltproj(K)
%
%    K -- SeDuMi type cone description structure.
%
%    V2U -- matrix sending SeDuMi variable to
%           new vector with each vec'd matrix stored
%           as an upper triangular part, offdiag weighted by 2.
%    U2V -- inverse.
    
    [n,nf,nl,nq,nr,ns] = spot_sdp_cone_dim(K);
    
    V2U = [ speye(nf+nl+nq+nr) ];
    U2V = [ speye(nf+nl+nq+nr) ];
    for i = 1:length(K.s)
        n = K.s(i);
        d = nchoosek(n+1,2);
        r = repmat((1:n)',1,n);
        c = r';
        %V2U = blkdiag(V2U,sparse((1:d)',find(r <= c),ones(d,1)));
        Blk = sparse((1:n^2)',vec(mss_v2s(1:d)),ones(n^2,1));
        Blk(vec(r ~= c),:) =  Blk(vec(r ~= c),:)/sqrt(2);
        U2V = blkdiag(U2V,Blk);
    end
end