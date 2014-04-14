function D = spot_vecBlkDiag(V)
%
% D = spot_vecBlkDiag(V)
%
%  V -- d-by-N matrix where d = (n+1)*n/2.
%
%  Creates a block diagonal matrix with the i-th block
%  being mss_v2s(V(:,i))
    
    d = size(V,1);
    N = size(V,2);
    [n,err] = spot_psdNoToDim(d);
    if err, error('First dimension not reasonable.'); end
    
    I = repmat((1:n)',1,n);
    r = repmat(mss_s2v(I),N,1);
    c = repmat(mss_s2v(I'),N,1);
    b = reshape(repmat(1:N,d,1),[],1);
    
    offdiag = r~=c;
    
    i = [ r + (b-1)*n
          c(offdiag) + (b(offdiag)-1)*n ];
    j = [ c + (b-1)*n
          r(offdiag) + (b(offdiag)-1)*n ];
    
    D = sparse(i,j,[V(:) ; V(offdiag)]);
    D = (D+D')/2;
end