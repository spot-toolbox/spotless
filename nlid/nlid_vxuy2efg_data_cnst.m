% 
%
%  prog = nlid_vxu2ef_data_cnst(prog,vars,Data,Hdata,rindx,rdata)
%
%  Helper function for adding per-data-point constraints.
%
%  prog  -- mssprog to add constriants to
%  vars  -- n-by-1 msspoly of variables to substitute.
%  Data  -- n-by-N matrix of values to substitute.
%  Hdata -- k-by-k matrix of polynomials.  This will be required to
%           be positive definite at each data point.
%  rindx -- 1-by-1 index into Hdata where the slack variable r
%           should be added.
%  rdata -- 1-by-N msspoly of slack variables for each data-point.
%
function prog = nlid_vxu2ef_data_cnst(prog,vars,Data,Hdata,rindx,rdata)
    D = msubs(Hdata(:),vars,Data);
    D(rindx,:) =  D(rindx,:) + rdata';
    
    if deg(D(:,1)) < 2
        N = nchoosek(size(Hdata,1),2)+size(Hdata,1);
        Q = msspoly('Q',N*size(D,2));
        Q = reshape(Q,N,size(D,2));
    
        prog.psd = Q;
    
        ii = mss_s2v(reshape(1:size(Hdata,1)*size(Hdata,1),size(Hdata,1),size(Hdata,1)));
        prog.eq = Q - D(ii(:),:);
    else
        prog.sos = D;
    end
end
