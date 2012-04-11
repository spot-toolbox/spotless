function [data,affine] = sim_data_normalize(data0,domain)
    data = struct();
    affine = struct();
    
    function [Sbar,nrm] = signal_stats(S)
        nrm = struct();
        nrm.mu = (max(S,[],2) + min(S,[],2))/2;

        Sbar = S - repmat(nrm.mu,1,size(S,2));
        
        nrm.T = diag((max(S,[],2) - min(S,[],2))/2);
        Sbar = inv(nrm.T)*Sbar;
    end
    
    [data.Y,affine.Y] = signal_stats(data0.Y);
    [data.U,affine.U] = signal_stats(data0.U);

    if strcmp(domain,'CT')
            [data.X,affine.X] = signal_stats(data0.X);
            affine.V.mu = zeros(size(data0.V,1),1);
            affine.V.T = affine.X.T;
            data.V = inv(affine.V.T)*data0.V;
            tscl = 1;% max(abs(data.V(:)));
            affine.V.T = affine.V.T*tscl;
            data.V = data.V/tscl;
    elseif strcmp(domain,'DT')
            XV = [data0.X data0.V];
            [XVbar,affine.X] = signal_stats(XV);
            affine.V = affine.X;
            data.X = XVbar(:,1:size(data0.X,2));
            data.V = XVbar(:,size(data0.X,2)+(1:size(data0.V,2)));
    else
        error('Unknown domain');
    end

end