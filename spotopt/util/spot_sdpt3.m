function [x,y,z,info] = spot_sedumi(A,b,c,K,options)
    if nargin < 5, options = spot_sdp_default_options(); end
    
    sdpt3_options = struct('printlevel', 3*options.verbose);
    
    [blk,At,C,b] = read_sedumi(A,b,c,K);
    [obj,X,y,Z,sdpt3info] = sqlp(blk,At,C,b,sdpt3_options);
    x = [];
    for i = 1:size(blk,1)
        dim = blk{i,2};
        switch blk{i,1}
          case 's',
            xi = zeros(sum(dim.^2),1);
            m = size(X{i},1);
            for k = 1:length(dim)
                n = dim(k);
                fromOff = (m+1)*sum(dim(1:k-1));
                iFrom = fromOff + ...
                        repmat((1:n)',1,n) + ...
                        repmat((0:n-1)*m,n,1);
                
                toOff = sum(dim(1:k-1).^2);
                iTo = toOff + (1:n^2)';
                
                xi(iTo) = X{i}(iFrom);
            end
          otherwise,
            xi = X{i}(:);
        end
        x = [ x ; xi];
    end
    
    z = spot_sdp_cone_symm(c - A'*y,K);
    
    info.primalInfeasible = sdpt3info.termcode == 1;
    info.dualInfeasible   = sdpt3info.termcode == 2;
    info.solverName = 'sdpt3';
    info.solverInfo = sdpt3info;
    info.dimacs = spot_sdp_dimacs(A,b,c,K,x,y,z);
end
    