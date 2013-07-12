function [x,y,z,info] = spot_sedumi(A,b,c,K,options)
    if nargin < 5, options = spot_sdp_default_options(); end
    
    sedumi_options = struct('fid',options.verbose,'errors',1);
    
    [x,y,s_info] = sedumi(A,b,c,K,sedumi_options);
    z = spot_sdp_cone_symm(c - A'*y,K);
    
    info.primalInfeasible = s_info.pinf;
    info.dualInfeasible   = s_info.dinf;
    info.solverName = 'sedumi';
    info.solverInfo = s_info;
    info.dimacs = spot_sdp_dimacs(A,b,c,K,x,y,z);
end