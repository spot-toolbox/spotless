 function [x,y,z,info] = spot_sedumi_fr(A,b,c,K,options)
    if nargin < 5, options = spot_sdp_default_options(); end
    
    info = struct(); 
    
    sedumi_options = struct('fid',options.verbose,'errors',1);
    
    [A,c] = makeSymmetric(A,c,K);
    prg = frlibPrg(A,b,c,K);
    prgR = prg.ReducePrimal('dd');
    
    
    start = spot_now();
    [xFr,~,s_info] = prgR.Solve();
    % [x,y,s_info] = sedumi(A,b,c,K,sedumi_options);
    [info.ctime,info.wtime] = spot_etime(spot_now(),start);
    x = prgR.RecoverPrimal(xFr); 
    
    y = [];
    z = []; %spot_sdp_cone_symm(c - A'*y,K);
    
    info.primalInfeasible = s_info.pinf;
    info.dualInfeasible   = s_info.dinf;
    info.solverName = 'sedumi';
    info.solverInfo = s_info;
    % info.dimacs = spot_sdp_dimacs(A,b,c,K,x,y,z);
end