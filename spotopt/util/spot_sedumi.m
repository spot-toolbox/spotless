function [x,y,z,info] = spot_sedumi(A,b,c,K)
    [x,y,s_info] = sedumi(A,b,c,K);
    z = c - A'*y;
    
    info.primalInfeasible = s_info.pinf;
    info.dualInfeasible   = s_info.dinf;
    info.solverName = 'sedumi';
    info.solverInfo = s_info;
    info.dimacs = spot_sdp_dimacs(A,b,c,K,x,y,z);
end