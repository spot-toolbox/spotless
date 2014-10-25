function [x,y,z,info] = spot_sedumi(A,b,c,K,options)
    if nargin < 5, options = spot_sdp_default_options(); end
    
    info = struct();
    
    if isfield(options.solver_options,'sedumi')
        sedumi_options = options.solver_options.sedumi;
    else
        sedumi_options = struct('errors',0);
    end    
    if ~isfield(sedumi_options,'fid')
        sedumi_options.fid = options.verbose;
    end
    
    start = spot_now();
    [x,y,s_info] = sedumi(A,b,c,K,sedumi_options);
    [info.ctime,info.wtime] = spot_etime(spot_now(),start);
    
    z = spot_sdp_cone_symm(c - A'*y,K);
    
    if s_info.pinf & s_info.dinf
        status = spotsolstatus.STATUS_PRIMAL_AND_DUAL_INFEASIBLE;
    elseif s_info.pinf
        status = spotsolstatus.STATUS_PRIMAL_INFEASIBLE;
    elseif s_info.dinf
        status = spotsolstatus.STATUS_DUAL_INFEASIBLE;
    else
        status = spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE;
    end
    
    if s_info.numerr > 0
        status = spotsolstatus.STATUS_NUMERICAL_PROBLEMS;
    end

    info.status = status;
    info.solverName = 'sedumi';
    info.solverInfo = s_info;
    %info.dimacs = spot_sdp_dimacs(A,b,c,K,x,y,z);
end
