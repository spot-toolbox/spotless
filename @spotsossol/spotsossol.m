classdef spotsossol < spotsqlsol
    properties 
        gramMatrices = {};
        decompBasis = {};
    end
    
    methods
        function sol = spotsossol(sqlsol,Q,phi)
            sol@spotsqlsol(sqlsol.program,sqlsol.info,...
                           sqlsol.objective,sqlsol.dualObjective,...
                           sqlsol.primalSolution,sqlsol.dualSolution,...
                           sqlsol.dualSlack);
            sol.gramMatrices = Q;
            sol.decompBasis = phi;
        end
        
        function [Q,phi] = getSOSDecomp(sol,token)
            if ~spot_hasSize(token,[1 1]) & ~spot_isIntGe(token,0) ...
                    & token <= length(sol.gramMatrices)
                error('Invalid token for SOS decomposition.');
            end
            
            Q = sol.gramMatrices{token};
            phi = sol.decompBasis{token};
        end
    end
end