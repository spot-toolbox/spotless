classdef spotsossol < spotsqlsol
    properties 
        gramMatrices = {};
        decompBasis = {};
        decompGramBasis = {};
        decompDual = {};
    end
    
    methods
        function sol = spotsossol(sqlsol,Q,phi,y,basis)
            sol@spotsqlsol(sqlsol.program,sqlsol.info,...
                           sqlsol.objective,sqlsol.dualObjective,...
                           sqlsol.primalSolution,sqlsol.dualSolution,...
                           sqlsol.dualSlack);
            sol.gramMatrices = Q;
            sol.decompGramBasis = phi;
            sol.decompBasis = basis;
            sol.decompDual = y;
        end
        
        function [Q,phi,y,basis] = getSOSDecomp(sol,token)
            if ~spot_hasSize(token,[1 1]) & ~spot_isIntGe(token,0) ...
                    & token <= length(sol.gramMatrices)
                error('Invalid token for SOS decomposition.');
            end
            
            Q = sol.gramMatrices{token};
            phi = sol.decompGramBasis{token};
            basis = sol.decompBasis{token};
            y = sol.decompDual{token};
        end
    end
end