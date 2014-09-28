% TODO:
%   Really we want a solution to encode the entire process of
%   producing the solution (i.e. what program was /actually/ solved
%   and the series of linear transformations that get us back to
%   this problem).
%
%
classdef spotprogsol
    properties
        x = [];
        y = [];
        z = [];
        info = struct();
        status = spotsolstatus.STATUS_UNSOLVED
        prog = [];
        objective = [];
        dualize = 0;
        variables = [];
        gramMatrices = {};
        gramMonomials = {};
    end

    methods (Static)
         function feas = statusIsPrimalFeasible(status)
            feas = (status == spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE ...
                    | status == spotsolstatus.STATUS_DUAL_INFEASIBLE | status == spotsolstatus.STATUS_NUMERICAL_PROBLEMS);
        end

        function feas = statusIsDualFeasible(status)
            feas = (status == spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE ...
                    | status == spotsolstatus.STATUS_PRIMAL_INFEASIBLE | status == spotsolstatus.STATUS_NUMERICAL_PROBLEMS);
        end
    end

    methods
        function sol = spotprogsol(prog,objective,x,y,z,info,dualize)
            if ~prog.isPrimalWithFree()
                error(['Solutions must come from primal form with ' ...
                       'free variables right now.']);
            end
            sol.prog = prog;
            sol.x = x;
            sol.y = y;
            sol.z = z;
            sol.info = info;
            sol.status = info.status;
            sol.objective = objective;
            if nargin < 7, dualize = 0; end
            sol.dualize = dualize;
            
            sol.gramMatrices = prog.gramMatrices;
            sol.gramMonomials = prog.gramMonomials;
            
        end

        function feas = isPrimalFeasible(sol)
            feas = spotprogsol.statusIsPrimalFeasible(sol.status);
        end

        function feas = isDualFeasible(sol)
            feas = spotprogsol.statusIsDualFeasible(sol.status);
        end

        function err = dimacs(sol)
            [P,A,b,c,K,d] = sol.prog.toSedumi(sol.objective);
            nf = sol.prog.numFree;
            err = spot_sdp_dimacs(A,b,c,K,P'*sol.x,sol.y,...
                                  spot_sdp_cone_symm([ zeros(nf,1) ; P(nf+1:end,nf+1:end)'*sol.z],K));
        end

        function e = eval(sol,expr)
            if ~sol.isPrimalFeasible(),
                error('Cannot evaluate: primal infeasible.'); 
            end
            if sol.dualize,
                e = sol.solDualEval(expr);
            else, 
                e = sol.solPrimalEval(expr);
           end
        end

        function e = dualEval(sol,expr)
            if sol.dualize, e = sol.solPrimalEval(expr);
            else, e = sol.solDualEval(expr);
            end
        end
    end
    methods (Access = private)
        function e = solPrimalEval(sol,expr)
            if ~sol.isPrimalFeasible(),
                error('Cannot evaluate: primal infeasible.'); 
            end

            e = subs(expr,...
                     sol.prog.variables,...
                     sol.prog.decToVar(sol.x));
        end

        function e = solDualEval(sol,expr)
            if ~sol.isDualFeasible,
                error('Cannot evaluate: dual infeasible.'); 
            end

            nf = sol.prog.numFree;
            e = subs(expr,...
                     [sol.prog.variables
                      sol.prog.dualEqVariables],...
                     [sol.prog.decToVar([zeros(nf,1);sol.z])
                      sol.y]);
        end
    end
end