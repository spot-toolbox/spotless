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
        primalInfeasible = [];
        dualInfeasible = [];
        prog = [];
        objective = [];
        dualize = 0;
        variables = [];
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
            sol.primalInfeasible = info.primalInfeasible;
            sol.dualInfeasible = info.dualInfeasible;
            sol.objective = objective;
            if nargin < 7, dualize = 0; end
            sol.dualize = dualize;
        end
        
        function err = dimacs(sol)
            [P,A,b,c,K,d] = sol.prog.toSedumi(sol.objective);
            nf = sol.prog.numFree;
            err = spot_sdp_dimacs(A,b,c,K,P'*sol.x,sol.y,...
                                  spot_sdp_cone_symm([ zeros(nf,1) ; P(nf+1:end,nf+1:end)'*sol.z],K));
        end

        function e = eval(sol,expr)
            if sol.primalInfeasible,
                error('Cannot evaluate: primal infeasible.'); 
            end
           if sol.dualize, e = sol.solDualEval(expr);
           else, e = sol.solPrimalEval(expr);
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
            if sol.primalInfeasible,
                error('Cannot evaluate: primal infeasible.'); 
            end

            e = subs(expr,...
                     sol.prog.variables,...
                     sol.prog.decToVar(sol.x));
        end
        
        function e = solDualEval(sol,expr)
            if sol.dualInfeasible,
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