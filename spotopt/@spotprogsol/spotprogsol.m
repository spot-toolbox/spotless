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
        prog = [];
        objective = [];
        
        variables = [];
    end
    
    methods
        function sol = spotprogsol(prog,objective,x,y,z,info)
            if ~prog.isPrimalWithFree()
                error(['Solutions must come from primal form with ' ...
                       'free variables right now.']);
            end
            sol.prog = prog;
            sol.x = x;
            sol.y = y;
            sol.z = z;
            sol.info = info;
            sol.objective = objective;
        end
        
        function err = dimacs(sol)
            [P,A,b,c,K,d] = sol.prog.toSedumi(sol.objective);
            err = spot_sdp_dimacs(A,b,c,K,P'*sol.x,sol.y,P'*sol.z);
        end
        
        function e = eval(sol,expr)
            e = subs(expr,...
                     sol.prog.variables,...
                     sol.prog.decToVar(sol.x));
        end
        
        function e = dualEval(sol,expr)
            nf = sol.prog.numFree;
            e = subs(expr,...
                     [sol.prog.variables
                      sol.prog.dualEqVariables],...
                     [sol.prog.decToVar([zeros(nf,1);sol.z])
                      sol.y]);
        end
    end
end