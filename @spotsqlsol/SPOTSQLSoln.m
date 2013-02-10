classdef SPOTSQLSoln
    properties 
        program = [];
        user_variables = [];
        primalSolution = [];
        G = [];
        h = [];
        info = [];
    end
    
    methods
        function f = primalInfeasible(sol)
            f = sol.info.pinf;
        end
        function f = dualInfeasible(sol)
            f = sol.info.dinf;
        end
        
        function sol = SPOTSQLSoln(prog,info,uservars,psol)
            sol.program = prog;
            sol.info = info;
            sol.user_variables = uservars;
            sol.primalSolution = psol;
        end
        
        function ev = eval(sol,exp)
            if sol.primalInfeasible
                error('Primal Infeasible, cannot evalutate primal variables.');
            end
            ev = subs(exp,[sol.user_variables],...
                      [sol.primalSolution]);
        end
    end
end