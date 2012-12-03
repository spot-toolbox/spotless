classdef spotsqlsol
    properties 
        program = [];
        primalSolution = [];
        dualSolution = [];
        dualSlack = [];
        objective = [];
        info = [];
    end
    
    methods
        function f = primalInfeasible(sol)
            f = sol.info.pinf;
        end
        function f = dualInfeasible(sol)
            f = sol.info.dinf;
        end
        
        function sol = spotsqlsol(prog,info,obj,psol,dsol,dslack)
            sol.program = prog;
            sol.info = info;
            sol.primalSolution = psol;
            sol.dualSolution = dsol;
            sol.dualSlack = dslack;
            sol.objective = obj;
        end
        
        function ev = eval(sol,exp)
            if sol.primalInfeasible
                error('Primal Infeasible, cannot evalutate primal variables.');
            end
            ev = double(subs(exp,sol.program.variables,sol.primalSolution));
        end
        function ev = dualEval(sol,exp)
            if sol.dualInfeasible
                error('Dual Infeasible, cannot evaluate dual variables.');
            end
            ev = double(subs(exp,sol.program.variables,sol.dualSlack));
        end
    end
end