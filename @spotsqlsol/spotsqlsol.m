classdef spotsqlsol
    properties 
        program = [];
        primalSolution = [];
        dualSolution = [];
        dualSlack = [];
        objective = [];
        dualObjective = [];
        info = [];
    end
    
    methods
        function f = primalInfeasible(sol)
            f = sol.info.pinf;
        end
        function f = dualInfeasible(sol)
            f = sol.info.dinf;
        end
        
        function sol = spotsqlsol(prog,info,obj,dobj,psol,dsol,dslack)
            sol.program = prog;
            sol.info = info;
            sol.primalSolution = psol;
            sol.dualSolution = dsol;
            sol.dualSlack = dslack;
            sol.objective = obj;
            sol.dualObjective = dobj;
        end
        
        function ev = eval(sol,exp)
            if sol.primalInfeasible
                error('Primal Infeasible, cannot evalutate primal variables.');
            end
            ev = subs(exp,[sol.program.variables],...
                      [sol.primalSolution]);
        end
        function ev = dualEval(sol,exp)
            if sol.dualInfeasible
                error('Dual Infeasible, cannot evalutate dual variables.');
            end
            ev = subs(exp,[sol.program.dualVariables],...
                      [sol.dualSolution]);
        end
        function ev = dualSlackEval(sol,exp)
            if sol.dualInfeasible
                error('Dual Infeasible, cannot evaluate dual variables.');
            end
            ev = subs(exp,sol.program.variables,sol.dualSlack);
        end
    end
end