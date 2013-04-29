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
        
        variables = [];
    end
    
    methods
        function sol = spotprogsol(prog,x,y,z,info)
            if ~prog.isPrimalWithFree()
                error(['Solutions must come from primal form with ' ...
                       'free variables right now.']);
            end
            sol.prog = prog;
            sol.x = x;
            sol.y = y;
            sol.z = z;
            sol.info = info;
        end
        
        function e = eval(sol,expr)
            e = subs(expr,sol.prog.variables,sol.prog.decToVar(sol.x));
        end
    end
end