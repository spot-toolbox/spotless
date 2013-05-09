function [score,sol,dsol] = testCstrNamed(varargin)
    solver = @spot_sedumi;
    names = varargin;
    pr = spotsosprog;
    
    n = length(names);
    zero = msspoly(zeros(n,1));
    obj = msspoly(zeros(n,1));
    tol = zeros(n,1);
    
    for i = 1:n
        fn = str2func([names{i}]);
        [pr,obj(i),zero(i),tol(i)] = fn(pr);
    end
    
    pobj = sum(obj);
    sol = pr.minimize(pobj,solver);
    
    dsol = pr.minimize(pobj,solver,struct('dualize',1));
    
    % [dl,dobj] = pr.toDual(pobj);
    % dsol = dl.minimize(-dobj,solver,struct('dualize',1));

    score=double([ sol.eval(zero) dsol.eval(zero) tol]);
    
end