function [] = testCstrNamed(varargin)
    names = varargin;
    pr = spotsosprog;
    
    n = length(names);
    zero = msspoly(zeros(n,1));
    obj = msspoly(zeros(n,1));
    tol = zeros(n,1);
    
    for i = 1:n
        fn = str2func(['testCstr/' names{i}]);
        [pr,obj(i),zero(i),tol(i)] = fn(pr);
    end
    sol = pr.minimize(sum(obj));
    [ sol.eval(zero) tol]
    
end