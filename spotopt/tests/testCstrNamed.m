function score = testCstrNamed(varargin)
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
    sol = pr.minimize(sum(obj));
    score=double([ sol.eval(zero) tol]);
    
end