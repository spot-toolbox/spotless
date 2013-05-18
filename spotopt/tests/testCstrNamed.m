function [score,sol,dsol,err] = testCstrNamed(varargin)
    solver = @spot_mosek;
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
    
    score=double([ sol.eval(zero) tol]);
    % opt = spot_sdp_default_options();
    % opt.dualize = 1;
    % dsol = pr.minimize(pobj,solver,opt);

    %score=double([ sol.eval(zero) dsol.eval(zero) tol]);
    
end